"""
pipeline.py — Core processing classes for ZTF-Tools

AstroTools: Full image processing pipeline (alignment, stacking,
            subtraction, photometry, transient detection)
aphoto_modified: Aperture photometry utilities

Author: Florian Devender-Dauge
Co-developer (image alignment & reference pipeline): Quentin Arvois
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import shutil
import glob
import warnings
import math

from astropy.io import fits
from astropy.wcs import WCS
from astropy import table, coordinates, units, wcs, time
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
from astropy.stats import SigmaClip, sigma_clipped_stats
from astropy.time import Time
from astropy.utils.exceptions import AstropyWarning

from reproject import reproject_interp
from scipy.stats import norm
from scipy.ndimage import gaussian_filter, zoom, median_filter
from scipy.ndimage import zoom as scipy_zoom
from datetime import datetime, timedelta
from matplotlib.patches import Rectangle

from photutils.background import Background2D, MedianBackground
from photutils.detection import DAOStarFinder


# ---------------------------------------------------------------------------
# Aperture photometry (legacy class, kept for compatibility)
# ---------------------------------------------------------------------------

class aphoto_modified:
    """Simple aperture photometry on a single FITS image.

    Parameters
    ----------
    fdata_r : HDUList
        Opened FITS file (astropy).
    """

    def __init__(self, fdata_r):
        img_r = fdata_r[0].data
        self.img_r_bkg = img_r - np.nanmedian(img_r)
        self.magzp_r = fdata_r[0].header['MAGZP']
        self.wcs_img_r = wcs.WCS(header=fdata_r[0].header)

    def star_photometry(self, x0, y0, radius):
        """Compute aperture flux at pixel position (x0, y0)."""
        sky = self.wcs_img_r.pixel_to_world(x0, y0)
        xr, yr = self.wcs_img_r.world_to_pixel(sky)
        x0r = round(xr.mean())
        y0r = round(yr.mean())
        flux_r = 0
        for xi in range(x0r - radius, x0r + radius):
            for yi in range(y0r - radius, y0r + radius):
                ri = np.sqrt((xi - x0r) ** 2 + (yi - y0r) ** 2)
                if ri < radius:
                    flux_r += self.img_r_bkg[yi, xi]
        flux_r = max(flux_r, 1)
        mag_r = -2.5 * np.log10(flux_r) + self.magzp_r
        ra = sky.ra.degree
        dec = sky.dec.degree
        return (ra, dec, flux_r, mag_r)

    def photometry(self, list_star):
        """Run aperture photometry on a list of (x, y, radius) tuples."""
        ra, dec, flux_r, mag_r = [], [], [], []
        for x, y, r in list_star:
            data = self.star_photometry(x, y, r)
            ra.append(data[0])
            dec.append(data[1])
            flux_r.append(data[2])
            mag_r.append(data[3])
        return Table([ra, dec, flux_r, mag_r], names=('ra', 'dec', 'flux_r', 'mag_r'))


# ---------------------------------------------------------------------------
# Main pipeline class
# ---------------------------------------------------------------------------

class AstroTools:
    """Full ZTF image processing pipeline.

    Handles alignment, reference image construction, temporal stacking,
    image subtraction, automatic transient detection, and light curve
    extraction.

    Parameters
    ----------
    fits_files : list of str
        Sorted list of paths to raw FITS science images.
    file_name : str
        Identifier string used for output directory naming.
        Typical format: '000664_zr_c06_o_q1'

    Examples
    --------
    >>> tool = AstroTools(fits_files, file_name="000664_zr_c06_o_q1")
    >>> tool.distrib_seeing()
    >>> tool.align_images()
    >>> tool.get_refimg(aligned_files, ref_method="median", seeing_max=6)
    """

    def __init__(self, fits_files, file_name='00Field_Filter_cCcd_o_qQuadrant'):
        self.fits_files = fits_files
        self.first_file = fits_files[0]
        self.first_hdu = fits.open(self.first_file)[0]
        self.first_wcs = WCS(self.first_hdu.header)
        self.first_shape = self.first_hdu.data.shape
        self.file_name = file_name

    # ------------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------------

    def show_img(self, data, title='', vmin=None, vmax=None, cmap=None):
        """Display a 2D image with colorbar."""
        plt.figure(figsize=(7, 7))
        im = plt.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax)
        plt.colorbar(im, fraction=0.046, pad=0.04).set_label('Flux (ADU)', fontsize=14)
        plt.title(title, fontsize=16, pad=20)
        plt.xlabel('X pixel', fontsize=14)
        plt.ylabel('Y pixel', fontsize=14)
        plt.tight_layout()
        plt.show()

    def radec_to_pixel(self, ra, dec, image_wcs):
        """Convert sky coordinates (RA, Dec) to pixel coordinates."""
        skycoord = SkyCoord(ra, dec, unit='deg')
        x, y = image_wcs.world_to_pixel(skycoord)
        return x, y

    def pixel_to_radec(self, x, y, image_wcs):
        """Convert pixel coordinates to sky coordinates (RA, Dec)."""
        skycoord = image_wcs.pixel_to_world(x, y)
        return skycoord.ra.degree, skycoord.dec.degree

    def get_zoom(self, image, ra, dec, image_wcs, zoom=50):
        """Return pixel bounding box for a zoom region centred on (RA, Dec)."""
        x, y = self.radec_to_pixel(ra, dec, image_wcs)
        x_min = max(int(x - zoom), 0)
        x_max = min(int(x + zoom), image.shape[1])
        y_min = max(int(y - zoom), 0)
        y_max = min(int(y + zoom), image.shape[0])
        return x_min, x_max, y_min, y_max

    def reproject_astropy(self, source_data, source_wcs, target_wcs, target_shape):
        """Reproject source image to target WCS using bilinear interpolation."""
        array, _ = reproject_interp(
            (source_data, source_wcs), target_wcs,
            shape_out=target_shape, order='bilinear'
        )
        return array

    def degrade_psf(self, image, seeing_from, seeing_to, pixel_scale=1.01):
        """Homogenize PSF by applying a Gaussian blur.

        Parameters
        ----------
        seeing_from : float
            Current image seeing (arcsec).
        seeing_to : float
            Target seeing (arcsec, must be >= seeing_from).
        pixel_scale : float
            Arcsec per pixel (ZTF default: 1.01).
        """
        if seeing_to <= seeing_from:
            return image
        fwhm_diff = np.sqrt(seeing_to ** 2 - seeing_from ** 2)
        sigma_pix = fwhm_diff / (2.3548 * pixel_scale)
        return gaussian_filter(image, sigma=sigma_pix)

    # ------------------------------------------------------------------
    # Background subtraction
    # ------------------------------------------------------------------

    def remove_background_photutils(self, image, box_size=64, filter_size=(3, 3)):
        """Remove sky background using photutils 2D background estimator."""
        sigma_clip = SigmaClip(sigma=3.0)
        bkg_estimator = MedianBackground()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            bkg = Background2D(image, box_size=box_size, filter_size=filter_size,
                               sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
        return image - bkg.background, bkg.background

    def remove_local_background(self, image, block_size=64):
        """Remove sky background using a local median block map."""
        ny, nx = image.shape
        ny_b, nx_b = ny // block_size, nx // block_size
        bkg_map = np.zeros((ny_b, nx_b))
        for i in range(ny_b):
            for j in range(nx_b):
                patch = image[i * block_size:(i+1) * block_size,
                              j * block_size:(j+1) * block_size]
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    bkg_map[i, j] = np.nanmedian(patch)
        zy = ny / bkg_map.shape[0]
        zx = nx / bkg_map.shape[1]
        bkg_full = scipy_zoom(bkg_map, (zy, zx), order=1)
        return image - bkg_full, bkg_full

    # ------------------------------------------------------------------
    # Seeing statistics
    # ------------------------------------------------------------------

    def distrib_seeing(self, fits_file=None, plot=True):
        """Compute and optionally plot the seeing distribution.

        Parameters
        ----------
        fits_file : list of str, optional
            Override default file list.
        plot : bool
            Display histogram if True.

        Returns
        -------
        tuple : (mean, median, sigma, max) seeing values
        """
        files = fits_file if fits_file is not None else self.fits_files
        seeing_values = []
        for f in files:
            try:
                with fits.open(f) as h:
                    s = h[0].header['SEEING']
                    if s is not None:
                        seeing_values.append(s)
            except Exception as e:
                print(f"Error reading {f}: {e}")

        mean_s = np.mean(seeing_values)
        median_s = np.median(seeing_values)
        sigma_s = np.std(seeing_values)
        max_s = np.max(seeing_values)

        if plot:
            plt.figure(figsize=(8, 5))
            plt.hist(seeing_values, bins=50, color='darkblue', edgecolor='orange')
            plt.axvline(mean_s, color='r', linewidth=2.5,
                        label=f'Mean = {mean_s:.3f}')
            plt.axvline(median_s, color='g', linewidth=2.5,
                        label=f'Median = {median_s:.3f}')
            plt.axvline(mean_s + sigma_s, color='r', linestyle='--', linewidth=1.5)
            plt.axvline(mean_s - sigma_s, color='r', linestyle='--', linewidth=1.5,
                        label=f'σ = ±{sigma_s:.3f}')
            plt.title(f'Seeing distribution — {self.file_name} ({len(seeing_values)} images)',
                      fontsize=14)
            plt.xlabel('Seeing (arcsec)', fontsize=13)
            plt.ylabel('N images', fontsize=13)
            plt.legend(fontsize=13)
            plt.tight_layout()
            plt.show()

        print(f"Seeing — mean: {mean_s:.3f}  median: {median_s:.3f}  "
              f"σ: {sigma_s:.3f}  max: {max_s:.3f}")
        return mean_s, median_s, sigma_s, max_s

    # ------------------------------------------------------------------
    # Image alignment
    # ------------------------------------------------------------------

    def align_images(self, output_dir=None, show=False, max_plots=None):
        """Align all images to the WCS of the first file via reprojection.

        Parameters
        ----------
        output_dir : str, optional
            Output directory for aligned FITS files.
        show : bool
            Display each aligned image.
        max_plots : int, optional
            Maximum number of images to display.

        Notes
        -----
        Processing time: ~12 s/image on a standard workstation.
        All aligned images are saved with the prefix ``aligned_``.
        """
        if output_dir is None:
            output_dir = os.path.join("ztf_data_aligned", self.file_name)
        os.makedirs(output_dir, exist_ok=True)

        with fits.open(self.first_file) as h:
            self.first_hdu = h[0]
            self.first_wcs = WCS(self.first_hdu.header)
            self.first_shape = self.first_hdu.data.shape
            self.first_wcs_header = self.first_wcs.to_header()

        for n, fits_file in enumerate(self.fits_files):
            fname = os.path.basename(fits_file)
            with fits.open(fits_file) as h:
                hdu = h[0]
                src_wcs = WCS(hdu.header)
                reprojected = self.reproject_astropy(
                    hdu.data, src_wcs, self.first_wcs, self.first_shape)
                new_hdr = hdu.header.copy()
                new_hdr.update(self.first_wcs_header)
                new_hdr['NAXIS'] = 2
                new_hdr['NAXIS1'] = reprojected.shape[1]
                new_hdr['NAXIS2'] = reprojected.shape[0]
                out_path = os.path.join(output_dir, f"aligned_{fname}")
                fits.PrimaryHDU(
                    data=reprojected.astype(np.float32), header=new_hdr
                ).writeto(out_path, overwrite=True)
            print(f"[{n+1}/{len(self.fits_files)}] Aligned → {out_path}")
            if show and (max_plots is None or n < max_plots):
                vmed = np.nanmedian(reprojected)
                self.show_img(reprojected, f"Aligned: {fname}",
                              vmin=vmed, vmax=vmed + 100)

        print("✓ Alignment complete.")

    # ------------------------------------------------------------------
    # Reference image construction
    # ------------------------------------------------------------------

    def get_refimg(self, aligned_files, ref_method='median', seeing_max=6,
                   zp_ref=None, limit_bkg=400, Imshow_ref=True, Imshow=False,
                   use_photutils_bkg=False, block_size=64, max_plots=None,
                   stack_number=100, output_file='ztf_data_ref', n=''):
        """Build a reference (sky background) image from aligned frames.

        The reference is constructed as a pixel-wise median (or mean) of
        images passing seeing and background quality cuts. The median is
        computed **line by line** from memory-mapped arrays to avoid
        loading all images simultaneously into RAM.

        Parameters
        ----------
        aligned_files : list of str
            Paths to aligned FITS files.
        ref_method : {'median', 'mean'}
            Stacking method.
        seeing_max : float
            Maximum seeing (arcsec) for a frame to be included.
        zp_ref : float, optional
            Zero-point for photometric normalization. Defaults to first image MAGZP.
        limit_bkg : float
            Maximum median sky background (ADU) — rejects cloudy frames.
        output_file : str
            Base output directory.
        n : str
            Suffix appended to output filename (useful for batch stacking).
        """
        first_hdu_a = fits.open(aligned_files[0])[0]
        num_images = 0
        valid_paths = []

        zp_ref = zp_ref if zp_ref is not None else first_hdu_a.header['MAGZP']
        print(f"Zero-point reference: {zp_ref}")

        if ref_method == 'mean':
            image_sum = np.zeros(first_hdu_a.data.shape, dtype=np.float32)
        elif ref_method == 'median':
            temp_dir = f".temp_refimg_{self.file_name}"
            os.makedirs(temp_dir, exist_ok=True)

        for fits_file in aligned_files:
            with fits.open(fits_file) as h:
                hdu = h[0]
                img = hdu.data.copy()
                zp_img = hdu.header['MAGZP']
                seeing = hdu.header['SEEING']
                saturate = hdu.header['SATURATE']
                sky_bkg = np.nanmedian(img)

                if sky_bkg > limit_bkg and not use_photutils_bkg:
                    print(f"  SKIP (cloudy, bkg={sky_bkg:.0f}): {os.path.basename(fits_file)}")
                    continue
                if seeing > seeing_max:
                    print(f"  SKIP (seeing={seeing:.2f}): {os.path.basename(fits_file)}")
                    continue

                num_images += 1
                img[img >= saturate] = np.nan

                if use_photutils_bkg and sky_bkg > limit_bkg:
                    img_sub, _ = self.remove_background_photutils(img, box_size=block_size)
                else:
                    img_sub, _ = self.remove_local_background(img, block_size=block_size)

                img_corr = img_sub * 10 ** ((zp_ref - zp_img) / 2.5)
                img_corr = self.degrade_psf(img_corr, seeing, seeing_max)
                print(f"  [{num_images}] seeing={seeing:.2f}  MAGZP={zp_img:.2f}  "
                      f"file={os.path.basename(fits_file)}")

                if ref_method == 'median':
                    tmp = os.path.join(temp_dir, f"img_{num_images:04d}.npy")
                    np.save(tmp, img_corr.astype(np.float32))
                    valid_paths.append(tmp)
                elif ref_method == 'mean':
                    image_sum += img_corr.astype(np.float32)

                if num_images >= stack_number:
                    break

        if num_images == 0:
            raise ValueError("No valid images found for reference construction.")

        out_path = os.path.join(output_file, self.file_name)
        os.makedirs(out_path, exist_ok=True)

        if ref_method == 'median':
            H, W = np.load(valid_paths[0]).shape
            image_median = np.empty((H, W), dtype=np.float32)
            print("Computing pixel-wise median (line by line)…")
            for i in range(H):
                rows = [np.load(p, mmap_mode='r')[i, :] for p in valid_paths]
                image_median[i, :] = np.median(np.vstack(rows), axis=0)
            out_file = os.path.join(out_path, f"refimg_median_{self.file_name}{n}.fits")
            hdr = first_hdu_a.header.copy()
            hdr['SEEING'] = (round(float(seeing_max), 3), 'PSF equalization target seeing')
            hdr['MAGZP'] = zp_ref
            hdr['NFRAMES'] = (num_images, 'Number of frames stacked')
            fits.PrimaryHDU(data=image_median, header=hdr).writeto(out_file, overwrite=True)
            print(f"✓ Median reference saved → {out_file}  ({num_images} frames)")
            if Imshow_ref:
                self.show_img(image_median, "Reference Image (Median)", vmin=0, vmax=30)
            shutil.rmtree(temp_dir)

        elif ref_method == 'mean':
            image_avg = image_sum / num_images
            out_file = os.path.join(out_path, f"refimg_mean_{self.file_name}{n}.fits")
            hdr = first_hdu_a.header.copy()
            hdr['SEEING'] = (round(float(seeing_max), 3), 'PSF equalization target seeing')
            hdr['NFRAMES'] = (num_images, 'Number of frames averaged')
            fits.PrimaryHDU(data=image_avg, header=hdr).writeto(out_file, overwrite=True)
            print(f"✓ Mean reference saved → {out_file}  ({num_images} frames)")
            if Imshow_ref:
                self.show_img(image_avg, "Reference Image (Mean)", vmin=0, vmax=100)

    # ------------------------------------------------------------------
    # Aperture photometry (standalone)
    # ------------------------------------------------------------------

    def star_photometry(self, x0, y0, radius, img, magzp):
        """Circular aperture photometry at pixel (x0, y0)."""
        x0 = round(float(x0))
        y0 = round(float(y0))
        flux = 0
        for xi in range(x0 - radius, x0 + radius):
            for yi in range(y0 - radius, y0 + radius):
                if np.sqrt((xi - x0) ** 2 + (yi - y0) ** 2) < radius:
                    flux += img[yi, xi]
        flux = max(flux, 1)
        mag = -2.5 * np.log10(flux) + magzp
        return Table(rows=[(flux, mag)], names=('flux', 'mag'))

    # ------------------------------------------------------------------
    # Transient detection
    # ------------------------------------------------------------------

    def detect_transients(self, diff_img, fwhm=3.0, threshold_sigma=5.0):
        """Detect point sources in a difference image using DAOStarFinder.

        Parameters
        ----------
        diff_img : 2D array
            Difference image (science − reference).
        fwhm : float
            Expected PSF FWHM in pixels.
        threshold_sigma : float
            Detection threshold in units of sigma.

        Returns
        -------
        astropy.table.Table or None
        """
        mean, median, std = sigma_clipped_stats(diff_img, sigma=3.0)
        daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold_sigma * std)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", AstropyWarning)
            sources = daofind(diff_img - median)
        return sources

    def group_transients_by_position(self, table, tol=3.0):
        """Spatially cluster transient detections across epochs."""
        positions = np.column_stack([table['x'], table['y']])
        used = np.zeros(len(positions), dtype=bool)
        clusters = []
        for i in range(len(positions)):
            if used[i]:
                continue
            group = [i]
            used[i] = True
            for j in range(i + 1, len(positions)):
                if not used[j] and np.linalg.norm(positions[i] - positions[j]) < tol:
                    group.append(j)
                    used[j] = True
            clusters.append(group)
        return clusters

    # ------------------------------------------------------------------
    # Image subtraction pipeline
    # ------------------------------------------------------------------

    def subtraction(self, aligned_files, ref_filename, ra, dec,
                    Start_Year, Start_Month, Start_Day,
                    End_Year, End_Month, End_Day,
                    zp=None, limit_bkg=500, seeing_reference=None,
                    Imshow=True, zoom=50, aperture=15,
                    use_photutils_bkg=False, block_size=64, max_plots=None):
        """Run image subtraction and extract light curve.

        For each science image in the date range:
        1. Applies photometric normalization and PSF equalization.
        2. Subtracts the reference image.
        3. Detects transient sources via DAOStarFinder.
        4. Runs aperture photometry on detected sources.

        Parameters
        ----------
        aligned_files : list of str
            Aligned FITS science images.
        ref_filename : str
            Path to the reference FITS image.
        ra, dec : float
            Sky coordinates (degrees) of the target.
        Start_Year, Start_Month, Start_Day : int
            Start of the analysis window.
        End_Year, End_Month, End_Day : int
            End of the analysis window.
        aperture : int
            Aperture radius in pixels for photometry.

        Returns
        -------
        astropy.table.Table or None
        """
        ref_hdu = fits.open(ref_filename)[0]
        image_ref = ref_hdu.data - np.nanmedian(ref_hdu.data)
        zp_ref = ref_hdu.header['MAGZP']
        seeing_ref = seeing_reference if seeing_reference else ref_hdu.header['SEEING']
        wcs_ref = WCS(ref_hdu.header)
        x_min_r, x_max_r, y_min_r, y_max_r = self.get_zoom(image_ref, ra, dec, wcs_ref, zoom)

        if zp is not None:
            image_ref = image_ref * 10 ** ((zp - zp_ref) / 2.5)
            zp_ref = zp

        start_mjd = Time(datetime(Start_Year, Start_Month, Start_Day)).mjd
        end_mjd = Time(datetime(End_Year, End_Month, End_Day)).mjd
        filtered_files = []
        for f in aligned_files:
            try:
                with fits.open(f) as h:
                    mjd = h[0].header['OBSMJD']
                    if start_mjd <= mjd <= end_mjd:
                        filtered_files.append(f)
            except Exception:
                continue

        tab_data = []
        num_images = 0

        for fits_file in filtered_files:
            with fits.open(fits_file) as h:
                hdu = h[0]
                img = hdu.data.copy()
                img_wcs = WCS(hdu.header)
                zp_img = hdu.header['MAGZP']
                seeing_img = hdu.header['SEEING']
                saturate = hdu.header['SATURATE']

                if np.nanmedian(img) > limit_bkg and not use_photutils_bkg:
                    continue
                if seeing_img > seeing_ref:
                    continue

                mjd = hdu.header['OBSMJD']
                date_obs = Time(mjd, format='mjd').to_datetime()
                date_str = date_obs.strftime('%Y-%m-%d')
                fracday = f"{date_obs.hour:02}{date_obs.minute:02}{date_obs.second:02}"

                img[img >= saturate] = np.nan
                if use_photutils_bkg and np.nanmedian(img) > limit_bkg:
                    img_sub, _ = self.remove_background_photutils(img, box_size=block_size)
                else:
                    img_sub, _ = self.remove_local_background(img, block_size=block_size)

                img_corr = img_sub * 10 ** ((zp_ref - zp_img) / 2.5)
                if seeing_img < seeing_ref:
                    img_corr = self.degrade_psf(img_corr, seeing_img, seeing_ref)

                diff = img_corr - image_ref
                diff -= np.nanmedian(diff)
                num_images += 1

                transients = self.detect_transients(diff, fwhm=seeing_ref)
                if transients is not None and len(transients) > 0:
                    print(f"  {len(transients)} transient(s) — {date_str}")
                    for idx, src in enumerate(transients):
                        x, y = src['xcentroid'], src['ycentroid']
                        phot = self.star_photometry(x, y, aperture, diff, zp_ref)
                        phot['date'] = date_str
                        phot['fracday'] = fracday
                        phot['x'] = x
                        phot['y'] = y
                        phot['id'] = idx
                        tab_data.append(phot)
                else:
                    print(f"  No transient — {date_str}")

                print(f"  Image {num_images}  seeing={seeing_img:.2f}  {date_str}")

        if not tab_data:
            print("No transients detected in the selected date range.")
            return None

        final_table = vstack(tab_data)
        clusters = self.group_transients_by_position(final_table)
        print(f"\n{len(final_table)} total detections → {len(clusters)} unique transients.")
        return final_table
