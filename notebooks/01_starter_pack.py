"""
notebooks/01_starter_pack.py
=============================
Full pipeline walkthrough — equivalent to Astro_StarterPack.ipynb

Run as a Jupyter notebook or convert with:
    jupyter nbconvert --to notebook --execute 01_starter_pack.py

Author: Florian Devender-Dauge
"""

# ─────────────────────────────────────────────
# 0. Imports
# ─────────────────────────────────────────────

import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astrotools import AstroTools, ZTFDownloader

# ─────────────────────────────────────────────
# 1. Settings — choose your target
# ─────────────────────────────────────────────

Field    = "0664"      # ZTF field (0001–1800)
Filter   = "zr"        # zr | zg | zi
Ccd      = "06"        # 01–16
Quadrant = "1"         # 1–4
START    = "20200101"  # earliest: 20180321
END      = "20200501"  # latest:   20210108

FFCQ = f"00{Field}_{Filter}_c{Ccd}_o_q{Quadrant}"

# ─────────────────────────────────────────────
# 2. Download images (requires fichiers_par_field.csv)
# ─────────────────────────────────────────────

downloader = ZTFDownloader(data_type="sciimg", output_base_dir="ztf_data")
downloader.download_from_csv(
    csv_path="fichiers_par_field.csv",
    target_field=Field,
    filter_band=Filter,
    ccd=Ccd,
    quadrant=Quadrant,
    start_date=START,
    end_date=END,
)

# For parallel multi-quadrant download (16 threads):
# downloader.download_parallel(
#     csv_path="fichiers_par_field.csv",
#     target_field=Field, filter_band=Filter,
#     ccd_range=range(1, 5), quad_range=range(1, 5),
#     start_date=START, end_date=END,
#     max_workers=16
# )

# ─────────────────────────────────────────────
# 3. Load files & initialize AstroTools
# ─────────────────────────────────────────────

fits_files = sorted(glob.glob(f"ztf_data/{FFCQ}/*.fits"))
print(f"Found {len(fits_files)} FITS files")

tool = AstroTools(fits_files, file_name=FFCQ)

# ─────────────────────────────────────────────
# 4. Seeing statistics
# ─────────────────────────────────────────────

mean_s, median_s, sigma_s, max_s = tool.distrib_seeing(plot=True)
print(f"Seeing: mean={mean_s:.2f}  median={median_s:.2f}  max={max_s:.2f}")

# ─────────────────────────────────────────────
# 5. Align all images to common WCS
# ─────────────────────────────────────────────
# ~12 seconds per image on a standard workstation

tool.align_images()

aligned_files = sorted(glob.glob(f"ztf_data_aligned/{FFCQ}/*.fits"))
print(f"Aligned: {len(aligned_files)} images")

# ─────────────────────────────────────────────
# 6. Build reference image (pixel-wise median)
# ─────────────────────────────────────────────
# Memory-efficient: computed line by line from memory-mapped arrays

tool.get_refimg(
    aligned_files,
    ref_method="median",
    seeing_max=6,
    zp_ref=27.5,
    limit_bkg=300,
    Imshow_ref=True,
)

ref_filename = f"ztf_data_ref/{FFCQ}/refimg_median_{FFCQ}.fits"

# ─────────────────────────────────────────────
# 7. Temporal stacking (groups of 3)
# ─────────────────────────────────────────────
# SNR gain: √3 ≈ +0.6 mag in limiting magnitude

valid_files = []
for f in aligned_files:
    with fits.open(f) as h:
        seeing = h[0].header["SEEING"]
        bkg = np.nanmedian(h[0].data)
        if seeing <= 6 and bkg < 400:
            valid_files.append(f)

print(f"{len(valid_files)} images pass quality cuts")

for i in range(0, len(valid_files), 3):
    batch = valid_files[i:i+3]
    if len(batch) >= 2:
        tool.get_refimg(
            batch,
            ref_method="median",
            seeing_max=6,
            zp_ref=27.5,
            limit_bkg=400,
            Imshow_ref=False,
            output_file="ztf_data_stack",
            n=str(i),
        )

# ─────────────────────────────────────────────
# 8. Image subtraction + transient detection
# ─────────────────────────────────────────────

stack_files = sorted(glob.glob(f"ztf_data_stack/{FFCQ}/*.fits"))

# Target coordinates (example: field centre)
with fits.open(aligned_files[0]) as h:
    wcs = WCS(h[0].header)
ra, dec = tool.pixel_to_radec(1400, 600, wcs)
print(f"Target: RA={ra:.4f}  Dec={dec:.4f}")

results = tool.subtraction(
    stack_files,
    ref_filename,
    ra, dec,
    Start_Year=2020, Start_Month=1, Start_Day=1,
    End_Year=2020, End_Month=5, End_Day=1,
    zp=None,
    limit_bkg=400,
    Imshow=True,
    zoom=50,
    aperture=15,
)

print(results)

# ─────────────────────────────────────────────
# 9. Plot light curve
# ─────────────────────────────────────────────

if results is not None:
    plt.figure(figsize=(10, 5))
    plt.scatter(results["date"], results["flux"], s=30, color="steelblue", label="Flux (ADU)")
    plt.xlabel("Date")
    plt.ylabel("Flux (ADU)")
    plt.title(f"Light curve — {FFCQ}")
    plt.xticks(rotation=30)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"lightcurve_{FFCQ}.png", dpi=150)
    plt.show()
