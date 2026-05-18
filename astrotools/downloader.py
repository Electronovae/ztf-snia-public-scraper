"""
downloader.py — ZTF public archive access utilities

ZTFDownloader: Download science images from the IRSA public archive
               using a pre-built CSV index.

Author: Florian Devender-Dauge
"""

import csv
import os
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm


class ZTFDownloader:
    """Download ZTF science images from the IRSA public archive.

    Images are identified via a CSV index built from IRSA HTML pages
    (see ``notebooks/00_build_index.ipynb``).

    Parameters
    ----------
    data_type : str
        FITS file type suffix. Use ``'sciimg'`` for science images
        or ``'mskimg'`` for mask images.
    output_base_dir : str
        Root directory where downloaded FITS files are stored.

    Examples
    --------
    Single quadrant download::

        downloader = ZTFDownloader(data_type="sciimg")
        downloader.download_from_csv(
            csv_path="fichiers_par_field.csv",
            target_field="0664",
            filter_band="zr",
            ccd="06",
            quadrant="1",
            start_date="20200101",
            end_date="20200501"
        )

    Parallel multi-quadrant download::

        downloader.download_parallel(
            csv_path="fichiers_par_field.csv",
            target_field="0664",
            filter_band="zr",
            ccd_range=range(1, 5),
            quad_range=range(1, 5),
            start_date="20200101",
            end_date="20200501",
            max_workers=16
        )
    """

    BASE_URL = 'https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/'

    def __init__(self, data_type="sciimg", output_base_dir='ztf_data'):
        self.data_type = data_type
        self.output_base_dir = output_base_dir
        os.makedirs(self.output_base_dir, exist_ok=True)

    def _build_filename(self, year, month, day, fractime, field, filter_band, ccd, quadrant):
        """Build ZTF FITS filename from components."""
        datefull = f"{year}{month}{day}{fractime}"
        prefix = f"ztf_{datefull}_00{field}_{filter_band}_c{ccd}_o_q"
        return f"{prefix}{quadrant}_{self.data_type}.fits"

    def get_one_file(self, year, month, day, fractime, field, filter_band, ccd, quadrant):
        """Download a single FITS file from IRSA.

        Parameters
        ----------
        year, month, day : str
            Observation date components (e.g., '2020', '01', '15').
        fractime : str
            Fractional time of observation (6-digit string).
        field : str
            ZTF field number (4-digit string, e.g., '0664').
        filter_band : str
            Photometric band ('zr', 'zg', or 'zi').
        ccd : str
            CCD number ('01' to '16').
        quadrant : str
            Readout quadrant ('1' to '4').
        """
        filename = self._build_filename(
            year, month, day, fractime, field, filter_band, ccd, quadrant)
        output_dir = os.path.join(
            self.output_base_dir, f"00{field}_{filter_band}_c{ccd}_o_q{quadrant}")
        os.makedirs(output_dir, exist_ok=True)
        filepath = os.path.join(output_dir, filename)

        if os.path.exists(filepath):
            return  # skip already downloaded files

        url = f"{self.BASE_URL}/{year}/{month}{day}/{fractime}/{filename}"
        response = requests.get(url, stream=True)
        if response.status_code == 200:
            total = int(response.headers.get('content-length', 0))
            with open(filepath, 'wb') as f:
                for chunk in response.iter_content(1024):
                    f.write(chunk)
        elif response.status_code == 404:
            pass  # file not available for this observation
        else:
            print(f"HTTP {response.status_code}: {url}")

    def download_from_csv(self, csv_path, target_field, filter_band, ccd,
                          quadrant, start_date, end_date):
        """Download images for one quadrant from the CSV index.

        Parameters
        ----------
        csv_path : str
            Path to the pre-built ``fichiers_par_field.csv`` index.
        target_field : str
            ZTF field number (4-digit string).
        filter_band : str
            Photometric band.
        ccd : str
            CCD number.
        quadrant : str
            Readout quadrant.
        start_date, end_date : str
            Date range in ``'YYYYMMDD'`` format.
        """
        with open(csv_path, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            if target_field not in reader.fieldnames:
                print(f"Field '{target_field}' not found in CSV.")
                return
            for row in reader:
                entry = row[target_field]
                if entry and start_date <= entry[:8] <= end_date:
                    year, month, day = entry[:4], entry[4:6], entry[6:8]
                    fractime = entry[8:]
                    self.get_one_file(
                        year, month, day, fractime,
                        target_field, filter_band, ccd, quadrant)

    def download_parallel(self, csv_path, target_field, filter_band,
                          ccd_range, quad_range, start_date, end_date,
                          max_workers=16):
        """Download multiple quadrants in parallel using ThreadPoolExecutor.

        Parameters
        ----------
        ccd_range : iterable of int
            CCD numbers to process (e.g., ``range(1, 5)``).
        quad_range : iterable of int
            Quadrant numbers to process (e.g., ``range(1, 5)``).
        max_workers : int
            Number of parallel download threads.

        Examples
        --------
        >>> downloader.download_parallel(
        ...     csv_path="fichiers_par_field.csv",
        ...     target_field="0664", filter_band="zr",
        ...     ccd_range=range(1, 5), quad_range=range(1, 5),
        ...     start_date="20200101", end_date="20200501",
        ...     max_workers=16
        ... )
        """
        args_list = [
            (f"{ccd:02d}", str(quad))
            for ccd in ccd_range
            for quad in quad_range
        ]

        def _task(args):
            ccd_str, quad_str = args
            try:
                self.download_from_csv(
                    csv_path, target_field, filter_band,
                    ccd_str, quad_str, start_date, end_date)
                return f"✅ CCD {ccd_str} Q{quad_str}"
            except Exception as e:
                return f"❌ CCD {ccd_str} Q{quad_str}: {e}"

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(_task, a): a for a in args_list}
            for future in tqdm(as_completed(futures), total=len(futures),
                               desc="Downloading quadrants"):
                print(future.result())
