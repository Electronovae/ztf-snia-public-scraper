"""
ZTF-Tools: Open-source Python framework for ZTF public data processing
and Type Ia supernova detection.

Author: Florian Devender-Dauge
Institution: EUPI — Université Clermont Auvergne
Supervisors: Philippe Rosnet, Marie Aubert (LPC Clermont)
"""

from .pipeline import AstroTools, aphoto_modified
from .downloader import ZTFDownloader

__version__ = "0.2.0"
__author__ = "Florian Devender-Dauge"
__email__ = "fdevenderdauge@protonmail.com"
__license__ = "MIT"

__all__ = [
    "AstroTools",
    "aphoto_modified",
    "ZTFDownloader",
]
