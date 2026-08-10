"""
gut_ibm_tools – Analysis toolkit for GutIBM simulation output.
Compatible with nufeb_tools HDF5 format.
"""

from . import (
    analysis,
    colony,
    fish_observation,
    spatial_stats,
    validation,
    visualization,
)
from .hdf5_reader import GutIBMData

__version__ = "0.1.0"

__all__ = [
    "GutIBMData",
    "__version__",
    "analysis",
    "colony",
    "fish_observation",
    "spatial_stats",
    "validation",
    "visualization",
]
