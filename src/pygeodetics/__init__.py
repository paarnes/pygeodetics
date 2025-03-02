"""
PyGeodetics - A Python library for geodetic calculations.

Author: Per Helge Aarnes
Email: per.helge.aarnes@gmail.com
"""

# Define package version
__version__ = "0.1.0"

# Import necessary modules
from .Ellipsoid import Ellipsoid, WGS84
from .geodetics.geod2ECEF import geod2ECEF
from .geodetics.ECEF2geod import ECEF2geodb
from .geodetics.ECEF2enu import ECEF2enu
from .geodetics.geodetic_inverse_problem import geodetic_inverse_problem
from .geodetics.geodetic_direct_problem import geodetic_direct_problem
from .geodetics.radius_of_curvature_azimuth import radius_of_curvature_azimuth

# Import main Geodetic class
from Geodetic import Geodetic


# Expose public API
__all__ = [
    "Ellipsoid",
    "WGS84",
    "geod2ECEF",
    "ECEF2geodb",
    "ECEF2enu",
    "geodetic_inverse_problem",
    "geodetic_direct_problem",
    "radius_of_curvature_azimuth",
    "Geodetic",
]
