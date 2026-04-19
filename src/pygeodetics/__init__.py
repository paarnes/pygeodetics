"""
PyGeodetics - A Python library for geodetic calculations.

Author: Per Helge Aarnes
Email: per.helge.aarnes@gmail.com
"""

# Define package version
__version__ = "0.1.0"

# Import necessary modules
from .Ellipsoid import (
    Ellipsoid,
    WGS84,
    GRS80,
    International1924,
    Clarke1866,
    BesselModified,
    Bessel1841,
)
from .geodetics.geod2ECEF import geod2ECEF
from .geodetics.ECEF2geod import ECEF2geodb, ECEF2geodv, ECEF2geod
from .geodetics.ECEF2enu import ECEF2enu
from .geodetics.geodetic_inverse_problem import geodetic_inverse_problem
from .geodetics.geodetic_direct_problem import geodetic_direct_problem
from .geodetics.radius_of_curvature_azimuth import radius_of_curvature_azimuth
from .geodetics.Vincenty import vincenty_distance
from .geodetics.Mrad import Mrad
from .geodetics.Nrad import Nrad
from .geodetics.footpoint_latitude import footpoint_latitude
from .geodetics.meridional_arc_dist import meridional_arc_dist
from .projections.MercatorVariantC import MercatorVariantC
from .projections.TransverseMercator import TransverseMercator
from .projections.GridConvergence import (
    tm_grid_convergence_geographic,
    tm_grid_convergence_projected,
    tm_sphere_grid_conv_projected,
)
from .projections.ScaleFactor import (
    tm_point_scale_factor_geographic,
    tm_point_scale_factor_projected,
    tm_sphere_point_scale_factor,
)

# Import main Geodetic class
from .Geodetic import Geodetic


# Expose public API
__all__ = [
    "Ellipsoid",
    "WGS84",
    "GRS80",
    "International1924",
    "Clarke1866",
    "BesselModified",
    "Bessel1841",
    "geod2ECEF",
    "ECEF2geod",
    "ECEF2geodb",
    "ECEF2geodv",
    "ECEF2enu",
    "geodetic_inverse_problem",
    "geodetic_direct_problem",
    "radius_of_curvature_azimuth",
    "vincenty_distance",
    "Mrad",
    "Nrad",
    "footpoint_latitude",
    "meridional_arc_dist",
    "MercatorVariantC",
    "TransverseMercator",
    "tm_grid_convergence_geographic",
    "tm_grid_convergence_projected",
    "tm_sphere_grid_conv_projected",
    "tm_point_scale_factor_geographic",
    "tm_point_scale_factor_projected",
    "tm_sphere_point_scale_factor",
    "Geodetic",
]
