"""
PyGeodetics - A Python library for geodetic calculations.

Author: Per Helge Aarnes
Email: per.helge.aarnes@gmail.com
"""

from importlib.metadata import PackageNotFoundError, version as _pkg_version

# Derive package version from installed distribution metadata so it stays
# in sync with pyproject.toml without manual edits on each release.
try:
    __version__ = _pkg_version("pygeodetics")
except PackageNotFoundError:  # package is not installed (e.g. running from source tree)
    __version__ = "0.0.0+unknown"

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
from .geodetics.local_frames import (
    geodetic2enu,
    enu2geodetic,
    enu2ecef,
    geodetic2ned,
    ned2geodetic,
    ned2ecef,
    enu2aer,
    aer2enu,
    ned2aer,
    aer2ned,
    ecef2aer,
    aer2ecef,
    geodetic2aer,
    aer2geodetic,
)
from .geodetics.geodetic_inverse_problem import geodetic_inverse_problem
from .geodetics.geodetic_direct_problem import geodetic_direct_problem
from .geodetics.geodesic_polygon import (
    polygon_perimeter,
    polygon_area,
    polygon_centroid,
    polygon_bounds,
    polygon_densify,
    geodesic_interpolate,
)
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
from .projections.UTM import (
    geodetic_to_utm,
    utm_to_geodetic,
    mgrs_band_letter,
)
from .projections.PolarStereographic import PolarStereographic
from .projections.UPS import (
    geodetic_to_ups,
    ups_to_geodetic,
)
from .projections.MGRS import (
    to_mgrs,
    from_mgrs,
    utm_epsg,
    ups_epsg,
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
    "geodetic2enu",
    "enu2geodetic",
    "enu2ecef",
    "geodetic2ned",
    "ned2geodetic",
    "ned2ecef",
    "enu2aer",
    "aer2enu",
    "ned2aer",
    "aer2ned",
    "ecef2aer",
    "aer2ecef",
    "geodetic2aer",
    "aer2geodetic",
    "geodetic_inverse_problem",
    "geodetic_direct_problem",
    "polygon_perimeter",
    "polygon_area",
    "polygon_centroid",
    "polygon_bounds",
    "polygon_densify",
    "geodesic_interpolate",
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
    "geodetic_to_utm",
    "utm_to_geodetic",
    "mgrs_band_letter",
    "PolarStereographic",
    "geodetic_to_ups",
    "ups_to_geodetic",
    "to_mgrs",
    "from_mgrs",
    "utm_epsg",
    "ups_epsg",
    "Geodetic",
]
