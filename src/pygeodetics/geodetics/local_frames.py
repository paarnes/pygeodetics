"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com


Local-tangent-plane (ENU / NED) and AER coordinate conversions.

Implements bidirectional transformations between geodetic coordinates,
ECEF, the local East-North-Up (ENU) frame, the local North-East-Down
(NED) frame, and observer-target Azimuth-Elevation-Range (AER) tuples.

All functions take ``radians=False`` by default (latitude / longitude /
azimuth / elevation in degrees) for consistency with the rest of the
package, and accept either scalar or NumPy-array inputs.

The core ECEF<->ENU rotation reuses the existing
:func:`pygeodetics.geodetics.ECEF2enu.ECEF2enu` helper. Geodetic<->ECEF
reuses :func:`pygeodetics.geodetics.geod2ECEF.geod2ECEF` and
:func:`pygeodetics.geodetics.ECEF2geod.ECEF2geodb`.

Azimuth convention
------------------
Azimuth is measured clockwise from the local north axis, in the range
``[0°, 360°)``. Elevation is measured from the local horizontal plane,
positive upward, in the range ``[-90°, 90°]``. Slant range is the
straight-line 3D distance in metres.
"""

from __future__ import annotations

from typing import Tuple, Union

import numpy as np

from ..Ellipsoid import Ellipsoid, WGS84
from .ECEF2enu import ECEF2enu
from .ECEF2geod import ECEF2geodb
from .geod2ECEF import geod2ECEF


ArrayLike = Union[float, int, np.ndarray, list, tuple]
_TripleArr = Tuple[Union[float, np.ndarray],
                   Union[float, np.ndarray],
                   Union[float, np.ndarray]]


def _resolve_ellipsoid(ellipsoid):
    if ellipsoid is None:
        return WGS84()
    if not isinstance(ellipsoid, Ellipsoid):
        raise TypeError(
            "ellipsoid must be a pygeodetics.Ellipsoid instance "
            f"(got {type(ellipsoid).__name__})"
        )
    return ellipsoid


def _scalarize(*arrays, scalar):
    if scalar:
        return tuple(float(np.asarray(a).reshape(())) for a in arrays)
    return arrays


def _geod_to_ecef_arr(lat, lon, h, ellipsoid, radians):
    if not radians:
        lat = np.radians(lat)
        lon = np.radians(lon)
    return geod2ECEF(lat, lon, h, ellipsoid=ellipsoid, radians=True)


def _ecef_to_geod_arr(X, Y, Z, ellipsoid, radians):
    """Vectorised wrapper around the Bowring (ECEF2geodb) implementation."""
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    Z = np.asarray(Z, dtype=float)
    angle_unit = "rad" if radians else "deg"
    return ECEF2geodb(ellipsoid.a, ellipsoid.b, X, Y, Z, angle_unit=angle_unit)


def geodetic2enu(lat: ArrayLike, lon: ArrayLike, h: ArrayLike,
                 lat0: ArrayLike, lon0: ArrayLike, h0: ArrayLike,
                 ellipsoid: Ellipsoid = None,
                 radians: bool = False) -> _TripleArr:
    """
    Convert geodetic target coordinates ``(lat, lon, h)`` into the local
    ENU frame anchored at the reference point ``(lat0, lon0, h0)``.

    Parameters
    ----------
    lat, lon, h : float or array_like
        Target geodetic coordinates (height in metres).
    lat0, lon0, h0 : float or array_like
        Reference (origin) geodetic coordinates.
    ellipsoid : Ellipsoid, optional
        Reference ellipsoid (defaults to WGS84).
    radians : bool, default False
        If ``True``, latitude / longitude inputs are in radians.

    Returns
    -------
    e, n, u : float or np.ndarray
        East, North, Up offsets in metres.
    """
    ellipsoid = _resolve_ellipsoid(ellipsoid)
    X, Y, Z = _geod_to_ecef_arr(lat, lon, h, ellipsoid, radians)
    return ECEF2enu(X, Y, Z, lat0, lon0, h0,
                    ellipsoid=ellipsoid, radians=radians)


def enu2ecef(e: ArrayLike, n: ArrayLike, u: ArrayLike,
             lat0: ArrayLike, lon0: ArrayLike, h0: ArrayLike,
             ellipsoid: Ellipsoid = None,
             radians: bool = False) -> _TripleArr:
    """
    Convert local ENU offsets back to absolute ECEF coordinates given the
    reference point ``(lat0, lon0, h0)``.
    """
    ellipsoid = _resolve_ellipsoid(ellipsoid)
    e_arr = np.asarray(e, dtype=float)
    n_arr = np.asarray(n, dtype=float)
    u_arr = np.asarray(u, dtype=float)
    scalar = (e_arr.ndim == 0 and n_arr.ndim == 0 and u_arr.ndim == 0
              and np.ndim(lat0) == 0 and np.ndim(lon0) == 0
              and np.ndim(h0) == 0)

    lat0_r = np.asarray(lat0, dtype=float)
    lon0_r = np.asarray(lon0, dtype=float)
    if not radians:
        lat0_r = np.radians(lat0_r)
        lon0_r = np.radians(lon0_r)

    sin_lat = np.sin(lat0_r)
    cos_lat = np.cos(lat0_r)
    sin_lon = np.sin(lon0_r)
    cos_lon = np.cos(lon0_r)

    # Inverse of the ENU rotation (its transpose).
    dx = -sin_lon * e_arr - sin_lat * cos_lon * n_arr + cos_lat * cos_lon * u_arr
    dy =  cos_lon * e_arr - sin_lat * sin_lon * n_arr + cos_lat * sin_lon * u_arr
    dz =                     cos_lat            * n_arr + sin_lat            * u_arr

    x0, y0, z0 = geod2ECEF(lat0_r, lon0_r, h0,
                           ellipsoid=ellipsoid, radians=True)
    X = dx + x0
    Y = dy + y0
    Z = dz + z0
    return _scalarize(X, Y, Z, scalar=scalar)


def enu2geodetic(e: ArrayLike, n: ArrayLike, u: ArrayLike,
                 lat0: ArrayLike, lon0: ArrayLike, h0: ArrayLike,
                 ellipsoid: Ellipsoid = None,
                 radians: bool = False) -> _TripleArr:
    """
    Convert local ENU offsets to absolute geodetic coordinates.
    """
    ellipsoid = _resolve_ellipsoid(ellipsoid)
    X, Y, Z = enu2ecef(e, n, u, lat0, lon0, h0,
                       ellipsoid=ellipsoid, radians=radians)
    return _ecef_to_geod_arr(X, Y, Z, ellipsoid, radians)


def geodetic2ned(lat, lon, h, lat0, lon0, h0,
                 ellipsoid: Ellipsoid = None,
                 radians: bool = False) -> _TripleArr:
    """Geodetic -> local NED (North, East, Down) frame."""
    e, n, u = geodetic2enu(lat, lon, h, lat0, lon0, h0,
                           ellipsoid=ellipsoid, radians=radians)
    return n, e, _neg(u)


def ned2ecef(n, e, d, lat0, lon0, h0,
             ellipsoid: Ellipsoid = None,
             radians: bool = False) -> _TripleArr:
    """Local NED -> absolute ECEF."""
    return enu2ecef(e, n, _neg(d), lat0, lon0, h0,
                    ellipsoid=ellipsoid, radians=radians)


def ned2geodetic(n, e, d, lat0, lon0, h0,
                 ellipsoid: Ellipsoid = None,
                 radians: bool = False) -> _TripleArr:
    """Local NED -> absolute geodetic."""
    return enu2geodetic(e, n, _neg(d), lat0, lon0, h0,
                        ellipsoid=ellipsoid, radians=radians)


def _neg(x):
    if isinstance(x, np.ndarray):
        return -x
    return -np.asarray(x) if not np.isscalar(x) else -x


#
# Azimuth: clockwise from north, [0, 360) deg (or [0, 2*pi) rad).
# Elevation: from local horizontal plane, [-90, 90] deg.
# Range: 3-D Euclidean slant distance in metres.

def enu2aer(e: ArrayLike, n: ArrayLike, u: ArrayLike,
            radians: bool = False) -> _TripleArr:
    """Convert ENU offsets to (azimuth, elevation, slant range)."""
    e_arr = np.asarray(e, dtype=float)
    n_arr = np.asarray(n, dtype=float)
    u_arr = np.asarray(u, dtype=float)
    scalar = e_arr.ndim == 0 and n_arr.ndim == 0 and u_arr.ndim == 0

    horiz = np.hypot(e_arr, n_arr)
    srange = np.hypot(horiz, u_arr)
    az = np.arctan2(e_arr, n_arr)
    az = np.mod(az, 2.0 * np.pi)
    # Avoid 0/0 -> NaN at the origin: elevation defaults to 0.
    with np.errstate(divide="ignore", invalid="ignore"):
        el = np.where(srange > 0,
                      np.arctan2(u_arr, horiz),
                      np.zeros_like(srange))

    if not radians:
        az = np.degrees(az)
        el = np.degrees(el)

    return _scalarize(az, el, srange, scalar=scalar)


def aer2enu(az: ArrayLike, el: ArrayLike, srange: ArrayLike,
            radians: bool = False) -> _TripleArr:
    """Convert (azimuth, elevation, slant range) to ENU offsets."""
    az_arr = np.asarray(az, dtype=float)
    el_arr = np.asarray(el, dtype=float)
    r_arr = np.asarray(srange, dtype=float)
    scalar = az_arr.ndim == 0 and el_arr.ndim == 0 and r_arr.ndim == 0

    if np.any(r_arr < 0):
        raise ValueError("Slant range must be non-negative.")

    if not radians:
        az_arr = np.radians(az_arr)
        el_arr = np.radians(el_arr)

    cos_el = np.cos(el_arr)
    e = r_arr * cos_el * np.sin(az_arr)
    n = r_arr * cos_el * np.cos(az_arr)
    u = r_arr * np.sin(el_arr)
    return _scalarize(e, n, u, scalar=scalar)


def ned2aer(n, e, d, radians: bool = False) -> _TripleArr:
    """Convert NED offsets to (azimuth, elevation, slant range)."""
    return enu2aer(e, n, _neg(d), radians=radians)


def aer2ned(az, el, srange, radians: bool = False) -> _TripleArr:
    """Convert (azimuth, elevation, slant range) to NED offsets."""
    e, n, u = aer2enu(az, el, srange, radians=radians)
    return n, e, _neg(u)


def ecef2aer(X, Y, Z, lat0, lon0, h0,
             ellipsoid: Ellipsoid = None,
             radians: bool = False) -> _TripleArr:
    """Compute (azimuth, elevation, slant range) of an ECEF point relative
    to a geodetic observer ``(lat0, lon0, h0)``."""
    ellipsoid = _resolve_ellipsoid(ellipsoid)
    e, n, u = ECEF2enu(X, Y, Z, lat0, lon0, h0,
                       ellipsoid=ellipsoid, radians=radians)
    return enu2aer(e, n, u, radians=radians)


def aer2ecef(az, el, srange, lat0, lon0, h0,
             ellipsoid: Ellipsoid = None,
             radians: bool = False) -> _TripleArr:
    """Convert AER (relative to a geodetic observer) to absolute ECEF."""
    ellipsoid = _resolve_ellipsoid(ellipsoid)
    e, n, u = aer2enu(az, el, srange, radians=radians)
    return enu2ecef(e, n, u, lat0, lon0, h0,
                    ellipsoid=ellipsoid, radians=radians)


def geodetic2aer(lat, lon, h, lat0, lon0, h0,
                 ellipsoid: Ellipsoid = None,
                 radians: bool = False) -> _TripleArr:
    """Compute (azimuth, elevation, slant range) of a geodetic target
    relative to a geodetic observer."""
    ellipsoid = _resolve_ellipsoid(ellipsoid)
    e, n, u = geodetic2enu(lat, lon, h, lat0, lon0, h0,
                           ellipsoid=ellipsoid, radians=radians)
    return enu2aer(e, n, u, radians=radians)


def aer2geodetic(az, el, srange, lat0, lon0, h0,
                 ellipsoid: Ellipsoid = None,
                 radians: bool = False) -> _TripleArr:
    """Convert AER (relative to a geodetic observer) to absolute geodetic
    coordinates."""
    ellipsoid = _resolve_ellipsoid(ellipsoid)
    X, Y, Z = aer2ecef(az, el, srange, lat0, lon0, h0,
                       ellipsoid=ellipsoid, radians=radians)
    return _ecef_to_geod_arr(X, Y, Z, ellipsoid, radians)
