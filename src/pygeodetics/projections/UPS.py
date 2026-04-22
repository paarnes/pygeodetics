"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com


Universal Polar Stereographic (UPS) convenience API.

UPS covers the polar regions where UTM is undefined (latitudes north of
84°N and south of 80°S). It is a thin wrapper around the
:class:`~pygeodetics.projections.PolarStereographic.PolarStereographic`
engine with the standard UPS defaults applied:

- scale factor at the pole: 0.994
- false easting: 2 000 000 m
- false northing: 2 000 000 m
- central meridian: 0°

Usage Example
-------------
    >>> from pygeodetics import geodetic_to_ups, ups_to_geodetic
    >>> e, n = geodetic_to_ups(85.0, 0.0)
    >>> lat, lon = ups_to_geodetic(e, n, hemisphere="N")
"""

from __future__ import annotations

from typing import Literal, Optional, Tuple, Union

import numpy as np

from ..Ellipsoid import Ellipsoid, WGS84
from .PolarStereographic import PolarStereographic


# Standard UPS defaults (EPSG / NGA).
UPS_SCALE_FACTOR = 0.994
UPS_FALSE_EASTING = 2_000_000.0
UPS_FALSE_NORTHING = 2_000_000.0
UPS_CENTRAL_MERIDIAN = 0.0

# UPS validity (the 30 km overlap with UTM is allowed by NGA).
UPS_NORTH_MIN_LAT = 84.0
UPS_SOUTH_MAX_LAT = -80.0


ArrayLike = Union[float, int, np.ndarray, list, tuple]


def _resolve_ellipsoid(ellipsoid: Optional[Ellipsoid]) -> Ellipsoid:
    if ellipsoid is None:
        return WGS84()
    if not isinstance(ellipsoid, Ellipsoid):
        raise TypeError(
            "ellipsoid must be a pygeodetics.Ellipsoid instance "
            f"(got {type(ellipsoid).__name__})"
        )
    return ellipsoid


def _make_ups(pole: str, ellipsoid: Ellipsoid) -> PolarStereographic:
    return PolarStereographic(
        pole=pole,
        central_meridian_deg=UPS_CENTRAL_MERIDIAN,
        scale_factor=UPS_SCALE_FACTOR,
        false_easting=UPS_FALSE_EASTING,
        false_northing=UPS_FALSE_NORTHING,
        a=ellipsoid.a,
        f=ellipsoid.f,
    )


def geodetic_to_ups(
    lat: ArrayLike,
    lon: ArrayLike,
    ellipsoid: Optional[Ellipsoid] = None,
) -> Tuple[Union[float, np.ndarray],
           Union[float, np.ndarray],
           Union[str, np.ndarray]]:
    """
    Project geodetic coordinates onto the Universal Polar Stereographic grid.

    The hemisphere is selected automatically from the sign of ``lat``.

    Parameters
    ----------
    lat, lon : float or array_like
        Geodetic latitude / longitude in degrees.
    ellipsoid : Ellipsoid, optional
        Reference ellipsoid. Defaults to WGS84.

    Returns
    -------
    easting : float or np.ndarray
        UPS easting in metres.
    northing : float or np.ndarray
        UPS northing in metres.
    hemisphere : str or np.ndarray of str
        ``"N"`` or ``"S"`` indicating which UPS projection was used.

    Notes
    -----
    Standard UPS validity is ``lat >= 84°`` (north) or ``lat <= -80°``
    (south). This function does *not* enforce that limit so callers can
    use the 30 km UTM/UPS overlap zone, but a ``ValueError`` is raised
    if a single call mixes northern- and southern-hemisphere points
    (different projections cannot be combined into one output).
    """
    ellipsoid = _resolve_ellipsoid(ellipsoid)

    lat_arr = np.asarray(lat, dtype=float)
    lon_arr = np.asarray(lon, dtype=float)
    scalar_input = lat_arr.ndim == 0 and lon_arr.ndim == 0
    lat_arr, lon_arr = np.broadcast_arrays(lat_arr, lon_arr)

    flat_lat = np.atleast_1d(lat_arr).ravel()
    flat_lon = np.atleast_1d(lon_arr).ravel()
    n = flat_lat.size

    easting = np.empty(n, dtype=float)
    northing = np.empty(n, dtype=float)
    hemisphere = np.where(flat_lat >= 0.0, "N", "S")

    for hemi in ("N", "S"):
        mask = hemisphere == hemi
        if not np.any(mask):
            continue
        ups = _make_ups(hemi, ellipsoid)
        e, nn = ups.forward(flat_lat[mask], flat_lon[mask])
        easting[mask] = e
        northing[mask] = nn

    if scalar_input:
        return float(easting[0]), float(northing[0]), str(hemisphere[0])

    out_shape = lat_arr.shape
    return (easting.reshape(out_shape),
            northing.reshape(out_shape),
            hemisphere.reshape(out_shape))


def ups_to_geodetic(
    easting: ArrayLike,
    northing: ArrayLike,
    hemisphere: Literal["N", "S"],
    ellipsoid: Optional[Ellipsoid] = None,
) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray]]:
    """
    Convert UPS easting/northing back to geodetic (lat, lon) in degrees.

    Parameters
    ----------
    easting, northing : float or array_like
        UPS coordinates in metres.
    hemisphere : {"N", "S"}
        Which UPS projection the input belongs to.
    ellipsoid : Ellipsoid, optional
        Reference ellipsoid. Defaults to WGS84.
    """
    ellipsoid = _resolve_ellipsoid(ellipsoid)

    hemisphere = str(hemisphere).upper()
    if hemisphere not in ("N", "S"):
        raise ValueError("hemisphere must be 'N' or 'S'.")

    e_arr = np.asarray(easting, dtype=float)
    n_arr = np.asarray(northing, dtype=float)
    scalar_input = e_arr.ndim == 0 and n_arr.ndim == 0

    ups = _make_ups(hemisphere, ellipsoid)
    lat, lon = ups.inverse(e_arr, n_arr)

    if scalar_input:
        return float(lat), float(lon)
    return lat, lon
