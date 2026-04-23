"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com


Universal Transverse Mercator (UTM) convenience API.

This module is a thin, user-facing wrapper around the existing
:class:`pygeodetics.projections.TransverseMercator.TransverseMercator`
projection engine and the
:class:`pygeodetics.utils.ProjectionZoneCalc.ProjectionZoneCalc` zone
helper. It applies the standard UTM defaults so the user can convert
between geodetic coordinates and UTM grid coordinates in a single call:

- scale factor at the central meridian: 0.9996
- false easting: 500 000 m
- false northing: 0 m (northern hemisphere) / 10 000 000 m (southern hemisphere)
- 6° wide longitudinal zones, numbered 1..60 starting at 180°W
- latitude of natural origin: equator (0°)

No new geodetic mathematics is implemented here.

Usage Example
-------------
    >>> from pygeodetics import geodetic_to_utm, utm_to_geodetic
    >>> easting, northing, zone, band = geodetic_to_utm(60.0, 10.75)
    >>> lat, lon = utm_to_geodetic(easting, northing, zone=zone, hemisphere="N")
"""

from __future__ import annotations

from typing import Literal, Optional, Tuple, Union

import numpy as np

from ..Ellipsoid import Ellipsoid, WGS84
from ..utils.ProjectionZoneCalc import ProjectionZoneCalc
from .TransverseMercator import TransverseMercator


# Standard UTM defaults (EPSG / NGA).
UTM_SCALE_FACTOR = 0.9996
UTM_FALSE_EASTING = 500_000.0
UTM_FALSE_NORTHING_NORTH = 0.0
UTM_FALSE_NORTHING_SOUTH = 10_000_000.0

# MGRS latitude band letters, 8° tall, from 80°S northward. The final band
# (``X``) covers 72°N..84°N (12° tall). The letters ``I`` and ``O`` are
# intentionally skipped to avoid confusion with the digits 1 and 0.
_MGRS_BAND_LETTERS = "CDEFGHJKLMNPQRSTUVWX"

ArrayLike = Union[float, int, np.ndarray, list, tuple]


def mgrs_band_letter(lat: ArrayLike) -> Union[str, np.ndarray]:
    """
    Return the MGRS latitude band letter (C..X) for the given geodetic
    latitude(s) in degrees.

    Bands are 8° tall, starting at 80°S. The northernmost band ``X`` covers
    72°N..84°N (12° tall). The letters ``I`` and ``O`` are skipped.

    Parameters
    ----------
    lat : float or array_like
        Geodetic latitude(s) in degrees. Must satisfy ``-80 <= lat <= 84``.

    Returns
    -------
    str or np.ndarray of str
        MGRS band letter(s). Scalar input yields a scalar string; array
        input yields a NumPy array of single-character strings.

    Raises
    ------
    ValueError
        If any latitude lies outside the MGRS validity range
        (``-80 <= lat <= 84``). Polar regions outside this range are
        covered by UPS (Universal Polar Stereographic) and not by UTM.
    """
    lat_arr = np.asarray(lat, dtype=float)

    if np.any(lat_arr < -80.0) or np.any(lat_arr > 84.0):
        raise ValueError(
            "MGRS latitude band letters are only defined for "
            "-80° <= lat <= 84°. Use UPS for polar regions."
        )

    # Clamp index to the last band (X) so latitudes in [72, 84] all map to X.
    idx = np.minimum(np.floor((lat_arr + 80.0) / 8.0).astype(int), 19)
    letters = np.array(list(_MGRS_BAND_LETTERS))[idx]

    if lat_arr.ndim == 0:
        return str(letters)
    return letters


def _hemisphere_from_lat(lat: np.ndarray) -> np.ndarray:
    """Return an array of 'N'/'S' chars matching ``lat``."""
    return np.where(lat >= 0.0, "N", "S")


def _resolve_ellipsoid(ellipsoid: Optional[Ellipsoid]) -> Ellipsoid:
    if ellipsoid is None:
        return WGS84()
    if not isinstance(ellipsoid, Ellipsoid):
        raise TypeError(
            "ellipsoid must be a pygeodetics.Ellipsoid instance "
            f"(got {type(ellipsoid).__name__})"
        )
    return ellipsoid


def _make_tm(central_meridian_deg: float, false_northing: float,
             ellipsoid: Ellipsoid) -> TransverseMercator:
    """Build a Transverse Mercator engine configured for UTM."""
    return TransverseMercator(
        lat_origin=0.0,
        lon_origin=np.radians(central_meridian_deg),
        scale_factor=UTM_SCALE_FACTOR,
        false_easting=UTM_FALSE_EASTING,
        false_northing=false_northing,
        a=ellipsoid.a,
        f=ellipsoid.f,
    )


def geodetic_to_utm(
    lat: ArrayLike,
    lon: ArrayLike,
    ellipsoid: Optional[Ellipsoid] = None,
    force_zone: Optional[int] = None,
) -> Tuple[Union[float, np.ndarray],
           Union[float, np.ndarray],
           Union[int, np.ndarray],
           Union[str, np.ndarray]]:
    """
    Convert geodetic coordinates to UTM grid coordinates.

    The UTM zone is selected automatically from the longitude unless
    ``force_zone`` is supplied. The hemisphere (and therefore the false
    northing) is selected automatically from the latitude.

    Parameters
    ----------
    lat : float or array_like
        Geodetic latitude(s) in degrees. Valid range: -80°..84°.
    lon : float or array_like
        Geodetic longitude(s) in degrees.
    ellipsoid : Ellipsoid, optional
        Reference ellipsoid. Defaults to WGS84.
    force_zone : int, optional
        Override the automatically computed UTM zone (1..60). Useful near
        zone borders or for the Norway / Svalbard exceptions, which this
        wrapper does not apply automatically.

    Returns
    -------
    easting : float or np.ndarray
        UTM easting in metres.
    northing : float or np.ndarray
        UTM northing in metres (false northing 10 000 000 m applied for
        southern hemisphere points).
    zone : int or np.ndarray of int
        UTM zone number(s) (1..60).
    band : str or np.ndarray of str
        MGRS latitude band letter(s).

    Notes
    -----
    - Scalar inputs return Python scalars; array inputs return NumPy arrays.
    - When ``lat`` / ``lon`` are arrays spanning multiple zones or both
      hemispheres, the points are grouped by ``(zone, hemisphere)`` and
      projected with the correct central meridian / false northing per
      group.
    """
    ellipsoid = _resolve_ellipsoid(ellipsoid)

    lat_arr = np.asarray(lat, dtype=float)
    lon_arr = np.asarray(lon, dtype=float)
    scalar_input = lat_arr.ndim == 0 and lon_arr.ndim == 0

    lat_arr, lon_arr = np.broadcast_arrays(lat_arr, lon_arr)
    flat_lat = np.atleast_1d(lat_arr).ravel()
    flat_lon = np.atleast_1d(lon_arr).ravel()
    n = flat_lat.size

    zone_calc = ProjectionZoneCalc()
    if force_zone is not None:
        if not (1 <= int(force_zone) <= 60):
            raise ValueError("force_zone must be an integer in [1, 60].")
        zones = np.full(n, int(force_zone), dtype=int)
    else:
        zones = zone_calc.calculate_zone(flat_lon).astype(int)

    bands = mgrs_band_letter(flat_lat) if n > 0 else np.array([], dtype="<U1")
    bands = np.atleast_1d(bands)
    hemispheres = _hemisphere_from_lat(flat_lat)

    easting = np.empty(n, dtype=float)
    northing = np.empty(n, dtype=float)

    # Group by (zone, hemisphere) so we instantiate the TM engine once
    # per group rather than once per point.
    keys = np.array([f"{z}_{h}" for z, h in zip(zones, hemispheres)])
    for key in np.unique(keys):
        mask = keys == key
        zone = int(zones[mask][0])
        hemi = hemispheres[mask][0]
        false_n = (UTM_FALSE_NORTHING_NORTH if hemi == "N"
                   else UTM_FALSE_NORTHING_SOUTH)
        cm = zone_calc.calculate_central_meridian(zone)
        tm = _make_tm(cm, false_n, ellipsoid)
        coords = np.column_stack((flat_lon[mask], flat_lat[mask]))
        e, n_ = tm.geog_to_projected(coords, unit="deg")
        easting[mask] = e
        northing[mask] = n_

    if scalar_input:
        return (float(easting[0]), float(northing[0]),
                int(zones[0]), str(bands[0]))

    out_shape = lat_arr.shape
    return (easting.reshape(out_shape),
            northing.reshape(out_shape),
            zones.reshape(out_shape),
            bands.reshape(out_shape))


def utm_to_geodetic(
    easting: ArrayLike,
    northing: ArrayLike,
    zone: int,
    hemisphere: Literal["N", "S"] = "N",
    ellipsoid: Optional[Ellipsoid] = None,
) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray]]:
    """
    Convert UTM grid coordinates back to geodetic latitude and longitude.

    Parameters
    ----------
    easting : float or array_like
        UTM easting(s) in metres.
    northing : float or array_like
        UTM northing(s) in metres.
    zone : int
        UTM zone number (1..60). A single zone applies to all input points.
    hemisphere : {"N", "S"}, default "N"
        Hemisphere of the input points. Selects the false northing
        (0 m for "N", 10 000 000 m for "S").
    ellipsoid : Ellipsoid, optional
        Reference ellipsoid. Defaults to WGS84.

    Returns
    -------
    lat : float or np.ndarray
        Geodetic latitude(s) in degrees.
    lon : float or np.ndarray
        Geodetic longitude(s) in degrees.
    """
    ellipsoid = _resolve_ellipsoid(ellipsoid)

    if not (1 <= int(zone) <= 60):
        raise ValueError("zone must be an integer in [1, 60].")

    hemisphere = str(hemisphere).upper()
    if hemisphere not in ("N", "S"):
        raise ValueError("hemisphere must be 'N' or 'S'.")

    easting_arr = np.asarray(easting, dtype=float)
    northing_arr = np.asarray(northing, dtype=float)
    scalar_input = easting_arr.ndim == 0 and northing_arr.ndim == 0

    easting_arr, northing_arr = np.broadcast_arrays(easting_arr, northing_arr)
    flat_e = np.atleast_1d(easting_arr).ravel()
    flat_n = np.atleast_1d(northing_arr).ravel()

    false_n = (UTM_FALSE_NORTHING_NORTH if hemisphere == "N"
               else UTM_FALSE_NORTHING_SOUTH)
    cm = ProjectionZoneCalc().calculate_central_meridian(int(zone))
    tm = _make_tm(cm, false_n, ellipsoid)

    coords = np.column_stack((flat_e, flat_n))
    lon_deg, lat_deg = tm.projected_to_geog(coords, unit="deg")

    if scalar_input:
        return float(lat_deg[0]), float(lon_deg[0])

    out_shape = easting_arr.shape
    return lat_deg.reshape(out_shape), lon_deg.reshape(out_shape)
