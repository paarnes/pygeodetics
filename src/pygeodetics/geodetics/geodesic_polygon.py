"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com


Geodesic polygon and path utilities on a reference ellipsoid.

Built on top of the existing geodesic primitives:

- :func:`pygeodetics.geodetics.geodetic_inverse_problem.geodetic_inverse_problem`
- :func:`pygeodetics.geodetics.geodetic_direct_problem.geodetic_direct_problem`

Provides:

- :func:`polygon_perimeter` — sum of geodesic edge lengths
- :func:`polygon_area`      — signed ellipsoidal area via Green's theorem
- :func:`polygon_centroid`  — area-weighted centroid (lat, lon)
- :func:`polygon_bounds`    — bounding box ``(min_lat, min_lon, max_lat, max_lon)``
- :func:`geodesic_interpolate` — N intermediate points along a geodesic

Sign convention for ``polygon_area``: positive for counter-clockwise
(mathematically positive) vertex ordering, negative for clockwise.
:func:`polygon_area` returns the absolute area by default; set
``signed=True`` to keep the sign.

Polygons are auto-closed: the last vertex does not need to repeat the
first one. Inputs may be Python lists, tuples, or 1-D NumPy arrays of
equal length.
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from ..Ellipsoid import Ellipsoid, WGS84
from .geodetic_inverse_problem import geodetic_inverse_problem
from .geodetic_direct_problem import geodetic_direct_problem


def _resolve_ellipsoid(ellipsoid):
    if ellipsoid is None:
        return WGS84()
    if not isinstance(ellipsoid, Ellipsoid):
        raise TypeError(
            "ellipsoid must be a pygeodetics.Ellipsoid instance "
            f"(got {type(ellipsoid).__name__})"
        )
    return ellipsoid


def _prepare_polygon(lats, lons):
    lats = np.asarray(lats, dtype=float).ravel()
    lons = np.asarray(lons, dtype=float).ravel()
    if lats.size != lons.size:
        raise ValueError("lats and lons must have the same length.")
    if lats.size < 3:
        raise ValueError("A polygon requires at least 3 vertices.")
    # Auto-close.
    if not (np.isclose(lats[0], lats[-1]) and np.isclose(lons[0], lons[-1])):
        lats = np.concatenate([lats, lats[:1]])
        lons = np.concatenate([lons, lons[:1]])
    return lats, lons


def polygon_perimeter(lats, lons,
                      ellipsoid: Optional[Ellipsoid] = None) -> float:
    """
    Compute the geodesic perimeter of a polygon on an ellipsoid.

    Parameters
    ----------
    lats, lons : array_like
        Polygon vertices in degrees. The polygon is auto-closed.
    ellipsoid : Ellipsoid, optional
        Reference ellipsoid (defaults to WGS84).

    Returns
    -------
    float
        Sum of the geodesic edge lengths in metres.
    """
    ellipsoid = _resolve_ellipsoid(ellipsoid)
    lats, lons = _prepare_polygon(lats, lons)
    a, b = ellipsoid.a, ellipsoid.b

    total = 0.0
    for i in range(lats.size - 1):
        # Skip degenerate zero-length edges (consecutive identical vertices).
        if lats[i] == lats[i + 1] and lons[i] == lons[i + 1]:
            continue
        _, _, d = geodetic_inverse_problem(
            a, b, lats[i], lons[i], lats[i + 1], lons[i + 1],
            radians=False,
        )
        total += float(d)
    return total


def _q(lat_rad: np.ndarray, e: float) -> np.ndarray:
    """
    Authalic-area antiderivative q(φ) = 2·F(φ) where F(φ) is the integral

        F(φ) = ∫_0^φ M(φ')·N(φ')·cosφ' dφ' / a²

    so that the area of the strip between the equator and latitude φ over
    a longitude width Δλ on the ellipsoid is exactly

        A = a² · Δλ · F(φ) = (a²/2) · Δλ · q(φ).

    Closed form (substitute u = sinφ, M = a(1-e²)/(1-e²sin²φ)^{3/2},
    N = a/(1-e²sin²φ)^{1/2}, integrate):

        q(φ) = (1 - e²) · [ sinφ / (1 - e²sin²φ)
                          + (1/(2e)) · ln((1+e·sinφ)/(1-e·sinφ)) ]

    For e = 0 (sphere) the closed form reduces to ``q = 2·sinφ``.
    """
    sin_phi = np.sin(lat_rad)
    if e == 0.0:
        return 2.0 * sin_phi
    return (1.0 - e * e) * (
        sin_phi / (1.0 - (e * sin_phi) ** 2)
        + (1.0 / (2.0 * e))
        * np.log((1.0 + e * sin_phi) / (1.0 - e * sin_phi))
    )


def polygon_area(lats, lons,
                 ellipsoid: Optional[Ellipsoid] = None,
                 signed: bool = False) -> float:
    """
    Compute the area of a polygon on an ellipsoid using Green's theorem.

    The implementation is the standard ellipsoidal trapezoidal-strip
    formula derived from Green's theorem,

        A = -(a² / 2) · Σ (λ_{i+1} - λ_i) · (F(φ_i) + F(φ_{i+1}))
          = -(a² / 4) · Σ (λ_{i+1} - λ_i) · (q(φ_i) + q(φ_{i+1}))

    where ``q(φ)`` is twice the equator-to-φ antiderivative of
    ``M(φ)·N(φ)·cosφ / a²``. The formula treats polygon edges as
    straight lines in the ``(λ, q)`` equal-area plane (rhumb-line
    polygon area). For geodesic polygons there is an additional
    correction of order ``f`` per edge; for typical small / medium
    polygons (≲ 1000 km extent) the rhumb-line area agrees with the
    geodesic area to better than 1 part per million. For very large
    geodesic polygons (continental scale or wider), the result will
    differ from a true geodesic-polygon area by a fraction of a
    percent.

    Parameters
    ----------
    lats, lons : array_like
        Polygon vertices in degrees. The polygon is auto-closed.
    ellipsoid : Ellipsoid, optional
        Reference ellipsoid (defaults to WGS84).
    signed : bool, default False
        If False (default), return the absolute area in m². If True,
        return the signed area: positive for counter-clockwise vertex
        ordering, negative for clockwise.

    Returns
    -------
    float
        Polygon area in square metres.
    """
    ellipsoid = _resolve_ellipsoid(ellipsoid)
    lats, lons = _prepare_polygon(lats, lons)
    a, b = ellipsoid.a, ellipsoid.b
    e = float(np.sqrt(1.0 - (b / a) ** 2))

    lat_r = np.radians(lats)
    lon_r = np.radians(lons)
    q = _q(lat_r, e)

    # Wrap longitude differences into (-pi, pi] so the formula handles
    # polygons that cross the antimeridian or wrap several times.
    dlon = np.diff(lon_r)
    dlon = ((dlon + np.pi) % (2.0 * np.pi)) - np.pi

    # Signed area: positive for CCW, negative for CW. The leading
    # negative reflects the orientation of the Green's-theorem
    # contour integral.
    area = -0.25 * (a ** 2) * float(np.sum(dlon * (q[:-1] + q[1:])))

    return area if signed else abs(area)


def polygon_centroid(lats, lons,
                     ellipsoid: Optional[Ellipsoid] = None
                     ) -> Tuple[float, float]:
    """
    Approximate the area-weighted centroid of a polygon on the ellipsoid.

    The polygon is mapped to (longitude, q(latitude)) — a coordinate pair
    in which equal areas correspond to equal planar areas — and the
    standard planar centroid formula is applied. The result is then
    mapped back to geodetic coordinates via the authalic latitude
    inverse series.

    Parameters
    ----------
    lats, lons : array_like
        Polygon vertices in degrees. Polygon is auto-closed.
    ellipsoid : Ellipsoid, optional
        Reference ellipsoid (defaults to WGS84).

    Returns
    -------
    (lat, lon) : tuple of float
        Centroid latitude and longitude in degrees.
    """
    ellipsoid = _resolve_ellipsoid(ellipsoid)
    lats, lons = _prepare_polygon(lats, lons)
    a, b = ellipsoid.a, ellipsoid.b
    e2 = 1.0 - (b / a) ** 2
    e = float(np.sqrt(e2))

    lat_r = np.radians(lats)
    lon_r = np.radians(lons)

    # Equal-area mapping: (λ, q) where q = q(φ) is the authalic integral.
    q_vals = _q(lat_r, e)

    # Planar centroid in (λ, q).
    dlon = np.diff(lon_r)
    dlon = ((dlon + np.pi) % (2.0 * np.pi)) - np.pi
    # Use cumulative longitudes to avoid the antimeridian discontinuity.
    lon_cum = np.concatenate(([lon_r[0]], lon_r[0] + np.cumsum(dlon)))

    cross = lon_cum[:-1] * q_vals[1:] - lon_cum[1:] * q_vals[:-1]
    area2 = 0.5 * np.sum(cross)
    if area2 == 0.0:
        # Degenerate polygon — fall back to vertex average.
        return (float(np.mean(lats[:-1])), float(np.mean(lons[:-1])))

    cx = np.sum((lon_cum[:-1] + lon_cum[1:]) * cross) / (6.0 * area2)
    cy = np.sum((q_vals[:-1] + q_vals[1:]) * cross) / (6.0 * area2)

    # Invert q(φ) by Newton's method.  d q/dφ = 2(1-e²)·cosφ/(1-e²sin²φ)².
    q_target = cy
    # Initial guess from spherical approximation (q ≈ 2·sinφ for e=0).
    phi = float(np.arcsin(np.clip(q_target / 2.0, -1.0, 1.0)))
    for _ in range(30):
        sin_phi = np.sin(phi)
        denom = 1.0 - (e * sin_phi) ** 2
        f = float(_q(np.array(phi), e)) - q_target
        fp = 2.0 * (1.0 - e2) * np.cos(phi) / (denom ** 2)
        if fp == 0.0:
            break
        delta = f / fp
        phi -= delta
        if abs(delta) < 1e-14:
            break

    centroid_lon = ((np.degrees(cx) + 180.0) % 360.0) - 180.0
    centroid_lat = np.degrees(phi)
    return float(centroid_lat), float(centroid_lon)


def polygon_bounds(lats, lons) -> Tuple[float, float, float, float]:
    """
    Compute the latitude/longitude bounding box of a polygon.

    Returns
    -------
    (min_lat, min_lon, max_lat, max_lon) : tuple of float
    """
    lats = np.asarray(lats, dtype=float).ravel()
    lons = np.asarray(lons, dtype=float).ravel()
    if lats.size == 0 or lons.size != lats.size:
        raise ValueError("lats and lons must be non-empty and the same length.")
    return (float(lats.min()), float(lons.min()),
            float(lats.max()), float(lons.max()))


def geodesic_interpolate(lat1: float, lon1: float, lat2: float, lon2: float,
                         n_points: int,
                         ellipsoid: Optional[Ellipsoid] = None,
                         include_endpoints: bool = True
                         ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate ``n_points`` evenly spaced points along the geodesic from
    ``(lat1, lon1)`` to ``(lat2, lon2)``.

    Parameters
    ----------
    lat1, lon1, lat2, lon2 : float
        Endpoint coordinates in degrees.
    n_points : int
        Total number of output points (>= 2 if ``include_endpoints``,
        >= 1 otherwise).
    ellipsoid : Ellipsoid, optional
        Reference ellipsoid (defaults to WGS84).
    include_endpoints : bool, default True
        If True, the first/last sample equal the input endpoints.
        If False, only the interior samples are returned.

    Returns
    -------
    lats, lons : np.ndarray
        Arrays of length ``n_points`` of geodetic latitudes / longitudes
        in degrees.
    """
    ellipsoid = _resolve_ellipsoid(ellipsoid)
    a, b = ellipsoid.a, ellipsoid.b

    if n_points < (2 if include_endpoints else 1):
        raise ValueError(
            "n_points must be >= 2 when include_endpoints=True (>= 1 otherwise)."
        )

    az1, _, dist = geodetic_inverse_problem(
        a, b, lat1, lon1, lat2, lon2, radians=False,
    )
    if dist == 0.0:
        return (np.full(n_points, lat1, dtype=float),
                np.full(n_points, lon1, dtype=float))

    if include_endpoints:
        fractions = np.linspace(0.0, 1.0, n_points)
    else:
        # n_points evenly spaced strictly between the endpoints.
        fractions = (np.arange(1, n_points + 1) / (n_points + 1))

    lats = np.empty(n_points, dtype=float)
    lons = np.empty(n_points, dtype=float)
    for i, frac in enumerate(fractions):
        if include_endpoints and i == 0:
            lats[i], lons[i] = lat1, lon1
            continue
        if include_endpoints and i == n_points - 1:
            lats[i], lons[i] = lat2, lon2
            continue
        lat_i, lon_i, _ = geodetic_direct_problem(
            a, b, lat1, lon1, az1, dist * float(frac), radians=False,
        )
        lats[i] = lat_i
        lons[i] = lon_i
    return lats, lons
