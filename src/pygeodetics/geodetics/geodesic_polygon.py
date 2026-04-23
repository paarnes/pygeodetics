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
- :func:`polygon_densify`   — add geodesic vertices so no edge exceeds a max length
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
                 signed: bool = False,
                 geodesic: bool = False,
                 max_segment_length: float = 10_000.0) -> float:
    """
    Compute the area of a polygon on an ellipsoid using Green's theorem.

    The implementation is the standard ellipsoidal trapezoidal-strip
    formula derived from Green's theorem,

        A = -(a² / 2) · Σ (λ_{i+1} - λ_i) · (F(φ_i) + F(φ_{i+1}))
          = -(a² / 4) · Σ (λ_{i+1} - λ_i) · (q(φ_i) + q(φ_{i+1}))

    where ``q(φ)`` is twice the equator-to-φ antiderivative of
    ``M(φ)·N(φ)·cosφ / a²``. The formula treats polygon edges as
    straight lines in the ``(λ, q)`` equal-area plane, which is
    equivalent to treating each edge as a *rhumb line*.

    For most polygons this matches the true geodesic-polygon area to
    better than 1 ppm. For very large polygons (continental scale or
    wider), the rhumb / geodesic discrepancy can reach a fraction of a
    percent — set ``geodesic=True`` to get the true geodesic-polygon
    area instead.

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
    geodesic : bool, default False
        If True, treat polygon edges as true geodesics (great-ellipse
        arcs) instead of rhumb lines. The polygon is internally
        densified with :func:`polygon_densify` so that no edge segment
        is longer than ``max_segment_length`` metres before applying
        the Green's-theorem strip formula. As segments shrink, the
        rhumb-edge approximation converges to the geodesic area; with
        the default 10 km cap the residual error is < 1 ppm even for
        hemispheric polygons.
    max_segment_length : float, default 10 000
        Maximum geodesic segment length (in metres) used for the
        internal densification when ``geodesic=True``. Ignored when
        ``geodesic=False``. Smaller values give a tighter match to the
        true geodesic area but increase the cost roughly linearly.

    Returns
    -------
    float
        Polygon area in square metres.
    """
    ellipsoid = _resolve_ellipsoid(ellipsoid)

    if geodesic:
        # Densify to small enough segments that the rhumb-line edges of
        # the densified polygon are an excellent local approximation to
        # the true geodesic edges. The total area then converges to the
        # geodesic-polygon area.
        lats, lons = polygon_densify(
            lats, lons,
            max_segment_length=max_segment_length,
            ellipsoid=ellipsoid,
            close=True,
        )

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


def polygon_densify(lats, lons,
                    max_segment_length: float,
                    ellipsoid: Optional[Ellipsoid] = None,
                    close: bool = True,
                    ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Densify a polygon by inserting geodesic vertices along every edge so
    that no two consecutive output vertices are farther apart than
    ``max_segment_length`` metres along the ellipsoid.

    The added vertices lie exactly on the geodesic connecting the two
    original vertices of each edge — they are obtained by sampling the
    geodesic at evenly spaced fractions of its arc length using the
    Vincenty direct problem. This preserves the true ellipsoidal shape
    of the boundary, which is important when projecting the polygon to a
    flat map (where straight cartographic edges between widely spaced
    vertices distort the original geometry).

    Parameters
    ----------
    lats, lons : array_like
        Polygon vertices in degrees. Order is preserved. The polygon is
        treated as closed: a closing edge from the last to the first
        vertex is also densified (and its closing vertex is included in
        the output if ``close=True``).
    max_segment_length : float
        Maximum allowed geodesic length (in metres) between consecutive
        output vertices. Must be > 0. Existing edges shorter than this
        are left untouched.
    ellipsoid : Ellipsoid, optional
        Reference ellipsoid (defaults to WGS84).
    close : bool, default True
        If True, the returned ring is explicitly closed (last vertex
        equals the first). If False, the returned arrays end at the
        densified version of the last input vertex (open ring).

    Returns
    -------
    lats_out, lons_out : np.ndarray
        Densified vertex arrays (1-D, dtype ``float``) in degrees.

    Raises
    ------
    ValueError
        If ``max_segment_length`` is not strictly positive, or if fewer
        than three input vertices are supplied.

    Notes
    -----
    - The number of intermediate points inserted on an edge of length
      ``d`` is ``ceil(d / max_segment_length) - 1``, so each resulting
      sub-segment has length ``d / ceil(d / max_segment_length)``,
      which is ``<= max_segment_length`` by construction.
    - The original input vertices are always preserved exactly; only
      new points are inserted between them.
    - For an open polyline (not a polygon), use
      :func:`geodesic_interpolate` per segment instead.
    """
    ellipsoid = _resolve_ellipsoid(ellipsoid)
    a, b = ellipsoid.a, ellipsoid.b

    if not np.isfinite(max_segment_length) or max_segment_length <= 0.0:
        raise ValueError("max_segment_length must be a positive, finite number of metres.")

    lats_c, lons_c = _prepare_polygon(lats, lons)  # auto-closed ring
    n_edges = len(lats_c) - 1

    # First pass: compute every edge's geodesic length and starting azimuth,
    # and how many sub-segments each edge needs. We need all of this up front
    # so we can pre-allocate the output and avoid any Python-list growth.
    az1s = np.empty(n_edges, dtype=float)
    dists = np.empty(n_edges, dtype=float)
    n_segs = np.empty(n_edges, dtype=np.int64)
    for i in range(n_edges):
        az1, _, dist = geodetic_inverse_problem(
            a, b,
            float(lats_c[i]), float(lons_c[i]),
            float(lats_c[i + 1]), float(lons_c[i + 1]),
            radians=False,
        )
        az1s[i] = az1
        dists[i] = dist
        n_segs[i] = max(1, int(np.ceil(dist / max_segment_length)))

    # Total output length: each edge contributes n_seg points (n_seg-1
    # interior + 1 end vertex), plus the very first start vertex.
    total = int(n_segs.sum()) + 1
    out_lats = np.empty(total, dtype=float)
    out_lons = np.empty(total, dtype=float)
    out_lats[0] = float(lats_c[0])
    out_lons[0] = float(lons_c[0])

    # Second pass: for each edge, compute all interior points in a single
    # vectorised call to the Vincenty direct solver, then write the slice
    # into the pre-allocated output (interior points + the original end
    # vertex of the edge).
    write = 1
    for i in range(n_edges):
        n_seg = int(n_segs[i])
        end_lat = float(lats_c[i + 1])
        end_lon = float(lons_c[i + 1])

        if n_seg > 1:
            # Distances at fractions 1/n_seg ... (n_seg-1)/n_seg of the edge.
            sample_dists = dists[i] * (np.arange(1, n_seg, dtype=float) / n_seg)
            lat_int, lon_int, _ = geodetic_direct_problem(
                a, b,
                float(lats_c[i]), float(lons_c[i]), float(az1s[i]),
                sample_dists,
                radians=False,
            )
            n_int = n_seg - 1
            out_lats[write:write + n_int] = lat_int
            out_lons[write:write + n_int] = lon_int
            write += n_int

        # Always preserve the edge's end vertex exactly (the original input).
        out_lats[write] = end_lat
        out_lons[write] = end_lon
        write += 1

    if not close:
        # Drop the duplicated closing vertex.
        out_lats = out_lats[:-1]
        out_lons = out_lons[:-1]

    return out_lats, out_lons
