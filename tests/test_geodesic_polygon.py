"""
Tests for geodesic polygon utilities.
"""

import numpy as np
import pytest

from pygeodetics import (
    polygon_perimeter,
    polygon_area,
    polygon_centroid,
    polygon_bounds,
    polygon_densify,
    geodesic_interpolate,
    geodetic_inverse_problem,
)
from pygeodetics.Ellipsoid import WGS84


PERIM_CASES = [
    ([60.0, 60.0, 60.5, 60.5], [10.0, 11.0, 11.0, 10.0], 222171.32504913595),
    ([0.0, 0.0, 0.001, 0.001], [0.0, 0.001, 0.001, 0.0], 443.7875332131202),
    ([45.0, 45.0, 50.0, 50.0], [-75.0, -70.0, -70.0, -75.0], 1864388.9349665223),
    ([-22.0, -22.0, -23.0, -23.0], [-43.0, -42.0, -42.0, -43.0], 427258.50430609565),
]


@pytest.mark.parametrize("lats, lons, ref", PERIM_CASES)
def test_polygon_perimeter_reference(lats, lons, ref):
    assert polygon_perimeter(lats, lons) == pytest.approx(ref, rel=1e-9, abs=1e-3)


def test_perimeter_auto_close():
    lats_open = [60.0, 60.0, 60.5, 60.5]
    lons_open = [10.0, 11.0, 11.0, 10.0]
    lats_closed = lats_open + [lats_open[0]]
    lons_closed = lons_open + [lons_open[0]]
    assert polygon_perimeter(lats_open, lons_open) == pytest.approx(
        polygon_perimeter(lats_closed, lons_closed), abs=1e-9)


# (lats, lons, expected_area_m2, rel_tol)
AREA_CASES = [
    ([60.0, 60.0, 60.5, 60.5], [10.0, 11.0, 11.0, 10.0], 3084929054.2945557, 1e-4),
    ([0.0, 0.0, 0.01, 0.01], [0.0, 0.01, 0.01, 0.0], 1230907.2049932536, 1e-7),
    ([45.0, 45.0, 45.5, 45.5], [-75.0, -74.5, -74.5, -75.0], 2181131784.213135, 1e-4),
    ([-22.0, -22.0, -22.5, -22.5], [-43.0, -42.5, -42.5, -43.0], 2853612847.0240784, 1e-4),
]


@pytest.mark.parametrize("lats, lons, ref, rel_tol", AREA_CASES)
def test_polygon_area_reference(lats, lons, ref, rel_tol):
    assert polygon_area(lats, lons) == pytest.approx(ref, rel=rel_tol)


def test_polygon_area_signed():
    lats = [0.0, 0.0, 0.01, 0.01]
    lons = [0.0, 0.01, 0.01, 0.0]
    a_signed = polygon_area(lats, lons, signed=True)
    a_unsigned = polygon_area(lats, lons, signed=False)
    assert a_unsigned > 0
    assert abs(a_signed) == pytest.approx(a_unsigned, rel=1e-12)
    a_rev = polygon_area(list(reversed(lats)), list(reversed(lons)), signed=True)
    assert a_rev == pytest.approx(-a_signed, rel=1e-12)


def test_polygon_area_too_few_points():
    with pytest.raises(ValueError):
        polygon_area([0.0, 1.0], [0.0, 1.0])


def test_centroid_small_square_at_equator():
    lats = [0.0, 0.0, 0.01, 0.01]
    lons = [0.0, 0.01, 0.01, 0.0]
    lat_c, lon_c = polygon_centroid(lats, lons)
    assert lat_c == pytest.approx(0.005, abs=1e-6)
    assert lon_c == pytest.approx(0.005, abs=1e-6)


def test_centroid_norway_box_inside_polygon():
    lats = [60.0, 60.0, 60.5, 60.5]
    lons = [10.0, 11.0, 11.0, 10.0]
    lat_c, lon_c = polygon_centroid(lats, lons)
    assert 60.0 < lat_c < 60.5
    assert 10.0 < lon_c < 11.0
    # Equal-area centroid pulls slightly south of geometric mid-latitude.
    assert lat_c == pytest.approx(60.25, abs=0.01)
    assert lon_c == pytest.approx(10.5, abs=1e-6)


def test_polygon_bounds():
    lats = [60.0, 60.5, 61.0, 60.25]
    lons = [10.0, 11.0, 10.5, 9.5]
    lat_min, lon_min, lat_max, lon_max = polygon_bounds(lats, lons)
    assert lat_min == 60.0
    assert lat_max == 61.0
    assert lon_min == 9.5
    assert lon_max == 11.0


def test_geodesic_interpolate_endpoints():
    lat1, lon1 = 60.0, 10.0
    lat2, lon2 = 61.0, 11.0
    lats, lons = geodesic_interpolate(lat1, lon1, lat2, lon2, n_points=5)
    assert len(lats) == 5
    assert lats[0] == pytest.approx(lat1, abs=1e-12)
    assert lons[0] == pytest.approx(lon1, abs=1e-12)
    assert lats[-1] == pytest.approx(lat2, abs=1e-6)
    assert lons[-1] == pytest.approx(lon2, abs=1e-6)


def test_geodesic_interpolate_lies_on_geodesic():
    """Cumulative segment lengths should equal the total geodesic distance."""
    ell = WGS84()
    lat1, lon1 = 60.0, 10.0
    lat2, lon2 = 60.5, 11.5
    lats, lons = geodesic_interpolate(lat1, lon1, lat2, lon2, n_points=11)
    _, _, total = geodetic_inverse_problem(ell.a, ell.b, lat1, lon1, lat2, lon2)
    seg = 0.0
    for i in range(len(lats) - 1):
        _, _, d = geodetic_inverse_problem(
            ell.a, ell.b, lats[i], lons[i], lats[i + 1], lons[i + 1])
        seg += d
    assert seg == pytest.approx(total, rel=1e-9)


def test_geodesic_interpolate_invalid_n():
    with pytest.raises(ValueError):
        geodesic_interpolate(60.0, 10.0, 61.0, 11.0, n_points=1)


def test_polygon_array_inputs():
    lats = np.array([60.0, 60.0, 60.5, 60.5])
    lons = np.array([10.0, 11.0, 11.0, 10.0])
    a = polygon_area(lats, lons)
    p = polygon_perimeter(lats, lons)
    assert a > 0 and p > 0


def _max_segment_length(lats, lons):
    ell = WGS84()
    longest = 0.0
    for i in range(len(lats) - 1):
        _, _, d = geodetic_inverse_problem(
            ell.a, ell.b, lats[i], lons[i], lats[i + 1], lons[i + 1]
        )
        longest = max(longest, d)
    return longest


def test_densify_max_segment_respected():
    lats = [60.0, 60.0, 61.0, 61.0]
    lons = [10.0, 11.0, 11.0, 10.0]
    cap = 5000.0
    dlats, dlons = polygon_densify(lats, lons, max_segment_length=cap)
    # Add a small numerical tolerance for Vincenty round-off.
    assert _max_segment_length(dlats, dlons) <= cap + 1e-6


def test_densify_preserves_original_vertices_and_perimeter():
    lats = [60.0, 60.0, 61.0, 61.0]
    lons = [10.0, 11.0, 11.0, 10.0]
    cap = 5000.0
    dlats, dlons = polygon_densify(lats, lons, max_segment_length=cap)

    # Perimeter must be preserved (same edges, just sampled densely).
    p_orig = polygon_perimeter(lats, lons)
    p_dens = polygon_perimeter(dlats, dlons)
    assert p_dens == pytest.approx(p_orig, rel=0, abs=1e-3)

    # Each original vertex must appear exactly in the densified output.
    out_pts = set(zip(map(float, dlats), map(float, dlons)))
    for la, lo in zip(lats, lons):
        assert (float(la), float(lo)) in out_pts


def test_densify_closed_vs_open_ring():
    lats = [60.0, 60.0, 60.5, 60.5]
    lons = [10.0, 11.0, 11.0, 10.0]
    closed_lats, closed_lons = polygon_densify(lats, lons, 4000.0, close=True)
    open_lats, open_lons = polygon_densify(lats, lons, 4000.0, close=False)
    # Open is exactly one (closing) vertex shorter than closed.
    assert len(closed_lats) == len(open_lats) + 1
    # Closed ring's last point equals the first.
    assert closed_lats[0] == closed_lats[-1]
    assert closed_lons[0] == closed_lons[-1]


def test_densify_no_op_when_cap_exceeds_longest_edge():
    lats = [60.0, 60.0, 60.001, 60.001]
    lons = [10.0, 10.001, 10.001, 10.0]
    # Edges are ~110 m; cap is much larger -> no inserted vertices.
    dlats, dlons = polygon_densify(lats, lons, max_segment_length=10_000.0)
    # Output is the auto-closed input (4 vertices + 1 closing).
    assert len(dlats) == 5
    np.testing.assert_allclose(dlats[:4], lats)
    np.testing.assert_allclose(dlons[:4], lons)


def test_densify_inserted_points_lie_on_geodesic():
    """
    Each inserted vertex must lie on the geodesic of the edge it
    subdivides — i.e. the sum of geodesic distances from edge start to
    inserted point and from inserted point to edge end equals the
    original edge length.
    """
    lats = [59.0, 59.0, 60.0]
    lons = [10.0, 12.0, 12.0]
    cap = 25_000.0
    dlats, dlons = polygon_densify(lats, lons, max_segment_length=cap, close=False)

    ell = WGS84()
    # Find the first inserted point along edge 0 (between vertex 0 and 1)
    # and verify it lies on that geodesic.
    _, _, edge_len = geodetic_inverse_problem(
        ell.a, ell.b, lats[0], lons[0], lats[1], lons[1]
    )
    # The first inserted point is at index 1 in the output (index 0 is the
    # original start vertex). It must split the edge: dist(start, p) +
    # dist(p, end) == edge_len.
    _, _, d_start_to_p = geodetic_inverse_problem(
        ell.a, ell.b, lats[0], lons[0], dlats[1], dlons[1]
    )
    _, _, d_p_to_end = geodetic_inverse_problem(
        ell.a, ell.b, dlats[1], dlons[1], lats[1], lons[1]
    )
    assert d_start_to_p + d_p_to_end == pytest.approx(edge_len, rel=0, abs=1e-6)


def test_densify_invalid_inputs():
    lats = [60.0, 60.0, 61.0, 61.0]
    lons = [10.0, 11.0, 11.0, 10.0]
    with pytest.raises(ValueError):
        polygon_densify(lats, lons, max_segment_length=0.0)
    with pytest.raises(ValueError):
        polygon_densify(lats, lons, max_segment_length=-100.0)
    with pytest.raises(ValueError):
        polygon_densify(lats, lons, max_segment_length=float("nan"))


def test_densify_accepts_numpy_inputs():
    lats = np.array([60.0, 60.0, 60.5, 60.5])
    lons = np.array([10.0, 11.0, 11.0, 10.0])
    dlats, dlons = polygon_densify(lats, lons, 10_000.0)
    assert isinstance(dlats, np.ndarray)
    assert isinstance(dlons, np.ndarray)
    assert dlats.dtype == float and dlons.dtype == float


# Reference values for a WGS84 triangle [(59,10), (60,10), (60,11.5)]
# densified with max_segment_length = 50 000 m, generated independently
# with a different geodesic library (GeographicLib via pyproj.Geod).
# Our implementation must reproduce these to better than 1e-9 degrees
# (~0.1 mm on the ground), i.e. agree with the reference engine within
# floating-point round-off of the Vincenty iteration.
_REF_DENSIFY_LATS = [
    59.0,
    59.333350468120514,
    59.66668376426481,
    60.0,
    60.00212914232574,
    60.0,
    59.66860531053436,
    59.335251085945316,
    59.0,
]
_REF_DENSIFY_LONS = [
    10.0,
    10.0,
    10.0,
    10.0,
    10.75,
    11.5,
    10.990094730652336,
    10.490184377826507,
    10.0,
]


def test_densify_matches_independent_geodesic_reference():
    lats = [59.0, 60.0, 60.0]
    lons = [10.0, 10.0, 11.5]
    dlats, dlons = polygon_densify(lats, lons, max_segment_length=50_000.0)
    np.testing.assert_allclose(dlats, _REF_DENSIFY_LATS, atol=1e-9)
    np.testing.assert_allclose(dlons, _REF_DENSIFY_LONS, atol=1e-9)


# Reference values for the geodesic-area mode of polygon_area, generated
# independently with a different geodesic library (GeographicLib via
# pyproj.Geod(ellps="WGS84").geometry_area_perimeter on a shapely Polygon).
# Each tuple: (lats, lons, expected_geodesic_area_m2).
GEODESIC_AREA_CASES = [
    # Continental US bounding box (~58 deg wide, 24 deg tall): rhumb error
    # is ~1.2% here, geodesic mode must match the reference to ~1 ppm.
    ([49.0, 49.0, 25.0, 25.0], [-125.0, -67.0, -67.0, -125.0],
     13490379239003.164),
    # South America bounding box (~48 deg wide, 67 deg tall).
    ([12.0, 12.0, -55.0, -55.0], [-82.0, -34.0, -34.0, -82.0],
     35816818445838.03),
    # Hemispheric strip (160 deg wide, 160 deg tall).
    ([80.0, 80.0, -80.0, -80.0], [0.0, 90.0, 90.0, 0.0],
     126262770415437.17),
    # Small box where rhumb and geodesic effectively coincide; geodesic
    # mode must not regress on small polygons.
    ([60.0, 60.0, 60.5, 60.5], [10.0, 11.0, 11.0, 10.0],
     3084929054.2945557),
]


@pytest.mark.parametrize("lats, lons, ref", GEODESIC_AREA_CASES)
def test_polygon_area_geodesic_mode_matches_independent_reference(lats, lons, ref):
    area = polygon_area(lats, lons, geodesic=True)
    assert area == pytest.approx(ref, rel=1e-5)


def test_polygon_area_geodesic_mode_improves_continental_polygon():
    # Continental US bbox: rhumb formula is off by ~1%, geodesic mode
    # must agree with the GeographicLib reference far more tightly.
    lats = [49.0, 49.0, 25.0, 25.0]
    lons = [-125.0, -67.0, -67.0, -125.0]
    ref = 13490379239003.164
    err_rhumb = abs(polygon_area(lats, lons, geodesic=False) - ref) / ref
    err_geod = abs(polygon_area(lats, lons, geodesic=True) - ref) / ref
    assert err_rhumb > 1e-3
    assert err_geod < 1e-5
    assert err_geod < err_rhumb / 100.0


def test_polygon_area_geodesic_mode_signed_orientation():
    lats = [49.0, 49.0, 25.0, 25.0]
    lons = [-125.0, -67.0, -67.0, -125.0]
    a_unsigned = polygon_area(lats, lons, geodesic=True)
    a_ccw = polygon_area(list(reversed(lats)), list(reversed(lons)),
                        geodesic=True, signed=True)
    a_cw = polygon_area(lats, lons, geodesic=True, signed=True)
    assert a_unsigned > 0
    assert a_cw < 0
    assert a_ccw > 0
    assert abs(a_cw) == pytest.approx(a_unsigned, rel=1e-12)
    assert a_ccw == pytest.approx(-a_cw, rel=1e-12)


def test_polygon_area_geodesic_mode_segment_length_convergence():
    # Tighter densification should not move the answer materially once
    # we're already at the 10 km default; both should match the reference.
    lats = [12.0, 12.0, -55.0, -55.0]
    lons = [-82.0, -34.0, -34.0, -82.0]
    ref = 35816818445838.03
    a_default = polygon_area(lats, lons, geodesic=True)
    a_tight = polygon_area(lats, lons, geodesic=True, max_segment_length=2_000.0)
    assert a_default == pytest.approx(ref, rel=1e-5)
    assert a_tight == pytest.approx(ref, rel=1e-5)
    assert a_tight == pytest.approx(a_default, rel=1e-6)
