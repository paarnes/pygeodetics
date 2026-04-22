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
