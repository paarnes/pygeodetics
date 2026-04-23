"""
Tests for the Universal Polar Stereographic (UPS) wrapper.
"""

import numpy as np
import pytest

from pygeodetics import geodetic_to_ups, ups_to_geodetic
from pygeodetics.Ellipsoid import WGS84


# Reference values from pyproj (EPSG:32661 north / EPSG:32761 south).
NORTH_REFERENCE = [
    # (lat, lon, easting, northing)
    (85.0, 0.0,   2_000_000.0,        1_444_542.6086173225),
    (85.0, 90.0,  2_555_457.3913826775, 2_000_000.0),
    (85.0, 45.0,  2_392_767.6881068815, 1_607_232.3118931185),
]
SOUTH_REFERENCE = [
    (-85.0, 0.0,   2_000_000.0,        2_555_457.3913826775),
    (-85.0, 45.0,  2_392_767.6881068815, 2_392_767.6881068815),
    (-85.0, -90.0, 1_444_542.6086173225, 2_000_000.0),
]


@pytest.mark.parametrize("lat, lon, e_ref, n_ref", NORTH_REFERENCE)
def test_ups_north_forward_matches_pyproj(lat, lon, e_ref, n_ref):
    e, n, hemi = geodetic_to_ups(lat, lon)
    assert hemi == "N"
    assert e == pytest.approx(e_ref, abs=1e-3)
    assert n == pytest.approx(n_ref, abs=1e-3)


@pytest.mark.parametrize("lat, lon, e_ref, n_ref", SOUTH_REFERENCE)
def test_ups_south_forward_matches_pyproj(lat, lon, e_ref, n_ref):
    e, n, hemi = geodetic_to_ups(lat, lon)
    assert hemi == "S"
    assert e == pytest.approx(e_ref, abs=1e-3)
    assert n == pytest.approx(n_ref, abs=1e-3)


@pytest.mark.parametrize("lat, lon", [(85.0, 0.0), (87.5, 23.0),
                                      (-85.0, 45.0), (-89.0, -120.0)])
def test_ups_round_trip_sub_mm(lat, lon):
    e, n, hemi = geodetic_to_ups(lat, lon)
    lat_back, lon_back = ups_to_geodetic(e, n, hemisphere=hemi)
    assert lat_back == pytest.approx(lat, abs=1e-10)
    assert lon_back == pytest.approx(lon, abs=1e-9)


def test_ups_north_pole_is_false_origin():
    e, n, hemi = geodetic_to_ups(90.0, 0.0)
    assert hemi == "N"
    assert e == pytest.approx(2_000_000.0, abs=1e-6)
    assert n == pytest.approx(2_000_000.0, abs=1e-6)


def test_ups_south_pole_is_false_origin():
    e, n, hemi = geodetic_to_ups(-90.0, 0.0)
    assert hemi == "S"
    assert e == pytest.approx(2_000_000.0, abs=1e-6)
    assert n == pytest.approx(2_000_000.0, abs=1e-6)


def test_ups_pole_inverse_returns_pole():
    lat, lon = ups_to_geodetic(2_000_000.0, 2_000_000.0, hemisphere="N")
    assert lat == pytest.approx(90.0, abs=1e-9)
    # longitude is undefined at the pole; we just require no crash / NaN.
    assert np.isfinite(lon)


def test_ups_array_input_round_trip():
    lats = np.array([85.0, 86.5, 88.0])
    lons = np.array([0.0, 30.0, 90.0])
    e, n, hemi = geodetic_to_ups(lats, lons)
    assert e.shape == lats.shape
    assert np.all(hemi == "N")
    lat_back, lon_back = ups_to_geodetic(e, n, hemisphere="N")
    np.testing.assert_allclose(lat_back, lats, atol=1e-10)
    np.testing.assert_allclose(lon_back, lons, atol=1e-9)


def test_ups_invalid_hemisphere():
    with pytest.raises(ValueError):
        ups_to_geodetic(2_000_000.0, 2_000_000.0, hemisphere="X")


def test_ups_custom_ellipsoid():
    e1, n1, _ = geodetic_to_ups(85.0, 30.0)
    e2, n2, _ = geodetic_to_ups(85.0, 30.0, ellipsoid=WGS84())
    assert e1 == pytest.approx(e2, abs=1e-9)
    assert n1 == pytest.approx(n2, abs=1e-9)
