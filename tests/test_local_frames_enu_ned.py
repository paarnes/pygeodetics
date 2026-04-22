"""
Tests for bidirectional ENU / NED conversions.
"""

import numpy as np
import pytest

from pygeodetics import (
    geodetic2enu,
    enu2geodetic,
    enu2ecef,
    geodetic2ned,
    ned2geodetic,
    ned2ecef,
)
from pygeodetics.Ellipsoid import WGS84, GRS80


# (lat0, lon0, h0, lat, lon, h)
CASES = [
    (60.0, 10.75, 100.0,  60.1, 10.9, 250.0),
    (-33.8688, 151.2093, 50.0,  -34.0, 151.3, 200.0),
    (0.0, 0.0, 0.0,  0.001, 0.001, 10.0),
    (45.0, -75.0, 200.0,  46.0, -74.0, 1000.0),
    (-22.9068, -43.1729, 5.0,  -22.95, -43.20, 800.0),
    (78.2232, 15.6267, 12.0,  78.5, 16.0, 100.0),
]

ENU_REF = [
    (8345.0441139962, 11151.2044941953, 134.8151449579),
    (8379.5590755459, -14557.0031728094, 127.8302755237),
    (111.3196653037, 110.5744503490, 9.9980636081),
    (77471.4898437192, 111631.4213016773, -648.1122752632),
    (-2779.7313267054, -4784.9566248038, 792.5906659039),
    (8311.6762240596, 30931.0473310455, 7.8222511830),
]

ECEF_REF = [
    (3130067.722908, 602756.037498, 5506256.092066),
    (-4643106.679524, 2542026.655197, -3546558.402361),
    (6378146.998064, 111.319665, 110.574450),
    (1223558.292002, -4267054.860409, 4565966.880633),
    (4284153.442049, -4023087.866597, -2471933.325044),
    (1226303.504313, 351636.871401, 6228402.467347),
]

NED_REF = [
    (11151.2044941953, 8345.0441139962, -134.8151449579),
    (-14557.0031728094, 8379.5590755459, -127.8302755237),
    (110.5744503490, 111.3196653037, -9.9980636081),
    (111631.4213016773, 77471.4898437192, 648.1122752632),
    (-4784.9566248038, -2779.7313267054, -792.5906659039),
    (30931.0473310455, 8311.6762240596, -7.8222511830),
]


@pytest.mark.parametrize("case, ref", list(zip(CASES, ENU_REF)))
def test_geodetic2enu_reference(case, ref):
    lat0, lon0, h0, lat, lon, h = case
    e, n, u = geodetic2enu(lat, lon, h, lat0, lon0, h0)
    assert e == pytest.approx(ref[0], abs=1e-6)
    assert n == pytest.approx(ref[1], abs=1e-6)
    assert u == pytest.approx(ref[2], abs=1e-6)


@pytest.mark.parametrize("case, ref", list(zip(CASES, ECEF_REF)))
def test_enu2ecef_reference(case, ref):
    lat0, lon0, h0, lat, lon, h = case
    e, n, u = geodetic2enu(lat, lon, h, lat0, lon0, h0)
    X, Y, Z = enu2ecef(e, n, u, lat0, lon0, h0)
    assert X == pytest.approx(ref[0], abs=1e-3)
    assert Y == pytest.approx(ref[1], abs=1e-3)
    assert Z == pytest.approx(ref[2], abs=1e-3)


@pytest.mark.parametrize("case", CASES)
def test_geodetic_enu_round_trip(case):
    lat0, lon0, h0, lat, lon, h = case
    e, n, u = geodetic2enu(lat, lon, h, lat0, lon0, h0)
    lat_b, lon_b, h_b = enu2geodetic(e, n, u, lat0, lon0, h0)
    assert lat_b == pytest.approx(lat, abs=1e-10)
    assert lon_b == pytest.approx(lon, abs=1e-10)
    assert h_b == pytest.approx(h, abs=1e-6)


@pytest.mark.parametrize("case, ref", list(zip(CASES, NED_REF)))
def test_geodetic2ned_reference(case, ref):
    lat0, lon0, h0, lat, lon, h = case
    n_, e_, d_ = geodetic2ned(lat, lon, h, lat0, lon0, h0)
    assert n_ == pytest.approx(ref[0], abs=1e-6)
    assert e_ == pytest.approx(ref[1], abs=1e-6)
    assert d_ == pytest.approx(ref[2], abs=1e-6)


@pytest.mark.parametrize("case", CASES)
def test_geodetic_ned_round_trip(case):
    lat0, lon0, h0, lat, lon, h = case
    n_, e_, d_ = geodetic2ned(lat, lon, h, lat0, lon0, h0)
    lat_b, lon_b, h_b = ned2geodetic(n_, e_, d_, lat0, lon0, h0)
    assert lat_b == pytest.approx(lat, abs=1e-10)
    assert lon_b == pytest.approx(lon, abs=1e-10)
    assert h_b == pytest.approx(h, abs=1e-6)


def test_ned2ecef_reference():
    lat0, lon0, h0 = 50.0, 8.0, 0.0
    n_, e_, d_ = 1234.5, -567.8, 90.1
    X, Y, Z = ned2ecef(n_, e_, d_, lat0, lon0, h0)
    assert X == pytest.approx(4066971.832648, abs=1e-3)
    assert Y == pytest.approx(571002.235976, abs=1e-3)
    assert Z == pytest.approx(4863513.538406, abs=1e-3)


def test_array_inputs_round_trip():
    lat0, lon0, h0 = 60.0, 10.0, 0.0
    lats = np.array([60.0, 60.1, 60.2, 60.3])
    lons = np.array([10.0, 10.05, 10.10, 10.15])
    hs = np.array([0.0, 50.0, 100.0, 150.0])
    e, n, u = geodetic2enu(lats, lons, hs, lat0, lon0, h0)
    assert isinstance(e, np.ndarray)
    assert e.shape == lats.shape
    lat_b, lon_b, h_b = enu2geodetic(e, n, u, lat0, lon0, h0)
    np.testing.assert_allclose(lat_b, lats, atol=1e-10)
    np.testing.assert_allclose(lon_b, lons, atol=1e-10)
    np.testing.assert_allclose(h_b, hs, atol=1e-6)


def test_radians_input():
    lat0, lon0, h0 = np.radians(60.0), np.radians(10.0), 0.0
    lat, lon, h = np.radians(60.001), np.radians(10.001), 50.0
    e, n, u = geodetic2enu(lat, lon, h, lat0, lon0, h0, radians=True)
    assert e == pytest.approx(55.7987538841, abs=1e-6)
    assert n == pytest.approx(111.4135903132, abs=1e-6)
    assert u == pytest.approx(49.9987842681, abs=1e-6)


def test_custom_ellipsoid_grs80():
    lat0, lon0, h0 = 60.0, 10.0, 0.0
    lat, lon, h = 60.5, 10.5, 100.0
    e1, n1, u1 = geodetic2enu(lat, lon, h, lat0, lon0, h0, ellipsoid=WGS84())
    e2, n2, u2 = geodetic2enu(lat, lon, h, lat0, lon0, h0, ellipsoid=GRS80())
    # WGS84 and GRS80 differ only at sub-millimetre scale for this geometry.
    assert e1 == pytest.approx(e2, abs=1e-3)
    assert n1 == pytest.approx(n2, abs=1e-3)
    assert u1 == pytest.approx(u2, abs=1e-3)


def test_invalid_ellipsoid():
    with pytest.raises(TypeError):
        geodetic2enu(60.1, 10.1, 0.0, 60.0, 10.0, 0.0, ellipsoid="WGS84")
