"""
Tests for AER (azimuth / elevation / slant-range) conversions.
"""

import numpy as np
import pytest

from pygeodetics import (
    geodetic2aer,
    aer2geodetic,
    enu2aer,
    aer2enu,
    ned2aer,
    aer2ned,
    ecef2aer,
    aer2ecef,
    geodetic2enu,
)


# (lat0, lon0, h0, lat, lon, h)
CASES = [
    (60.0, 10.75, 100.0, 60.1, 10.9, 250.0),
    (-33.8688, 151.2093, 50.0, -34.0, 151.3, 200.0),
    (45.0, -75.0, 200.0, 46.0, -74.0, 1000.0),
    (-22.9068, -43.1729, 5.0, -22.95, -43.20, 800.0),
    (78.2232, 15.6267, 12.0, 78.5, 16.0, 100.0),
]

AER_REF = [
    (36.8094778709, 0.5545734359, 13928.650260),
    (150.0736490298, 0.4360421331, 16797.014385),
    (34.7603997832, -0.2732836671, 135881.661787),
    (210.1536471243, 8.1509095097, 5590.251883),
    (15.0409992898, 0.0139933000, 32028.326710),
]


@pytest.mark.parametrize("case, ref", list(zip(CASES, AER_REF)))
def test_geodetic2aer_reference(case, ref):
    lat0, lon0, h0, lat, lon, h = case
    az, el, rng = geodetic2aer(lat, lon, h, lat0, lon0, h0)
    # Use unit-circle comparison for azimuth to avoid 0/360 wrap issues.
    assert np.cos(np.radians(az)) == pytest.approx(np.cos(np.radians(ref[0])), abs=1e-9)
    assert np.sin(np.radians(az)) == pytest.approx(np.sin(np.radians(ref[0])), abs=1e-9)
    assert el == pytest.approx(ref[1], abs=1e-9)
    assert rng == pytest.approx(ref[2], abs=1e-3)


@pytest.mark.parametrize("case", CASES)
def test_geodetic_aer_round_trip(case):
    lat0, lon0, h0, lat, lon, h = case
    az, el, rng = geodetic2aer(lat, lon, h, lat0, lon0, h0)
    lat_b, lon_b, h_b = aer2geodetic(az, el, rng, lat0, lon0, h0)
    assert lat_b == pytest.approx(lat, abs=1e-9)
    assert lon_b == pytest.approx(lon, abs=1e-9)
    assert h_b == pytest.approx(h, abs=1e-6)


def test_azimuth_in_zero_to_360():
    # Point straight west of an equatorial observer -> azimuth 270.
    lat0, lon0, h0 = 0.0, 0.0, 0.0
    az, el, _ = geodetic2aer(0.0, -0.001, 0.0, lat0, lon0, h0)
    assert az == pytest.approx(270.0, abs=1e-6)
    assert 0.0 <= az < 360.0
    assert el == pytest.approx(0.0, abs=1e-3)


def test_zero_range():
    az, el, rng = geodetic2aer(60.0, 10.0, 100.0, 60.0, 10.0, 100.0)
    assert rng == pytest.approx(0.0, abs=1e-9)


def test_enu_aer_round_trip():
    e, n, u = 1234.5, -567.8, 90.1
    az, el, rng = enu2aer(e, n, u)
    e2, n2, u2 = aer2enu(az, el, rng)
    assert e2 == pytest.approx(e, abs=1e-9)
    assert n2 == pytest.approx(n, abs=1e-9)
    assert u2 == pytest.approx(u, abs=1e-9)


def test_ned_aer_round_trip():
    n, e, d = 100.0, 200.0, -50.0
    az, el, rng = ned2aer(n, e, d)
    n2, e2, d2 = aer2ned(az, el, rng)
    assert n2 == pytest.approx(n, abs=1e-9)
    assert e2 == pytest.approx(e, abs=1e-9)
    assert d2 == pytest.approx(d, abs=1e-9)


def test_ecef_aer_round_trip():
    lat0, lon0, h0 = 60.0, 10.0, 0.0
    X, Y, Z = 3100000.0, 700000.0, 5500000.0
    az, el, rng = ecef2aer(X, Y, Z, lat0, lon0, h0)
    X2, Y2, Z2 = aer2ecef(az, el, rng, lat0, lon0, h0)
    assert X2 == pytest.approx(X, abs=1e-3)
    assert Y2 == pytest.approx(Y, abs=1e-3)
    assert Z2 == pytest.approx(Z, abs=1e-3)


def test_aer_array_inputs():
    lat0, lon0, h0 = 60.0, 10.0, 0.0
    lats = np.array([60.1, 60.2, 60.3])
    lons = np.array([10.1, 10.2, 10.3])
    hs = np.array([100.0, 200.0, 300.0])
    az, el, rng = geodetic2aer(lats, lons, hs, lat0, lon0, h0)
    assert az.shape == lats.shape
    az_ref = np.array([26.525525863413968, 26.447334545262503, 26.369078618330274])
    el_ref = np.array([0.4040581842290203, 0.3483394827961724, 0.29265435452716515])
    rng_ref = np.array([12457.272878108812, 24907.327069485927, 37350.13780203923])
    np.testing.assert_allclose(np.cos(np.radians(az)), np.cos(np.radians(az_ref)), atol=1e-9)
    np.testing.assert_allclose(np.sin(np.radians(az)), np.sin(np.radians(az_ref)), atol=1e-9)
    np.testing.assert_allclose(el, el_ref, atol=1e-9)
    np.testing.assert_allclose(rng, rng_ref, atol=1e-3)


def test_enu2aer_reference():
    az, el, rng = enu2aer(1500.0, -700.0, 50.0)
    assert az == pytest.approx(115.0168934781, abs=1e-9)
    assert el == pytest.approx(1.7301562377, abs=1e-9)
    assert rng == pytest.approx(1656.049516, abs=1e-3)


def test_geodetic2aer_vs_geodetic2enu_consistency():
    lat0, lon0, h0 = 50.0, 8.0, 0.0
    lat, lon, h = 50.1, 8.2, 300.0
    az, el, rng = geodetic2aer(lat, lon, h, lat0, lon0, h0)
    e, n, u = geodetic2enu(lat, lon, h, lat0, lon0, h0)
    az2, el2, rng2 = enu2aer(e, n, u)
    assert az == pytest.approx(az2, abs=1e-9)
    assert el == pytest.approx(el2, abs=1e-9)
    assert rng == pytest.approx(rng2, abs=1e-9)
