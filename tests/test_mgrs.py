"""
Tests for the MGRS encoder / parser and EPSG helpers.

Reference values are cross-validated against the third-party ``mgrs``
package where it is available; otherwise hard-coded reference strings
from authoritative sources are used.
"""

import numpy as np
import pytest

from pygeodetics import (
    to_mgrs,
    from_mgrs,
    geodetic_to_utm,
    utm_to_geodetic,
    geodetic_to_ups,
    ups_to_geodetic,
    utm_epsg,
    ups_epsg,
)



REFERENCE_POINTS = [
    # UTM, northern hemisphere
    (60.0, 10.75, "32VNM9760352702"),
    # UTM, southern hemisphere
    (-33.8688, 151.2093, "56HLH3436850948"),
    # Equator / prime-meridian neighbourhood
    (0.0, 0.0, "31NAA6602100000"),
    # North America
    (45.0, -75.0, "18TWQ0000082950"),
    # UPS north (zone Z, lon >= 0)
    (85.0, 0.0, "ZAB0000044542"),
    # UPS north (zone Y, lon < 0)
    (85.0, -45.0, "YUD0723207232"),
    # UPS south (zone B, lon >= 0)
    (-85.0, 0.0, "BAT0000055457"),
    (-85.0, 45.0, "BFR9276792767"),
    # Exact poles (mgrs convention sends them to zone Z / B respectively)
    (90.0, 0.0, "ZAH0000000000"),
    (-90.0, 0.0, "BAN0000000000"),
]


@pytest.mark.parametrize("lat, lon, expected", REFERENCE_POINTS)
def test_to_mgrs_matches_reference(lat, lon, expected):
    assert to_mgrs(lat, lon, precision=5) == expected


@pytest.mark.parametrize("lat, lon, mgrs_str", REFERENCE_POINTS)
def test_from_mgrs_recovers_position(lat, lon, mgrs_str):
    """
    Parsing returns the SW corner of the 1 m MGRS cell. Recovered position
    must be within ~1 m of the original.
    """
    lat_back, lon_back = from_mgrs(mgrs_str)
    # 1e-4 deg ≈ 11 m, but we expect <2 m -> use 2e-5 deg ≈ 2.2 m.
    assert lat_back == pytest.approx(lat, abs=2e-5)
    # Longitude tolerance must scale with cos(lat); use a generous bound
    # since at high latitudes 1 m → many degrees.
    cos_lat = max(np.cos(np.radians(lat)), 0.05)
    assert lon_back == pytest.approx(lon, abs=2e-5 / cos_lat)



@pytest.mark.parametrize("lat, lon", [
    (60.0, 10.75),
    (-33.8688, 151.2093),
    (0.0, 0.0),
    (45.0, -75.0),
    (-22.9, -43.2),
    (-80.0 + 0.001, 0.0),
    (84.0 - 0.001, 0.0),
])
def test_round_trip_utm_precision5(lat, lon):
    mgrs_str = to_mgrs(lat, lon, precision=5)
    lat_back, lon_back = from_mgrs(mgrs_str)
    # SW-corner of the 1 m cell -> position diff is bounded by sqrt(2) m.
    assert lat_back == pytest.approx(lat, abs=2e-5)
    cos_lat = max(np.cos(np.radians(lat)), 0.05)
    assert lon_back == pytest.approx(lon, abs=2e-5 / cos_lat)



def test_precision_levels_format_and_size():
    lat, lon = 60.0, 10.75
    expected_lengths = {
        0: 3,   # GZD + 100km square only
        1: 5,
        2: 7,
        3: 9,
        4: 11,
        5: 13,
    }
    for p, length in expected_lengths.items():
        s = to_mgrs(lat, lon, precision=p)
        # 32V (3) + NM (2) + 2*p digits = 5 + 2p
        assert len(s) == 5 + 2 * p, f"precision={p} -> '{s}'"


@pytest.mark.parametrize("bad", [-1, 6, 10])
def test_invalid_precision(bad):
    with pytest.raises(ValueError):
        to_mgrs(60.0, 10.75, precision=bad)


def test_precision_3_string_known():
    # 100m precision => first 3 digits of east 97603 -> 976, of north 52702 -> 527
    s = to_mgrs(60.0, 10.75, precision=3)
    assert s == "32VNM976527"



@pytest.mark.parametrize("bad_str", [
    "",                       # empty
    "32V",                    # missing 100km square
    "32VNM1",                 # odd-length numeric part
    "32VNM123",               # odd-length numeric part
    "61VNM9760352702",        # zone out of range
    "32VIM9760352702",        # invalid 'I' in 100km square
    "32VOM9760352702",        # invalid 'O'
    "32IM",                   # invalid band 'I'
    "32YNM",                  # zone+band but Y is north UPS letter
    "garbage",
])
def test_from_mgrs_rejects_malformed(bad_str):
    with pytest.raises(ValueError):
        from_mgrs(bad_str)


def test_from_mgrs_accepts_spaces_and_lowercase():
    a = from_mgrs("32VNM9760352702")
    b = from_mgrs("32V NM 97603 52702")
    c = from_mgrs("32vnm9760352702")
    np.testing.assert_allclose(a, b, atol=1e-12)
    np.testing.assert_allclose(a, c, atol=1e-12)


def test_from_mgrs_type_check():
    with pytest.raises(TypeError):
        from_mgrs(12345)



mgrs_pkg = pytest.importorskip("mgrs")


@pytest.mark.parametrize("lat, lon", [
    (60.0, 10.75),
    (-33.8688, 151.2093),
    (0.0, 0.0),
    (45.0, -75.0),
    (-22.9, -43.2),
    (40.7128, -74.0060),    # New York
    (35.6762, 139.6503),    # Tokyo
    (-77.8419, 166.6863),   # McMurdo (just inside UTM band C edge)
    (85.0, 0.0),
    (85.0, -45.0),
    (-85.0, 0.0),
    (-85.0, 45.0),
])
def test_cross_check_against_mgrs_package(lat, lon):
    ref = mgrs_pkg.MGRS().toMGRS(lat, lon, MGRSPrecision=5)
    ours = to_mgrs(lat, lon, precision=5)
    assert ours == ref



@pytest.mark.parametrize("zone, hemi, expected", [
    (32, "N", 32632),
    (33, "N", 32633),
    (56, "S", 32756),
    (1,  "N", 32601),
    (60, "S", 32760),
])
def test_utm_epsg_wgs84(zone, hemi, expected):
    assert utm_epsg(zone, hemi) == expected


def test_utm_epsg_nad83():
    assert utm_epsg(15, "N", datum="NAD83") == 26915
    with pytest.raises(ValueError):
        utm_epsg(15, "S", datum="NAD83")


@pytest.mark.parametrize("bad_zone", [0, 61, -1])
def test_utm_epsg_validation(bad_zone):
    with pytest.raises(ValueError):
        utm_epsg(bad_zone, "N")


def test_ups_epsg():
    assert ups_epsg("N") == 32661
    assert ups_epsg("S") == 32761
    with pytest.raises(ValueError):
        ups_epsg("X")
