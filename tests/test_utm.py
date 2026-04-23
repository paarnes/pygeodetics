"""
Tests for the high-level UTM wrapper API
(:func:`pygeodetics.geodetic_to_utm` / :func:`pygeodetics.utm_to_geodetic`).
"""

import numpy as np
import pytest

from pygeodetics import (
    geodetic_to_utm,
    utm_to_geodetic,
    mgrs_band_letter,
)
from pygeodetics.Ellipsoid import WGS84, GRS80



@pytest.mark.parametrize(
    "lat, expected",
    [
        (-80.0, "C"),
        (-72.0, "D"),
        (0.0, "N"),
        (8.0, "P"),       # 'I' and 'O' skipped
        (60.0, "V"),
        (72.0, "X"),
        (83.9, "X"),
        (84.0, "X"),
    ],
)
def test_mgrs_band_letter_scalar(lat, expected):
    assert mgrs_band_letter(lat) == expected


def test_mgrs_band_letter_array():
    lats = np.array([-80.0, 0.0, 60.0, 84.0])
    bands = mgrs_band_letter(lats)
    assert list(bands) == ["C", "N", "V", "X"]


@pytest.mark.parametrize("lat", [-80.1, 84.1, -90.0, 90.0])
def test_mgrs_band_letter_out_of_range(lat):
    with pytest.raises(ValueError):
        mgrs_band_letter(lat)



ROUND_TRIP_POINTS = [
    # (lat, lon, expected_zone, expected_band)
    (60.0, 10.75, 32, "V"),    # Oslo, Norway (N hemisphere)
    (-33.8688, 151.2093, 56, "H"),  # Sydney, Australia (S hemisphere)
    (0.0, 0.0, 31, "N"),       # Equator / prime meridian
    (-80.0, -179.999, 1, "C"), # Far SW corner of UTM validity
    (84.0, 179.999, 60, "X"),  # Far NE corner of UTM validity
    (45.0, -75.0, 18, "T"),    # Ottawa-ish
    (-22.9, -43.2, 23, "K"),   # Rio de Janeiro
]


@pytest.mark.parametrize("lat, lon, expected_zone, expected_band",
                         ROUND_TRIP_POINTS)
def test_round_trip_sub_mm(lat, lon, expected_zone, expected_band):
    e, n, zone, band = geodetic_to_utm(lat, lon)
    assert zone == expected_zone
    assert band == expected_band

    hemi = "N" if lat >= 0 else "S"
    lat_back, lon_back = utm_to_geodetic(e, n, zone=zone, hemisphere=hemi)

    # Latitude / longitude tolerances: 1e-10 deg ~ 1.1e-5 m on Earth, so
    # round-trip error in metres is well below 1 mm.
    assert lat_back == pytest.approx(lat, abs=1e-10)
    assert lon_back == pytest.approx(lon, abs=1e-10)



def test_known_point_oslo_zone32_north():
    """
    Oslo (60°N, 10.75°E) on WGS84 in UTM zone 32N.

    Reference values were computed independently with the GeographicLib /
    pyproj pipeline to ~1e-3 m accuracy.
    """
    e, n, zone, band = geodetic_to_utm(60.0, 10.75)
    assert zone == 32
    assert band == "V"
    # Reference: pyproj EPSG:4326 -> EPSG:32632
    assert e == pytest.approx(597603.3590, abs=1e-3)
    assert n == pytest.approx(6652702.2062, abs=1e-3)


def test_known_point_sydney_zone56_south():
    """
    Sydney (-33.8688°, 151.2093°) on WGS84 in UTM zone 56H.

    Reference values from pyproj (EPSG:4326 -> EPSG:32756).
    """
    e, n, zone, band = geodetic_to_utm(-33.8688, 151.2093)
    assert zone == 56
    assert band == "H"
    assert e == pytest.approx(334368.6336, abs=1e-3)
    assert n == pytest.approx(6250948.3454, abs=1e-3)
    # Southern hemisphere => false northing 10 000 000 must be applied.
    assert 0.0 < n < 10_000_000.0



def test_southern_hemisphere_false_northing_applied():
    """A point just south of the equator must yield N close to 1e7."""
    _, n, _, _ = geodetic_to_utm(-0.0001, 10.0)
    assert n == pytest.approx(10_000_000.0, abs=50.0)
    assert n < 10_000_000.0


def test_northern_hemisphere_false_northing_zero_at_equator():
    _, n, _, _ = geodetic_to_utm(0.0001, 10.0)
    assert n == pytest.approx(0.0, abs=50.0)
    assert n > 0.0



def test_force_zone_overrides_auto_zone():
    # 10.75°E would normally fall in zone 32; force zone 33 instead.
    e_auto, _, zone_auto, _ = geodetic_to_utm(60.0, 10.75)
    e_forced, _, zone_forced, _ = geodetic_to_utm(60.0, 10.75, force_zone=33)
    assert zone_auto == 32
    assert zone_forced == 33
    assert not np.isclose(e_auto, e_forced)

    # Round trip with the forced zone must still recover the input.
    lat_back, lon_back = utm_to_geodetic(
        e_forced,
        geodetic_to_utm(60.0, 10.75, force_zone=33)[1],
        zone=33,
        hemisphere="N",
    )
    assert lat_back == pytest.approx(60.0, abs=1e-10)
    assert lon_back == pytest.approx(10.75, abs=1e-10)


@pytest.mark.parametrize("bad_zone", [0, 61, -1, 100])
def test_force_zone_validation(bad_zone):
    with pytest.raises(ValueError):
        geodetic_to_utm(60.0, 10.0, force_zone=bad_zone)


@pytest.mark.parametrize("bad_zone", [0, 61, -3])
def test_inverse_zone_validation(bad_zone):
    with pytest.raises(ValueError):
        utm_to_geodetic(500000.0, 6_651_000.0, zone=bad_zone, hemisphere="N")


def test_inverse_hemisphere_validation():
    with pytest.raises(ValueError):
        utm_to_geodetic(500000.0, 6_651_000.0, zone=32, hemisphere="X")



def test_array_inputs_round_trip():
    lats = np.array([60.0, -33.8688, 45.0, -22.9])
    lons = np.array([10.75, 151.2093, -75.0, -43.2])
    e, n, zones, bands = geodetic_to_utm(lats, lons)

    assert e.shape == lats.shape
    assert n.shape == lats.shape
    assert zones.shape == lats.shape
    assert bands.shape == lats.shape
    assert list(zones) == [32, 56, 18, 23]
    assert list(bands) == ["V", "H", "T", "K"]

    # Round trip each point through its own (zone, hemisphere).
    for i in range(lats.size):
        hemi = "N" if lats[i] >= 0 else "S"
        lat_back, lon_back = utm_to_geodetic(
            e[i], n[i], zone=int(zones[i]), hemisphere=hemi
        )
        assert lat_back == pytest.approx(lats[i], abs=1e-10)
        assert lon_back == pytest.approx(lons[i], abs=1e-10)


def test_array_inputs_inverse_single_zone():
    """Inverse should accept arrays for a single zone/hemisphere."""
    lats = np.array([59.0, 60.0, 61.0])
    lons = np.array([10.0, 10.5, 11.0])
    e, n, zones, _ = geodetic_to_utm(lats, lons)
    assert np.all(zones == 32)

    lat_back, lon_back = utm_to_geodetic(e, n, zone=32, hemisphere="N")
    np.testing.assert_allclose(lat_back, lats, atol=1e-10)
    np.testing.assert_allclose(lon_back, lons, atol=1e-10)



def test_custom_ellipsoid_grs80():
    """GRS80 should give nearly identical results to WGS84 (a equal,
    f differs at ~1e-9 level)."""
    lat, lon = 60.0, 10.75
    e_wgs, n_wgs, _, _ = geodetic_to_utm(lat, lon, ellipsoid=WGS84())
    e_grs, n_grs, _, _ = geodetic_to_utm(lat, lon, ellipsoid=GRS80())
    assert e_wgs == pytest.approx(e_grs, abs=1e-3)
    assert n_wgs == pytest.approx(n_grs, abs=1e-3)


def test_invalid_ellipsoid_type():
    with pytest.raises(TypeError):
        geodetic_to_utm(60.0, 10.0, ellipsoid="WGS84")
