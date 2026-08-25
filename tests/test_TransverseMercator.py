"""
Tests for the TransverseMercator class.

The forward and reverse series are Kruger's expansions in the third flattening ``n``,
carried to ``n**6`` (Karney 2011). Since ``n`` is about ``1/599``, the first neglected
term is of order ``1e-19``, so the projection is exact to double precision and the
tolerances below are at the nanometre level rather than the millimetre level.

All reference values are taken from pyproj/PROJ. pyproj is intentionally not a test
dependency, so the reference values are hard-coded.
"""

import numpy as np
import pytest
from pygeodetics.projections.TransverseMercator import TransverseMercator

WGS84_A = 6378137.0
WGS84_F = 1 / 298.257223563

# (lon, lat, easting, northing) for WGS84 / UTM zone 32N (EPSG:32632).
# The last three points are far outside the 6 degree zone, where a series truncated at
# n**4 is no longer exact.
UTM32N_CASES = [
    (3.0, 60.0, 165640.3321077152, 6666593.572146888),
    (10.0, 59.0, 557450.9412736214, 6540481.778558198),
    (14.0, 70.0, 690670.6937077054, 7773697.111538875),
    (6.0, 45.0, 263553.9738987921, 4987329.504698914),
    (3.5, 80.0, 393532.3345924397, 8886622.364645604),
    (9.0, 0.0, 500000.0, 0.0),
    (16.0, 20.0, 1233644.4158883446, 2226862.8936689002),
    (25.0, 10.0, 2275474.0483386135, 1149317.9780958903),
]

# WGS84 / UTM zone 32S (EPSG:32732), which needs the 10 000 000 m false northing.
UTM32S_CASES = [
    (9.5, -30.0, 548224.1512265273, 6681109.436954567),
    (11.0, -45.0, 657630.6407299498, 5015103.828728243),
    (12.0, -70.0, 614473.7147457979, 2231309.891224767),
]


@pytest.fixture
def utm32n():
    return TransverseMercator(0.0, np.radians(9.0), 0.9996, 500000.0, 0.0, WGS84_A, WGS84_F)


@pytest.fixture
def utm32s():
    return TransverseMercator(0.0, np.radians(9.0), 0.9996, 500000.0, 10000000.0,
                              WGS84_A, WGS84_F)


@pytest.mark.parametrize("lon, lat, east, north", UTM32N_CASES)
def test_geog_to_projected_utm32n(utm32n, lon, lat, east, north):
    E, N = utm32n.geog_to_projected(np.array([[lon, lat]]))
    assert E[0] == pytest.approx(east, abs=1e-6)
    assert N[0] == pytest.approx(north, abs=1e-6)


@pytest.mark.parametrize("lon, lat, east, north", UTM32S_CASES)
def test_geog_to_projected_utm32s(utm32s, lon, lat, east, north):
    E, N = utm32s.geog_to_projected(np.array([[lon, lat]]))
    assert E[0] == pytest.approx(east, abs=1e-6)
    assert N[0] == pytest.approx(north, abs=1e-6)


@pytest.mark.parametrize("lon, lat, east, north", UTM32N_CASES)
def test_projected_to_geog_utm32n(utm32n, lon, lat, east, north):
    result = utm32n.projected_to_geog(np.array([[east, north]]))
    assert result[0][0] == pytest.approx(lon, abs=1e-11)
    assert result[1][0] == pytest.approx(lat, abs=1e-11)


def test_round_trip_is_exact(utm32n):
    """Forward followed by reverse must return the original coordinates."""
    coords = np.array([[lon, lat] for lon, lat, _, _ in UTM32N_CASES])
    projected = np.array(utm32n.geog_to_projected(coords)).T
    back = np.array(utm32n.projected_to_geog(projected)).T

    np.testing.assert_allclose(back, coords, rtol=0, atol=1e-11)


def test_far_field_beats_the_n4_truncation(utm32n):
    """
    At 16 degrees from the central meridian a series truncated at n**4 is off by roughly
    a millimetre. The n**6 series must do far better than that.
    """
    lon, lat, east, north = 25.0, 10.0, 2275474.0483386135, 1149317.9780958903
    E, N = utm32n.geog_to_projected(np.array([[lon, lat]]))

    assert abs(E[0] - east) < 1e-7
    assert abs(N[0] - north) < 1e-7


def test_vectorised_matches_pointwise(utm32n):
    coords = np.array([[lon, lat] for lon, lat, _, _ in UTM32N_CASES])
    expected = np.array([[e, n] for _, _, e, n in UTM32N_CASES])

    projected = np.array(utm32n.geog_to_projected(coords)).T

    np.testing.assert_allclose(projected, expected, rtol=0, atol=1e-6)


def test_height_is_passed_through(utm32n):
    coords = np.array([[10.0, 59.0, 123.456]])
    result = utm32n.geog_to_projected(coords)
    assert result.shape[0] == 3
    assert result[2][0] == pytest.approx(123.456)


def test_radian_input(utm32n):
    coords = np.array([[np.radians(10.0), np.radians(59.0)]])
    E, N = utm32n.geog_to_projected(coords, unit="rad")
    assert E[0] == pytest.approx(557450.9412736214, abs=1e-6)
    assert N[0] == pytest.approx(6540481.778558198, abs=1e-6)


def test_non_zero_latitude_of_origin():
    """
    OSGB36 / British National Grid (EPSG:27700) exercises a non-zero latitude of natural
    origin together with a negative false northing.
    """
    airy_a = 6377563.396
    airy_f = 1 / 299.3249646
    bng = TransverseMercator(np.radians(49.0), np.radians(-2.0), 0.9996012717,
                             400000.0, -100000.0, airy_a, airy_f)

    E, N = bng.geog_to_projected(np.array([[0.5, 53.0]]))
    assert E[0] == pytest.approx(567745.131459005, abs=1e-6)
    assert N[0] == pytest.approx(347699.8539085306, abs=1e-6)

    back = bng.projected_to_geog(np.array([[E[0], N[0]]]))
    assert back[0][0] == pytest.approx(0.5, abs=1e-11)
    assert back[1][0] == pytest.approx(53.0, abs=1e-11)


def test_invalid_shape_raises(utm32n):
    with pytest.raises(ValueError):
        utm32n.geog_to_projected(np.array([[1.0]]))
    with pytest.raises(ValueError):
        utm32n.projected_to_geog(np.array([[1.0]]))


def test_invalid_type_raises(utm32n):
    with pytest.raises(TypeError):
        utm32n.geog_to_projected("not coordinates")
