"""
Tests for the Transverse Mercator point scale factor functions.

All reference values are taken from pyproj/PROJ (`Proj.get_factors(...).meridional_scale`),
which uses a rigorous (non-truncated) formulation. pyproj is intentionally not a test
dependency, so the reference values are hard-coded.
"""

import doctest

import numpy as np
import pytest
from pygeodetics.Ellipsoid import WGS84, GRS80
from pygeodetics.projections import ScaleFactor
from pygeodetics.projections.ScaleFactor import (
    tm_point_scale_factor_geographic,
    tm_point_scale_factor_projected,
)

# Airy 1830, used by the British National Grid (EPSG:27700)
AIRY1830_A = 6377563.396
AIRY1830_B = 6356256.909237285


### TEST 1 - scale factor from geographic coordinates

# Issue #15: the 4th-order coefficient of the old series was wrong, which showed up as a
# 2.5e-05 error 8 degrees off the central meridian. Reference values from pyproj, GRS80,
# pure Transverse Mercator (k0 = 1), lon0 = 9E, lat = 59N.
test_cases_geog = [
    (10.0, 1.0000404725867253, "1 degree from the central meridian"),
    (11.0, 1.000161870124904, "2 degrees from the central meridian"),
    (13.0, 1.000647156435287, "4 degrees from the central meridian"),
    (17.0, 1.0025834089311234, "8 degrees from the central meridian"),
]


@pytest.mark.parametrize("method", ["kruger", "snyder"])
@pytest.mark.parametrize("lon, k_true, description", test_cases_geog,
                         ids=[c[2] for c in test_cases_geog])
def test_tm_point_scale_factor_geographic(lon, k_true, description, method):
    a, b = GRS80().a, GRS80().b
    k = tm_point_scale_factor_geographic(a, b, lon, 59.0, 9.0, method=method)

    assert np.isclose(k, k_true, rtol=0, atol=1e-9), (
        f"Test failed for case: {description} (method={method})\n"
        f"Computed scale factor: {k}\n"
        f"Expected scale factor: {k_true}"
    )


def test_tm_point_scale_factor_geographic_k0():
    """Issue #17: `k0` scales the returned point scale factor."""
    a, b = GRS80().a, GRS80().b
    pure = tm_point_scale_factor_geographic(a, b, 10.0, 59.0, 9.0)
    utm = tm_point_scale_factor_geographic(a, b, 10.0, 59.0, 9.0, k0=0.9996)

    assert np.isclose(utm, 0.9996 * pure, rtol=0, atol=0)
    # ETRS89 / UTM zone 32N at 10E 59N, pyproj reference
    assert np.isclose(utm, 0.9996404563973336, rtol=0, atol=1e-9)


def test_tm_point_scale_factor_geographic_degrees_and_radians_agree():
    a, b = GRS80().a, GRS80().b
    deg = tm_point_scale_factor_geographic(a, b, 10.0, 59.0, 9.0)
    rad = tm_point_scale_factor_geographic(a, b, np.deg2rad(10.0), np.deg2rad(59.0),
                                           np.deg2rad(9.0), radians=True)
    assert np.isclose(deg, rad, rtol=0, atol=1e-15)


def test_tm_point_scale_factor_geographic_vectorised():
    a, b = GRS80().a, GRS80().b
    lon = np.array([10.0, 17.0])
    lat = np.array([59.0, 59.0])
    k = tm_point_scale_factor_geographic(a, b, lon, lat, 9.0)

    assert isinstance(k, np.ndarray)
    assert np.allclose(k, [1.0000404725867253, 1.0025834089311234], rtol=0, atol=1e-9)


def test_tm_point_scale_factor_geographic_is_one_on_central_meridian():
    a, b = GRS80().a, GRS80().b
    assert np.isclose(tm_point_scale_factor_geographic(a, b, 9.0, 59.0, 9.0),
                      1.0, rtol=0, atol=1e-15)
    assert np.isclose(tm_point_scale_factor_geographic(a, b, 9.0, 59.0, 9.0, k0=0.9996),
                      0.9996, rtol=0, atol=1e-15)


### TEST 2 - scale factor from projected coordinates

test_cases_proj = [
    {
        "a": WGS84().a,
        "b": WGS84().b,
        "x": 555776.2667516097,
        "y": 6651832.735433666,
        "false_easting": 500000.0,
        "false_northing": 0.0,
        "k0": 0.9996,
        "lat0": 0.0,
        "description": "WGS84 / UTM zone 32N (EPSG:32632), 10E 60N",
        "k_true": 0.9996381243277405,
    },
    {
        "a": GRS80().a,
        "b": GRS80().b,
        "x": 557450.9412743163,
        "y": 6540481.778434776,
        "false_easting": 500000.0,
        "false_northing": 0.0,
        "k0": 0.9996,
        "lat0": 0.0,
        "description": "ETRS89 / UTM zone 32N (EPSG:25832), 10E 59N",
        "k_true": 0.9996404563973336,
    },
    {
        "a": GRS80().a,
        "b": GRS80().b,
        "x": 280841.3814,
        "y": 6571720.2565,
        "false_easting": 500000.0,
        "false_northing": 0.0,
        "k0": 0.9996,
        "lat0": 0.0,
        "description": "ETRS89 / UTM zone 32N, west of the central meridian (issue #14)",
        "k_true": 1.0001887417671618,
    },
    {
        "a": WGS84().a,
        "b": WGS84().b,
        "x": 548224.1512265273,
        "y": 6681109.436954567,
        "false_easting": 500000.0,
        "false_northing": 10000000.0,
        "k0": 0.9996,
        "lat0": 0.0,
        "description": "WGS84 / UTM zone 32S (EPSG:32732), 9.5E 30S - needs false_northing",
        "k_true": 0.9996286912772635,
    },
    {
        "a": GRS80().a,
        "b": GRS80().b,
        "x": 55798.58618677315,
        "y": 6654494.533123902,
        "false_easting": 0.0,
        "false_northing": 0.0,
        "k0": 1.0,
        "lat0": 0.0,
        "description": "Pure Transverse Mercator on GRS80, lon0 = 9E, 10E 60N",
        "k_true": 1.0000381395763323,
    },
    {
        "a": AIRY1830_A,
        "b": AIRY1830_B,
        "x": 567745.131459005,
        "y": 347699.8539085306,
        "false_easting": 400000.0,
        "false_northing": -100000.0,
        "k0": 0.9996012717,
        "lat0": 49.0,
        "description": "OSGB36 / British National Grid (EPSG:27700), 0.5E 53N",
        "k_true": 0.9999467012262234,
    },
]


@pytest.mark.parametrize("method", ["kruger", "snyder"])
@pytest.mark.parametrize("case", test_cases_proj, ids=lambda c: c["description"])
def test_tm_point_scale_factor_projected(case, method):
    k = tm_point_scale_factor_projected(
        case["a"], case["b"], case["x"], case["y"],
        lat0=case["lat0"],
        false_easting=case["false_easting"],
        false_northing=case["false_northing"],
        k0=case["k0"],
        method=method,
    )

    assert np.isclose(k, case["k_true"], rtol=0, atol=1e-9), (
        f"Test failed for case: {case['description']} (method={method})\n"
        f"Computed scale factor: {k}\n"
        f"Expected scale factor: {case['k_true']}"
    )


def test_tm_point_scale_factor_projected_does_not_mutate_input():
    """Issue #16: the caller's NumPy arrays must not be modified in place."""
    a, b = GRS80().a, GRS80().b
    x = np.array([600000.0, 700000.0])
    y = np.array([6.6e6, 6.6e6])
    x_copy, y_copy = x.copy(), y.copy()

    first = tm_point_scale_factor_projected(a, b, x, y, lat0=0.0, false_easting=500000.0)

    assert np.array_equal(x, x_copy), "x was modified in place"
    assert np.array_equal(y, y_copy), "y was modified in place"

    second = tm_point_scale_factor_projected(a, b, x, y, lat0=0.0, false_easting=500000.0)
    assert np.array_equal(first, second), "repeated calls gave different results"


def test_tm_point_scale_factor_projected_vectorised():
    a, b = GRS80().a, GRS80().b
    x = np.array([557450.9412743163, 280841.3814])
    y = np.array([6540481.778434776, 6571720.2565])

    k = tm_point_scale_factor_projected(a, b, x, y, lat0=0.0,
                                        false_easting=500000.0, k0=0.9996)

    assert isinstance(k, np.ndarray)
    assert np.allclose(k, [0.9996404563973336, 1.0001887417671618], rtol=0, atol=1e-9)


def test_tm_point_scale_factor_projected_requires_keyword_arguments():
    """Issue #22: lat0/false_easting/k0 are keyword-only."""
    a, b = GRS80().a, GRS80().b
    with pytest.raises(TypeError):
        tm_point_scale_factor_projected(a, b, 557450.9, 6540481.8, 0.0, 500000.0)


def test_tm_point_scale_factor_projected_matches_geographic():
    """The two formulations must agree for the same physical point."""
    a, b = GRS80().a, GRS80().b
    projected = tm_point_scale_factor_projected(a, b, 557450.9412743163, 6540481.778434776,
                                                lat0=0.0, false_easting=500000.0, k0=0.9996)
    geographic = tm_point_scale_factor_geographic(a, b, 10.0, 59.0, 9.0, k0=0.9996)
    assert np.isclose(projected, geographic, rtol=0, atol=1e-9)


def test_docstring_examples():
    """Issue #20: keep the documented examples from drifting away from the signatures."""
    results = doctest.testmod(ScaleFactor, verbose=False)
    assert results.failed == 0
    assert results.attempted > 0


### TEST 3 - the Kruger series (default) versus the truncated Snyder series

# Rigorous reference values (pyproj), pure Transverse Mercator on GRS80 with lon0 = 9E,
# at lat = 59N, far outside a UTM zone.
far_field_cases = [
    (16.0, 1.0102482198388723),
    (24.0, 1.0227228513216862),
    (30.0, 1.0349437565300619),
]


@pytest.mark.parametrize("dlon, k_true", far_field_cases)
def test_kruger_is_accurate_far_from_the_central_meridian(dlon, k_true):
    """
    The Kruger series is in the third flattening, not in `dlon`, so it stays exact to
    double precision arbitrarily far from the central meridian.
    """
    a, b = GRS80().a, GRS80().b
    k = tm_point_scale_factor_geographic(a, b, 9.0 + dlon, 59.0, 9.0)
    assert np.isclose(k, k_true, rtol=0, atol=1e-9)


def test_snyder_degrades_far_from_the_central_meridian():
    """
    Documented limitation of `method="snyder"`: it is accurate close to the central
    meridian but loses accuracy quickly beyond a UTM zone width.
    """
    a, b = GRS80().a, GRS80().b
    near = tm_point_scale_factor_geographic(a, b, 10.0, 59.0, 9.0, method="snyder")
    far = tm_point_scale_factor_geographic(a, b, 39.0, 59.0, 9.0, method="snyder")

    assert np.isclose(near, 1.0000404725867253, rtol=0, atol=1e-9)
    assert not np.isclose(far, 1.0349437565300619, rtol=0, atol=1e-6)


@pytest.mark.parametrize("func, args, kwargs", [
    (tm_point_scale_factor_geographic, (GRS80().a, GRS80().b, 10.0, 59.0, 9.0), {}),
    (tm_point_scale_factor_projected, (GRS80().a, GRS80().b, 557450.94, 6540481.78),
     {"false_easting": 500000.0, "k0": 0.9996}),
])
def test_unknown_method_raises(func, args, kwargs):
    with pytest.raises(ValueError, match="method must be"):
        func(*args, method="bogus", **kwargs)
