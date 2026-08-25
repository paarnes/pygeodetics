"""
Tests for the Transverse Mercator grid convergence functions.

All reference values are taken from pyproj/PROJ (`Proj.get_factors(...).meridian_convergence`),
which uses a rigorous (non-truncated) formulation. pyproj is intentionally not a test
dependency, so the reference values are hard-coded.
"""

import doctest

import numpy as np
import pytest
from pygeodetics.Ellipsoid import WGS84, GRS80
from pygeodetics.projections import GridConvergence
from pygeodetics.projections.GridConvergence import (
    tm_grid_convergence_projected,
    tm_grid_convergence_geographic,
)

# Airy 1830, used by the British National Grid (EPSG:27700)
AIRY1830_A = 6377563.396
AIRY1830_B = 6356256.909237285


### TEST 1 - grid convergence from projected coordinates

test_cases_tm_proj = [
    {
        "a": WGS84().a,
        "b": WGS84().b,
        "x": 555776.2667516097,
        "y": 6651832.735433666,
        "false_easting": 500000.0,
        "false_northing": 0.0,
        "k0": 0.9996,
        "lat0": 0.0,
        "radians": False,
        "description": "WGS84 / UTM zone 32N (EPSG:32632), 10E 60N",
        "gamma_true": 0.866047498561436,
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
        "radians": False,
        "description": "ETRS89 / UTM zone 32N (EPSG:25832), 10E 59N",
        "gamma_true": 0.8571905119316398,
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
        "radians": False,
        "description": "ETRS89 / UTM zone 32N, west of the central meridian (issue #14)",
        "gamma_true": -3.3018337424508464,
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
        "radians": False,
        "description": "WGS84 / UTM zone 32S (EPSG:32732), 9.5E 30S - needs false_northing",
        "gamma_true": -0.2500048321934752,
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
        "radians": False,
        "description": "Pure Transverse Mercator on GRS80, lon0 = 9E, 10E 60N",
        "gamma_true": 0.8660474985775397,
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
        "radians": False,
        "description": "OSGB36 / British National Grid (EPSG:27700), 0.5E 53N",
        "gamma_true": 1.9970510531036239,
    },
]


@pytest.mark.parametrize("method", ["kruger", "snyder"])
@pytest.mark.parametrize("case", test_cases_tm_proj, ids=lambda c: c["description"])
def test_tm_grid_convergence_proj(case, method):
    """Grid convergence from projected coordinates against rigorous reference values."""
    gamma = tm_grid_convergence_projected(
        case["a"], case["b"], case["x"], case["y"],
        lat0=case["lat0"],
        false_easting=case["false_easting"],
        false_northing=case["false_northing"],
        k0=case["k0"],
        radians=case["radians"],
        method=method,
    )

    assert np.isclose(gamma, case["gamma_true"], rtol=0, atol=1e-6), (
        f"Test failed for case: {case['description']} (method={method})\n"
        f"Computed Grid Convergence: {gamma} degrees\n"
        f"Expected Grid Convergence: {case['gamma_true']} degrees"
    )


def test_tm_grid_convergence_proj_radians_matches_degrees():
    """The `radians` flag must only change the units, not the value."""
    a, b = GRS80().a, GRS80().b
    deg = tm_grid_convergence_projected(a, b, 557450.9412743163, 6540481.778434776,
                                        lat0=0.0, false_easting=500000.0, k0=0.9996)
    rad = tm_grid_convergence_projected(a, b, 557450.9412743163, 6540481.778434776,
                                        lat0=0.0, false_easting=500000.0, k0=0.9996,
                                        radians=True)
    assert np.isclose(np.degrees(rad), deg, rtol=0, atol=1e-12)


def test_tm_grid_convergence_proj_does_not_mutate_input():
    """Issue #16: the caller's NumPy arrays must not be modified in place."""
    a, b = GRS80().a, GRS80().b
    x = np.array([600000.0, 700000.0])
    y = np.array([6.6e6, 6.6e6])
    x_copy, y_copy = x.copy(), y.copy()

    first = tm_grid_convergence_projected(a, b, x, y, lat0=0.0, false_easting=500000.0)

    assert np.array_equal(x, x_copy), "x was modified in place"
    assert np.array_equal(y, y_copy), "y was modified in place"

    second = tm_grid_convergence_projected(a, b, x, y, lat0=0.0, false_easting=500000.0)
    assert np.array_equal(first, second), "repeated calls gave different results"


def test_tm_grid_convergence_proj_vectorised():
    """The projected variant must vectorise over NumPy arrays."""
    a, b = GRS80().a, GRS80().b
    x = np.array([557450.9412743163, 280841.3814])
    y = np.array([6540481.778434776, 6571720.2565])
    expected = np.array([0.8571905119316398, -3.3018337424508464])

    gamma = tm_grid_convergence_projected(a, b, x, y, lat0=0.0,
                                          false_easting=500000.0, k0=0.9996)

    assert isinstance(gamma, np.ndarray)
    assert np.allclose(gamma, expected, rtol=0, atol=1e-6)


def test_tm_grid_convergence_proj_requires_keyword_arguments():
    """Issue #22: lat0/false_easting/k0 are keyword-only."""
    a, b = GRS80().a, GRS80().b
    with pytest.raises(TypeError):
        tm_grid_convergence_projected(a, b, 557450.9, 6540481.8, 0.0, 500000.0)


def test_tm_grid_convergence_proj_defaults_are_pure_tm():
    """Defaults must correspond to pure TM: k0 = 1, no false easting/northing, lat0 = 0."""
    a, b = GRS80().a, GRS80().b
    x, y = 55798.58618677315, 6654494.533123902
    assert tm_grid_convergence_projected(a, b, x, y) == tm_grid_convergence_projected(
        a, b, x, y, lat0=0.0, false_easting=0.0, false_northing=0.0, k0=1.0
    )


### TEST 2 - grid convergence from geographic coordinates

test_cases_geog = [
    {
        "a": WGS84().a,
        "b": WGS84().b,
        "lon": np.deg2rad(10),
        "lat": np.deg2rad(60),
        "central_meridian": np.deg2rad(9),
        "radians": True,
        "description": "WGS84, lon0 = 9E, 10E 60N",
        "gamma_true": 0.8660474985710462,
    },
    {
        "a": GRS80().a,
        "b": GRS80().b,
        "lon": np.deg2rad(10),
        "lat": np.deg2rad(59),
        "central_meridian": np.deg2rad(9),
        "radians": True,
        "description": "GRS80, lon0 = 9E, 10E 59N",
        "gamma_true": 0.8571905119311873,
    },
    {
        "a": GRS80().a,
        "b": GRS80().b,
        "lon": np.deg2rad(17),
        "lat": np.deg2rad(59),
        "central_meridian": np.deg2rad(9),
        "radians": True,
        "description": "GRS80, lon0 = 9E, 17E 59N (8 degrees from the central meridian)",
        "gamma_true": 6.869212584047117,
    },
]


@pytest.mark.parametrize("method", ["kruger", "snyder"])
@pytest.mark.parametrize("case", test_cases_geog, ids=lambda c: c["description"])
def test_tm_grid_convergence_geog(case, method):
    """Grid convergence from geographic coordinates against rigorous reference values."""
    gamma = np.degrees(
        tm_grid_convergence_geographic(
            case["a"], case["b"], case["lon"], case["lat"], case["central_meridian"],
            radians=case["radians"], method=method,
        )
    )

    assert np.isclose(gamma, case["gamma_true"], rtol=0, atol=1e-6), (
        f"Test failed for case: {case['description']} (method={method})\n"
        f"Computed Grid Convergence: {gamma} degrees\n"
        f"Expected Grid Convergence: {case['gamma_true']} degrees"
    )


def test_tm_grid_convergence_geog_degrees_and_radians_agree():
    a, b = GRS80().a, GRS80().b
    deg = tm_grid_convergence_geographic(a, b, 10.0, 59.0, 9.0)
    rad = tm_grid_convergence_geographic(a, b, np.deg2rad(10.0), np.deg2rad(59.0),
                                         np.deg2rad(9.0), radians=True)
    assert np.isclose(np.degrees(rad), deg, rtol=0, atol=1e-12)


def test_tm_grid_convergence_geog_vectorised():
    a, b = GRS80().a, GRS80().b
    lon = np.array([10.0, 17.0])
    lat = np.array([59.0, 59.0])
    gamma = tm_grid_convergence_geographic(a, b, lon, lat, 9.0)

    assert isinstance(gamma, np.ndarray)
    assert np.allclose(gamma, [0.8571905119311873, 6.869212584047117], rtol=0, atol=1e-6)


def test_tm_grid_convergence_geog_is_independent_of_k0():
    """
    Issue #17: grid convergence does not depend on k0. The geographic variant (which has
    no k0 argument) must agree with the projected variant evaluated on a k0-scaled grid.
    """
    a, b = GRS80().a, GRS80().b
    geographic = tm_grid_convergence_geographic(a, b, 10.0, 59.0, 9.0)
    projected = tm_grid_convergence_projected(a, b, 557450.9412743163, 6540481.778434776,
                                              lat0=0.0, false_easting=500000.0, k0=0.9996)
    assert np.isclose(geographic, projected, rtol=0, atol=1e-6)


def test_tm_grid_convergence_geog_prime_meridian_frames_agree():
    """
    Issue #19: `lon` and `lon0` must be in the same prime-meridian frame. NGO 1948
    zone I (EPSG:27391) has its central meridian 4 deg 40' west of Oslo.
    """
    a, b = GRS80().a, GRS80().b
    oslo_from_greenwich = 10.0 + 43 / 60 + 22.5 / 3600  # 10.7229166667 deg
    lon_greenwich, lat = 11.0, 60.0

    greenwich_frame = tm_grid_convergence_geographic(
        a, b, lon_greenwich, lat, oslo_from_greenwich - 4 - 40 / 60)
    oslo_frame = tm_grid_convergence_geographic(
        a, b, lon_greenwich - oslo_from_greenwich, lat, -4 - 40 / 60)

    assert np.isclose(greenwich_frame, oslo_frame, rtol=0, atol=1e-9)


def test_docstring_examples():
    """Issue #20: keep the documented examples from drifting away from the signatures."""
    results = doctest.testmod(GridConvergence, verbose=False)
    assert results.failed == 0
    assert results.attempted > 0


### TEST 3 - the Kruger series (default) versus the truncated Snyder series

# Rigorous reference values (pyproj), pure Transverse Mercator on GRS80 with lon0 = 9E,
# at lat = 59N, for increasing distance from the central meridian.
far_field_cases = [
    (1.0, 0.8571905119311873),
    (2.0, 1.7145202838493365),
    (4.0, 3.4301544266421655),
    (8.0, 6.869212584047117),
    (16.0, 13.809415701653856),
    (30.0, 26.33304621386647),
]


@pytest.mark.parametrize("dlon, gamma_true", far_field_cases)
def test_kruger_is_accurate_far_from_the_central_meridian(dlon, gamma_true):
    """
    The Kruger series is in the third flattening, not in `dlon`, so it stays exact to
    double precision arbitrarily far from the central meridian.
    """
    a, b = GRS80().a, GRS80().b
    gamma = tm_grid_convergence_geographic(a, b, 9.0 + dlon, 59.0, 9.0)
    assert np.isclose(gamma, gamma_true, rtol=0, atol=1e-9)


def test_snyder_degrades_far_from_the_central_meridian():
    """
    Documented limitation of `method="snyder"`: it is accurate close to the central
    meridian but loses accuracy quickly beyond a UTM zone width.
    """
    a, b = GRS80().a, GRS80().b
    near = tm_grid_convergence_geographic(a, b, 10.0, 59.0, 9.0, method="snyder")
    far = tm_grid_convergence_geographic(a, b, 39.0, 59.0, 9.0, method="snyder")

    assert np.isclose(near, 0.8571905119311873, rtol=0, atol=1e-9)
    assert not np.isclose(far, 26.33304621386647, rtol=0, atol=1e-3)


@pytest.mark.parametrize("func, args, kwargs", [
    (tm_grid_convergence_geographic, (GRS80().a, GRS80().b, 10.0, 59.0, 9.0), {}),
    (tm_grid_convergence_projected, (GRS80().a, GRS80().b, 557450.94, 6540481.78),
     {"false_easting": 500000.0, "k0": 0.9996}),
])
def test_unknown_method_raises(func, args, kwargs):
    with pytest.raises(ValueError, match="method must be"):
        func(*args, method="bogus", **kwargs)
