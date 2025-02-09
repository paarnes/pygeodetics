import numpy as np
import pytest
from Ellipsoid import WGS84, GRS80
from projections.GridConvergence import tm_grid_convergence_projected, tm_grid_convergence_geographic



# Test cases for different locations
test_cases_tm_proj = [
    {   "a": GRS80().a,
        "b": GRS80().b,
        "x": 555776.2667516097,
        "y": 6651832.735433666,
        "false_easting": 500000,
        "lat0": 0,
        "radians": False,
        "description": "Test grid convergences projected",
        "gamma_true": 0.864869193938,  # True meridian convergence
    }
]

@pytest.mark.parametrize("case", test_cases_tm_proj)
def test_tm_grid_convergence_proj(case):
    """
    Test the custom tms_grid_convergence function against known true values.
    """
    a, b = case["a"], case["b"]
    x, y = case["x"], case["y"]
    lat0 = case["lat0"]
    false_easting = case["false_easting"]
    radians = case["radians"]
    gamma_true = case["gamma_true"]
    description = case["description"]

    # Compute grid convergence
    gamma = tm_grid_convergence_projected(a, b, x, y, lat0, false_easting, radians)

    # Assert the grid convergence is close to the expected value
    assert np.isclose(gamma, gamma_true, atol=1e-6), (
        f"Test failed for case: {description}\n"
        f"Computed Grid Convergence: {gamma} degrees\n"
        f"Expected Grid Convergence: {gamma_true} degrees"
    )


### TEST 2

# Test cases for different locations
test_cases_geog = [
    {
        "a": WGS84().a,
        "b": WGS84().b,
        "lon": np.deg2rad(10),
        "lat": np.deg2rad(60),
        "false_easting": 500000,
        "central_meridian": np.deg2rad(9),
        "radians": True,
        "description": "Test 2",
        "gamma_true": 0.8660474985710462,  # True meridian convergence
    }
]

@pytest.mark.parametrize("case", test_cases_geog)
def test_tm_grid_convergence_geog(case):
    """
    Test the custom tms_grid_convergence function against known true values.
    """
    a, b = case["a"], case["b"]
    lon, lat = case["lon"], case["lat"]
    false_easting = case["false_easting"]
    lon0 = case["central_meridian"]
    radians = case["radians"]
    gamma_true = case["gamma_true"]
    description = case["description"]

    # Compute grid convergence
    gamma = np.degrees(tm_grid_convergence_geographic(a, b, lon, lat, lon0, radians=radians))

    # Assert the grid convergence is close to the expected value
    assert np.isclose(gamma, gamma_true, atol=1e-12), (
        f"Test failed for case: {description}\n"
        f"Computed Grid Convergence: {gamma} degrees\n"
        f"Expected Grid Convergence: {gamma_true} degrees"
    )