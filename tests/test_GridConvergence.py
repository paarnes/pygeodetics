import numpy as np
import pytest
from Ellipsoid import WGS84
from projections.GridConvergence import tms_grid_convergence_projected, tm_grid_convergence_geographic

# Define WGS84 ellipsoid parameters
ellip = WGS84()
a = ellip.a
b = ellip.b

# Test cases for different locations
test_cases_tms = [
    {
        "x": 555776.266751609742641,
        "y": 6651832.735433666035533,
        "false_easting": 500000,
        "angle_unit": "deg",
        "description": "Test 1",
        "gamma_true": 0.8660474985710462,  # True meridian convergence
    }
]

@pytest.mark.parametrize("case", test_cases_tms)
def test_tms_grid_convergence(case):
    """
    Test the custom tms_grid_convergence function against known true values.
    """
    x, y = case["x"], case["y"]
    false_easting = case["false_easting"]
    angle_unit = case["angle_unit"]
    gamma_true = case["gamma_true"]
    description = case["description"]

    # Compute grid convergence
    gamma = tms_grid_convergence_projected(x, y, false_easting, angle_unit=angle_unit)

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