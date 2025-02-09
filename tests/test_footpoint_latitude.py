import numpy as np
import pytest
from Ellipsoid import WGS84
from geodetics.footpoint_latitude import footpoint_latitude


# Test cases for different locations
test_cases_1 = [
    {
        "a": WGS84().a,
        "b": WGS84().b,
        "y": 6651832.735433666035533,
        "lat0": 0,
        "radians": False,
        "description": "Test 1",
        "footlat_true": 59.979893712987,  # True footpoint lat (degree)
    }
]

@pytest.mark.parametrize("case", test_cases_1)
def test_calc_footpoint_latitude(case):
    """
    Test the custom tms_grid_convergence function against known true values.
    """
    a, b = case["a"], case["b"]
    y = case["y"]
    lat0 = case["lat0"]
    radians = case["radians"]
    footlat_true = case["footlat_true"]
    description = case["description"]

    # Compute grid convergence
    footlat_deg = footpoint_latitude(a, b, y, lat0, radians)

    # Assert the grid convergence is close to the expected value
    assert np.isclose(footlat_deg, footlat_true, atol=1e-12), (
        f"Test failed for case: {description}\n"
        f"Computed Foot point latitude: {footlat_deg} degrees\n"
        f"Expected Foot point latitude: {footlat_true} degrees"
    )