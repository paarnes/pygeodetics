"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""

import numpy as np
import pytest
from Ellipsoid import WGS84
from geodetics.Mrad import Mrad


# Test cases for different locations
test_cases_1 = [
    {
        "a": WGS84().a,
        "b": WGS84().b,
        "lat": 60,
        "radians": False,
        "description": "Test 1: Meridional radius of curvature (M) at latitude 60",
        "Mrad_true": 6383453.85723,
    },
    {
        "a": WGS84().a,
        "b": WGS84().b,
        "lat": 45,
        "radians": False,
        "description": "Test 2: Meridional radius of curvature (M) at latitude 45",
        "Mrad_true": 6367381.81562,
    },
    {
        "a": WGS84().a,
        "b": WGS84().b,
        "lat": 0,
        "radians": False,
        "description": "Test 3: Meridional radius of curvature (M) at latitude 0",
        "Mrad_true": 6335439.32729,
    }
]

@pytest.mark.parametrize("case", test_cases_1)
def test_Mrad(case):
    """
    Test the custom tms_grid_convergence function against known true values.
    """
    a, b = case["a"], case["b"]
    lat = case["lat"]
    radians = case["radians"]
    mrad_true = case["Mrad_true"]
    description = case["description"]

    # Compute grid convergence
    mrad = Mrad(a, b, lat, radians)

    # Assert the grid convergence is close to the expected value
    assert np.isclose(mrad, mrad_true, atol=1e-12), (
        f"Test failed for case: {description}\n"
        f"Computed normal radius of curvature: {mrad} degrees\n"
        f"Expected normal radius of curvature: {mrad_true} degrees"
    )