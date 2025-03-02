"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""

import numpy as np
import pytest
from Ellipsoid import WGS84
from geodetics.Nrad import Nrad


# Test cases for different locations
test_cases_1 = [
    {
        "a": WGS84().a,
        "b": WGS84().b,
        "lat": 60,
        "radians": False,
        "description": "Test 1: Normal radius of curvature (N) at latitude 60",
        "Nrad_true": 6394209.17385,
    }
]

@pytest.mark.parametrize("case", test_cases_1)
def test_Nrad(case):
    """
    Test the custom tms_grid_convergence function against known true values.
    """
    a, b = case["a"], case["b"]
    lat = case["lat"]
    radians = case["radians"]
    nrad_true = case["Nrad_true"]
    description = case["description"]

    # Compute grid convergence
    nrad = Nrad(a, b, lat, radians)

    # Assert the grid convergence is close to the expected value
    assert np.isclose(nrad, nrad_true, atol=1e-12), (
        f"Test failed for case: {description}\n"
        f"Computed normal radius of curvature: {nrad} degrees\n"
        f"Expected normal radius of curvature: {nrad_true} degrees"
    )