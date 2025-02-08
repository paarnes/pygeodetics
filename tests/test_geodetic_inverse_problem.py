"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""


import numpy as np
import pytest
from Ellipsoid import WGS84
from geodetics.geodetic_inverse_problem import geodetic_inverse_problem


# Define WGS84 ellipsoid parameters
ellip = WGS84()
a = ellip.a
b = ellip.b

# Test case: Warsaw to Rome
test_cases = [
    {
        "a": a,
        "b": b,
        "lat1": np.radians(52.2296756),  # Point 1 latitude (in radians)
        "lon1": np.radians(21.0122287),  # Point 1 longitude (in radians)
        "lat2": np.radians(41.8919300),  # Point 2 latitude (in radians)
        "lon2": np.radians(12.5113300),  # Point 2 longitude (in radians)
        "quadrant_correction": False,  # Default behavior
        "description": "Test case 1",
        "az1_true": np.radians(-147.4628043168),  # Expected Forward Azimuth (radians)
        "az2_true": np.radians(-153.7168672619),  # Expected Reverse Azimuth (radians)
        "s_true": 1316208.0833061778,  # Expected Distance (meters)
    }
]

@pytest.mark.parametrize("case", test_cases)
def test_geodetic_inverse_problem(case):
    """
    Test the geodetic inverse problem function against known true values.
    """
    a, b = case["a"], case["b"]
    lat1, lon1, lat2, lon2 = case["lat1"], case["lon1"], case["lat2"], case["lon2"]
    quadrant_correction = case["quadrant_correction"]
    az1_true, az2_true, s_true = case["az1_true"], case["az2_true"], case["s_true"]
    description = case["description"]

    # Compute geodetic inverse solution
    az1, az2, s = geodetic_inverse_problem(a, b, lat1, lon1, lat2, lon2, quadrant_correction)

    # Assert the computed values are close to the expected values
    assert np.isclose(az1, az1_true, atol=1e-10), (
        f"Test failed for case: {description}\n"
        f"Computed Forward Azimuth: {np.degrees(az1):.10f}°\n"
        f"Expected Forward Azimuth: {np.degrees(az1_true):.10f}°"
    )

    assert np.isclose(az2, az2_true, atol=1e-10), (
        f"Test failed for case: {description}\n"
        f"Computed Reverse Azimuth: {np.degrees(az2):.10f}°\n"
        f"Expected Reverse Azimuth: {np.degrees(az2_true):.10f}°"
    )

    assert np.isclose(s, s_true, atol=1e-6), (
        f"Test failed for case: {description}\n"
        f"Computed Distance: {s:.10f} meters\n"
        f"Expected Distance: {s_true:.10f} meters"
    )
