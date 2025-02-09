"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""

import numpy as np
import pytest
from Ellipsoid import WGS84
from geodetics.geodetic_direct_problem import geodetic_direct_problem

# Define WGS84 ellipsoid parameters
ellip = WGS84()
a = ellip.a
b = ellip.b

# Test cases with known values
test_cases = [
    {
        "lat1": np.radians(52.2296756),  # Start latitude (radians)
        "lon1": np.radians(21.0122287),  # Start longitude (radians)
        "az1": np.radians(-147.4628043168),  # Initial azimuth (radians)
        "d": 1316208.08334,  # Distance (meters)
        "quadrant_correction": False,  # Apply quadrant correction
        "lat2_true": np.radians(41.8919300),  # True destination latitude (radians)
        "lon2_true": np.radians(12.5113300),  # True destination longitude (radians)
        "az2_true": np.radians(-153.7168672618),  # True final azimuth (radians)
        "description": "Test case 1",
    },
    {
        "lat1": np.radians(52.2296756),  # Start latitude (radians)
        "lon1": np.radians(21.0122287),  # Start longitude (radians)
        "az1": np.radians(-147.4628043168),  # Initial azimuth (radians)
        "d": 1316208.08334,  # Distance (meters)
        "quadrant_correction": True,  # Apply quadrant correction
        "lat2_true": np.radians(41.8919300),  # True destination latitude (radians)
        "lon2_true": np.radians(12.5113300),  # True destination longitude (radians)
        "az2_true": np.radians(206.2831327381),  # True final azimuth (radians)
        "description": "Test case 2",
    }
]


@pytest.mark.parametrize("case", test_cases)
def test_geodetic_direct_problem(case):
    """
    Test the geodetic direct problem function against known true values.
    """
    lat1, lon1, az1, d = case["lat1"], case["lon1"], case["az1"], case["d"]
    lat2_true, lon2_true, az2_true = case["lat2_true"], case["lon2_true"], case["az2_true"]
    quadrant_correction = case["quadrant_correction"]
    description = case["description"]

    # Compute the destination point using the function
    lat2, lon2, az2 = geodetic_direct_problem(a, b, lat1, lon1, az1, d, quadrant_correction)

    # Assert that the computed values are close to the true values
    assert np.isclose(lat2, lat2_true, atol=1e-6), (
        f"Test failed for case: {description}\n"
        f"Computed latitude: {np.degrees(lat2):.10f} degrees\n"
        f"Expected latitude: {np.degrees(lat2_true):.10f} degrees"
    )

    assert np.isclose(lon2, lon2_true, atol=1e-6), (
        f"Test failed for case: {description}\n"
        f"Computed longitude: {np.degrees(lon2):.10f} degrees\n"
        f"Expected longitude: {np.degrees(lon2_true):.10f} degrees"
    )

    assert np.isclose(az2, az2_true, atol=1e-6), (
        f"Test failed for case: {description}\n"
        f"Computed final azimuth: {np.degrees(az2):.10f} degrees\n"
        f"Expected final azimuth: {np.degrees(az2_true):.10f} degrees"
    )
