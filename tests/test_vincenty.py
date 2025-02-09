"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""


import numpy as np
import pytest
from pygeodesy.ellipsoidalVincenty import LatLon
from geodetics.Vincenty import vincenty_distance  # Replace `your_module` with the actual module name
import Ellipsoid as Ellipsoid



# Define WGS84 ellipsoid parameters
ellip = Ellipsoid.WGS84()
a = ellip.a
b = ellip.b

# Test cases
test_cases = [
    # Crossing the equator
    {"lat1": 1.0, "lon1": 30.0, "lat2": -1.0, "lon2": 30.0, "description": "Crossing the equator"},
    # Crossing the anti-meridian
    {"lat1": 10.0, "lon1": 179.9, "lat2": 10.0, "lon2": -179.9, "description": "Crossing the anti-meridian"},
    # Over the North Pole
    {"lat1": 89.9, "lon1": 45.0, "lat2": 89.9, "lon2": -135.0, "description": "Over the North Pole"},
    # Over the South Pole
    {"lat1": -89.9, "lon1": 45.0, "lat2": -89.9, "lon2": -135.0, "description": "Over the South Pole"},
    # Arbitrary locations
    {"lat1": 52.2296756, "lon1": 21.0122287, "lat2": 41.8919300, "lon2": 12.5113300, "description": "Warsaw to Rome"},
]

@pytest.mark.parametrize("case", test_cases)
def test_vincenty(case):
    """
    Test the custom vincenty_distance function against pygeodesy's implementation.
    """
    lat1, lon1, lat2, lon2 = case["lat1"], case["lon1"], case["lat2"], case["lon2"]
    description = case["description"]

    # Convert to radians for the custom function
    lat1_rad = np.radians(lat1)
    lon1_rad = np.radians(lon1)
    lat2_rad = np.radians(lat2)
    lon2_rad = np.radians(lon2)

    # Calculate distance using the custom vincenty_distance function
    custom_distance = vincenty_distance(lon1_rad, lat1_rad, lon2_rad, lat2_rad, a=a, b=b)

    # Calculate distance using pygeodesy
    p1 = LatLon(lat1, lon1)
    p2 = LatLon(lat2, lon2)
    pygeodesy_distance = p1.distanceTo(p2)

    # Assert the distances match within a small tolerance
    assert np.isclose(custom_distance, pygeodesy_distance, atol=1e-6), (
        f"Test failed for case: {description}\n"
        f"Custom distance: {custom_distance} m\n"
        f"PyGeodesy distance: {pygeodesy_distance} m"
    )
