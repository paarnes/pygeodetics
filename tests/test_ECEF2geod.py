"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""

import numpy as np
import pytest
from Ellipsoid import WGS84
from geodetics.ECEF2geod import ECEF2geod, ECEF2geodb, ECEF2geodv

# Define WGS84 ellipsoid parameters
ellip = WGS84()
a = ellip.a
b = ellip.b

# Define Cartesian ECEF coordinates and true geodetic values
X, Y, Z = 3149785.9652, 598260.8822, 5495348.4927
true_lat = 59.907072474276958  # in degrees
true_lon = 10.754482924017791  # in degrees
true_h = 63.8281  # in meters


@pytest.mark.parametrize("func", [ECEF2geod, ECEF2geodb, ECEF2geodv])
def test_ecef2geod(func):
    """
    Test all ECEF2geod functions against the true values.
    """
    lat, lon, h = func(a, b, X, Y, Z, angle_unit="deg")

    # Compare against true values
    assert np.isclose(lat, true_lat, atol=1e-6), (
        f"{func.__name__}: Latitude mismatch\nExpected: {true_lat}\nGot: {lat}"
    )
    assert np.isclose(lon, true_lon, atol=1e-6), (
        f"{func.__name__}: Longitude mismatch\nExpected: {true_lon}\nGot: {lon}"
    )
    assert np.isclose(h, true_h, atol=1e-3), (
        f"{func.__name__}: Height mismatch\nExpected: {true_h}\nGot: {h}"
    )

