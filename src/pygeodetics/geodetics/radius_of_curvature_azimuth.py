"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""

import numpy as np

from .Mrad import Mrad
from .Nrad import Nrad

def radius_of_curvature_azimuth(a: float, b: float, lat: float, az: float, radians: bool = False) -> float:
    """
    Compute the radius of curvature in a given direction is calculated using Euler's equation.
    This considers both the meridional (M) and normal (N) radius of curvature.

    Parameters
    ----------
    a : float. Semi-major axis of the ellipsoid (meters).
    b : float. Semi-minor axis of the ellipsoid (meters).
    lat : float. Geodetic latitude in degrees (if radians=False) or in radians (if radians=True).
    az : float. Azimuth angle in degrees (if radians=False) or in radians (if radians=True).
    radians : bool, optional. If `True`, input and output are in radians. If `False`, input is in degrees and output is converted to degrees. Defaults to `False`.

    Returns
    -------
    float. Radius of curvature for the given azimuth (meters).
    """
    if not radians:
        lat, az = np.radians(lat), np.radians(az)

    # lat is now in radians; pass radians=True to avoid a second conversion
    M = Mrad(a, b, lat, radians=True)  # Compute meridional radius of curvature
    N = Nrad(a, b, lat, radians=True)  # Compute normal radius of curvature

    return (M * N) / (M * np.sin(az)**2 + N * np.cos(az)**2)




if __name__ == '__main__':
    from pygeodetics.Ellipsoid import WGS84

    # Example usage
    ellip = WGS84()
    a, b = ellip.a, ellip.b

    lat = 45.0  # Latitude in degrees
    az = 30.0   # Azimuth in degrees

    radius = radius_of_curvature_azimuth(a, b, lat, az, radians=False)

    print(f"Radius of Curvature at lat {lat}°, az {az}°: {radius:.3f} meters")
