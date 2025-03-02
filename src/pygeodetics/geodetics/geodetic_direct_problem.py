"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""

from typing import Tuple
import numpy as np

def geodetic_direct_problem(
    a: float,
    b: float,
    lat1: float,
    lon1: float,
    az1: float,
    d: float,
    quadrant_correction: bool = False,
    radians: bool = False
) -> Tuple[float, float, float]:
    """
    Solve the geodetic direct problem: Compute the destination coordinates and final azimuth
    given an initial point, azimuth, and distance.

    Notes
    -----
    This function implements Vincenty's direct method to solve the geodetic direct problem.
    It calculates the destination latitude, longitude, and final azimuth after traveling
    a given distance along a given azimuth on an ellipsoid.

    Parameters
    ----------
    a : float. Semi-major axis of the ellipsoid (meters).
    b : float. Semi-minor axis of the ellipsoid (meters).
    lat1 : float. Latitude of the initial point in radians.
    lon1 : float. Longitude of the initial point in radians.
    az1 : float. Forward azimuth at the initial point in radians.
    d : float. Distance to travel along the geodesic (meters).
    quadrant_correction : bool, optional. If True, ensures azimuths
        are in the range [0, 2π]. Default is False.

    Returns
    -------
    lat2 : float. Latitude of the destination point
        (degrees if `radians=False`, radians if `radians=True`).
    lon2 : float. Longitude of the destination point
        (degrees if `radians=False`, radians if `radians=True`).
    az2 : float. Reverse azimuth at the destination point
        (degrees if `radians=False`, radians if `radians=True`).
    """

    # Convert inputs to radians if they are in degrees
    if not radians:
        lat1, lon1, az1 = np.radians(lat1), np.radians(lon1), np.radians(az1)

    # Compute flattening and second eccentricity squared
    f = (a - b) / a
    e2m = (a**2 - b**2) / b**2

    # Compute reduced latitude
    beta1 = np.arctan((b / a) * np.tan(lat1))
    az0 = np.arcsin(np.sin(az1) * np.cos(beta1))
    sigma1 = np.arctan(np.tan(beta1) / np.cos(az1))

    # Compute auxiliary variables
    g = e2m * np.cos(az0)**2
    H = (1 / 8) * g - (1 / 16) * g**2 + (37 / 1024) * g**3
    b0 = b * (1 + (1 / 4) * g - (3 / 64) * g**2 + (5 / 256) * g**3)

    # Compute initial geodesic arc length
    d1 = b0 * (sigma1 - H * np.sin(2 * sigma1) - (H**2 / 4) * np.sin(4 * sigma1) - (H**3 / 6) * np.sin(6 * sigma1))
    d2 = d1 + d

    # Compute new geodesic arc length
    sigma2 = (d2 / b0) + (H - (3 / 4) * H**3) * np.sin(2 * d2 / b0) + (5 / 4) * H**2 * np.sin(4 * d2 / b0) + (29 / 12) * H**3 * np.sin(6 * d2 / b0)
    sigma = sigma2 - sigma1

    # Compute destination point
    X = np.cos(beta1) * np.cos(sigma) - np.sin(beta1) * np.sin(sigma) * np.cos(az1)
    Y = np.sin(sigma) * np.sin(az1)
    Z = np.sin(beta1) * np.cos(sigma) + np.cos(beta1) * np.sin(sigma) * np.cos(az1)

    beta2 = np.arctan(Z / np.sqrt(X**2 + Y**2))
    dlon = np.arctan2(Y, X)

    # Compute correction for longitude
    K = ((f + f**2) / 4) * np.cos(az0)**2 - (f**2 / 4) * np.cos(az0)**4
    dlon = dlon - f * np.sin(az0) * (
        (1 - K - K**2) * sigma + K * np.sin(sigma) * np.cos(sigma1 + sigma2)
        + K**2 * np.sin(sigma) * np.cos(sigma) * np.cos(2 * (sigma1 + sigma2))
    )

    lat2 = np.arctan((a / b) * np.tan(beta2))
    lon2 = lon1 + dlon

    # Compute reverse azimuth
    az2 = np.arctan2(np.sin(az1) * np.cos(beta1), (np.cos(beta1) * np.cos(sigma) * np.cos(az1) - np.sin(beta1) * np.sin(sigma)))

    # Adjust azimuth to be in range [0, 2π] if requested
    if quadrant_correction:
        az2 = np.mod(az2, 2 * np.pi)

    # Convert outputs to degrees if input was in degrees
    if not radians:
        lat2, lon2, az2 = np.degrees(lat2), np.degrees(lon2), np.degrees(az2)

    return lat2, lon2, az2





if __name__ == "__main__":
    import sys, os
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
    from Ellipsoid import WGS84

    ellip = WGS84()
    a = ellip.a
    b = ellip.b

    lat1 = np.radians(52.2296756)  # Start latitude (in radians)
    lon1 = np.radians(21.0122287)  # Start longitude (in radians)
    az1 = np.radians(-147.4628043168) # Initial azimuth (in radians)
    d = 1316208.08334 # Distance in meters

    lat2, lon2, az2 = geodetic_direct_problem(a, b, lat1, lon1, az1, d, quadrant_correction=True, radians=True)

    print(f"Destination Latitude: {np.degrees(lat2):.10f} degrees")
    print(f"Destination Longitude: {np.degrees(lon2):.10f} degrees")
    print(f"Final Azimuth: {np.degrees(az2):.10f} degrees")
