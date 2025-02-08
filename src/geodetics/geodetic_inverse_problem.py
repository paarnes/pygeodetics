"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""


import numpy as np


def geodetic_inverse_problem(a: float, b: float, lat1: float, lon1: float,
                             lat2: float, lon2: float, quadrant_correction: bool = False):
    """
    Solve the geodetic indirect problem.

    Parameters
    ----------
    a : float. Semi-major axis of the ellipsoid (meters).
    b : float. Semi-minor axis of the ellipsoid (meters).
    lat1 : float. Latitude of the first point in radians.
    lon1 : float. Longitude of the first point in radians.
    lat2 : float. Latitude of the second point in radians.
    lon2 : float. Longitude of the second point in radians.
    quadrant_correction : bool, optional. If True, ensures
        azimuths are in the range [0, 2π]. Default is False.

    Returns
    -------
    az1 : float. Forward azimuth at point 1 (radians).
    az2 : float. Reverse azimuth at point 2 (radians).
    d : float. Geodesic distance between the two points (meters).
    """
    f = (a - b) / a
    e2m = (a**2 - b**2) / b**2

    beta1 = np.arctan(b / a * np.tan(lat1))
    beta2 = np.arctan(b / a * np.tan(lat2))

    epsilon = 1e-10
    dlon_new = lon2 - lon1
    dlon = 0

    while np.abs(dlon_new - dlon) > epsilon:
        dlon = dlon_new

        X = np.cos(beta1) * np.sin(beta2) - np.sin(beta1) * np.cos(beta2) * np.cos(dlon)
        Y = np.cos(beta2) * np.sin(dlon)
        Z = np.sin(beta1) * np.sin(beta2) + np.cos(beta1) * np.cos(beta2) * np.cos(dlon)

        sigma = np.arctan(np.sqrt(X**2 + Y**2) / Z)
        az1 = np.arctan2(Y, X)
        az0 = np.arcsin(np.sin(az1) * np.cos(beta1))

        sigma1 = np.arctan(np.tan(beta1) / np.cos(az1))
        sigma2 = sigma1 + sigma

        K = (f + f**2) / 4 * np.cos(az0)**2 - f**2 / 4 * np.cos(az0)**4

        dlon_new = (lon2 - lon1) + f * np.sin(az0) * (
            (1 - K - K**2) * sigma +
            K * np.sin(sigma) * np.cos(sigma1 + sigma2) +
            K**2 * np.sin(sigma) * np.cos(sigma) * np.cos(2 * (sigma1 + sigma2))
        )

    dlon = dlon_new
    az2 = np.arctan2(np.cos(beta1) * np.sin(dlon),
                     np.cos(beta1) * np.sin(beta2) * np.cos(dlon) - np.sin(beta1) * np.cos(beta2))

    # Apply quadrant correction if requested
    if quadrant_correction:
        az1 = np.mod(az1, 2 * np.pi)
        az2 = np.mod(az2, 2 * np.pi)

    g = e2m * np.cos(az0)**2
    H = 1/8 * g - 1/16 * g**2 + 37/1024 * g**3
    b0 = b * (1 + 1/4 * g - 3/64 * g**2 + 5/256 * g**3)

    d = b0 * (
        sigma - 2 * H * np.sin(sigma) * np.cos(sigma1 + sigma2)
        - H**2 / 2 * np.sin(2 * sigma) * np.cos(2 * (sigma1 + sigma2))
        - H**3 / 3 * np.sin(3 * sigma) * np.cos(3 * (sigma1 + sigma2))
    )

    return az1, az2, d




if __name__ == "__main__":
    import sys, os
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
    from Ellipsoid import WGS84

    ellip = WGS84()
    a = ellip.a
    b = ellip.b
    f = ellip.f
    lat1 = np.radians(52.2296756)  # Point 1 latitude (in radians)
    lon1 = np.radians(21.0122287)  # Point 1 longitude (in radians)
    lat2 = np.radians(41.8919300)  # Point 2 latitude (in radians)
    lon2 = np.radians(12.5113300)  # Point 2 longitude (in radians)

    az1, az2, distance = geodetic_inverse_problem(a, b, lat1, lon1, lat2, lon2)

    # Validate using pyproj.Geod
    from pyproj import Geod
    geod = Geod(a=a, b=b)
    az1_pyproj, az2_pyproj, s_pyproj = geod.inv(np.degrees(lon1), np.degrees(lat1), np.degrees(lon2), np.degrees(lat2))

    # Geographlib
    from geographiclib.geodesic import Geodesic
    geod_lib = Geodesic(a, f)
    data_geoglib = geod_lib.Inverse(np.degrees(lat1),np.degrees(lon1), np.degrees(lat2), np.degrees(lon2))

    print("Custom geod2_indir results:")
    print(f"Forward Azimuth: {np.degrees(az1):.10f} degrees")
    print(f"Reverse Azimuth: {np.degrees(az2):.10f} degrees")
    print(f"Distance: {distance:.10f} meters\n")

    print("Validation using pyproj.Geod:")
    print(f"Forward Azimuth: {az1_pyproj:.10f} degrees")
    print(f"Reverse Azimuth: {az2_pyproj:.10f} degrees")
    print(f"Distance: {s_pyproj:.10f} meters\n")

    print("Validation using GEOGLIB:")
    print(f"Forward Azimuth: {data_geoglib["azi1"]:.10f} degrees")
    print(f"Reverse Azimuth: {data_geoglib["azi2"]:.10f} degrees")
    print(f"Forward Azimuth: {data_geoglib["s12"]:.10f} degrees")
