import numpy as np
from Ellipsoid import WGS84


def geod2_indir(a: float, b: float, lat1: float, lon1: float, lat2: float, lon2: float):
    """
    Solve the indirect geodetic problem: compute the forward and reverse azimuths and the geodesic distance
    between two points on an ellipsoid.

    DONT GIVE CORRECT RESULTS.

    Parameters
    ----------
    a : float. Semi-major axis of the ellipsoid (meters).
    b : float. Semi-minor axis of the ellipsoid (meters).
    lat1 : float. Latitude of point 1 in radians.
    lon1 : float. Longitude of point 1 in radians.
    lat2 : float. Latitude of point 2 in radians.
    lon2 : float. Longitude of point 2 in radians.

    Returns
    -------
    az1 : float. Forward azimuth at point 1 (radians).
    az2 : float. Reverse azimuth at point 2 (radians).
    s : float. Geodesic distance between the two points (meters).

    Notes
    -----
    This function implements the iterative solution of Vincenty's formula for the geodesic on an ellipsoid.
    """
    epsilon = 1e-10
    dlon_old = -1  # Initialize to a value that guarantees iteration

    # Flattening and eccentricity squared
    f = (a - b) / a
    e2 = (a**2 - b**2) / a**2

    # Reduced latitudes
    beta1 = np.arctan((b / a) * np.tan(lat1))
    beta2 = np.arctan((b / a) * np.tan(lat2))

    # Initial difference in longitude
    dlon = lon2 - lon1
    dlon_new = dlon

    while np.abs(dlon_new - dlon_old) > epsilon:
        dlon_old = dlon_new

        # Compute intermediate values
        X = np.cos(beta1) * np.sin(beta2) - np.sin(beta1) * np.cos(beta2) * np.cos(dlon_old)
        Y = np.cos(beta2) * np.sin(dlon_old)
        Z = np.sin(beta1) * np.sin(beta2) + np.cos(beta1) * np.cos(beta2) * np.cos(dlon_old)

        sigma = np.arctan2(np.sqrt(X**2 + Y**2), Z)
        az1 = np.arctan2(Y, X)
        az0 = np.arcsin(np.sin(az1) * np.cos(beta1))

        sigma1 = np.arctan2(np.tan(beta1), np.cos(az1))
        sigma2 = sigma1 + sigma

        K = (f + f**2 / 4) * np.cos(az0)**2 - (f**2 / 4) * np.cos(az0)**4

        # Update longitude difference
        dlon_new = dlon + f * np.sin(az0) * (
            (1 - K - K**2) * sigma +
            K * np.sin(sigma) * np.cos(sigma1 + sigma2) +
            K**2 * np.sin(sigma) * np.cos(sigma) * np.cos(2 * (sigma1 + sigma2))
        )

    dlon = dlon_new

    # Compute reverse azimuth
    az2 = np.arctan2(np.cos(beta1) * np.sin(dlon),
                     np.cos(beta1) * np.sin(beta2) * np.cos(dlon) - np.sin(beta1) * np.cos(beta2))

    # Ensure az2 lies in the range [0, 2π]
    az2 = np.mod(az2, 2 * np.pi)

    # Compute auxiliary values for distance calculation
    g = e2 * np.cos(az0)**2
    H = (1 / 8) * g - (1 / 16) * g**2 + (37 / 1024) * g**3
    b0 = b * (1 + (1 / 4) * g - (3 / 64) * g**2 + (5 / 256) * g**3)

    # Compute geodesic distance
    s = b0 * (
        sigma -
        2 * H * np.sin(sigma) * np.cos(sigma1 + sigma2) -
        (H**2 / 2) * np.sin(2 * sigma) * np.cos(2 * (sigma1 + sigma2)) -
        (H**3 / 3) * np.sin(3 * sigma) * np.cos(3 * (sigma1 + sigma2))
    )

    return az1, az2, s



if __name__ == "__main__":
    ellip = WGS84()
    a = ellip.a
    b = ellip.b
    f = ellip.f
    lat1 = np.radians(52.2296756)  # Point 1 latitude (in radians)
    lon1 = np.radians(21.0122287)  # Point 1 longitude (in radians)
    lat2 = np.radians(41.8919300)  # Point 2 latitude (in radians)
    lon2 = np.radians(12.5113300)  # Point 2 longitude (in radians)

    az1, az2, distance = geod2_indir(a, b, lat1, lon1, lat2, lon2)
    # print(f"Forward Azimuth: {np.degrees(az1):.6f} degrees")
    # print(f"Reverse Azimuth: {np.degrees(az2):.6f} degrees")
    # print(f"Distance: {distance:.3f} meters")

    # Validate using pyproj.Geod
    from pyproj import Geod
    geod = Geod(a=a, b=b)
    az1_pyproj, az2_pyproj, s_pyproj = geod.inv(np.degrees(lon1), np.degrees(lat1), np.degrees(lon2), np.degrees(lat2))

    # Geographlib
    from geographiclib.geodesic import Geodesic
    geod_lib = Geodesic(a, f)
    data_geoglib = geod_lib.Inverse(np.degrees(lat1),np.degrees(lon1), np.degrees(lat2), np.degrees(lon2))

    print("Custom geod2_indir results:")
    print(f"Forward Azimuth: {np.degrees(az1):.6f} degrees")
    print(f"Reverse Azimuth: {np.degrees(az2):.6f} degrees")
    print(f"Distance: {distance:.3f} meters\n")

    print("Validation using pyproj.Geod:")
    print(f"Forward Azimuth: {az1_pyproj:.6f} degrees")
    print(f"Reverse Azimuth: {az2_pyproj:.6f} degrees")
    print(f"Distance: {s_pyproj:.3f} meters\n")

    print("Validation using GEOGLIB:")
    print(f"Forward Azimuth: {data_geoglib["azi1"]:.6f} degrees")
    print(f"Reverse Azimuth: {data_geoglib["azi2"]:.6f} degrees")
    print(f"Forward Azimuth: {data_geoglib["s12"]:.6f} degrees")

    # # Compare results
    # assert np.isclose(np.degrees(az1), az1_pyproj, atol=1e-6), "Forward azimuth mismatch!"
    # assert np.isclose(np.degrees(az2), az2_pyproj, atol=1e-6), "Reverse azimuth mismatch!"
    # assert np.isclose(distance, s_pyproj, atol=1e-6), "Distance mismatch!"
    # print("Validation passed: Results match pyproj.")
