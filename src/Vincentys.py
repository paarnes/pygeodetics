
from typing import Union
import numpy as np
import Ellipsoid as Ellipsoid

def vincenty_distance(lat1: Union[float, np.ndarray], lon1: Union[float, np.ndarray],
                      lat2: Union[float, np.ndarray], lon2: Union[float, np.ndarray],
                      a: float, b: float, tol: float = 1e-12, max_iter: int = 200) -> Union[float, np.ndarray]:
    """
    Compute the geodesic distance between two points on an ellipsoid using Vincenty's formulae.

    Parameters
    ----------
    lat1 : float or np.ndarray
        Geodetic latitude(s) of the first point(s) in radians.
    lon1 : float or np.ndarray
        Geodetic longitude(s) of the first point(s) in radians.
    lat2 : float or np.ndarray
        Geodetic latitude(s) of the second point(s) in radians.
    lon2 : float or np.ndarray
        Geodetic longitude(s) of the second point(s) in radians.
    a : float
        Semi-major axis of the ellipsoid (meters).
    b : float
        Semi-minor axis of the ellipsoid (meters).
    tol : float, optional
        Convergence tolerance (default: 1e-12).
    max_iter : int, optional
        Maximum number of iterations for convergence (default: 200).

    Returns
    -------
    distance : float or np.ndarray
        The geodesic distance(s) between the two points (meters).

    Examples
    --------
    >>> a = 6378137.0  # Semi-major axis (meters, WGS84)
    >>> b = 6356752.314245  # Semi-minor axis (meters, WGS84)
    >>> lat1, lon1 = np.radians(52.2296756), np.radians(21.0122287)
    >>> lat2, lon2 = np.radians(41.8919300), np.radians(12.5113300)
    >>> distance = vincenty_distance(lat1, lon1, lat2, lon2, a, b)
    >>> print(f"Distance: {distance:.3f} meters")
    """
    # Flattening
    f = (a - b) / a

    # Ensure inputs are numpy arrays
    lat1 = np.asarray(lat1,dtype=np.float64)
    lon1 = np.asarray(lon1,dtype=np.float64)
    lat2 = np.asarray(lat2,dtype=np.float64)
    lon2 = np.asarray(lon2,dtype=np.float64)

    # Difference in longitudes
    L = lon2 - lon1

    # Reduced latitudes
    U1 = np.arctan((1 - f) * np.tan(lat1))
    U2 = np.arctan((1 - f) * np.tan(lat2))

    sinU1, cosU1 = np.sin(U1), np.cos(U1)
    sinU2, cosU2 = np.sin(U2), np.cos(U2)

    # Initialize iteration variables
    lamb = L  # Initial lambda = difference in longitude
    iter_count = 0

    while True:
        sin_lambda = np.sin(lamb)
        cos_lambda = np.cos(lamb)

        sin_sigma = np.sqrt((cosU2 * sin_lambda)**2 +
                            (cosU1 * sinU2 - sinU1 * cosU2 * cos_lambda)**2)
        cos_sigma = sinU1 * sinU2 + cosU1 * cosU2 * cos_lambda
        sigma = np.arctan2(sin_sigma, cos_sigma)

        sin_alpha = cosU1 * cosU2 * sin_lambda / sin_sigma
        cos2_alpha = 1 - sin_alpha**2
        cos2_sigma_m = cos_sigma - 2 * sinU1 * sinU2 / cos2_alpha

        C = f / 16 * cos2_alpha * (4 + f * (4 - 3 * cos2_alpha))

        lamb_prev = lamb
        lamb = L + (1 - C) * f * sin_alpha * (
            sigma + C * sin_sigma * (cos2_sigma_m + C * cos_sigma *
                                     (-1 + 2 * cos2_sigma_m**2)))

        # Convergence check
        if np.all(np.abs(lamb - lamb_prev) < tol) or iter_count >= max_iter:
            break
        iter_count += 1

    if iter_count >= max_iter:
        raise RuntimeError("Vincenty's formula failed to converge")

    # Compute distance
    u2 = cos2_alpha * (a**2 - b**2) / b**2
    A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))

    delta_sigma = B * sin_sigma * (
        cos2_sigma_m + B / 4 *
        (cos_sigma * (-1 + 2 * cos2_sigma_m**2) -
         B / 6 * cos2_sigma_m * (-3 + 4 * sin_sigma**2) *
         (-3 + 4 * cos2_sigma_m**2)))

    distance = b * A * (sigma - delta_sigma)
    return distance


if __name__ == "__main__":
    # # Example usage
    # a = 6378137.0  # Semi-major axis (meters, WGS84)
    # b = 6356752.314245  # Semi-minor axis (meters, WGS84)

    # # Example 1: Single pair of coordinates
    # lat1 = np.radians(52.2296756)
    # lon1 = np.radians(21.0122287)
    # lat2 = np.radians(41.8919300)
    # lon2 = np.radians(12.5113300)

    # distance = vincenty_distance(lat1, lon1, lat2, lon2, a, b)
    # print(f"Distance (single point): {distance:.3f} meters")

    # # Example 2: Multiple coordinates
    # lat1 = np.radians([52.2296756, 48.856614])
    # lon1 = np.radians([21.0122287, 2.3522219])
    # lat2 = np.radians([41.8919300, 51.507222])
    # lon2 = np.radians([12.5113300, -0.1275])

    # distances = vincenty_distance(lat1, lon1, lat2, lon2, a, b)
    # print(f"Distances (multiple points): {distances}")


    from pygeodesy.ellipsoidalVincenty import LatLon

    ellip = Ellipsoid.WGS84()
    a = ellip.a
    b = ellip.b
    # Example 1: Single pair of coordinates
    lat1 = 52.2296756
    lon1 = 21.0122287
    lat2 = 41.8919300
    lon2 = 12.5113300

    # Convert to radians for vincenty_distance
    lat1_rad = np.radians(lat1)
    lon1_rad = np.radians(lon1)
    lat2_rad = np.radians(lat2)
    lon2_rad = np.radians(lon2)

    # Calculate distance using your vincenty_distance function
    distance_custom = vincenty_distance(lat1_rad, lon1_rad, lat2_rad, lon2_rad, a, b)
    print(f"Custom Vincenty Distance: {distance_custom:.12f} meters")

    # Calculate distance using pygeodesy
    p1 = LatLon(lat1, lon1)
    p2 = LatLon(lat2, lon2)
    distance_pygeodesy = p1.distanceTo(p2)
    print(f"PyGeodesy Vincenty Distance: {distance_pygeodesy:.12f} meters")
