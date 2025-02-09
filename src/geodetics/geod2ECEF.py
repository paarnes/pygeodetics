"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com

Assumes earth is a perfect sphere.

Based on the "THE MERCATOR PROJECTIONS" book from  Peter Osborne, 2013
See "The scale factor for the TMS projection" section at page 63.

"""


from typing import Tuple
import numpy as np
from Ellipsoid import Ellipsoid, WGS84




def geod2ECEF(lat: float, lon: float, h: float, ellipsoid: Ellipsoid = WGS84(), a: float = None, b: float = None) -> Tuple:
    """
    Convert geodetic coordinates (latitude, longitude, height) to ECEF (Earth-Centered, Earth-Fixed)
    Cartesian coordinates.

    Parameters
    ----------
    lat : float. Geodetic latitude in radians.
    lon : float. Geodetic longitude in radians.
    h : float. Height above the ellipsoid in meters.
    ellipsoid : Ellipsoid, optional
        An instance of the Ellipsoid class defining the ellipsoid. Defaults to WGS84.
    a : float, optional. Semi-major axis of the ellipsoid (meters). Overrides `ellipsoid` if provided.
    b : float, optional. Semi-minor axis of the ellipsoid (meters). Overrides `ellipsoid` if provided.

    Returns
    -------
    A tuple containing the ECEF coordinates (X, Y, Z) in meters:
        - X : float. X-coordinate in meters.
        - Y : float. Y-coordinate in meters.
        - Z : float. Z-coordinate in meters.

    Examples
    --------
    >>> lat = np.radians(59.907072474276958)  # Latitude in radians
    >>> lon = np.radians(10.754482924017791)  # Longitude in radians
    >>> h = 63.8281  # Height in meters
    >>> # Example using WGS84 ellipsoid
    >>> X, Y, Z = geod2ECEF(lat, lon, h)
    >>> print(f"X: {X:.4f}, Y: {Y:.4f}, Z: {Z:.4f}")
    """
    # Determine ellipsoid parameters
    if a is None and b is None:
        a, b = ellipsoid.a, ellipsoid.b

    # Compute the prime vertical radius of curvature
    e2 = (a**2 - b**2) / a**2  # Square of the first eccentricity
    N = a / np.sqrt(1 - e2 * np.sin(lat)**2)

    # Calculate ECEF coordinates
    X = (N + h) * np.cos(lat) * np.cos(lon)
    Y = (N + h) * np.cos(lat) * np.sin(lon)
    Z = ((b**2 / a**2) * N + h) * np.sin(lat)

    return X, Y, Z


if __name__ == "__main__":
    # Example usage with WGS84 ellipsoid
    lat = np.radians(59.907072474276958)  # Latitude in radians
    lon = np.radians(10.754482924017791)  # Longitude in radians
    h = 63.8281  # Height in meters

    # Using WGS84 ellipsoid (default)
    X, Y, Z = geod2ECEF(lat, lon, h)
    print(f"Using WGS84 Ellipsoid:\nX: {X:.4f} m\nY: {Y:.4f} m\nZ: {Z:.4f} m")

    # Example with manually specified ellipsoid parameters
    a = 6378137.0  # Semi-major axis
    b = 6356752.314245  # Semi-minor axis
    X, Y, Z = geod2ECEF(lat, lon, h, a=a, b=b)
    print(f"Using Manual Ellipsoid Parameters:\nX: {X:.4f} m\nY: {Y:.4f} m\nZ: {Z:.4f} m")
