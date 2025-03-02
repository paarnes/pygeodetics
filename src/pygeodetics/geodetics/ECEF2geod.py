
"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com

Three differnt function/approcaches to convert Cartesian ECEF coordinates to geodetic coordinates.

Recommended to use ECEF2geodb or ECEF2geodv, they are faster and slightly more accurate than the iterative method.
Both these methods give the same results as pyproj. ECEF2geodv seems to be slightly faster than ECEF2geodb.

"""
from typing import Literal, Tuple
import numpy as np


def ECEF2geod(a: float, b: float, X: float, Y: float,
              Z: float, angle_unit: Literal["deg","rad"]= "deg") -> Tuple[float, float, float]:
    """
    Convert Cartesian ECEF (Earth-Centered, Earth-Fixed) coordinates to geodetic coordinates using iteration.

    Parameters
    ----------
    a : float. Semi-major axis of the ellipsoid (meters).
    b : float. Semi-minor axis of the ellipsoid (meters).
    X : float. ECEF X coordinate (meters).
    Y : float. ECEF Y coordinate (meters).
    Z : float. ECEF Z coordinate (meters).

    Returns
    -------
    lat : float. Latitude in radians.
    lon : float. Longitude in radians.
    h : float. Height above the ellipsoid (meters).

    Examples
    --------
    >>> a = 6378137.0
    >>> b = 6356752.314245
    >>> X, Y, Z = 4510731.0, 4510731.0, 4510731.0
    >>> lat, lon, h = ECEF2geod(a, b, X, Y, Z)
    >>> print(lat, lon, h)
    """
    e2 = (a**2 - b**2) / a**2  # Square of the first eccentricity
    p = np.sqrt(X**2 + Y**2)
    epsilon = 1e-10  # Convergence threshold

    # Initial latitude estimate
    lat = np.arctan2(Z, p)
    lat_new = 0

    # Iterative process
    while np.abs(lat_new - lat) > epsilon:
        lat = lat_new or lat
        N = a / np.sqrt(1 - e2 * np.sin(lat)**2)  # Prime vertical radius of curvature
        lat_new = np.arctan2(Z + N * e2 * np.sin(lat), p)

    lat = lat_new
    lon = np.arctan2(Y, X)  # Longitude
    N = a / np.sqrt(1 - e2 * np.sin(lat)**2)
    h = p / np.cos(lat) - N  # Height above ellipsoid

    if angle_unit == "deg":
        lat = np.degrees(lat)
        lon = np.degrees(lon)

    return lat, lon, h


def ECEF2geodb(a: float, b: float, X: float, Y: float,
              Z: float, angle_unit: Literal["deg","rad"]= "deg") -> Tuple[float, float, float]:
    """
    Convert Cartesian ECEF coordinates to geodetic coordinates using Bowring's method.
    Faster and slightly more accurate than the iterative method.

    Parameters
    ----------
    a : float
        Semi-major axis of the ellipsoid (meters).
    b : float
        Semi-minor axis of the ellipsoid (meters).
    X : float
        ECEF X coordinate (meters).
    Y : float
        ECEF Y coordinate (meters).
    Z : float
        ECEF Z coordinate (meters).

    Returns
    -------
    lat : float. Latitude in radians.
    lon : float. Longitude in radians.
    h : float. Height above the ellipsoid (meters).

    Examples
    --------
    >>> a = 6378137.0
    >>> b = 6356752.314245
    >>> X, Y, Z = 4510731.0, 4510731.0, 4510731.0
    >>> lat, lon, h = ECEF2geodb(a, b, X, Y, Z)
    >>> print(lat, lon, h)
    """
    e2m = (a**2 - b**2) / b**2  # Second eccentricity squared
    e2 = (a**2 - b**2) / a**2  # First eccentricity squared
    rho = np.sqrt(X**2 + Y**2)

    my = np.arctan2(Z * a, rho * b)
    lat = np.arctan2(Z + e2m * b * np.sin(my)**3, rho - e2 * a * np.cos(my)**3)
    lon = np.arctan2(Y, X)
    N = a / np.sqrt(1 - e2 * np.sin(lat)**2)
    h = rho * np.cos(lat) + Z * np.sin(lat) - N * (1 - e2 * np.sin(lat)**2)

    if angle_unit == "deg":
        lat = np.degrees(lat)
        lon = np.degrees(lon)
    return lat, lon, h


def ECEF2geodv(a: float, b: float, X: float, Y: float,
              Z: float, angle_unit: Literal["deg","rad"]= "deg") -> Tuple[float, float, float]:
    """
    Convert Cartesian ECEF coordinates to geodetic coordinates using Vermeille's method.
    Faster and slightly more accurate than the iterative method. Gives excact same results as pyproj.
    Generates the same results as ECEF2geodb, but seems to be slightly faster.

    Parameters
    ----------
    a : float. Semi-major axis of the ellipsoid (meters).
    b : float. Semi-minor axis of the ellipsoid (meters).
    X : float. ECEF X coordinate (meters).
    Y : float. ECEF Y coordinate (meters).
    Z : float. ECEF Z coordinate (meters).

    Returns
    -------
    lat : float
        Latitude in radians.
    lon : float
        Longitude in radians.
    h : float
        Height above the ellipsoid (meters).

    Examples
    --------
    >>> a = 6378137.0
    >>> b = 6356752.314245
    >>> X, Y, Z = 4510731.0, 4510731.0, 4510731.0
    >>> lat, lon, h = ECEF2geodv(a, b, X, Y, Z)
    >>> print(lat, lon, h)
    """
    e2 = (a**2 - b**2) / a**2  # Square of the first eccentricity
    p = (X**2 + Y**2) / a**2
    q = ((1 - e2) * Z**2) / a**2
    r = (p + q - e2**2) / 6

    s = e2**2 * p * q / (4 * r**3)
    t = (1 + s + np.sqrt(s * (2 + s)))**(1 / 3)
    u = r * (1 + t + 1 / t)
    v = np.sqrt(u**2 + e2**2 * q)
    w = e2 * (u + v - q) / (2 * v)
    k = np.sqrt(u + v + w**2) - w
    D = k * np.sqrt(X**2 + Y**2) / (k + e2)

    lat = 2 * np.arctan2(Z, D + np.sqrt(D**2 + Z**2))
    lon = np.arctan2(Y, X)
    h = ((k + e2 - 1) / k) * np.sqrt(D**2 + Z**2)

    if angle_unit == "deg":
        lat = np.degrees(lat)
        lon = np.degrees(lon)

    return lat, lon, h



if __name__ == "__main__":
    import sys
    import os
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
    from Ellipsoid import WGS84
    from pyproj import Transformer

    # Get the WGS84 ellipsoid parameters
    a = WGS84().a
    b = WGS84().b

    # Define Cartesian ECEF coordinates
    X, Y, Z = 3149785.9652, 598260.8822, 5495348.4927

    # Use ECEF2geod
    lat, lon, h = ECEF2geod(a, b, X, Y, Z)
    print(f"ECEF2geod:\nLatitude={lat:.15f}\nLongitude={lon:.15f}, Height={h:.4f}\n")

    # Use ECEF2geodb
    lat, lon, h = ECEF2geodb(a, b, X, Y, Z)
    print(f"ECEF2geodb:\nLatitude={lat:.15f}\nLongitude={lon:.15f}, Height={h:.4f}\n")

    # Use ECEF2geodv
    lat, lon, h = ECEF2geodv(a, b, X, Y, Z)
    print(f"ECEF2geodv:\nLatitude={lat:.15f}\nLongitude={lon:.15f}, Height={h:.4f}\n")

    # Use pyproj to verify the results
    transformer = Transformer.from_crs("epsg:4978", "epsg:4326", always_xy=True)
    lon, lat, h = transformer.transform(X, Y, Z)
    print(f"pyproj:\nLatitude={lat:.15f}\nLongitude={lon:.15f}, Height={h:.4f}\n")