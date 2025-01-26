"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com


https://en.wikipedia.org/wiki/Transverse_Mercator_projection
https://en.wikipedia.org/wiki/Scale_(map)

"""

import sys
import os
from typing import Tuple
import numpy as np
from pyproj import Proj
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
from Ellipsoid import Ellipsoid


def scale_factor(lat: float, lon: float, central_lon: float, ellipsoid_a: float, ellipsoid_b: float) -> float:
    """
    Compute the scale factor (distortion parameter) of a map projection at a given point.

    Parameters
    ----------
    lat : float. Geodetic latitude in radians.
    lon : float. Geodetic longitude in radians.
    central_lon : float. Central meridian of the projection in radians.
    ellipsoid_a : float. Semi-major axis of the ellipsoid (meters).
    ellipsoid_b : float. Semi-minor axis of the ellipsoid (meters).

    Returns
    -------
    k : float. The scale factor at the given point.

    Examples
    --------
    >>> a = 6378137.0  # Semi-major axis of WGS84 (meters)
    >>> b = 6356752.314245  # Semi-minor axis of WGS84 (meters)
    >>> lat, lon = np.radians(59.907072474276958), np.radians(10.754482924017791)
    >>> central_lon = np.radians(10.0)
    >>> k = scale_factor(lat, lon, central_lon, a, b)
    >>> print(f"Scale factor: {k:.10f}")
    """
    # Compute eccentricity squared
    e2 = (ellipsoid_a**2 - ellipsoid_b**2) / ellipsoid_a**2

    # Radius of curvature in the meridian (M) and prime vertical (N)
    N = ellipsoid_a / np.sqrt(1 - e2 * np.sin(lat)**2)

    # Projection coordinates (basic transverse Mercator equations)
    delta_lon = lon - central_lon
    t = np.tan(lat)
    eta2 = (ellipsoid_a**2 / ellipsoid_b**2 - 1) * np.cos(lat)**2

    # Scaling distortion factor (simplified for transverse Mercator projection)
    k = np.sqrt(1 + eta2 * np.cos(delta_lon)**2)

    return k



def tm_point_scale_factor_geog(lon: float, lat: float, k0: float = 0.9996) -> float:
    """
    Compute the point scale factor of a Transverse Mercator projection in terms of geographical coordinates.

    Parameters
    ----------
    lon : float
        Longitude in radians.
    lat : float
        Latitude in radians.
    k0 : float, optional
        Scale factor at the central meridian (default is 0.9996).

    Returns
    -------
    float
        Point scale factor at the given geographical coordinates.

    Examples
    --------
    >>> lon = np.radians(10)
    >>> lat = np.radians(50)
    >>> k = tm_point_scale_factor_geog(lon, lat)
    >>> print(k)
    """
    sin_lon = np.sin(lon)
    cos_lat = np.cos(lat)
    k = k0 / np.sqrt(1 - (sin_lon**2 * cos_lat**2))
    return k


def tm_point_scale_factor_proj(x: float, k0: float = 0.9996, a: float = 6378137.0) -> float:
    """
    Compute the point scale factor of a Transverse Mercator projection in terms of projection coordinates.

    Parameters
    ----------
    x : float
        X-coordinate (distance from the central meridian) in meters.
    k0 : float, optional
        Scale factor at the central meridian (default is 0.9996).
    a : float, optional
        Semi-major axis of the ellipsoid in meters (default is WGS84: 6378137.0).

    Returns
    -------
    float
        Point scale factor at the given projection coordinates.

    Examples
    --------
    >>> x = 180000  # 180 km from the central meridian
    >>> k = tm_point_scale_factor_proj(x)
    >>> print(k)
    """
    k = k0 * np.cosh(x / (k0 * a))
    return k



if __name__ == "__main__":

    from Ellipsoid import WGS84
    ellip = WGS84()
    a = ellip.a
    b = ellip.b
    lat = np.radians(59.907072474276958)  # Latitude in radians
    lon = np.radians(10.754482924017791)  # Longitude in radians
    central_lon = np.radians(9.0)  # Central meridian in radians

    # Compute scale factor using the custom function
    custom_k = scale_factor(lat, lon, central_lon, a, b)
    print(f"Custom scale factor: {custom_k:.10f}")

    # Compute scale factor using pyproj
    # projection = Proj(proj="tmerc", lon_0=np.degrees(central_lon), ellps="WGS84")
    projection = Proj(proj='utm',zone=32,ellps='WGS84', preserve_units=False)
    # _, _, pyproj_k = projection(lon, lat, inverse=False, radians=True, errcheck=True, return_scale=True)
    pyproj_k = projection.get_factors(lon, lat)
    print(f"pyproj scale factor: {pyproj_k}")
    # print(f"pyproj scale factor: {pyproj_k:.10f}")

    # # Compare results
    # assert np.isclose(custom_k, pyproj_k, atol=1e-10), "Scale factors do not match!"
    # print("Scale factors match.")

    lon = np.radians(10.754482924017791)  # 10 degrees longitude
    lat = np.radians(59.907072474276958)    # 50 degrees latitude
    k0 = 0.9996  # Central meridian scale factor
    x = 180000  # Distance from the central meridian in meters

    # Point scale factor in geographical coordinates
    k_geog = tm_point_scale_factor_geog(lon, lat, k0)
    print(f"Scale factor (geographical): {k_geog:.6f}")

    # # Point scale factor in projection coordinates
    # k_proj = tm_point_scale_factor_proj(x, k0, a)
    # print(f"Scale factor (projection): {k_proj:.6f}")
