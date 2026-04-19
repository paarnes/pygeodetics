"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""


import numpy as np
from typing import Union, Tuple
from .Nrad import Nrad
from .geod2ECEF import geod2ECEF
from ..Ellipsoid import Ellipsoid, WGS84


def ECEF2enu_dvec(
    dX: Union[float, np.ndarray],
    dY: Union[float, np.ndarray],
    dZ: Union[float, np.ndarray],
    lat: Union[float, np.ndarray],
    lon: Union[float, np.ndarray],
    radians: bool = False
    ) -> Tuple[Union[np.ndarray, float], Union[np.ndarray, float], Union[np.ndarray, float]]:
    """
    Converts ECEF displacement vectors to ENU displacement vectors which
    is a local topocentric ENU (East-North-Up) coordinate system. It computes the
    displacement (ΔX, ΔY, ΔZ) relative to the observer.

    Parameters
    ----------
    lat : float or np.ndarray. Geodetic latitude(s) of the reference point(s).
    lon : float or np.ndarray. Geodetic longitude(s) of the reference point(s).
    dX : float or np.ndarray. Difference (∆X) in X-coordinate (ECEF in meters).
    dY : float or np.ndarray. Difference (∆Y) in X-coordinate(ECEF in meters).
    dZ : float or np.ndarray. Difference (∆Z) in X-coordinate (ECEF in meters).
    radians : bool, optional. If `False`, assumes `lat` and `lon` are in degrees and converts them to radians.
              Defaults to `True` (assumes input is already in radians).

    Returns
    -------
    e : float|np.ndarray. East coordinate in the ENU coordinate system (meters).
    n : float|np.ndarray. North coordinate in the ENU coordinate system (meters).
    u : float|np.ndarray. Up coordinate in the ENU coordinate system (meters).

    Examples
    --------
    >>> lat = [45.0, 46.0]
    >>> lon = [9.0, 10.0]
    >>> dX = np.array([100.0, 200.0])
    >>> dY = np.array([200.0, 300.0])
    >>> dZ = np.array([50.0, 75.0])
    >>> e, n, u = ECEF2enu_dvec(dX, dY, dZ, lat, lon, radians=False)
    >>> print(e, n, u)
    """
    # Ensure inputs are numpy arrays and at least 1D to handle scalars
    is_float = isinstance(lat, float)

    lat = np.atleast_1d(lat)
    lon = np.atleast_1d(lon)
    dX = np.atleast_1d(dX)
    dY = np.atleast_1d(dY)
    dZ = np.atleast_1d(dZ)

    # Convert lat/lon to radians if needed
    if not radians:
        lat = np.radians(lat)
        lon = np.radians(lon)

    # ECEF displacement vector
    dP_ECEF = np.stack((dX, dY, dZ), axis=-1)  # Shape (N, 3) for multiple inputs

    # Create rotation matrix for ENU transformation
    M = np.array([
        [-np.sin(lon), np.cos(lon), np.zeros_like(lon)],
        [-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)],
        [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)],
    ]).transpose(2, 0, 1)  # Shape (N, 3, 3)

    # Perform matrix multiplication (vectorized)
    dP_ENU = np.einsum('nij,nj->ni', M, dP_ECEF)  # Shape (N, 3)

    # Extract East, North, Up components
    e = dP_ENU[:, 0]  # East
    n = dP_ENU[:, 1]  # North
    u = dP_ENU[:, 2]  # Up

    if is_float:
        return e.item(), n.item(), u.item()  # Ensures single float output

    return e, n, u



def ECEF2enu(
    X: Union[float, np.ndarray],
    Y: Union[float, np.ndarray],
    Z: Union[float, np.ndarray],
    lat0: Union[float, np.ndarray],
    lon0: Union[float, np.ndarray],
    h0: Union[float, np.ndarray],
    ellipsoid: Ellipsoid = WGS84(),
    radians: bool = False) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Converts absolute ECEF (Earth-Centered, Earth-Fixed) coordinates to the
    local topocentric ENU (East-North-Up) coordinate system.

    Parameters
    ----------
    X, Y, Z : float or np.ndarray. ECEF coordinates of the target point(s) in meters.
    lat0, lon0 : float or np.ndarray.Geodetic latitude and longitude of the reference point(s).
    h0 : float or np.ndarray. Altitude of the reference point(s) in meters.
    radians : bool, optional. If False (default), assumes `lat0` and `lon0` are in degrees and converts them to radians.

    Returns
    -------
    e, n, u : np.ndarray. ENU (East, North, Up) coordinates in meters.

    Examples
    --------
    >>> X, Y, Z = 12738186.63827794, -15447555.301322976, 10385003.518329535
    >>> lat0, lon0, h0 = 60.0, 10.0, 100.0
    >>> e, n, u = ecef_to_enu(X, Y, Z, lat0, lon0, h0, radians=False)
    >>> print(e, n, u)
    """

    # Ensure inputs are numpy arrays
    is_float = isinstance(X, float)
    X, Y, Z = np.atleast_1d(X), np.atleast_1d(Y), np.atleast_1d(Z)
    lat0, lon0, h0 = np.atleast_1d(lat0), np.atleast_1d(lon0), np.atleast_1d(h0)

    # Convert lat/lon to radians if needed
    if not radians:
        lat0 = np.radians(lat0)
        lon0 = np.radians(lon0)

    # Precompute trigonometric values
    sin_lat = np.sin(lat0)
    cos_lat = np.cos(lat0)
    sin_lon = np.sin(lon0)
    cos_lon = np.cos(lon0)

    # Compute ECEF coordinates of the reference point (lat0, lon0, h0)
    x0, y0, z0 = geod2ECEF(lat0, lon0, h0, ellipsoid=ellipsoid, radians=True)

    # Compute relative ECEF coordinates
    xd, yd, zd = X - x0, Y - y0, Z - z0

    # Compute ENU coordinates using rotation matrix
    e = -sin_lon * xd + cos_lon * yd
    n = -cos_lon * sin_lat * xd - sin_lat * sin_lon * yd + cos_lat * zd
    u = cos_lat * cos_lon * xd + cos_lat * sin_lon * yd + sin_lat * zd

    # If inputs were scalars, return scalars instead of arrays
    if is_float:
        return e.item(), n.item(), u.item()  # Ensures single float output

    return e, n, u




if __name__ == "__main__":
    # Simple demo
    lat0 = 59.907072474276958
    lon0 = 10.754482924017791
    h0 = 63.8281
    X, Y, Z = 3149785.9652, 598260.8822, 5495348.4927

    e, n, u = ECEF2enu(X, Y, Z, lat0, lon0, h0, radians=False)
    print("ECEF to ENU results:")
    print(f"  East : {e:.6f} m")
    print(f"  North: {n:.6f} m")
    print(f"  Up   : {u:.6f} m")
