"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""


import numpy as np
from typing import Union, Tuple

def ECEF2enu(lat: Union[float, np.ndarray], lon: Union[float, np.ndarray],
             dX: Union[float, np.ndarray], dY: Union[float, np.ndarray],
             dZ: Union[float, np.ndarray]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Converts a displacement vector in ECEF coordinates to the local topocentric ENU (East-North-Up) coordinate system.

    Parameters
    ----------
    lat : float or np.ndarray
        Geodetic latitude(s) of the reference point(s) in radians.
    lon : float or np.ndarray
        Geodetic longitude(s) of the reference point(s) in radians.
    dX : float or np.ndarray
        X displacement(s) in ECEF coordinates (meters).
    dY : float or np.ndarray
        Y displacement(s) in ECEF coordinates (meters).
    dZ : float or np.ndarray
        Z displacement(s) in ECEF coordinates (meters).

    Returns
    -------
    e : np.ndarray
        East displacement(s) in the ENU coordinate system (meters).
    n : np.ndarray
        North displacement(s) in the ENU coordinate system (meters).
    u : np.ndarray
        Up displacement(s) in the ENU coordinate system (meters).

    Examples
    --------
    >>> lat = np.radians([45.0, 46.0])
    >>> lon = np.radians([9.0, 10.0])
    >>> dX = np.array([100.0, 200.0])
    >>> dY = np.array([200.0, 300.0])
    >>> dZ = np.array([50.0, 75.0])
    >>> e, n, u = ECEF2enu(lat, lon, dX, dY, dZ)
    >>> print(e, n, u)
    """
    # Ensure inputs are numpy arrays
    lat = np.asarray(lat)
    lon = np.asarray(lon)
    dX = np.asarray(dX)
    dY = np.asarray(dY)
    dZ = np.asarray(dZ)

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

    return e, n, u




if __name__ == "__main__":

    e_true = -17553188.505429048
    n_true = -3126586.21287199
    u_true = 13814727.418270092

    # Define reference point geodetic coordinates (latitude and longitude in degrees)
    lat_rad = np.array([1.0455756599070236])
    lon_rad =np.array([0.18770113637361763])
    dX = np.array([12738186.63827794])
    dY = np.array([-15447555.301322976])
    dZ = np.array([10385003.518329535])

    # Perform ECEF to ENU conversion
    e, n, u = ECEF2enu(lat_rad, lon_rad, dX, dY, dZ)

    # Print results
    print(f"ECEF to ENU results:")
    for i in range(len(lat_rad)):
        print(f"Point {i+1}:")
        print(f"  East displacement (e): {e[i]:.6f} meters")
        print(f"  North displacement (n): {n[i]:.6f} meters")
        print(f"  Up displacement (u): {u[i]:.6f} meters")
