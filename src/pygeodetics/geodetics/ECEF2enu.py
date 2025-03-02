"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""


import numpy as np
from typing import Union, Tuple
from geodetics.Nrad import Nrad
from geodetics.geod2ECEF import geod2ECEF
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
from Ellipsoid import Ellipsoid, WGS84


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
    >>> X = np.array([100.0, 200.0])
    >>> Y = np.array([200.0, 300.0])
    >>> Z = np.array([50.0, 75.0])
    >>> e, n, u = ECEF2enu(lat, lon, X, Y, Z, radians=False)  # Convert degrees to radians
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

    e_true = -17553188.505429048
    n_true = -3126586.21287199
    u_true = 13814727.418270092

    # Define reference point geodetic coordinates (latitude and longitude in degrees)
    lat_rad = np.array([1.0455756599070236]) # 59.90707247427837
    lon_rad =np.array([0.18770113637361763]) # 10.754482924017791
    X = np.array([12738186.63827794])
    Y = np.array([-15447555.301322976])
    Z = np.array([10385003.518329535])
    h0 = np.array([0])
    # Perform ECEF to ENU conversion
    e, n, u = ECEF2enu(X, Y, Z, lat_rad, lon_rad, h0, radians=True)

    # Print results
    print(f"ECEF to ENU results:")
    for i in range(len(lat_rad)):
        print(f"Point {i+1}:")
        print(f"  East displacement (e): {e[i]:.6f} meters")
        print(f"  North displacement (n): {n[i]:.6f} meters")
        print(f"  Up displacement (u): {u[i]:.6f} meters")



    # Define reference point geodetic coordinates (latitude and longitude in degrees)
    lat = np.degrees(1.0455756599070236)
    lon = np.degrees(0.18770113637361763)
    X = 12738186.63827794
    Y = -15447555.301322976
    Z = 10385003.518329535
    h0 = 0

    # Perform ECEF to ENU conversion
    e, n, u = ECEF2enu(X, Y, Z, lat, lon, h0, radians=False)

    # Print results
    print(f"ECEF to ENU results:")
    print(f"  East displacement (e): {e:.6f} meters")
    print(f"  North displacement (n): {n:.6f} meters")
    print(f"  Up displacement (u): {u:.6f} meters")


    # ECEF2enu2
    h = 0
    e, n, u =ECEF2enu_dvec(X, Y, Z, lat, lon, radians=False)
    # Print results
    print(f"ECEF to ENU results (INTERNETT):")
    print(f"  East displacement (e): {e} meters")
    print(f"  North displacement (n): {n} meters")
    print(f"  Up displacement (u): {u} meters")

    # ECEF2enu PYMAP
    from pymap3d import ecef2enu, ecef2enuv, Ellipsoid as ELL, geodetic2ecef
    h = 0
    e, n, u =ECEF2enu_dvec(X, Y, Z, lat, lon, radians=False)
    print(f"ECEF to ENU results (PYMAP):")
    print(f"  East displacement (e): {e} meters")
    print(f"  North displacement (n): {n} meters")
    print(f"  Up displacement (u): {u} meters")

    ev, nv,uv = ecef2enuv(X, Y, Z, lat, lon, deg=True)
    print(f"ECEF to ENU results (PYMAP ENUV):")
    print(f"  East displacement (e): {ev} meters")
    print(f"  North displacement (n): {nv} meters")
    print(f"  Up displacement (u): {uv} meters")


    # ECEF2enu make fasit
    import sys
    from pprint import pprint
    import os
    from geod2ECEF import geod2ECEF
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
    from Ellipsoid import Ellipsoid, WGS84



    lat0=59.907072474276958
    lon0=10.754482924017791
    h0 = 63.8281
    X, Y, Z = 3149785.9652, 598260.8822, 5495348.4927

    lat0=60
    lon0=5
    h0 = 0
    X,Y,Z =12738186.63827794, -15447555.301322976, 10385003.518329535


    lat0=-49.907072474276958
    lon0=-115.754482924017791
    h0 = 10
    X,Y,Z =-2949198.7309,-3422172.1149, -4487393.5421



    a_wgs84 = WGS84().a
    b_wgs84 = WGS84().b
    ell = ELL(a_wgs84, b_wgs84)

    eceff = geodetic2ecef(lat0, lon0, h0, ell=ell, deg=True)
    min_ecef = geod2ECEF(np.radians(lat0), np.radians(lon0), h0)

    # e, n, u = ecef2enu(X, Y, Z, lat0, lon0, h0, ell=ell, deg= True )
    ee, nn, uu = ECEF2enu(X, Y, Z, lat0, lon0, h0, radians=False)
    ei,ni,ui = ecef_to_enu_internett(X, Y, Z, lat0, lon0, h0)
    e, n, u = ecef2enu(X, Y, Z, lat0, lon0, h0 )

    print(f"\n\nECEF to ENU results (MIN):"
        f"\nE: {ee}\nN: {nn}\nU: {uu}")

    print(f"\n\nECEF to ENU results (INTERNETT):"
        f"\nE: {ei}\nN: {ni}\nU: {ui}")

    print(f"\n\nECEF to ENU results (PYMAP ENUV):"
          f"\nE: {e}\nN: {n}\nU: {u}")



