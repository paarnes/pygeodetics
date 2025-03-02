"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com

"""

import sys
import os
from typing import Literal
import numpy as np
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
from geodetics.footpoint_latitude import footpoint_latitude
from geodetics.Nrad import Nrad



def tm_grid_convergence_geographic(a: float, b: float, lon: float, lat: float,
                                   lon0: float, radians: bool=False) -> float:
    """
    Compute the Transverse Mercator meridian convergence (grid convergence) based on geographic coordinates (latitude, longitude).

    Notes
    -----
    The meridian convergence (gamma) is the angular difference between grid north and true north in the Transverse Mercator projection.

    Parameters
    ----------
    a : float. Semi-major axis of the ellipsoid (meters).
    b : float. Semi-minor axis of the ellipsoid (meters).
    lon : float. Geodetic longitude in radians.
    lat : float. Geodetic latitude in radians.
    lon0 : float. Central meridian of the projection in radians. The tangent of the Meridian (Ex 9° for UTM Z32N)
    radians : bool, optional. If `True`, input and output are in radians. If False, input is in degrees and output is converted to degrees. Defaults to `False`.


    Returns
    -------
    gamma : float. Meridian convergence in degrees (if radians=False) or radians (if radians=True).

    """
    if not radians:
        lon, lat, lon0 = np.radians([lon, lat, lon0])

    dlon = lon - lon0  # Difference in longitude from the central meridian
    e2 = (a**2 - b**2) / a**2
    eps2 = e2 / (1 - e2) * np.cos(lat)**2

    gamma = (
        dlon * np.sin(lat)
        + (dlon**3 / 3) * np.sin(lat) * np.cos(lat)**2 * (1 + 3 * eps2 + 2 * eps2**2)
        + (dlon**5 / 15) * np.sin(lat) * np.cos(lat)**4 * (2 - np.tan(lat)**2)
    )

    if not radians:
        gamma = np.degrees(gamma)

    return gamma


def tm_grid_convergence_projected(a: float, b: float, x: float, y: float, lat0: float, false_easting: float, radians: bool = False) -> float:
    """
    Compute the Transverse Mercator meridian convergence (grid convergence) based on projected coordinates (x, y).

    Notes
    -----
    The meridian convergence (gamma) is the angular difference between grid north and true north.
    It is computed using projected coordinates (eastings and northings) rather than geographic coordinates.

    Parameters
    ----------
    a : float. Semi-major axis of the ellipsoid (meters).
    b : float. Semi-minor axis of the ellipsoid (meters).
    x : float. Easting coordinate (meters).
    y : float. Northing coordinate (meters).
    lat0 : float. Latitude of natural origin (degrees if radians=False, radians if radians=True).
    false_easting : float. False easting value (meters) to adjust the easting coordinate.
                    UTM uses 500 000 meters here. Set to 0 (zero) if not applicable.
    radians : bool, optional. If `True`, input latitude is in radians. If `False`, input latitude is in degrees. Defaults to `False`.

    Returns
    -------
    gamma : float. Meridian convergence in degrees (if radians=False) or radians (if radians=True).
    """
    if not radians:
        lat0 = np.radians(lat0)  # Convert latitude of origin to radians

    # Adjust x-coordinate by subtracting false easting
    x -= false_easting

    # Compute footpoint latitude
    latf = footpoint_latitude(a, b, y, lat0, radians=True)  # Footpoint latitude in radians

    # Compute first eccentricity squared
    e2 = (a**2 - b**2) / a**2

    # Compute normal radius of curvature at footpoint latitude
    Nf = Nrad(a, b, latf, radians=True)

    # Compute second eccentricity squared at footpoint latitude
    epsf2 = e2 / (1 - e2) * np.cos(latf)**2

    # Compute meridian convergence (gamma)
    gamma = (
        (x * np.tan(latf) / Nf) -
        (x**3 * np.tan(latf) / (3 * Nf**3) * (1 + np.tan(latf)**2 - epsf2 - 2 * epsf2**2))
    )

    if not radians:
        gamma = np.degrees(gamma)

    return gamma




def tm_sphere_grid_conv_projected(x: float, y: float, false_easting: float, R: float=6371000, angle_unit: Literal["deg", "rad"]="deg") -> float:
    """
    Compute the grid convergence of a Transverse Mercator projection on a sphere (TMS)
    by using projection coordinates.

    Note: Assumes earth is a perfect sphere. Based on the "THE MERCATOR PROJECTIONS" book from  Peter Osborne, 2013
    See "The scale factor for the TMS projection" section at page 63.

    Parameters
    ----------
    x : float. X-coordinate (distance from the central meridian) in meters.
    y : float. Y-coordinate (distance from the equator) in meters.
    false_easting: float. The false easting value in meters that is used by the projection
    R : float, optional. Radius of the sphere in meters. Default is 6371000 meters.
    angle_unit : str, optional. Unit of the grid convergence angle. Default is degrees

    Returns
    -------
    float. Grid convergence angle (gamma) in radians.

    """
    x = x - false_easting
    gamma = np.arctan(np.tanh(x / (R)) * np.tan(y / (R)))
    if angle_unit == "deg":
        gamma = np.degrees(gamma)
    return gamma





if __name__ == "__main__":
    from pyproj import Proj, Transformer
    import sys, os
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
    from Ellipsoid import WGS84

    ellip = WGS84()
    a = ellip.a
    b = ellip.b
    lon = np.radians(10)  # Longitude in radians
    lat = np.radians(60)  # Latitude in radians
    x, y = Transformer.from_crs("EPSG:4326", "EPSG:32632", always_xy=True).transform(lon, lat, radians=True)
    false_easting=500000 # false easting in meters
    lon0 = np.deg2rad(9) # central meridian/tangent meridian
    lat0 = np.deg2rad(0)


    # Grid convergence in projection coordinates
    gamma_proj = tm_sphere_grid_conv_projected(x, y, false_easting)
    print(f"Grid convergence (projection): {gamma_proj:.6f} degrees\n")

    gamma_geog = tm_grid_convergence_geographic(a, b, lon, lat,lon0, radians=True)
    print(f"Grid convergence (geographic): {np.rad2deg(gamma_geog):.12f} degrees\n")

    # gamma_proj = tm_grid_convergence_projected(a, b, x, y, lat0, false_easting=false_easting, radians=True)
    gamma_proj = tm_grid_convergence_projected(a, b, x, y, lat0, false_easting=false_easting, radians=True)
    print(f"Grid convergence (projected NEW): {np.rad2deg(gamma_proj):.12f} degrees\n")

    # projection = Proj(proj="tmerc", lon_0=9, ellps="WGS84")
    projection = Proj(proj='utm',zone=32,ellps='WGS84')
    pyproj_k = projection.get_factors(longitude=lon, latitude=lat, radians=True)
    print(f"Pyproj projection factors: \n - Meridian convergence: {pyproj_k.meridian_convergence}")

    # print(f"pyproj scale factor: {np.degrees(pyproj_k)}")
