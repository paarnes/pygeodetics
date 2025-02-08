"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""
import sys
import os
import numpy as np
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
from geodetics.Nrad import Nrad
from geodetics.Mrad import Mrad
from geodetics.footpoint_latitude import footpoint_latitude


def tm_point_scale_factor_geographic(a: float, b: float, lon: float, lat: float,
                                     lon0: float, radians: bool = False) -> float:
    """
    Compute the Transverse Mercator point scale factor based on geographic coordinates (latitude, longitude).

    Notes
    -----
    The point scale factor describes how much a unit distance on the ellipsoid is distorted in the projection.
    This is computed based on geographic coordinates rather than projected coordinates.

    Parameters
    ----------
    a : float. Semi-major axis of the ellipsoid (meters).
    b : float. Semi-minor axis of the ellipsoid (meters).
    lon : float. Geodetic longitude.
    lat : float. Geodetic latitude.
    lon0 : float. Central meridian of the projection.
    radians : bool, optional. If `True`, input and output are in radians. If `False`, input is in degrees and output is converted to degrees. Defaults to `False`.

    Returns
    -------
    scale : float. Scale factor at the given geographic location.

    """
    if not radians:
        lat, lon, lon0 = np.radians([lat, lon, lon0])

    # Compute longitude difference
    dlon = lon - lon0

    # Compute first eccentricity squared
    e2 = (a**2 - b**2) / a**2

    # Compute second eccentricity squared at given latitude
    eps2 = e2 / (1 - e2) * np.cos(lat)**2

    # Compute scale factor
    scale = (
        1 + (dlon**2 / 2) * np.cos(lat)**2 * (1 + eps2) +
        (dlon**4 / 24) * np.cos(lat)**4 * (5 + 4 * np.tan(lat)**2)
    )

    return scale


def tm_point_scale_factor_projected(a: float, b: float, x: float, y: float, lat0: float, false_easting: float, radians: bool = False) -> float:
    """
    Compute the Transverse Mercator scale factor based on projected coordinates (x, y).

    Notes
    -----
    The scale factor represents the distortion introduced by the projection.
    It is computed using the projected coordinates (eastings and northings) rather than geographic coordinates.

    Parameters
    ----------
    a : float. Semi-major axis of the ellipsoid (meters).
    b : float. Semi-minor axis of the ellipsoid (meters).
    x : float. Easting coordinate (meters).
    y : float. Northing coordinate (meters).
    lat0 : float. Latitude of natural origin (degrees if radians=False, radians if radians=True).
    false_easting : float. False easting value (meters) to adjust the easting coordinate. UTM uses 500 000 meters here.
    radians : bool, optional. If `True`, input latitude is in radians. If `False`, input latitude is in degrees. Defaults to `False`.

    Returns
    -------
    scale : float. Scale factor at the given projected coordinates.
    """
    if not radians:
        lat0 = np.radians(lat0)

    # Adjust x-coordinate by subtracting false easting
    x -= false_easting

    # Compute footpoint latitude
    latf = footpoint_latitude(a, b, y, lat0, radians=True)

    # Compute normal and meridional radii of curvature
    Nf = Nrad(a, b, latf, radians=True)
    Mf = Mrad(a, b, latf, radians=True)

    # Compute Transverse Mercator scale factor
    scale = 1 + (x**2 / (2 * Mf * Nf)) + (x**4 / (24 * Nf**4))

    return scale


def tm_sphere_point_scale_factor(x: float, false_esting: float, R: float = 6371000.0) -> float:
    """
    Compute the point scale factor of a Transverse Mercator projection for a sphere (TMS) using projected coordinates.
    Note: Assumes earth is a perfect sphere.

    Based on the "THE MERCATOR PROJECTIONS" book from  Peter Osborne, 2013
    See "The scale factor for the TMS projection" section at page 61.

    Parameters
    ----------
    x : float. X-coordinate (distance from the central meridian) in meters.
    false_easting : float, optional. False easting value for the projection in meters.
    R : float, optional. Radius of the sphere in meters.

    Returns
    -------
    float. Point scale factor at the given projection coordinates.
    """
    # Subtract false easting to get distance from the central meridian
    x = x - false_esting
    k = np.cosh(x / (R))
    return k




if __name__ == "__main__":

    from pyproj import Proj, Transformer
    import sys, os
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
    from Ellipsoid import WGS84, GRS80
    # lat = np.radians(59.907072474276958)  # Latitude in radians
    # lon = np.radians(10.754482924017791)  # Longitude in radians
    ellip = WGS84()
    a = ellip.a
    b = ellip.b
    lat = np.radians(60)  # Latitude in radians
    lon = np.radians(10)  # Longitude in radians
    central_lon = np.radians(9.0)  # Central meridian in radians
    x, y = Transformer.from_crs("EPSG:4326", "EPSG:32632", always_xy=True).transform(lon, lat, radians=True)
    false_easting=500000.0
    lon0 = np.deg2rad(9) # UTM Z32N
    lat0 = np.deg2rad(0)


    # Compute scale factor using the custom function
    custom_k_sphere = tm_sphere_point_scale_factor(x, false_esting=false_easting)
    print(f"Custom scale factor projected (sphere): {custom_k_sphere:.12f}")

    custom_k = tm_point_scale_factor_geographic(a, b, lon, lat, lon0, radians=True)
    print(f"Custom scale factor geograpic: {custom_k:.12f}")

    custom_k_proj = tm_point_scale_factor_projected(a, b, x, y, lat0, false_easting, radians=True)
    print(f"Custom scale factor projected: {custom_k_proj:.12f}")

    # Compute scale factor using pyproj
    projection = Proj(proj="tmerc", lon_0=np.degrees(central_lon), ellps="WGS84")
    # projection = Proj(proj='utm',zone=32,ellps='WGS84')
    pyproj_k = projection.get_factors(longitude=lon, latitude=lat, radians=True)
    print(f"pyproj scale factor: {pyproj_k.meridional_scale:.12f}")

