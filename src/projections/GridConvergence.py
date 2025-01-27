"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com

Assumes earth is a perfect sphere.

Based on the "THE MERCATOR PROJECTIONS" book from  Peter Osborne, 2013
See "The scale factor for the TMS projection" section at page 63.

"""


from typing import Literal
import numpy as np


def tms_grid_convergence(x: float, y: float, false_easting: float, R: float=6371000, angle_unit: Literal["deg", "rad"]="deg") -> float:
    """
    Compute the grid convergence of a Transverse Mercator projection on a sphere (TMS)
    by using projection coordinates.

    Note: Assumes earth is a perfect sphere.

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
    lon = np.radians(10)  # Longitude in radians
    lat = np.radians(60)  # Latitude in radians
    x, y = Transformer.from_crs("EPSG:4326", "EPSG:32632", always_xy=True).transform(lon, lat, radians=True)
    false_esting=500000 # false easting in meters


    # Grid convergence in projection coordinates
    gamma_proj = tms_grid_convergence(x, y, false_esting)
    print(f"Grid convergence (projection): {gamma_proj:.6f} degrees\n")

    # projection = Proj(proj="tmerc", lon_0=9, ellps="WGS84")
    projection = Proj(proj='utm',zone=32,ellps='WGS84')
    pyproj_k = projection.get_factors(longitude=lon, latitude=lat, radians=True)
    print(f"Pyproj projection factors: \n - Meridian convergence: {pyproj_k.meridian_convergence}")

    # print(f"pyproj scale factor: {np.degrees(pyproj_k)}")
