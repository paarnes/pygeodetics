"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com

Assumes earth is a perfect sphere.

Based on the "THE MERCATOR PROJECTIONS" book from  Peter Osborne, 2013
See "The scale factor for the TMS projection" section at page 61.

"""

import numpy as np


def tms_point_scale_factor(x: float, false_esting: float, R: float = 6371000.0) -> float:
    """
    Compute the point scale factor of a Transverse Mercator projection for a sphere (TMS) using projected coordinates.
    Note: Assumes earth is a perfect sphere.

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
    lat = np.radians(59.907072474276958)  # Latitude in radians
    lon = np.radians(10.754482924017791)  # Longitude in radians
    central_lon = np.radians(9.0)  # Central meridian in radians
    x, _ = Transformer.from_crs("EPSG:4326", "EPSG:32632", always_xy=True).transform(lon, lat, radians=True)
    false_easting=500000.0


    # Compute scale factor using the custom function
    custom_k = tms_point_scale_factor(x, false_esting=false_easting)
    print(f"Custom scale factor projected: {custom_k:.10f}")

    # Compute scale factor using pyproj
    projection = Proj(proj="tmerc", lon_0=np.degrees(central_lon), ellps="WGS84")
    # projection = Proj(proj='utm',zone=32,ellps='WGS84')
    pyproj_k = projection.get_factors(longitude=lon, latitude=lat, radians=True)
    print(f"pyproj scale factor: {pyproj_k.meridional_scale:.10f}")

