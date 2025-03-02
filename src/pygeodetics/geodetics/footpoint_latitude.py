"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com

"""

import sys
import os
import numpy as np
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
from geodetics.meridional_arc_dist import meridional_arc_dist

def footpoint_latitude(a: float, b: float, y: float, lat0: float, radians: bool = False) -> float:
    """
    Compute the footpoint latitude (latitude of the origin of meridian arc) for a given northing coordinate.

    Notes
    -----
    The footpoint latitude is an intermediate latitude used in the inverse Transverse Mercator projection
    to convert northing (y) values back to geodetic latitude.

    Parameters
    ----------
    a : float. Semi-major axis of the ellipsoid (meters).
    b : float. Semi-minor axis of the ellipsoid (meters).
    y : float. Northing coordinate (meters).
    lat0 : float. Latitude of the origin of meridian arc in degrees (if radians=False) or in radians (if radians=True).
    radians : bool, optional. If `True`, input and output are in radians. If `False`, input is in degrees and output is converted to degrees. Defaults to `False`.

    Returns
    -------
    latf : float. Footpoint latitude in degrees (if radians=False) or radians (if radians=True).
    """
    if not radians:
        lat0 = np.radians(lat0)  # Convert latitude of origin to radians if needed

    # Compute flattening
    f = (a - b) / a

    # Compute base coefficient
    b0 = a * (1 - 0.5 * f + (1 / 16) * f**2 + (1 / 32) * f**3)

    # Compute meridional arc length
    B = meridional_arc_dist(a, b, lat0, radians=True) + y  # Using fixed function for meridional arc distance

    # Compute footpoint latitude
    latf = (
        B / b0
        + (3 / 4 * f + 3 / 8 * f**2 + 21 / 256 * f**3) * np.sin(2 * B / b0)
        + (21 / 64 * f**2 + 21 / 64 * f**3) * np.sin(4 * B / b0)
        + (151 / 768 * f**3) * np.sin(6 * B / b0)
    )

    if not radians:
        latf = np.degrees(latf)

    return latf




if __name__ == '__main__':
    import sys, os
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
    from Ellipsoid import WGS84
    # Compute the arc length
    ellip = WGS84()
    a, b = ellip.a, ellip.b
    y = 6651832.735433666035533 # y coordinate
    lat0 = 0
    mer_arc_dist = meridional_arc_dist(a, b, 60)
    print(f"Meridional arc length: {mer_arc_dist:.12f} meters")

    fot_lat = footpoint_latitude(a, b, y, lat0)
    print(f"Foot point latitude: {fot_lat:.12f} degrees")