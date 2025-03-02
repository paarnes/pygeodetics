"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com

"""

import numpy as np

def meridional_arc_dist(a: float, b: float, lat: float, radians: bool = False) -> float:
    """
    Compute the meridional arc length from the equator to a given latitude.

    Notes
    -----
    The meridional arc distance is the north-south distance along a meridian between the equator and a given latitude.

    Parameters
    ----------
    a : float. Semi-major axis of the ellipsoid (meters).
    b : float. Semi-minor axis of the ellipsoid (meters).
    lat : float. Geodetic latitude in degrees (if radians=False) or in radians (if radians=True).
    radians : bool, optional. If `True`, input latitude is in radians. If `False`, input is in degrees and will be converted. Defaults to `False`.

    Returns
    -------
    B : float. Meridional arc length in meters.
    """
    if not radians:
        lat = np.radians(lat)  # Convert latitude to radians if needed

    # Compute flattening
    f = (a - b) / a

    # Compute base coefficient
    b0 = a * (1 - (1 / 2) * f + (1 / 16) * f**2 + (1 / 32) * f**3)

    # Compute meridional arc distance
    B = (
        b0 * (
            lat
            - (3 / 4 * f + 3 / 8 * f**2 + 15 / 128 * f**3) * np.sin(2 * lat)
            + (15 / 64 * f**2 + 15 / 64 * f**3) * np.sin(4 * lat)
            - (35 / 384 * f**3) * np.sin(6 * lat)
        )
    )

    return B



if __name__ == '__main__':
    import sys, os
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
    from Ellipsoid import WGS84
    # Compute the arc length
    ellip = WGS84()
    a, b = ellip.a, ellip.b
    mer_arc_dist = meridional_arc_dist(a, b, 60)
    print(f"Meridional arc length: {mer_arc_dist:.12f} meters")