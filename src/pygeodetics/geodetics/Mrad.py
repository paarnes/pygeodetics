"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com

"""

import numpy as np

def Mrad(a: float, b: float, lat: float, radians: bool = False) -> float:
    """
    Compute the meridional radius of curvature (M) at a given latitude.

    The meridional radius of curvature (M) represents the radius of curvature
    in the north-south direction along a meridian.

    Parameters
    ----------
    a : float. Semi-major axis of the ellipsoid (meters).
    b : float. Semi-minor axis of the ellipsoid (meters).
    lat : float. Geodetic latitude.
    radians : bool, optional. If `True`, input latitude is in radians. If `False`, input latitude is in degrees. Defaults to `False`.

    Returns
    -------
    M : float. Meridional radius of curvature at the given latitude (meters).

    """
    if not radians:
        lat = np.radians(lat)

    # Compute first eccentricity squared
    e2 = (a**2 - b**2) / a**2

    # Compute meridional radius of curvature (M)
    M = a * (1 - e2) / (1 - e2 * np.sin(lat)**2)**(3/2)

    return M




if __name__ == '__main__':
    import sys, os
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
    from Ellipsoid import WGS84
    ellip = WGS84()
    a, b = ellip.a, ellip.b
    lat = 60
    mrad = Mrad(a, b, lat)
    print(f"Mrad: {mrad:.12f} meters")
