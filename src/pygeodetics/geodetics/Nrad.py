"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com

"""


import numpy as np

def Nrad(a: float, b: float, lat: float, radians: bool = False) -> float:
    """
    Compute the normal radius of curvature (N) at a given latitude.

    The normal radius of curvature (N) is perpendicular to the meridian
    and represents the radius of curvature in the east-west direction.

    Parameters
    ----------
    a : float. Semi-major axis of the ellipsoid (meters).
    b : float. Semi-minor axis of the ellipsoid (meters).
    lat : float. Geodetic latitude in degrees (if radians=False) or in radians (if radians=True).
    radians : bool, optional. If `True`, input latitude is in radians. If `False`, input is in degrees. Defaults to `False`.

    Returns
    -------
    N : float. Normal radius of curvature (meters).
    """
    if not radians:
        lat = np.radians(lat)  # Convert latitude to radians if needed

    e2 = (a**2 - b**2) / a**2  # First eccentricity squared
    N = a / np.sqrt(1 - e2 * np.sin(lat)**2)

    return N



if __name__ == '__main__':
    import sys, os
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
    from Ellipsoid import WGS84
    ellip = WGS84()
    a, b = ellip.a, ellip.b
    lat = 60
    nrad = Nrad(a,b,lat)
    print(f"Nrad: {nrad:.12f} meters")