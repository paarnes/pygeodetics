import numpy as np

def tm_grid_convergence_geog(lon: float, lat: float) -> float:
    """
    Compute the grid convergence of a Transverse Mercator projection in terms of geographical coordinates.

    Parameters
    ----------
    lon : float
        Longitude in radians.
    lat : float
        Latitude in radians.

    Returns
    -------
    float
        Grid convergence angle (gamma) in radians.

    Examples
    --------
    >>> lon = np.radians(10)
    >>> lat = np.radians(50)
    >>> gamma = tm_grid_convergence_geog(lon, lat)
    >>> print(np.degrees(gamma))  # Convert to degrees
    """
    gamma = np.arctan(np.tan(lon) * np.sin(lat))
    return gamma


def tm_grid_convergence_proj(x: float, y: float, k0: float = 0.9996, a: float = 6378137.0) -> float:
    """
    Compute the grid convergence of a Transverse Mercator projection in terms of projection coordinates.

    Parameters
    ----------
    x : float
        X-coordinate (distance from the central meridian) in meters.
    y : float
        Y-coordinate (distance from the equator) in meters.
    k0 : float, optional
        Scale factor at the central meridian (default is 0.9996).
    a : float, optional
        Semi-major axis of the ellipsoid in meters (default is WGS84: 6378137.0).

    Returns
    -------
    float
        Grid convergence angle (gamma) in radians.

    Examples
    --------
    >>> x = 180000  # 180 km from the central meridian
    >>> y = 50000  # 50 km north of the equator
    >>> gamma = tm_grid_convergence_proj(x, y)
    >>> print(np.degrees(gamma))  # Convert to degrees
    """
    gamma = np.arctan(np.tanh(x / (k0 * a)) * np.tan(y / (k0 * a)))
    return gamma


if __name__ == "__main__":
    from pyproj import Proj
    lon = np.radians(10)  # Longitude in radians
    lat = np.radians(60)  # Latitude in radians
    k0 = 0.9996  # Central meridian scale factor
    a = 6378137.0  # Semi-major axis of WGS84 ellipsoid in meters
    x = 180000  # Distance from the central meridian in meters
    y = 50000  # Distance from the equator in meters

    # Grid convergence in geographical coordinates
    gamma_geog = tm_grid_convergence_geog(lon, lat)
    print(f"Grid convergence (geographical): {np.degrees(gamma_geog):.6f} degrees")

    # Grid convergence in projection coordinates
    gamma_proj = tm_grid_convergence_proj(x, y, k0, a)
    print(f"Grid convergence (projection): {np.degrees(gamma_proj):.6f} degrees")


    projection = Proj(proj='utm',zone=32,ellps='WGS84', preserve_units=False)
    # _, _, pyproj_k = projection(lon, lat, inverse=False, radians=True, errcheck=True, return_scale=True)
    pyproj_k = projection.get_factors(lon, lat)
    print(f"pyproj scale factor: {pyproj_k}")
