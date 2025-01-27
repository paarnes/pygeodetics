import numpy as np

def geod2ECEF(a: float, b: float, lat: float, lon: float, h: float) -> tuple:
    """
    Convert geodetic coordinates (latitude, longitude, height) to ECEF (Earth-Centered, Earth-Fixed)
    Cartesian coordinates.

    Parameters
    ----------
    a : float
        Semi-major axis of the ellipsoid (meters).
    b : float
        Semi-minor axis of the ellipsoid (meters).
    lat : float
        Geodetic latitude in radians.
    lon : float
        Geodetic longitude in radians.
    h : float
        Height above the ellipsoid in meters.

    Returns
    -------
    tuple
        A tuple containing the ECEF coordinates (X, Y, Z) in meters:
        - X : float. X-coordinate in meters.
        - Y : float. Y-coordinate in meters.
        - Z : float. Z-coordinate in meters.

    Examples
    --------
    >>> a = 6378137.0  # Semi-major axis of WGS84 ellipsoid
    >>> b = 6356752.314245  # Semi-minor axis of WGS84 ellipsoid
    >>> lat = np.radians(59.907072474276958)  # Latitude in radians
    >>> lon = np.radians(10.754482924017791)  # Longitude in radians
    >>> h = 63.8281  # Height in meters
    >>> X, Y, Z = geod2ECEF(a, b, lat, lon, h)
    >>> print(f"X: {X:.3f}, Y: {Y:.3f}, Z: {Z:.3f}")
    """
    # Compute the prime vertical radius of curvature
    e2 = (a**2 - b**2) / a**2  # Square of the first eccentricity
    N = a / np.sqrt(1 - e2 * np.sin(lat)**2)

    # Calculate ECEF coordinates
    X = (N + h) * np.cos(lat) * np.cos(lon)
    Y = (N + h) * np.cos(lat) * np.sin(lon)
    Z = ((b**2 / a**2) * N + h) * np.sin(lat)

    return X, Y, Z


if __name__ == "__main__":
    # Example usage
    a = 6378137.0  # Semi-major axis of WGS84 ellipsoid
    b = 6356752.314245  # Semi-minor axis of WGS84 ellipsoid
    lat = np.radians(59.907072474276958)  # Latitude in radians
    lon = np.radians(10.754482924017791)  # Longitude in radians
    h = 63.8281  # Height in meters
    # Should get X, Y, Z = 3149785.9652, 598260.8822, 5495348.4927

    X, Y, Z = geod2ECEF(a, b, lat, lon, h)
    print(f"ECEF coordinates:\nX: {X:.4f} m\nY: {Y:.4f} m\nZ: {Z:.4f} m")
