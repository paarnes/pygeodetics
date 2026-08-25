"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""
from typing import Literal, Union
import numpy as np
from ._kruger import tm_factors_geographic, tm_factors_projected
from ..geodetics.Nrad import Nrad
from ..geodetics.Mrad import Mrad
from ..geodetics.footpoint_latitude import footpoint_latitude

Method = Literal["kruger", "snyder"]


def _check_method(method: str) -> None:
    if method not in ("kruger", "snyder"):
        raise ValueError(f"method must be 'kruger' or 'snyder', got {method!r}")


def tm_point_scale_factor_geographic(a: float, b: float, lon, lat,
                                     lon0: float, *, k0: float = 1.0,
                                     radians: bool = False,
                                     method: Method = "kruger") -> Union[float, np.ndarray]:
    """
    Compute the Transverse Mercator point scale factor based on geographic coordinates (latitude, longitude).

    Notes
    -----
    The point scale factor describes how much a unit distance on the ellipsoid is distorted in the projection.
    This is computed based on geographic coordinates rather than projected coordinates.
    The scale factor is dimensionless, so the `radians` flag only controls the unit of the
    input angles - it has no effect on the returned value.

    The series used is the standard one from Snyder (1987), *Map Projections - A Working
    Manual*, Transverse Mercator, with ``A = dlon * cos(lat)``, ``T = tan(lat)**2`` and
    ``C = e'**2 * cos(lat)**2``::

        k / k0 = 1 + (1 + C) * A**2 / 2
                   + (5 - 4T + 42C + 13C**2 - 28e'**2) * A**4 / 24
                   + (61 - 148T + 16T**2) * A**6 / 720

    That series is selected with ``method="snyder"``. The default ``method="kruger"`` uses
    the expansion in the third flattening (Karney 2011) instead, which is exact to double
    precision. Both are multiplied by `k0`, so pass ``k0=0.9996`` for UTM. Forgetting it is
    a silent 400 ppm error (40 cm per kilometre).

    ``lon`` and ``lon0`` must be referenced to the **same prime meridian**. Mixing frames is a
    silent error that still returns a plausible-looking value. This is a real trap for grids
    whose EPSG definition uses a non-Greenwich prime meridian, e.g. the Norwegian NGO 1948
    axes: EPSG:27391 (NGO zone I) lists a longitude of natural origin of ``4°40'00"W``, which
    is ``-4.666667°`` **relative to Oslo**, not to Greenwich. Either convert the central
    meridian to Greenwich, or convert the point longitude to the local frame; both of the
    following are equivalent::

        # Greenwich throughout: Oslo (10°43'22.5"E) - 4°40' = 6°03'22.5"E
        tm_point_scale_factor_geographic(a, b, lon_greenwich, lat, lon0=6.05625)

        # Oslo-relative throughout:
        tm_point_scale_factor_geographic(a, b, lon_greenwich - 10.7229166667, lat, lon0=-4.6666667)

    Accuracy: with ``method="kruger"`` the error against a rigorous solution (pyproj/PROJ
    ``tmerc``) is below ``1e-10`` everywhere from the equator to the poles and out to
    ``±45°`` from the central meridian. With ``method="snyder"`` the ``A**6`` truncation
    keeps the error below ``1e-10`` within ``±4°`` of the central meridian on GRS80 at
    ``lat = 59°``, rising to about ``6e-10`` at ``±8°``.

    Parameters
    ----------
    a : float. Semi-major axis of the ellipsoid (meters).
    b : float. Semi-minor axis of the ellipsoid (meters).
    lon : float or numpy.ndarray. Geodetic longitude in degrees (if `radians=False`) or radians (if `radians=True`).
    lat : float or numpy.ndarray. Geodetic latitude in degrees (if `radians=False`) or radians (if `radians=True`).
    lon0 : float. Central meridian of the projection, in degrees (if `radians=False`) or radians (if `radians=True`).
    k0 : float, keyword-only, optional. Scale factor at the central meridian. UTM uses `0.9996`.
         Defaults to `1.0` (pure Transverse Mercator).
    radians : bool, keyword-only, optional. If `True`, the input angles are in radians.
              If `False`, they are in degrees. Defaults to `False`. The returned scale factor
              is dimensionless and unaffected by this flag.
    method : {"kruger", "snyder"}, keyword-only, optional. Series to use. `"kruger"` is the
             Kruger expansion in the third flattening (Karney 2011), accurate to double
             precision. `"snyder"` is the classic expansion in `dlon` truncated at the 6th
             order. Defaults to `"kruger"`.

    Returns
    -------
    scale : float or numpy.ndarray. Dimensionless point scale factor at the given geographic location.

    Examples
    --------
    >>> from pygeodetics.Ellipsoid import GRS80
    >>> ellip = GRS80()
    >>> k = tm_point_scale_factor_geographic(ellip.a, ellip.b, 10.0, 59.0, 9.0)
    >>> round(float(k), 9)
    1.000040473

    Same point, but as a UTM zone 32N scale factor:

    >>> k_utm = tm_point_scale_factor_geographic(ellip.a, ellip.b, 10.0, 59.0, 9.0, k0=0.9996)
    >>> round(float(k_utm), 9)
    0.999640456
    """
    _check_method(method)
    if not radians:
        lat, lon, lon0 = np.radians(lat), np.radians(lon), np.radians(lon0)

    # Compute longitude difference
    dlon = lon - lon0

    if method == "kruger":
        _, scale = tm_factors_geographic(a, b, dlon, lat, k0)
        return scale

    # Compute first and second eccentricity squared
    e2 = (a**2 - b**2) / a**2
    ep2 = e2 / (1 - e2)

    T = np.tan(lat)**2
    C = ep2 * np.cos(lat)**2
    A = dlon * np.cos(lat)

    # Compute scale factor (Snyder 1987, Transverse Mercator)
    scale = k0 * (
        1
        + (1 + C) * A**2 / 2
        + (5 - 4 * T + 42 * C + 13 * C**2 - 28 * ep2) * A**4 / 24
        + (61 - 148 * T + 16 * T**2) * A**6 / 720
    )

    return scale


def tm_point_scale_factor_projected(a: float, b: float, x, y, *,
                                    lat0: float = 0.0,
                                    false_easting: float = 0.0,
                                    false_northing: float = 0.0,
                                    k0: float = 1.0,
                                    radians: bool = False,
                                    method: Method = "kruger") -> Union[float, np.ndarray]:
    """
    Compute the Transverse Mercator scale factor based on projected coordinates (x, y).

    Notes
    -----
    The scale factor represents the distortion introduced by the projection.
    It is computed using the projected coordinates (eastings and northings) rather than geographic coordinates.
    The scale factor is dimensionless, so the `radians` flag only controls the unit of `lat0`.

    The series is defined on *unscaled* coordinates measured from the natural origin, so the
    grid coordinates are first reduced with::

        x' = (E - false_easting) / k0
        y' = (N - false_northing) / k0

    and the resulting pure-TM scale factor is multiplied by `k0` at the end. Ignoring `k0`
    for UTM gives an error of about 400 ppm, and ignoring `false_northing` breaks the
    southern hemisphere completely (UTM uses ``false_northing = 10 000 000`` m there).

    Accuracy: with ``method="kruger"`` (the default) the error against a rigorous solution
    (pyproj/PROJ ``tmerc``) is below ``1e-10`` everywhere from the equator to the poles and
    out to ``±45°`` from the central meridian. With ``method="snyder"`` the series is
    truncated after the ``x**4`` term, giving about ``2e-11`` at ``±2°`` from the central
    meridian, ``7e-10`` at ``±4°`` and ``1e-8`` at ``±8°`` on GRS80 at ``lat = 59°``.

    Parameters
    ----------
    a : float. Semi-major axis of the ellipsoid (meters).
    b : float. Semi-minor axis of the ellipsoid (meters).
    x : float or numpy.ndarray. Easting coordinate (meters), including the false easting.
    y : float or numpy.ndarray. Northing coordinate (meters), including the false northing.
    lat0 : float, keyword-only, optional. Latitude of natural origin (degrees if radians=False,
           radians if radians=True). Defaults to `0.0`, which is the correct value for UTM.
    false_easting : float, keyword-only, optional. False easting of the projection (meters).
                    UTM uses 500 000 meters here. Defaults to `0.0`.
    false_northing : float, keyword-only, optional. False northing of the projection (meters).
                     UTM uses 0 on the northern hemisphere and 10 000 000 meters on the
                     southern hemisphere. Defaults to `0.0`.
    k0 : float, keyword-only, optional. Scale factor at the central meridian. UTM uses `0.9996`.
         Defaults to `1.0` (pure Transverse Mercator).
    radians : bool, keyword-only, optional. If `True`, `lat0` is in radians. If `False`, `lat0`
              is in degrees. Defaults to `False`.
    method : {"kruger", "snyder"}, keyword-only, optional. Series to use. `"kruger"` is the
             Kruger expansion in the third flattening (Karney 2011), accurate to double
             precision. `"snyder"` is the classic expansion in `x` truncated at the 4th
             order. Defaults to `"kruger"`.

    Returns
    -------
    scale : float or numpy.ndarray. Dimensionless point scale factor at the given projected coordinates.

    Examples
    --------
    ETRS89 / UTM zone 32N (EPSG:25832):

    >>> from pygeodetics.Ellipsoid import GRS80
    >>> ellip = GRS80()
    >>> k = tm_point_scale_factor_projected(
    ...     ellip.a, ellip.b, 280841.3814, 6571720.2565,
    ...     lat0=0.0, false_easting=500000.0, k0=0.9996)
    >>> round(float(k), 9)
    1.000188742
    """
    _check_method(method)
    if not radians:
        lat0 = np.radians(lat0)

    if method == "kruger":
        _, scale = tm_factors_projected(a, b, x, y, lat0, false_easting, false_northing, k0)
        return scale

    # Reduce the grid coordinates to unscaled coordinates relative to the natural origin.
    # Deliberately not `x -= ...`, which would mutate a caller-supplied NumPy array in place.
    x = (x - false_easting) / k0
    y = (y - false_northing) / k0

    # Compute footpoint latitude
    latf = footpoint_latitude(a, b, y, lat0, radians=True)

    # Compute normal and meridional radii of curvature
    Nf = Nrad(a, b, latf, radians=True)
    Mf = Mrad(a, b, latf, radians=True)

    # Compute Transverse Mercator scale factor
    scale = k0 * (1 + (x**2 / (2 * Mf * Nf)) + (x**4 / (24 * Nf**4)))

    return scale


def tm_sphere_point_scale_factor(x, false_esting: float, R: float = 6371000.0) -> Union[float, np.ndarray]:
    """
    Compute the point scale factor of a Transverse Mercator projection for a sphere (TMS) using projected coordinates.
    Note: Assumes earth is a perfect sphere.

    Based on the "THE MERCATOR PROJECTIONS" book from  Peter Osborne, 2013
    See "The scale factor for the TMS projection" section at page 61.

    Parameters
    ----------
    x : float or numpy.ndarray. X-coordinate (distance from the central meridian) in meters.
    false_esting : float. False easting value for the projection in meters.
    R : float, optional. Radius of the sphere in meters.

    Returns
    -------
    float or numpy.ndarray. Point scale factor at the given projection coordinates.
    """
    # Subtract false easting to get distance from the central meridian
    x = x - false_esting
    k = np.cosh(x / (R))
    return k




if __name__ == "__main__":

    from pyproj import Proj, Transformer
    from pygeodetics.Ellipsoid import WGS84, GRS80
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
    k0 = 0.9996  # UTM scale factor at the central meridian


    # Compute scale factor using the custom function
    custom_k_sphere = tm_sphere_point_scale_factor(x, false_esting=false_easting)
    print(f"Custom scale factor projected (sphere): {custom_k_sphere:.12f}")

    custom_k = tm_point_scale_factor_geographic(a, b, lon, lat, lon0, k0=k0, radians=True)
    print(f"Custom scale factor geograpic: {custom_k:.12f}")

    custom_k_proj = tm_point_scale_factor_projected(a, b, x, y, lat0=lat0, false_easting=false_easting,
                                                    k0=k0, radians=True)
    print(f"Custom scale factor projected: {custom_k_proj:.12f}")

    # Compute scale factor using pyproj
    projection = Proj(proj='utm', zone=32, ellps="WGS84")
    pyproj_k = projection.get_factors(longitude=lon, latitude=lat, radians=True)
    print(f"pyproj scale factor: {pyproj_k.meridional_scale:.12f}")

