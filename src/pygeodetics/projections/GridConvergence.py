"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com

"""

from typing import Literal, Union
import numpy as np
from ._kruger import tm_factors_geographic, tm_factors_projected
from ..geodetics.footpoint_latitude import footpoint_latitude
from ..geodetics.Nrad import Nrad

Method = Literal["kruger", "snyder"]


def _check_method(method: str) -> None:
    if method not in ("kruger", "snyder"):
        raise ValueError(f"method must be 'kruger' or 'snyder', got {method!r}")


def tm_grid_convergence_geographic(a: float, b: float, lon, lat,
                                   lon0: float, *, radians: bool = False,
                                   method: Method = "kruger") -> Union[float, np.ndarray]:
    """
    Compute the Transverse Mercator meridian convergence (grid convergence) based on geographic coordinates (latitude, longitude).

    Notes
    -----
    The meridian convergence (gamma) is the angular difference between grid north and true north in the Transverse Mercator projection.


    ``lon`` and ``lon0`` must be referenced to the **same prime meridian**. Mixing frames is a
    silent error that still returns a plausible-looking angle. This is a real trap for grids
    whose EPSG definition uses a non-Greenwich prime meridian, e.g. the Norwegian NGO 1948
    axes: EPSG:27391 (NGO zone I) lists a longitude of natural origin of ``4°40'00"W``, which
    is ``-4.666667°`` **relative to Oslo**, not to Greenwich. Either convert the central
    meridian to Greenwich, or convert the point longitude to the local frame; both of the
    following are equivalent::

        # Greenwich throughout: Oslo (10°43'22.5"E) - 4°40' = 6°03'22.5"E
        tm_grid_convergence_geographic(a, b, lon_greenwich, lat, lon0=6.05625)

        # Oslo-relative throughout:
        tm_grid_convergence_geographic(a, b, lon_greenwich - 10.7229166667, lat, lon0=-4.6666667)

    Accuracy: with ``method="kruger"`` (the default) the series is in the third flattening
    ``n``, truncated at ``n**6``. Since ``n`` is about ``1/599``, the neglected term is of
    order ``1e-19`` and the result is exact to double precision. Measured against a rigorous
    solution (pyproj/PROJ ``tmerc``) the error is below ``3e-10°`` (0.000001") everywhere
    from the equator to the poles and out to ``±45°`` from the central meridian.

    With ``method="snyder"`` the series is in ``dlon``, truncated after the ``dlon**5`` term.
    On GRS80 at ``lat = 59°`` the error is then ``5e-10°`` at ``±2°`` from the central
    meridian, ``2e-8°`` (0.00007") at ``±4°`` and ``7e-7°`` (0.003") at ``±8°``, and it
    degrades rapidly beyond that.

    Parameters
    ----------
    a : float. Semi-major axis of the ellipsoid (meters).
    b : float. Semi-minor axis of the ellipsoid (meters).
    lon : float or numpy.ndarray. Geodetic longitude in degrees (if `radians=False`) or radians (if `radians=True`).
    lat : float or numpy.ndarray. Geodetic latitude in degrees (if `radians=False`) or radians (if `radians=True`).
    lon0 : float. Central meridian of the projection in degrees (if `radians=False`) or radians
           (if `radians=True`). The tangent meridian (Ex 9° for UTM Z32N).
    radians : bool, keyword-only, optional. If `True`, input and output are in radians. If False, input is in degrees and output is converted to degrees. Defaults to `False`.
    method : {"kruger", "snyder"}, keyword-only, optional. Series to use. `"kruger"` is the
             Kruger expansion in the third flattening (Karney 2011), accurate to double
             precision. `"snyder"` is the classic expansion in `dlon` truncated at the 5th
             order. Defaults to `"kruger"`.


    Returns
    -------
    gamma : float or numpy.ndarray. Meridian convergence in degrees (if radians=False) or radians (if radians=True).

    Examples
    --------
    >>> from pygeodetics.Ellipsoid import WGS84
    >>> ellip = WGS84()
    >>> gamma = tm_grid_convergence_geographic(ellip.a, ellip.b, 10.0, 60.0, 9.0)
    >>> round(float(gamma), 9)
    0.866047499
    """
    _check_method(method)
    if not radians:
        lon, lat, lon0 = np.radians(lon), np.radians(lat), np.radians(lon0)

    dlon = lon - lon0  # Difference in longitude from the central meridian

    if method == "kruger":
        gamma, _ = tm_factors_geographic(a, b, dlon, lat)
        return gamma if radians else np.degrees(gamma)

    e2 = (a**2 - b**2) / a**2
    eps2 = e2 / (1 - e2) * np.cos(lat)**2

    gamma = (
        dlon * np.sin(lat)
        + (dlon**3 / 3) * np.sin(lat) * np.cos(lat)**2 * (1 + 3 * eps2 + 2 * eps2**2)
        + (dlon**5 / 15) * np.sin(lat) * np.cos(lat)**4 * (2 - np.tan(lat)**2)
    )

    if not radians:
        gamma = np.degrees(gamma)

    return gamma


def tm_grid_convergence_projected(a: float, b: float, x, y, *,
                                  lat0: float = 0.0,
                                  false_easting: float = 0.0,
                                  false_northing: float = 0.0,
                                  k0: float = 1.0,
                                  radians: bool = False,
                                  method: Method = "kruger") -> Union[float, np.ndarray]:
    """
    Compute the Transverse Mercator meridian convergence (grid convergence) based on projected coordinates (x, y).

    Notes
    -----
    The meridian convergence (gamma) is the angular difference between grid north and true north.
    It is computed using projected coordinates (eastings and northings) rather than geographic coordinates.

    The series is defined on *unscaled* coordinates measured from the natural origin, so the
    grid coordinates are first reduced with::

        x' = (E - false_easting) / k0
        y' = (N - false_northing) / k0

    Grid convergence is an angle, so there is **no** final multiplication by ``k0`` - only the
    unscaling above. Ignoring ``k0`` for UTM (``k0 = 0.9996``) gives an error of roughly 16"
    near the zone edge, and ignoring ``false_northing`` breaks the southern hemisphere
    completely (UTM uses ``false_northing = 10 000 000`` m there).

    Accuracy: with ``method="kruger"`` (the default) the series is in the third flattening
    ``n``, truncated at ``n**6``, which is exact to double precision. Measured against a
    rigorous solution (pyproj/PROJ ``tmerc``) the error is below ``3e-10°`` (0.000001")
    everywhere from the equator to the poles and out to ``±45°`` from the central meridian.

    With ``method="snyder"`` the series is in ``x``, truncated after the ``x**5`` term. On
    GRS80 at ``lat = 59°`` the error is then ``2e-10°`` at ``±2°`` from the central meridian,
    ``4e-8°`` (0.0002") at ``±4°`` and ``6e-6°`` (0.02") at ``±8°``. That variant degrades
    faster than the geographic one because its series is in ``x`` rather than in ``dlon``.

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
    radians : bool, keyword-only, optional. If `True`, `lat0` and the returned angle are in radians.
              If `False`, both are in degrees. Defaults to `False`.
    method : {"kruger", "snyder"}, keyword-only, optional. Series to use. `"kruger"` is the
             Kruger expansion in the third flattening (Karney 2011), accurate to double
             precision. `"snyder"` is the classic expansion in `x` truncated at the 5th
             order. Defaults to `"kruger"`.

    Returns
    -------
    gamma : float or numpy.ndarray. Meridian convergence in degrees (if radians=False) or radians (if radians=True).

    Examples
    --------
    ETRS89 / UTM zone 32N (EPSG:25832):

    >>> from pygeodetics.Ellipsoid import GRS80
    >>> ellip = GRS80()
    >>> gamma = tm_grid_convergence_projected(
    ...     ellip.a, ellip.b, 280841.3814, 6571720.2565,
    ...     lat0=0.0, false_easting=500000.0, k0=0.9996)
    >>> round(float(gamma), 6)
    -3.301834
    """
    _check_method(method)
    if not radians:
        lat0 = np.radians(lat0)  # Convert latitude of origin to radians

    if method == "kruger":
        gamma, _ = tm_factors_projected(a, b, x, y, lat0, false_easting, false_northing, k0)
        return gamma if radians else np.degrees(gamma)

    # Reduce the grid coordinates to unscaled coordinates relative to the natural origin.
    # Deliberately not `x -= ...`, which would mutate a caller-supplied NumPy array in place.
    x = (x - false_easting) / k0
    y = (y - false_northing) / k0

    # Compute footpoint latitude
    latf = footpoint_latitude(a, b, y, lat0, radians=True)  # Footpoint latitude in radians

    # Compute first eccentricity squared
    e2 = (a**2 - b**2) / a**2

    # Compute normal radius of curvature at footpoint latitude
    Nf = Nrad(a, b, latf, radians=True)

    # Compute second eccentricity squared at footpoint latitude
    epsf2 = e2 / (1 - e2) * np.cos(latf)**2
    Tf = np.tan(latf)**2

    # Compute meridian convergence (gamma)
    gamma = (x * np.tan(latf) / Nf) * (
        1
        - (x**2 / (3 * Nf**2)) * (1 + Tf - epsf2 - 2 * epsf2**2)
        + (x**4 / (15 * Nf**4)) * (2 + 5 * Tf + 3 * Tf**2)
    )

    if not radians:
        gamma = np.degrees(gamma)

    return gamma




def tm_sphere_grid_conv_projected(x, y, false_easting: float, R: float = 6371000,
                                  angle_unit: Literal["deg", "rad"] = "deg") -> Union[float, np.ndarray]:
    """
    Compute the grid convergence of a Transverse Mercator projection on a sphere (TMS)
    by using projection coordinates.

    Note: Assumes earth is a perfect sphere. Based on the "THE MERCATOR PROJECTIONS" book from  Peter Osborne, 2013
    See "The scale factor for the TMS projection" section at page 63.

    Parameters
    ----------
    x : float or numpy.ndarray. X-coordinate (distance from the central meridian) in meters.
    y : float or numpy.ndarray. Y-coordinate (distance from the equator) in meters.
    false_easting: float. The false easting value in meters that is used by the projection
    R : float, optional. Radius of the sphere in meters. Default is 6371000 meters.
    angle_unit : str, optional. Unit of the grid convergence angle. Default is degrees

    Returns
    -------
    float or numpy.ndarray. Grid convergence angle (gamma), in degrees if `angle_unit="deg"`
    and in radians if `angle_unit="rad"`.

    """
    x = x - false_easting
    gamma = np.arctan(np.tanh(x / (R)) * np.tan(y / (R)))
    if angle_unit == "deg":
        gamma = np.degrees(gamma)
    return gamma





if __name__ == "__main__":
    from pyproj import Proj, Transformer
    from pygeodetics.Ellipsoid import WGS84

    ellip = WGS84()
    a = ellip.a
    b = ellip.b
    lon = np.radians(10)  # Longitude in radians
    lat = np.radians(60)  # Latitude in radians
    x, y = Transformer.from_crs("EPSG:4326", "EPSG:32632", always_xy=True).transform(lon, lat, radians=True)
    false_easting = 500000  # false easting in meters
    lon0 = np.deg2rad(9)  # central meridian/tangent meridian
    lat0 = np.deg2rad(0)
    k0 = 0.9996  # UTM scale factor at the central meridian


    # Grid convergence in projection coordinates
    gamma_proj = tm_sphere_grid_conv_projected(x, y, false_easting)
    print(f"Grid convergence (projection): {gamma_proj:.6f} degrees\n")

    gamma_geog = tm_grid_convergence_geographic(a, b, lon, lat, lon0, radians=True)
    print(f"Grid convergence (geographic): {np.rad2deg(gamma_geog):.12f} degrees\n")

    gamma_proj = tm_grid_convergence_projected(a, b, x, y, lat0=lat0, false_easting=false_easting,
                                               k0=k0, radians=True)
    print(f"Grid convergence (projected NEW): {np.rad2deg(gamma_proj):.12f} degrees\n")

    # projection = Proj(proj="tmerc", lon_0=9, ellps="WGS84")
    projection = Proj(proj='utm',zone=32,ellps='WGS84')
    pyproj_k = projection.get_factors(longitude=lon, latitude=lat, radians=True)
    print(f"Pyproj projection factors: \n - Meridian convergence: {pyproj_k.meridian_convergence}")

    # print(f"pyproj scale factor: {np.degrees(pyproj_k)}")
