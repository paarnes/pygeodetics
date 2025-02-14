"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""

from typing import Literal, Optional, Tuple
import numpy as np
from Ellipsoid import Ellipsoid, WGS84
from geodetics.ECEF2enu import ECEF2enu
from geodetics.ECEF2geod import ECEF2geodb
from geodetics.geod2ECEF import geod2ECEF
from geodetics.geodetic_inverse_problem import geodetic_inverse_problem
from geodetics.geodetic_direct_problem import geodetic_direct_problem
from geodetics.radius_of_curvature_azimuth import radius_of_curvature_azimuth


class Geodetic(Ellipsoid):
    """
    Geodetic class for performing various geodetic calculations on an ellipsoid.

    Inherits from the `Ellipsoid` class, allowing geodetic computations on different reference ellipsoids.

    Attributes
    ----------
    a : float. Semi-major axis of the ellipsoid (meters).
    b : float. Semi-minor axis of the ellipsoid (meters).
    """

    def __init__(self, ellipsoid: Optional[Ellipsoid] = None, a: Optional[float] = None,
                 b: Optional[float] = None, f: Optional[float] = None):
        """
        Initialize the Geodetic class with a specific ellipsoid or custom ellipsoid parameters.

        Parameters
        ----------
        ellipsoid : Ellipsoid, optional. Reference ellipsoid to use. Default is WGS84.
        a : float, optional. Semi-major axis of the ellipsoid (meters).
        b : float, optional. Semi-minor axis of the ellipsoid (meters). If None, will be computed from `a` and `f`.
        f : float, optional. Flattening of the ellipsoid. Used only if `b` is not provided.
        """
        if a is not None and (b is not None or f is not None):
            super().__init__(a=a, b=b, f=f)  # Use Ellipsoid constructor logic
            self.ellipsoid_name = "CustomEllipsoid"
            self.ellipsoid = ellipsoid
        else:
            # Use the provided ellipsoid, defaulting to WGS84 if None
            self.ellipsoid = ellipsoid if ellipsoid else WGS84()
            super().__init__(a=self.ellipsoid.a, b=self.ellipsoid.b)
            self.ellipsoid_name = self.ellipsoid.__class__.__name__

    def __repr__(self):
        return f"Geodetic(ellipsoid={self.ellipsoid_name}, a={self.a}, b={self.b}, f={self.f})"

    def ecef2enu(self, lat: float, lon: float, X: float, Y: float, Z: float, radians: bool = True):
        """
        Converts from ECEF coordinates to the local topocentric ENU (East-North-Up) coordinate system.

        Parameters
        ----------
        lat : float. Geodetic latitude of the reference point in radians.
        lon : float. Geodetic longitude of the reference point in radians.
        X : float. X coordinate (ECEF in meters).
        Y : float. Y coordinate (ECEF in meters).
        Z : float. Z coordinate (ECEF in meters).
        radians : bool, optional. If `False`, assumes `lat` and `lon` are in degrees and converts them to radians.
            Defaults to `True` (assumes input is already in radians).

        Returns
        -------
        e : float. East coordinate in the ENU coordinate system (meters).
        n : float. North coordinate in the ENU coordinate system (meters).
        u : float. Up coordinate in the ENU coordinate system (meters).
        """
        return ECEF2enu(lat, lon, X, Y, Z)

    def ecef2ned(self, lat: float, lon: float, X: float, Y: float, Z: float):
        """
        Converts from ECEF coordinates to the local topocentric NED (North-East-Down) coordinate system.

        Parameters
        ----------
        lat : float. Geodetic latitude of the reference point in radians.
        lon : float. Geodetic longitude of the reference point in radians.
        X : float. X coordinate (ECEF in meters).
        Y : float. Y coordinate (ECEF in meters).
        Z : float. Z coordinate (ECEF in meters).

        Returns
        -------
        n : float. North coordinate in the NED coordinate system (meters).
        e : float. East coordinate in the NED coordinate system (meters).
        d : float. Down coordinate in the NED coordinate system (meters).

        Examples
        --------
        >>> geod = Geodetic(WGS84())
        >>> n, e, d = geod.ecef2ned(np.radians(45.0), np.radians(9.0), 1111321.0, 1234567.0, 5678901.0)
        >>> print(f"N: {n:.6f}, E: {e:.6f}, D: {d:.6f}")
        """
        e, n, u = self.ecef2enu(lat, lon, X, Y, Z)
        d = -u
        return n, e, d

    def inverse_problem(self, lat1: float, lon1: float, lat2: float, lon2: float, quadrant_correction: bool = False):
        """
        Compute the geodetic inverse problem: given two points (lat1, lon1) and (lat2, lon2),
        determine the forward azimuth, reverse azimuth, and geodesic distance.

        Parameters
        ----------
        lat1, lon1 : float. Latitude and longitude of point 1 in radians.
        lat2, lon2 : float. Latitude and longitude of point 2 in radians.
        quadrant_correction : bool, optional. If True, ensures azimuths are in range [0, 2π]. Default is False.

        Returns
        -------
        az1 : float. Forward azimuth at point 1 (radians).
        az2 : float. Reverse azimuth at point 2 (radians).
        d : float. Geodesic distance between the two points (meters).
        """
        return geodetic_inverse_problem(self.a, self.b, lat1, lon1, lat2, lon2, quadrant_correction)

    def direct_problem(self, lat1: float, lon1: float, az1: float, d: float, quadrant_correction: bool = False):
        """
        Compute the geodetic direct problem: given an initial point (lat1, lon1), an initial azimuth (az1),
        and a geodesic distance (d), determine the final point coordinates and azimuth.

        Parameters
        ----------
        lat1, lon1 : float. Latitude and longitude of starting point in radians.
        az1 : float. Initial azimuth in radians.
        d : float. Distance to travel along the geodesic (meters).
        quadrant_correction : bool, optional. If True, ensures azimuths are in range [0, 2π]. Default is False.

        Returns
        -------
        lat2 : float. Latitude of the destination point in radians.
        lon2 : float. Longitude of the destination point in radians.
        az2 : float. Reverse azimuth at the destination point in radians.
        """
        return geodetic_direct_problem(self.a, self.b, lat1, lon1, az1, d, quadrant_correction)

    def radius_of_curvature(self, lat: float, az: float, radians: bool = False) -> float:
        """
        Compute the radius of curvature for a given azimuth using Euler's equation.

        Parameters
        ----------
        lat : float. Geodetic latitude in degrees (if radians=False) or radians (if radians=True).
        az : float. Azimuth angle in degrees (if radians=False) or radians (if radians=True).
        radians : bool, optional. If `True`, input and output are in radians. Default is `False`.

        Returns
        -------
        float. Radius of curvature for the given azimuth (meters).
        """
        return radius_of_curvature_azimuth(self.a, self.b, lat, az, radians)

    def ecef2geod(self, X: float, Y: float, Z: float, angle_unit: Literal['deg', 'rad'] = "deg"):
        """
        Convert Earth-Centered, Earth-Fixed (ECEF) coordinates to geodetic coordinates (latitude, longitude, height).

        Notes
        -----
        This function computes geodetic latitude, longitude, and height above the ellipsoid
        from ECEF coordinates using an Bowrings method.

        Parameters
        ----------
        a : float. Semi-major axis of the reference ellipsoid (meters).
        b : float. Semi-minor axis of the reference ellipsoid (meters).
        X : float. ECEF X coordinate (meters).
        Y : float. ECEF Y coordinate (meters).
        Z : float. ECEF Z coordinate (meters).
        angle_unit : {'deg', 'rad'}, optional
            Unit of the output latitude and longitude. If 'deg' (default), results are in degrees.
            If 'rad', results are in radians.

        Returns
        -------
        lat : float. Geodetic latitude (degrees if `angle_unit='deg'`, radians if `angle_unit='rad'`).
        lon : float. Geodetic longitude (degrees if `angle_unit='deg'`, radians if `angle_unit='rad'`).
        h : float. Height above the ellipsoid (meters).

        """
        return ECEF2geodb(self.a, self.b, X, Y, Z, angle_unit)

    def geod2ecef(self, lat: float, lon: float, h: float, radians: bool = False) -> Tuple[float, float, float]:
        """
        Convert geodetic coordinates (latitude, longitude, height) to ECEF (Earth-Centered, Earth-Fixed) coordinates.

        Parameters
        ----------
        lat : float. Geodetic latitude in degrees (if radians=False) or radians (if radians=True).
        lon : float. Geodetic longitude in degrees (if radians=False) or radians (if radians=True).
        h : float. Height above the ellipsoid (meters).
        radians : bool, optional. If `False`, assumes `lat` and `lon` are in degrees and converts them to radians.

        Returns
        -------
        X : float. ECEF X coordinate (meters).
        Y : float. ECEF Y coordinate (meters).
        Z : float. ECEF Z coordinate (meters).

        Examples
        --------
        >>> geod = Geodetic(WGS84())
        >>> X, Y, Z = geod.geod2ECEF(60.0, 10.0, 100.0, radians=False)
        >>> print(f"X: {X:.3f}, Y: {Y:.3f}, Z: {Z:.3f}")
        """
        return geod2ECEF(lat, lon, h, ellipsoid=self.ellipsoid, radians=radians)









if __name__ == "__main__":
    # Example usage
    geod = Geodetic(WGS84())
    print(geod)

    # Geodetic inverse problem
    lat1 = np.radians(52.2296756)
    lon1 = np.radians(21.0122287)
    lat2 = np.radians(41.8919300)
    lon2 = np.radians(12.5113300)

    az1, az2, distance = geod.inverse_problem(lat1, lon1, lat2, lon2, quadrant_correction=True)

    print(f"Forward Azimuth: {np.degrees(az1):.10f} degrees")
    print(f"Reverse Azimuth: {np.degrees(az2):.10f} degrees")
    print(f"Distance: {distance:.10f} meters\n")

    # Geodetic direct problem
    az1 = np.radians(-147.4628043168)
    d = 1316208.08334

    lat2, lon2, az2 = geod.direct_problem(lat1, lon1, az1, d, quadrant_correction=True)

    print(f"Destination Latitude: {np.degrees(lat2):.10f} degrees")
    print(f"Destination Longitude: {np.degrees(lon2):.10f} degrees")
    print(f"Final Azimuth: {np.degrees(az2):.10f} degrees")

    # Radius of curvature for given azimuth
    radius = geod.radius_of_curvature(45.0, 30.0, radians=False)
    print(f"Radius of Curvature: {radius:.3f} meters")

    # Get mean radius
    mean_radius = geod.get_mean_radius()
    print(f"Mean radius: {mean_radius:.3f} meters")


    # ECEF to Geodetic example
    X, Y, Z = 3149785.9652, 598260.8822, 5495348.4927
    lat, lon, h = geod.ecef2geod(X, Y, Z, angle_unit='deg')

    print(f"ECEF to Geodetic:")
    print(f"Latitude: {lat:.10f} degrees")
    print(f"Longitude: {lon:.10f} degrees")
    print(f"Height: {h:.3f} meters\n")

    # Geodetic to ECEF example
    lat_geo = 59.907072474276958  # Latitude in degrees
    lon_geo = 10.754482924017791  # Longitude in degrees
    h_geo = 63.8281  # Height in meters


    X, Y, Z = geod.geod2ecef(lat_geo, lon_geo, h_geo, radians=False)
    print(f"Geodetic to ECEF:")
    print(f"X: {X:.4f} meters")
    print(f"Y: {Y:.4f} meters")
    print(f"Z: {Z:.4f} meters\n")

