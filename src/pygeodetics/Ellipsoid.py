"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""


from typing import Literal
import numpy as np


class Ellipsoid:
    """
    Class for defining an ellipsoid and performing
    calculations related to it.


    Attributes
    ----------
    a : float. Semi-major axis of the ellipsoid.
    b : float. Semi-minor axis of the ellipsoid.
    f : float. Flattening of the ellipsoid.
    inv_f : float. Inverse flattening of the ellipsoid.
    e2 : float. Square of the first eccentricity e^2.
    e : float. First eccentricity.
    e2_prime : float. Square of the second eccentricity e'^2.

    """

    def __init__(self, a: float, b: float = None, f: float = None):
        if b is None and f is None:
            raise ValueError("Either 'b' or 'f' must be defined.")
        self.a = a
        if b is not None and f is None:
            self.b = b
            self.f = (self.a - self.b) / self.a
        elif f is not None and b is None:
            self.f = f
            self.b = self.a * (1 - self.f)
        elif b is not None and f is not None:
            self.b = b
            self.f = f

        self.inv_f = 1 / self.f
        self.e2 = self.eccentricity_squared() # e^2
        self.e = np.sqrt(self.e2) # e
        self.e2_prime = self.e2 / (1 - self.e2)  # e'^2
        self.mean_radius = self.get_mean_radius()

    def __repr__(self):
        return f"Ellipsoid(a={self.a}, b={self.b}, f={self.f})"

    def eccentricity_squared(self) -> float:
        """Calculate the eccentricity of the ellipsoid."""
        return (self.a**2 - self.b**2) / self.a**2

    def calc_b(self) -> float:
        """Calculate the semi-minor axis of the ellipsoid."""
        return self.a * (1 - self.f)

    def flattening(self) -> float:
        """Calculate the flattening of the ellipsoid."""
        return (self.a - self.b) / self.a

    def get_mean_radius(self) -> float:
        """Calculate the mean radius of the ellipsoid."""
        return (2 * self.a + self.b) / 3

    def surface_area(self) -> float:
        """Calculate the surface area [km^2] of the ellipsoid."""
        surface_area_m2 = 2 * np.pi * self.a**2 + np.pi * (self.b**2 / self.e) * np.log((1 + self.e) / (1 - self.e))
        return surface_area_m2 / 1e6  # Convert from m^2 to km^2

    def marc(self, lat: float, angle_unit: Literal["deg", "rad"] = "rad") -> float:
        """
        Calculate the meridian arc length.

        Parameters
        ----------
        lat : float
            Latitude in degrees or radians.
        angle_unit : Literal["deg", "rad"], optional
            Unit of the latitude angle, either 'deg' or 'rad'. Defaults to 'rad'.

        Returns
        -------
        float
            Meridian arc length.
        """
        if angle_unit == 'deg':
            lat = np.radians(lat)

        b0 = self.a * (1 - 1/2 * self.f + 1/16 * self.f**2 + 1/32 * self.f**3)

        B = b0 * (lat - (3/4 * self.f + 3/8 * self.f**2 + 15/128 * self.f**3) * np.sin(2 * lat)
                + (15/64 * self.f**2 + 15/64 * self.f**3) * np.sin(4 * lat)
                - 35/384 * self.f**3 * np.sin(6 * lat))
        return B

    def mrad(self, lat: float, angle_unit: Literal["deg", "rad"] = "rad") -> float:
        """
        Calculate the Earth's meridional radius of curvature at a given latitude (north-south direction).

        Parameters
        ----------
        lat: float. Latitude in degrees or radians.
        angle_unit: Literal["deg", "rad"], optional
            Unit of the latitude angle, either 'deg' or 'rad'. Defaults to 'rad'.

        Returns
        -------
        float: Meridional radius of curvature (M).
        """
        if angle_unit == "deg":
            lat = np.radians(lat)

        M = self.a * (1 - self.e2) / (1 - self.e2 * np.sin(lat)**2)**(3/2)
        return M

    def nrad(self, lat: float, angle_unit: Literal["deg", "rad"] = "rad") -> float:
        """
        Calculate the Earth's transverse radius of curvature at a given latitude (east-west direction).

        Parameters
        ----------
        lat: float. Latitude in degrees or radians.
        angle_unit: Literal["deg", "rad"], optional
            Unit of the latitude angle, either 'deg' or 'rad'. Defaults to 'rad'.

        Returns
        -------
        float: Transverse radius of curvature (N).
        """
        if angle_unit == "deg":
            lat = np.radians(lat)

        N = self.a / (1 - self.e2 * np.sin(lat)**2)**(1/2)
        return N

    def footlat(self, y: float, lat0: float, a: float = None, b: float = None) -> float:
        """
        Compute the footpoint latitude (latitude of the origin of meridian arc).

        Parameters
        ----------
        y : float. Northing (meters).
        lat0 : float. Latitude of the natural origin (radians).
        a : float, optional. Semi-major axis of the ellipsoid (meters). If not provided, uses instance attribute.
        b : float, optional. Semi-minor axis of the ellipsoid (meters). If not provided, uses instance attribute.

        Returns
        -------
        float. Footpoint latitude (radians).
        """
        # Use provided a and b or default to instance attributes
        a = a if a is not None else self.a
        b = b if b is not None else self.b

        f = (a - b) / a
        b0 = a * (1 - 1 / 2 * f + 1 / 16 * f**2 + 1 / 32 * f**3)

        # Compute meridian arc length
        B = self.marc(lat0, angle_unit="rad") + y

        # Compute footpoint latitude
        latf = (
            B / b0
            + (3 / 4 * f + 3 / 8 * f**2 + 21 / 256 * f**3) * np.sin(2 * B / b0)
            + (21 / 64 * f**2 + 21 / 64 * f**3) * np.sin(4 * B / b0)
            + (151 / 768 * f**3) * np.sin(6 * B / b0)
        )
        return latf

    def mean_radis_for_latitude(self, lat: float, angle_unit: Literal["deg", "rad"] = "rad") -> float:
        """Calculate the mean radius of curvature for a given latitude"""
        N = self.nrad(lat, angle_unit)
        M = self.mrad(lat, angle_unit)
        return np.sqrt(M*N)

    def radius_of_curvature_azimuth(self, lat: float, az: float, angle_unit: Literal["deg", "rad"] = "deg") -> float:
        """
        This function calculates the radius of curvature in a given direction using Euler's equation.
        It considers the meridional and normal radius of curvature to determine the radius for a specific azimuth.

        Parameters
        ----------
        lat : float. Latitude in degrees or radians.
        az : float. Azimuth angle in degrees or radians.
        angle_unit : Literal["deg", "rad"], optional. Specifies whether input angles are in degrees or radians.
                    Defaults to "rad".

        Returns
        -------
        float. Radius of curvature for the given azimuth (meters).
        """
        if angle_unit == "deg":
            lat, az = np.radians(lat), np.radians(az)

        M = self.mrad(lat, angle_unit="rad")  # Meridional radius of curvature
        N = self.nrad(lat, angle_unit="rad")  # Transverse radius of curvature

        return (M * N) / (M * np.sin(az)**2 + N * np.cos(az)**2)

    def tm_conv(self, x: float, y: float, lat0: float, false_easting: float,
                a: float = None, b: float = None) -> float:
        """
        Compute the meridian convergence (gamma) at a point (x, y).

        Parameters
        ----------
        x : float. Easting (meters).
        y : float. Northing (meters).
        lat0 : float. Latitude of the natural origin (radians).
        false_easting: float. False easting to be subtracted from the easting coordinates
        a : float, optional. Semi-major axis of the ellipsoid (meters). If not provided, uses instance attribute.
        b : float, optional. Semi-minor axis of the ellipsoid (meters). If not provided, uses instance attribute.

        Returns
        -------
        float. Meridian convergence (gamma) in radians.
        """
        a = a if a is not None else self.a
        b = b if b is not None else self.b

        x = x - false_easting

        e2 = (a**2 - b**2) / a**2
        latf = self.footlat(y, lat0, a, b)
        N = a / np.sqrt(1 - e2 * np.sin(latf) ** 2)
        eps2 = (e2 / (1 - e2)) * np.cos(latf) ** 2
        gamma = (
            (x * np.tan(latf) / N)
            - ((x**3 * np.tan(latf)) / (3 * N**3)) * (1 + np.tan(latf) ** 2 - eps2 - 2 * eps2**2)
        )
        return gamma




class WGS84(Ellipsoid):
    def __init__(self):
        super().__init__(a=6378137, f=1/298.257223563)

class GRS80(Ellipsoid):
    def __init__(self):
        super().__init__(a=6378137, f=1/298.257222101)

class International1924(Ellipsoid):
    def __init__(self):
        super().__init__(a=6378388, f=1/297)

class Clarke1866(Ellipsoid):
    def __init__(self):
        super().__init__(a=6378206.4, f=1/294.9786982)

class BesselModified(Ellipsoid):
    def __init__(self):
        super().__init__(a=6377492.018, f=1/299.1528128)

class Bessel1841(Ellipsoid):
    def __init__(self):
        super().__init__(a=6377397.155, f=1/299.1528128)



if __name__ == "__main__":
    # ellipsoid = Ellipsoid(a=6378137, f=1/298.257223563)
    ellipsoid = WGS84()
    print(ellipsoid)
    print("b =", ellipsoid.calc_b())
    print("Flattening =", ellipsoid.flattening())
    print("Eccentricity Squared (e^2) =", ellipsoid.eccentricity_squared())
    print("Eccentricity (e) =", ellipsoid.e)
    print("Eccentricity Prime Squared (e'^2) =", ellipsoid.e2_prime)
    print("The arc length along the meridian at 45° latitude:", ellipsoid.marc(45, angle_unit='deg'))
    print("Meridional Radius (M) at 45° latitude:", ellipsoid.mrad(45, angle_unit='deg'))
    print("Transverse Radius (N) at 45° latitude:", ellipsoid.nrad(45, angle_unit='deg'))
    print("Mean radius:", ellipsoid.mean_radius)
    print("Surface area:", ellipsoid.surface_area())
    print("Mean radius for latitude:", ellipsoid.mean_radis_for_latitude(45, angle_unit='deg'))
    print(f"Footpoint Latitude: {np.degrees(ellipsoid.footlat(1000000.0, np.radians(0.0))):.6f} degrees")


    x,y = 555776.2667516097, 6651832.735433666
    lat0 = 0
    false_easting = 500000
    gamma = ellipsoid.tm_conv(x, y, lat0, false_easting)
    print(f"Meridian Convergence (gamma): {np.degrees(gamma):.6f} degrees")


