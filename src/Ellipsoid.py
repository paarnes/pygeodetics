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
        if b is not None:
            self.b = b
            self.f = (self.a - self.b) / self.a
        elif f is not None:
            self.f = f
            self.b = self.a * (1 - self.f)
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

    def mean_radis_for_latitude(self, lat: float, angle_unit: Literal["deg", "rad"] = "rad") -> float:
        """Calculate the mean radius of curvature for a given latitude"""
        N = self.nrad(lat, angle_unit)
        M = self.mrad(lat, angle_unit)
        return np.sqrt(M*N)




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


