"""
Transverse Mercator Projection Module
-------------------------------------

This module provides an implementation of the Transverse Mercator projection,
including forward and inverse transformations. The implementation is based on
the JHS formulas and the document "Geomatics Guidance Note 7, part 2
Coordinate Conversions & Transformations including Formulas" from IOGP.


Usage Example:
--------------
    # Define projection parameters
    lat_origin = 0  # Latitude of natural origin in radians
    lon_origin = math.radians(9)  # Longitude of natural origin in radians
    scale_factor = 0.9996  # Scale factor at the natural origin
    false_easting = 500000  # False easting in meters
    false_northing = 0  # False northing in meters

    # Ellipsoid parameters (WGS84)
    a = 6378137.0  # Semi-major axis in meters
    f = 1 / 298.257223563  # Flattening

    # Create an instance of the TransverseMercator class
    tm = TransverseMercator(lat_origin, lon_origin, scale_factor, false_easting, false_northing, a, f)

    # Example geographic coordinates (latitude, longitude in radians)
    lat = math.radians(60.0)  # Latitude in radians
    lon = math.radians(3.0)  # Longitude in radians

    # Perform forward projection
    easting, northing = tm.forward(lat, lon)
    print(f"Projected Coordinates:\nEasting = {round(easting,4)}\nNorthing = {round(northing,4)}")

    # Perform inverse projection
    lat_back, lon_back = tm.inverse(easting, northing)
    print(f"Geographic Coordinates: Latitude = {math.degrees(lat_back)}, Longitude = {math.degrees(lon_back)}")
"""

from typing import Tuple
import math



class TransverseMercator:
    """
    A class to handle Transverse Mercator projections including forward and reverse transformations.

    Attributes:
    - lat_origin (float): Latitude of natural origin (in radians).
    - lon_origin (float): Longitude of natural origin (in radians).
    - scale_factor (float): Scale factor at the natural origin.
    - false_easting (float): False easting.
    - false_northing (float): False northing.
    - a (float): Semi-major axis of the ellipsoid.
    - f (float): Flattening of the ellipsoid.
    """

    def __init__(self, lat_origin, lon_origin, scale_factor, false_easting, false_northing, a, f):
        self.lat_origin = lat_origin
        self.lon_origin = lon_origin
        self.scale_factor = scale_factor
        self.false_easting = false_easting
        self.false_northing = false_northing
        self.close_to_poles = self.lat_orig_is_close_to_poles()
        self.a = a
        self.f = f
        self.n = f / (2 - f)
        self.B = a / (1 + self.n) * (1 + self.n**2 / 4 + self.n**4 / 64)
        self.h_coeffs = self._calculate_h_coefficients()

    def _calculate_h_coefficients(self):
        """
        Calculate the coefficients h1, h2, h3, h4 for the series expansion.
        """
        n = self.n
        return {
            "h1": n / 2 - (2 / 3) * n**2 + (5 / 16) * n**3 + (41 / 180) * n**4,
            "h2": (13 / 48) * n**2 - (3 / 5) * n**3 + (557 / 1440) * n**4,
            "h3": (61 / 240) * n**3 - (103 / 140) * n**4,
            "h4": (49561 / 161280) * n**4
        }

    def lat_orig_is_close_to_poles(self) -> bool:
        """
        Latitude of the origin is closer than 2" to a pole,
        but not on at a pole.
        """
        is_close_to_poles = False
        if round(self.lat_origin, 8) == 0:
            return False
        elif round(self.lat_origin, 8) == abs(round(math.pi/2,8)):
            return False
        if abs(self.lat_origin - math.pi / 2) < 2 / 3600 * math.pi / 180:
            is_close_to_poles = True
        return is_close_to_poles

    def _calculate_meridional_arc_pole_safe(self, lat) -> float:
        """
        Calculate the meridional arc distance using a series expansion.
        """
        e2 = self.f * (2 - self.f)
        a0 = 1 - e2 / 4 - 3 * e2**2 / 64 - 5 * e2**3 / 256
        a2 = (3 / 8) * (e2 + e2**2 / 4 + 15 * e2**3 / 128)
        a4 = (15 / 256) * (e2**2 + 3 * e2**3 / 4)
        a6 = (35 / 3072) * e2**3

        M = self.a * (a0 * lat - a2 * math.sin(2 * lat) + a4 * math.sin(4 * lat) - a6 * math.sin(6 * lat))
        return M

    def _calculate_meridional_arc(self) -> float:
        """
        Calculate the meridional arc distance from the equator to the origin latitude
        using the alternative method for latitudes close to the poles.

        Returns:
        - M0 (float): Meridional arc distance in meters.
        """
        # Handle special cases for latitude of origin
        if abs(self.lat_origin) < 1e-8:  # Latitude of origin is 0
            return 0
        elif abs(self.lat_origin - math.pi / 2) < 1e-8:  # Latitude of origin is +90°N
            return self.B * (math.pi / 2)
        elif abs(self.lat_origin + math.pi / 2) < 1e-8:  # Latitude of origin is -90°S
            return self.B * (-math.pi / 2)

        # General case for latitude of origin
        e2 = self.f * (2 - self.f)  # Eccentricity squared
        e = math.sqrt(e2)  # Eccentricity

        # Compute Q0
        Qo = math.asinh(math.tan(self.lat_origin)) - e * math.atanh(e * math.sin(self.lat_origin))
        beta_o = math.atan(math.sinh(Qo))

        # Compute xi_o0
        xi_o0 = math.asin(math.sin(beta_o))

        # Compute xi_o using the series expansion
        xi_o1 = self.h_coeffs[f"h1"] * math.sin(2*xi_o0)
        xi_o2 = self.h_coeffs[f"h2"] * math.sin(4*xi_o0)
        xi_o3 = self.h_coeffs[f"h3"] * math.sin(6*xi_o0)
        xi_o4 = self.h_coeffs[f"h4"] * math.sin(8*xi_o0)
        xi_o = xi_o0 + xi_o1 + xi_o2 + xi_o3 + xi_o4

        # Compute M0
        M0 = self.B * xi_o

        return M0

    def calculate_xi_eta_pole_safe(self, eta0, xi0) -> Tuple[float, float]:
        # Compute full xi and eta
        xi = xi0
        eta = eta0
        for i in range(1, 5):
            xi += self.h_coeffs[f"h{i}"] * math.sin(2 * i * xi0) * math.cosh(2 * i * eta0)
            eta += self.h_coeffs[f"h{i}"] * math.cos(2 * i * xi0) * math.sinh(2 * i * eta0)
        return xi, eta

    def geog_to_projected(self, lat, lon):
        """
        Convert geographic coordinates (latitude, longitude) to projected coordinates (easting, northing).

        Args:
        - lat (float): Latitude in radians.
        - lon (float): Longitude in radians.

        Returns:
        - (float, float): Easting (E) and Northing (N).
        """
        # First, calculate Q and beta
        e2 = self.f * (2 - self.f)
        e = math.sqrt(e2) # Take the square root to get the eccentricity
        Q = math.asinh(math.tan(lat)) - e * math.atanh(e * math.sin(lat))
        beta = math.atan(math.sinh(Q))

        # Compute initial xi0 and eta0
        eta0 = math.atanh(math.cos(beta) * math.sin(lon - self.lon_origin))
        xi0 = math.asin(math.sin(beta) * math.cosh(eta0))

        # Compute full xi and eta
        xi, eta = self.calculate_xi_eta_pole_safe(eta0, xi0)

        # Compute M0 (Meridional Arc at Origin)
        if self.close_to_poles:
            M0 = self._calculate_meridional_arc_pole_safe(self.lat_origin)
        else:
            M0 = self._calculate_meridional_arc()
        # M0 = self.B * xi # test = ξO0 + ξO1+ ξO2 + ξO3 + ξO4 (test*self.B = 5429228.602101349)
        # M0 = 5429228.602
        # Compute Easting and Northing
        E = self.false_easting + self.scale_factor * self.B * eta
        N = self.false_northing + self.scale_factor * (self.B * xi - M0)

        return E, N

    def projected_to_geog(self, E, N):
        """
        Convert projected coordinates (easting, northing) back to geographic coordinates (latitude, longitude).

        Args:
        - E (float): Easting.
        - N (float): Northing.

        Returns:
        - (float, float): Latitude (lat) and Longitude (lon) in radians.
        """
        h_prime_coeffs = {
            "h1": self.n / 2 - (2 / 3) * self.n**2 + (37 / 96) * self.n**3 - (1 / 360) * self.n**4,
            "h2": (1 / 48) * self.n**2 + (1 / 15) * self.n**3 - (437 / 1440) * self.n**4,
            "h3": (17 / 480) * self.n**3 - (37 / 840) * self.n**4,
            "h4": (4397 / 161280) * self.n**4
        }

        eta_prime = (E - self.false_easting) / (self.B * self.scale_factor)
        if self.close_to_poles:
            xi_prime = (N - self.false_northing + self.scale_factor * self._calculate_meridional_arc_pole_safe(self.lat_origin)) / (self.B * self.scale_factor)
        else:
            xi_prime = (N - self.false_northing + self.scale_factor * self._calculate_meridional_arc()) / (self.B * self.scale_factor)

        xi0_prime = xi_prime
        eta0_prime = eta_prime

        for i in range(1, 5):
            xi0_prime -= h_prime_coeffs[f"h{i}"] * math.sin(2 * i * xi_prime) * math.cosh(2 * i * eta_prime)
            eta0_prime -= h_prime_coeffs[f"h{i}"] * math.cos(2 * i * xi_prime) * math.sinh(2 * i * eta_prime)

        beta_prime = math.asin(math.sin(xi0_prime) / math.cosh(eta0_prime))
        Q_prime = math.asinh(math.tan(beta_prime))

        e2 = self.f * (2 - self.f)
        e = math.sqrt(e2)
        Q_double_prime = Q_prime
        while True:
            Q_new = Q_prime + e * math.atanh(e * math.tanh(Q_double_prime))
            if abs(Q_new - Q_double_prime) < 1e-12:
                break
            Q_double_prime = Q_new

        lat = math.atan(math.sinh(Q_double_prime))
        lon = self.lon_origin + math.asin(math.tanh(eta0_prime) / math.cos(beta_prime))

        return lat, lon




if __name__ == "__main__":
    from pyproj import Transformer, CRS
    from TransversMercator import TransverseMercator
    proj_true_values = (555776.2668, 6651832.7354) # east, north

    # Define projection parameters
    lat_origin = 0  # Latitude of natural origin in radians
    lon_origin = math.radians(9)  # Longitude of natural origin in radians
    scale_factor = 0.9996  # Scale factor at the natural origin
    false_easting = 500000  # False easting in meters
    false_northing = 0  # False northing in meters

    # Ellipsoid parameters (WGS84)
    a = 6378137.0  # Semi-major axis in meters
    f = 1 / 298.257223563  # Flattening

    # Create an instance of the TransverseMercator class
    tm = TransverseMercator(lat_origin, lon_origin, scale_factor, false_easting, false_northing, a, f)

    # Example geographic coordinates (latitude, longitude in radians)
    lat = math.radians(60.0)  # Latitude in radians
    lon = math.radians(3.0)  # Longitude in radians

    # Perform forward projection
    easting, northing = tm.geog_to_projected(lat, lon)
    easting_true, northing_true = Transformer.from_crs("EPSG:4326", "EPSG:32632", always_xy=True).transform(lon,lat,radians= True)
    print(f"Projected Coordinates:\nEasting = {round(easting,4)}\nNorthing = {round(northing,4)}")
    print(f"\nProjected Coordinates Pyproj:\nEasting = {round(easting_true,4)}\nNorthing = {round(northing_true,4)}")

    # Perform inverse projection
    lat_back, lon_back = tm.projected_to_geog(easting, northing)
    print(f"Geographic Coordinates: Latitude = {math.degrees(lat_back)}, Longitude = {math.degrees(lon_back)}")








