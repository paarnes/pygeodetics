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

    # Example geographic coordinates (longitude, latidtude in degrees) (Order: [[x1,y2], [x2,y2]...[xn,yn]])
    coords = np.array([[3.44958679, 60.83740565], [3.52343265, 61.43578934]])

    # Get projected coordinates from the latitude and longitude
    easting, northing = tm.geog_to_projected(coords, unit="deg")

"""

import numpy as np
from typing import Tuple, Literal, Union, List


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

    Note:
    -----
    Order of input coordinates: First longitude (x), then latitude (y).
        - [[x1,y2], [x2,y2]...[xn,yn]]
        - [[lon1,lat1], [lon2,lat2]...[lonn,latn]]
        - [[east1,north1], [east2,north2]...[eastn,northn]]

    Example:
    --------
    .. code-block:: python

            # Define projection parameters
            lat_origin = 0              # Latitude of natural origin in radians
            lon_origin = np.radians(9)  # Longitude of natural origin in radians
            scale_factor = 0.9996       # Scale factor at the natural origin
            false_easting = 500000      # False easting in meters
            false_northing = 0          # False northing in meters

            # Ellipsoid parameters (WGS84)
            a = 6378137.0               # Semi-major axis in meters
            f = 1 / 298.257223563       # Flattening

            # Create an instance of the TransverseMercator class
            tm = TransverseMercator(lat_origin, lon_origin, scale_factor, false_easting, false_northing, a, f)

            # Example geographic coordinates (longitude, latidtude in degrees) (Order: [[x1,y2], [x2,y2]...[xn,yn]])
            coords = np.array([[3.44958679, 60.83740565], [3.52343265, 61.43578934]])

            # Get projected coordinates from the latitude and longitude
            easting, northing = tm.geog_to_projected(coords, unit="deg")


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
        if np.isclose(self.lat_origin, 0, atol=1e-8):
            return False
        elif np.isclose(self.lat_origin, np.pi / 2, atol=1e-8) or np.isclose(self.lat_origin, -np.pi / 2, atol=1e-8):
            return False
        if np.abs(self.lat_origin - np.pi / 2) < 2 / 3600 * np.pi / 180:
            return True
        return False

    def _calculate_meridional_arc_pole_safe(self, lat) -> float:
        """
        Calculate the meridional arc distance using a series expansion.
        """
        e2 = self.f * (2 - self.f)
        a0 = 1 - e2 / 4 - 3 * e2**2 / 64 - 5 * e2**3 / 256
        a2 = (3 / 8) * (e2 + e2**2 / 4 + 15 * e2**3 / 128)
        a4 = (15 / 256) * (e2**2 + 3 * e2**3 / 4)
        a6 = (35 / 3072) * e2**3

        M = self.a * (a0 * lat - a2 * np.sin(2 * lat) + a4 * np.sin(4 * lat) - a6 * np.sin(6 * lat))
        return M

    def _calculate_meridional_arc(self) -> float:
        """
        Calculate the meridional arc distance from the equator to the origin latitude
        using the alternative method for latitudes close to the poles.

        Returns:
        - M0 (float): Meridional arc distance in meters.
        """
        # Handle special cases for latitude of origin
        if np.abs(self.lat_origin) < 1e-8:
            return 0
        elif np.isclose(self.lat_origin, np.pi / 2, atol=1e-8):
            return self.B * (np.pi / 2)
        elif np.isclose(self.lat_origin, -np.pi / 2, atol=1e-8):
            return self.B * (-np.pi / 2)

        # General case for latitude of origin
        e2 = self.f * (2 - self.f)  # Eccentricity squared
        e = np.sqrt(e2)  # Eccentricity

        # Compute Q0
        Qo = np.arcsinh(np.tan(self.lat_origin)) - e * np.arctanh(e * np.sin(self.lat_origin))
        beta_o = np.arctan(np.sinh(Qo))

        # Compute xi_o0
        xi_o0 = np.arcsin(np.sin(beta_o))

        # Compute xi_o using the series expansion
        xi_o1 = self.h_coeffs[f"h1"] * np.sin(2 * xi_o0)
        xi_o2 = self.h_coeffs[f"h2"] * np.sin(4 * xi_o0)
        xi_o3 = self.h_coeffs[f"h3"] * np.sin(6 * xi_o0)
        xi_o4 = self.h_coeffs[f"h4"] * np.sin(8 * xi_o0)
        xi_o = xi_o0 + xi_o1 + xi_o2 + xi_o3 + xi_o4

        # Compute M0
        M0 = self.B * xi_o

        return M0

    def calculate_xi_eta(self, eta0, xi0) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute full xi and eta using the series expansion for an array of inputs.

        Args:
        - eta0 (np.ndarray): Initial eta values.
        - xi0 (np.ndarray): Initial xi values.

        Returns:
        - Tuple[np.ndarray, np.ndarray]: Full xi and eta values.
        """
        xi = xi0.copy()
        eta = eta0.copy()
        for i in range(1, 5):
            xi += self.h_coeffs[f"h{i}"] * np.sin(2 * i * xi0) * np.cosh(2 * i * eta0)
            eta += self.h_coeffs[f"h{i}"] * np.cos(2 * i * xi0) * np.sinh(2 * i * eta0)
        return xi, eta

    def convert_to_np_array(self, coords: Union[List[List[float]], np.ndarray]):
        """
        Convert input list of lists or NumPy array to a NumPy array if not already one.
        """
        if isinstance(coords, list):
            return np.array(coords)
        elif isinstance(coords, np.ndarray):
            return coords
        else:
            raise TypeError("Input coordinates must be a list of lists or a NumPy array.")

    def geog_to_projected(self, coordinates: np.ndarray, unit: Literal["deg", "rad"]="deg") -> np.ndarray:
        """
        Convert geographic coordinates (longitude, latitude) to projected coordinates (easting, northing).

        Note: The order must be first longitude (x), then latitude (y).

        Args:
        - coordinates (np.ndarray): Array of shape (n, 2) or (n, 3) containing longitude and latitude in degrees or radians.
        - unit (str): The unit of the input coordinates. Either "deg" or "rad".

        Returns:
        - np.ndarray: Array of shape (n, 2) containing Easting (E) and Northing (N).
        """
        coordinates = self.convert_to_np_array(coordinates)

        if coordinates.shape[1] not in [2, 3]:
            raise ValueError("Input coordinates must have shape (n, 2) or (n, 3)")

        lon, lat = coordinates[:, 0], coordinates[:, 1]
        height = coordinates[:, 2] if coordinates.shape[1] == 3 else None

        if unit == "deg":
            lat = np.radians(lat)
            lon = np.radians(lon)

        # First, calculate Q and beta
        e2 = self.f * (2 - self.f)
        e = np.sqrt(e2)
        Q = np.arcsinh(np.tan(lat)) - e * np.arctanh(e * np.sin(lat))
        beta = np.arctan(np.sinh(Q))

        # Compute initial xi0 and eta0
        eta0 = np.arctanh(np.cos(beta) * np.sin(lon - self.lon_origin))
        xi0 = np.arcsin(np.sin(beta) * np.cosh(eta0))

        # Compute full xi and eta
        xi, eta = self.calculate_xi_eta(eta0, xi0)

        # Compute M0 (Meridional Arc at Origin)
        if self.close_to_poles:
            M0 = self._calculate_meridional_arc_pole_safe(self.lat_origin)
        else:
            M0 = self._calculate_meridional_arc()

        # Compute Easting and Northing
        E = self.false_easting + self.scale_factor * self.B * eta
        N = self.false_northing + self.scale_factor * (self.B * xi - M0)

        # Include height if present
        if height is not None:
            return np.vstack((E, N, height))
        return np.vstack((E, N))

    def projected_to_geog(self, proj_coordinates: np.ndarray, unit: Literal["deg", "rad"] = "deg") -> np.ndarray:
        """
        Convert projected coordinates (easting, northing) to geographic coordinates (longitude, latitude).

        Args:
        - proj_coordinates (np.ndarray): Array of shape (n, 2) containing easting and northing coordinates.
        - unit (str): The desired output unit. Either "deg" (degrees) or "rad" (radians).

        Returns:
        - np.ndarray: Array of shape (n, 2) containing longitude and latitude in the specified unit.
        """
        proj_coordinates = self.convert_to_np_array(proj_coordinates)

        if proj_coordinates.shape[1] not in [2, 3]:
            raise ValueError("Input coordinates must have shape (n, 2) or (n, 3)")

        E, N = proj_coordinates[:, 0], proj_coordinates[:, 1]
        height = proj_coordinates[:, 2] if proj_coordinates.shape[1] == 3 else None

        # Define reverse series coefficients
        h_prime_coeffs = {
            "h1": self.n / 2 - (2 / 3) * self.n**2 + (37 / 96) * self.n**3 - (1 / 360) * self.n**4,
            "h2": (1 / 48) * self.n**2 + (1 / 15) * self.n**3 - (437 / 1440) * self.n**4,
            "h3": (17 / 480) * self.n**3 - (37 / 840) * self.n**4,
            "h4": (4397 / 161280) * self.n**4
        }

        # Compute eta_prime and xi_prime
        eta_prime = (E - self.false_easting) / (self.B * self.scale_factor)
        if self.close_to_poles:
            xi_prime = (N - self.false_northing + self.scale_factor * self._calculate_meridional_arc_pole_safe(self.lat_origin)) / (self.B * self.scale_factor)
        else:
            xi_prime = (N - self.false_northing + self.scale_factor * self._calculate_meridional_arc()) / (self.B * self.scale_factor)

        # Backward series expansion for xi0_prime and eta0_prime
        xi0_prime = xi_prime.copy()
        eta0_prime = eta_prime.copy()
        for i in range(1, 5):
            xi0_prime -= h_prime_coeffs[f"h{i}"] * np.sin(2 * i * xi_prime) * np.cosh(2 * i * eta_prime)
            eta0_prime -= h_prime_coeffs[f"h{i}"] * np.cos(2 * i * xi_prime) * np.sinh(2 * i * eta_prime)

        # Compute beta_prime and Q_prime
        beta_prime = np.arcsin(np.sin(xi0_prime) / np.cosh(eta0_prime))
        Q_prime = np.arcsinh(np.tan(beta_prime))

        # Iteratively compute Q_double_prime for latitude
        e2 = self.f * (2 - self.f)
        e = np.sqrt(e2)
        Q_double_prime = Q_prime.copy()
        while True:
            Q_new = Q_prime + e * np.arctanh(e * np.tanh(Q_double_prime))
            if np.max(np.abs(Q_new - Q_double_prime)) < 1e-12:
                break
            Q_double_prime = Q_new

        # Compute latitude and longitude
        lat = np.arctan(np.sinh(Q_double_prime))
        lon = self.lon_origin + np.arcsin(np.tanh(eta0_prime) / np.cos(beta_prime))

        # Combine into a single array and convert units if necessary
        geog_coords = np.vstack((lon, lat))  # Longitude first, then latitude
        if unit == "deg":
            geog_coords = np.degrees(geog_coords)

        # Include height if present
        if height is not None:
            return np.vstack((geog_coords, height))

        return geog_coords



if __name__ == "__main__":
    from pyproj import Transformer, CRS
    np.set_printoptions(precision=8, suppress=True)
    proj_true_values = (555776.2668, 6651832.7354)  # east, north

    # Define projection parameters
    lat_origin = 0  # Latitude of natural origin in radians
    lon_origin = np.radians(9)  # Longitude of natural origin in radians (zone 32)
    # lon_origin = np.radians(15)  # Longitude of natural origin in radians (zone 33)
    scale_factor = 0.9996  # Scale factor at the natural origin
    false_easting = 500000  # False easting in meters
    false_northing = 0  # False northing in meters

    # Ellipsoid parameters (WGS84)
    a = 6378137.0  # Semi-major axis in meters
    f = 1 / 298.257223563  # Flattening

    # Create an instance of the TransverseMercator class
    tm = TransverseMercator(lat_origin, lon_origin, scale_factor, false_easting, false_northing, a, f)

    # Example geographic coordinates (latitude, longitude in radians)
    lat = np.radians([60.0, 60.1, 60.2])  # Latitude in radians
    lon = np.radians([3.0, 3.1, 3.2])  # Longitude in radians

    coords = np.array([[3, 60], [3.2, 61]])
    coords = np.array([[3, 60, 100], [3.2, 61, 102]])
    # coords = [[3, 60], [3.2, 61]] # list of lists

    # Perform forward projection
    easting, northing, *height = tm.geog_to_projected(coords, unit="deg")
    results = np.vstack((easting, northing, height)).T if height else np.vstack((easting, northing)).T
    print(f"\nProjected Coordinates TM class:\n{results}")

    # # # Perform inverse projection
    proj_coordinates = np.vstack([easting, northing, height]).T if height else np.vstack([easting, northing]).T
    lon_back, lat_back, *height = tm.projected_to_geog(proj_coordinates, unit="deg")
    results = np.vstack((lon_back, lat_back, height)).T if height else np.vstack((lon_back, lat_back)).T
    print(f"\nGeographic Coordinates TM class:\n{results}")

    # # Test using Pyproj
    easting_true, northing_true, *height_true = Transformer.from_crs("EPSG:4326", "EPSG:32632", always_xy=True).transform(coords[:, 0], coords[:, 1], coords[:, 1])
    results_pyproj = np.vstack((easting_true, northing_true, height_true)).T if height else np.vstack((easting_true, northing_true)).T
    print(f"\nProjected Coordinates Pyproj:\n{results_pyproj}")










