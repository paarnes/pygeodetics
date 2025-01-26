"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com

MercatorVariantC class for performing Mercator Variant C projection conversions.
This class provides methods to initialize projection parameters, convert geographic coordinates to projected coordinates,
and convert projected coordinates back to geographic coordinates.
"""

import numpy as np
from typing import Literal, Union, List




class MercatorVariantC:
    """
    Initialize the parameters for the Mercator Variant C projection.

    Instance attributes:
    --------------------
    - a = semi-major axis of the ellipsoid
    - f = flattening of the ellipsoid
    - latSP1 = latitude of the first standard parallel (in radians)
    - latFO = latitude of false origin (in radians)
    - lonFO = longitude of false origin (in radians)
    - EFO = eastings at false origin
    - NFO = northings at false origin
    """

    proj_params_dict = {
        'Example_CRS_1' : {
            'a': 6378137,
            'f': 1/298.257223563,
            'Latitude of the first standard parallel': 0,
            'Latitude of false origin': 11,
            'Longitude of false origin': -58,
            'Easting at false origin': 0,
            'Northing at false origin': 0
            },
        'Example_CRS_2': {
            'a': 6378137,
            'f': 1/298.257223563,
            'Latitude of the first standard parallel': 0,
            'Latitude of false origin': 10,
            'Longitude of false origin': 117,
            'Easting at false origin': 3000000,
            'Northing at false origin': 2500000
            }
        }


    def __init__(self, proj_crs: Literal["Example_CRS_1", "Example_CRS_2"]=None,
                 a=None, f=None, latSP1=None, latFO=None, lonFO=None, EFO=None, NFO=None):
        """
        Initialize the parameters for the Mercator Variant C projection.
        """
        if proj_crs:
            params = MercatorVariantC.proj_params_dict[proj_crs]
            a = params['a']
            f = params['f']
            latSP1 = np.deg2rad(params['Latitude of the first standard parallel'])
            latFO = np.deg2rad(params['Latitude of false origin'])
            lonFO = np.deg2rad(params['Longitude of false origin'])
            EFO = params['Easting at false origin']
            NFO = params['Northing at false origin']

        self.a = a
        self.f = f
        self.latSP1 = abs(latSP1)  # Ensure the latitude is positive as per instructions
        self.latFO = latFO
        self.lonFO = lonFO
        self.EFO = EFO
        self.NFO = NFO

        # Derived parameters
        self.e = np.sqrt(2 * f - f ** 2)  # Eccentricity
        self.esq = self.e ** 2  # Square of eccentricity

    def calc_ko(self):
        """
        Calculate the projection constant ko based on the latitude of the first standard parallel.
        """
        return np.cos(self.latSP1) / ((1 - self.esq * np.sin(self.latSP1) ** 2) ** 0.5)

    def calc_M(self, lat, ko):
        """
        Calculate the projection constant M based on the latitude of the first standard parallel.
        """
        return self.a * ko * np.log(np.tan(np.pi/4 + lat/2) *
                                     ((1 - self.e * np.sin(lat)) / (1 + self.e * np.sin(lat))) ** (self.e/2))

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

    def geog_to_projected(self, coords, unit: Literal["deg", "rad"]="deg"):
        """
        Convert a list or array of geographic coordinates (longitude, latitude) to projected easting and northing coordinates.
        """
        coords = self.convert_to_np_array(coords)
        if  unit == "deg":
            coords = np.deg2rad(coords)
        ko = self.calc_ko()
        # Calculate Easting and Northing for each coordinate pair
        E = self.EFO + self.a * ko * (coords[:, 0] - self.lonFO)
        M_latFO = self.calc_M(self.latFO, ko)
        N = (self.NFO - M_latFO) + np.apply_along_axis(lambda lat: self.calc_M(lat, ko), 0, coords[:, 1])

        return np.column_stack((E, N))  # Combine Easting and Northing into one array

    def projected_to_geog(self, coords, unit: Literal["deg","rad"]="deg"):
        """
        Convert a list or array of projected coordinates (Easting, Northing) to geographic coordinates (latitude and longitude).
        """
        coords = self.convert_to_np_array(coords)
        ko = self.calc_ko()

        # Calculate longitude for each pair
        lon = ((coords[:, 0] - self.EFO) / (self.a * ko)) + self.lonFO

        # Calculate chi for each pair
        B = np.exp(1)
        t = B ** ((self.NFO - self.calc_M(self.latFO, ko) - coords[:, 1]) / (self.a * ko))
        chi = np.pi/2 - 2 * np.arctan(t)

        # Calculate latitude using iterative formula for each pair
        lat = chi + (self.esq/2 + 5*self.e**4/24 + self.e**6/12 + 13*self.e**8/360) * np.sin(2*chi) \
                  + (7*self.e**4/48 + 29*self.e**6/240 + 811*self.e**8/11520) * np.sin(4*chi) \
                  + (7*self.e**6/120 + 81*self.e**8/1120) * np.sin(6*chi) \
                  + (4279*self.e**8/161280) * np.sin(8*chi)

       # Combine longitude and latitude into one array
        geog_coords = np.column_stack((lon, lat))

        if unit == "deg":
            geog_coords = np.rad2deg(geog_coords)

        return geog_coords






if __name__=="__main__":

    # Example usage with provided parameters and values
    a = 6378245.00  # Semi-major axis
    f = 1/298.300  # Flattening
    latSP1 = np.deg2rad(42)  # Latitude of the first standard parallel
    latFO = np.deg2rad(42)  # Latitude of false origin
    lonFO = np.deg2rad(51)  # Longitude of false origin
    EF = 0.00  # Eastings at false origin
    NF = 0.00  # Northings at false origin

    converter = MercatorVariantC(a=a, f=f, latSP1=latSP1, latFO=latFO, lonFO=lonFO, EFO=EF, NFO=NF)
    lat = np.deg2rad(53)  # Latitude
    lon = np.deg2rad(53)  # Longitude
    coords = [[lon, lat]]
    proj_coords = converter.geog_to_projected(coords, unit="rad")
    reversed_geog = converter.projected_to_geog(proj_coords)

    E, N = proj_coords.flatten()
    reversed_lon, reversed_lat = reversed_geog.flatten()
    print(f"Easting E: {E:.2f} m")
    print(f"Northing N: {N:.2f} m")
    print(f"Reversed Latitude: {reversed_lat:.8f}°")
    print(f"Reversed Longitude: {reversed_lon:.8f}°")

    # #%% Example_CRS_1
    # a = 6378137
    # f = 1/298.257223563
    # latFO = np.deg2rad(11) # Latitude of false origin
    # latSP1 = np.deg2rad(0)  # Latitude of the first standard parallel
    # lonFO = np.deg2rad(-58)  # Longitude of false origin
    # EF = 0.00  # Eastings at false origin
    # NF = 0.00  # Northings at false origin



    # converter = MercatorVariantC(a=a, f=f, latSP1=latSP1, latFO=latFO, lonFO=lonFO, EFO=EF, NFO=NF)
    # lat = np.deg2rad(14)# Latitude
    # lon = np.deg2rad(-58)  # Longitude
    # coords = [[lon, lat]]
    # proj_coords = converter.geog_to_projected(coords, unit="rad")
    # reversed_geog = converter.projected_to_geog(proj_coords)

    # E, N = proj_coords.flatten()
    # reversed_lon, reversed_lat = reversed_geog.flatten()
    # print(f"Easting E: {E:.2f} m")
    # print(f"Northing N: {N:.2f} m")
    # print(f"Reversed Latitude: {reversed_lat:.8f}°")
    # print(f"Reversed Longitude: {reversed_lon:.8f}°")


    # #%%
    # converter = MercatorVariantC(proj_crs='Example_CRS_1')
    # lat = np.deg2rad(14)# Latitude
    # lon = np.deg2rad(-58)  # Longitude
    # coords = [[lon, lat]]
    # proj_coords = converter.geog_to_projected(coords, unit="rad")
    # reversed_geog = converter.projected_to_geog(proj_coords)

    # E, N = proj_coords.flatten()
    # reversed_lon, reversed_lat = reversed_geog.flatten()
    # print(f"Easting E: {E:.2f} m")
    # print(f"Northing N: {N:.2f} m")
    # print(f"Reversed Latitude: {reversed_lat:.8f}°")
    # print(f"Reversed Longitude: {reversed_lon:.8f}°")


    # #%%
    # converter = MercatorVariantC(proj_crs='Example_CRS_2')
    # lat = 10 # Latitude
    # lon = 115 # Longitude
    # coords = [[lon, lat]]
    # proj_coords = converter.geog_to_projected(coords)
    # reversed_geog = converter.projected_to_geog(proj_coords)

    # E, N = proj_coords.flatten()
    # reversed_lon, reversed_lat = reversed_geog.flatten()
    # print(f"Easting E: {E:.2f} m")
    # print(f"Northing N: {N:.2f} m")
    # print(f"Reversed Latitude: {reversed_lat:.8f}°")
    # print(f"Reversed Longitude: {reversed_lon:.8f}°")


    # #%%
    # converter = MercatorVariantC(proj_crs='Example_CRS_2')

    # # Convert list of lists
    # list_of_coords = np.array([[115, 10], [116, 11]])  # Each inner list is [longitude, latitude]
    # # list_of_coords = np.deg2rad(list_of_coords)
    # projected_coords = converter.geog_to_projected(list_of_coords)
    # print("Projected coordinates (list of lists):")
    # print(projected_coords)


    # # Convert list of lists
    # list_of_proj_coords = [[2777361.018, 2500000], [2888680.509, 2612483.937]]  # Each inner list is [Easting, Northing]
    # geog_coords = converter.projected_to_geog(list_of_proj_coords)
    # print("Geographic coordinates (list of lists):")
    # print(geog_coords)