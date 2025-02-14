"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""

import numpy as np

class ProjectionZoneCalc:
    """
    Class to calculate UTM projection zones based on longitude.
    """
    def __init__(self, initial_longitude=-180.0, zone_width=6.0, false_easting=500000.0, scale_factor=0.9996):
        """
        Initialize the projection zone calculator.

        Parameters
        ----------
        initial_longitude : float. Western limit of zone 1 in degrees. Default is -180.0.
        zone_width : float. Width of each zone in degrees. Default is 6.0 degrees
        false_easting : float. False easting value for the projection. Default is 500000.0.meters
        scale_factor : float. Scale factor at the natural origin. Default is 0.9996.
        """
        self.initial_longitude = initial_longitude
        self.zone_width = zone_width
        self.false_easting = false_easting
        self.scale_factor = scale_factor

    def calculate_zone(self, longitude):
        """
        Calculate the zone number based on the given longitude.

        Parameters
        ----------
        longitude : float or np.ndarray. Longitude(s) in degrees.

        Returns
        -------
        int or np.ndarray. Zone number(s).
        """
        zones = np.floor((longitude - self.initial_longitude) / self.zone_width) % (360 / self.zone_width) + 1
        return zones.astype(int)

    def calculate_central_meridian(self, zone):
        """
        Calculate the central meridian for a given zone.

        Parameters
        ----------
        zone : int. Zone number.

        Returns
        -------
        float. Central meridian for the given zone in degrees.
        """
        return self.initial_longitude + zone * self.zone_width - self.zone_width / 2



if __name__ == "__main__":

    # Create an instance of the ProjectionZoneCalc
    calculator = ProjectionZoneCalc()

    # Example longitude
    longitude = 10.9857
    # Calculate zone
    zone = calculator.calculate_zone(longitude)
    print(f"Zone for longitude {longitude}: {zone}")

    # Calculate central meridian
    central_meridian = calculator.calculate_central_meridian(zone)
    print(f"Central meridian for zone {zone}: {central_meridian}")

