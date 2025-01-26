"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""

from decimal import Decimal, getcontext
from round_to_significant_digits import round_to_significant_digits

# Set precision for Decimal operations
getcontext().prec = 15

def deg2dms(deg: float, sig_figs: int = 9) -> tuple:
    """
    Convert from decimal degrees to degrees, minutes, and seconds (DMS).

    Parameters
    ----------
    deg : float. Angle in decimal degrees.
    sig_figs : int. Number of significant figures to round seconds to (default is 15).

    Returns
    -------
    A tuple (degrees, minutes, seconds) where:
        - degrees (int): Integer degrees.
        - minutes (int): Integer minutes (0 to 59).
        - seconds (float): Fractional seconds (0 to 59.999...).
    """
    # Use Decimal for precise calculations
    deg = Decimal(deg)

    # Extract sign for handling negative values
    sign = -1 if deg < 0 else 1

    # Take absolute value for calculation
    abs_deg = abs(deg)

    # Calculate degrees, minutes, and seconds
    degrees = int(abs_deg)  # Integer part of degrees
    fractional_degrees = abs_deg - degrees
    minutes = int(fractional_degrees * 60)
    seconds = (fractional_degrees * 60 - minutes) * 60

    # Apply the sign to the degrees
    degrees *= sign
    seconds = round_to_significant_digits(seconds, sig_figs=sig_figs)
    return degrees, minutes, seconds


if __name__ == "__main__":
    print(deg2dms(30.2638888889))  # (30, 15, 50.0)
    print(deg2dms(-45.5))        # (-45, 30, 0.0)
    print(deg2dms(0.0))          # (0, 0, 0.0)
