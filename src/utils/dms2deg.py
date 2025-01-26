"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""

from decimal import Decimal, getcontext
import numpy as np
from round_to_significant_digits import round_to_significant_digits

# Set precision for Decimal calculations
getcontext().prec = 25

def dms2deg(degrees: float, minutes: float, seconds: float, sig_figs: int = 15) -> Decimal:
    """
    Convert from degrees, minutes, and seconds (DMS) to decimal degrees (DD).

    Parameters
    ----------
    degrees : float. The degree component of the angle (integer or float). Positive for north/east, negative for south/west.
    minutes : float. The minute component of the angle. Must be between 0 and 59.
    seconds : float. The second component of the angle. Must be between 0 and 59.999... .
    sig_figs : int. Number of significant figures to round the result to (default is 15).

    Returns
    -------
    decimal_degrees : Decimal. The equivalent angle in decimal degrees.

    Raises
    ------
    ValueError
        If minutes or seconds are out of valid ranges.

    Examples
    --------
    >>> dms2deg(30, 15, 50)
    Decimal('30.26388888888889')
    >>> dms2deg(-45, 30, 0)
    Decimal('-45.5')
    """
    # Input validation
    if not (0 <= minutes < 60):
        raise ValueError("Minutes must be in the range [0, 60).")
    if not (0 <= seconds < 60):
        raise ValueError("Seconds must be in the range [0, 60).")

    # Calculate decimal degrees
    sign = np.sign(degrees)
    degrees = Decimal(abs(degrees))
    minutes = Decimal(minutes)
    seconds = Decimal(seconds)

    decimal_degrees = degrees + minutes / 60 + seconds / 3600
    decimal_degrees = round_to_significant_digits(decimal_degrees, sig_figs=sig_figs)
    return decimal_degrees * sign


# Example usage
if __name__ == "__main__":
    print(dms2deg(30, 15, 50))   # 30.26388888888889
    print(dms2deg(-45, 30, 0))  # -45.5
