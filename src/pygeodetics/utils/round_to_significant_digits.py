"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""

import numpy as np
from decimal import Decimal, getcontext

def round_to_significant_digits(value: float, sig_figs: int = 15) -> float:
    """
    Rounds a value to a specified number of significant figures.

    Parameters
    ----------
    value : float. The value to round.
    sig_figs : int, optional. The number of significant figures to round to (default is 15).

    Returns
    -------
    float: The value rounded to the specified number of significant figures.

    Raises
    ------
    ValueError: If sig_figs is not a positive integer.
    """
    if sig_figs <= 0:
        raise ValueError("Number of significant figures (sig_figs) must be a positive integer.")

    if value == 0:
        return 0.0

    # Use Decimal for precision
    getcontext().prec = sig_figs + 2  # Extra precision to ensure accuracy
    value_decimal = Decimal(value)
    magnitude = np.floor(np.log10(abs(value)))  # Magnitude of the value
    rounded_value = round(value_decimal, -int(magnitude - (sig_figs - 1)))

    return float(rounded_value)




if __name__ == "__main__":
    rounded_num = round_to_significant_digits(30.26388888888888888888889, 10)
    print(rounded_num)  # 30.26388889
