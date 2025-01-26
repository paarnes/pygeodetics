"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""


from typing import Tuple
import numpy as np


def tissots_indicatrix(h: float, k: float) -> Tuple[float, float, float]:
    """
    Compute Tissot's Indicatrix parameters: semimajor axis, semiminor axis, and angular distortion.

    Parameters
    ----------
    h : float
        Meridional scale factor (north-south direction).
    k : float
        Parallel scale factor (east-west direction).

    Returns
    -------
    tissot_a : float
        Semimajor axis of the indicatrix.
    tissot_b : float
        Semiminor axis of the indicatrix.
    omega : float
        Angular distortion (in radians).

    Examples
    --------
    >>> h, k = 1.0002, 0.9998
    >>> tissot_a, tissot_b, omega = tissots_indicatrix(h, k)
    >>> print(f"Tissot's Semimajor Axis: {tissot_a:.6f}")
    >>> print(f"Tissot's Semiminor Axis: {tissot_b:.6f}")
    >>> print(f"Angular Distortion: {np.degrees(omega):.6f} degrees")
    """

    # Semimajor axis (a) and semiminor axis (b) of Tissot's indicatrix
    tissot_a = np.sqrt((h**2 + k**2 + np.sqrt((h**2 - k**2)**2)) / 2)
    tissot_b = np.sqrt((h**2 + k**2 - np.sqrt((h**2 - k**2)**2)) / 2)

    # Angular distortion (omega)
    if h != 0 and k != 0:
        omega = np.arctan2(2 * h * k, h**2 - k**2) / 2
    else:
        omega = 0.0  # No distortion if one of the scales is zero

    return tissot_a, tissot_b, omega


if __name__ == "__main__":
    # Example usage
    h = 1.0002  # Example meridional scale factor
    k = 0.9998  # Example parallel scale factor

    tissot_a, tissot_b, omega = tissots_indicatrix(h, k)

    print(f"Tissot's Semimajor Axis: {tissot_a:.6f}")
    print(f"Tissot's Semiminor Axis: {tissot_b:.6f}")
    print(f"Angular Distortion: {np.degrees(omega):.6f} degrees")
