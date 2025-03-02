"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""


import sys
import os
import numpy as np
import pytest
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
from geodetics.ECEF2enu import ECEF2enu

# Define test cases with expected results
test_cases = [
    {
        "X": 12738186.63827794,
        "Y": -15447555.301322976,
        "Z": 10385003.518329535,
        "lat0": np.radians(60.0),
        "lon0": np.radians(5.0),
        "h0": 0,
        "radians": True,
        "e_true": -16498978.807374395,
        "n_true": -4612610.325871408,
        "u_true": 8303257.059092532,
        "description": "Case 1: ECEF to ENU using radians"
    },
    {
        "X": 3149785.9652,
        "Y": 598260.8822,
        "Z": 5495348.4927,
        "lat0": 59.907072474276958,
        "lon0": 10.754482924017791,
        "h0": 63.8281,
        "radians": False,
        "e_true": 0,
        "n_true": 0,
        "u_true": 0,
        "description": "Case 2: ECEF to ENU using degrees"
    },
   {
        "X": 4438287.3531,
        "Y": 842994.9643,
        "Z": -4487393.5421,
        "lat0": 49.907072474276958,
        "lon0": 5.754482924017791,
        "h0": 10.0,
        "radians": False,
        "e_true": 393737.9220936171,
        "n_true": -6311779.278360579,
        "u_true": -6900082.913685788,
        "description": "Case 2: ECEF to ENU using degrees"
    },
   {
        "X": -2949198.730,
        "Y": -3422172.1149,
        "Z": -4487393.5421,
        "lat0": -49.907072474276958,
        "lon0": -115.754482924017791,
        "h0": 10.0,
        "radians": False,
        "e_true": -1169250.2387980018,
        "n_true": 427145.20298326924,
        "u_true": -122429.64300887837,
        "description": "Case 2: ECEF to ENU using degrees"
    }
]


@pytest.mark.parametrize("case", test_cases)
def test_ECEF2enu(case):
    """
    Test the ECEF2enu function using various test cases.
    """

    X, Y, Z = case["X"], case["Y"], case["Z"]
    lat0, lon0, h0  = case["lat0"], case["lon0"], case["h0"]
    radians = case["radians"]
    e_true, n_true, u_true = case["e_true"], case["n_true"], case["u_true"]
    description = case["description"]

    # Compute ENU coordinates
    e, n, u = ECEF2enu(X, Y, Z, lat0, lon0, h0, radians=radians)

    # Check results
    assert np.isclose(e, e_true, atol=1e-4), (
        f"Test failed for {description}\nComputed East: {e} | Expected East: {e_true}"
    )
    assert np.isclose(n, n_true, atol=1e-4), (
        f"Test failed for {description}\nComputed North: {n} | Expected North: {n_true}"
    )
    assert np.isclose(u, u_true, atol=1e-4), (
        f"Test failed for {description}\nComputed Up: {u} | Expected Up: {u_true}"
    )

