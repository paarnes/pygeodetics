"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""

import numpy as np
import pytest
import Ellipsoid as Ellipsoid
from Geodetic import Geodetic

# Define WGS84 ellipsoid parameters
ellip = Ellipsoid.WGS84()
a = ellip.a
b = ellip.b
geod = Geodetic(ellip)

# Test cases for Vincenty distance function
vincenty_test_cases = [
    {
        "lat1": 1.0,
        "lon1": 30.0,
        "lat2": -1.0,
        "lon2": 30.0,
        "expected_distance": 221148.7771,
        "description": "Crossing the equator",
    },
    {
        "lat1": 10.0,
        "lon1": 179.9,
        "lat2": 10.0,
        "lon2": -179.9,
        "expected_distance": 21927.8725,
        "description": "Crossing the anti-meridian",
    },
    {
        "lat1": 89.9,
        "lon1": 45.0,
        "lat2": 89.9,
        "lon2": -135.0,
        "expected_distance": 22338.7957,
        "description": "Over the North Pole",
    },
    {
        "lat1": -89.9,
        "lon1": 45.0,
        "lat2": -89.9,
        "lon2": -135.0,
        "expected_distance": 22338.7957,
        "description": "Over the South Pole",
    },
    {
        "lat1": 52.2296756,
        "lon1": 21.0122287,
        "lat2": 41.8919300,
        "lon2": 12.5113300,
        "expected_distance": 1316208.0833,
        "description": "Warsaw to Rome",
    },
]

# Test cases for geodetic inverse problem
inverse_test_cases = [
    {
        "lat1": 52.2296756,
        "lon1": 21.0122287,
        "lat2": 41.8919300,
        "lon2": 12.5113300,
        "quadrant_correction": False,
        "radians": False,
        "expected_az1": -147.4628043168,
        "expected_az2": -153.7168672619,
        "expected_distance": 1316208.0833,
        "description": "Geodetic inverse problem 1",
    },
]

# Test cases for geodetic direct problem
direct_test_cases = [
    {
        "lat1": 52.2296756,
        "lon1": 21.0122287,
        "az1": -147.4628,
        "d": 1316208.08334,
        "quadrant_correction": False,
        "radians": False,
        "expected_lat2": 41.8919300,
        "expected_lon2": 12.5113300,
        "expected_az2": -153.7168672618,
        "description": "Geodetic direct problem 1",
    },
]

# Test cases for radius of curvature
radius_test_cases = [
    {
        "lat": 45.0,
        "az": 30.0,
        "expected": 6346070.049,
        "description": "Radius at 45° latitude",
    },
]

# Test cases for geodetic to ECEF conversion
geod2ecef_test_cases = [
    {
        "lat": 59.9070724743,
        "lon": 10.754482924,
        "h": 63.8281,
        "expected_X": 3149785.9652,
        "expected_Y": 598260.8822,
        "expected_Z": 5495348.4927,
        "description": "ECEF Conversion Test",
    },
]

# Test cases for distance between two points
distance_test_cases = [
    {
        "lat1": 52.2296756,
        "lon1": 21.0122287,
        "lat2": 41.8919300,
        "lon2": 12.5113300,
        "expected_distance": 1316208.08334,
        "description": "Warsaw to Rome",
    },
]


# Test cases for ECEF to ENU conversion
ecef2enu_test_cases = [
    {
        "X": 3149785.9652,
        "Y": 598260.8822,
        "Z": 5495348.4927,
        "lat0": 59.9070724743,
        "lon0": 10.754482924,
        "h0": 63.8281,
        "expected_e": 0.0,
        "expected_n": 0.0,
        "expected_u": 0.0,
        "description": "ECEF to ENU",
    },
]



# Test cases for nrad method (normal radius of curvature)
nrad_test_cases = [
    {
        "lat": 60,
        "radians": False,
        "expected_nrad": 6394209.1738,
        "description": "Test 1: Normal radius of curvature at 60° latitude",
    },
    {
        "lat": 45,
        "radians": False,
        "expected_nrad": 6388838.2901,
        "description": "Test 2: Normal radius of curvature at 45° latitude",
    },
    {
        "lat": 0,
        "radians": False,
        "expected_nrad": 6378137.0000,
        "description": "Test 3: Normal radius of curvature at 0° latitude",
    },
]


# Test cases for Mrad method (meridional radius of curvature)
mrad_test_cases = [
    {
        "lat": 60,
        "radians": False,
        "expected_mrad": 6383453.85723,
        "description": "Test 1: Meridional radius of curvature at 60° latitude",
    },
    {
        "lat": 45,
        "radians": False,
        "expected_mrad": 6367381.81562,
        "description": "Test 2: Meridional radius of curvature at 45° latitude",
    },
    {
        "lat": 0,
        "radians": False,
        "expected_mrad": 6335439.32729,
        "description": "Test 3: Meridional radius of curvature at 0° latitude",
    },
]





# Run the tests
@pytest.mark.parametrize("case", vincenty_test_cases)
def test_vincenty_degrees(case):
    """
    Test the custom vincenty_distance function against pygeodesy's implementation.
    """
    lat1_rad, lon1_rad, lat2_rad, lon2_rad = [
        case["lat1"],
        case["lon1"],
        case["lat2"],
        case["lon2"],
    ]
    custom_distance = geod.distance_between_two_points(
        lon1_rad, lat1_rad, lon2_rad, lat2_rad, radians=False
    )
    assert np.isclose(custom_distance, case["expected_distance"], atol=1e-6)

@pytest.mark.parametrize("case", vincenty_test_cases)
def test_vincenty_radians(case):
    """
    Test the custom vincenty_distance function against pygeodesy's implementation.
    """
    lat1_rad, lon1_rad, lat2_rad, lon2_rad = map(
        np.radians, [case["lat1"], case["lon1"], case["lat2"], case["lon2"]]
    )
    custom_distance = geod.distance_between_two_points(
        lon1_rad, lat1_rad, lon2_rad, lat2_rad, radians=True
    )
    assert np.isclose(custom_distance, case["expected_distance"], atol=1e-6)

@pytest.mark.parametrize("case", inverse_test_cases)
def test_geodetic_inverse_problem(case):
    """
    Test geodetic inverse problem method.
    """
    lat1, lon1, lat2, lon2 = case["lat1"], case["lon1"], case["lat2"], case["lon2"]
    az1, az2, distance = geod.inverse_problem(
        lat1, lon1, lat2, lon2, quadrant_correction=case["quadrant_correction"], radians=case["radians"]
    )

    assert np.isclose(az1, case["expected_az1"], atol=1e-2)
    assert np.isclose(az2, case["expected_az2"], atol=1e-2)
    assert np.isclose(distance, case["expected_distance"], atol=1e-2)


@pytest.mark.parametrize("case", inverse_test_cases)
def test_geodetic_inverse_problem_radians(case):
    """
    Test geodetic inverse problem method.
    """
    lat1, lon1, lat2, lon2 = map(
        np.radians, [case["lat1"], case["lon1"], case["lat2"], case["lon2"]]
    )
    az1, az2, distance = geod.inverse_problem(
        lat1, lon1, lat2, lon2, quadrant_correction=case["quadrant_correction"], radians=True
    )

    assert np.isclose(np.degrees(az1), case["expected_az1"], atol=1e-2)
    assert np.isclose(np.degrees(az2), case["expected_az2"], atol=1e-2)
    assert np.isclose(distance, case["expected_distance"], atol=1e-2)

@pytest.mark.parametrize("case", direct_test_cases)
def test_geodetic_direct_problem(case):
    """
    Test geodetic direct problem method.
    """
    lat1, lon1, az1, d = case["lat1"], case["lon1"], case["az1"], case["d"]
    lat2, lon2, az2 = geod.direct_problem(
        lat1, lon1, az1, d, quadrant_correction=case["quadrant_correction"], radians=case["radians"]
    )

    assert np.isclose(lat2, case["expected_lat2"], atol=1e-4)
    assert np.isclose(lon2, case["expected_lon2"], atol=1e-4)
    assert np.isclose(az2, case["expected_az2"], atol=1e-4)

@pytest.mark.parametrize("case", direct_test_cases)
def test_geodetic_direct_problem_radians(case):
    """
    Test geodetic direct problem method.
    """
    lat1, lon1, az1 = map(np.radians, [case["lat1"], case["lon1"], case["az1"]])
    lat2, lon2, az2 = geod.direct_problem(
        lat1, lon1, az1, case["d"], quadrant_correction=case["quadrant_correction"], radians=True
    )

    assert np.isclose(np.degrees(lat2), case["expected_lat2"], atol=1e-4)
    assert np.isclose(np.degrees(lon2), case["expected_lon2"], atol=1e-4)
    assert np.isclose(np.degrees(az2), case["expected_az2"], atol=1e-4)


@pytest.mark.parametrize("case", radius_test_cases)
def test_radius_of_curvature(case):
    """
    Test radius of curvature computation.
    """
    radius = geod.radius_of_curvature(case["lat"], case["az"], radians=False)
    assert np.isclose(radius, case["expected"], atol=1e-2)


@pytest.mark.parametrize("case", geod2ecef_test_cases)
def test_geod2ecef(case):
    """
    Test geodetic to ECEF conversion.
    """
    X, Y, Z = geod.geod2ecef(case["lat"], case["lon"], case["h"], radians=False)

    assert np.isclose(X, case["expected_X"], atol=1e-2)
    assert np.isclose(Y, case["expected_Y"], atol=1e-2)
    assert np.isclose(Z, case["expected_Z"], atol=1e-2)


@pytest.mark.parametrize("case", distance_test_cases)
def test_distance_between_two_points(case):
    """
    Test distance between two points computation.
    """
    lat1, lon1, lat2, lon2 = map(
        np.radians, [case["lat1"], case["lon1"], case["lat2"], case["lon2"]]
    )
    distance = geod.distance_between_two_points(lon1, lat1, lon2, lat2, radians=True)

    assert np.isclose(distance, case["expected_distance"], atol=1e-2)


@pytest.mark.parametrize("case", ecef2enu_test_cases)
def test_ecef2enu(case):
    """
    Test ECEF to ENU conversion.
    """
    e, n, u = geod.ecef2enu(
        case["X"],
        case["Y"],
        case["Z"],
        case["lat0"],
        case["lon0"],
        case["h0"],
        radians=False,
    )

    assert np.isclose(e, case["expected_e"], atol=1e-4)
    assert np.isclose(n, case["expected_n"], atol=1e-4)
    assert np.isclose(u, case["expected_u"], atol=1e-4)


@pytest.mark.parametrize("case", nrad_test_cases)
def test_nrad(case):
    """
    Test normal radius of curvature computation.
    """
    nrad = geod.nrad(case["lat"], radians=case["radians"])
    assert np.isclose(nrad, case["expected_nrad"], atol=1e-6), (
        f"Test failed for case: {case['description']}\n"
        f"Computed Nrad: {nrad}\n"
        f"Expected Nrad: {case['expected_nrad']}"
    )

@pytest.mark.parametrize("case", mrad_test_cases)
def test_mrad(case):
    """
    Test meridional radius of curvature computation.
    """
    nrad = geod.mrad(case["lat"], radians=case["radians"])
    assert np.isclose(nrad, case["expected_mrad"], atol=1e-6), (
        f"Test failed for case: {case['description']}\n"
        f"Computed Mrad: {nrad}\n"
        f"Expected Mrad: {case['expected_mrad']}"
    )




if __name__ == "__main__":
    pytest.main()


