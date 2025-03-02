
<p align="center">
    <img src="docs/icon/pygeodetics.png" alt="PyGeodetics Logo" width="700">
</p>


## Introduction
PyGeodetics is a Python library for performing geodetic computations like geodetic inverse and direct problems, conversions between different reference systems like ECEF to ENU, ECEF to geographic etc.

## Features
- Convert geodetic coordinates (latitude, longitude, height) to ECEF (Earth-Centered, Earth-Fixed)
- Convert ECEF to geodetic coordinates
- Transform ECEF to local ENU (East-North-Up) and NED (North-East-Down) systems
- Solve geodetic inverse and direct problems
- Distance between two points along the ellipsoid
- Compute radius of curvature and mean radius of the reference ellipsoid
- Support for different reference ellipsoids

## Installation
```sh
pip install pygeodetics
```

## TOC Code examples
- [Geodetic to ECEF](#geodetic-to-ecef)
- [ECEF to Geodetic](#ecef-to-geodetic)
- [ECEF to ENU](#ecef-to-enu)
- [ECEF to NED](#ecef-to-ned)
- [Geodetic Inverse Problem on the GRS80 ellipsoid](#geodetic-inverse-problem-on-the-grs80-ellipsoid)
- [Geodetic Direct Problem on the GRS80 ellipsoid](#geodetic-direct-problem-on-the-grs80-ellipsoid)
- [Radius of Curvature for a given Azimuth using Euler's equation.](#radius-of-curvature-for-a-given-azimuth-using-eulers-equation)
- [Calculate the mean Radius of the International1924 Ellipsoid](#calculate-the-mean-radius-of-the-international1924-ellipsoid)
- [Calculate the distance between two points on the ellipsoid (Vincenty formula)](#calculate-the-distance-between-two-points-on-the-ellipsoid-vincenty-formula)
- [Calculate the meridional radius of curvature (M) at a given latitude](#calculate-the-meridional-radius-of-curvature-m-at-a-given-latitude)
- [Calculate the normal radius of curvature (N) at a given latitude.](#calculate-the-normal-radius-of-curvature-n-at-a-given-latitude)

## Usage Examples

### Geodetic to ECEF
```python
from pygeodetics import Geodetic

# Initialize Geodetic class with WGS84 ellipsoid
geod = Geodetic()

lat = 59.907072474276958  # Latitude in degrees
lon = 10.754482924017791  # Longitude in degrees
h = 63.8281  # Height in meters

X, Y, Z = geod.geod2ecef(lat, lon, h)
print(f"Geodetic to ECEF:\nX: {X:.4f} m\nY: {Y:.4f} m\nZ: {Z:.4f} m")

```

### ECEF to Geodetic
```python
from pygeodetics import Geodetic

X, Y, Z = 3149785.9652, 598260.8822, 5495348.4927
geod = Geodetic()
lat, lon, h = geod.ecef2geod(X, Y, Z, angle_unit='deg')
print(f"ECEF to Geodetic:\nLatitude: {lat:.6f}°\nLongitude: {lon:.6f}°\nHeight: {h:.3f} m")

```

### ECEF to ENU
```python
from pygeodetics import Geodetic

X, Y, Z = 3149785.9652, 598260.8822, 5495348.4927
lat0, lon0, h0 = 58.907072, 9.75448, 63.8281

e, n, u = Geodetic().ecef2enu(X, Y, Z, lat0, lon0, h0, radians=False)
print(f"ECEF to ENU:\nEast: {e:.6f} m\nNorth: {n:.6f} m\nUp: {u:.6f} m")

```

### ECEF to NED
```python
from pygeodetics import Geodetic

X, Y, Z = 3149785.9652, 598260.8822, 5495348.4927
lat0, lon0, h0 = 58.907072, 9.75448, 63.8281

n, e, d = Geodetic().ecef2ned(X, Y, Z, lat0, lon0, h0)
print(f"ECEF to NED:\nNorth: {n:.6f} m\nEast: {e:.6f} m\nDown: {d:.6f} m")

```

### Geodetic Inverse Problem on the GRS80 ellipsoid
```python
from pygeodetics import Geodetic
from pygeodetics.Ellipsoid import GRS80

geod = Geodetic(GRS80())

lat1, lon1 = 52.2296756, 21.0122287
lat2, lon2 = 41.8919300, 12.5113300

az1, az2, distance = geod.inverse_problem(lat1, lon1, lat2, lon2, quadrant_correction=False)
print(f"Geodetic Inverse Problem:\nForward Azimuth: {az1:.6f}°\nReverse Azimuth: {az2:.6f}°\nDistance: {distance:.3f} m")

```

### Geodetic Direct Problem on the GRS80 ellipsoid
```python
from pygeodetics import Geodetic
from pygeodetics.Ellipsoid import GRS80

geod = Geodetic(GRS80())

az1 = -147.4628043168
d = 1316208.08334

lat2, lon2, az2 = geod.direct_problem(lat1, lon1, az1, d, quadrant_correction=True)
print(f"Geodetic Direct Problem:\nDestination Latitude: {lat2:.6f}°\nDestination Longitude: {lon2:.6f}°\nFinal Azimuth: {az2:.6f}°")

```

### Radius of Curvature for a given Azimuth using Euler's equation.
```python
from pygeodetics import Geodetic

lat = 45
azimuth = 30

radius = Geodetic().radius_of_curvature(lat, azimuth, radians=False)
print(f"Radius of Curvature:\n{radius:.3f} meters")

```

### Calculate the mean Radius of the International1924 Ellipsoid
```python
from pygeodetics import Geodetic
from pygeodetics.Ellipsoid import International1924

geod = Geodetic(International1924())

mean_radius = geod.get_mean_radius()
print(f"Mean Radius of the Ellipsoid:\n{mean_radius:.3f} meters")

```

### Calculate the distance between two points on the ellipsoid (Vincenty formula)
```python
from pygeodetics import Geodetic

# Define the coordinates of the first point
lat1 = 52.2296756
lon1 = 21.0122287

# Define the coordinates of the second point
lat2 = 41.8919300
lon2 = 12.5113300

distances = Geodetic().distance_between_two_points(lon1, lat1, lon2, lat2, radians=False)
print(f"Distances between the two points: {distances}")

```

### Calculate the meridional radius of curvature (M) at a given latitude

```python
from pygeodetics import Geodetic

# Compute the mean radius of the ellipsoid at a given latitude
lat = 61.456121547 # Latitude in degrees
mradius = Geodetic().mrad(lat)
print(f"Mean Radius of the Ellipsoid at Latitude {lat}°: {mradius:.3f} meters")
```


### Calculate the normal radius of curvature (N) at a given latitude.

```python
from pygeodetics import Geodetic

# Compute the normal radius of the ellipsoid at a given latitude
lat = 61.456121547 # Latitude in degrees
mradius = Geodetic().nrad(lat)
print(f"Normal Radius of the Ellipsoid at Latitude {lat}°:\n{mradius:.3f} meters")
```



## License
This project is licensed under the MIT License.

