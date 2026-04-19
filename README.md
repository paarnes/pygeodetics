<p align="center">
    <img src="https://github.com/paarnes/pygeodetics/blob/main/docs/icon/pygeodetics.png?raw=true" alt="PyGeodetics Logo" width="700">
</p>

<p align="center">
    <a href="https://github.com/paarnes/pygeodetics/actions/workflows/tests.yml"><img alt="Tests" src="https://github.com/paarnes/pygeodetics/actions/workflows/tests.yml/badge.svg"></a>
    <a href="https://pypi.org/project/pygeodetics/"><img alt="PyPI" src="https://img.shields.io/pypi/v/pygeodetics?color=blue"></a>
    <a href="https://pypi.org/project/pygeodetics/"><img alt="Python versions" src="https://img.shields.io/pypi/pyversions/pygeodetics"></a>
    <a href="https://github.com/paarnes/pygeodetics/blob/main/LICENSE"><img alt="License: MIT" src="https://img.shields.io/github/license/paarnes/pygeodetics"></a>
    <a href="https://github.com/paarnes/pygeodetics/stargazers"><img alt="GitHub stars" src="https://img.shields.io/github/stars/paarnes/pygeodetics"></a>
    <a href="https://github.com/paarnes/pygeodetics/issues"><img alt="GitHub issues" src="https://img.shields.io/github/issues/paarnes/pygeodetics"></a>
    <a href="https://pypi.org/project/pygeodetics/"><img alt="Downloads" src="https://img.shields.io/pypi/dm/pygeodetics"></a>
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
- [Use of Mercator Variant C projection](#use-of-mercator-variant-c-projection)
- [Use of Transverse Mercator projection](#use-of-transverse-mercator-projection)

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
print(f"Geodetic Inverse Problem:\nForward Azimuth: {az1:.6f}°\nFinal Azimuth at P2: {az2:.6f}°\nDistance: {distance:.3f} m")

```

### Geodetic Direct Problem on the GRS80 ellipsoid
```python
from pygeodetics import Geodetic
from pygeodetics.Ellipsoid import GRS80

geod = Geodetic(GRS80())

lat1, lon1 = 52.2296756, 21.0122287
az1 = -147.4628043168
d = 1316208.08334

lat2, lon2, az2 = geod.direct_problem(lat1, lon1, az1, d, quadrant_correction=True)
print(f"Geodetic Direct Problem:\nDestination Latitude: {lat2:.6f}°\nDestination Longitude: {lon2:.6f}°\nFinal Azimuth at destination: {az2:.6f}°")

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

### Use of Mercator Variant C projection
```python

from pygeodetics import MercatorVariantC
import numpy as np

conv = MercatorVariantC(
    a=6378245.0, f=1/298.3,
    latSP1=np.deg2rad(42), latFO=np.deg2rad(42), lonFO=np.deg2rad(51),
    EFO=0.0, NFO=0.0
)

lon, lat = 53, 53
E, N = conv.geog_to_projected([[lon, lat]], unit="deg").ravel()
rlon, rlat = conv.projected_to_geog([[E, N]]).ravel()

print(f"Easting = {E:.2f} m\nNorthing = {N:.2f} m")
print(f"Reversed lon = {rlon:.8f}°\nReversed lat = {rlat:.8f}°")
```


### Use of Transverse Mercator projection

```python
import numpy as np
from pygeodetics import TransverseMercator
np.set_printoptions(precision=8, suppress=True)

# Projection parameters (radians for origins)
lat_origin = 0.0                         # φ0 (rad)
lon_origin = np.radians(9.0)             # λ0 (rad), e.g. UTM zone 32
scale_factor = 0.9996                    # k0
false_easting = 500000.0                 # FE (m)
false_northing = 0.0                     # FN (m)

# Ellipsoid (WGS84)
a = 6378137.0
f = 1 / 298.257223563

tm = TransverseMercator(
    lat_origin, lon_origin, scale_factor,
    false_easting, false_northing, a, f
)

# Input coordinates (lon, lat, h) 
coords = np.array([
    [3.0, 60.0, 100.0],
    [3.2, 61.0, 102.0],
])

# Perform forward projection
easting, northing, height = tm.geog_to_projected(coords, unit="deg")
results = np.vstack((easting, northing, height)).T
print(f"\nProjected Coordinates TM:\n{results}")

# Perform inverse projection
proj_coordinates = np.vstack([easting, northing, height]).T
lon_back, lat_back, height = tm.projected_to_geog(proj_coordinates, unit="deg")
results = np.vstack((lon_back, lat_back, height)).T
print(f"\nGeographic Coordinates TM:\n{results}")

```


## Math and the Theory Basis

This section provides the mathematical foundation for the computations performed in PyGeodetics. 

### 1. Geodetic to ECEF Conversion

The conversion from geodetic coordinates $(\phi, \lambda, h)$ to Earth-Centered, Earth-Fixed (ECEF) coordinates $(X, Y, Z)$ is given by:

$$
\begin{aligned}
    X &= (N + h) \cos\phi \cos\lambda, \\
    Y &= (N + h) \cos\phi \sin\lambda, \\
    Z &= \left( \frac{b^2}{a^2} N + h \right) \sin\phi,
\end{aligned}
$$

where:
- $N = \frac{a}{\sqrt{1 - e^2 \sin^2\phi}}$ is the prime vertical radius of curvature,
- $a$ is the semi-major axis,
- $b$ is the semi-minor axis,
- $e^2 = \frac{a^2 - b^2}{a^2}$ is the first eccentricity squared,
- $\phi$ is the geodetic latitude,
- $\lambda$ is the geodetic longitude,
- $h$ is the height above the ellipsoid.

### 2. ECEF to Geodetic Conversion

PyGeodetics implements three solvers for the inverse problem $(X, Y, Z) \to (\phi, \lambda, h)$ — an iterative Heiskanen–Moritz solver, Bowring's closed-form method, and Vermeille's closed-form method. The recommended (default) solver is Bowring's method, summarised below:

1. Compute the longitude:

   $$\lambda=\arctan2(Y, X)$$

2. Compute the intermediate values:

   $$p=\sqrt{X^2 + Y^2}, \quad \theta = \arctan2(Z\,a,\; p\,b)$$

3. Compute the latitude:

   $$\phi=\arctan2\!\left(Z + e'^2 b \sin^3\theta,\; p - e^2 a \cos^3\theta\right)$$
   ,where $e'^2=\frac{a^2 - b^2}{b^2}$ is the second eccentricity squared.

4. Compute the height:

   $$h=\frac{p}{\cos\phi}-N$$

### 3. Geodetic Inverse Problem

The geodetic inverse problem calculates the geodesic distance $s$ and azimuths $(\alpha_1, \alpha_2)$ between two points $(\phi_1, \lambda_1)$ and $(\phi_2, \lambda_2)$. PyGeodetics solves it with Vincenty's formulae, which iterate on the auxiliary sphere using the *reduced* latitudes:

1. Compute the reduced latitudes (also called parametric latitudes):

   $$U_1 = \arctan\left((1 - f) \tan\phi_1\right), \quad U_2 = \arctan\left((1 - f) \tan\phi_2\right)$$
   ,where $f = \frac{a - b}{a}$ is the flattening.

2. Set $L = \lambda_2 - \lambda_1$ and initialise $\lambda = L$. Iterate $\lambda$ (longitude on the auxiliary sphere) together with the angular distance $\sigma$:

   $$\sigma = \arctan2\!\left(\sqrt{(\cos U_2 \sin\lambda)^2 + (\cos U_1 \sin U_2 - \sin U_1 \cos U_2 \cos\lambda)^2},\; \sin U_1 \sin U_2 + \cos U_1 \cos U_2 \cos\lambda\right)$$

   The full update equation for $\lambda$ and the auxiliary quantities $\sin\alpha$, $\cos 2\sigma_m$, $C$ are given in [Section 7](#7-vincentys-formula-for-geodesic-distance).

3. Once converged, compute the forward and reverse azimuths:

   $$\alpha_1 = \arctan2\left(\cos U_2 \sin\lambda,\; \cos U_1 \sin U_2 - \sin U_1 \cos U_2 \cos\lambda\right)$$

   $$\alpha_2 = \arctan2\left(\cos U_1 \sin\lambda,\; -\sin U_1 \cos U_2 + \cos U_1 \sin U_2 \cos\lambda\right)$$

4. Compute the geodesic distance using the series expansion in $u^2 = \cos^2\alpha\,(a^2 - b^2)/b^2$:

   $$s = b\,A\,(\sigma - \Delta\sigma)$$

   with $A$, $B$ and $\Delta\sigma$ defined in [Section 7](#7-vincentys-formula-for-geodesic-distance). The simpler one-term form $s \approx b(\sigma - \tfrac{f}{4} \sin\sigma \cos(2\sigma_m + \sigma))$ is **not** used here — it loses several kilometres of accuracy on continental baselines.

### 4. Radius of Curvature

The radius of curvature in the meridian ($M$) and the prime vertical ($N$) are given by:

$$
M = \frac{a(1 - e^2)}{(1 - e^2 \sin^2\phi)^{3/2}}, \quad N = \frac{a}{\sqrt{1 - e^2 \sin^2\phi}}.
$$

### 5. Local ENU Coordinates

The transformation from ECEF to local East-North-Up (ENU) coordinates is given by:

$$
\begin{bmatrix}
e \\
n \\
u
\end{bmatrix}=\begin{bmatrix}
-\sin\lambda & \cos\lambda & 0 \\
-\sin\phi\cos\lambda & -\sin\phi\sin\lambda & \cos\phi \\
\cos\phi\cos\lambda & \cos\phi\sin\lambda & \sin\phi
\end{bmatrix}
\begin{bmatrix}
X - X_0 \\
Y - Y_0 \\
Z - Z_0
\end{bmatrix}
$$



where $(X_0, Y_0, Z_0)$ are the ECEF coordinates of the reference point.

### 6. Mean Radius of the Ellipsoid

PyGeodetics returns the IUGG arithmetic mean radius:

$$
R_1 = \frac{2a + b}{3}.
$$

This is one of several conventional definitions (others include the authalic and volumetric mean radii); it is the simple arithmetic mean of the three semi-axes $(a, a, b)$.

### 7. Vincenty's Formula for Geodesic Distance

Vincenty's formula is used to calculate the geodesic distance between two points on an ellipsoid. The method iteratively solves for the distance $s$ and azimuths $(\alpha_1, \alpha_2)$ between two points $(\phi_1, \lambda_1)$ and $(\phi_2, \lambda_2)$.

#### Steps:

1. **Reduced Latitude**:
   Compute the reduced latitudes:

   $$U_1=\arctan\left((1 - f) \cdot \tan\phi_1\right), \quad U_2=\arctan\left((1 - f) \cdot\tan\phi_2\right)$$

   ,where $f = \frac{a - b}{a}$ is the flattening.

2. **Longitude Difference**:
   Compute the difference in longitudes:

   $$L=\lambda_2 - \lambda_1$$

3. **Iterative Solution**:
   Initialize $\lambda = L$ and iteratively solve for $\lambda$ using:
   
   $$\lambda=L + (1 - C) f \sin\alpha \left[\sigma + C \sin\sigma \left(\cos2\sigma_m + C \cos\sigma \left(-1 + 2 \cos^2 2\sigma_m\right)\right)\right]$$

   where:
   - $\sigma = \arctan2\left(\sqrt{(\cos U_2 \sin\lambda)^2 + (\cos U_1 \sin U_2 - \sin U_1 \cos U_2 \cos\lambda)^2}, \sin U_1 \sin U_2 + \cos U_1 \cos U_2 \cos\lambda\right)$,
   - $\sin\alpha = \frac{\cos U_1 \cos U_2 \sin\lambda}{\sin\sigma}$,
   - $\cos2\sigma_m = \cos\sigma - \frac{2 \sin U_1 \sin U_2}{\cos^2\alpha}$,
   - $C = \frac{f}{16} \cos^2\alpha \left(4 + f \left(4 - 3 \cos^2\alpha\right)\right)$.

   The iteration stops when $|\lambda - \lambda_{\text{prev}}| < \text{tolerance}$.

4. **Geodesic Distance**:
   Compute the geodesic distance:

   $$s=b A \left[\sigma - \delta\sigma\right]$$
   ,where:
   - $u^2=\frac{\cos^2\alpha (a^2 - b^2)}{b^2}$,
   - $A=1 + \frac{u^2}{16384} \left(4096 + u^2 \left(-768 + u^2 (320 - 175 u^2)\right)\right)$,
   - $B=\frac{u^2}{1024} \left(256 + u^2 \left(-128 + u^2 (74 - 47 u^2)\right)\right)$,
   - $\delta\sigma=B \sin\sigma \left[\cos2\sigma_m + \frac{B}{4} \left(\cos\sigma \left(-1 + 2 \cos^2 2\sigma_m\right) - \frac{B}{6} \cos2\sigma_m \left(-3 + 4 \sin^2\sigma\right) \left(-3 + 4 \cos^2 2\sigma_m\right)\right)\right]$.

5. **Azimuths**:
   Compute the forward and reverse azimuths:

   $$\alpha_1 = \arctan2\left(\cos U_2 \sin\lambda, \cos U_1 \sin U_2 - \sin U_1 \cos U_2 \cos\lambda\right)$$

   $$\alpha_2 = \arctan2\left(\cos U_1 \sin\lambda, -\sin U_1 \cos U_2 + \cos U_1 \sin U_2 \cos\lambda\right)$$

Vincenty's formula is highly accurate for most geodetic calculations but may fail to converge for nearly antipodal points.

### 8. Meridional Radius of Curvature (M)

The meridional radius of curvature, denoted as $M$, represents the radius of curvature in the north-south direction along a meridian. It is computed as:

$$M=\frac{a(1 - e^2)}{(1 - e^2 \sin^2\phi)^{3/2}}$$

where:
- $a$ is the semi-major axis of the ellipsoid,
- $e^2 = \frac{a^2 - b^2}{a^2}$ is the first eccentricity squared,
- $\phi$ is the geodetic latitude.

This formula accounts for the flattening of the ellipsoid and the variation in curvature with latitude.

### 9. Normal Radius of Curvature (N)

The normal radius of curvature, denoted as $N$, represents the radius of curvature in the east-west direction perpendicular to the meridian. It is computed as:

$$
N = \frac{a}{\sqrt{1 - e^2 \sin^2\phi}},
$$

where:
- $a$ is the semi-major axis of the ellipsoid,
- $e^2 = \frac{a^2 - b^2}{a^2}$ is the first eccentricity squared,
- $\phi$ is the geodetic latitude.

The value of $N$ varies with latitude due to the ellipsoidal shape of the Earth, being largest at the equator and smallest at the poles.

### 10. Grid Convergence in the Transverse Mercator Projection

The grid convergence, denoted as $\gamma$, is the angular difference between grid north and true north in the Transverse Mercator (TM) projection. It can be computed using either geographic coordinates or projected coordinates.

#### 10.1 Grid Convergence Using Geographic Coordinates

The grid convergence $\gamma$ at a point $(\phi, \lambda)$ is given by:

$$
\gamma = \Delta\lambda \sin\phi + \frac{\Delta\lambda^3}{3} \sin\phi \cos^2\phi (1 + 3\epsilon^2 + 2\epsilon^4) + \frac{\Delta\lambda^5}{15} \sin\phi \cos^4\phi (2 - \tan^2\phi),
$$

where:
- $\Delta\lambda = \lambda - \lambda_0$ is the longitude difference from the central meridian,
- $\epsilon^2 = \frac{e^2}{1 - e^2} \cos^2\phi$ is the second eccentricity squared,
- $e^2 = \frac{a^2 - b^2}{a^2}$ is the first eccentricity squared,
- $a$ is the semi-major axis,
- $b$ is the semi-minor axis,
- $\phi$ is the geodetic latitude,
- $\lambda$ is the geodetic longitude,
- $\lambda_0$ is the central meridian.

#### 10.2 Grid Convergence Using Projected Coordinates

The grid convergence $\gamma$ at a point $(x, y)$ in projected coordinates is given by:

$$
\gamma = \frac{x \tan\phi_f}{N_f} - \frac{x^3 \tan\phi_f}{3 N_f^3} \left(1 + \tan^2\phi_f - \epsilon_f^2 - 2\epsilon_f^4\right),
$$

where:
- $\phi_f$ is the footpoint latitude, computed iteratively,
- $N_f = \frac{a}{\sqrt{1 - e^2 \sin^2\phi_f}}$ is the normal radius of curvature at the footpoint latitude,
- $\epsilon_f^2 = \frac{e^2}{1 - e^2} \cos^2\phi_f$ is the second eccentricity squared at the footpoint latitude,
- $x$ is the easting coordinate (adjusted for false easting),
- $y$ is the northing coordinate.

### 11. Scale Factor in the Transverse Mercator Projection

The scale factor, denoted as $k$, describes the distortion of distances in the Transverse Mercator projection. It can be computed using either geographic coordinates or projected coordinates.

#### 11.1 Scale Factor Using Geographic Coordinates

The scale factor $k$ at a point $(\phi, \lambda)$ is given by:

$$
k = 1 + \frac{\Delta\lambda^2}{2} \cos^2\phi (1 + \epsilon^2) + \frac{\Delta\lambda^4}{24} \cos^4\phi (5 + 4\tan^2\phi),
$$

where:
- $\Delta\lambda = \lambda - \lambda_0$ is the longitude difference from the central meridian,
- $\epsilon^2 = \frac{e^2}{1 - e^2} \cos^2\phi$ is the second eccentricity squared,
- $e^2 = \frac{a^2 - b^2}{a^2}$ is the first eccentricity squared,
- $a$ is the semi-major axis,
- $b$ is the semi-minor axis,
- $\phi$ is the geodetic latitude,
- $\lambda$ is the geodetic longitude,
- $\lambda_0$ is the central meridian.

#### 11.2 Scale Factor Using Projected Coordinates

The scale factor $k$ at a point $(x, y)$ in projected coordinates is given by:

$$
k = 1 + \frac{x^2}{2 M_f N_f} + \frac{x^4}{24 N_f^4},
$$

where:
- $M_f = \frac{a(1 - e^2)}{(1 - e^2 \sin^2\phi_f)^{3/2}}$ is the meridional radius of curvature at the footpoint latitude,
- $N_f = \frac{a}{\sqrt{1 - e^2 \sin^2\phi_f}}$ is the normal radius of curvature at the footpoint latitude,
- $\phi_f$ is the footpoint latitude, computed iteratively,
- $x$ is the easting coordinate (adjusted for false easting),
- $y$ is the northing coordinate.

#### 11.3 Scale Factor for a Sphere

For a spherical Earth model, the scale factor $k$ is given by:

$$
k = \cosh\left(\frac{x}{R}\right),
$$

where:
- $x$ is the easting coordinate (adjusted for false easting),
- $R$ is the radius of the sphere.

### 12. Projections

This section explains the mathematical foundation for the projections method supported by the library. In this release supported methods are:
 - Mercator Variant C 
 - Transverse Mercator projections.

#### 12.1 Mercator Variant C Projection

The Mercator Variant C projection is a cylindrical map projection that preserves angles, making it conformal.
This projection is a variant of the Mercator projection, with specific parameters defined for its implementation.

##### Parameters of the Mercator Variant C Projection

The following table summarizes the key parameters of the Mercator Variant C projection as defined by EPSG:

| **Parameter Name**                  | **Parameter EPSG Code** | **Sign Reversal** | **Description**                                                                                                                                                                                                 |
|-------------------------------------|--------------------|-------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Latitude of 1st standard parallel   | 8823               | No                | For a conic projection with two standard parallels, this is the latitude of one of the parallels of intersection of the cone with the ellipsoid. Scale is true along this parallel.                              |
| Longitude of natural origin         | 8802               | No                | The longitude of the point from which the values of both the geographical coordinates on the ellipsoid and the grid coordinates on the projection are deemed to increment or decrement for computational purposes. |
| Latitude of false origin            | 8821               | No                | The latitude of the point which is not the natural origin and at which grid coordinate values false easting and false northing are defined.                                                                     |
| Easting at false origin             | 8826               | No                | The easting value assigned to the false origin.                                                                                                                                                                 |
| Northing at false origin            | 8827               | No                | The northing value assigned to the false origin.                                                                                                                   

### 12.1.1 Forward Projection (Geographic → Projected)

The forward projection transforms geographic coordinates $(\lambda, \phi)$ (longitude, latitude) to projected coordinates $(E, N)$ (easting, northing).  

- Constants:
  - $k_0 = \frac{\cos \phi_1}{\sqrt{1 - e^2 \sin^2 \phi_1}}$
  - All logarithms below are **natural** logarithms (base $\mathrm{e}$). To avoid clashing with the eccentricity $e$, the natural-log base is denoted $\mathrm{e}$ in this section.

1. **Compute the projection constant $k_0$:**

   $$k_0 = \frac{\cos \phi_1}{\sqrt{1 - e^2 \sin^2 \phi_1}}$$

   where:
   - $\phi_1$ = latitude of the first standard parallel (absolute value, positive),
   - $e^2 = \tfrac{a^2 - b^2}{a^2}$ = eccentricity squared,
   - $a$ = semi-major axis, $b$ = semi-minor axis.

2. **Compute the meridional arc at latitude $\phi$:**

   $$M(\phi) = a k_0 \,\ln \left[
   \tan\!\left(\tfrac{\pi}{4} + \tfrac{\phi}{2}\right) 
   \left(\frac{1 - e \sin\phi}{1 + e \sin\phi}\right)^{\tfrac{e}{2}}
   \right]$$

3. **Compute Easting and Northing:**

   $$E = E_0 + a k_0 (\lambda - \lambda_0)$$

   $$N = N_0 - M(\phi_0) + M(\phi)$$

   where:
   - $(\lambda_0, \phi_0)$ = longitude/latitude of the false origin,
   - $(E_0, N_0)$ = false easting and northing.

---

### 12.1.2 Inverse Projection (Projected → Geographic)

The inverse projection transforms projected coordinates $(E, N)$ back to geographic $(\lambda, \phi)$.

1. **Longitude:**

   $$\lambda = \lambda_0 + \frac{E - E_0}{a k_0}$$

2. **Latitude (iterative series):**

   Define:
   $$\chi = \frac{\pi}{2} - 2 \arctan \!\left( 
   \exp\!\left[-\frac{N - N_0 + M(\phi_0)}{a k_0}\right] 
   \right)$$

   Then:

   $$\phi=\chi+\left(\tfrac{e^2}{2} + \tfrac{5 e^4}{24} + \tfrac{e^6}{12} + \tfrac{13 e^8}{360}\right) \sin(2\chi)+\left(\tfrac{7 e^4}{48} +\tfrac{29 e^6}{240} + \tfrac{811 e^8}{11520}\right) \sin(4\chi)+\left(\tfrac{7 e^6}{120} + \tfrac{81 e^8}{1120}\right) \sin(6\chi)+\left(\tfrac{4279 e^8}{161280}\right) \sin(8\chi)$$

---


#### 12.2 Transverse Mercator Projection

The Transverse Mercator (TM) projection is a conformal map projection widely used for large-scale mapping, such as the Universal Transverse Mercator (UTM) system.

#### Parameters of the Transverse Mercator Projection

The following table summarizes the key parameters of the Transverse Mercator projection as defined by EPSG:

| **Parameter Name**               | **Parameter EPSG Code** | **Sign Reversal** | **Description**                                                                                                                                                                                                 |
|----------------------------------|--------------------|-------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Latitude of natural origin       | 8801               | No                | The latitude of the point from which the values of both the geographical coordinates on the ellipsoid and the grid coordinates on the projection are deemed to increment or decrement for computational purposes. |
| Longitude of natural origin      | 8802               | No                | The longitude of the point from which the values of both the geographical coordinates on the ellipsoid and the grid coordinates on the projection are deemed to increment or decrement for computational purposes. |
| Scale factor at natural origin   | 8805               | No                | The factor by which the map grid is reduced or enlarged during the projection process, defined by its value at the natural origin.                                                                              |
| False easting                    | 8806               | No                | The value assigned to the abscissa (east or west) axis of the projection grid at the natural origin to avoid negative coordinates over parts of the mapped area.                                                |
| False northing                   | 8807               | No                | The value assigned to the ordinate (north or south) axis of the projection grid at the natural origin to avoid negative coordinates over parts of the mapped area.                                              |


##### 12.2.1 Forward Projection (Geographic to Projected)

The forward projection transforms geographic coordinates $(\lambda, \phi)$ to projected coordinates $(E, N)$ as follows:

1. **Calculate constants for the projection**:

   $$n = \frac{f}{2 - f}$$

   $$B = \frac{a}{1 + n} \left(1 + \frac{n^2}{4} + \frac{n^4}{64}\right)$$

   $$h_1=\frac{n}{2} - \frac{2}{3}n^2 + \frac{5}{16}n^3 + \frac{41}{180}n^4$$

   $$h_2=\frac{13}{48}n^2 - \frac{3}{5}n^3 + \frac{557}{1440}n^4$$

   $$h_3=\frac{61}{240}n^3 - \frac{103}{140}n^4$$

   $$h_4=\frac{49561}{161280}n^4$$

2. **Compute the meridional arc distance from the equator to the projection origin**:

   If $\phi_0 = 0$

      $$M_0 = 0$$

   or if $\phi_0 = 90^\circ , (\pi/2 , \text{radians})$
   
      $$M_0 = B \cdot \frac{\pi}{2}$$

   or if $\phi_0 = -90^\circ , (-\pi/2 , \text{radians})$ 
   
      $$M_0 = B \cdot \left(-\frac{\pi}{2}\right)$$
   
   Otherwise:

   $$Q_0 = \sinh^{-1}(\tan\phi_0) - e \tanh^{-1}(e \sin\phi_0)$$

   $$\beta_0 = \tan^{-1}(\sinh Q_0)$$

   $$\xi_{0,0} = \sin^{-1}(\sin\beta_0)$$

   $$\xi_0 = \xi_{0,0} + h_1 \sin(2\xi_{0,0}) + h_2 \sin(4\xi_{0,0}) + h_3 \sin(6\xi_{0,0}) + h_4 \sin(8\xi_{0,0})$$

   $$M_0 = B \xi_0$$

3. **Compute intermediate values for the given latitude $\phi$**:

   $$Q = \sinh^{-1}(\tan\phi) - e \tanh^{-1}(e \sin\phi)$$

   $$\beta = \tan^{-1}(\sinh Q)$$

   $$\eta_0 = \tanh^{-1}(\cos\beta \sin(\lambda - \lambda_0))$$

   $$\xi_0 = \sin^{-1}(\sin\beta \cosh\eta_0)$$

   $$\xi = \xi_0 + h_1 \sin(2\xi_0) \cosh(2\eta_0) + h_2 \sin(4\xi_0) \cosh(4\eta_0) + h_3 \sin(6\xi_0) \cosh(6\eta_0) + h_4 \sin(8\xi_0) \cosh(8\eta_0)$$

   $$\eta = \eta_0 + h_1 \cos(2\xi_0) \sinh(2\eta_0) + h_2 \cos(4\xi_0) \sinh(4\eta_0) + h_3 \cos(6\xi_0) \sinh(6\eta_0) + h_4 \cos(8\xi_0) \sinh(8\eta_0)$$

4. **Compute the easting and northing**:

   $$E = E_0 + k_0 B \eta$$

   $$N = N_0 + k_0 \left(B \xi - M_0\right)$$

##### 12.2.2 Inverse Projection (Projected to Geographic)

The inverse projection transforms projected coordinates $(E, N)$ back to geographic coordinates $(\lambda, \phi)$ as follows:

1. **Compute intermediate values**:

   $$\eta' = \frac{E - E_0}{k_0 B}, \quad \xi' = \frac{(N - N_0) + k_0 M_0}{k_0 B}$$

2. **Iteratively compute $\xi_0'$ and $\eta_0'$**:

   $$\xi_0' = \xi' - \left(h_1 \sin(2\xi') \cosh(2\eta') + h_2 \sin(4\xi') \cosh(4\eta') + h_3 \sin(6\xi') \cosh(6\eta') + h_4 \sin(8\xi') \cosh(8\eta')\right)$$

   $$\eta_0' = \eta' - \left(h_1 \cos(2\xi') \sinh(2\eta') + h_2 \cos(4\xi') \sinh(4\eta') + h_3 \cos(6\xi') \sinh(6\eta') + h_4 \cos(8\xi') \sinh(8\eta')\right)$$

3. **Compute $\beta'$ and $Q'$**:

   $$\beta' = \sin^{-1}\left(\frac{\sin\xi_0'}{\cosh\eta_0'}\right)$$

   $$Q' = \sinh^{-1}(\tan\beta')$$

4. **Iteratively compute latitude $\phi$**:

   $$Q'' = Q' + e \tanh^{-1}(e \tanh Q')$$

   Repeat until the change in $Q''$ is insignificant. Then:

   $$\phi = \tan^{-1}(\sinh Q'')$$

5. **Compute longitude $\lambda$**:

   $$\lambda = \lambda_0 + \sin^{-1}\left(\frac{\tanh\eta_0'}{\cos\beta'}\right)$$


### EPSG Reference

- These formulas follow **EPSG Guidance Note 7-2**.  

## License
This project is licensed under the MIT License.

