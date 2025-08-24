<p align="center">
    <img src="https://github.com/paarnes/pygeodetics/blob/main/docs/icon/pygeodetics.png?raw=true" alt="PyGeodetics Logo" width="700">
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

The iterative method for converting ECEF coordinates $(X, Y, Z)$ to geodetic coordinates $(\phi, \lambda, h)$ involves:

1. Compute the longitude:

   $$\lambda=\arctan2(Y, X)$$

2. Compute the intermediate values:

   $$p=\sqrt{X^2 + Y^2}, \quad \theta = \arctan\left(\frac{Z a}{p b}\right)$$

3. Compute the latitude:

   $$\phi=\arctan\left(\frac{Z + e'^2 b \sin^3\theta}{p - e^2 a \cos^3\theta}\right)$$
   ,where $e'^2=\frac{a^2 - b^2}{b^2}$ is the second eccentricity squared.

4. Compute the height:

   $$h=\frac{p}{\cos\phi}-N$$

### 3. Geodetic Inverse Problem

The geodetic inverse problem calculates the geodesic distance $s$ and azimuths $(\alpha_1, \alpha_2)$ between two points $(\phi_1, \lambda_1)$ and $(\phi_2, \lambda_2)$. Using Vincenty's formulae:

1. Compute the reduced latitude:

   $$\beta=\arctan\left((1 - f) \tan\phi\right)$$
   ,where $f = \frac{a - b}{a}$ is the flattening.

2. Iteratively solve for the longitude difference $\Delta\lambda$ and the spherical arc $\sigma$:

   $$\sigma = \arctan\left(\frac{\sqrt{(\cos\beta_2 \sin\Delta\lambda)^2 + (\cos\beta_1 \sin\beta_2 - \sin\beta_1 \cos\beta_2 \cos\Delta\lambda)^2}}{\sin\beta_1 \sin\beta_2 + \cos\beta_1 \cos\beta_2 \cos\Delta\lambda}\right)$$

3. Compute the azimuths:

   $$\alpha_1 = \arctan2\left(\cos\beta_2 \sin\Delta\lambda, \cos\beta_1 \sin\beta_2 - \sin\beta_1 \cos\beta_2 \cos\Delta\lambda\right)$$

   $$\alpha_2 = \arctan2\left(\cos\beta_1 \sin\Delta\lambda, -\sin\beta_1 \cos\beta_2 + \cos\beta_1 \sin\beta_2 \cos\Delta\lambda\right)$$

4. Compute the geodesic distance:

   $$s=b \left(\sigma - \frac{f}{4} \sin\sigma \cos(2\sigma_m + \sigma)\right)$$
   ,where $\sigma_m = \frac{\sigma_1 + \sigma_2}{2}$.

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

The mean radius of the ellipsoid is computed as:

$$
R = \frac{2a + b}{3}.
$$

This formula provides an average radius considering both the equatorial and polar radii.

### 7. Vincenty's Formula for Geodesic Distance

Vincenty's formula is used to calculate the geodesic distance between two points on an ellipsoid. The method iteratively solves for the distance $s$ and azimuths $(\alpha_1, \alpha_2)$ between two points $(\phi_1, \lambda_1)$ and $(\phi_2, \lambda_2)$.

#### Steps:

1. **Reduced Latitude**:
   Compute the reduced latitudes:

   $$U_1=\arctan\left((1 - f) \tan\phi_1\right), \quad U_2=\arctan\left((1 - f) \tan\phi_2\right)$$

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

##### 12.1.1 Forward Projection (Geographic to Projected)

The forward projection transforms geographic coordinates $(\lambda, \phi)$ (longitude, latitude) to projected coordinates $(E, N)$ (easting, northing) as follows:

1. Compute the projection constant $k_0$:

   $$k_0=\frac{\cos\phi_1}{\sqrt{1 - e^2 \sin^2\phi_1}}$$
   
   where:
   - $\phi_1$ is the latitude of the first standard parallel,
   - $e^2 = \frac{a^2 - b^2}{a^2}$ is the first eccentricity squared.

2. Compute the meridional arc $M$:

   $$M(\phi)=a k_0 \ln\left(\tan\left(\frac{\pi}{4} + \frac{\phi}{2}\right) \cdot \left(\frac{1 - e \sin\phi}{1 + e \sin\phi}\right)^{\frac{e}{2}}\right),
   $$

   where:
   - $a$ is the semi-major axis,
   - $e$ is the eccentricity.

3. Compute the easting and northing:

   $$E=E_0 + a k_0 (\lambda - \lambda_0)$$

   $$N = N_0 + M(\phi) - M(\phi_0)$$

   where:
   - $\lambda_0$ is the longitude of the false origin,
   - $\phi_0$ is the latitude of the false origin,
   - $E_0$ and $N_0$ are the false easting and northing.

##### 12.1.2 Inverse Projection (Projected to Geographic)

The inverse projection transforms projected coordinates $(E, N)$ back to geographic coordinates $(\lambda, \phi)$:

1. Compute the longitude:

   $$\lambda = \lambda_0 + \frac{E - E_0}{a k_0}$$

2. Compute the latitude iteratively:
   
   $$\phi = \chi + \sum_{n=1}^4 c_n \sin(2n\chi)$$

   where:
   - $\chi = \frac{\pi}{2} - 2 \arctan\left(e^{-\frac{N - N_0 + M(\phi_0)}{a k_0}}\right)$,
   - $c_n$ are coefficients derived from the eccentricity.

#### 12.2 Transverse Mercator Projection

The Transverse Mercator (TM) projection is a conformal map projection widely used for large-scale mapping, such as the Universal Transverse Mercator (UTM) system.

##### 12.2.1 Forward Projection (Geographic to Projected)

The forward projection transforms geographic coordinates $(\lambda, \phi)$ to projected coordinates $(E, N)$:

1. Compute the meridional arc $M$:

   $$M=a \left(A_0 \phi - A_2 \sin(2\phi) + A_4 \sin(4\phi) - A_6 \sin(6\phi)\right)$$
   
   where:
   - $A_0$, $A_2$, $A_4$, $A_6$ are coefficients derived from the eccentricity.

2. Compute the easting and northing:
   
   $$E=E_0 + k_0 N \eta$$
   
   $$N=N_0 + k_0 \left(M - M_0 + \frac{\nu}{2} \sin\phi \cos\phi \eta^2 + \frac{\nu}{24} \sin\phi \cos^3\phi (5 - \tan^2\phi + 9\epsilon^2) \eta^4\right)$$

   where:
   - $\eta = \lambda - \lambda_0$,
   - $\nu = \frac{a}{\sqrt{1 - e^2 \sin^2\phi}}$ is the radius of curvature in the prime vertical,
   - $\epsilon^2 = \frac{e^2}{1 - e^2} \cos^2\phi$ is the second eccentricity squared.

##### 12.2.2 Inverse Projection (Projected to Geographic)

The inverse projection transforms projected coordinates $(E, N)$ back to geographic coordinates $(\lambda, \phi)$:

1. Compute the footpoint latitude $\phi_f$:

   $$\phi_f = \frac{N - N_0}{a A_0}$$

2. Compute the latitude iteratively:
   
   $$\phi = \phi_f + \sum_{n=1}^4 b_n \sin(2n\phi_f)$$
   
   where:
   - $b_n$ are coefficients derived from the eccentricity.

3. Compute the longitude:
   
   $$\lambda = \lambda_0 + \frac{E - E_0}{k_0 \nu}$$



## License
This project is licensed under the MIT License.

