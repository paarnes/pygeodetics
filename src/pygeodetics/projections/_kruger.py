"""
Kruger series for the ellipsoidal Transverse Mercator projection.

This is the expansion in the third flattening ``n = f / (2 - f)`` described in

    C. F. F. Karney (2011), *Transverse Mercator with an accuracy of a few nanometers*,
    Journal of Geodesy 85(8), 475-485,

which is also the algorithm used by PROJ and GeographicLib. The series is carried to
``n**6``. Because ``n`` is about ``1/599`` for terrestrial ellipsoids, the neglected
``n**7`` term is of order ``1e-19`` - below double precision - so the result is exact
for every practical purpose, and it stays exact far outside a single 6-degree zone.

Contrast this with the classic Snyder expansions in ``dlon`` or ``x``, which are
truncated in a quantity that reaches ~0.14 rad at a UTM zone edge and therefore lose
accuracy rapidly with distance from the central meridian.

The module is private; it exists so that the grid convergence, point scale factor and
the projection itself share one verified implementation.

author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com
"""

from functools import lru_cache
from typing import Tuple, Union
import numpy as np

# Series order. Karney's coefficients are given to n**6.
_ORDER = 6
_J = np.arange(1, _ORDER + 1)


def _series(coeffs: np.ndarray, zeta):
    """Sum ``coeffs[j] * sin(2 (j+1) zeta)`` over the series, for scalar or array `zeta`."""
    z = np.asarray(zeta)[..., np.newaxis]
    return np.sum(coeffs * np.sin(2 * _J * z), axis=-1)


def _series_derivative(coeffs: np.ndarray, zeta):
    """Derivative of :func:`_series` with respect to `zeta`, plus the identity term."""
    z = np.asarray(zeta)[..., np.newaxis]
    return 1 + np.sum(2 * _J * coeffs * np.cos(2 * _J * z), axis=-1)


def _alpha(n: float) -> np.ndarray:
    """Forward Kruger coefficients (Gauss-Schreiber -> Gauss-Kruger), Karney (2011) eq. 35."""
    return np.array([
        n / 2 - 2 * n**2 / 3 + 5 * n**3 / 16 + 41 * n**4 / 180 - 127 * n**5 / 288
        + 7891 * n**6 / 37800,
        13 * n**2 / 48 - 3 * n**3 / 5 + 557 * n**4 / 1440 + 281 * n**5 / 630
        - 1983433 * n**6 / 1935360,
        61 * n**3 / 240 - 103 * n**4 / 140 + 15061 * n**5 / 26880 + 167603 * n**6 / 181440,
        49561 * n**4 / 161280 - 179 * n**5 / 168 + 6601661 * n**6 / 7257600,
        34729 * n**5 / 80640 - 3418889 * n**6 / 1995840,
        212378941 * n**6 / 319334400,
    ])


def _beta(n: float) -> np.ndarray:
    """Reverse Kruger coefficients (Gauss-Kruger -> Gauss-Schreiber), Karney (2011) eq. 38."""
    return np.array([
        n / 2 - 2 * n**2 / 3 + 37 * n**3 / 96 - n**4 / 360 - 81 * n**5 / 512
        + 96199 * n**6 / 604800,
        n**2 / 48 + n**3 / 15 - 437 * n**4 / 1440 + 46 * n**5 / 105
        - 1118711 * n**6 / 3870720,
        17 * n**3 / 480 - 37 * n**4 / 840 - 209 * n**5 / 4480 + 5569 * n**6 / 90720,
        4397 * n**4 / 161280 - 11 * n**5 / 504 - 830251 * n**6 / 7257600,
        4583 * n**5 / 161280 - 108847 * n**6 / 3991680,
        20648693 * n**6 / 638668800,
    ])


class KrugerTM:
    """
    Ellipsoid-dependent constants and primitives for the Kruger Transverse Mercator series.

    Parameters
    ----------
    a : float. Semi-major axis of the ellipsoid (meters).
    b : float. Semi-minor axis of the ellipsoid (meters).
    """

    def __init__(self, a: float, b: float):
        self.a = float(a)
        self.b = float(b)
        self.f = (self.a - self.b) / self.a
        self.e2 = self.f * (2 - self.f)
        self.e = np.sqrt(self.e2)
        self.n = self.f / (2 - self.f)
        self.alpha = _alpha(self.n)
        self.beta = _beta(self.n)
        # Rectifying radius: meridian quadrant = A1 * pi / 2. Karney (2011) eq. 14.
        self.A1 = self.a / (1 + self.n) * (
            1 + self.n**2 / 4 + self.n**4 / 64 + self.n**6 / 256
        )

    def taupf(self, tau):
        """tan(conformal latitude) from tan(geodetic latitude). Karney (2011) eq. 7."""
        tau1 = np.hypot(1.0, tau)
        sig = np.sinh(self.e * np.arctanh(self.e * tau / tau1))
        return np.hypot(1.0, sig) * tau - sig * tau1

    def tauf(self, taup):
        """
        tan(geodetic latitude) from tan(conformal latitude), by Newton's method.

        Converges quadratically from the starting guess below; five iterations reach
        double precision for any terrestrial ellipsoid.
        """
        e2m = 1 - self.e2
        tau = taup / e2m
        for _ in range(5):
            taupa = self.taupf(tau)
            tau = tau + (taup - taupa) * (1 + e2m * tau**2) / (
                e2m * np.hypot(1.0, tau) * np.hypot(1.0, taupa)
            )
        return tau

    def _xi0(self, lat0_rad: float) -> float:
        """Rectifying-plane ordinate of the natural origin (the point lat = lat0 on the CM)."""
        xip0 = np.arctan(self.taupf(np.tan(lat0_rad)))
        return float(xip0 + _series(self.alpha, xip0))

    def forward(self, dlon_rad, lat_rad) -> Tuple[np.ndarray, np.ndarray,
                                                  np.ndarray, np.ndarray]:
        """
        Forward projection factors for a point given by geodetic coordinates.

        Parameters
        ----------
        dlon_rad : float or numpy.ndarray. Longitude relative to the central meridian (radians).
        lat_rad : float or numpy.ndarray. Geodetic latitude (radians).

        Returns
        -------
        xi, eta : Rectifying-plane coordinates, normalised by `A1`.
        gamma : Meridian convergence (radians), for `k0 = 1` (it is independent of `k0`).
        k : Point scale factor for `k0 = 1`.
        """
        tau = np.tan(lat_rad)
        taup = self.taupf(tau)
        sin_lam, cos_lam = np.sin(dlon_rad), np.cos(dlon_rad)
        denom = np.hypot(taup, cos_lam)

        # Gauss-Schreiber projection: transverse Mercator of the conformal sphere.
        zeta_p = np.arctan2(taup, cos_lam) + 1j * np.arcsinh(sin_lam / denom)
        gamma_p = np.arctan2(sin_lam * taup, cos_lam * np.hypot(1.0, taup))
        k_p = np.sqrt(1 - self.e2 * np.sin(lat_rad)**2) * np.hypot(1.0, tau) / denom

        # Kruger correction from Gauss-Schreiber to Gauss-Kruger, plus its derivative.
        zeta = zeta_p + _series(self.alpha, zeta_p)
        dzeta = _series_derivative(self.alpha, zeta_p)

        gamma = gamma_p - np.angle(dzeta)
        k = k_p * (self.A1 / self.a) * np.abs(dzeta)
        return np.real(zeta), np.imag(zeta), gamma, k

    def reverse(self, xi, eta) -> Tuple[np.ndarray, np.ndarray]:
        """
        Geodetic coordinates from rectifying-plane coordinates.

        Parameters
        ----------
        xi, eta : float or numpy.ndarray. Rectifying-plane coordinates, normalised by `A1`.

        Returns
        -------
        dlon_rad : Longitude relative to the central meridian (radians).
        lat_rad : Geodetic latitude (radians).
        """
        zeta = xi + 1j * eta
        zeta_p = zeta - _series(self.beta, zeta)
        xi_p, eta_p = np.real(zeta_p), np.imag(zeta_p)

        # Inverse Gauss-Schreiber.
        s, c = np.sinh(eta_p), np.cos(xi_p)
        r = np.hypot(s, c)
        dlon_rad = np.arctan2(s, c)
        lat_rad = np.arctan(self.tauf(np.sin(xi_p) / r))
        return dlon_rad, lat_rad


@lru_cache(maxsize=32)
def _tm(a: float, b: float) -> KrugerTM:
    """Cached :class:`KrugerTM`, so repeated calls do not rebuild the coefficients."""
    return KrugerTM(a, b)


def tm_factors_geographic(a: float, b: float, dlon_rad, lat_rad,
                          k0: float = 1.0) -> Tuple[Union[float, np.ndarray],
                                                    Union[float, np.ndarray]]:
    """
    Meridian convergence (radians) and point scale factor from geodetic coordinates.

    `dlon_rad` is the longitude relative to the central meridian.
    """
    tm = _tm(float(a), float(b))
    _, _, gamma, k = tm.forward(dlon_rad, lat_rad)
    return gamma, k0 * k


def tm_factors_projected(a: float, b: float, x, y, lat0_rad: float = 0.0,
                         false_easting: float = 0.0, false_northing: float = 0.0,
                         k0: float = 1.0) -> Tuple[Union[float, np.ndarray],
                                                   Union[float, np.ndarray]]:
    """
    Meridian convergence (radians) and point scale factor from projected coordinates.

    The grid coordinates are unprojected with the reverse Kruger series and the factors
    are then evaluated with the forward series, so the two directions cannot disagree.
    """
    tm = _tm(float(a), float(b))
    xi = (y - false_northing) / (k0 * tm.A1) + tm._xi0(lat0_rad)
    eta = (x - false_easting) / (k0 * tm.A1)
    dlon_rad, lat_rad = tm.reverse(xi, eta)
    _, _, gamma, k = tm.forward(dlon_rad, lat_rad)
    return gamma, k0 * k
