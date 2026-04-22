"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com


Polar Stereographic projection (EPSG 9810, Variant A).

Implements the forward and inverse Polar Stereographic projection used
by UPS (Universal Polar Stereographic). For UPS itself, see
:mod:`pygeodetics.projections.UPS`, which configures this engine with
the standard UPS defaults.

Formulas follow IOGP "Geomatics Guidance Note 7-2" and Snyder
"Map Projections - A Working Manual" (USGS 1395).
"""

from __future__ import annotations

from typing import Literal

import numpy as np


class PolarStereographic:
    """
    Polar Stereographic projection (Variant A) about the geographic pole.

    Parameters
    ----------
    pole : {"N", "S"}
        Which pole the projection is centred on. ``"N"`` corresponds to
        latitude of natural origin = +90°, ``"S"`` to -90°.
    central_meridian_deg : float
        Longitude of natural origin in degrees.
    scale_factor : float
        Scale factor at the pole (k0). For UPS this is 0.994.
    false_easting : float
        False easting in metres.
    false_northing : float
        False northing in metres.
    a : float
        Semi-major axis of the ellipsoid in metres.
    f : float
        Flattening of the ellipsoid.
    """

    def __init__(self, pole: Literal["N", "S"], central_meridian_deg: float,
                 scale_factor: float, false_easting: float,
                 false_northing: float, a: float, f: float):
        pole = str(pole).upper()
        if pole not in ("N", "S"):
            raise ValueError("pole must be 'N' or 'S'.")
        self.pole = pole
        self.lon0 = np.radians(central_meridian_deg)
        self.k0 = float(scale_factor)
        self.fe = float(false_easting)
        self.fn = float(false_northing)
        self.a = float(a)
        self.f = float(f)
        self.e2 = self.f * (2.0 - self.f)
        self.e = np.sqrt(self.e2)
        # Common scaling constant: sqrt((1+e)^(1+e) * (1-e)^(1-e))
        self._k_const = np.sqrt(
            (1.0 + self.e) ** (1.0 + self.e)
            * (1.0 - self.e) ** (1.0 - self.e)
        )

    def forward(self, lat_deg, lon_deg):
        """
        Project geodetic (lat, lon) in degrees to (easting, northing).

        Accepts scalars or NumPy arrays; output matches input shape.
        """
        lat = np.radians(np.asarray(lat_deg, dtype=float))
        lon = np.radians(np.asarray(lon_deg, dtype=float))

        if self.pole == "N":
            # EPSG 9810 north-pole formula
            esp = self.e * np.sin(lat)
            t = np.tan(np.pi / 4.0 - lat / 2.0) * np.power(
                (1.0 + esp) / (1.0 - esp), self.e / 2.0
            )
            rho = 2.0 * self.a * self.k0 * t / self._k_const
            theta = lon - self.lon0
            dE = rho * np.sin(theta)
            dN = -rho * np.cos(theta)
        else:
            # EPSG 9810 south-pole formula
            esp = self.e * np.sin(lat)
            t = np.tan(np.pi / 4.0 + lat / 2.0) / np.power(
                (1.0 + esp) / (1.0 - esp), self.e / 2.0
            )
            rho = 2.0 * self.a * self.k0 * t / self._k_const
            theta = lon - self.lon0
            dE = rho * np.sin(theta)
            dN = rho * np.cos(theta)

        return self.fe + dE, self.fn + dN

    def inverse(self, easting, northing):
        """
        Convert (easting, northing) back to geodetic (lat, lon) in degrees.
        """
        dE = np.asarray(easting, dtype=float) - self.fe
        dN = np.asarray(northing, dtype=float) - self.fn

        rho = np.hypot(dE, dN)
        # Avoid division by zero exactly at the pole.
        with np.errstate(divide="ignore", invalid="ignore"):
            t = rho * self._k_const / (2.0 * self.a * self.k0)

            # Conformal latitude χ.
            if self.pole == "N":
                chi = np.pi / 2.0 - 2.0 * np.arctan(t)
            else:
                chi = 2.0 * np.arctan(t) - np.pi / 2.0

            # Series inverse from conformal to geodetic latitude.
            e2 = self.e2
            e4 = e2 * e2
            e6 = e4 * e2
            e8 = e4 * e4
            lat = (
                chi
                + (e2 / 2.0 + 5.0 * e4 / 24.0 + e6 / 12.0 + 13.0 * e8 / 360.0)
                * np.sin(2.0 * chi)
                + (7.0 * e4 / 48.0 + 29.0 * e6 / 240.0 + 811.0 * e8 / 11520.0)
                * np.sin(4.0 * chi)
                + (7.0 * e6 / 120.0 + 81.0 * e8 / 1120.0) * np.sin(6.0 * chi)
                + (4279.0 * e8 / 161280.0) * np.sin(8.0 * chi)
            )

            if self.pole == "N":
                lon = self.lon0 + np.arctan2(dE, -dN)
            else:
                lon = self.lon0 + np.arctan2(dE, dN)

            # At the pole rho == 0 → longitude is undefined; force to lon0.
            at_pole = rho == 0.0
            if np.ndim(lon) == 0:
                if at_pole:
                    lon = self.lon0
            else:
                lon = np.where(at_pole, self.lon0, lon)

        # Wrap longitude into (-180, 180].
        lon_deg = np.degrees(lon)
        lon_deg = ((lon_deg + 180.0) % 360.0) - 180.0
        return np.degrees(lat), lon_deg
