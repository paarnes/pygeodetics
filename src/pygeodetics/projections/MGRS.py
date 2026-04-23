"""
author: Per Helge Aarnes
email: per.helge.aarnes@gmail.com


Military Grid Reference System (MGRS) encoding and parsing.

This module wires together the existing UTM and UPS wrappers to provide
turn-key MGRS string encoding and parsing. No new geodetic math is
implemented here:

- UTM regions (latitudes -80°..84°) are handled via
  :func:`pygeodetics.geodetic_to_utm` / :func:`pygeodetics.utm_to_geodetic`.
- Polar regions (latitudes <-80° or >84°) are handled via
  :func:`pygeodetics.geodetic_to_ups` / :func:`pygeodetics.ups_to_geodetic`.

MGRS string layout
------------------
For UTM regions::

    <zone:1-2 digits><band:1 letter><col:1 letter><row:1 letter><E><N>

For UPS (polar) regions::

    <pole-zone:1 letter A/B/Y/Z><col:1 letter><row:1 letter><E><N>

``<E>`` and ``<N>`` are zero-padded numeric strings whose length equals
the requested precision (1..5 digits per axis):

- precision=5 → 1 m
- precision=4 → 10 m
- precision=3 → 100 m
- precision=2 → 1 km
- precision=1 → 10 km
- precision=0 → 100 km square only (no easting/northing digits)
"""

from __future__ import annotations

import re
from typing import Optional, Tuple, Union

import numpy as np

from ..Ellipsoid import Ellipsoid, WGS84
from .UTM import geodetic_to_utm, utm_to_geodetic, mgrs_band_letter
from .UPS import geodetic_to_ups, ups_to_geodetic


#
# UTM column letters depend on (zone-1) % 3 ("set"):
#   set 0 (zones 1,4,7,...): A..H  for easting columns 1..8 (100k..800k)
#   set 1 (zones 2,5,8,...): J..R  (skipping I and O)
#   set 2 (zones 3,6,9,...): S..Z
_UTM_COL_ALPHABETS = ("ABCDEFGH", "JKLMNPQR", "STUVWXYZ")

# UTM row letters: 20 letters cycling every 2 000 000 m of northing.
# Odd zones start at A (northing 0 → A); even zones start at F.
_UTM_ROW_ALPHABET = "ABCDEFGHJKLMNPQRSTUV"

# UPS 100 km column / row alphabets (NGA / GeographicLib convention).
# Polar zones:
#   south: A (lon < 0), B (lon >= 0)
#   north: Y (lon < 0), Z (lon >= 0)
_UPS_COL_LEFT = "JKLPQRSTUXYZ"   # zones A and Y, easting starts at  800 000
_UPS_COL_RIGHT = "ABCFGHJKLPQR"  # zones B and Z, easting starts at 2 000 000

# Row alphabet for SOUTH polar (24 letters), starting at northing 800 000.
_UPS_ROW_SOUTH = "ABCDEFGHJKLMNPQRSTUVWXYZ"
# Row alphabet for NORTH polar (14 letters), starting at northing 1 300 000.
_UPS_ROW_NORTH = "ABCDEFGHJKLMNP"


_MGRS_RE_UTM = re.compile(
    r"^\s*(?P<zone>\d{1,2})(?P<band>[C-HJ-NP-X])"
    r"(?P<col>[A-HJ-NP-Z])(?P<row>[A-HJ-NP-V])"
    r"(?P<digits>\d*)\s*$"
)
_MGRS_RE_UPS = re.compile(
    r"^\s*(?P<pole>[ABYZ])"
    r"(?P<col>[A-HJ-NP-Z])(?P<row>[A-HJ-NP-Z])"
    r"(?P<digits>\d*)\s*$"
)



def _utm_col_letter(zone: int, easting: float) -> str:
    alphabet = _UTM_COL_ALPHABETS[(zone - 1) % 3]
    idx = int(easting // 100_000) - 1
    if idx < 0 or idx >= len(alphabet):
        raise ValueError(
            f"UTM easting {easting:.3f} out of MGRS column range "
            f"for zone {zone}."
        )
    return alphabet[idx]


def _utm_row_letter(zone: int, northing: float) -> str:
    # For southern hemisphere, the false-northing offset of 10 000 000 m
    # is exactly 5 cycles of 2 000 000 m, so the same modulo works.
    row_idx = int((northing % 2_000_000) // 100_000)
    if zone % 2 == 0:
        row_idx = (row_idx + 5) % 20
    return _UTM_ROW_ALPHABET[row_idx]


def _utm_easting_from_col(zone: int, col_letter: str) -> float:
    alphabet = _UTM_COL_ALPHABETS[(zone - 1) % 3]
    if col_letter not in alphabet:
        raise ValueError(
            f"Invalid MGRS column letter '{col_letter}' for zone {zone}."
        )
    return (alphabet.index(col_letter) + 1) * 100_000.0


def _utm_northing_from_row(zone: int, row_letter: str,
                           band_letter: str) -> float:
    if row_letter not in _UTM_ROW_ALPHABET:
        raise ValueError(f"Invalid MGRS row letter '{row_letter}'.")
    row_idx = _UTM_ROW_ALPHABET.index(row_letter)
    if zone % 2 == 0:
        row_idx = (row_idx - 5) % 20
    base = row_idx * 100_000.0  # northing within 0..2 000 000 cycle

    if band_letter not in "CDEFGHJKLMNPQRSTUVWX":
        raise ValueError(f"Invalid MGRS band letter '{band_letter}'.")
    band_idx = "CDEFGHJKLMNPQRSTUVWX".index(band_letter)
    band_lat_lo = -80.0 + 8.0 * band_idx
    # Band X is 12° tall (72°..84°); all other bands are 8° tall.
    band_lat_hi = 84.0 if band_letter == "X" else band_lat_lo + 8.0
    hemi = "N" if band_lat_lo >= 0.0 else "S"
    cm = -177.0 + 6.0 * (zone - 1)

    # Try each candidate cycle of 2 000 000 m and pick the one whose
    # recovered latitude lies inside the band. This is exact even for
    # the 12°-tall band X.
    n_cycles = 5 if hemi == "N" else 5  # northings cycle 0..1e7
    best_n = None
    best_score = float("inf")
    for k in range(n_cycles + 1):
        candidate = base + k * 2_000_000.0
        if candidate < 0.0 or candidate > 10_000_000.0:
            continue
        try:
            cand_lat, _ = utm_to_geodetic(500_000.0, candidate,
                                          zone=zone, hemisphere=hemi)
        except Exception:
            continue
        if band_lat_lo - 1e-6 <= cand_lat <= band_lat_hi + 1e-6:
            return candidate
        # Fall-back score: distance to band centre.
        score = abs(cand_lat - 0.5 * (band_lat_lo + band_lat_hi))
        if score < best_score:
            best_score = score
            best_n = candidate
    if best_n is None:
        raise ValueError(
            f"Could not recover northing for zone {zone}, "
            f"row '{row_letter}', band '{band_letter}'."
        )
    return best_n



def _ups_pole_zone(lat: float, lon: float) -> str:
    if lat >= 0.0:
        return "Z" if lon >= 0.0 else "Y"
    return "B" if lon >= 0.0 else "A"


def _ups_col_letter(pole_zone: str, easting: float) -> str:
    if pole_zone in ("A", "Y"):
        alphabet = _UPS_COL_LEFT
        base = 800_000.0
    else:
        alphabet = _UPS_COL_RIGHT
        base = 2_000_000.0
    idx = int((easting - base) // 100_000)
    if idx < 0 or idx >= len(alphabet):
        raise ValueError(
            f"UPS easting {easting:.3f} out of MGRS column range "
            f"for zone {pole_zone}."
        )
    return alphabet[idx]


def _ups_row_letter(pole_zone: str, northing: float) -> str:
    if pole_zone in ("A", "B"):
        alphabet = _UPS_ROW_SOUTH
        base = 800_000.0
    else:
        alphabet = _UPS_ROW_NORTH
        base = 1_300_000.0
    idx = int((northing - base) // 100_000)
    if idx < 0 or idx >= len(alphabet):
        raise ValueError(
            f"UPS northing {northing:.3f} out of MGRS row range "
            f"for zone {pole_zone}."
        )
    return alphabet[idx]


def _ups_easting_from_col(pole_zone: str, col_letter: str) -> float:
    if pole_zone in ("A", "Y"):
        alphabet = _UPS_COL_LEFT
        base = 800_000.0
    else:
        alphabet = _UPS_COL_RIGHT
        base = 2_000_000.0
    if col_letter not in alphabet:
        raise ValueError(
            f"Invalid UPS column letter '{col_letter}' for zone {pole_zone}."
        )
    return base + alphabet.index(col_letter) * 100_000.0


def _ups_northing_from_row(pole_zone: str, row_letter: str) -> float:
    if pole_zone in ("A", "B"):
        alphabet = _UPS_ROW_SOUTH
        base = 800_000.0
    else:
        alphabet = _UPS_ROW_NORTH
        base = 1_300_000.0
    if row_letter not in alphabet:
        raise ValueError(
            f"Invalid UPS row letter '{row_letter}' for zone {pole_zone}."
        )
    return base + alphabet.index(row_letter) * 100_000.0



def to_mgrs(
    lat: float,
    lon: float,
    precision: int = 5,
    ellipsoid: Optional[Ellipsoid] = None,
) -> str:
    """
    Encode a geodetic position as an MGRS string.

    Parameters
    ----------
    lat, lon : float
        Geodetic latitude / longitude in degrees.
    precision : int, default 5
        Number of digits per axis (0..5):

        - 5 → 1 m
        - 4 → 10 m
        - 3 → 100 m
        - 2 → 1 km
        - 1 → 10 km
        - 0 → 100 km square only

    ellipsoid : Ellipsoid, optional
        Reference ellipsoid. Defaults to WGS84.

    Returns
    -------
    str
        MGRS string with no internal spaces, e.g. ``"32VNM9760352702"``.
    """
    if not (0 <= int(precision) <= 5):
        raise ValueError("precision must be an integer in [0, 5].")
    precision = int(precision)

    lat_f = float(lat)
    lon_f = float(lon)

    if -80.0 <= lat_f <= 84.0:
        # UTM region.
        e, n, zone, band = geodetic_to_utm(lat_f, lon_f, ellipsoid=ellipsoid)
        col = _utm_col_letter(zone, e)
        row = _utm_row_letter(zone, n)
        gzd = f"{zone:d}{band}"
    else:
        # UPS (polar) region.
        e, n, hemi = geodetic_to_ups(lat_f, lon_f, ellipsoid=ellipsoid)
        pole_zone = _ups_pole_zone(lat_f, lon_f)
        col = _ups_col_letter(pole_zone, e)
        row = _ups_row_letter(pole_zone, n)
        gzd = pole_zone

    if precision == 0:
        return f"{gzd}{col}{row}"

    scale = 10 ** (5 - precision)
    sub_e = int(e) % 100_000 // scale
    sub_n = int(n) % 100_000 // scale
    width = precision
    return f"{gzd}{col}{row}{sub_e:0{width}d}{sub_n:0{width}d}"


def from_mgrs(
    mgrs_string: str,
    ellipsoid: Optional[Ellipsoid] = None,
) -> Tuple[float, float]:
    """
    Parse an MGRS string into geodetic (lat, lon) in degrees.

    The returned position is the south-west corner of the MGRS cell at
    the given precision. For a center-of-cell position, add half the
    cell size to easting and northing before re-projecting.

    Parameters
    ----------
    mgrs_string : str
        MGRS reference. Spaces inside the string are ignored.
    ellipsoid : Ellipsoid, optional
        Reference ellipsoid. Defaults to WGS84.

    Returns
    -------
    lat, lon : float
        Geodetic latitude and longitude in degrees.

    Raises
    ------
    ValueError
        If the string is not a syntactically valid MGRS reference.
    """
    if not isinstance(mgrs_string, str):
        raise TypeError("mgrs_string must be a string.")

    cleaned = mgrs_string.replace(" ", "").upper()
    if not cleaned:
        raise ValueError("Empty MGRS string.")

    m_utm = _MGRS_RE_UTM.match(cleaned)
    m_ups = _MGRS_RE_UPS.match(cleaned)

    if m_utm:
        digits = m_utm.group("digits")
        if len(digits) % 2 != 0 or len(digits) > 10:
            raise ValueError(
                f"MGRS numeric part must have an even length 0..10 "
                f"(got {len(digits)}: '{digits}')."
            )
        precision = len(digits) // 2
        scale = 10 ** (5 - precision) if precision > 0 else 100_000

        zone = int(m_utm.group("zone"))
        band = m_utm.group("band")
        col = m_utm.group("col")
        row = m_utm.group("row")

        if not (1 <= zone <= 60):
            raise ValueError(f"UTM zone {zone} out of range 1..60.")

        e_base = _utm_easting_from_col(zone, col)
        n_base = _utm_northing_from_row(zone, row, band)

        if precision > 0:
            sub_e = int(digits[:precision]) * scale
            sub_n = int(digits[precision:]) * scale
        else:
            sub_e = sub_n = 0

        easting = e_base + sub_e
        northing = n_base + sub_n
        hemi = "N" if band >= "N" else "S"
        return utm_to_geodetic(easting, northing, zone=zone,
                               hemisphere=hemi, ellipsoid=ellipsoid)

    if m_ups:
        digits = m_ups.group("digits")
        if len(digits) % 2 != 0 or len(digits) > 10:
            raise ValueError(
                f"MGRS numeric part must have an even length 0..10 "
                f"(got {len(digits)}: '{digits}')."
            )
        precision = len(digits) // 2
        scale = 10 ** (5 - precision) if precision > 0 else 100_000

        pole = m_ups.group("pole")
        col = m_ups.group("col")
        row = m_ups.group("row")

        e_base = _ups_easting_from_col(pole, col)
        n_base = _ups_northing_from_row(pole, row)

        if precision > 0:
            sub_e = int(digits[:precision]) * scale
            sub_n = int(digits[precision:]) * scale
        else:
            sub_e = sub_n = 0

        easting = e_base + sub_e
        northing = n_base + sub_n
        hemi = "N" if pole in ("Y", "Z") else "S"
        return ups_to_geodetic(easting, northing, hemisphere=hemi,
                               ellipsoid=ellipsoid)

    raise ValueError(f"Malformed MGRS string: '{mgrs_string}'.")



def utm_epsg(zone: int, hemisphere: str, datum: str = "WGS84") -> int:
    """
    Return the EPSG code for a WGS84 / NAD83 UTM zone.

    Parameters
    ----------
    zone : int
        UTM zone number (1..60).
    hemisphere : {"N", "S"}
        Hemisphere of the zone.
    datum : {"WGS84", "NAD83"}, default "WGS84"
        Geodetic datum. NAD83 is northern hemisphere only.
    """
    if not (1 <= int(zone) <= 60):
        raise ValueError("zone must be in 1..60.")
    hemisphere = str(hemisphere).upper()
    if hemisphere not in ("N", "S"):
        raise ValueError("hemisphere must be 'N' or 'S'.")
    datum = datum.upper()
    if datum == "WGS84":
        return (32600 if hemisphere == "N" else 32700) + int(zone)
    if datum == "NAD83":
        if hemisphere != "N":
            raise ValueError("NAD83 UTM zones are northern hemisphere only.")
        return 26900 + int(zone)
    raise ValueError(f"Unsupported datum '{datum}' (use 'WGS84' or 'NAD83').")


def ups_epsg(hemisphere: str, datum: str = "WGS84") -> int:
    """
    Return the EPSG code for a WGS84 UPS hemisphere.

    Parameters
    ----------
    hemisphere : {"N", "S"}
    datum : {"WGS84"}, default "WGS84"
    """
    hemisphere = str(hemisphere).upper()
    if hemisphere not in ("N", "S"):
        raise ValueError("hemisphere must be 'N' or 'S'.")
    if datum.upper() != "WGS84":
        raise ValueError("Only WGS84 UPS EPSG codes are defined here.")
    return 32661 if hemisphere == "N" else 32761
