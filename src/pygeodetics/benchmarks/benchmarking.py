import timeit

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../geodetics')))
from ECEF2geod import ECEF2geod, ECEF2geodb, ECEF2geodv

# Define parameters
a = 6378137.0  # Semi-major axis (meters)
b = 6356752.314245  # Semi-minor axis (meters)
X, Y, Z = 3149785.9652, 598260.8822, 5495348.4927
n_sim = 10000 # Number of iterations

# Benchmark ECEF2geod
time_geod = timeit.timeit(lambda: ECEF2geod(a, b, X, Y, Z), number=n_sim)

# Benchmark ECEF2geodb
time_geodb = timeit.timeit(lambda: ECEF2geodb(a, b, X, Y, Z), number=n_sim)

# Benchmark ECEF2geodv
time_geodv = timeit.timeit(lambda: ECEF2geodv(a, b, X, Y, Z), number=n_sim)

print(f"ECEF2geod (iteration): {time_geod:.6f} seconds for {n_sim} iterations")
print(f"ECEF2geodb (Bowring): {time_geodb:.6f} seconds for {n_sim} iterations")
print(f"ECEF2geodv (Vermeille): {time_geodv:.6f} seconds for {n_sim} iterations")
