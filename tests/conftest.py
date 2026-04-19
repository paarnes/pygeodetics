"""
Author: Per Helge Aarnes
Email: per.helge.aarnes@gmail.com

Pytest configuration. Adds the `src/` directory to `sys.path` so the
`pygeodetics` package is importable when tests are run from the repository
root without first installing the package.
"""

import sys
import os

sys.path.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")),
)
