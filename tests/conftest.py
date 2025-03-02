"""
Author: Per Helge Aarnes
Email: per.helge.aarnes@gmail.com

This file is automatically loaded by pytest and is used to configure the test environment.
In this project, it ensures that the `src` directory is added to the Python module search path (`sys.path`),
allowing tests to import modules from the `src` folder without individual path modifications.

Notes:
- The `src` directory is appended to `sys.path` dynamically so all tests can access modules like `Vincentys`.
- Avoid hardcoding imports or making changes to individual test files for path resolution.
"""

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/pygeodetics')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/pygeodetics/geodetics')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/pygeodetics/projections')))