#!/usr/bin/env bash
# Post-create setup for the pygeodetics devcontainer.
# Installs the package in editable mode plus everything needed for
# development, testing, and running the example notebooks.

set -euxo pipefail

python -m pip install --upgrade pip setuptools wheel

# Install the project in editable mode.
python -m pip install -e .

# Test dependencies.
python -m pip install -r tests/requirements.txt

# Notebook + scientific stack used by examples.ipynb and benchmarks.
# ipykernel is pinned to <7 to avoid the traitlets `_control_lock`
# regression that affects ipykernel 7.x.
python -m pip install \
    "ipykernel<7" \
    jupyterlab \
    notebook \
    matplotlib \
    pandas

# Optional cross-check libraries used during development (NOT imported
# by the test suite -- tests use hard-coded reference values).
python -m pip install \
    pyproj \
    pymap3d \
    pygeodesy \
    mgrs

# Developer tooling.
python -m pip install \
    ruff \
    build \
    twine

echo "pygeodetics devcontainer ready."
