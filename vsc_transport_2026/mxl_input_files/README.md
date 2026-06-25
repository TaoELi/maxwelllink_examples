# MXL Input Files

This folder contains input templates and helper scripts for running MXL simulations with LAMMPS and Maxwell-link style cavity coupling.

## Contents

- `1d_mxl/`: Input and submission files for one-dimensional nonequilibrium MXL simulations.
- `1d_mxl_eq/`: Input and submission files for one-dimensional equilibrium simulations.
- `1d_mxl_ndependence/`: Input and analysis files for studying system-size or molecule-number dependence in 1D simulations.
- `2d_mxl/`: Input and submission files for two-dimensional MXL simulations.

Each case directory contains job submission scripts, launch scripts, MXL driver files, and LAMMPS input files under `files/`.  Some directories also include analysis helpers such as `coenergy.py`, `analysis.sh`, and `run.sh`.
