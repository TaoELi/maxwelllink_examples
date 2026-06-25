# Plotting Scripts

This folder contains the plotting and post-processing scripts used to generate the figures for the MXL simulation results.  The scripts read precomputed simulation data from `data/`, generate publication-style plots and animations, and write the resulting files to `figures/`.

## Contents

- `plotting.py`: Main figure-generation script for 1D and 2D electric-field plots, co-energy plots, mean-square-displacement analysis, velocity-group analysis, nonequilibrium spectra, and animations.
- `simulate_1d_symmetry_breaking.py`: Standalone model for a one-dimensional polariton spectrum with spatial translational-symmetry breaking.  It can generate both spectrum plots and the corresponding model data.
- `columnplots.py`: Local plotting utilities for figure layout, axis styling, legends, and image placement.
- `data/`: Input simulation and model data, organized by system size, scan parameters, co-energy data, and number-dependence studies.
- `figures/`: Generated figures, movies, and presentation-ready outputs.

Most scripts assume they are run from inside this directory so that relative paths such as `./data/...` and `./figures/...` resolve correctly.
