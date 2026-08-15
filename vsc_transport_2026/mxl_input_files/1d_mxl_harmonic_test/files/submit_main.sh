#!/bin/bash
# following 3 commands required by OMP(OpenMP) only within one node
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=mxl_meso
#SBATCH --partition=shared
#SBATCH --time=1-00:00:00
#SBATCH --output=./nve_%A_%a.out
#SBATCH --error=./nve_%A_%a.err

ulimit -n 30000
python -u multmodes_vsc_hpc.py

