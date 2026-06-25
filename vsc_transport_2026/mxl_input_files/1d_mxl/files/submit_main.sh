#!/bin/bash
# following 3 commands required by OMP(OpenMP) only within one node
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=mxl_meso
#SBATCH --partition=standard
#SBATCH --time=7-00:00:00
#SBATCH --output=./nve_%A_%a.out
#SBATCH --error=./nve_%A_%a.err

ulimit -n 30000
python -u multmodes_vsc_hpc.py

