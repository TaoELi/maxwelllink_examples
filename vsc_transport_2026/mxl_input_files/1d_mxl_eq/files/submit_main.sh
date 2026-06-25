#!/bin/bash
# following 3 commands required by OMP(OpenMP) only within one node
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=mxl_meso
#SBATCH --partition=taoeli
#SBATCH --time=7-00:00:00
#SBATCH --output=./nve_%A_%a.out
#SBATCH --error=./nve_%A_%a.err


python -u multmodes_vsc_hpc.py

