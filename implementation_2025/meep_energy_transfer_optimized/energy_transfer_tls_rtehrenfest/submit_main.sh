#!/bin/bash

#SBATCH --nodes=1         # Total # of nodes (must be 1 for serial job)
#SBATCH --ntasks=16       # Total # of MPI tasks (should be 1 for serial job)
#SBATCH --time=2-00:00:00    # Total run time limit (hh:mm:ss)
#SBATCH -J etm_rteh     # Job name
#SBATCH -o main.o%j      # Name of stdout output file
#SBATCH -e main.e%j      # Name of stderr error file
#SBATCH -p shared  # Queue (partition) name


mpirun -np 16 python meep_energy_transfer_3d.py


