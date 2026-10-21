#!/bin/bash
#SBATCH --job-name=r-1-56
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --time=144:00:00
#SBATCH --mem-per-cpu=5G
##SBATCH --mem=120G
#SBATCH --partition=taoeli
#SBATCH --output=job.out
#SBATCH --error=job.err

source ~/opt/conda/bin/activate
conda activate mxl

unset PYTHONPATH
export OMP_NUM_THREADS=1
export MPIR_CVAR_CH4_SHM_ENABLE=1

echo "Node: $(hostname)"
echo "SLURM_NTASKS=$SLURM_NTASKS"
echo "SLURM_CPUS_ON_NODE=$SLURM_CPUS_ON_NODE"

echo "CPU affinity:"
$CONDA_PREFIX/bin/mpirun -bind-to core -np "$SLURM_NTASKS" bash -c 'echo "rank=$PMI_RANK pid=$$"; taskset -pc $$'
echo "Start Meep:"
$CONDA_PREFIX/bin/mpirun -bind-to core -np "$SLURM_NTASKS" $CONDA_PREFIX/bin/python npom_mb_scattering_spectrum.py -u > out-maxwelllink
