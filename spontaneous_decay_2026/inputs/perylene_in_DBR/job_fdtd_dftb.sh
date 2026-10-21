#!/bin/bash
#SBATCH --partition=taoeli
#SBATCH --job-name=c-2.634eV
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=100:00:00
##SBATCH --mem-per-cpu=7G
#SBATCH --mem=32G
#SBATCH --output=job.out
#SBATCH --error=job.err

source ~/opt/conda/bin/activate
conda activate mxl
export PATH="$HOME/.local/bin:$PATH"
export DFTBPLUS_BIN="$HOME/.local/bin/dftb+"
export LD_LIBRARY_PATH="$HOME/opt/conda/envs/dftb-build/lib:$LD_LIBRARY_PATH"

python -u fdtd.py > out-maxwelllink
python -u fit_energy_decay_rate.py
