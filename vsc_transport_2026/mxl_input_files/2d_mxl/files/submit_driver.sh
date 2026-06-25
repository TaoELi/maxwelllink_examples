#!/bin/bash
# following 3 commands required by OMP(OpenMP) only within one node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=mxl_meso
#SBATCH --partition=standard
#SBATCH --time=7-00:00:00
#SBATCH --array=0-35
#SBATCH --output=./driver_%A_%a.out
#SBATCH --error=./driver_%A_%a.err

# 0. wait for the main job to write the bridge manifest, then start this bridge node
traj=$SLURM_ARRAY_TASK_ID
export OMP_NUM_THREADS=1

# Wait for the main job to write the shared manifest.
until [[ -f aggregation.json ]]; do sleep 10; done

# Start this node's bridge: one upstream TCP link, one local UNIX hub.
mxl_bridge --info aggregation.json --idx ${traj} &

HOST="bridge_${traj}"
PORT=32435
SOCKET_PATH="/tmp/socketmxl_${HOST}"
until [[ -S "$SOCKET_PATH" ]]; do sleep 10; done

# launch LAMMPS molecular drivers
n_molecule_start=$((576*traj + 1))
n_molecule_end=$((576*traj + 576))
# launch many drivers in the background if n_molecule > 1
for ((i=n_molecule_start; i<n_molecule_end; i++)); do
    echo "Launching LAMMPS driver for molecule $i"
    sh launch_lmp_xml.sh $HOST $PORT $i &
    sleep 0.1s  # slight delay to avoid overwhelming the hub
done

# launch the last driver in the foreground
sh launch_lmp_xml.sh $HOST $PORT $n_molecule_end
#sh launch_lmp_xml.sh $HOST $PORT $traj
