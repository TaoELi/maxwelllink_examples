#!/bin/bash
# following 3 commands required by OMP(OpenMP) only within one node
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=mxl_meso
#SBATCH --partition=taoeli
#SBATCH --time=7-00:00:00
#SBATCH --array=1-36
#SBATCH --output=./driver_%A_%a.out
#SBATCH --error=./driver_%A_%a.err

# 0. wait for a few seconds to ensure the main job is running and the tcp host and port info file is created
traj=$SLURM_ARRAY_TASK_ID
sleep $(echo "10 + 1 * $traj" | bc)s
HOST_PORT_FILE="tcp_host_port_info.txt"
if [[ ! -f "$HOST_PORT_FILE" ]]; then
    echo "Error: Host and port info file '$HOST_PORT_FILE' not found!"
    exit 1
fi
HOST=$(sed -n '1p' "$HOST_PORT_FILE")
PORT=$(sed -n '2p' "$HOST_PORT_FILE")

# launch tls driver
n_molecule_start=$((4*traj - 3))
n_molecule_end=$((4*traj))
# launch many drivers in the background if n_molecule > 1
for ((i=n_molecule_start; i<n_molecule_end; i++)); do
    echo "Launching tls_driver for molecule $i"
    sh launch_lmp_xml.sh $HOST $PORT $i &
    sleep 1s  # slight delay to avoid overwhelming the hub
done

# launch the last driver in the foreground
sh launch_lmp_xml.sh $HOST $PORT $n_molecule_end
#sh launch_lmp_xml.sh $HOST $PORT $traj
