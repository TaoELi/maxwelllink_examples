#!/bin/bash
# this job should be run after submit_main.sh and after the tcp host and port info file is created in file

#SBATCH --nodes=1         # Total # of nodes (must be 1 for serial job)
#SBATCH --ntasks=1        # Total # of MPI tasks (should be 1 for serial job)
#SBATCH --time=2-00:00:00    # Total run time limit (hh:mm:ss)
#SBATCH -J tls_se      # Job name
#SBATCH -o driver.o%j      # Name of stdout output file
#SBATCH -e driver.e%j      # Name of stderr error file
#SBATCH -p standard  # Queue (partition) name

# 0. wait for a few seconds to ensure the main job is running and the tcp host and port info file is created
sleep 10s

# 1. read tcp host and port info from file
HOST_PORT_FILE="tcp_host_port_info.txt"
if [[ ! -f "$HOST_PORT_FILE" ]]; then
    echo "Error: Host and port info file '$HOST_PORT_FILE' not found!"
    exit 1
fi
HOST=$(sed -n '1p' "$HOST_PORT_FILE")
PORT=$(sed -n '2p' "$HOST_PORT_FILE")

# launch tls driver
n_molecule=1
# launch many drivers in the background if n_molecule > 1
if (( n_molecule > 1 )); then
    for i in $(seq 1 $((n_molecule-1))); do
        echo "Launching tls_driver for molecule $i"
        mxl_driver --model tls --address $HOST --port $PORT --param "omega=0.242, mu12=187, orientation=2, pe_initial=1e-4" &
        sleep 0.01s  # slight delay to avoid overwhelming the hub
    done
fi

# launch the last driver in the foreground
mxl_driver --model tls --address $HOST --port $PORT --param "omega=0.242, mu12=187, orientation=2, pe_initial=1e-4"
