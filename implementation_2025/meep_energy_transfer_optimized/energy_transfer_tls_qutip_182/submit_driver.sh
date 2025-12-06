#!/bin/bash
# this job should be run after submit_main.sh and after the tcp host and port info file is created in file

#SBATCH --nodes=1         # Total # of nodes (must be 1 for serial job)
#SBATCH --ntasks=2        # Total # of MPI tasks (should be 1 for serial job)
#SBATCH --time=2-00:00:00    # Total run time limit (hh:mm:ss)
#SBATCH -J etd_qutip      # Job name
#SBATCH -o driver.o%j      # Name of stdout output file
#SBATCH -e driver.e%j      # Name of stderr error file
#SBATCH -p shared  # Queue (partition) name

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

# launch tls donor driver in the background
mxl_driver --model tls --address $HOST --port $PORT --param "omega=0.491865, mu12=170.3072, orientation=2, pe_initial=0.5" &
sleep 10s

# launch acceptor driver in the foreground
mxl_driver --model qutip --address $HOST --port $PORT --param "preset=custom, module=hcn_qutip.py"  --verbose
