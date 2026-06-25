#!/bin/bash
traj="$1"
up_freq="$2"
amp="$3"
ngrid="$4"
# change bridge prefix to avoid conflicting
bridge_prefix="rtx4090_"
rm -rf /tmp/"$bridge_prefix"*
# the driver code should be submitted after the main code is submitted and running
cp -r files zrun_"$amp"_"$traj"_"$ngrid"/
cd zrun_"$amp"_"$traj"_"$ngrid"
sed -i 's/unix_prefix="bridge_"/unix_prefix="'"$bridge_prefix"'"/' ./multmodes_vsc_hpc.py
sed -i 's/HOST="bridge_${traj}"/HOST="'"$bridge_prefix"'${traj}"/' ./submit_driver.sh
sed -i "s/omega_au=2413.82/omega_au=$up_freq/" ./multmodes_vsc_hpc.py
sed -i "s/k_parallel_au=1/k_parallel_au=$traj/" ./multmodes_vsc_hpc.py
sed -i "s/amplitude_au=0.002/amplitude_au=$amp/" ./multmodes_vsc_hpc.py
sed -i "s/N_grid = 144/N_grid = $ngrid/" ./multmodes_vsc_hpc.py
sed -i "s|#SBATCH --array=0-3|#SBATCH --array=0-$((ngrid / 36 -1))|" submit_driver.sh
job_main_id=$(sbatch submit_main.sh | awk '{print $4}')
echo "Submitted main job with Job ID: ${job_main_id}"

sbatch --dependency=after:${job_main_id} submit_driver.sh
