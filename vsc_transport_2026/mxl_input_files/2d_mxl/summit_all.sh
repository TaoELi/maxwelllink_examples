#!/bin/bash
traj=36_0.002
# change bridge prefix to avoid conflicting
bridge_prefix="rtx5080_"
rm -rf /tmp/"$bridge_prefix"*
# the driver code should be submitted after the main code is submitted and running
cp -r files zrun_"$traj"/
cd zrun_"$traj"
sed -i 's/unix_prefix="bridge_"/unix_prefix="'"$bridge_prefix"'"/' ./multmodes_vsc_hpc.py
sed -i 's/HOST="bridge_${traj}"/HOST="'"$bridge_prefix"'${traj}"/' ./submit_driver.sh

job_main_id=$(sbatch submit_main.sh | awk '{print $4}')
echo "Submitted main job with Job ID: ${job_main_id}"

sbatch --dependency=after:${job_main_id} submit_driver.sh
