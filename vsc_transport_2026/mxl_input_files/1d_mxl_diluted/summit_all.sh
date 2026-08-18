#!/bin/bash
traj="$1"
up_freq="$2"
amp="$3"
nframe="$4"
tau="$5"
qc="$6"

# change bridge prefix to avoid conflicting
bridge_prefix="rtx4090_"
rm -rf /tmp/"$bridge_prefix"*
# the driver code should be submitted after the main code is submitted and running
cp -r files zrun_"$tau"_"$nframe"/
cd zrun_"$tau"_"$nframe"
sed -i 's/unix_prefix="bridge_"/unix_prefix="'"$bridge_prefix"'"/' ./multmodes_vsc_hpc.py
sed -i 's/HOST="bridge_${traj}"/HOST="'"$bridge_prefix"'${traj}"/' ./submit_driver.sh
sed -i "s/omega_au=2413.82/omega_au=$up_freq/" ./multmodes_vsc_hpc.py
sed -i "s/k_parallel_au=1/k_parallel_au=$traj/" ./multmodes_vsc_hpc.py
sed -i "s/amplitude_au=0.002/amplitude_au=$amp/" ./multmodes_vsc_hpc.py
sed -i "s/tau_au=5000/tau_au=$tau/" ./multmodes_vsc_hpc.py
sed -i "s/photon_partial_charge=0.066/photon_partial_charge=$qc/" ./multmodes_vsc_hpc.py

new_seed=$((114514 + nframe))
sed -i "s/seed=114514/seed=${new_seed}/" multmodes_vsc_hpc.py
sed -i "s/sample_idx=-1/sample_idx=${nframe}/" launch_lmp_xml.sh

job_main_id=$(sbatch submit_main.sh | awk '{print $4}')
echo "Submitted main job with Job ID: ${job_main_id}"

sbatch --dependency=after:${job_main_id} submit_driver.sh
