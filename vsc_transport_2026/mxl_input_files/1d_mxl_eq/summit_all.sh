#!/bin/bash
# the driver code should be submitted after the main code is submitted and running
for nframe in {0..9}; do
    cp -r files zrun_eq_"$nframe"/
    cd zrun_eq_"$nframe"

    echo "Changed to directory: zrun_eq_$nframe"
    new_seed=$((114514 + nframe))
    sed -i "s/seed=114514/seed=${new_seed}/" multmodes_vsc_hpc.py
    sed -i "s/sample_idx=-1/sample_idx=${nframe}/" launch_lmp_xml.sh

    job_main_id=$(sbatch submit_main.sh | awk '{print $4}')
    echo "Submitted main job with Job ID: ${job_main_id}"

    sbatch --dependency=after:${job_main_id} submit_driver.sh
    cd ..
done