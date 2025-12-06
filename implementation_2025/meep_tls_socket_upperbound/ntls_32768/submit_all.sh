#!/bin/bash

# the driver code should be submitted after the main code is submitted and running

job_main_id=$(sbatch submit_main.sh | awk '{print $4}')
echo "Submitted main job with Job ID: ${job_main_id}"

copy_driver=16
# submit multiple drivers if copy_driver > 1
for (( i=1; i<=copy_driver; i++ ))
do
    sbatch --dependency=after:${job_main_id} submit_driver.sh
done