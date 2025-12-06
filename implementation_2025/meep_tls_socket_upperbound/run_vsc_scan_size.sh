#!/bin/bash

for size in 1 16 128 256 512 1024 2048 4096 8192 16384 32768 65536
do

    DIR="ntls_${size}"
    mkdir -p ${DIR}
    cd ${DIR}

    # copy input files from template
    cp ../template/* .

    # modify input parameters to accomodate different system sizes
    # 1. modify python file to adjust ntls accordingly
    sed -i -e "s/n_molecule=2/n_molecule=${size}/g" meep_tls_spontaneous_emission.py
    # for driver, the maximal molecular number is set to 2048: if size > 2048, we will launch multiple drivers
    size_driver=$size
    if (( size_driver > 2048 )); then
        size_driver=2048
    fi
    ratio_driver=$(( size / size_driver ))
    echo "Total size: ${size}, size per driver: ${size_driver}, number of drivers launched: ${ratio_driver}"

    sed -i -e "s/n_molecule=2/n_molecule=${size_driver}/g" submit_driver.sh
    sed -i -e "s/-J tls_se/-J sed_${size}/g" submit_driver.sh
    sed -i -e "s/-J tls_se/-J sem_${size}/g" submit_main.sh

    if (( size > 2048 )); then
        echo "changing copy_driver from 1 to ${ratio_driver} in submit_all.sh"
        sed -i -e "s/copy_driver=1/copy_driver=${ratio_driver}/" submit_all.sh
    fi

    # 4. submit jobs
    ./submit_all.sh

    cd ..

done
