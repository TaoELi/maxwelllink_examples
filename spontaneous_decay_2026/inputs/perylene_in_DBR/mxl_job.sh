#!/bin/bash

#for ((ii=40; ii <= 64; ii++))
for ((ii=41; ii <= 64; ii++))
#for ii in $heights;
do
    frequency=$(echo "scale=2; $ii/20" | bc)
    if [ ! -d $frequency ]; then
        mkdir $frequency
    #else
       #rm -r $frequency/*
    fi

    cp -r fit_energy_decay_rate.py dftbplus_input $frequency
    sed -e "s/c-2.634eV/$frequency/g" job_fdtd_dftb.sh > $frequency/job_fdtd_dftb.sh
    sed -e "s/CAVITY_FREQUENCY/$frequency/g" fdtd.py > $frequency/fdtd.py
    #sed -e "s/adi-1/adi-$ii/g" job_cp2k.sh > $ii/job_cp2k.sh
    cd $frequency
    sbatch job_fdtd_dftb.sh
    cd ../
    sleep 0.25
done
