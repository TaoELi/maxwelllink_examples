#!/bin/bash

# nmol should be the square of an int 
for nmol in 256
do
    # copy template files
    cp run_only_mol.sh run_only_mol_nmol${nmol}.sh
    cp run_both.sh run_both_nmol${nmol}.sh
    # modify nmol in the copied scripts
    sed -i "s/nmol=1;/nmol=${nmol};/" run_only_mol_nmol${nmol}.sh
    sed -i "s/nmol=1;/nmol=${nmol};/" run_both_nmol${nmol}.sh
    # submit the jobs
    # bash run_only_mol_nmol${nmol}.sh
    bash run_both_nmol${nmol}.sh
done
