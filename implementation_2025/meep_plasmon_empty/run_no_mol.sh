#!/bin/bash

# submit the job without molecule to generate spectrum
dir="no_mol_with_dielectric"
mkdir -p ${dir}
cd ${dir}

cp ../template/* .

cp ../vac/emitter-refl-flux.h5 .

sed -i "s/plasmon_nomol/plasmon_nomol/" submit_main.sh

sbatch submit_main.sh

cd ..
