#!/bin/bash

# submit the job without molecule to generate spectrum
dir="vac"
mkdir -p ${dir}
cd ${dir}

cp ../template/* .

sed -i "s/plasmon_nomol/plasmon_vac/" submit_vac.sh

sbatch submit_vac.sh

cd ..
