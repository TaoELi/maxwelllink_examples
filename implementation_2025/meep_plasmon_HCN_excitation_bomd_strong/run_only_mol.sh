#!/bin/bash

# submit the job without dielectric medium and with molecule
nmol=1;
dir="nmol_${nmol}_no_dielectric"
mkdir -p ${dir}
cd ${dir}

cp ../template/* .

sed -i 's/-dielectric/-mol/g' submit_main.sh
sed -i "s/plasmon_nomol/mol_${nmol}/" submit_main.sh
sed -i "s/nmol=1/nmol=${nmol}/" submit_main.sh
# if nmol <= 128, we can set ntasks = nmol; else ntasks = 128 and we will submit multiple driver jobs later
if [ ${nmol} -le 128 ]; then
    sed -i "s/--ntasks=128/--ntasks=${nmol}/" submit_driver.sh
    sed -i "s/nmol=128/nmol=${nmol}/" submit_driver.sh
fi
sed -i "s/nmol=128;/nmol=${nmol};/" submit_all.sh


./submit_all.sh

cd ..
