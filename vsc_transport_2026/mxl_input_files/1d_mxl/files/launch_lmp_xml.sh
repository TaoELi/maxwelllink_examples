#!/bin/bash

cd lmp_input/

# This script prepares LAMMPS input files for MaxwellLink tutorials.
HOST=$1
PORT=$2
BOX_ID=$3
echo "Preparing LAMMPS input files with host $HOST, port $PORT and box ID $BOX_ID..."

# modify in_mxl.lmp to set the correct port number
cp in_mxl.lmp in_$BOX_ID.lmp
sed -i -e "s/HOST/$HOST/g" in_$BOX_ID.lmp
sed -i -e "s/PORT/$PORT/g" in_$BOX_ID.lmp
sed -i -e "s/BOX_ID/$BOX_ID/g" in_$BOX_ID.lmp
randomseed=$((100 * BOX_ID))
sed -i -e "s/velocity all create 300.0 2345187/velocity all create 300.0 $randomseed/g" in_$BOX_ID.lmp

# check whether lmp_mxl command is available
if ! command -v lmp_mxl &> /dev/null
then
    echo "lmp_mxl could not be found. Please ensure LAMMPS with MaxwellLink is installed and lmp_mxl is in your PATH."
    exit 1
fi

# launch lammps from bash script
#lmp_mxl < in_$BOX_ID.lmp
lmp_mxl -in in_$BOX_ID.lmp -screen out_$BOX_ID.lmp -log none
rm in_$BOX_ID.lmp
#rm out_$BOX_ID.lmp
cd ..
