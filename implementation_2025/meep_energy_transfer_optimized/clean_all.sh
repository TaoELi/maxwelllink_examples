#!/bin/bash

for dir in energy_transfer_*
do

echo "Running in directory: $dir"
cd $dir

./clean.sh
rm psi4* *.npz

cd ..

done
