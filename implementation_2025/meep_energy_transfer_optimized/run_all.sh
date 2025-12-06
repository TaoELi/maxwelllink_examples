#!/bin/bash

for dir in energy_transfer_*
do

echo "Running in directory: $dir"
cd $dir

./submit_all.sh

cd ..

done