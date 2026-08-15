#!/bin/bash
# following 3 commands required by OMP(OpenMP) only within one node
#SBATCH -A che250091
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=mxl_meso
#SBATCH --partition=shared
#SBATCH --time=1-00:00:00
#SBATCH --output=./nve_%A_%a.out

if [ ! -d data ] ; then
    mkdir data
fi

amp_list=(0.001 0.003 0.005 0.007)
up_list=(24 27 30 36 42 48 54 60 66)
CURRENTFOLDER=$(pwd)
for aidx in {0..3}; do
    for traj in {0..8}; do
        for nframe in {0..9}; do
        python coenergy.py $CURRENTFOLDER/zrun_"${amp_list[$aidx]}"_"${up_list[$traj]}"_"$nframe"/lmp_input
        #mv $CURRENTFOLDER/zrun_"${amp_list[$aidx]}"_"${up_list[$traj]}"_"$nframe"/multimode_cavmd_4t4_neq.h5 $CURRENTFOLDER/data/multimode_cavmd_"$nframe"_"${up_list[$traj]}"_"${amp_list[$aidx]}"_neq.h5
        mv $CURRENTFOLDER/zrun_"${amp_list[$aidx]}"_"${up_list[$traj]}"_"$nframe"/lmp_input/coenergy.npy $CURRENTFOLDER/data/coenergy_"$nframe"_"${up_list[$traj]}"_"${amp_list[$aidx]}".npy
        done
    done
done

if [ ! -d scpdata ] ; then
    mkdir scpdata
fi

python get_average_data.py
