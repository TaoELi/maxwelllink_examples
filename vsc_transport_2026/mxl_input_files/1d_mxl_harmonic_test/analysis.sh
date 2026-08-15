#!/bin/bash
# following 3 commands required by OMP(OpenMP) only within one node
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=analysis
#SBATCH --partition=standard
#SBATCH --time=7-00:00:00
#SBATCH --output=./nve_%A_%a.out

if [ ! -d data ] ; then
    mkdir data
fi
ngrid=144
amp_list=(0.001 0.007)
up_list=(96)
CURRENTFOLDER=$(pwd)
for aidx in {0..1}; do
    for nframe in {0..9}; do
        traj=0
        python coenergy.py $CURRENTFOLDER/zrun_"${amp_list[$aidx]}"_"${up_list[$traj]}"_"$ngrid"_"$nframe"/lmp_input
        mv $CURRENTFOLDER/zrun_"${amp_list[$aidx]}"_"${up_list[$traj]}"_"$ngrid"_${nframe}/multimode_cavmd_4t4_neq.h5 $CURRENTFOLDER/data/multimode_cavmd_"$nframe"_"${up_list[$traj]}"_"${amp_list[$aidx]}"_neq.h5
        mv $CURRENTFOLDER/zrun_"${amp_list[$aidx]}"_"${up_list[$traj]}"_"$ngrid"_${nframe}/lmp_input/coenergy.npy $CURRENTFOLDER/data/coenergy_"$nframe"_"${up_list[$traj]}"_"${amp_list[$aidx]}".npy
    done
done

if [ ! -d scpdata ] ; then
    mkdir scpdata
fi

python get_noneq_data.py
