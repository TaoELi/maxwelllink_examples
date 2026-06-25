if [ ! -d data ] ; then
    mkdir data
fi
ngrid=144
amp_list=(0.001 0.003 0.005 0.007)
up_list=(24 27 30 36 42 48 54 60 66)
CURRENTFOLDER=$(pwd)
for aidx in {0,3}; do
    for traj in {0..8}; do
        python coenergy.py $CURRENTFOLDER/zrun_"${amp_list[$aidx]}"_"${up_list[$traj]}"_"$ngrid"/lmp_input
        mv $CURRENTFOLDER/zrun_"${amp_list[$aidx]}"_"${up_list[$traj]}"_"$ngrid"/multimode_cavmd_4t4_neq.h5 $CURRENTFOLDER/data/multimode_cavmd_"$ngrid"_"${up_list[$traj]}"_"${amp_list[$aidx]}"_neq.h5
        mv $CURRENTFOLDER/zrun_"${amp_list[$aidx]}"_"${up_list[$traj]}"_"$ngrid"/lmp_input/coenergy.npy $CURRENTFOLDER/data/coenergy_"$ngrid"_"${up_list[$traj]}"_"${amp_list[$aidx]}".npy
        rm -r $CURRENTFOLDER/zrun_"${amp_list[$aidx]}"_"${up_list[$traj]}"_"$ngrid"
    done
done
