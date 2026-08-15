if [ ! -d data ] ; then
    mkdir data
fi
ngrid=144
amp_list=(0.001 0.0016 0.003 0.005 0.007)
up_list=(24 27 30 36 42 48 54 60 66)
CURRENTFOLDER=$(pwd)
for aidx in {0..4}; do
    for nframe in {0..9}; do
        traj=0
        mv $CURRENTFOLDER/zrun_"${amp_list[$aidx]}"_"${up_list[$traj]}"_"$ngrid"_"$nframe"/multimode_cavmd_4t4_neq.h5 $CURRENTFOLDER/data/multimode_cavmd_"$ngrid"_"${up_list[$traj]}"_"${amp_list[$aidx]}"_"$nframe"_neq.h5
        rm -r $CURRENTFOLDER/zrun_"${amp_list[$aidx]}"_"${up_list[$traj]}"_"$ngrid"_"$nframe"
    done
done

if [ ! -d scpdata ] ; then
    mkdir scpdata
fi

python get_final_data.py