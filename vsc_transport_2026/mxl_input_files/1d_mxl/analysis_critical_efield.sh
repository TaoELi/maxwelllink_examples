if [ ! -d data ] ; then
    mkdir data
fi

amp_list=(0.0016 0.0024)
up_list=(24 36)
CURRENTFOLDER=$(pwd)
for aidx in {0..1}; do
    for nframe in {0..9}; do
        mv $CURRENTFOLDER/zrun_"${amp_list[$aidx]}"_"${up_list[$aidx]}"_"$nframe"/multimode_cavmd_4t4_neq.h5 $CURRENTFOLDER/data/multimode_cavmd_"$nframe"_"${up_list[$aidx]}"_"${amp_list[$aidx]}"_neq.h5
    done
done

if [ ! -d scpdata ] ; then
    mkdir scpdata
fi

python get_average_data_critical_efield.py
