if [ ! -d data ] ; then
    mkdir data
fi

amp_list=(5000 2000 1000 800 500 100)
CURRENTFOLDER=$(pwd)
for aidx in {0..0}; do
    for nframe in {0..9}; do
        mv $CURRENTFOLDER/zrun_"${amp_list[$aidx]}"_"$nframe"/multimode_cavmd_4t4_neq.h5 $CURRENTFOLDER/data/multimode_cavmd_"$nframe"_"${amp_list[$aidx]}"_neq.h5
    done
done

if [ ! -d scpdata ] ; then
    mkdir scpdata
fi

python get_average_data.py
