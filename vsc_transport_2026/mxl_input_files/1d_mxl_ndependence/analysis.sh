if [ ! -d data ] ; then
    mkdir data
fi
ngrid_list=(144 288 576 1152)
amp=0.007
up=36
CURRENTFOLDER=$(pwd)
for nidx in {0..3}; do
    for nframe in {0..9}; do
        mv $CURRENTFOLDER/zrun_"$amp"_"$up"_"${ngrid_list[$nidx]}"_"$nframe"/multimode_cavmd_4t4_neq.h5 $CURRENTFOLDER/data/multimode_cavmd_"${ngrid_list[$nidx]}"_"$nframe"_"$up"_"$amp"_neq.h5
    done
done

if [ ! -d scpdata ] ; then
    mkdir scpdata
fi

python get_average_data.py