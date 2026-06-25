if [ ! -d data ] ; then
    mkdir data
fi
ngrid_list=(144 288 576 1152)
amp=0.007
up=36
CURRENTFOLDER=$(pwd)
for nidx in {0..3}; do
    #python coenergy.py $CURRENTFOLDER/zrun_"$amp"_"$up"_"${ngrid_list[$nidx]}"/lmp_input
    mv $CURRENTFOLDER/zrun_"$amp"_"$up"_"${ngrid_list[$nidx]}"/multimode_cavmd_4t4_neq.h5 $CURRENTFOLDER/data/multimode_cavmd_"${ngrid_list[$nidx]}"_"$up"_"$amp"_neq.h5
    #mv $CURRENTFOLDER/zrun_"$amp"_"$up"_"${ngrid_list[$nidx]}"/lmp_input/coenergy.npy $CURRENTFOLDER/data/coenergy_"${ngrid_list[$nidx]}"_"$up"_"$amp".npy
    rm -r $CURRENTFOLDER/zrun_"$amp"_"$up"_"${ngrid_list[$nidx]}"
done
