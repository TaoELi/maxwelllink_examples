if [ ! -d data ] ; then
    mkdir data
fi

CURRENTFOLDER=$(pwd)
for nframe in {0..9}; do
    mv $CURRENTFOLDER/zrun_eq_"$nframe"/multimode_cavmd_4t4_eq.h5 $CURRENTFOLDER/data/multimode_cavmd_144_"$nframe"_eq.h5
    rm -r $CURRENTFOLDER/zrun_eq_"$nframe"
done

python get_final_data.py
