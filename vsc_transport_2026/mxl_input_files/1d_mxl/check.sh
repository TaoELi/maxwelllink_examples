amp_list=(0.0016 0.0024)
list=(24 36)
freq=(2423.83 2438.84)
CURRENTFOLDER=$(pwd)
for aidx in {0..1}; do
    for nframe in {0..9}; do
    cd zrun_"${amp_list[$aidx]}"_"${list[$aidx]}"_"$nframe"
    tail -n 4 nve_*.out
    cd ..
    done
done
