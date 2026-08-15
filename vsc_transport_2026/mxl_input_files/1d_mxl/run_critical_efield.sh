amp_list=(0.0016 0.0024)
list=(24 36)
freq=(2423.83 2438.84)
CURRENTFOLDER=$(pwd)
for aidx in {0..0}; do
    for nframe in {0..0}; do
        cp summit_all.sh summit_all_"$aidx"_"$aidx"_"$nframe".sh
        sed -i "s|bridge_prefix=\"rtx4090_\"|bridge_prefix=\"rtx1080_${aidx}_${aidx}_${nframe}\"|" "summit_all_${aidx}_${aidx}_${nframe}.sh"
        sh summit_all_"$aidx"_"$aidx"_"$nframe".sh "${list[$aidx]}" "${freq[$aidx]}" "${amp_list[$aidx]}" "$nframe"
        sleep 1s
        rm summit_all_"$aidx"_"$aidx"_"$nframe".sh
    done
done
