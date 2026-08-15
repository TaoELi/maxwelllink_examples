amp_list=(0.001 0.003 0.005 0.007)
list=(24 27 30 36 42 48 54 60 66)
freq=(2423.83 2427.16 2430.50 2438.84 2450.52 2462.20 2477.21 2495.56 2515.58)
CURRENTFOLDER=$(pwd)
for aidx in {2..2}; do
    for traj in {3..3}; do
        for nframe in {4..4}; do
            cp summit_all.sh summit_all_"$aidx"_"$traj"_"$nframe".sh
            sed -i "s|bridge_prefix=\"rtx4090_\"|bridge_prefix=\"rtx1080_${aidx}_${traj}_${nframe}\"|" "summit_all_${aidx}_${traj}_${nframe}.sh"
            sh summit_all_"$aidx"_"$traj"_"$nframe".sh "${list[$traj]}" "${freq[$traj]}" "${amp_list[$aidx]}" "$nframe"
            sleep 0.1s
            rm summit_all_"$aidx"_"$traj"_"$nframe".sh
        done
    done
done
