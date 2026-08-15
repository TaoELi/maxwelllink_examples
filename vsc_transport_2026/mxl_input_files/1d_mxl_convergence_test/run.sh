ngrid=144
amp_list=(0.001 0.0016 0.003 0.005 0.007)
list=(24 27 30 36 42 48 54 60 66)
freq=(2423.83 2427.16 2430.50 2438.84 2450.52 2462.20 2477.21 2495.56 2515.58)
CURRENTFOLDER=$(pwd)
for aidx in {0..4}; do
    for nframe in {0..9}; do
        traj=0
        cp summit_all.sh summit_all_"$aidx"_"$nframe"_"$ngrid".sh
        sed -i "s|bridge_prefix=\"rtx4090_\"|bridge_prefix=\"rtx1080_${aidx}_${nframe}_\"|" "summit_all_${aidx}_${nframe}_${ngrid}.sh"
        sh summit_all_"$aidx"_"$nframe"_"$ngrid".sh "${list[$traj]}" "${freq[$traj]}" "${amp_list[$aidx]}" "$ngrid" "$nframe"
        sleep 1s
        rm summit_all_"$aidx"_"$nframe"_"$ngrid".sh
    done
done
