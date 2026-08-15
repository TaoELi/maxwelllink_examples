ngrid=144
#amp_list=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8)
#amp_list=(0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08)
#list=(48)
#freq=(290.26)
list=(96)
freq=(925.82)
amp_list=(0.001 0.007)
CURRENTFOLDER=$(pwd)
for aidx in {0..0}; do
    for nframe in {3,4,5,7}; do
        traj=0
        cp summit_all.sh summit_all_"$aidx"_"$traj"_"$nframe".sh
        sed -i "s|bridge_prefix=\"rtx4090_\"|bridge_prefix=\"rtx1080_${aidx}_${traj}_\"|" "summit_all_${aidx}_${traj}_${nframe}.sh"
        sh summit_all_"$aidx"_"$traj"_"$nframe".sh "${list[$traj]}" "${freq[$traj]}" "${amp_list[$aidx]}" "$ngrid" "$nframe"
        sleep 1s
        rm summit_all_"$aidx"_"$traj"_"$nframe".sh
    done
done
