cavityloss=(5000 2000 1000 800 500 100)
photoncharge=(0.066 0.105 0.151 0.171 0.210 0.480)
CURRENTFOLDER=$(pwd)
for aidx in {0..0}; do
    for nframe in {0..9}; do
        cp summit_all.sh summit_all_"$aidx"_"$aidx"_"$nframe".sh
        sed -i "s|bridge_prefix=\"rtx4090_\"|bridge_prefix=\"rtx1080_${aidx}_${aidx}_${nframe}\"|" "summit_all_${aidx}_${aidx}_${nframe}.sh"
        sh summit_all_"$aidx"_"$aidx"_"$nframe".sh 36 2438.84 0.007 "$nframe" "${cavityloss[$aidx]}" "${photoncharge[$aidx]}"
        sleep 1s
        rm summit_all_"$aidx"_"$aidx"_"$nframe".sh
    done
done
