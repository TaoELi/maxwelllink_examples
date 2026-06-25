ngrid_list=(144 288 576 1152)
amp=0.007
up=36
freq=2438.84

CURRENTFOLDER=$(pwd)
for nidx in {0..3}; do
    cp summit_all.sh summit_all_"$nidx".sh
    sed -i "s|bridge_prefix=\"rtx4090_\"|bridge_prefix=\"rtx4060_${nidx}\"|" "summit_all_${nidx}.sh"
    sh summit_all_"$nidx".sh "$up" "$freq" "$amp" "${ngrid_list[$nidx]}"
    sleep 1s
    rm summit_all_"$nidx".sh
done
