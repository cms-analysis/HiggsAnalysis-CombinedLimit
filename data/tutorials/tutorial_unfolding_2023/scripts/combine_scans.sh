
for p in r_zh_75_150 r_zh_150_250noj r_zh_150_250wj r_zh_250_400 r_zh_gt400
do
    hadd -k -f scan_blinded_${p}.root higgsCombine.scans_blinded.${p}.POINTS.*.MultiDimFit.mH120.root
done
