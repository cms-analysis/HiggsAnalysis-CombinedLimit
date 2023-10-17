
for p in r_zh_75_150 r_zh_150_250noj r_zh_150_250wj r_zh_250_400 r_zh_gt400
do
    python scripts/plot1DScan.py scan_blinded_${p}.root -o scan_plot_${p}_blinded --POI ${p} --json summary_zh_stxs_blinded.json  
done
