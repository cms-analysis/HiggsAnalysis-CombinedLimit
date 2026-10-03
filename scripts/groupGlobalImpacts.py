#!/usr/bin/env python3
import math
import json
import argparse

def IsConstrained(param_info):
    return param_info["type"] != "Unconstrained"

parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", help="input json file")
args = parser.parse_args()

# Load the json output of combineTool.py -M Impacts
data = {}
with open(args.input) as jsonfile:
    data = json.load(jsonfile)

# # Set the global plotting style
# plot.ModTDRStyle(l=args.left_margin, b=0.10, width=(900 if args.checkboxes else 700), height=args.height)

# # We will assume the first POI is the one to plot
POIs = [ele["name"] for ele in data["POIs"]]
POIdata = {ele["name"]: ele for ele in data["POIs"]}
proto = {POI: 0. for POI in POIs}
fit_val = {}
err_tot_hi = {}
err_tot_lo = {}

# Get the list of groups
groups = set(["TOTAL"])
for ele in data["params"]:
    groups.update(ele["groups"])

err_hi = {grp : dict(proto) for grp in groups}
err_lo = {grp : dict(proto) for grp in groups}


for POI in POIs:
    fit_val[POI] = POIdata[POI]["fit"][1]
    err_tot_hi[POI] = POIdata[POI]["fit"][2] - POIdata[POI]["fit"][1]
    err_tot_lo[POI] = POIdata[POI]["fit"][1] - POIdata[POI]["fit"][0]


for ele in data["params"]:
    if IsConstrained(ele):
        name = ele["name"]
        # We expect to find "global_[POI]"
        for POI in POIs:
            impact = ele[f"global_{POI}"]
            impact_hi = impact[2] - impact[1] # impact from shifting global obs up
            impact_lo = impact[0] - impact[1] # impact from shifting global obs down
            # Check different situations
            contrib_hi = 0.
            contrib_lo = 0.
            if impact_hi >= 0. and impact_lo <= 0.:
                contrib_hi = impact_hi
                contrib_lo = impact_lo
            elif impact_hi <=0. and impact_lo >= 0.:
                contrib_lo = impact_hi
                contrib_hi = impact_lo
            elif impact_hi >= 0. and impact_lo >= 0.:
                print(f"Warning, parameter {name} has global impact on {POI} that are both positive: ({impact_hi},{impact_lo}), we will take the max of the two")
                contrib_hi = max(impact_hi, impact_lo)
            elif impact_lo <= 0. and impact_lo <= 0.:
                print(f"Warning, parameter {name} has global impact on {POI} that are both negative: ({impact_hi},{impact_lo}), we will take the min of the two")

                contrib_lo = min(impact_hi, impact_lo)
            for grp in ["TOTAL"] + ele["groups"]:
                err_hi[grp][POI] += pow(contrib_hi, 2)
                err_lo[grp][POI] += pow(contrib_lo, 2)
            # print(ele["name"], POI, impact_hi, impact_lo, contrib_hi, contrib_lo)

err_hi["STAT"] = dict(proto)
err_lo["STAT"] = dict(proto)
for POI in POIs:
    err_hi["STAT"][POI] = math.sqrt(pow(err_tot_hi[POI], 2) - err_hi["TOTAL"][POI])
    err_lo["STAT"][POI] = math.sqrt(pow(err_tot_lo[POI], 2) - err_lo["TOTAL"][POI])
for grp in groups:
    for POI in POIs:
        err_hi[grp][POI] = math.sqrt(err_hi[grp][POI])
        err_lo[grp][POI] = math.sqrt(err_lo[grp][POI])

for POI in POIs:
    print(f'POI: {POI:<20} {"Best-fit":>20} : {fit_val[POI]:<10.3f}')
    print(f'POI: {POI:<20} {"Total uncertainty":>20} : +{err_tot_hi[POI]:<10.3f} -{err_tot_lo[POI]:<10.3f}')
    print(f'POI: {POI:<20} {"Stat":>20} : +{err_hi["STAT"][POI]:<10.3f} -{err_lo["STAT"][POI]:<10.3f}')
    print(f'POI: {POI:<20} {"Syst":>20} : +{err_hi["TOTAL"][POI]:<10.3f} -{err_lo["TOTAL"][POI]:<10.3f}')
    for grp in groups:
        if grp in ["TOTAL", "STAT"]:
            continue
        print(f'POI: {POI:<20} {grp:>20} : +{err_hi[grp][POI]:<10.3f} -{err_lo[grp][POI]:<10.3f}')
 