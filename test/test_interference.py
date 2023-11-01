#!/usr/bin/env python3
import ROOT
import numpy as np
import subprocess
import json


def array2vector_2d(array):
    assert len(array.shape) == 2
    out = ROOT.std.vector["std::vector<double>"]()
    out.reserve(len(array))
    for row in array:
        out.push_back(row)
    return out


def to_TH1(name, array):
    h = ROOT.TH1D(name, name, 20, 250, 1200)
    for i, val in enumerate(array):
        h.SetBinContent(i + 1, val)
    return h


# make some shapes
fout = ROOT.TFile("shapes.root", "recreate")
nom = 100 * np.array(
    [
        0.14253603,
        0.21781641,
        0.22698837,
        0.19603483,
        0.19259561,
        0.15552859,
        0.13909682,
        0.09438712,
        0.08521593,
        0.06878416,
        0.06419854,
        0.04318116,
        0.04776676,
        0.03057073,
        0.02866007,
        0.02292805,
        0.02139951,
        0.02063524,
        0.01222829,
        0.01337466,
    ]
)
obs = nom.copy()

hnom = to_TH1("VBFHH", nom)
for i, val in enumerate(nom):
    hnom.SetBinError(i + 1, np.sqrt(val) * 0.1)
hnom.Write()
to_TH1("VBFHH_jesUp", nom * np.linspace(0.95, 1.05, 20)).Write()
to_TH1("VBFHH_jesDown", nom / np.linspace(0.95, 1.05, 20)).Write()

nom = 100 * np.exp(-np.linspace(0, 1, 20))
obs += nom
hnom = to_TH1("background", nom)
for i, val in enumerate(nom):
    hnom.SetBinError(i + 1, np.sqrt(val) * 0.1)
hnom.Write()
to_TH1("background_jesUp", nom * np.linspace(0.95, 1.05, 20)).Write()
to_TH1("background_jesDown", nom / np.linspace(0.95, 1.05, 20)).Write()

to_TH1("data_obs", np.round(obs)).Write()

# write a card
with open("card.txt", "w") as fout:
    fout.write(
        """\
Combination of card.txt
imax 1 number of bins
jmax 1 number of processes minus 1
kmax 3 number of nuisance parameters
----------------------------------------------------------------------------------------------------------------------------------
shapes *    ch1  shapes.root $PROCESS $PROCESS_$SYSTEMATIC
----------------------------------------------------------------------------------------------------------------------------------
bin          ch1
observation  -1
----------------------------------------------------------------------------------------------------------------------------------
bin                               ch1         ch1
process                           VBFHH       background
process                           0           1
rate                              -1          -1
----------------------------------------------------------------------------------------------------------------------------------
bgnorm                  lnN       -           1.3
jes                     shape     1.0         1.0
lumi                    lnN       1.02        1.02
* autoMCStats 0
"""
    )

# write the scaling data
scaling = [
    {
        "channel": "ch1",
        "process": "VBFHH",
        "parameters": [
            "expr::a0('@0*@1', kv[1,0,2], kl[1,0,2])",
            "expr::a1('@0*@0', kv[1,0,2])",
            "k2v[1,0,2]",
        ],
        "scaling": [
            [3.30353674666415, -8.54170982038222, 22.96464188467882, 4.2353483207128, -11.07996258835088, 5.504469544697623],
            [2.20644332142891, -7.076836641962523, 23.50989689214267, 4.053185685866683, -13.08569222837996, 7.502346155380032],
            [2.323314512827915, -9.040565356058327, 35.97836755817654, 6.135410429832603, -23.87500375686657, 16.25863529518014],
            [0.5925805332888091, -3.139640204484572, 17.13589230378096, 1.976422909655008, -10.515009347222, 6.627980447033357],
            [0.3505170703269003, -2.237551781735236, 14.84170437843385, 1.483246407331564, -9.493221195626012, 6.302831691298613],
            [0.2334740816002986, -1.79981530542472, 14.40275907010786, 1.22573469124098, -9.473328823485497, 6.458585723630323],
            [0.1959374052725543, -1.757624939190617, 16.26453309470772, 1.257071800464856, -11.27433100208976, 8.089297781650776],
            [0.186504029519429, -1.822069966644771, 18.49884116334009, 1.391201452068932, -13.60944960795888, 10.39529105220993],
            [0.0972664880869715, -1.169327097322687, 14.59008882076805, 0.863528399754145, -10.3792682464399, 7.682778579161866],
            [0.0878055201507951, -1.156467907936071, 15.89008007570035, 0.9032722960981284, -11.89269332401088, 9.31389227584649],
            [0.2006114827320551, -2.776793529232688, 39.12099794229794, 2.335437756631023, -32.32054661687562, 27.20219535392458],
            [0.05179066270217598, -0.8550833654170061, 14.82948643635128, 0.6753041305320768, -11.17147343869588, 8.821228248108168],
            [0.05095452842967559, -0.948210749952708, 18.25692153880631, 0.7561787601252029, -14.08060600559859, 11.23739992361621],
            [0.01801729907958822, -0.4102710173962256, 9.828612567782274, 0.2950547079487041, -6.724075992501888, 4.831954737036956],
            [0.02762233787081839, -0.6400852678490596, 15.46123935330864, 0.5127706114146487, -11.88156027745066, 9.528888176590677],
            [0.031654990010153, -0.7790988142691099, 19.86193598270219, 0.648607157613526, -15.97020317079789, 13.30779868219462],
            [0.01316620750610005, -0.3821583465978776, 11.67177892863827, 0.2914985575363581, -8.480579644763935, 6.457533731506537],
            [0.02886802887767344, -0.8369930524359994, 24.91187028333814, 0.7103796160970196, -20.57801184441048, 17.4685122492831],
            [0.02930989281275648, -0.9242683589392606, 29.81526833971985, 0.7918145320034853, -25.01001674094063, 21.4403629032202],
            [0.0516036089235802, -1.673006024492507, 55.06880032083917, 1.526329404862879, -49.49344845658728, 45.15984622267106],
        ],
    },
]

with open("scaling.json", "w") as fout:
    json.dump(scaling, fout)

t2wcmd = [
    "text2workspace.py",
    "card.txt",
    "-P",
    "HiggsAnalysis.CombinedLimit.InterferenceModels:interferenceModel",
    "--PO",
    "verbose",
    "--PO",
    "scalingData=scaling.json",
    "--PO",
    "POIs=kl[1,0,2]:kv[1,0,2]:k2v[1,0,2]",
]

ret = subprocess.call(t2wcmd)
assert ret == 0

histsum_args = ["--for-fits", "--no-wrappers", "--use-histsum", "-o", "card_histsum.root"]
ret = subprocess.call(t2wcmd + histsum_args)
assert ret == 0

fws = ROOT.TFile.Open("card.root")
w = fws.Get("w")


def setvars(x, kl, kv, k2v):
    w.var("CMS_th1x").setVal(x)
    w.var("kl").setVal(kl)
    w.var("kv").setVal(kv)
    w.var("k2v").setVal(k2v)


func = w.function("shapeSig_ch1_VBFHH_morph_externalMorph")
assert func

setvars(0, 1, 1, 1)
assert abs(func.getVal() - 1.0) < 1e-14, func.getVal()
setvars(1, 1.1, 1, 1)
assert abs(func.getVal() - 0.8586229062809139) < 1e-14, func.getVal()
setvars(2, 1.1, 0.9, 1.3)
assert abs(func.getVal() - 4.372110974178483) < 1e-14, func.getVal()

# toy generation is different between the histsum and histfunc models, somehow
ntoys = 10
ret = subprocess.call("combine -M GenerateOnly card.root -t {ntoys} --saveToys".format(ntoys=ntoys).split(" "))
assert ret == 0

ret = subprocess.call(
    "combine -M MultiDimFit card.root -t {ntoys} --toysFile higgsCombineTest.GenerateOnly.mH120.123456.root -n HistFunc".format(ntoys=ntoys).split(" ")
)
assert ret == 0

ret = subprocess.call(
    "combine -M MultiDimFit card_histsum.root -t {ntoys} --toysFile higgsCombineTest.GenerateOnly.mH120.123456.root -n HistSum".format(ntoys=ntoys).split(" ")
)
assert ret == 0

f_histfunc = ROOT.TFile.Open("higgsCombineHistFunc.MultiDimFit.mH120.123456.root")
f_histsum = ROOT.TFile.Open("higgsCombineHistSum.MultiDimFit.mH120.123456.root")

ndiff = {"kl": 0, "kv": 0, "k2v": 0}
for row1, row2 in zip(f_histfunc.Get("limit"), f_histsum.Get("limit")):
    if abs(row1.kl - row2.kl) > 1e-4:
        ndiff["kl"] += 1
    if abs(row1.kv - row2.kv) > 1e-4:
        ndiff["kv"] += 1
    if abs(row1.k2v - row2.k2v) > 1e-4:
        ndiff["k2v"] += 1

print("Out of {ntoys} toys, {ndiff} are not matching (tolerance: 1e-4) between CMSHistFunc and CMSHistSum".format(ntoys=ntoys, ndiff=ndiff))
