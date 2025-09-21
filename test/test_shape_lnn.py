#!/usr/bin/env python3
# 1ch, 1sig,1bkg; 2 shape systs; sys1 uses 'shape?' with asym lnN on sig and shape-only on bkg; autoMCStats on

import numpy as np, ROOT
ROOT.gROOT.SetBatch(); np.random.seed(123)

# config
nb, lo, hi = 40, 0.0, 10.0
Nsig, Nbkg = 600, 5000
mu, sig, lam = 4.0, 0.6, 0.45
outroot, card = "shapes.root", "datacard.txt"

# helpers
def H(name): h=ROOT.TH1F(name,name,nb,lo,hi); h.Sumw2(); return h
def fill(h, a): [h.Fill(x) for x in a]
def tnorm(n,m,s,lo,hi):
    out=[]
    while len(out)<n:
        x=np.random.normal(m,s,size=n)
        x=x[(x>=lo)&(x<hi)]; out.extend(x.tolist())
    return np.array(out[:n])
def texp(n,l,lo,hi):
    L=hi-lo; u=np.random.uniform(0,1,size=n); return lo - np.log(1-u*(1-np.exp(-l*L)))/l

# nominal
hs=H("sig"); fill(hs, tnorm(Nsig, mu, sig, lo, hi))
hb=H("bkg"); fill(hb, texp(Nbkg, lam, lo, hi))
hd=H("data_obs")
for i in range(1,nb+1):
    m=hs.GetBinContent(i)+hb.GetBinContent(i)
    hd.SetBinContent(i, np.random.poisson(m)); hd.SetBinError(i,0.0)

# define two shape systematics (both have Up/Down for both processes)
# sys1: small mean/sigma tweaks on sig, lambda tweak on bkg
hs_sys1Up   = H("sig_sys1Up");   fill(hs_sys1Up,   tnorm(Nsig, mu+0.20, sig*1.00, lo, hi))
hs_sys1Down = H("sig_sys1Down"); fill(hs_sys1Down, tnorm(Nsig, mu-0.20, sig/1.00, lo, hi))
hb_sys1Up   = H("bkg_sys1Up");   fill(hb_sys1Up,   texp(Nbkg, lam*1.10, lo, hi))
hb_sys1Down = H("bkg_sys1Down"); fill(hb_sys1Down, texp(Nbkg, lam/1.10, lo, hi))

# sys2: different tweaks
hs_sys2Up   = H("sig_sys2Up");   fill(hs_sys2Up,   tnorm(Nsig, mu-0.15, sig*1.20, lo, hi))
hs_sys2Down = H("sig_sys2Down"); fill(hs_sys2Down, tnorm(Nsig, mu+0.15, sig/1.20, lo, hi))
hb_sys2Up   = H("bkg_sys2Up");   fill(hb_sys2Up,   texp(Nbkg, lam*0.85, lo, hi))
hb_sys2Down = H("bkg_sys2Down"); fill(hb_sys2Down, texp(Nbkg, lam/0.85, lo, hi))

# write shapes
f=ROOT.TFile(outroot,"RECREATE")
for h in [hd,hs,hb,hs_sys1Up,hs_sys1Down,hb_sys1Up,hb_sys1Down,hs_sys2Up,hs_sys2Down,hb_sys2Up,hb_sys2Down]: h.Write()
f.Close()

# datacard:
# - sys1 uses 'shape?' with asymmetric lnN factor on sig (1.07/0.93) and shape-only (1) on bkg
# - sys2 uses standard 'shape' for both
dc=[]
dc+=["imax 1","jmax 1","kmax *",""]
dc+=[f"shapes * ch1 {outroot} $PROCESS $PROCESS_$SYSTEMATIC",
     f"shapes data_obs ch1 {outroot} data_obs",""]
dc+=["bin ch1","observation -1",""]
dc+=["bin ch1 ch1","process sig bkg","process 0 1","rate -1 -1",""]
dc+=["sys1 shape? 1.07/0.93 1","sys2 shape 1 1",""]  # <- key test line for the bugfix
dc+=["* autoMCStats 10",""]                           # per-bin MC stat from Sumw2
open(card,"w").write("\n".join(dc))

print(f"Wrote {outroot} and {card}")
