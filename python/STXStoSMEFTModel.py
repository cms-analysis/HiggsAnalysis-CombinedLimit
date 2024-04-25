# Author: Jonathon Langford (ICL)
# Date: 10/2022
# Description: Model to describe how bins in STXS stage 1.2 scale using full set of dimension-6 EFT parameters
#              Equations calculated using nanoAOD reweighting using Madgraph reweighting modules
#              SMEFTsim

from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder
from math import exp
import ROOT
import os
import re
import sys
import json
import yaml
from collections import OrderedDict as od

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Ordered dicts crucial: e.g. choose WH_had before WH

MAP_HIGGS_DECAY_SMEFT = od()
MAP_HIGGS_DECAY_SMEFT["hgg"] = "gamgam"
MAP_HIGGS_DECAY_SMEFT["hzz"] = "ZZ"
MAP_HIGGS_DECAY_SMEFT["hbb"] = "bb"
MAP_HIGGS_DECAY_SMEFT["hww"] = "WW"
MAP_HIGGS_DECAY_SMEFT["htt"] = "tautau"
MAP_HIGGS_DECAY_SMEFT["hmm"] = "mumu"
MAP_HIGGS_DECAY_SMEFT["hzg"] = "Zgam"

MAP_HIGGS_PROD_SMEFT = od()
MAP_HIGGS_PROD_SMEFT["ggH"] = "GG2H"
MAP_HIGGS_PROD_SMEFT["qqH"] = "QQ2HQQ"
MAP_HIGGS_PROD_SMEFT["WH_had"] = "QQ2HQQ"
MAP_HIGGS_PROD_SMEFT["ZH_had"] = "QQ2HQQ"
MAP_HIGGS_PROD_SMEFT["ggZH_had"] = "GG2H"
MAP_HIGGS_PROD_SMEFT["ggZH_qq"] = "GG2H"
MAP_HIGGS_PROD_SMEFT["WH_lep"] = "QQ2HLNU"
MAP_HIGGS_PROD_SMEFT["ZH_lep"] = "QQ2HLL"
MAP_HIGGS_PROD_SMEFT["ggZH_lep"] = "GG2HLL"
MAP_HIGGS_PROD_SMEFT["ggZH_ll"] = "GG2HLL"
MAP_HIGGS_PROD_SMEFT["ggZH_nunu"] = "GG2HLL"
MAP_HIGGS_PROD_SMEFT["ttH"] = "TTH"
MAP_HIGGS_PROD_SMEFT["tHq"] = "THQ"
MAP_HIGGS_PROD_SMEFT["tHW"] = "THW"
MAP_HIGGS_PROD_SMEFT["bbH"] = "BBH"

MAP_HIGGS_PROD_SMEFT_ATLAS = od()
MAP_HIGGS_PROD_SMEFT_ATLAS["ggH"] = "GG2H"
MAP_HIGGS_PROD_SMEFT_ATLAS["qqH"] = "QQ2HQQ"
MAP_HIGGS_PROD_SMEFT_ATLAS["WH_had"] = "QQ2HQQ"
MAP_HIGGS_PROD_SMEFT_ATLAS["ZH_had"] = "QQ2HQQ"
MAP_HIGGS_PROD_SMEFT_ATLAS["ggZH_had"] = "GG2H"
MAP_HIGGS_PROD_SMEFT_ATLAS["ggZH_qq"] = "GG2H"
MAP_HIGGS_PROD_SMEFT_ATLAS["WH_lep"] = "QQ2HLNU"
MAP_HIGGS_PROD_SMEFT_ATLAS["ZH_lep"] = "QQ2HLL"
MAP_HIGGS_PROD_SMEFT_ATLAS["ggZH_lep"] = "GG2HLL"
MAP_HIGGS_PROD_SMEFT_ATLAS["ggZH_ll"] = "GG2HLL"
MAP_HIGGS_PROD_SMEFT_ATLAS["ggZH_nunu"] = "GG2HLL"
MAP_HIGGS_PROD_SMEFT_ATLAS["ttH"] = "TTH"
MAP_HIGGS_PROD_SMEFT_ATLAS["tHq"] = "TH"
MAP_HIGGS_PROD_SMEFT_ATLAS["tHW"] = "TH"
MAP_HIGGS_PROD_SMEFT_ATLAS["bbH"] = "BBH"


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Global function to extract reco category, STXS bin, decay mode and energy from process name
def getProcessInfo(bin, process):
    foundRecoCategory = bin
    foundSTXSBin = process
    foundDecay = None
    foundEnergy = "13TeV"
    # Iterate over Higgs decays
    matchedDecayString = False
    for D in ALL_HIGGS_DECAYS:
        if matchedDecayString:
            continue
        if "_%s" % D in foundSTXSBin:
            foundSTXSBin = re.sub("_%s" % D, "", foundSTXSBin)
            foundDecay = D
            matchedDecayString = True
    # Also drop year tag in STXS bin name if present
    for Y in ["2016", "2017", "2018"]:
        if "_%s" % Y in foundSTXSBin:
            foundSTXSBin = re.sub("_%s" % Y, "", foundSTXSBin)

    # Catch for H->Zgam
    if (foundDecay == "hzg") | ("bkg" in foundSTXSBin):
        foundSTXSBin = foundSTXSBin.split("_")[0]

    if not matchedDecayString:
        raise RuntimeError("Validation Error: no supported decay found in process")

    return (foundRecoCategory, foundSTXSBin, foundDecay, foundEnergy)


#################################################################################################################
# STXS to EFT abstract base class: inherited classes for different stages
class STXStoSMEFTBaseModel(SMLikeHiggsModel):
    def __init__(self, fixProcesses=[]):
        SMLikeHiggsModel.__init__(self)
        self.PROCESSES = None
        self.DECAYS = None
        # Dicts to store pois + scaling functions
        self.pois = None
        self.poiNameMap = {}  # To account for exponents in poi names, which aren't in input json files
        self.STXSScalingTerms = None
        self.DecayScalingTerms = None
        self.map_prod = MAP_HIGGS_PROD_SMEFT
        self.map_decay = MAP_HIGGS_DECAY_SMEFT
        # Options
        self.floatMass = False
        self.fixProcesses = fixProcesses  # Option to fix certain STXS bins: comma separated list of STXS bins
        self.stage0 = False
        self.parametrisation = "CMS-prelim-SMEFT-topU3l_22_05_05_AccCorr_0p01"
        self.eigenvalueThreshold = -1
        # Configuration of model: taylor expansion vs full, save only linear/linquad
        self.linear_only = False
        self.linquad_only = False
        self.expand_equations = False  # Option for Taylor-expanded equations

    def setPhysicsOptionsBase(self, physOptions):
        for po in physOptions:
            if po.startswith("higgsMassRange="):
                self.floatMass = True
                self.mHRange = po.replace("higgsMassRange=", "").split(",")
                if len(self.mHRange) != 2:
                    raise RuntimeError("Higgs mass range definition requires two extrema")
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError("Extrema for Higgs mass range defined with inverterd order. Second must be larger the first")
            if po.startswith("fixProcesses="):
                self.fixProcesses = (po.replace("fixProcesses=", "")).split(",")
            if po.startswith("stage0="):
                self.stage0 = po.replace("stage0=", "") in ["yes", "1", "Yes", "True", "true"]
            if po.startswith("parametrisation="):
                self.parametrisation = po.replace("parametrisation=", "")
                if "ATLAS" in self.parametrisation:
                    self.map_prod = MAP_HIGGS_PROD_SMEFT_ATLAS
            if po.startswith("eigenvalueThreshold="):
                self.eigenvalueThreshold = po.replace("eigenvalueThreshold=", "")
            if po.startswith("linear_only="):
                self.linear_only = po.replace("linear_only=", "") in ["yes", "1", "Yes", "True", "true"]
            if po.startswith("linquad_only="):
                self.linquad_only = po.replace("linquad_only=", "") in ["yes", "1", "Yes", "True", "true"]
            if po.startswith("expand_equations="):
                self.expand_equations = po.replace("expand_equations=", "") in ["yes", "1", "Yes", "True", "true"]

        # Output options to screen
        print(" --> [STXStoSMEFT] Using (%s) parametrisation" % self.parametrisation)
        if len(self.fixProcesses) > 0:
            print(" --> [STXStoSMEFT] Fixing following processes to SM: %s" % self.fixProcesses)

    def doMH(self):
        if self.floatMass:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]), float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
            else:
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0], self.mHRange[1]))
        else:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
                self.modelBuilder.out.var("MH").setConstant(True)
            else:
                self.modelBuilder.doVar("MH[%g]" % self.options.mass)

    # Overwrite getYieldScale to extract (RECO-category,STXS bin,decay,energy)
    def getYieldScale(self, bin, process):
        if not self.DC.isSignal[process]:
            return 1.0

        # Extract process line info
        (recocat, stxsbin, decay, energy) = getProcessInfo(bin, process)

        # Return 1 (no scaling) for fixed processes and scaling for non-fixed
        if stxsbin in self.fixProcesses:
            return 1.0
        else:
            procStr = stxsbin
            return self.getHiggsSignalYieldScale(procStr, decay, energy)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Extract pois from yaml file
    def extractPOIs(self, filename):
        with open(filename, "r") as fpois:
            try:
                self.pois = yaml.safe_load(fpois)
                # Apply eigenvector threshold if set
                if self.eigenvalueThreshold != -1.0:
                    pois_to_keep = {}
                    for poi, v in self.pois.items():
                        if "eigenvalue" in v:
                            if v["eigenvalue"] > float(self.eigenvalueThreshold):
                                pois_to_keep[poi] = v
                        else:
                            pois_to_keep[poi] = v
                    self.pois = pois_to_keep

            except yaml.YAMLERROR as exc:
                print(exc)

    # Function to extract STXS scaling terms from json file
    def extractSTXSScalingTerms(self, filename=""):
        if filename != "":
            with open(filename, "r") as jf:
                self.STXSScalingTerms = json.load(jf)
        else:
            self.STXSScalingTerms = {}

    # Function to extract decay scaling functions from file
    def extractDecayScalingTerms(self, filename=""):
        if filename != "":
            with open(filename, "r") as jf:
                self.DecayScalingTerms = json.load(jf)
        else:
            self.DecayScalingTerms = {}

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Function to make individual scaling functions in workspace
    def makeScalingFunction(self, what, is_decay=False, mode=""):
        mode_str = "_%s" % mode if mode != "" else ""

        # Apply mapping of production mode/decay to match inputs in json file
        k = what
        if is_decay:
            for D in self.map_decay.keys():
                if what == D:
                    k = self.map_decay[D]
        else:
            for P in self.map_prod.keys():
                if P in what:
                    k = re.sub(P, self.map_prod[P], what)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Pre-processing fix:
        # Fix for ttH multilepton: missing lep label in VH bins
        if "WH_PTV" in k:
            k = re.sub("WH", "QQ2HLNU", k)
        if "ZH_PTV" in k:
            k = re.sub("ZH", "QQ2HLL", k)
        if "ggZH_PTV" in k:
            k = re.sub("ggZH", "GG2HLL", k)

        # Fix for ttH multilepton: duplicate of proc names
        for P in self.map_prod.values():
            if "%s_%s" % (P, P) in k:
                k = re.sub("%s_%s" % (P, P), P, k)

        # Fix for VH procs without had/lep label: use leptonic scaling function. Is this accurate?
        if k == "WH":
            k = "QQ2HLNU"
        if k == "ZH":
            k = "QQ2HLL"
        if k == "ggZH":
            k = "GG2HLL"
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # Extract terms for dict
        if k in self.STXSScalingTerms:
            terms = self.STXSScalingTerms[k]
        elif k in self.DecayScalingTerms:
            terms = self.DecayScalingTerms[k]
        else:
            print(" --> [WARNING] Scaling terms for %s do not exist in input json. Setting to 1" % k)
            terms = {}

        # Loop over pois and extract the terms from scaling function, stored in C++ map
        coeffs = ROOT.std.map("string", "double")()
        list_pois = []

        for j, jpoi in enumerate(self.pois):
            if mode == "linear":
                jpoi_name = "l%s" % self.poiNameMap[jpoi]
            else:
                jpoi_name = self.poiNameMap[jpoi]
            e_jpoi = 10 ** (-1 * self.pois[jpoi]["exponent"])

            # Interference terms: Aj
            if "A_%s" % jpoi in terms:
                coeffs[jpoi_name] = e_jpoi * terms["A_%s" % jpoi]
                if jpoi_name not in list_pois:
                    list_pois.append(jpoi_name)

            if mode != "linear":
                # BSM-only terms: Bjk
                if "B_%s_2" % jpoi in terms:
                    coeffs["%s_2" % jpoi_name] = e_jpoi * e_jpoi * terms["B_%s_2" % jpoi]
                    if jpoi_name not in list_pois:
                        list_pois.append(jpoi_name)

                # Cross terms
                for k, kpoi in enumerate(self.pois):
                    if k > j:
                        if mode == "linear":
                            kpoi_name = "l%s" % self.poiNameMap[kpoi]
                        else:
                            kpoi_name = self.poiNameMap[kpoi]
                        e_kpoi = 10 ** (-1 * self.pois[kpoi]["exponent"])

                        if "B_%s_%s" % (jpoi, kpoi) in terms:
                            coeffs["%s_%s" % (jpoi_name, kpoi_name)] = e_jpoi * e_kpoi * terms["B_%s_%s" % (jpoi, kpoi)]
                        elif "B_%s_%s" % (kpoi, jpoi) in terms:
                            coeffs["%s_%s" % (jpoi_name, kpoi_name)] = e_jpoi * e_kpoi * terms["B_%s_%s" % (kpoi, jpoi)]

                        if kpoi_name not in list_pois:
                            list_pois.append(kpoi_name)

        # Make RooArgList of pois in equation
        arglist_pois = ROOT.RooArgList()
        for jpoi in list_pois:
            arglist_pois.add(self.modelBuilder.out.var(jpoi))

        # Make RooEFTScalingFunction
        if is_decay:
            if what != "tot":
                name = "scaling%s_partial_%s" % (mode_str, what)
            else:
                name = "scaling%s_%s" % (mode_str, what)
        else:
            name = "scaling%s_XS_%s" % (mode_str, what)
        eft_scaling = ROOT.RooEFTScalingFunction(name, name, coeffs, arglist_pois)

        # Add scaling function as RooAddition into model
        self.modelBuilder.out.Import(eft_scaling)

    # Function to make BR scaling functions: partial width/total width
    def makeBRScalingFunction(self, what, mode=""):
        mode_str = "_%s" % mode if mode != "" else ""
        self.modelBuilder.factory_('expr::scaling%s_BR_%s("@0/@1", scaling%s_partial_%s, scaling%s_tot)' % (mode_str, what, mode_str, what, mode_str))

    # Function to make Taylor-expanded scaling functions in workspace
    def makeScalingFunction_expand(self, what_production, what_decay):
        # Apply mapping of production mode/decay to match inputs in json file
        k_production = what_production
        for P in self.map_prod.keys():
            if P in what_production:
                k_production = re.sub(P, self.map_prod[P], what_production)

        k_decay = what_decay
        for D in self.map_decay.keys():
            if what_decay == D:
                k_decay = self.map_decay[D]

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Pre-processing fix:
        # Fix for ttH multilepton: missing lep label in VH bins
        if "WH_PTV" in k_production:
            k_production = re.sub("WH", "QQ2HLNU", k_production)
        if "ZH_PTV" in k_production:
            k_production = re.sub("ZH", "QQ2HLL", k_production)
        if "ggZH_PTV" in k_production:
            k_production = re.sub("ggZH", "GG2HLL", k_production)

        # Fix for ttH multilepton: duplicate of proc names
        for P in self.map_prod.values():
            if "%s_%s" % (P, P) in k_production:
                k_production = re.sub("%s_%s" % (P, P), P, k_production)

        # Fix for VH procs without had/lep label: use leptonic scaling function. Is this accurate?
        if k_production == "WH":
            k_production = "QQ2HLNU"
        if k_production == "ZH":
            k_production = "QQ2HLL"
        if k_production == "ggZH":
            k_production = "GG2HLL"
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # Extract terms from dict
        if k_production in self.STXSScalingTerms:
            terms_production = self.STXSScalingTerms[k_production]
        else:
            print(" --> [WARNING] Scaling terms for %s do not exist in input json. Setting to 1" % k_production)
            terms_production = {}

        if k_decay in self.DecayScalingTerms:
            terms_decay = self.DecayScalingTerms[k_decay]
        else:
            print(" --> [WARNING] Scaling terms for %s do not exist in input json. Setting to 1" % k_decay)
            terms_decay = {}

        if "tot" in self.DecayScalingTerms:
            terms_tot = self.DecayScalingTerms["tot"]
        else:
            print(" --> [WARNING] Scaling terms for Higgs total decay width (tot) do not exist in input json. Setting to 1")
            terms_tot = {}

        # Loop over pois and extract the terms from scaling function, stored in C++ map
        coeffs = ROOT.std.map("string", "double")()
        list_pois = []
        # Use different pois for linearised model so can store both in same workspace
        lcoeffs = ROOT.std.map("string", "double")()
        list_lpois = []

        for j, jpoi in enumerate(self.pois):
            Aj_sum = 0
            Bjj_sum = 0
            jpoi_name = self.poiNameMap[jpoi]
            jpoi_lname = "l%s" % self.poiNameMap[jpoi]
            e_jpoi = 10 ** (-1 * self.pois[jpoi]["exponent"])

            # Interference terms: O(c_j)
            if "A_%s" % jpoi in terms_production:
                Aj_sum += terms_production["A_%s" % jpoi]
            if "A_%s" % jpoi in terms_decay:
                Aj_sum += terms_decay["A_%s" % jpoi]
            if "A_%s" % jpoi in terms_tot:
                Aj_sum -= terms_tot["A_%s" % jpoi]
            # Multiply by exponent
            Aj_sum *= e_jpoi
            # If non-zero then add to coeffs
            if Aj_sum != 0.0:
                if jpoi_name not in list_pois:
                    list_pois.append(jpoi_name)
                coeffs[jpoi_name] = Aj_sum
                if jpoi_lname not in list_lpois:
                    list_lpois.append(jpoi_lname)
                lcoeffs[jpoi_lname] = Aj_sum

            if not self.linear_only:
                # Squared term: O(c_j^2)
                if "B_%s_2" % jpoi in terms_production:
                    Bjj_sum += terms_production["B_%s_2" % jpoi]
                if "B_%s_2" % jpoi in terms_decay:
                    Bjj_sum += terms_decay["B_%s_2" % jpoi]
                if "B_%s_2" % jpoi in terms_tot:
                    Bjj_sum -= terms_tot["B_%s_2" % jpoi]
                if ("A_%s" % jpoi in terms_production) & ("A_%s" % jpoi in terms_decay):
                    Bjj_sum += terms_production["A_%s" % jpoi] * terms_decay["A_%s" % jpoi]
                if ("A_%s" % jpoi in terms_production) & ("A_%s" % jpoi in terms_tot):
                    Bjj_sum -= terms_production["A_%s" % jpoi] * terms_tot["A_%s" % jpoi]
                if ("A_%s" % jpoi in terms_decay) & ("A_%s" % jpoi in terms_tot):
                    Bjj_sum -= terms_decay["A_%s" % jpoi] * terms_tot["A_%s" % jpoi]
                if "A_%s" % jpoi in terms_tot:
                    Bjj_sum += terms_tot["A_%s" % jpoi] * terms_tot["A_%s" % jpoi]
                # Multiply by exponent
                Bjj_sum *= e_jpoi * e_jpoi
                # If non zero then add to coeffs
                if Bjj_sum != 0.0:
                    if jpoi_name not in list_pois:
                        list_pois.append(jpoi_name)
                    coeffs["%s_2" % jpoi_name] = Bjj_sum

                # Cross terms
                for k, kpoi in enumerate(self.pois):
                    if k > j:
                        Bjk_sum = 0
                        kpoi_name = self.poiNameMap[kpoi]
                        e_kpoi = 10 ** (-1 * self.pois[kpoi]["exponent"])
                        if "B_%s_%s" % (jpoi, kpoi) in terms_production:
                            Bjk_sum += terms_production["B_%s_%s" % (jpoi, kpoi)]
                        elif "B_%s_%s" % (kpoi, jpoi) in terms_production:
                            Bjk_sum += terms_production["B_%s_%s" % (kpoi, jpoi)]
                        if "B_%s_%s" % (jpoi, kpoi) in terms_decay:
                            Bjk_sum += terms_decay["B_%s_%s" % (jpoi, kpoi)]
                        elif "B_%s_%s" % (kpoi, jpoi) in terms_decay:
                            Bjk_sum += terms_decay["B_%s_%s" % (kpoi, jpoi)]
                        if "B_%s_%s" % (jpoi, kpoi) in terms_tot:
                            Bjk_sum -= terms_tot["B_%s_%s" % (jpoi, kpoi)]
                        elif "B_%s_%s" % (kpoi, jpoi) in terms_tot:
                            Bjk_sum -= terms_tot["B_%s_%s" % (kpoi, jpoi)]
                        if ("A_%s" % jpoi in terms_production) & ("A_%s" % kpoi in terms_decay):
                            Bjk_sum += terms_production["A_%s" % jpoi] * terms_decay["A_%s" % kpoi]
                        if ("A_%s" % kpoi in terms_production) & ("A_%s" % jpoi in terms_decay):
                            Bjk_sum += terms_production["A_%s" % kpoi] * terms_decay["A_%s" % jpoi]
                        if ("A_%s" % jpoi in terms_production) & ("A_%s" % kpoi in terms_tot):
                            Bjk_sum -= terms_production["A_%s" % jpoi] * terms_tot["A_%s" % kpoi]
                        if ("A_%s" % kpoi in terms_production) & ("A_%s" % jpoi in terms_tot):
                            Bjk_sum -= terms_production["A_%s" % kpoi] * terms_tot["A_%s" % jpoi]
                        if ("A_%s" % jpoi in terms_decay) & ("A_%s" % kpoi in terms_tot):
                            Bjk_sum -= terms_decay["A_%s" % jpoi] * terms_tot["A_%s" % kpoi]
                        if ("A_%s" % kpoi in terms_decay) & ("A_%s" % jpoi in terms_tot):
                            Bjk_sum -= terms_decay["A_%s" % kpoi] * terms_tot["A_%s" % jpoi]
                        if ("A_%s" % jpoi in terms_tot) & ("A_%s" % kpoi in terms_tot):
                            Bjk_sum += 2 * terms_tot["A_%s" % jpoi] * terms_tot["A_%s" % kpoi]
                        # Multiply by exponent
                        Bjk_sum *= e_jpoi * e_kpoi
                        if Bjk_sum != 0.0:
                            if jpoi_name not in list_pois:
                                list_pois.append(jpoi_name)
                            if kpoi_name not in list_pois:
                                list_pois.append(kpoi_name)
                            coeffs["%s_%s" % (jpoi_name, kpoi_name)] = Bjk_sum

        if not self.linquad_only:
            # Make RooArgList of pois in equation
            arglist_lpois = ROOT.RooArgList()
            for jpoi in list_lpois:
                arglist_lpois.add(self.modelBuilder.out.var(jpoi))
            # Make RooEFTScalingFunction
            name = "scaling_linear_expand_XS_%s_BR_%s" % (what_production, what_decay)
            eft_scaling_linear = ROOT.RooEFTScalingFunction(name, name, lcoeffs, arglist_lpois)
            # Add scaling function into model
            self.modelBuilder.out.Import(eft_scaling_linear)

        if not self.linear_only:
            # Make RooArgList of pois in equation
            arglist_pois = ROOT.RooArgList()
            for jpoi in list_pois:
                arglist_pois.add(self.modelBuilder.out.var(jpoi))
            # Make RooEFTScalingFunction
            name = "scaling_linquad_expand_XS_%s_BR_%s" % (what_production, what_decay)
            eft_scaling_linquad = ROOT.RooEFTScalingFunction(name, name, coeffs, arglist_pois)
            # Add scaling function into model
            self.modelBuilder.out.Import(eft_scaling_linquad)


#################################################################################################################
class STXSToSMEFTModel(STXStoSMEFTBaseModel):
    def __init__(self):
        STXStoSMEFTBaseModel.__init__(self)

    def setPhysicsOptions(self, physOptions):
        self.setPhysicsOptionsBase(physOptions)

    def doParametersOfInterest(self):
        if self.floatMass:
            print(" --> [WARNING] Floating Higgs mass selected. STXStoSMEFT model assumes MH=125.0 GeV")
        self.doMH()
        self.SMH = SMHiggsBuilder(self.modelBuilder)

        # Read in parameters of interest from yaml file
        self.extractPOIs("%s/src/HiggsAnalysis/CombinedLimit/data/eft/STXStoSMEFT/%s/pois.yaml" % (os.environ["CMSSW_BASE"], self.parametrisation))

        # Create list of pois and build RooRealVars
        POIs, lPOIs = [], []
        for poi in self.pois:
            if self.pois[poi]["exponent"] < 0:
                poi_name = "%sXEm%g" % (poi, abs(self.pois[poi]["exponent"]))
            elif self.pois[poi]["exponent"] >= 1:
                poi_name = "%sXE%g" % (poi, self.pois[poi]["exponent"])
            else:
                poi_name = poi
            self.poiNameMap[poi] = poi_name

            if not self.linear_only:
                POIs.append(poi_name)
                self.modelBuilder.doVar("%s[%g,%g,%g]" % (poi_name, self.pois[poi]["val"], self.pois[poi]["min"], self.pois[poi]["max"]))
                self.modelBuilder.out.var(poi_name).setConstant(True)

            if not self.linquad_only:
                lpoi_name = "l%s" % poi_name
                lPOIs.append(lpoi_name)
                self.modelBuilder.doVar("%s[%g,%g,%g]" % (lpoi_name, self.pois[poi]["val"], self.pois[poi]["min"], self.pois[poi]["max"]))
                self.modelBuilder.out.var(lpoi_name).setConstant(True)

        # Create POI sets
        if self.linear_only:
            self.modelBuilder.doSet("POI", ",".join(lPOIs))
            self.POIs = ROOT.RooArgList(self.modelBuilder.out.set("POI"))
        elif self.linquad_only:
            self.modelBuilder.doSet("POI", ",".join(POIs))
            self.POIs = ROOT.RooArgList(self.modelBuilder.out.set("POI"))
        else:
            self.modelBuilder.doSet("POI", ",".join(POIs))
            self.POIs = ROOT.RooArgList(self.modelBuilder.out.set("POI"))
            self.modelBuilder.doSet("lPOI", ",".join(lPOIs))
            self.lPOIs = ROOT.RooArgList(self.modelBuilder.out.set("lPOI"))

        # set up model
        self.setup()

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def setup(self):
        # Extract scaling terms from json files: inclusive vs reco-level
        self.extractSTXSScalingTerms(
            filename="%s/src/HiggsAnalysis/CombinedLimit/data/eft/STXStoSMEFT/%s/prod.json" % (os.environ["CMSSW_BASE"], self.parametrisation)
        )
        self.extractDecayScalingTerms(
            filename="%s/src/HiggsAnalysis/CombinedLimit/data/eft/STXStoSMEFT/%s/decay.json" % (os.environ["CMSSW_BASE"], self.parametrisation)
        )

        if not self.expand_equations:
            if not self.linquad_only:
                self.makeScalingFunction("tot", is_decay=True, mode="linear")
            if not self.linear_only:
                self.makeScalingFunction("tot", is_decay=True, mode="linquad")

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getHiggsSignalYieldScale(self, production, decay, energy):
        # Function to convert troublesome procs into viable one for HC combination
        production = convert_to_STXS(production, decay)

        # Stage 0 option: use stage 0 bin scaling
        if self.stage0:
            for P in self.map_prod.keys():
                if P in production:
                    production = P

        # For model with Taylor-expanded equations
        if self.expand_equations:
            if self.linear_only:
                name = "scaling_linear_expand_XS_%s_BR_%s" % (production, decay)
                if self.modelBuilder.out.function(name) == None:
                    print(" --> [STXStoSMEFT] Making linearised model for (STXS bin,decay): (%s,%s)" % (production, decay))
                    self.makeScalingFunction_expand(production, decay)

            elif self.linquad_only:
                name = "scaling_linquad_expand_XS_%s_BR_%s" % (production, decay)
                if self.modelBuilder.out.function(name) == None:
                    print(" --> [STXStoSMEFT] Making linear+quad model for (STXS bin,decay): (%s,%s)" % (production, decay))
                    self.makeScalingFunction_expand(production, decay)

            else:
                name = "stxstosmeft_scaling_expand_XS_%s_BR_%s" % (production, decay)
                name_lin = "scaling_linear_expand_XS_%s_BR_%s" % (production, decay)
                name_linquad = "scaling_linquad_expand_XS_%s_BR_%s" % (production, decay)
                if (self.modelBuilder.out.function(name_lin) == None) | (self.modelBuilder.out.function(name_linquad) == None):
                    print(" --> [STXStoSMEFT] Making linearised and linear+quadratic model for (STXS bin,decay): (%s,%s)" % (production, decay))
                    self.makeScalingFunction_expand(production, decay)
                # Combine linear and linear+quadratic models into same scaling (use different POIs)
                self.modelBuilder.factory_("prod::%s(%s)" % (name, ",".join([name_lin, name_linquad])))

        # For model with full equations: xs*partial/total
        else:
            if self.linear_only:
                name_lin_xs = "scaling_linear_XS_%s" % production
                if self.modelBuilder.out.function(name_lin_xs) == None:
                    print(" --> [STXStoSMEFT] Making linear scaling function for STXS bin: %s" % production)
                    self.makeScalingFunction(production, mode="linear")
                name_lin_br = "scaling_linear_BR_%s" % decay
                if self.modelBuilder.out.function(name_lin_br) == None:
                    print(" --> [STXStoSMEFT] Making linear scaling function for decay: %s" % decay)
                    self.makeScalingFunction(decay, is_decay=True, mode="linear")
                    self.makeBRScalingFunction(decay, mode="linear")
                name = "scaling_linear_XS_%s_BR_%s" % (production, decay)
                if self.modelBuilder.out.function(name) == None:
                    self.modelBuilder.factory_("prod::%s(%s)" % (name, ",".join([name_lin_xs, name_lin_br])))

            elif self.linquad_only:
                name_linquad_xs = "scaling_linquad_XS_%s" % production
                if self.modelBuilder.out.function(name_linquad_xs) == None:
                    print(" --> [STXStoSMEFT] Making linear+quad scaling function for STXS bin: %s" % production)
                    self.makeScalingFunction(production, mode="linquad")
                name_linquad_br = "scaling_linquad_BR_%s" % decay
                if self.modelBuilder.out.function(name_linquad_br) == None:
                    print(" --> [STXStoSMEFT] Making linear+quad scaling function for decay: %s" % decay)
                    self.makeScalingFunction(decay, is_decay=True, mode="linquad")
                    self.makeBRScalingFunction(decay, mode="linquad")
                name = "scaling_linquad_XS_%s_BR_%s" % (production, decay)
                if self.modelBuilder.out.function(name) == None:
                    self.modelBuilder.factory_("prod::%s(%s)" % (name, ",".join([name_linquad_xs, name_linquad_br])))

            # Combine linear and linear+quadratic models into same scaling (use different POIs) as product
            else:
                # Linear part
                name_lin_xs = "scaling_linear_XS_%s" % production
                if self.modelBuilder.out.function(name_lin_xs) == None:
                    print(" --> [STXStoSMEFT] Making linear scaling function for STXS bin: %s" % production)
                    self.makeScalingFunction(production, mode="linear")
                name_lin_br = "scaling_linear_BR_%s" % decay
                if self.modelBuilder.out.function(name_lin_br) == None:
                    print(" --> [STXStoSMEFT] Making linear scaling function for decay: %s" % decay)
                    self.makeScalingFunction(decay, is_decay=True, mode="linear")
                    self.makeBRScalingFunction(decay, mode="linear")
                name_lin = "scaling_linear_XS_%s_BR_%s" % (production, decay)
                if self.modelBuilder.out.function(name_lin) == None:
                    self.modelBuilder.factory_("prod::%s(%s)" % (name_lin, ",".join([name_lin_xs, name_lin_br])))

                # Linear+quadratic part
                name_linquad_xs = "scaling_linquad_XS_%s" % production
                if self.modelBuilder.out.function(name_linquad_xs) == None:
                    print(" --> [STXStoSMEFT] Making linear+quad scaling function for STXS bin: %s" % production)
                    self.makeScalingFunction(production, mode="linquad")
                name_linquad_br = "scaling_linquad_BR_%s" % decay
                if self.modelBuilder.out.function(name_linquad_br) == None:
                    print(" --> [STXStoSMEFT] Making linear+quad scaling function for decay: %s" % decay)
                    self.makeScalingFunction(decay, is_decay=True, mode="linquad")
                    self.makeBRScalingFunction(decay, mode="linquad")
                name_linquad = "scaling_linquad_XS_%s_BR_%s" % (production, decay)
                if self.modelBuilder.out.function(name_linquad) == None:
                    self.modelBuilder.factory_("prod::%s(%s)" % (name_linquad, ",".join([name_linquad_xs, name_linquad_br])))

                # Combination
                name = "stxstosmeft_scaling_XS_%s_BR_%s" % (production, decay)
                if self.modelBuilder.out.function(name) == None:
                    self.modelBuilder.factory_("prod::%s(%s)" % (name, ",".join([name_lin, name_linquad])))

        return name


#################################################################################################################
# Function to convert troublesome procs to the names in the json files
def convert_to_STXS(_production, _decay):
    # Add string replace functions
    return _production


#################################################################################################################
# Instantiation of STXStoSMEFT model
STXStoSMEFT = STXSToSMEFTModel()
