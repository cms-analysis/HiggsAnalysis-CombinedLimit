import gzip
import json
import pickle
import re

import numpy as np
import ROOT
from HiggsAnalysis.CombinedLimit.PhysicsModel import PhysicsModelBase_NiceSubclasses


def read_scaling(path):
    if path.endswith(".json"):
        with open(path) as fin:
            out = json.load(fin)
    elif path.endswith(".json.gz"):
        with gzip.open(path, "rt") as fin:
            out = json.load(fin)
    elif path.endswith(".pkl.gz"):
        with gzip.open(path, "rb") as fin:
            out = pickle.load(fin)
    else:
        raise RuntimeError("Unrecognized scaling data path {path}; must be either .json, .json.gz, or .pkl.gz".format(path=path))
    # normalize
    if not isinstance(out, list):
        raise RuntimeError("Scaling data in invalid format: expected list")
    expected_fields = {"channel", "process", "parameters", "scaling"}
    for item in out:
        if not isinstance(item, dict):
            raise RuntimeError("Scaling data in invalid format: expected each item in list to be a dict")
        missing = expected_fields - set(item)
        if missing:
            raise RuntimeError("Missing fields in scaling item: {missing}".format(missing=missing))
        shortname = item["channel"] + "/" + item["process"]
        if not all(isinstance(par, str) for par in item["parameters"]):
            raise RuntimeError("Parameters must be a list of strings in {shortname}".format(shortname=shortname))
        try:
            # coerce into numpy with C-contiguous memory (needed for fast std::vector copy)
            item["scaling"] = np.ascontiguousarray(item["scaling"], dtype=float)
        except ValueError:
            # python3: raise from ex
            raise RuntimeError("Scaling data invalid: could not normalize array for {shortname}".format(shortname=shortname))
        if len(item["scaling"].shape) != 2:
            raise RuntimeError("Scaling data invalid: array shape incorrect for {shortname}".format(shortname=shortname))
        npar = len(item["parameters"])
        if item["scaling"].shape[1] != npar * (npar + 1) // 2:
            raise RuntimeError("Scaling data invalid: array has insufficent terms for parameters in {shortname}".format(shortname=shortname))
    return out


class InterferenceModel(PhysicsModelBase_NiceSubclasses):
    def __init__(self):
        self.verbose = False
        self.scaling = None
        self.scaling_map = {}
        self.explict_pois = None
        super(InterferenceModel, self).__init__()

    def processPhysicsOptions(self, physOptions):
        processed = []
        for po in physOptions:
            if po == "verbose":
                self.verbose = True
                processed.append(po)
            elif po.startswith("scalingData="):
                self.scaling = read_scaling(po[len("scalingData=") :])
                processed.append(po)
            elif po.startswith("POIs="):
                self.explict_pois = po[len("POIs=") :].split(":")
                processed.append(po)
        if self.scaling is None:
            raise RuntimeError("Must specify --PO=scalingData= physics option!")
        # make map for quick lookup
        for item in self.scaling:
            self.scaling_map[(item["channel"], item["process"])] = item
        return processed + super(InterferenceModel, self).processPhysicsOptions(physOptions)

    def getPOIList(self):
        poiNames = []
        if self.explict_pois:
            if self.verbose:
                print("Building explicitly requested POIs:")
            for poi in self.explict_pois:
                if self.verbose:
                    print(" - {poi}".format(poi=poi))
                self.modelBuilder.doVar(poi)
                poiname = re.sub(r"\[.*", "", poi)
                if not self.modelBuilder.out.var(poiname):
                    raise RuntimeError("Could not find POI {poiname} after factory call".format(poiname=poiname))
                poiNames.append(poiname)
        # RooFit would also detect duplicate params with mismatched factory, but this is faster
        constructed = {}
        for item in self.scaling:
            item["parameters_inworkspace"] = []
            for param in item["parameters"]:
                if param in constructed:
                    rooparam = constructed[param]
                else:
                    rooparam = self.modelBuilder.factory_(param)
                    if not rooparam:
                        raise RuntimeError("Failed to build {param}".format(param=param))
                    if not self.explict_pois and isinstance(rooparam, ROOT.RooRealVar) and not rooparam.isConstant():
                        if self.verbose:
                            print("Assuming parameter {param} is to be a POI (name: {name})".format(param=param, name=rooparam.GetName()))
                        poiNames.append(rooparam.GetName())
                    constructed[param] = rooparam
                # save for later use in making CMSInterferenceFunc
                item["parameters_inworkspace"].append(rooparam)
        if self.verbose:
            print("All POIs: {poiNames}".format(poiNames=poiNames))
        return poiNames

    def getYieldScale(self, channel, process):
        string = "%s/%s" % (channel, process)
        try:
            item = self.scaling_map[(channel, process)]
        except KeyError:
            if self.verbose:
                print("Will scale {string} by 1".format(string=string))
            return 1
        print("Will scale {string} using CMSInterferenceFunc dependent on parameters {parameters}".format(string=string, parameters=item["parameters"]))
        item["found"] = True
        # We don't actually scale the total via this mechanism
        # and we'll set up the shape effects in done() since shapes aren't available yet
        return 1

    def done(self):
        for item in self.scaling_map.values():
            channel = item["channel"]
            process = item["process"]
            string = "%s/%s" % (channel, process)
            if "found" not in item:
                if self.verbose:
                    print("Did not find {string} in workspace, even though it is in scaling data".format(string=string))
                continue

            hfname = "shapeSig_{channel}_{process}_morph".format(channel=channel, process=process)
            histfunc = self.modelBuilder.out.function(hfname)
            if not histfunc:
                # if there are no systematics, it ends up with a different name?
                hfname = "shapeSig_{process}_{channel}_rebinPdf".format(channel=channel, process=process)
                histfunc = self.modelBuilder.out.function(hfname)
            if not histfunc:
                # assume this is a CMSHistSum workspace
                histsum_name = "prop_bin{channel}".format(channel=channel)
                histfunc = self.modelBuilder.out.function(histsum_name)
                # in this workspace there is no object representing the process morphing
                # make one up so that funcname is still unique
                hfname = "prop_bin{channel}_process{process}".format(channel=channel, process=process)
            # TODO: support FastVerticalInterpHistPdf2
            if not isinstance(histfunc, (ROOT.CMSHistFunc, ROOT.CMSHistSum)):
                raise RuntimeError(
                    "Could not locate the CMSHistFunc or CMSHistSum for {string}.\n".format(string=string)
                    + "Note that CMSInterferenceFunc currently only supports workspaces that use these classes"
                )

            funcname = hfname + "_externalMorph"

            tpl = histfunc.cache()
            edges = [tpl.GetEdge(i) for i in range(tpl.size() + 1)]
            params = ROOT.RooArgList()
            for p in item["parameters_inworkspace"]:
                params.add(p)

            scaling_array = ROOT.std.vector["std::vector<double>"]()
            nbins, ncoef = item["scaling"].shape
            scaling_array.reserve(nbins)
            for sbin in item["scaling"]:
                scaling_array.push_back(sbin)

            self.modelBuilder.out.safe_import(ROOT.CMSInterferenceFunc(funcname, "", histfunc.getXVar(), edges, params, scaling_array))
            func = self.modelBuilder.out.function(funcname)
            if isinstance(histfunc, ROOT.CMSHistFunc):
                histfunc.injectExternalMorph(func)
            elif isinstance(histfunc, ROOT.CMSHistSum):
                postfix = "bin{channel}_proc_{process}".format(channel=channel, process=process)
                for idx, coef in enumerate(histfunc.coefList()):
                    if coef.GetName() in ("n_exp_" + postfix, "n_exp_final_" + postfix):
                        if self.verbose:
                            print("Injecting external morph for " + funcname)
                        histfunc.injectExternalMorph(idx, func)
                        break


interferenceModel = InterferenceModel()
