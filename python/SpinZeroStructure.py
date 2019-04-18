import collections, itertools, math
from HiggsAnalysis.CombinedLimit.PhysicsModel import *

### This is the base python class to study the SpinZero structure

class SpinZeroHiggsBase(PhysicsModelBase_NiceSubclasses):
    def __init__(self):
        super(SpinZeroHiggsBase, self).__init__()

        self.faiphiaistatus = collections.defaultdict(lambda: "fix")

        self.allowPMF = False

        self.HWWcombination = False

        self.offshell = False

    def setModelBuilder(self, modelBuilder):
        super(SpinZeroHiggsBase, self).setModelBuilder(modelBuilder)
        self.modelBuilder.doModelBOnly = False

    def getYieldScale(self,bin,process):
        "Split in production and decay, and call getHiggsSignalYieldScale; return 1 for backgrounds "
        result = super(SpinZeroHiggsBase, self).getYieldScale(bin, process)
        if result not in (0, 1): print "Process {0} will scale by {1}".format(process,result)
        return result

    def processPhysicsOptions(self,physOptions):
        processed = super(SpinZeroHiggsBase, self).processPhysicsOptions(physOptions)
        for po in physOptions:
            match = re.match("(f|phi)ai([0-9]+)(fixed|notpoi|floating|aspoi)$", po.lower())
            if match:
                parametertype = match.group(1)
                i = int(match.group(2))
                whattodo = match.group(3)

                if not 1 <= i <= self.numberoffais:
                    raise ValueError("There are only {} fais available, so can't do anything with {}".format(i, po))

                key = parametertype, i

                if key in self.faiphiaistatus:
                    raise ValueError("Specified multiple physics options for {}ai{}".format(parametertype, i))

                if whattodo == "fixed":
                    self.faiphiaistatus[key] = "fix"
                    if parametertype == "f" and i == 1:
                        print "Will fix CMS_zz4l_{}ai{} to 0".format(parametertype, i)
                elif whattodo == "floating" or whattodo == "notpoi":
                    self.faiphiaistatus[key] = "float"
                    if parametertype == "f" and i == 1:
                        print "CMS_zz4l_{}ai{} is NOT A POI".format(parametertype, i)
                    else:
                        print "Will float CMS_zz4l_{}ai{}".format(parametertype, i)
                elif whattodo == "aspoi":
                    self.faiphiaistatus[key] == "POI"
                    print "Will consider CMS_zz4l_{}ai{} as a parameter of interest".format(parametertype, i)
                else:
                    assert False, whattodo

                processed.append(po)

            if po.lower() == 'allowpmf':
                self.allowPMF = True
                processed.append(po)

            if po.lower() == 'hwwcombination':
                self.HWWcombination = True
                processed.append(po)

            if po.lower() == "offshell":
                self.offshell = True
                processed.append(po)

        return processed

    @property
    def numberoffais(self):
        return 2

    def getPOIList(self):

        poi = []
        poi += super(SpinZeroHiggsBase, self).getPOIList()

        if ("f", 1) not in self.faiphiaistatus: self.faiphiaistatus["f", 1] = "POI"

        for parametertype, parameterrange in ("f", (-1. if self.allowPMF else 0., 1.)), ("phi", (-math.pi, math.pi)):
            for i in xrange(1, self.numberoffais+1):
                varname = "CMS_zz4l_{}ai{}".format(parametertype, i)
                status = self.faiphiaistatus[parametertype, i]

                if status in ("float", "POI"):
                    if self.modelBuilder.out.var(varname):
                        self.modelBuilder.out.var(varname).setRange(*parameterrange)
                        self.modelBuilder.out.var(varname).setVal(0)
                    else:
                        self.modelBuilder.doVar("{}[0.,{},{}]".format(varname, *parameterrange))
                    self.modelBuilder.out.var(varname).setConstant(False)
                    if status == "POI":
                        print "Treating "+varname+" as a POI"
                        poi.append(varname)
                    else:
                        print "Floating "+varname
                        self.modelBuilder.out.var(varname).setAttribute("flatParam")
                    if parametertype == "f" and self.allowPMF: print "Allowing negative "+varname
                elif status == "fix":
                    if self.modelBuilder.out.var(varname):
                        self.modelBuilder.out.var(varname).setVal(0)
                        self.modelBuilder.out.var(varname).setConstant()
                    else:
                        self.modelBuilder.doVar("{}[0]".format(varname))
                    print "Fixing "+varname
                else:
                    assert False, status

        if self.HWWcombination:
            if self.modelBuilder.out.var("CMS_zz4l_alpha"):
                print "Found CMS_zz4l_alpha; setting to 0.5"
                self.modelBuilder.out.var("CMS_zz4l_alpha").setVal(0.5)
            else:
                print "Creating CMS_zz4l_alpha; setting to 0.5"
                self.modelBuilder.doVar("CMS_zz4l_alpha[0.5,-1,1]")
            print "Treating alpha (Rai1) as a POI"
            self.modelBuilder.out.var("CMS_zz4l_alpha").setConstant(False)
            poi.append("CMS_zz4l_alpha")
        else:
            if self.modelBuilder.out.var("CMS_zz4l_alpha"):
                print "Found CMS_zz4l_alpha; setting to constant 0"
                self.modelBuilder.out.var("CMS_zz4l_alpha").setVal(0)
                self.modelBuilder.out.var("CMS_zz4l_alpha").setConstant()

        # GGsm offshell variable
        if self.modelBuilder.out.var("CMS_zz4l_GGsm"):
            print "CMS_zz4l_GGsm is being renamed to GGsm"
            self.modelBuilder.out.var("CMS_zz4l_GGsm").SetName("GGsm")
        if self.offshell:
            if self.modelBuilder.out.var("GGsm"):
                print "Found GGsm, setting range to [1., 0., 50.]"
                self.modelBuilder.out.var("GGsm").setRange(0, 50)
                self.modelBuilder.out.var("GGsm").setVal(1)
            else:
                print "Creating GGsm; setting range to [1., 0., 50.]"
                self.modelBuilder.doVar("GGsm[1,0,50]")
            print "GGsm is a POI."
            self.modelBuilder.out.var("GGsm").setConstant(False)
            poi.append("GGsm")
        else:
            if self.modelBuilder.out.var("GGsm"):
                print "Found GGsm, fixing to 1"
                self.modelBuilder.out.var("GGsm").setVal(1)
                self.modelBuilder.out.var("GGsm").setConstant()

        return poi

class SpinZeroHiggs(SpinZeroHiggsBase):
    def __init__(self):
        super(SpinZeroHiggs, self).__init__()
        self.muFloating = True
        self.muAsPOI = False

    def processPhysicsOptions(self,physOptions):
        processed = super(SpinZeroHiggs, self).processPhysicsOptions(physOptions)
        for po in physOptions:
            if po.lower() == 'mufixed':
                print "Will consider the signal strength as a fixed parameter"
                self.muFloating = False
                processed.append(po)
            elif po.lower() == 'muaspoi':
                print "Will consider the signal strength as a parameter of interest"
                self.muAsPOI = True
                processed.append(po)

        if self.muAsPOI and not self.muFloating:
            raise ValueError("Specified both muFixed and muAsPOI!")

        if not self.muAsPOI and not self.fai1POI and not self.fai2POI and not self.phiai1POI and not self.phiai2POI:
            print "No POIs detected: Switching to default configuration: Floating nuisance mu, floating POI fai1, eveything else fixed"
            self.muFloating = True
            self.muAsPOI = False
            self.fai1Floating = True
            self.fai1POI = True
            self.phiai1Floating = False
            self.phiai1POI = False
            self.fai2Floating = False
            self.fai2POI = False
            self.phiai2Floating = False
            self.phiai2POI = False
            self.allowPMF = False

        return processed

    def getPOIList(self):
        poi = super(SpinZeroHiggs, self).getPOIList()
        if self.muFloating:
            if self.modelBuilder.out.var("r"):
                self.modelBuilder.out.var("r").setRange(0.,400.)
                self.modelBuilder.out.var("r").setVal(1)
            else:
                self.modelBuilder.doVar("r[1,0,400]")
            if self.HWWcombination:
                self.modelBuilder.out.var("r").removeMax()
                print "Removed maximum of r"

            self.modelBuilder.out.var("r").setConstant(False)

            if self.muAsPOI:
                print "Treating r as a POI"
                poi.append("r")
            else:
                self.modelBuilder.out.var("r").setAttribute("flatParam")
        else:
            if self.modelBuilder.out.var("r"):
                self.modelBuilder.out.var("r").setVal(1)
                self.modelBuilder.out.var("r").setConstant()
            else:
                self.modelBuilder.doVar("r[1]")
        return poi

class MultiSignalSpinZeroHiggs(SpinZeroHiggsBase,CanTurnOffBkgModel,MultiSignalModel):
    def __init__(self):
        super(MultiSignalSpinZeroHiggs, self).__init__()

        self.scalemuvfseparately = True
        self.scaledifferentsqrtsseparately = False
        self.uservoverrf = False
        self.sqrts = None
        self.fixed = []
        self.floated = []
        self.noRV = self.noRF = True
        #not doing muAsPOI or fixMu, there are too many permutations.
        #should just be set when running combine.

    def setPhysicsOptions(self, physOptions):
        if not any(po.startswith("map=") for po in physOptions):
            #no po started with map --> no manual overriding --> use the defaults
            #can still override with e.g. turnoff=ZH,WH
            physOptions = ["map=.*/(gg|qq|Z|W|tt|bb)H$:1"] + physOptions
        super(MultiSignalSpinZeroHiggs, self).setPhysicsOptions(physOptions)

    def processPhysicsOptions(self, physOptions):
        processed = super(MultiSignalSpinZeroHiggs, self).processPhysicsOptions(physOptions)
        for po in physOptions:
            if po.lower() == "scalemuvmuftogether":
                self.scalemuvfseparately = False
                processed.append(po)
            if po.lower() == "scaledifferentsqrtsseparately":
                self.scaledifferentsqrtsseparately = True
                processed.append(po)
            if po.lower().startswith("sqrts="):
                if self.sqrts is not None: raise ValueError("Duplicate physicsoption sqrts=?? provided")
                self.sqrts = [int(_) for _ in po.replace("sqrts=", "").split(",")]
                processed.append(po)
            if po.lower() == "uservoverrf":
                self.uservoverrf = True
                processed.append(po)
            if po.lower() == "norv":
                self.noRV = True
                processed.append(po)
            if po.lower() == "norf":
                self.noRF = True
                processed.append(po)

        if self.uservoverrf and not self.scalemuvfseparately:
            raise ValueError("can't specify both uservoverrf and scalemuvmuftogether")

        if self.sqrts is None:
            raise ValueError("PhysicsOption sqrts=?? is mandatory.  example: sqrts=7,8,13")

        if self.scaledifferentsqrtsseparately and self.scalemuvfseparately:
            if self.uservoverrf:
                self.fixed = ["RV", "RF", "R"] + ["RF_{}TeV".format(_) for _ in self.sqrts]
                self.floated = ["R{}_{}TeV".format(_1, _2) for _1 in ("V", "") for _2 in self.sqrts]
            else:
                self.fixed = ["RV", "RF", "R"] + ["R_{}TeV".format(_) for _ in self.sqrts]
                self.floated = ["R{}_{}TeV".format(_1, _2) for _1 in ("V", "F") for _2 in self.sqrts]
        elif self.scaledifferentsqrtsseparately and not self.scalemuvfseparately:
            self.fixed = ["RV", "RF", "R"] + ["R{}_{}TeV".format(_1, _2) for _1 in ("V", "F") for _2 in self.sqrts]
            self.floated = ["R_{}TeV".format(_) for _ in self.sqrts]
        elif not self.scaledifferentsqrtsseparately and self.scalemuvfseparately:
            if self.uservoverrf:
                self.fixed = ["RF"] + ["R{}_{}TeV".format(_1, _2) for _1 in ("V", "F", "") for _2 in self.sqrts]
                self.floated = ["RV", "R"]
            else:
                self.fixed = ["R"] + ["R{}_{}TeV".format(_1, _2) for _1 in ("V", "F", "") for _2 in self.sqrts]
                self.floated = ["RV", "RF"]
        elif not self.scaledifferentsqrtsseparately and not self.scalemuvfseparately:
            self.fixed = ["RV", "RF"] + ["R{}_{}TeV".format(_1, _2) for _1 in ("V", "F", "") for _2 in self.sqrts]
            self.floated = ["R"]
        else:
            assert False, "?????"

        if self.noRV:
            self.fixed = [_ for _ in self.fixed if "RV" not in _]
        if self.noRF:
            self.fixed = [_ for _ in self.fixed if "RF" not in _]

        return processed

    def getPOIList(self):
        result = super(MultiSignalSpinZeroHiggs, self).getPOIList()

        fixedorfloated = self.fixed+self.floated
        for variable in fixedorfloated:
            if not self.modelBuilder.out.var(variable):
                if variable in self.fixed: continue
                raise RuntimeError("{} does not exist in the workspace!  Check:\n - your datacard maker\n - your sqrts option".format(variable))
            else:
                if 'r' in variable.lower():
                    print "Setting {} range to [1.,0.,400.]".format(variable)
                    self.modelBuilder.out.var(variable).setRange(0.,400.)
                    self.modelBuilder.out.var(variable).setVal(1)
                elif "ggsm" in variable.lower():
                    print "Setting {} range to [1.,0.,50.]".format(variable)
                    self.modelBuilder.out.var(variable).setRange(0.,50.)
                    self.modelBuilder.out.var(variable).setVal(1)
                else:
                    print "Setting {} value to 0".format(variable)
                    self.modelBuilder.out.var(variable).setVal(0)
        for variable in self.fixed:
            if self.modelBuilder.out.var(variable):
                print "Fixing {}".format(variable)
                self.modelBuilder.out.var(variable).setConstant()
        for variable in self.floated:
            if self.modelBuilder.out.var(variable):
                print "Floating {} and assigning attribute flatParam".format(variable)
                self.modelBuilder.out.var(variable).setConstant(False)
                self.modelBuilder.out.var(variable).setAttribute("flatParam")

        return result

class HZZAnomalousCouplingsFromHistograms(MultiSignalSpinZeroHiggs):
    """
    General idea:
    This class expects histograms (which could be TH1 or RooDataHist) for each component of the PDF.
    They are supposed to be normalized to ai = 1.
    The pure components should be named ggH_0PM, qqH_0L1Zg, etc.

    The interference components should be split in two, one for positive bins and one for negative bins flipped.
    This way all bins are positive.
    They should be named ZH_g13g1prime21_positive or ttH_g11g41_negative, for example.
    And they should also be normalized to ai = 1, for all ais involved in that term.
    """

    aidecay = {
      "g2": 1.65684,
      "g4": 2.55052,
      "g1prime2": -12100.42,
      "ghzgs1prime2": -7613.351302119843,
    }

    def __init__(self):
        self.anomalouscouplings = []
        self.turnoff = []
        super(HZZAnomalousCouplingsFromHistograms, self).__init__()

    def setPhysicsOptions(self, physOptions):
        if not any(po.startswith("sqrts=") for po in physOptions):
            physOptions = physOptions + ["sqrts=13"]
        for po in physOptions:
            if po.startswith("turnoff="):
                self.turnoff += po.replace("turnoff=", "").split(",")
                #po gets removed in super
        super(MultiSignalSpinZeroHiggs, self).setPhysicsOptions(physOptions)
        if self.sqrts != [13]:
            raise ValueError("HZZAnomalousCouplingsFromHistograms is set up for 13 TeV only")
        if self.scaledifferentsqrtsseparately:
            raise ValueError("HZZAnomalousCouplingsFromHistograms is not set up for scaledifferentsqrtsseparately")
        if not self.scalemuvfseparately:
            raise ValueError("HZZAnomalousCouplingsFromHistograms is not set up for scalemuvmuftogether")

    def processPhysicsOptions(self,physOptions):
        processed = []
        for po in physOptions:
            if po in ("fa3", "fa2", "fL1", "fL1Zg"):
                if po in self.anomalouscouplings: raise ValueError("Provided physOption "+po+" twice")
                self.anomalouscouplings.append(po)
                processed.append(po)

        self.anomalouscouplings.sort(key="fa3 fa2 fL1 fL1Zg".index)

        for po in physOptions[:]:
            for i, fai in enumerate(self.anomalouscouplings, start=1):
                ai = fai[1:]
                if re.match("(f|phi){}(fixed|notpoi|floating|aspoi)$".format(ai).lower(), po.lower()):
                    physOptions.append(po.replace(ai, "ai{}".format(i)))
                    processed.append(po)

        processed += super(HZZAnomalousCouplingsFromHistograms, self).processPhysicsOptions(physOptions)

        if not self.anomalouscouplings: raise ValueError("Have to provide an anomalous coupling as a physOption (fa3, fa2, fL1, fL1Zg)")
        return processed

    @property
    def numberoffais(self):
        return len(self.anomalouscouplings)

    def getPOIList(self):
        self.modelBuilder.doVar("RF[1.0,0,10]")
        self.modelBuilder.doVar("RV[1.0,0,10]")
        self.modelBuilder.doVar("R[1.0,0,10]")

        pois = super(HZZAnomalousCouplingsFromHistograms, self).getPOIList()

        if not self.modelBuilder.out.var("CMS_zz4l_fa1"):
            expr = "-".join(["1"] + ["abs(@{})".format(i) for i in xrange(self.numberoffais)])
            fais = ", ".join("CMS_zz4l_fai{}".format(i) for i in xrange(1, self.numberoffais+1))
            self.modelBuilder.doVar('expr::CMS_zz4l_fa1("{}", {})'.format(expr, fais))

        if not self.modelBuilder.out.var("g1"):
            self.modelBuilder.doVar('expr::g1("sqrt(@0)", CMS_zz4l_fa1)')

        self.modelBuilder.doVar('expr::killswitch("@0>=0", CMS_zz4l_fa1)')

        couplings = ["g1"]
        i = 0
        for fai, ai in ("fa3", "g4"), ("fa2", "g2"), ("fL1", "g1prime2"), ("fL1Zg", "ghzgs1prime2"):
            if fai not in self.anomalouscouplings: continue
            i += 1

            kwargs = {
              "i": i,
              "ai": ai,
              "aidecay": self.aidecay[ai],
            }
            self.modelBuilder.doVar('expr::{ai}("(@0>0 ? 1 : -1) * sqrt(abs(@0))*{aidecay}", CMS_zz4l_fai{i})'.format(**kwargs))
            couplings.append(ai)

        if self.scaledifferentsqrtsseparately: raise ValueError("HZZAnomalousCouplingsFromHistograms is not compatible with scaledifferentsqrtsseparately")

        for g in couplings:
            self.modelBuilder.doVar('expr::ffH_{g}2("@0*@1*@2*@3*@3", killswitch, R, RF, {g})'.format(g=g))
            self.modelBuilder.doVar('expr::VVH_{g}4("@0*@1*@2*@3*@3*@3*@3", killswitch, R, RV, {g})'.format(g=g))

        kwargs = {}
        for kwargs["signname"], kwargs["sign"] in ("positive", ""), ("negative", "-"):
            for kwargs["g1"], kwargs["g2"] in itertools.combinations(couplings, 2):
                self.modelBuilder.doVar('expr::ffH_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@3*@4", killswitch, R, RF, {g1}, {g2})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVH_{g1}1{g2}3_{signname}("{sign}@0*@1*@2*@3*@4*@4*@4", killswitch, R, RV, {g1}, {g2})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVH_{g1}2{g2}2_{signname}("{sign}@0*@1*@2*@3*@3*@4*@4", killswitch, R, RV, {g1}, {g2})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVH_{g1}3{g2}1_{signname}("{sign}@0*@1*@2*@3*@3*@3*@4", killswitch, R, RV, {g1}, {g2})'.format(**kwargs))

            for kwargs["g1"], kwargs["g2"], kwargs["g3"] in itertools.combinations(couplings, 3):
                self.modelBuilder.doVar('expr::VVH_{g1}1{g2}1{g3}2_{signname}("{sign}@0*@1*@2*@3*@4*@5*@5", killswitch, R, RV, {g1}, {g2}, {g3})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVH_{g1}1{g2}2{g3}1_{signname}("{sign}@0*@1*@2*@3*@4*@4*@5", killswitch, R, RV, {g1}, {g2}, {g3})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVH_{g1}2{g2}1{g3}1_{signname}("{sign}@0*@1*@2*@3*@3*@4*@5", killswitch, R, RV, {g1}, {g2}, {g3})'.format(**kwargs))

            for kwargs["g1"], kwargs["g2"], kwargs["g3"], kwargs["g4"] in itertools.combinations(couplings, 4):
                self.modelBuilder.doVar('expr::VVH_{g1}1{g2}1{g3}1{g4}1_{signname}("{sign}@0*@1*@2*@3*@4*@5*@6", killswitch, R, RV, {g1}, {g2}, {g3}, {g4})'.format(**kwargs))

        return pois

    def getYieldScale(self,bin,process):
        match = re.match("(gg|tt|bb|qq|Z|W)H_(0(?:PM|M|PH|L1|L1Zg)|((?:g(?:1|2|4|1prime2|hzgs1prime2)[1234])*)_(positive|negative))$", process)
        if not match:
            return super(HZZAnomalousCouplingsFromHistograms, self).getYieldScale(bin, process)

        if match.group(1)+"H" in self.turnoff: return 0

        if match.group(1) in ("gg", "tt", "bb"): maxpower = 2; production = "ffH"
        elif match.group(1) in ("qq", "Z", "W"): maxpower = 4; production = "VVH"

        if match.group(2) == "0PM": powerdict = {"g1": maxpower}; sign = None
        elif match.group(2) == "0PH": powerdict = {"g2": maxpower}; sign = None
        elif match.group(2) == "0M": powerdict = {"g4": maxpower}; sign = None
        elif match.group(2) == "0L1": powerdict = {"g1prime2": maxpower}; sign = None
        elif match.group(2) == "0L1Zg": powerdict = {"ghzgs1prime2": maxpower}; sign = None
        else:

            powerdict = {coupling: int(power) for coupling, power in re.findall("(g(?:1|2|4|1prime2|hzgs1prime2))([1234])", match.group(3))}

            if sum(powerdict.values()) != maxpower:
                raise ValueError("power dict doesn't add up properly!  Sum should be {}\n{}\n{}".format(maxpower, process, powerdict))

            powerdict = collections.OrderedDict(
                sorted(powerdict.iteritems(), key = lambda x: "g1 g4 g2 g1prime2 ghzgs1prime2".index(x[0]))
            )

            sign = match.group(4)

        result = production + "_" + "".join("{}{}".format(k, v) for k, v in powerdict.iteritems())
        if sign is not None: result += "_" + sign

        print "Process {0} will scale by {1}".format(process,result)

        return result

spinZeroHiggs = SpinZeroHiggs()
multiSignalSpinZeroHiggs = MultiSignalSpinZeroHiggs()
hzzAnomalousCouplingsFromHistograms = HZZAnomalousCouplingsFromHistograms()
