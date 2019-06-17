import collections, itertools, math
from HiggsAnalysis.CombinedLimit.PhysicsModel import *

### This is the base python class to study the SpinZero structure

class SpinZeroHiggsBase(PhysicsModelBase_NiceSubclasses):
    def __init__(self):
        super(SpinZeroHiggsBase, self).__init__()

        self.faiphiaistatus = collections.defaultdict(lambda: "fix")
        self.fairelative = collections.defaultdict(lambda: False)
        self.allowPMF = collections.defaultdict(lambda: False)
        self.allowPMF[0] = False

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

        if {self.faidefinitionorder(i) for i in xrange(self.numberoffais)} != set(xrange(self.numberoffais)):
            raise ValueError("faidefinitionorder is not defined right, should go from 0 to {0} for inputs from 0 to {0}\n{1}\n{2}".format(self.numberoffais-1, {self.faidefinitionorder(i) for i in xrange(self.numberoffais)}, set(xrange(self.numberoffais))))

        for po in physOptions:
            newpo = po.lower().replace("a1", "ai0")
            match = re.match("(f|phi)ai([0-9]+)(fixed|notpoi|floating|aspoi)((?:relative)?)$", newpo)
            if match:
                parametertype = match.group(1)
                i = int(match.group(2))
                whattodo = match.group(3)
                relative = bool(match.group(4))

                if not 0 <= i < self.numberoffais:
                    raise ValueError("There are only {} fais available, so can't do anything with {}".format(self.numberoffais-1, po))

                if relative:
                    if whattodo == "fixed": raise ValueError("relative doesn't make sense for fixed fais")
                    if parametertype == "phi": raise ValueError("relative doesn't make sense for phi")

                key = parametertype, i

                if key in self.faiphiaistatus:
                    raise ValueError("Specified multiple physics options for {}ai{}".format(parametertype, i).replace("ai0", "a1"))
                if self.faidefinitionorder(i) == self.numberoffais-1:
                    if parametertype == "phi":
                        raise ValueError("fai{} is the last parameter to be defined, so it doesn't have a phi.  Can't do anything with {}".format(i, po).replace("ai0", "a1"))

                if whattodo == "fixed":
                    self.faiphiaistatus[key] = "fix"
                    if parametertype == "f" and i == 1:
                        print "Will fix {} to 0".format(self.parametername(parametertype, i, relative))
                elif whattodo == "floating" or whattodo == "notpoi":
                    self.faiphiaistatus[key] = "float"
                    if parametertype == "f" and i == 1:
                        print "{} is NOT A POI".format(self.parametername(parametertype, i, relative))
                    else:
                        print "Will float {}".format(self.parametername(parametertype, i, relative))
                elif whattodo == "aspoi":
                    self.faiphiaistatus[key] = "POI"
                    print "Will consider {} as a parameter of interest".format(self.parametername(parametertype, i, relative))
                else:
                    assert False, whattodo

                if relative:
                    self.fairelative[key] = True

                processed.append(po)

            if po.lower() == 'allowpmf':
                self.allowPMF.default_factory = lambda: True
                processed.append(po)

            match = re.match("(allow|forbid)pmfai([0-9]+)", newpo)
            if match:
                i = int(match.group(2))
                toset = {"allow": True, "forbid": False}[match.group(1)]
                if i == 0:
                    if toset == True:
                        raise ValueError("fa1 has to be fixed to positive")
                else:
                    if i in self.allowPMF: raise ValueError("Two different options set for allow/forbidPMfai{}".format(i).replace("ai0", "a1"))
                    self.allowPMF[i] = toset
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
        return 3 #including fa1

    def faidefinitionorder(self, i):
        if i == 0: return self.numberoffais-1
        return i-1
    def faidefinitionorderinverse(self, i):
        for _ in xrange(self.numberoffais):
            if self.faidefinitionorder(_) == i:
                return _
        raise ValueError("faidefinitionorder doesn't have anything that gives {}".format(i))
    def parametername(self, parametertype, i, relative):
      if parametertype == "phi" and self.faidefinitionorder(i) == self.numberoffais-1: raise ValueError("Can't have phiai for the last fai")
      if i >= self.numberoffais: raise ValueError("Only have {} fais".format(self.numberoffais-1))
      if i < 0: raise ValueError("Only have positive integer fais, plus 0 for fa1")
      return "CMS_zz4l_{}ai{:d}{}".format(parametertype, i, "_relative" if relative else "").replace("ai0", "a1")

    def getPOIList(self):
        poi = []
        poi += super(SpinZeroHiggsBase, self).getPOIList()

        for i in xrange(self.numberoffais):
            varname = self.parametername("f", i, False)
            if self.faidefinitionorder(i) == self.numberoffais-1:
                if self.faiphiaistatus["f", i] == "POI":
                    if not self.allowPMF[i]:
                        raise ValueError("fai{} is the last parameter to be defined, so it's defined as a function of the others.  It's sign is fixed.  Can't set it as a POI")
                    self.faiphiaistatus["f", i] = "lastPOI"
                else:
                    self.faiphiaistatus["f", i] = "last"
        if ("f", 1) not in self.faiphiaistatus: self.faiphiaistatus["f", 1] = "POI"

        for parametertype in "f", "phi":
            for i in xrange(self.numberoffais):
                if self.faidefinitionorder(i) == self.numberoffais-1: continue

                varname = self.parametername(parametertype, i, self.fairelative[parametertype, i])

                if not self.modelBuilder.out.var(varname):
                    self.modelBuilder.doVar(varname+"[0,0,1]") #will set the range later

                self.modelBuilder.out.var(varname).setVal(1 if i==0 else 0) #set fa1 to 1, anomalous couplings to 0

        parametertype = "f"
        done = {i: False for i in xrange(self.numberoffais)}
        while not all(done.values()):
            for i in xrange(self.numberoffais):
                if self.faidefinitionorder(i) == self.numberoffais-1: done[i] = True
                if done[i]: continue
                if self.fairelative[parametertype, i]:
                    otheris = [j for j in xrange(self.numberoffais) if self.faidefinitionorder(j) < self.faidefinitionorder(i)]
                    if not all(done[j] for j in otheris): continue
                    expr = "-".join(["1"] + ["abs(@{})".format(k) for k, j in enumerate(otheris)])
                    expr = "(" + expr + ")" + " * @{}".format(len(otheris))
                    fais = ", ".join(
                      [self.parametername("f", j, False) for j in otheris]
                      + [self.parametername("f", i, True)]
                    )
                    self.modelBuilder.doVar('expr::{}("{}", {})'.format(self.parametername("f", i, False), expr, fais))
                done[i] = True


        for parametertype in "f", "phi":
            for i in xrange(self.numberoffais):
                if self.faidefinitionorder(i) == self.numberoffais-1 and parametertype == "phi": continue
                relative = self.fairelative[parametertype, i]
                varname = self.parametername(parametertype, i, relative)
                status = self.faiphiaistatus[parametertype, i]
                if status in ("float", "POI"):
                    if parametertype == "f":
                        if self.faidefinitionorder(i) == 0 or relative:
                            parameterrange = (-1 if self.allowPMF[i] else 0), 1
                        else:
                            expr = "-".join(["1"] + ["abs(@{})".format(k) for k, j in enumerate(j for j in xrange(self.numberoffais) if self.faidefinitionorder(j) < self.faidefinitionorder(i))])
                            fais = ", ".join(self.parametername("f", j, False) for j in xrange(self.numberoffais) if self.faidefinitionorder(j) < self.faidefinitionorder(i))
                            self.modelBuilder.doVar('expr::max_'+varname+'("{}", {})'.format(expr, fais))
                            if self.allowPMF[i]:
                                self.modelBuilder.doVar('expr::min_{0}("-@0", max_{0})'.format(varname))
                            else:
                                self.modelBuilder.doVar('expr::min_{0}("0")'.format(varname))
                            parameterrange = (
                                self.modelBuilder.out.obj("min_"+varname),
                                self.modelBuilder.out.obj("max_"+varname),
                            )
                    elif parametertype == "phi":
                        parameterrange = -math.pi, math.pi
                    else:
                        assert False
                    self.modelBuilder.out.var(varname).setRange(*parameterrange)
                    self.modelBuilder.out.var(varname).setConstant(False)
                    if status == "POI":
                        print "Treating "+varname+" as a POI"
                        poi.append(varname)
                    else:
                        print "Floating "+varname
                        self.modelBuilder.out.var(varname).setAttribute("flatParam")
                    if parametertype == "f" and self.allowPMF[i]: print "Allowing negative "+varname
                elif status == "fix":
                    self.modelBuilder.out.var(varname).setConstant()
                    print "Fixing "+varname
                elif status in ("last", "lastPOI"):
                    expr = "-".join(["1"] + ["abs(@{})".format(j) for j in xrange(self.numberoffais-1)])
                    fais = ", ".join(self.parametername(parametertype, j, relative) for j in xrange(self.numberoffais) if j != i)

                    if self.allowPMF[i]:
                        self.modelBuilder.doVar('sgn{}[negative,positive]'.format(varname))
                        expr = "({}) * 2 * (@{}-0.5)".format(expr, self.numberoffais-1)
                        fais += ", sgn"+varname
                        self.modelBuilder.out.cat("sgn"+varname).setAttribute("flatParam")
                        self.modelBuilder.addDiscrete("sgn"+varname)
                        print "Floating the sign of "+varname

                    self.modelBuilder.doVar('expr::{}("{}", {})'.format(varname, expr, fais))
                    print "Setting "+varname+" to 1 - the sum of the other fais"
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
            self.allowPMF = defaultdict(lambda: False)

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
        self.scalegL1by10000 = False
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
            if po in ("fa3", "fa2", "fL1", "fL1Zg", "fa1"):
                if po in self.anomalouscouplings: raise ValueError("Provided physOption "+po+" twice")
                self.anomalouscouplings.append(po)
                processed.append(po)
            if po == "scalegL1by10000":
                self.scalegL1by10000 = True
                processed.append(po)

        if "fa1" not in self.anomalouscouplings: self.anomalouscouplings.append("fa1")

        for po in physOptions[:]:
            for i, fai in enumerate(self.anomalouscouplings):
                if fai == "fa1": continue #handled in the base class
                ai = fai[1:]
                if re.match("(f|phi){}(fixed|notpoi|floating|aspoi)(?:relative)?$".format(ai).lower(), po.lower()):
                    physOptions.append(po.replace(ai, "ai{}".format(self.faidefinitionorderinverse(i))))
                    processed.append(po)

        processed += super(HZZAnomalousCouplingsFromHistograms, self).processPhysicsOptions(physOptions)

        if not self.anomalouscouplings: raise ValueError("Have to provide an anomalous coupling as a physOption (fa3, fa2, fL1, fL1Zg)")
        return processed

    @property
    def numberoffais(self):
        return len(self.anomalouscouplings)  #including fa1

    def faidefinitionorder(self, i):
        #CMS_zz4l_fai1, CMS_zz4l_fai2, etc. correspond to fa3, fa2, fL1, fL1Zg in that order
        #However they might not be defined in that order, e.g. CMS_zz4l_fai1 might be restricted to (0, 1-CMS_zz4l_fai2)
        sortedcouplings = sorted(self.anomalouscouplings, key=["fa1", "fa3", "fa2", "fL1", "fL1Zg"].index)
        return self.anomalouscouplings.index(sortedcouplings[i])

    def getPOIList(self):
        self.modelBuilder.doVar("RF[1.0,0,10]")
        self.modelBuilder.doVar("RV[1.0,0,10]")
        self.modelBuilder.doVar("R[1.0,0,10]")

        pois = super(HZZAnomalousCouplingsFromHistograms, self).getPOIList()

        if not self.modelBuilder.out.var("g1"):
            self.modelBuilder.doVar('expr::g1("sqrt(@0)", CMS_zz4l_fa1)')

        couplings = ["g1"]
        i = 0
        for fai, ai in ("fa3", "g4"), ("fa2", "g2"), ("fL1", "g1prime2"), ("fL1Zg", "ghzgs1prime2"):
            if fai not in self.anomalouscouplings: continue
            i += 1

            if self.scalegL1by10000:
                divideby = {
                    "g4": 1,
                    "g2": 1,
                    "g1prime2": 10000,
                    "ghzgs1prime2": 10000,
                }[ai]
            else:
                divideby = 1

            kwargs = {
              "i": i,
              "ai": ai,
              "aidecay": self.aidecay[ai] / divideby,
            }
            self.modelBuilder.doVar('expr::{ai}("(@0>0 ? 1 : -1) * sqrt(abs(@0))*{aidecay}", CMS_zz4l_fai{i})'.format(**kwargs))
            couplings.append(ai)

        if self.scaledifferentsqrtsseparately: raise ValueError("HZZAnomalousCouplingsFromHistograms is not compatible with scaledifferentsqrtsseparately")

        for g in couplings:
            self.modelBuilder.doVar('expr::ffH_{g}2("@0*@1*@2*@2", R, RF, {g})'.format(g=g))
            self.modelBuilder.doVar('expr::VVH_{g}4("@0*@1*@2*@2*@2*@2", R, RV, {g})'.format(g=g))

        kwargs = {}
        for kwargs["signname"], kwargs["sign"] in ("positive", ""), ("negative", "-"):
            for kwargs["g1"], kwargs["g2"] in itertools.combinations(couplings, 2):
                self.modelBuilder.doVar('expr::ffH_{g1}1{g2}1_{signname}("{sign}@0*@1*@2*@3", R, RF, {g1}, {g2})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVH_{g1}1{g2}3_{signname}("{sign}@0*@1*@2*@3*@3*@3", R, RV, {g1}, {g2})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVH_{g1}2{g2}2_{signname}("{sign}@0*@1*@2*@2*@3*@3", R, RV, {g1}, {g2})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVH_{g1}3{g2}1_{signname}("{sign}@0*@1*@2*@2*@2*@3", R, RV, {g1}, {g2})'.format(**kwargs))

            for kwargs["g1"], kwargs["g2"], kwargs["g3"] in itertools.combinations(couplings, 3):
                self.modelBuilder.doVar('expr::VVH_{g1}1{g2}1{g3}2_{signname}("{sign}@0*@1*@2*@3*@4*@4", R, RV, {g1}, {g2}, {g3})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVH_{g1}1{g2}2{g3}1_{signname}("{sign}@0*@1*@2*@3*@3*@4", R, RV, {g1}, {g2}, {g3})'.format(**kwargs))
                self.modelBuilder.doVar('expr::VVH_{g1}2{g2}1{g3}1_{signname}("{sign}@0*@1*@2*@2*@3*@4", R, RV, {g1}, {g2}, {g3})'.format(**kwargs))

            for kwargs["g1"], kwargs["g2"], kwargs["g3"], kwargs["g4"] in itertools.combinations(couplings, 4):
                self.modelBuilder.doVar('expr::VVH_{g1}1{g2}1{g3}1{g4}1_{signname}("{sign}@0*@1*@2*@3*@4*@5", R, RV, {g1}, {g2}, {g3}, {g4})'.format(**kwargs))

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
