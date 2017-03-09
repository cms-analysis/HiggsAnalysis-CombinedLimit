import math
from HiggsAnalysis.CombinedLimit.PhysicsModel import *

### This is the base python class to study the SpinZero structure

class SpinZeroHiggsBase(PhysicsModelBase_NiceSubclasses):
    def __init__(self):
        super(SpinZeroHiggsBase, self).__init__()

        self.fai1Floating = True
        self.fai1POI = True

        self.fai2Floating = False
        self.fai2POI = False

        self.allowPMF = False

        self.phiai1Floating = False
        self.phiai1POI = False

        self.phiai2Floating = False
        self.phiai2POI = False

        self.HWWcombination = False

    def setModelBuilder(self, modelBuilder):
        super(SpinZeroHiggsBase, self).setModelBuilder(modelBuilder)
        self.modelBuilder.doModelBOnly = False

    def getYieldScale(self,bin,process):
        "Split in production and decay, and call getHiggsSignalYieldScale; return 1 for backgrounds "
        result = super(SpinZeroHiggsBase, self).getYieldScale(bin, process)
        if result not in (0, 1): print "Process {0} will scale by {1}".format(process,result)
        return result

    def setPhysicsOptions(self,physOptions):
        super(SpinZeroHiggsBase, self).setPhysicsOptions(physOptions)
        for po in physOptions:
            if 'fai1Fixed' in po or 'fai1NotPOI' in po:
                print "CMS_zz4l_fai1 is NOT A POI"
                self.fai1POI = False
                if 'fai1Fixed' in po:
                    print "Will fix CMS_zz4l_fai1 to 0"
                    self.fai1Floating = False

            if 'phiai1Floating' in po or 'phiai1AsPOI' in po:
                print "Will consider phase ai1 as a floating parameter"
                self.phiai1Floating = True
                if 'phiai1AsPOI' in po:
                    print "Will consider phase ai1 as a parameter of interest"
                    self.phiai1POI = True

            if 'fai2Floating' in po or 'fai2AsPOI' in po:
                print "Will float CMS_zz4l_fai2"
                self.fai2Floating = True
                if 'fai2AsPOI' in po:
                    self.fai2POI = True

            if 'phiai2Floating' in po or 'phiai2AsPOI' in po:
                print "Will consider phase ai2 as a floating parameter"
                self.phiai2Floating = True
                if 'phiai2AsPOI' in po:
                    print "Will consider phase ai2 as a parameter of interest"
                    self.phiai2POI = True

            if 'allowPMF' in po:
                self.allowPMF = True

            if 'HWWcombination' in po:
                self.HWWcombination = True

    def getPOIList(self):

        poi = []
        poi += super(SpinZeroHiggsBase, self).getPOIList()

        if self.fai1Floating:
            if self.modelBuilder.out.var("CMS_zz4l_fai1"):
                self.modelBuilder.out.var("CMS_zz4l_fai1").setRange(0.,1.)
                self.modelBuilder.out.var("CMS_zz4l_fai1").setVal(0)
            else:
                self.modelBuilder.doVar("CMS_zz4l_fai1[0.,0,1]")
            print "Floating CMS_zz4l_fai1"
            if self.allowPMF:
                self.modelBuilder.out.var("CMS_zz4l_fai1").setRange(-1.,1.)
                print "Allowing negative CMS_zz4l_fai1"
            if self.fai1POI:
                poi.append("CMS_zz4l_fai1")
        else:
            if self.modelBuilder.out.var("CMS_zz4l_fai1"):
                self.modelBuilder.out.var("CMS_zz4l_fai1").setVal(0)
                self.modelBuilder.out.var("CMS_zz4l_fai1").setConstant()
            else:
                self.modelBuilder.doVar("CMS_zz4l_fai1[0]")
            print "Fixing CMS_zz4l_fai1"

        if self.fai2Floating:
            if self.modelBuilder.out.var("CMS_zz4l_fai2"):
                self.modelBuilder.out.var("CMS_zz4l_fai2").setRange(0.,1.)
                self.modelBuilder.out.var("CMS_zz4l_fai2").setVal(0)
            else:
                self.modelBuilder.doVar("CMS_zz4l_fai2[0.,0,1]")
            print "Floating CMS_zz4l_fai2"
            if self.allowPMF:
                self.modelBuilder.out.var("CMS_zz4l_fai2").setRange(-1.,1.)
                print "Allowing negative CMS_zz4l_fai2"
            if self.fai2POI:
                poi.append("CMS_zz4l_fai2")
        else:
            if self.modelBuilder.out.var("CMS_zz4l_fai2"):
                self.modelBuilder.out.var("CMS_zz4l_fai2").setVal(0)
                self.modelBuilder.out.var("CMS_zz4l_fai2").setConstant()
            else:
                self.modelBuilder.doVar("CMS_zz4l_fai2[0]")
            print "Fixing CMS_zz4l_fai2"

        if self.phiai1Floating:
            if self.modelBuilder.out.var("CMS_zz4l_phiai1"):
                self.modelBuilder.out.var("CMS_zz4l_phiai1").setRange(-3.14159265359,3.14159265359)
                self.modelBuilder.out.var("CMS_zz4l_phiai1").setVal(0)
            else:
                self.modelBuilder.doVar("CMS_zz4l_phiai1[0.,-3.14159265359,3.14159265359]")
            print "Floating CMS_zz4l_phiai1"
            if self.phiai1POI:
                print "Treating phiai1 as a POI"
                poi.append("CMS_zz4l_phiai1")
        else:
            if self.modelBuilder.out.var("CMS_zz4l_phiai1"):
                self.modelBuilder.out.var("CMS_zz4l_phiai1").setVal(0)
                self.modelBuilder.out.var("CMS_zz4l_phiai1").setConstant()
            else:
                self.modelBuilder.doVar("CMS_zz4l_phiai1[0]")
            print "Fixing CMS_zz4l_phiai1"

        if self.phiai2Floating:
            if self.modelBuilder.out.var("CMS_zz4l_phiai2"):
                self.modelBuilder.out.var("CMS_zz4l_phiai2").setRange(-3.14159265359,3.14159265359)
                self.modelBuilder.out.var("CMS_zz4l_phiai2").setVal(0)
            else:
                self.modelBuilder.doVar("CMS_zz4l_phiai2[0.,-3.14159265359,3.14159265359]")
            print "Floating CMS_zz4l_phiai2"
            if self.phiai2POI:
                print "Treating phiai2 as a POI"
                poi.append("CMS_zz4l_phiai2")
        else:
            if self.modelBuilder.out.var("CMS_zz4l_phiai2"):
                self.modelBuilder.out.var("CMS_zz4l_phiai2").setVal(0)
                self.modelBuilder.out.var("CMS_zz4l_phiai2").setConstant()
            else:
                self.modelBuilder.doVar("CMS_zz4l_phiai2[0]")
            print "Fixing CMS_zz4l_phiai2"

        if self.HWWcombination:
            if self.modelBuilder.out.var("CMS_zz4l_alpha"):
                self.modelBuilder.out.var("CMS_zz4l_alpha").setVal(0.5)
                print "Found CMS_zz4l_alpha; setting to 0.5"
            else:
                self.modelBuilder.doVar("CMS_zz4l_alpha[0.5,-1,1]")
                print "Creating CMS_zz4l_alpha; setting to 0.5"
        else:
            if self.modelBuilder.out.var("CMS_zz4l_alpha"):
                self.modelBuilder.out.var("CMS_zz4l_alpha").setVal(0)
                self.modelBuilder.out.var("CMS_zz4l_alpha").setConstant()
                print "Found CMS_zz4l_alpha; setting to constant 0"

        return poi

class SpinZeroHiggs(SpinZeroHiggsBase):
    def __init__(self):
        super(SpinZeroHiggs, self).__init__()
        self.muFloating = True
        self.muAsPOI = False

    def setPhysicsOptions(self,physOptions):
        super(SpinZeroHiggs, self).setPhysicsOptions(physOptions)
        for po in physOptions:
            if 'muFixed' in po:
                print "Will consider the signal strength as a fixed parameter"
                self.muFloating = False
            elif 'muAsPOI' in po:
                print "Will consider the signal strength as a parameter of interest"
                self.muAsPOI = True

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
        self.sqrts = None
        self.fix = []
        self.float = []
        #not doing muAsPOI or fixMu, there are too many permutations.
        #should just be set when running combine.

    def setPhysicsOptions(self, physOptions):
        if not any(po.startswith("map=") for po in physOptions):
            #no po started with map --> no manual overriding --> use the defaults
            #can still override with e.g. turnoff=ZH,WH
            for po in physOptions[:]:
                if po.lower() == "scalemuvmuftogether":
                    self.scalemuvfseparately = False
                    physOptions.remove(po)
                if po.lower() == "scaledifferentsqrtsseparately":
                    self.scaledifferentsqrtsseparately = True
                    physOptions.remove(po)
                if po.lower().startswith("sqrts="):
                    if self.sqrts is not None: raise ValueError("Duplicate physicsoption sqrts=?? provided")
                    self.sqrts = [int(_) for _ in po.replace("sqrts=", "").split(",")]
                    physOptions.remove(po)

            if self.sqrts is None:
                raise ValueError("PhysicsOption sqrts=?? is mandatory.  example: sqrts=7,8,13")

            physOptions = ["map=.*/(gg|qq|Z|W|tt)H:1"] + physOptions

            physOptions.sort(key=lambda x: x.startswith("verbose"), reverse=True) #put verbose at the beginning

            if self.scaledifferentsqrtsseparately and self.scalemuvfseparately:
                self.fix = ["RV", "RF", "R"] + ["R_{}TeV".format(_) for _ in self.sqrts]
                self.float = ["R{}_{}TeV".format(_1, _2) for _1 in ("V", "F") for _2 in self.sqrts]
            elif self.scaledifferentsqrtsseparately and not self.scalemuvfseparately:
                self.fix = ["RV", "RF", "R"] + ["R{}_{}TeV".format(_1, _2) for _1 in ("V", "F") for _2 in self.sqrts]
                self.float = ["R_{}TeV".format(_) for _ in self.sqrts]
            elif not scaledifferentsqrtsseparately and self.scalemuvfseparately:
                self.fix = ["R"] + ["R{}_{}TeV".format(_1, _2) for _1 in ("V", "F", "") for _2 in self.sqrts]
                self.float = ["RV", "RF"]
            elif not scaledifferentsqrtsseparately and not self.scalemuvfseparately:
                self.fix = ["RV", "RF"] + ["R{}_{}TeV".format(_1, _2) for _1 in ("V", "F", "") for _2 in self.sqrts]
                self.float = ["R"]
            else:
                assert False, "?????"

        super(MultiSignalSpinZeroHiggs, self).setPhysicsOptions(physOptions)

    def getPOIList(self):
        result = super(SpinZeroHiggs, self).getPOIList()

        for variable in self.fix:
            if not self.out.var(variable):
                raise ValueError("{} does not exist in the workspace!  Check:\n - your datacard maker\n - your sqrts option")
            self.out.var(variable).setVal(1)
            self.out.var(variable).setConstant()
        for variable in self.float:
            if not self.out.var(variable):
                raise ValueError("{} does not exist in the workspace!  Check:\n - your datacard maker\n - your sqrts option")
            self.out.var(variable).setAttribute("flatParam")

        return result

    def getYieldScale(self,bin,process):
        result = super(MultiSignalSpinZeroHiggs, self).getYieldScale(bin,process)

        w = self.modelBuilder.out

        if w.function("muVratio") and w.var("r_VVH") and not w.function("muV_scaled"):
            self.modelBuilder.doExp("muV_scaled", "@0/@1", "r_VVH, muVratio")

        if w.function("mufratio") and w.var("r_ffH") and not w.function("muf_scaled"):
            self.modelBuilder.doExp("muf_scaled", "@0/@1", "r_ffH, mufratio")

        return result

spinZeroHiggs = SpinZeroHiggs()
multiSignalSpinZeroHiggs = MultiSignalSpinZeroHiggs()
