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
            if po.lower()=='fai1fixed' or po.lower()=='fai1notpoi':
                if not self.fai1POI: raise ValueError("Specified fai1Fixed and/or fai1NotPOI multiple times!\n{}".format(physOptions))
                print "CMS_zz4l_fai1 is NOT A POI"
                self.fai1POI = False
                if po.lower()=='fai1fixed':
                    print "Will fix CMS_zz4l_fai1 to 0"
                    self.fai1Floating = False
                processed.append(po)

            if po.lower()=='phiai1floating' or po.lower()=='phiai1aspoi':
                if self.phiai1Floating: raise ValueError("Specified phiai1Floating and/or phiai1AsPOI multiple times!\n{}".format(physOptions))
                print "Will consider phase ai1 as a floating parameter"
                self.phiai1Floating = True
                if po.lower()=='phiai1aspoi':
                    print "Will consider phase ai1 as a parameter of interest"
                    self.phiai1POI = True
                processed.append(po)

            if po.lower() == 'fai2floating' or po.lower() == 'fai2aspoi':
                if self.fai2Floating: raise ValueError("Specified fai2Floating and/or fai2AsPOI multiple times!\n{}".format(physOptions))
                print "Will float CMS_zz4l_fai2"
                self.fai2Floating = True
                if po.lower() == 'fai2aspoi':
                    self.fai2POI = True
                processed.append(po)

            if po.lower() == 'phiai2floating' or po.lower() == 'phiai2aspoi':
                if self.phiai2Floating: raise ValueError("Specified phiai2Floating and/or phiai2AsPOI multiple times!\n{}".format(physOptions))
                print "Will consider phase ai2 as a floating parameter"
                self.phiai2Floating = True
                if po.lower() == 'phiai2aspoi':
                    print "Will consider phase ai2 as a parameter of interest"
                    self.phiai2POI = True
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
            self.modelBuilder.out.var("CMS_zz4l_fai1").setConstant(False)
            if self.allowPMF:
                self.modelBuilder.out.var("CMS_zz4l_fai1").setRange(-1.,1.)
                print "Allowing negative CMS_zz4l_fai1"
            if self.fai1POI:
                poi.append("CMS_zz4l_fai1")
            else:
                self.modelBuilder.out.var("CMS_zz4l_fai1").setAttribute("flatParam")
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
            self.modelBuilder.out.var("CMS_zz4l_fai2").setConstant(False)
            if self.allowPMF:
                self.modelBuilder.out.var("CMS_zz4l_fai2").setRange(-1.,1.)
                print "Allowing negative CMS_zz4l_fai2"
            if self.fai2POI:
                poi.append("CMS_zz4l_fai2")
            else:
                self.modelBuilder.out.var("CMS_zz4l_fai2").setAttribute("flatParam")
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
            self.modelBuilder.out.var("CMS_zz4l_phiai1").setConstant(False)
            if self.phiai1POI:
                print "Treating phiai1 as a POI"
                poi.append("CMS_zz4l_phiai1")
            else:
                self.modelBuilder.out.var("CMS_zz4l_phiai1").setAttribute("flatParam")
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
            self.modelBuilder.out.var("CMS_zz4l_phiai2").setConstant(False)
            if self.phiai2POI:
                print "Treating phiai2 as a POI"
                poi.append("CMS_zz4l_phiai2")
            else:
                self.modelBuilder.out.var("CMS_zz4l_phiai2").setAttribute("flatParam")
        else:
            if self.modelBuilder.out.var("CMS_zz4l_phiai2"):
                self.modelBuilder.out.var("CMS_zz4l_phiai2").setVal(0)
                self.modelBuilder.out.var("CMS_zz4l_phiai2").setConstant()
            else:
                self.modelBuilder.doVar("CMS_zz4l_phiai2[0]")
            print "Fixing CMS_zz4l_phiai2"

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
        #not doing muAsPOI or fixMu, there are too many permutations.
        #should just be set when running combine.

    def setPhysicsOptions(self, physOptions):
        if not any(po.startswith("map=") for po in physOptions):
            #no po started with map --> no manual overriding --> use the defaults
            #can still override with e.g. turnoff=ZH,WH
            physOptions = ["map=.*/(gg|qq|Z|W|tt)H:1"] + physOptions
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

        return processed

    def getPOIList(self):
        result = super(MultiSignalSpinZeroHiggs, self).getPOIList()

        fixedorfloated = self.fixed+self.floated
        for variable in fixedorfloated:
            if not self.modelBuilder.out.var(variable):
                print "{} does not exist in the workspace!  Check:\n - your datacard maker\n - your sqrts option".format(variable)
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

spinZeroHiggs = SpinZeroHiggs()
multiSignalSpinZeroHiggs = MultiSignalSpinZeroHiggs()
