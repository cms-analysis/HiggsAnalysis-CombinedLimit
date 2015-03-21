import math
from HiggsAnalysis.CombinedLimit.PhysicsModel import *

### This is the base python class to study the SpinZeroHiggsCTau

class SpinZeroHiggsCTau(PhysicsModel):
    def __init__(self):
        self.mHRange = []

        self.muFloating = True
        self.muAsPOI = False

        self.ctauFloating = True
        self.ctauPOI = True

        self.poiMap = []
        self.pois = {}
        self.verbose = False

    def setModelBuilder(self, modelBuilder):
        PhysicsModel.setModelBuilder(self,modelBuilder)
        self.modelBuilder.doModelBOnly = False

    def getYieldScale(self,bin,process):
        "Split in production and decay, and call getHiggsSignalYieldScale; return 1 for backgrounds "
        #print "Bin ",bin
        #print "Process ",process
        if self.DC.isSignal[process]:
            self.my_norm = "r"

            print "Process {0} will scale by {1}".format(process,self.my_norm)
            return self.my_norm
        
        elif not self.DC.isSignal[process]: return 1
            

    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if 'muFixed' in po: 
                print "Will consider the signal strength as a fixed parameter"
                self.muFloating = False
            elif 'muAsPOI' in po: 
                print "Will consider the signal strength as a parameter of interest"
                self.muAsPOI = True

            if 'ctauFixed' in po or 'ctauNotPOI' in po: 
                print "CMS_zz4l_ctau is NOT A POI"
                self.ctauPOI = False
                if 'ctauFixed' in po:
                    print "Will fix CMS_zz4l_ctau to 0"
                    self.ctauFloating = False

            if not self.muAsPOI and not self.ctauPOI:
                print "No POIs detected: Switching to default configuration: Floating nuisance mu, floating POI ctau, eveything else fixed"
                self.muFloating = True
                self.muAsPOI = False
                self.ctauFloating = True
                self.ctauPOI = True
            
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        if self.ctauFloating:
            if self.modelBuilder.out.var("CMS_zz4l_ctau"):
                print "CTau variable is present, no need to rebuild"
                self.modelBuilder.out.var("CMS_zz4l_ctau").setVal(0)
            else:
                print "CTau variable is NOT present, rebuilding"
                self.modelBuilder.doVar("CMS_zz4l_ctau[0,0,1000]")
                self.modelBuilder.out.var("CMS_zz4l_ctau").setBins(50)
            print "Floating CMS_zz4l_ctau"
            if self.ctauPOI:
                poi = "CMS_zz4l_ctau"
            else:
                self.modelBuilder.out.var("CMS_zz4l_ctau").setAttribute("flatParam")

        else:
            if self.modelBuilder.out.var("CMS_zz4l_ctau"):
                self.modelBuilder.out.var("CMS_zz4l_ctau").setVal(0)
                self.modelBuilder.out.var("CMS_zz4l_ctau").setConstant()
            else:
                self.modelBuilder.doVar("CMS_zz4l_ctau[0]")
            print "Fixing CMS_zz4l_ctau"
                
        if self.muFloating:
            if self.modelBuilder.out.var("r"):
                self.modelBuilder.out.var("r").setRange(0.,40.)
                self.modelBuilder.out.var("r").setVal(1)
            else:
                self.modelBuilder.doVar("r[1,0,40]")

            if self.muAsPOI:
                print "Treating r as a POI"
                if self.ctauPOI:
                    poi += ",r"
                else:
                    poi = "r"
            else:
                self.modelBuilder.out.var("r").setAttribute("flatParam")
        else:
            if self.modelBuilder.out.var("r"):
                self.modelBuilder.out.var("r").setVal(1)
                self.modelBuilder.out.var("r").setConstant()
            else:
                self.modelBuilder.doVar("r[1]")

        if self.modelBuilder.out.var("CMS_zz4l_prod"):
            print "Found production systematic"
            self.modelBuilder.out.var("CMS_zz4l_prod").setAttribute("flatParam")

        self.modelBuilder.doSet("POI",poi)
        
spinZeroHiggsCTau = SpinZeroHiggsCTau()
