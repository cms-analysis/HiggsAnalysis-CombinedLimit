from HiggsAnalysis.CombinedLimit.PhysicsModel import *

### This is the base python class to study the SpinZero structure

class SpinZeroHiggs(PhysicsModel):
    def __init__(self):
        self.mHRange = []
        self.muAsPOI = False
        self.muFloating = False
        self.fg4fixed = False
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
            if self.muFloating:
                self.my_norm = "r"
            else:
                self.my_norm = 1
        
            print "Process {0} will scale by {1}".format(process,self.my_norm)
            return self.my_norm
        
        elif not self.DC.isSignal[process]: return 1
            

    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po == "muAsPOI":
                print "Will consider the signal strength as a parameter of interest"
                self.muAsPOI = True
                self.muFloating = True
            if po == "muFloating":
                print "Will consider the signal strength as a floating parameter (as a paraeter of interest if --PO muAsPOI is specified, as a nuisance otherwise)"
                self.muFloating = True
            if po == "fg4fixed":
                print "Will fix CMS_zz4l_fg4 to 0 and make signal strength as a parameter of interest"
                self.muAsPOI = True
                self.muFloating = True
                self.fg4fixed = True
            
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        if self.fg4fixed:
#            self.modelBuilder.doVar("CMS_zz4l_Gamma[1]")
            self.modelBuilder.doVar("CMS_zz4l_fg4[0]")
            self.modelBuilder.doVar("r[1,0,4]")
            print "Fixing CMS_zz4l_fg4"
            poi = "r"
        else:
            #self.modelBuilder.doVar("CMS_zz4l_Gamma[1.,0.001,25]")
            #poi = "CMS_zz4l_Gamma"
            self.modelBuilder.doVar("CMS_zz4l_fg4[-0.,-1,1]")
            poi = "CMS_zz4l_fg4"

            if self.muFloating:
                self.modelBuilder.doVar("r[1,0,4]")
                if self.muAsPOI:
                    print "Treating r as a POI"
                    poi += ",r"
                else:
                    self.modelBuilder.out.var("r").setAttribute("flatParam")
        
        self.modelBuilder.doSet("POI",poi)
        
spinZeroHiggs = SpinZeroHiggs()
