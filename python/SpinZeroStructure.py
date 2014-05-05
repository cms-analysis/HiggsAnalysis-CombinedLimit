import math
from HiggsAnalysis.CombinedLimit.PhysicsModel import *

### This is the base python class to study the SpinZero structure

class SpinZeroHiggs(PhysicsModel):
    def __init__(self):
        self.mHRange = []

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
            if 'muFixed' in po: 
                print "Will consider the signal strength as a fixed parameter"
                self.muFloating = False
            elif 'muAsPOI' in po: 
                print "Will consider the signal strength as a parameter of interest"
                self.muAsPOI = True

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
            
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        if self.fai1Floating:
            self.modelBuilder.doVar("CMS_zz4l_fai1[0.,0,1]")
            print "Floating CMS_zz4l_fai1"
            if self.fai1POI:
                poi = "CMS_zz4l_fai1"
        else:
            self.modelBuilder.doVar("CMS_zz4l_fai1[0]")
            print "Fixing CMS_zz4l_fai1"
                
        if self.fai2Floating:
            self.modelBuilder.doVar("CMS_zz4l_fai2[0.,0,1]")
            print "Floating CMS_zz4l_fai2"
            if self.fai2POI:
                if self.fai1POI:
                    poi += ",CMS_zz4l_fai2"
                else:
                    poi = "CMS_zz4l_fai2"
        else:
            self.modelBuilder.doVar("CMS_zz4l_fai2[0]")
            print "Fixing CMS_zz4l_fai2"

        if self.muFloating:
            self.modelBuilder.doVar("r[1,0,5]")
            if self.muAsPOI:
                print "Treating r as a POI"
                if self.fai1POI or self.fai2POI:
                    poi += ",r"
                else:
                    poi = "r"
            else:
                self.modelBuilder.out.var("r").setAttribute("flatParam")
        else:
            self.modelBuilder.doVar("r[1]")

        if self.phiai1Floating:
            self.modelBuilder.doVar("CMS_zz4l_phiai1[0.,-3.14159265359,3.14159265359]")
            print "Floating CMS_zz4l_phiai1"
            if self.phiai1POI:
                print "Treating phiai1 as a POI"
                if self.fai1POI or self.fai2POI or self.muAsPOI:
                    poi += ",CMS_zz4l_phiai1"
                else:
                    poi = "CMS_zz4l_phiai1"
        else:
            self.modelBuilder.doVar("CMS_zz4l_phiai1[0]")
            print "Fixing CMS_zz4l_phiai1"

        if self.phiai2Floating:
            self.modelBuilder.doVar("CMS_zz4l_phiai2[0.,-3.14159265359,3.14159265359]")
            print "Floating CMS_zz4l_phiai2"
            if self.phiai2POI:
                print "Treating phiai2 as a POI"
                if self.fai1POI or self.fai2POI or self.muAsPOI or self.phiai1POI:
                    poi += ",CMS_zz4l_phiai2"
                else:
                    poi = "CMS_zz4l_phiai2"
        else:
            self.modelBuilder.doVar("CMS_zz4l_phiai2[0]")
            print "Fixing CMS_zz4l_phiai2"
            
        self.modelBuilder.doSet("POI",poi)
        
spinZeroHiggs = SpinZeroHiggs()
