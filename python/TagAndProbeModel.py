from HiggsAnalysis.CombinedLimit.PhysicsModel import *
import re

class TagAndProbe(PhysicsModel):
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        self.modelBuilder.doVar("SF[1.0,0.0,2.0]");
        self.modelBuilder.doSet("POI","SF")
        if self.options.mass != 0:
            if self.modelBuilder.out.var("MH"):
              self.modelBuilder.out.var("MH").removeRange()
              self.modelBuilder.out.var("MH").setVal(self.options.mass)
            else:
              self.modelBuilder.doVar("MH[%g]" % self.options.mass);
        exp_pass=1
        exp_fail=1
        for b in self.DC.bins:
           for p in self.DC.exp[b].keys():
              if self.DC.isSignal[p]:
                 if re.search("pass",b):
                    exp_pass=self.DC.exp[b][p]
                 if re.search("fail",b):
                    exp_fail=self.DC.exp[b][p]
        self.modelBuilder.factory_('expr::fail_scale("(%f+%f-(%f*@0))/%f", SF)'%(exp_pass,exp_fail,exp_pass,exp_fail))
    def getYieldScale(self,bin,process):
        "Return the name of a RooAbsReal to scale this yield by or the two special values 1 and 0 (don't scale, and set to zero)"
        if self.DC.isSignal[process]:
                if re.search("pass",bin):
                        return "SF"
                elif re.search("fail",bin):
                        return "fail_scale"
        return 1

tagAndProbe = TagAndProbe()
