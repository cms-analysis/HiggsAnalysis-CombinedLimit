from HiggsAnalysis.CombinedLimit.PhysicsModel import PhysicsModel


class SMEFT_chg(PhysicsModel):
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        self.modelBuilder.doVar("A[39.54]")
        self.modelBuilder.doVar("B[245.32]")
        self.modelBuilder.out.var("A").setConstant(True)
        self.modelBuilder.out.var("B").setConstant(True)
        self.modelBuilder.doVar("chg[0,-1,1]")
        self.modelBuilder.factory_('expr::ggH_scaling_chg("1+@1*@0+@2*@0*@0", chg, A, B)')
        self.modelBuilder.doSet("POI", ",".join(["chg"]))

    def getYieldScale(self, bin, process):
        if process == "ggH":
            return "ggH_scaling_chg"
        else:
            return 1


smeft_chg_tutorial = SMEFT_chg()
