from HiggsAnalysis.CombinedLimit.PhysicsModel import PhysicsModel


class FractionModel(PhysicsModel):
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        self.modelBuilder.doVar("f[0,0,4]")
        self.modelBuilder.doVar("r[1,0,10]")
        self.modelBuilder.factory_('expr::scale_qqH("(1-@0)*@1", f,r)')
        self.modelBuilder.factory_('expr::scale_ggH("@0*@1", f,r)')
        self.modelBuilder.doSet("POI", ",".join(["f"]))
        self.modelBuilder.doSet("POI", ",".join(["r"]))

    def getYieldScale(self, bin, process):
        if process == "qqH":
            return "scale_qqH"
        elif process == "ggH":
            return "scale_ggH"
        else:
            return 1


Fraction_2signals = FractionModel()
