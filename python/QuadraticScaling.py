import numpy as np

from HiggsAnalysis.CombinedLimit.PhysicsModel import PhysicsModel


class Quadratic(PhysicsModel):
    """Apply process scaling due to EFT operators.

    This class takes a dictionary of quadratic fits describing how processes are
    scaled as a function of an EFT operator's Wilson coefficient and adds it to
    the workspace. For an example coefficient x, dictionary values should have
    the form `(a, b, c)` where `xsec_NP(x) / xsec_SM = a + bx + cx^2`.

    To produce an example dictionary, for coefficient `cuW`:
    >>> import numpy as np
    >>> scales = {'cuW': {'ttZ': (1, 0.322778, 653.371), 'ttW': (1, 1.20998, 205.528)}}
    >>> np.save('scales.npy', scales)

    Example for running:
    text2workspace.py ttV.txt -P HiggsAnalysis.CombinedLimit.EFT:quad --PO scaling=scales.npy --PO process=ttZ --PO process=ttW --PO coefficient=cuW -o cuW.root
    combine -M MultiDimFit cuW.root --setParameterRanges=-4,4
    """

    def setPhysicsOptions(self, options):
        self.coefficient = None
        self.processes = []
        for option, value in [x.split('=') for x in options]:
            if option == 'coefficient':
                if self.coefficient is not None:
                    raise NotImplementedError('only one coefficient currently supported')
                self.coefficient = value
            if option == 'process':  # processes which will be scaled
                self.processes.append(value)
            if option == 'scaling':
                self.scaling = value

    def setup(self):
        scaling = np.load(self.scaling)[()]
        for process in self.processes:
            self.modelBuilder.out.var(process)
            name = 'r_{0}_{1}'.format(process, self.coefficient)
            if not self.modelBuilder.out.function(name):
                template = "expr::{name}('{a0} + ({a1} * {c}) + ({a2} * {c} * {c})', {c})"
                a0, a1, a2 = scaling[self.coefficient][process]
                quadratic = self.modelBuilder.factory_(template.format(name=name, a0=a0, a1=a1, a2=a2, c=self.coefficient))
                self.modelBuilder.out._import(quadratic)

    def doParametersOfInterest(self):
        # user should call combine with `--setPhysicsModelParameterRanges` set to sensible ranges
        self.modelBuilder.doVar('{0}[0, -inf, inf]'.format(self.coefficient))
        self.modelBuilder.doSet('POI', self.coefficient)
        self.setup()

    def getYieldScale(self, bin, process):
        if process not in self.processes:
            return 1
        else:
            name = 'r_{0}_{1}'.format(process, self.coefficient)

            return name


quad = Quadratic()
