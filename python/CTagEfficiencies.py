from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from pdb import set_trace

class PhysOpts(object):
    def __init__(self):
       self.opts = set()
    def add(self, name, default=0.):
       setattr(self, name, default)
       self.opts.add(name)
    def parse(self, opt):
       key,value =tuple(opt.split('='))
       if key in self.opts:
          t = type(getattr(self, key))
          setattr(self, key, t(value))
       else:
          raise ValueError(
             "Model option %s not recognised!"
             " Available options are: %s" %(key, ','.join(self.opts))
             )
       

class CTagEfficiency(PhysicsModel):
    def __init__(self):
        PhysicsModel.__init__(self)
        self.opts = PhysOpts()
        self.opts.add('leadCharmFraction', -1.)
        self.opts.add('subCharmFraction' , -1.) 
        self.opts.add('mcCharmEff' , -1.) 
        self.opts.add('mcLightEff' , -1.) 
        self.opts.add('mcExpEvts' , -1.) 

    def setPhysicsOptions(self,physOptions):
        '''Receive a list of strings with the physics options from command line'''
        for po in physOptions:
           self.opts.parse(po)

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        #tt signal strenght 0-200% on over-all right combination ttbar scaling
        #self.modelBuilder.doVar('strength[4347,0,8000]') 
        #what we actually want to measure
        self.modelBuilder.doVar('charmTagScale[1,0,2]')
        self.modelBuilder.doVar('lightTagScale[1,0,2]')

        self.modelBuilder.doSet('POI','charmTagScale,lightTagScale')

        P_generic = '({ceff}*@0*{CFrac} + {leff}*@1*(1-{CFrac}))'
        P_lead = P_generic.format(
            CFrac=self.opts.leadCharmFraction,
            ceff =self.opts.mcCharmEff,
            leff =self.opts.mcLightEff 
            )
        P_sub = P_generic.format(
            CFrac=self.opts.subCharmFraction,
            ceff =self.opts.mcCharmEff,
            leff =self.opts.mcLightEff 
            )

        self.modelBuilder.factory_(
            'expr::Scaling_notag("{norm}*(1 -{F1} -{F2} +{F1}*{F2})",charmTagScale, lightTagScale)'.format(
                norm=self.opts.mcExpEvts,
                F1=P_lead,
                F2=P_sub
              ))
        self.modelBuilder.factory_(
            'expr::Scaling_leadtag("{norm}*{F1}*(1-{F2})",charmTagScale, lightTagScale)'.format(
                norm=self.opts.mcExpEvts,
                F1=P_lead,
                F2=P_sub
              ))
        self.modelBuilder.factory_(
            'expr::Scaling_subtag("{norm}*{F2}*(1-{F1})",charmTagScale, lightTagScale)'.format(
                norm=self.opts.mcExpEvts,
                F1=P_lead,
                F2=P_sub
              ))
        self.modelBuilder.factory_(
            'expr::Scaling_ditag("{norm}*{F2}*{F1}",charmTagScale, lightTagScale)'.format(
                norm=self.opts.mcExpEvts,
                F1=P_lead,
                F2=P_sub
              ))

        self.modelBuilder.out.Print()
        
    def getYieldScale(self,bin,process):
        if self.DC.isSignal[process]:
           return "Scaling_%s" % bin
        return 1


ctagEfficiency = CTagEfficiency()

