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
          setattr(self, key, t(eval(value)))
       else:
          raise ValueError(
             "Model option %s not recognised!"
             " Available options are: %s" %(key, ','.join(self.opts))
             )
       

class CTagEfficiency(PhysicsModel):
    def __init__(self):
        PhysicsModel.__init__(self)
        self.opts = PhysOpts()
        self.opts.add('fitLightEff', True)
        self.opts.add('verbose', False)
        self.opts.add('inclusive', False)
        self.paganini = True #avoid repetitions
        self.exprs = {}
        self.pars  = ['signal_norm', 'LCharmE', 'SLightE', 'lead_cfrac', 'LLightE', 'SCharmE', 'sub_cfrac']
        self.categories = set(['notag', 'leadtag', 'subtag', 'ditag'])

    def setPhysicsOptions(self,physOptions):
        '''Receive a list of strings with the physics options from command line'''
        for po in physOptions:
           self.opts.parse(po)
        if self.opts.inclusive:
            self.categories = set(['Inc_nolead', 'Inc_nosub', 'Inc_leadtag', 'Inc_subtag'])

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        #tt signal strenght 0-200% on over-all right combination ttbar scaling
        #self.modelBuilder.doVar('strength[4347,0,8000]') 
        #what we actually want to measure
        self.modelBuilder.doVar('charmSF[1,0,2]')
        if self.opts.fitLightEff:
            print 'PASSING HERE', self.opts.fitLightEff
            self.modelBuilder.doVar('lightSF[1,0,2]')
            self.modelBuilder.doSet('POI','charmSF,lightSF')
        else:
            self.modelBuilder.doSet('POI','charmSF')

        self.exprs = {
            'notag'   : 'expr::Scaling_notag(  "signal_norm*((1-LCharmE)*(1-SLightE)*lead_cfrac+(1-LLightE)*(1-SCharmE)*sub_cfrac+(1-LLightE)*(1-SLightE)*(1-lead_cfrac-sub_cfrac))", {PARS})',
            'leadtag' : 'expr::Scaling_leadtag("signal_norm*(LCharmE*(1-SLightE)*lead_cfrac+LLightE*(1-SCharmE)*sub_cfrac+LLightE*(1-SLightE)*(1-lead_cfrac-sub_cfrac))", {PARS})',
            'subtag'  : 'expr::Scaling_subtag( "signal_norm*((1-LCharmE)*SLightE*lead_cfrac+(1-LLightE)*SCharmE*sub_cfrac+(1-LLightE)*SLightE*(1-lead_cfrac-sub_cfrac))", {PARS})',
            'ditag'   : 'expr::Scaling_ditag(  "signal_norm*(LCharmE*SLightE*lead_cfrac+LLightE*SCharmE*sub_cfrac+LLightE*SLightE*(1-lead_cfrac-sub_cfrac))", {PARS})',
            ############################
            ##  INCLUSIVE CATEGORIES  ##
            ############################
            'Inc_nolead'  : 'expr::Scaling_Inc_nolead( "signal_norm*((1-LCharmE)*lead_cfrac+(1-LLightE)*(1-lead_cfrac))", {PARS})',
            'Inc_nosub'   : 'expr::Scaling_Inc_nosub(  "signal_norm*((1-SCharmE)*sub_cfrac+(1-SLightE)*(1-sub_cfrac))", {PARS})',
            'Inc_leadtag' : 'expr::Scaling_Inc_leadtag("signal_norm*(LCharmE*lead_cfrac+LLightE*(1-lead_cfrac))", {PARS})',
            'Inc_subtag'  : 'expr::Scaling_Inc_subtag( "signal_norm*(SCharmE*sub_cfrac+SLightE*(1-sub_cfrac))", {PARS})',
            }

    @staticmethod
    def replace_pars(expr, pars):
        idx = 0
        out_pars = []
        for par in pars:
            if par in expr:
                expr = expr.replace(par, '@%i' % idx)
                out_pars.append(par)
                idx += 1
        return expr, out_pars

    def getYieldScale(self,bin,process):
        #check category (bin) consistency
        if bin not in self.categories:
            raise ValueError(
                'Asking information for bin %s, but allowed bins are %s, '
                'did you forget the inclusive=True option?' % (bin, self.categories))

        if self.paganini:
            #make products of MC efficiency * SF to make category yields less painful
            #L(eading)/S(ubleading) Charm/Light E(fficiency)
            self.modelBuilder.factory_('expr::LCharmE("@0*@1", charmSF, mc_lead_charm_eff)')
            self.modelBuilder.factory_('expr::SCharmE("@0*@1", charmSF, mc_sub_charm_eff )')
            self.modelBuilder.factory_('expr::LLightE("@0*@1", lightSF, mc_lead_light_eff)') 
            self.modelBuilder.factory_('expr::SLightE("@0*@1", lightSF, mc_sub_light_eff )') 

        if self.DC.isSignal[process]:
            expr = self.exprs[bin]
            ## for idx, name in enumerate(self.pars):
            ##     expr = expr.replace(name, '@%i' % idx)
            ## expr = expr.format(PARS=(', '.join(self.pars)))
            expr, pars = CTagEfficiency.replace_pars(expr, self.pars)
            expr = expr.format(PARS=(', '.join(pars)))
            ## set_trace()
            try:
                self.modelBuilder.factory_(expr)
            except RuntimeError as e:
                if str(e) == 'Error in factory statement':
                    raise RuntimeError('%s: %s' % (str(e), expr))
                else:
                    raise e

        if self.paganini:
            if self.opts.verbose:
                self.modelBuilder.out.Print()
            self.paganini = False

        return 'Scaling_%s' % bin if self.DC.isSignal[process] else 1 


ctagEfficiency = CTagEfficiency()

