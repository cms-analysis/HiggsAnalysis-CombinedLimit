from math import *
from array import array
import os 
import ROOT
from HiggsAnalysis.CombinedLimit.PhysicsModel import ALL_HIGGS_DECAYS

class SMHiggsBuilder:
    def __init__(self,modelBuilder,datadir=None):
        self.modelBuilder = modelBuilder
        if datadir == None:
            datadir = os.environ['CMSSW_BASE']+"/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg"
        self.datadir = datadir
        self.brpath = os.path.join(self.datadir,'sm/br')
        self.coupPath = os.path.join(self.datadir,'couplings')

    def makeXS(self,process, energy='7TeV'):
        self.xspath = os.path.join(self.datadir, 'sm/xs', energy)
        if process == "ggH": self.textToSpline("SM_XS_ggH_"+energy, os.path.join(self.xspath, energy+"-ggH.txt") );
        if process == "qqH": self.textToSpline("SM_XS_qqH_"+energy, os.path.join(self.xspath, energy+"-vbfH.txt") );
        if process == "ttH": self.textToSpline("SM_XS_ttH_"+energy, os.path.join(self.xspath, energy+"-ttH.txt") );
        if process == "WH":  self.textToSpline("SM_XS_WH_"+energy,  os.path.join(self.xspath, energy+"-WH.txt") );
        if process == "ZH":  self.textToSpline("SM_XS_ZH_"+energy,  os.path.join(self.xspath, energy+"-ZH.txt") );
        if process == "bbH":  self.textToSpline("SM_XS_bbH_"+energy,  os.path.join(self.xspath, energy+"-bbH.txt") );
        if process == "tHq":  self.textToSpline("SM_XS_tHq_"+energy,  os.path.join(self.xspath, energy+"-tHq.txt") );
        if process == "tHW":  self.textToSpline("SM_XS_tHW_"+energy,  os.path.join(self.xspath, energy+"-tHW.txt") );
        if process == "VH":  
            makeXS("WH", energy); makeXS("ZH", energy);
            self.modelBuilder.factory_('sum::SM_XS_VH_'+energy+'(SM_XS_WH_'+energy+',SM_XS_ZH_'+energy+')')
    def makeTotalWidth(self):
        self.textToSpline("SM_GammaTot", os.path.join(self.brpath,"BR4.txt"), ycol=12);
    def makeBR(self,decay):
        if decay == "hww": self.textToSpline("SM_BR_hww", os.path.join(self.brpath, "BR4.txt"), ycol=10);
        if decay == "hzz": self.textToSpline("SM_BR_hzz", os.path.join(self.brpath, "BR4.txt"), ycol=11);
        if decay == "hgg": self.textToSpline("SM_BR_hgg", os.path.join(self.brpath, "BR4.txt"), ycol=8);
        if decay == "hzg": self.textToSpline("SM_BR_hzg", os.path.join(self.brpath, "BR4.txt"), ycol=9);
        if decay == "hbb": self.textToSpline("SM_BR_hbb", os.path.join(self.brpath, "BR4.txt"), ycol=1);
        if decay == "htt": self.textToSpline("SM_BR_htt", os.path.join(self.brpath, "BR4.txt"), ycol=2);
        if decay == "hmm": self.textToSpline("SM_BR_hmm", os.path.join(self.brpath, "BR4.txt"), ycol=3);
        if decay == "hss": self.textToSpline("SM_BR_hss", os.path.join(self.brpath, "BR4.txt"), ycol=5);
        if decay == "hcc": self.textToSpline("SM_BR_hcc", os.path.join(self.brpath, "BR4.txt"), ycol=4);
        if decay == "hgluglu": self.textToSpline("SM_BR_hgluglu", os.path.join(self.brpath, "BR4.txt"), ycol=7);
        if decay == "htoptop": self.textToSpline("SM_BR_htoptop", os.path.join(self.brpath, "BR4.txt"), ycol=6);
    def makePartialWidth(self,decay):
        self.makeTotalWidth(); 
        self.makeBR(decay);
        self.modelBuilder.factory_('prod::SM_Gamma_%s(SM_GammaTot,SM_BR_%s)' % (decay,decay))
    def makeScaling(self,what, Cb='Cb', Ctop='Ctop', CW='CW', CZ='CZ', Ctau='Ctau', Cc='Ctop', suffix=''):
        prefix = 'SM_%(what)s_' % locals()
        if suffix:
            suffix += '_'
            prefix += suffix
#        self.modelBuilder.doVar('One[1]')
#        self.modelBuilder.doVar('Zero[0]') 
        if what.startswith('qqH'):
            for sqrts in ('7TeV', '8TeV','13TeV','14TeV'):
                rooName = prefix+'RVBF_'+sqrts
                self.textToSpline(rooName, os.path.join(self.coupPath, 'R_VBF_%(sqrts)s.txt'%locals()), ycol=1 )
                scalingName = 'Scaling_'+what+'_'+sqrts
#                print 'Building '+ scalingName
                rooExpr = 'expr::%(scalingName)s(\
"(@0*@0 + @1*@1 * @2 )/(1+@2)",\
 %(CW)s, %(CZ)s,\
 %(rooName)s\
)'%locals()
#                print  rooExpr
                self.modelBuilder.factory_(rooExpr)
        elif what.startswith('ggH'):
            structure = {'c_kt2':1, 'c_kb2':2, 'c_ktkb':3, 'c_ktkc':4, 'c_kbkc':5, 'c_kc2':6}
            for sqrts in ('7TeV', '8TeV', '13TeV', '14TeV'):
                for qty, column in structure.iteritems():
                    rooName = prefix+qty+'_'+sqrts
                    self.textToSpline(rooName, os.path.join(self.coupPath, 'ggH_%(sqrts)s.txt'%locals()), ycol=column )
                scalingName = 'Scaling_'+what+'_'+sqrts
#                print 'Building '+scalingName
		coeffSum = 'expr::coeff_sum_%(scalingName)s(\
"@0+@1+@2+@3+@4+@5",\
 %(prefix)sc_kt2_%(sqrts)s, %(prefix)sc_kb2_%(sqrts)s, %(prefix)sc_ktkb_%(sqrts)s, %(prefix)sc_ktkc_%(sqrts)s, %(prefix)sc_kbkc_%(sqrts)s, %(prefix)sc_kc2_%(sqrts)s\
)'%locals()
                self.modelBuilder.factory_(coeffSum)
                rooExpr = 'expr::%(scalingName)s(\
"((@0*@0)*@3  + (@1*@1)*@4 + (@0*@1)*@5 + (@0*@2)*@6 + (@1*@2)*@7 + (@2*@2)*@8)/@9 ",\
 %(Ctop)s, %(Cb)s, %(Cc)s,\
 %(prefix)sc_kt2_%(sqrts)s, %(prefix)sc_kb2_%(sqrts)s, %(prefix)sc_ktkb_%(sqrts)s, %(prefix)sc_ktkc_%(sqrts)s, %(prefix)sc_kbkc_%(sqrts)s, %(prefix)sc_kc2_%(sqrts)s, coeff_sum_%(scalingName)s\
)'%locals()
#                print  rooExpr
                self.modelBuilder.factory_(rooExpr)
        elif what.startswith('hgluglu'):
            structure = {'Gamma_tt':2, 'Gamma_bb':3, 'Gamma_tb':4}
            for qty, column in structure.iteritems():
                rooName = prefix+qty
                self.textToSpline(rooName, os.path.join(self.coupPath, 'Gamma_Hgluongluon.txt'), ycol=column )
            scalingName = 'Scaling_'+what
#            print 'Building '+scalingName
            rooExpr = 'expr::%(scalingName)s(\
"(@0*@0)*@2  + (@1*@1)*@3 + (@0*@1)*@4",\
 %(Ctop)s, %(Cb)s,\
 %(prefix)sGamma_tt, %(prefix)sGamma_bb, %(prefix)sGamma_tb\
)'%locals()
#            print  rooExpr
            self.modelBuilder.factory_(rooExpr)
        elif what.startswith('hgg') or what.startswith('hzg'): #in ['hgg', 'hzg']:
            fileFor = {'hgg':'Gamma_Hgammagamma.txt',
                       'hzg':'Gamma_HZgamma.txt'}
            structure = {'Gamma_tt':2, 'Gamma_bb':3, 'Gamma_WW':4,
                         'Gamma_tb':5, 'Gamma_tW':6, 'Gamma_bW':7,
                         'Gamma_ll':8,
                         'Gamma_tl':9, 'Gamma_bl':10, 'Gamma_lW':11}
            for qty, column in structure.iteritems():
                rooName = prefix+qty
                self.textToSpline(rooName, os.path.join(self.coupPath, fileFor['hgg' if what.startswith('hgg') else 'hzg']), ycol=column )
            scalingName = 'Scaling_'+what
#            print 'Building '+scalingName
            rooExpr = 'expr::%(scalingName)s(\
"( (@0*@0)*@4 + (@1*@1)*@5 + (@2*@2)*@6 + (@0*@1)*@7 + (@0*@2)*@8 + (@1*@2)*@9 + (@3*@3)*@10 + (@0*@3)*@11 + (@1*@3)*@12 + (@2*@3)*@13 ) / (@4+@5+@6+@7+@8+@9+@10+@11+@12+@13)",\
 %(Ctop)s, %(Cb)s, %(CW)s, %(Ctau)s,\
 %(prefix)sGamma_tt, %(prefix)sGamma_bb, %(prefix)sGamma_WW,\
 %(prefix)sGamma_tb, %(prefix)sGamma_tW, %(prefix)sGamma_bW,\
 %(prefix)sGamma_ll,\
 %(prefix)sGamma_tl, %(prefix)sGamma_bl, %(prefix)sGamma_lW\
)'%locals()
#            print  rooExpr
            self.modelBuilder.factory_(rooExpr)
        elif what.startswith('ggZH'):
            structure = {'c_kt2':1, 'c_kb2':2, 'c_kZ2':3, 'c_ktkb':4, 'c_ktkZ':5, 'c_kbkZ':6}
            for sqrts in ('7TeV','8TeV','13TeV','14TeV'):
                for qty, column in structure.iteritems():
                    rooName = prefix+qty+'_'+sqrts
                    self.textToSpline(rooName, os.path.join(self.coupPath, 'ggZH_%(sqrts)s.txt'%locals()), ycol=column )
                scalingName = 'Scaling_'+what+'_'+sqrts
                coeffSum = 'expr::coeff_sum_%(scalingName)s( "@0+@1+@2+@3+@4+@5",\
		%(prefix)sc_kt2_%(sqrts)s, %(prefix)sc_kb2_%(sqrts)s, %(prefix)sc_kZ2_%(sqrts)s, %(prefix)sc_ktkb_%(sqrts)s, %(prefix)sc_ktkZ_%(sqrts)s, %(prefix)sc_kbkZ_%(sqrts)s)'%locals()
                self.modelBuilder.factory_(coeffSum)

                rooExpr = 'expr::%(scalingName)s( "( (@0*@0)*(@3)  + (@1*@1)*(@4) + (@2*@2)*(@5) + (@0*@1)*(@6)  + (@0*@2)*(@7) + (@1*@2)*(@8) )/@9",\
		%(Ctop)s, %(Cb)s, %(CZ)s, %(prefix)sc_kt2_%(sqrts)s, %(prefix)sc_kb2_%(sqrts)s, %(prefix)sc_kZ2_%(sqrts)s, %(prefix)sc_ktkb_%(sqrts)s, %(prefix)sc_ktkZ_%(sqrts)s, %(prefix)sc_kbkZ_%(sqrts)s, coeff_sum_%(scalingName)s)'%locals()
                self.modelBuilder.factory_(rooExpr)
        elif what.startswith('tHq'):
	    coeffs = {'7TeV':[3.099,3.980,-6.078], '8TeV':[2.984,3.886,-5.870],'13TeV':[2.633,3.578,-5.211],'14TeV':[2.582,3.538,-5.120]}  # coefficients for  kt^{2}, kW^{2}, ktkW @ MH = 125 GeV
            for sqrts in ('7TeV', '8TeV','13TeV','14TeV'):
                scalingName = 'Scaling_'+what+'_'+sqrts
                rooExpr = 'expr::%(scalingName)s'%locals()+'( "( (@0*@0)*(%g)  + (@1*@1)*(%g) + (@0*@1)*(%g) )/%g"'%tuple((coeffs[sqrts] + [sum(coeffs[sqrts])] ))+', %(Ctop)s, %(CW)s)'%locals()
                self.modelBuilder.factory_(rooExpr)
        elif what.startswith('tHW'):
	    coeffs = {'7TeV':[2.306,1.697,-3.003], '8TeV':[2.426,1.818,-3.244],'13TeV':[2.909,2.310,-4.220],'14TeV':[2.988,2.397,-4.385]}  # coefficients for  kt^{2}, kW^{2}, ktkW @ MH = 125 GeV
            for sqrts in ('7TeV', '8TeV','13TeV','14TeV'):
                scalingName = 'Scaling_'+what+'_'+sqrts
                rooExpr = 'expr::%(scalingName)s'%locals()+'( "( (@0*@0)*(%g)  + (@1*@1)*(%g) + (@0*@1)*(%g) )/%g"'%tuple((coeffs[sqrts] + [sum(coeffs[sqrts])] ))+', %(Ctop)s, %(CW)s,)'%locals()
                self.modelBuilder.factory_(rooExpr)
        else:
            raise RuntimeError, "There is no scaling defined for %(what)s" % locals()

                
    def makePartialWidthUncertainties(self):
        THU_GROUPS = [
           ('hvv' , [ 'hww', 'hzz' ] ),
           ('hqq' , [ 'hbb', 'hcc', 'hss' ] ),
           ('hll' , [ 'htt', 'hmm' ] ),
           ('hgg' , [ 'hgg' ] ),
           ('hzg' , [ 'hzg' ] ),
           ('hgluglu' , [ 'hgluglu' ] ),
        ]
        widthUncertainties = {}; widthUncertaintiesKeys = []
        for line in open(self.brpath+"/WidthUncertainties_125GeV_YR4.txt"):
            if widthUncertaintiesKeys == []:
                widthUncertaintiesKeys = line.split()[1:]
            else:
                fields = line.split()
                widthUncertainties[fields[0]] = dict([(k,0.01*float(v)) for (k,v) in zip(widthUncertaintiesKeys, fields[1:])]) 
        for K in widthUncertaintiesKeys[:-1]:
            self.modelBuilder.doVar("param_%s[-7,7]" % K)
        for K, DS in THU_GROUPS:
            self.modelBuilder.doVar("HiggsDecayWidthTHU_%s[-7,7]" % K)
        for D in ALL_HIGGS_DECAYS:
            #print "For decay %s: " % D,
            if D not in widthUncertainties:
                self.modelBuilder.doVar("HiggsDecayWidth_UncertaintyScaling_%s[1]" % D)
                #print " no uncertainties."
                continue
            pnorm = ROOT.ProcessNormalization("HiggsDecayWidth_UncertaintyScaling_%s" %D, "")
            for K in widthUncertaintiesKeys:
                if K == "thu":
                    var = None
                    for K2, DS in THU_GROUPS:
                        if D in DS:
                            var = self.modelBuilder.out.var("HiggsDecayWidthTHU_%s" % K2)
                            break
                    if var == None: continue
                else:
                    var = self.modelBuilder.out.var("param_%s" % K)
                #print " [ %+.1f%% from %s ]  " % (widthUncertainties[D][K]*100, var.GetName()),
                pnorm.addLogNormal(exp(widthUncertainties[D][K]), var)
            #print "."
            self.modelBuilder.out._import(pnorm)
    def dump(self,name,xvar,values,logfile):
        xv = self.modelBuilder.out.var(xvar)
        yf = self.modelBuilder.out.function(name)
        if yf == None: raise RuntimeError, "Missing "+name
        log = open(logfile, "w")
        for x in values:
            xv.setVal(x)
            log.write("%.3f\t%.7g\n" % (x, yf.getVal()) )
    def textToSpline(self,name,filename,xvar="MH",ycol=1,xcol=0,skipRows=1,algo="CSPLINE"):
        if (self.modelBuilder.out.function(name) != None): return
        x = []; y = []
        file = open(filename,'r')
        lines = [l for l in file]
        for line in lines[skipRows:]:
            if len(line.strip()) == 0: continue
            cols = line.split();
            x.append(float(cols[xcol]))
            y.append(float(cols[ycol]))
        xv = self.modelBuilder.out.var(xvar)
        spline = ROOT.RooSpline1D(name, "file %s, x=%d, y=%d" % (filename,xcol,ycol), xv, len(x), array('d', x), array('d', y), algo)
        self.modelBuilder.out._import(spline)

#if __name__ == "__main__":
#   sm = SMHiggsBuilder()
