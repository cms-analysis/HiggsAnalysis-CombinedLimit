#from HiggsAnalysis.CombinedLimit.LHCGModels import *
from HiggsAnalysis.CombinedLimit.PhysicsModel import *


### Conventions:
## https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsTemplateCrossSection

# For Stage0, the following convention for the production modes should be used:
# 
#
Stage0 = ['ggH', 'qqH', 'VH_had', 'WH_lep', 'ZH_lep', 'ggZH_lep', 'ttH', 'bbH', 'tH'] 

# For Stage1, the convention is to add the suffices in SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h. So for ggH the labels would be:
# Stage0->Stage1
Stage1 = {
        'ggH':
        ['ggH_VBFTOPO_JET3VETO', 'ggH_VBFTOPO_JET3',
            'ggH_0J',
            'ggH_1J_PTH_0_60',   'ggH_1J_PTH_60_120',
            'ggH_1J_PTH_120_200', 'ggH_1J_PTH_GT200',
            'ggH_GE2J_PTH_0_60',    'ggH_GE2J_PTH_60_120', 
            'ggH_GE2J_PTH_120_200', 'ggH_GE2J_PTH_GT200',
            ],
        'qqH':
        [  'qqH_fwd',
            'qqH_VBFTOPO_JET3VETO' , 'qqH_VBFTOPO_JET3' ,
            'qqH_VH2JET' , 'qqH_REST' , 'qqH_PTJET1_GT200' ,
            ],
        #// qq -> WH
        'WH_lep':
        [
            #WHlep_fwd ,
            'WH_lep_PTV_0_150' ,
            'WH_lep_PTV_150_250_0J' ,
            'WH_lep_PTV_150_250_GE1J' ,
            'WH_lep_PTV_GT250' ,
            ],
        #// qq -> ZH
        'ZH_lep':
        [
            #'ZH_lep_fwd' ,
            'ZH_lep_PTV_0_150' ,
            'ZH_lep_PTV_150_250_0J' ,
            'ZH_lep_PTV_150_250_GE1J' ,
            'ZH_lep_PTV_GT250' ,
            ],
        #// gg -> ZH
        'ggZH_lep':
        [
            #ggZH_lep_fwd ,
            'ggZH_lep_PTV_0_150' ,
            'ggZH_lep_PTV_GT150_0J' ,
            'ggZH_lep_PTV_GT150_GE1J' ,
            ],
        'VH_had':
        [
            'WH_had_VBFTOPO_JET3VETO',
            'WH_had_VBFTOPO_JET3',
            'WH_had_VH2JET',
            'WH_had_REST',
            'WH_had_PTJET1_GT200',
            'ZH_had_VBFTOPO_JET3VETO',
            'ZH_had_vbfTOPO_JET3',
            'ZH_had_VH2JET',
            'ZH_had_REST',
            'ZH_had_PtJET1_GT200',
        ],
        #// ttH
        'ttH':
        [
            #'TTH_fwd' , 
            'ttH' ,
            ],
        #// bbH -> these are not in the eft model
        #'bbH':['bbH'],
        #BBH_fwd , BBH ,
        #// tH
        #'tH':['tH',],
        #TH_fwd , TH 
}

ALL_STAGE1 = [ x for v in Stage1.itervalues() for x in v]

#https://twiki.cern.ch/twiki/bin/view/LHCPhysics/STXStoEFT

#class STXStoEFT(LHCHCGBaseModel):
class STXStoEFT_LO(SMLikeHiggsModel):
    def __init__(self):
        self.floatMass = False
        self.cWWcB=False ## cWW-cB -> sum and diff pois
        return 

    def readTxt(self):
        ''' Read the txt data and save them as dictionary'''
        import re
        import json # to write ?
        import os

        txt_bins={}
        self.txt_decays={}
        self.txt_crosssections={}
        self.txt_pois=[]
        
        data=os.environ["CMSSW_BASE"]+"/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg/eft"
        f1=open('/'.join([data,'bin.txt']))
        for line in f1:
            line=line.split('#')[0]
            if len(line.split()) <2: continue
            name=line.split()[0]
            binn=line.split()[1]

            parsed= re.sub('r_','',name)
            parsed= parsed.upper()
            parsed=re.sub('GG2H','ggH',parsed)
            parsed=re.sub('VBF_QQ2HQQ','qqH',parsed)
            parsed=re.sub('QQ2HLNU','WH_lep',parsed)
            parsed=re.sub('QQ2HLL','ZH_lep',parsed)
            parsed=re.sub('WH_QQ2HQQ','WH_had',parsed)
            parsed=re.sub('ZH_QQ2HQQ','ZH_had',parsed)
            parsed=re.sub('TTH','ttH',parsed)
            
            print "Bin:",binn,"orig",name,"->",parsed
            txt_bins[int(binn)]=parsed
        f1.close()
        f2=open('/'.join([data,"crosssections.txt"]))
        #self.txt_crosssections={}
        last_bin=None
        for line in f2:
            line=line.split('#')[0]
            if len(line.split()) ==0 : continue

            if line.startswith("Bin number"):
                last_bin=int( re.sub('Bin number: ','',line).split()[0])

            if line.startswith('1'): ## perturbation theory. I could put an else but in this way should be safer
                if last_bin==None: raise ValueError("Unable to parse. I should not get a cross section before a Bin specification")
                formula=line.replace('\n','')

                #print "DEBUG: looking for bin:",last_bin,"in",txt_bins
                #print "DEBUG: --> ",txt_bins[last_bin]
                self.txt_crosssections[txt_bins[last_bin]]= formula
                for p in formula.split():
                    if p=='+': continue
                    if p=='-': continue
                    if p=='*': continue
                    try:
                        num=float(p)
                    except ValueError:
                        self.txt_pois.append(p)
        self.txt_pois = list(set(self.txt_pois)) #uniq
        f2.close()
        
        f3=open('/'.join([data,"decay.txt"]))
        last_bin=None
        for line in f3:
            line=line.split('#')[0]
            if len(line.split()) ==0 : continue

            if line.startswith("Bin number"):
                name=line.split()[-1] ## this is a string not a number
                name=re.sub("hzzto4l",'hzz',name)
                last_bin=name

            if line.startswith('1'): ## perturbation theory. I could put an else but in this way should be safer
                if last_bin==None: raise ValueError("Unable to parse. I should not get a cross section before a Bin specification")
                formula=line.replace('\n','')
                self.txt_decays[last_bin]= formula
                for p in formula.split():
                    if p=='+': continue
                    if p=='-': continue
                    if p=='*': continue
                    try:
                        num=float(p)
                    except ValueError:
                        self.txt_pois.append(p)
        self.txt_pois = list(set(self.txt_pois)) #uniq
        
        print "--- POIS ---"
        print json.dumps(self.txt_pois)
        print "--- XS ---"
        print json.dumps(self.txt_crosssections)
        print "--- DECAYS ---"
        print json.dumps(self.txt_decays)
                
    def setPhysicsOptionsBase(self,physOptions):
        for po in physOptions:
            if po.startswith("higgsMassRange="):
                self.floatMass = True
                self.mHRange = po.replace("higgsMassRange=","").split(",")
                print 'The Higgs mass range:', self.mHRange
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema"
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Extrama for Higgs mass range defined with inverterd order. Second must be larger the first"
            if po.startswith("cWW-cB="):
                self.cWWcB=True
        return ## PO

    def doMH(self):
        if self.floatMass:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
            else:
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1])) 
        else:
            if self.modelBuilder.out.var("MH"):
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
                self.modelBuilder.out.var("MH").setConstant(True)
            else:
                self.modelBuilder.doVar("MH[%g]" % self.options.mass)
        return ## doMH

    def getHiggsSignalYieldScale(self, production, decay, energy):
        #if production not in SM_HIGGS_PROD: raise ValueError("production %s is not in the list of allowed production modes: %s"%(production,str(SM_HIGGS_PROD)))
        if decay not in SM_HIGG_DECAYS: raise ValueError("decay %s is not in the list of allowed decay modes: %s"%(decay, str(SM_HIGGS_DECAYS)))

        if decay not in self.txt_decays: raise ValueError("decay scalinig %s not available from the txt files."%decay)

        if production in Stage0:
            raise ValueError("Not implemented")
        elif production in ALL_STAGE1:
            return "scaling_%s_%s"%(production,decay)
        else:
            raise ValueError("production process '%s' not in Stage0/Stage1: %s"%(process, str(ALL_STAGE1)))

        #return "r"
        return 

    def doParametersOfInterest(self):
        if self.floatMass: 
            print "[WARNING] Not translate in EFT variations."
            self.doMH()
        self.readTxt()

        for poi in self.txt_pois:
            par_range="[0,-5,5]"
            if poi=='cG': par_range="0,-5e-4,5e-4"
            elif poi=='tcG': par_range="0,-5e-4,5e-4"
            elif poi=='cA': par_range="0,-2e-3,2e-3"
            elif poi=='tcA': par_range="0,-2e-3,2e-3"
            elif poi=='cH': par_range="0,-1,1"
            elif poi=='cHW': par_range="0,-0.05,0.05"
            elif poi=='tcHW': par_range="0,-0.1,0.1"
            elif poi=='cHB': par_range="0,-0.1,0.1"
            elif poi=='tcHB': par_range="0,-0.5,0.5"

            ## cWW - cB   -> -0.035, 0.005
            ##     +         -0.0033, 0.0018
            if self.cWWcB and (poi =='cWW' or poi=='cB'): continue
            self.modelBuilder.doVar(poi+par_range)


        POIs=','.join(self.txt_pois)
        if self.cWWcB: 
            self.modelBuilder.doVar('cWWPluscB[0,-0.005,0.005]')
            self.modelBuilder.doVar('cWWMinuscB[0,-0.05,0.05]')

            self.modelBuilder.factory_("expr::cWW(\"0.5*(@0+@1)\",{CWWPluscB,cWWMinusCB})")
            self.modelBuilder.factory_("expr::cB(\"0.5*(@0-@1)\",{CWWPluscB,cWWMinusCB})")
            POIs=','.join([ x for x in self.txt_pois if not x=='cWW' and not x=='cB' ]) + ',cWWPluscB,cWWMinuscB'
        self.modelBuilder.doSet("POI",POIs)
    

        dependence='{'+','.join(self.txt_pois)+'}'
        for key in self.txt_crosssections:
            self.modelBuilder.doExp(key, self.txt_crosssections[key],dependence)
        for key in self.txt_decays:
            self.modelBuilder.doExp(key, self.txt_decays[key],dependence)

        ## do scaling
        for D in self.decays:
            for P in self.txt_crossections:
                self.modelBulider.factory_("prod::scaling_%s_%s(%s)"%(P,D,",".join([P,D])))

stxsToEFT=STXStoEFT_LO()
if __name__ == "__main__":
    print " ------------------- "
    stxsToEFT.readTxt()
    print "--------------------"

