### Original Author: Andrea Carlo Marini
### Date: 9/11/2017

import ROOT 
from array import array
import re

# construct model and datacard
# type = cut & count / pdf
# bkg = yes / no

class Model():
    def __init__(self):
        self.type="c&c" ## pdf
        self.w = None
        ## truth bins
        self.nbins_x = 5
        ## reco bins
        self.nbins_y = 5 ## if different, TODO: likely not work the diagonal--off-diagonal matching, check for i==j
        ## ntruth events
        self.nevents=1000
        self.bkg=False ## non resonant bkg
        self.nbkg=100
        ## smearings, and efficiencies
        self.smearings = [0.75,.10,.025] ## diagonal, off-diagonola, 2nd ...
        self.eff = [.85,.85,.90,.90,.92,.92,.93,.94,.95,.99]  ## one per bin

        ###output file name
        self.fname = "workspace.root"
        self.dname = "datacard.txt"

        #Internal

    def __del__(self):
        return

    def info(self):
        print "-------------------------"
        print "bins truth:",self.nbins_x,"smeared",self.nbins_y
        print "n. events: sig:",self.nevents, ("bkg: %s"%self.nbkg) if self.bkg else "bkg: None"
        print "type:",self.type
        print "workspace:",self.fname
        print "datacard:",self.dname
        print "-------------------------"

    def run(self):
        self.doModel()
        self.doWorkspace()
        self.doDatacard()
        self.info()
        return self

    def doModel(self):
        xmin=0.
        xmax=10.
        self.x_cont = ROOT.RooRealVar("x_cont","x",xmin,xmax)
        self.pdf_cont = ROOT.RooGenericPdf("pt-sig","@0*@0*@0*TMath::Exp(-@0)",ROOT.RooArgList(self.x_cont))
        # this actually lives on the reco specturm
        self.bkg_cont = ROOT.RooGenericPdf("pt-bkg","Exp(-@0)",ROOT.RooArgList(self.x_cont))
        #self.cdf_pt = self.pt_cont.createCdf()

        ### 
        ## binned model
        self.x_th1 = ROOT.TH1D("x_th1","binned version of truth model",self.nbins_x, xmin,xmax)
        for i in range(0,self.nbins_x):
            rangename = "truth_range_%d"%i
            low  = i 
            high = i+1
            self.x_cont . setRange(rangename, low,high)
            count = self.pdf_cont.createIntegral(ROOT.RooArgSet(self.x_cont),ROOT.RooArgSet(self.x_cont),rangename).getVal() * self.nevents
            self.x_th1 .SetBinContent( self.x_th1.FindBin(i+0.5), count)
        ## construct smearing matrix
        ## y = M.x + b
        self.Mhat=ROOT.TH2D("Mhat","response probability matrix (X=truth, Y=Reco)",self.nbins_x,xmin,xmax,self.nbins_y,xmin,xmax)
        for i in range(0,self.nbins_x):
            for j in range(0,self.nbins_y):
                binX= self.Mhat.GetXaxis().FindBin(i+0.5)
                binY= self.Mhat.GetYaxis().FindBin(j+0.5)
                if len(self.smearings) >abs(i-j): self.Mhat.SetBinContent(binX,binY, self.smearings[abs(i-j)])
        # normalize Mhat to 1 per gen lines times efficiencies (are in truth dimension efficiencies)
        for i in range(0,self.nbins_x):
            binX= self.Mhat.GetXaxis().FindBin(i+0.5)
            S=0.
            for j in range(0,self.nbins_y):
                binY= self.Mhat.GetYaxis().FindBin(j+0.5)
                S+= self.Mhat.GetBinContent(binX,binY)
            for j in range(0,self.nbins_y):
                binY= self.Mhat.GetYaxis().FindBin(j+0.5)
                c= self.Mhat.GetBinContent(binX,binY)
                self.Mhat.SetBinContent(binX,binY, c * self.eff[i]/S )
        ## construct background
        self.b_th1 = ROOT.TH1D("b_th1","background count in each bin",self.nbins_y,xmin,xmax)
        if self.bkg:
            for i in range(0,self.nbins_y):
                rangename = "bkg_range_%d"%i
                low  = i 
                high = i+1
                self.x_cont . setRange(rangename, low,high)
                count = self.pdf_cont.createIntegral(ROOT.RooArgSet(self.x_cont),ROOT.RooArgSet(self.x_cont),rangename).getVal() * self.nbkg
                self.b_th1 .SetBinContent( self.b_th1.FindBin(i+0.5), count)
            
        ## construct reco y = Mhat . x + b
        self.y_th1 = ROOT.TH1D("y_th1","binned version of reco model",self.nbins_y, xmin,xmax)
        for j in range(0,self.nbins_y):
            binY= self.Mhat.GetYaxis().FindBin(j+0.5)
            S=0.
            for i in range(0,self.nbins_x):
                binX= self.Mhat.GetXaxis().FindBin(i+0.5)
                m = self.Mhat.GetBinContent(binX,binY)
                x = self.x_th1.GetBinContent(binX)
                S+= m*x
            b = self.b_th1.GetBinContent(binY)
            y= S+b
            self.y_th1.SetBinContent(binY,y)

        ## construct M for future references / RooUnfold way
        self.M=ROOT.TH2D("M","response matrix (X=truth, Y=Reco)",self.nbins_x,xmin,xmax,self.nbins_y,xmin,xmax)
        for i in range(0,self.nbins_x):
            for j in range(0,self.nbins_y):
                binX= self.Mhat.GetXaxis().FindBin(i+0.5)
                binY= self.Mhat.GetYaxis().FindBin(j+0.5)
                c = self.Mhat.GetBinContent(binX,binY)
                x = self.x_th1.GetBinContent(binX)
                self.M.SetBinContent(binX,binY,c * x) 

    def Import(self,obj):
        print "* importing",obj.GetName()
        #getattr(self.w,'import')(obj,ROOT.RooCmdArg())
        getattr(self.w,'import')(obj,ROOT.RooCmdArg(ROOT.RooFit.RecycleConflictNodes()))

    def doWorkspace(self):
        self.w = ROOT.RooWorkspace("w","w")
        ROOT.SetOwnership(self.w,False)
        ##self.pdf = ROOT.TF1("myfunc","x*x*x*TMath::Exp(-x)/5.9380",0,10)
        #create observable x (truth)
        print "-> construct real or fake pdfs for each entry in the matrix" 
        self.mgg = ROOT.RooRealVar("mgg","mgg",110,150)
        self.MH = ROOT.RooRealVar("MH","MH",125.)

        self.Import(self.mgg)
        self.Import(self.MH )

        ## construct data
        #self.y_th1 = ROOT.TH1D("y_th1","binned version of reco model",self.nbins_y, xmin,xmax)
        #self.dataObs=[]
        for j in range(0,self.nbins_y):
            nbins=200 ## every 0.2 GeV
            if self.type=='c&c':nbins=1
            h=ROOT.TH1D("h_RecoBin%d"%j,"Reco empty histogram",nbins,110,150) # in mgg
            dh = ROOT.RooDataHist("dataObs_RecoBin%d"%j,"Data Hist",ROOT.RooArgList(self.mgg),h.Clone())
            self.Import(dh)
            #self.dataObs.extend([h,dh])
            #ROOT.SetOwnership(h,False)
            #ROOT.SetOwnership(dh,False)
        
        ## construct pdfs
        if self.type =='pdf':
            self.sigmaD = ROOT.RooRealVar("sigmaD","sigma diagonal elements",1.)
            self.sigmaO = ROOT.RooRealVar("sigmaO","sigma off-diagonal elements",2.)

            self.Import(self.sigmaD)
            self.Import(self.sigmaO)

        for i in range(0,self.nbins_x):
            for j in range(0,self.nbins_y):
                name = "GenBin%d_RecoBin%d"%(i,j)

                binX= self.Mhat.GetXaxis().FindBin(i+0.5)
                binY= self.Mhat.GetYaxis().FindBin(j+0.5)

                norm = ROOT.RooRealVar(name +"_norm","normalization of xxx",self.M.GetBinContent(binX,binY))
                norm.setConstant()
                if self.type == 'c&c':
                    pdf = ROOT.RooUniform(name,"Fake pdf for "+name,ROOT.RooArgSet(self.mgg))
                elif self.type == 'pdf':
                    if i==j:pdf = ROOT.RooGaussian(name,"Gauss pdf for "+name,self.mgg,self.MH,self.sigmaD)
                    else:pdf = ROOT.RooGaussian(name,"Gauss pdf for "+name,self.mgg,self.MH,self.sigmaO)
                self.Import(norm)
                self.Import(pdf)

        if self.bkg:
            if self.type =='pdf':
                self.tau = ROOT.RooRealVar("tau","tau for bkg",-0.2)
                self.Import(self.sigmaD)

            for j in range(0,self.nbins_y):
                name = "Bkg_RecoBin%d"%(j)
                binY= self.Mhat.GetYaxis().FindBin(j+0.5)
                norm = ROOT.RooRealVar(name +"_norm","normalization of xxx",self.b_th1.GetBinContent(binY))
                norm.setConstant()
                
                if self.type=='c&c':
                    pdf = ROOT.RooUniform(name,"Fake pdf for "+name,ROOT.RooArgSet(self.mgg))
                elif self.type == 'pdf':
                    pdf = ROOT.RooExponential(name,"Exp pdf for "+name,self.mgg,self.tau)

                self.Import(norm)
                self.Import(pdf)
                    

        print "-> writing worspace to file",self.fname
        self.w.writeToFile(self.fname)
        print "-> adding  histograms to",self.fname
        self.fOut = ROOT.TFile.Open(self.fname,"UPDATE")
        self.fOut.cd()
        self.y_th1.Write()
        self.x_th1.Write()
        self.Mhat.Write()
        self.M.Write()
        self.fOut.Close()

    def doDatacard(self):
        self.datacard=open(self.dname,"w")
        self.datacard.write("------------------------\n")
        self.datacard.write("## Automatic created by makeModel\n")
        self.datacard.write("## Original Author: Andrea Carlo Marini\n")
        self.datacard.write("## Original Date: 9 Nov 2017\n")
        self.datacard.write("------------------------\n")
        self.datacard.write("* ibin\n")
        self.datacard.write("* jbin\n")
        self.datacard.write("* kbin\n")
        self.datacard.write("------------------------\n")
        ## shapes declaration
        #for i in range(0,self.nbins_x):
        #    for j in range(0,self.nbins_y):
        #name = "GenBin%d_RecoBin%d"%(i,j)
        #name = "Bkg_RecoBin%d"%(j)
        self.datacard.write("shapes data_obs * "+self.fname+" w:dataObs_$CHANNEL\n")
        self.datacard.write("shapes * * "+self.fname+" w:$PROCESS_$CHANNEL\n")
        self.datacard.write("------------------------\n")
        self.datacard.write("bin\t" + '\t'.join( ['RecoBin%d'%j for j in range(0,self.nbins_y) ] ) + '\n')
        self.datacard.write("observation\t" + '\t'.join( ['-1' for j in range(0,self.nbins_y) ] ) + '\n')
        ## observation
        self.datacard.write("------------------------\n")
        binline = "bin"
        procline = "process"
        procline2 = "process"
        rateline = "rate"

        for i in range(0,self.nbins_x):
            for j in range(0,self.nbins_y):
                binline +="\t"+'RecoBin%d'%j 
                procline += "\t"+'GenBin%d'%i
                procline2 += "\t"+'%d'%(-i)
                rateline += "\t" +"1"

        if self.bkg:
            for j in range(0,self.nbins_y):
                binline +="\t"+'RecoBin%d'%j 
                procline += "\t"+'Bkg'
                procline2 += "\t"+'1' ## >0: bkg ; <=0: sig
                rateline += "\t" +"1"

        self.datacard.write(binline +"\n"+ procline + "\n" + procline2+"\n" + rateline +"\n")
        ## expected
        self.datacard.write("------------------------\n")
        ## systematics

        ## write regularization lines
        self.datacard.write("#######################################\n")
        self.datacard.write("###### Regularization lines. Uncomment to enable regularization.\n")
        self.datacard.write("###### Requisites for regularization studies (PR 531):\n")
        self.datacard.write("######  git pull git@github.com:amarini/HiggsAnalysis-CombinedLimit/topic_regularization_102x\n")
        for i in range(0,self.nbins_x):
            if i==0 or i==self.nbins_x-1: continue
            ## different strength vars
            #self.datacard.write("## constr%d"%i+ " constr const%(bin)d_In[0.],RooFormulaVar::fconstr%(bin)d(\"r_Bin%(prev)d+r_Bin%(next)d-2*r_Bin%(bin)d\",{r_Bin%(prev)d,r_Bin%(bin)d,r_Bin%(next)d}),constr%(bin)d_S[%(delta)s]\n" %{"bin":i,"next":i+1,"prev":i-1,"delta":0.03} )
            ## one -> delta
            #self.datacard.write("## constr%d"%i+ " constr const%(bin)d_In[0.],RooFormulaVar::fconstr%(bin)d(\"r_Bin%(prev)d+r_Bin%(next)d-2*r_Bin%(bin)d\",{r_Bin%(prev)d,r_Bin%(bin)d,r_Bin%(next)d}),delta[%(delta)s]\n" %{"bin":i,"next":i+1,"prev":i-1,"delta":0.03} )
            #old -> new
            self.datacard.write("## constr%d"%i+ " constr r_Bin%(prev)d+r_Bin%(next)d-2*r_Bin%(bin)d delta[%(delta)s]\n" %{"bin":i,"next":i+1,"prev":i-1,"delta":0.03} )

        self.datacard.write("#######################################\n")
        self.datacard.write("###### TUnfold like\n")
        for i in range(0,self.nbins_x):
            if i==0 or i==self.nbins_x-1: continue
            d = {"bin":i,"next":i+1,"prev":i-1,"delta":0.03}
            lambdaMCList =[]
            lambdaMCPList =[]
            lambdaMCNList =[]
            varList=[]
            for j in range(0,self.nbins_y):
                #norm ='GenBin%d_RecoBin%d_norm'%(i,j )
                #p ='GenBin%d_RecoBin%d_norm'%(i-1,j )
                #n ='GenBin%d_RecoBin%d_norm'%(i+1,j )
                norm ='shapeSig_GenBin%d_RecoBin%d__norm'%(i,j )
                p ='shapeSig_GenBin%d_RecoBin%d__norm'%(i-1,j )
                n ='shapeSig_GenBin%d_RecoBin%d__norm'%(i+1,j )
                lambdaMCList .append(norm)
                lambdaMCPList .append(p)
                lambdaMCNList .append(n)
                varList.extend([norm,p,n])

            d['lambdaMC'] = "("+'+'.join(lambdaMCList)+")"
            d['lambdaMC_prev'] = "("+'+'.join(lambdaMCPList)+")"
            d['lambdaMC_next'] = "("+'+'.join(lambdaMCNList)+")"
            d['lambdaMC_vars'] = ','.join(varList)

            ##self.datacard.write("## constr%d"%i+ " constr const%(bin)d_In[0.],RooFormulaVar::fconstr%(bin)d(\"(r_Bin%(prev)d-1.)*%(lambdaMC_prev)s+r_Bin%(next)d*%(lambdaMC_next)s-2*r_Bin%(bin)d*%(lambdaMC)s\",{r_Bin%(prev)d,r_Bin%(bin)d,r_Bin%(next)d,%(lambdaMC_vars)s}),constr%(bin)d_S[%(delta)s]\n" % d )
            #self.datacard.write("## constr%d"%i+ " constr const%(bin)d_In[0.],RooFormulaVar::fconstr%(bin)d(\"(r_Bin%(prev)d-1.)*%(lambdaMC_prev)s+r_Bin%(next)d*%(lambdaMC_next)s-2*r_Bin%(bin)d*%(lambdaMC)s\",{r_Bin%(prev)d,r_Bin%(bin)d,r_Bin%(next)d,%(lambdaMC_vars)s}),delta[%(delta)s]\n" % d )
            self.datacard.write("## constr%d"%i+ " constr (r_Bin%(prev)d-1.)*%(lambdaMC_prev)s+(r_Bin%(next)d-1.)*%(lambdaMC_next)s-2*(r_Bin%(bin)d-1.)*%(lambdaMC)s {r_Bin%(prev)d,r_Bin%(bin)d,r_Bin%(next)d,%(lambdaMC_vars)s} delta[%(delta)s]\n" % d )
            #self.datacard.write("## constr%d"%i+ " constr (r_Bin%(prev)d-1.)*%(lambdaMC_prev)s+r_Bin%(next)d*%(lambdaMC_next)s-2*r_Bin%(bin)d*%(lambdaMC)s delta[%(delta)s]\n" % d )
        self.datacard.write("#######################################\n")

        self.datacard.write("\n")
        self.datacard.write("\n")
        self.datacard.write("\n")
        ## write commands

        #for i in range(nbinsG):
        #   cmd+= " --PO map='.*Bin%d.*':r_Bin%d[1,-1,20]"%(i,i)
        #   cmd+= " --PO map='.*Bkg.*':r_Bkg[1]"
        rootname = re.sub('.txt','.root',self.dname)
        self.datacard.write("###############################\n")
        self.datacard.write("##                             \n")
        self.datacard.write("## Commands:\n")
        self.datacard.write("## text2workspace.py -m 125 --X-allow-no-background -o "+rootname+" "+self.dname+"\n")
        pois_map = ' '.join(["--PO map='.*GenBin%d.*:r_Bin%d[1,-1,20]'"%(i,i) for i in range(0,self.nbins_x)])
        self.datacard.write("##   -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel "+pois_map+"\n")
        self.datacard.write("##\n")
        binbybin = ' '.join(["--PO map='.*RecoBin%d.*:r_Bin%d[1,-1,20]'"%(i,i) for i in range(0,self.nbins_y)])
        self.datacard.write("## for bin by bin use:\n")
        self.datacard.write("##   -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel "+binbybin+"\n")
        self.datacard.write("##\n")

        ## in the old combine version this was setPhysicsModelParameters
        expect="--setParameters="+','.join(["r_Bin%d=1"%i for i in range(0,self.nbins_x)])
        self.datacard.write("## combine -M MultiDimFit "+expect+" -t -1 -m 125 "+rootname+"\n")
        self.datacard.write("## combine -M MultiDimFit "+expect+" -t -1 -m 125 --algo=grid --points=100 -P r_Bin1 --setParameterRanges r_Bin1=0.5,1.5 --floatOtherPOIs=1 "+rootname+"\n")
        self.datacard.write("###############################\n")
        self.datacard.close()

########################
##        MAIN        ##
########################
if __name__=="__main__":
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("","--doCC",action='store_true',help="do Cut and Count [%default]",default=False)
    opts,args = parser.parse_args()

    model=Model()
    model.type='pdf'
    if opts.doCC : 
        model.type = 'c&c'
        model.fname = "workspace_cc.root"
        model.dname = "datacard_cc.txt"
    #model.type='c&c'
    model.run()
    print "---- DONE ----"
