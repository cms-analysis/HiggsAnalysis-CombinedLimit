import sys,os
from subprocess import call

### Original Author: Andrea Carlo Marini
### 21/08/2016

#global options
from optparse import OptionParser
parser=OptionParser()
#parser.add_option("-n","--name",help="NAME [%default]",default="")
parser.add_option("-o","--output",dest="output",help="output name [%default]",default="regularization.root")
parser.add_option("-w","--ws",dest="workspace",help="ws name [%default]",default="w")
parser.add_option("-d","--datacard",dest="datacard",help="datacard name [%default]",default="cms_datacard.txt")
parser.add_option("-r","--regularize",dest="reg",action='store_true',help="activate regularization [%default]",default=False)
opts,args= parser.parse_args()

# avoid that root reads argv
sys.argv=[]
import ROOT
# set batch mode
ROOT.gROOT.SetBatch()
ROOT.RooMsgService.instance().setSilentMode(True);

## create a workspace and a datacard
datacard=open(opts.datacard,"w")


ws=ROOT.RooWorkspace(opts.workspace,opts.workspace)

datacard.write("## Datacard automatically created by regularization/createWs.py\n")
datacard.write("## Original Author: Andrea Carlo Marini\n")
datacard.write("-------------------------------------\n")
datacard.write("imax *\n")
datacard.write("jmax *\n")
datacard.write("kmax *\n")
datacard.write("-------------------------------------\n")

# create smearing matrix -- squared for simplicity (mhat construction is not good for not squared ones)
nbinsR=4
nbinsG=4
mhat={}
for row in range(nbinsG):
   for col in range(nbinsR):
      #mhat[ (row,col) ] = ROOT.TMath.Gaus( row-col, 1)
      if col==row:
      	mhat[ (row,col) ] = .99
      else:
      	mhat[ (row,col) ] = 0.01 / abs(row-col)

#efficiency 95% +/- 0.02 random
for row in range(nbinsG):
   S=0.0
   for col in range(nbinsR):
      S+=mhat[ (row, col) ]
   for col in range(nbinsR):
      mhat[ (row,col) ] *= ROOT.TMath.Min(ROOT.TMath.Gaus(0.95,0.02),1.)/S

# create/inject gen level spectra
gen = ROOT.TH1D("gen","gen",nbinsG,0,nbinsG)
gen.SetBinContent(0,0) # underflow
gen.SetBinContent(nbinsG+1,0) #overflow
for i in range(nbinsG):
   gen.SetBinContent(i+1, 1000*ROOT.TMath.Exp(- i/0.5) ) 

bkg= ROOT.TH1D("bkg","bkg",nbinsR,0,nbinsR)
bkg.SetBinContent(0,0) # underflow
bkg.SetBinContent(nbinsR+1,0) #overflow
for i in range(nbinsR):
   bkg.SetBinContent(i+1, 2*ROOT.TMath.Exp(-i/1.0) )

reco = ROOT.TH1D("reco","reco",nbinsR,0,nbinsR)
reco.SetBinContent(0,0) # underflow
reco.SetBinContent(nbinsR+1,0) #overflow
for col in range(nbinsR):
   S=0.0
   for row in range(nbinsG):
      S+=mhat[ (row,col) ] * gen.GetBinContent(row+1)
   S+=bkg.GetBinContent(col+1) ## add bkg
   ## if you want smear it here
   ## S=ROOT.TMath.Poisson(S)
   ##
   ## BIAS
   if col == 2: 
   	   print "-> BIAS in Bin=2, to study reg"
   	   S *= 2
   reco.SetBinContent(col+1,S)

ws.factory("x[0,%d]"%(nbinsR+1))
x=ws.var("x")
data_obs= ROOT.RooDataHist("data_obs","data_obs",ROOT.RooArgList(x),reco)
getattr(ws,'import')(data_obs,ROOT.RooCmdArg()) 

datacard.write("shapes data_obs *\t" + opts.output + "\t" + opts.workspace+":data_obs\n")
## construct model

########################
def importPdf(h,name):
   ''' Import a Pdf based on hist with name. ws and x global vars used'''
   dh=ROOT.RooDataHist( "roohist_" + h.GetName(),h.GetName(),ROOT.RooArgList(x),h)
   pdf=ROOT.RooHistPdf(name,"dh_Bin%d",ROOT.RooArgSet(x),dh)
   getattr(ws,'import')(pdf,ROOT.RooCmdArg()) 
   ws.factory(name+"_norm[%f]"%h.Integral())
   return pdf,h.Integral()
#########################


if opts.reg:
   ## import POI in the ws
   for i in range(nbinsG):
      print "-> importing bin","r_Bin%d[1,-1,20]"%i
      ws.factory("r_Bin%d[1,-1,20]"%i)
   #ws.factory("delta[1.0]") # const, sigma=2/sqrt(delta), INF= NO REG, 0 = SUPER REG
   #ws.factory("delta[1e-5]") # const, sigma=2/sqrt(delta), INF= NO REG, 0 = SUPER REG
   #delta=ws.var("delta")
   delta=0.001

   #ws.factory("regIn[0.]")## regularization will pull to this value, const ?
   ##for i in range(nbinsG):
   ##	ws.factory("reg%d[0.,-4,4]"%i)## regularization will pull to this value, const ?
   ##	#regIn=ws.var("regIn")

constrains=[]
for row in range(nbinsG):
   model=ROOT.TH1D("h_Bin%d"%row,"model_Bin%d"%row,nbinsR,0,nbinsR)
   model.SetBinContent(0,0) #underflow
   model.SetBinContent(nbinsR+1,0) #overflow
   for col in range(nbinsR):
      model.SetBinContent(col+1, mhat[ (row,col) ] * gen.GetBinContent(row+1))
   ## import in Ws
   pdf,norm=importPdf(model,"model_Bin%d"%row)
   datacard.write("shapes Bin%d *\t"%row + opts.output + "\t" + opts.workspace+":model_Bin%d\n"%row)
   if not opts.reg or row==0 or row==nbinsG-1:
	   pass
   else:
	###
	#constrains.append( "constr%d"%row + " constr const%(bin)d_In[0.],RooFormulaVar::fconstr%(bin)d(\"r_Bin%(prev)d+r_Bin%(next)d-2*r_Bin%(bin)d\",{r_Bin%(prev)d,r_Bin%(bin)d,r_Bin%(next)d}),constr%(bin)d_S[%(reg)f]" %{"bin":row,"next":row+1,"prev":row-1,"reg":delta} )  ## REG STRENGTH 
	constrains.append( "constr%d"%row + " constr r_Bin%(prev)d+r_Bin%(next)d-2*r_Bin%(bin)d %(reg)f" %{"bin":row,"next":row+1,"prev":row-1,"reg":delta} )  ## REG STRENGTH 

	#constrains.append( "fconstr%d"%row ) 
	#ws.factory("RooFormulaVar::fconstr%(bin)d(\"r_Bin%(prev)d+r_Bin%(next)d-2*r_Bin%(bin)d\",{r_Bin%(prev)d,r_Bin%(next)d,r_Bin%(bin)d})"%{"bin":row,"next":row+1,"prev":row-1})
	#ws.factory("PROD::reg_Bin%(bin)d( model_Bin%(bin)d , Gaussian::constr_Bin%(bin)d( regIn, fconstr%(bin)d , delta ) )"%{"bin":row})
	#f= ws.function("fconstr%d"%row)
	#constr= ROOT.RooGaussian("constr_Bin%d","regularization",f,regIn,delta)
	#reg_pdf=ROOT.RooProdPdf("reg_Bin%d"%row,"reg pdf",pdf,constr)
	#getattr(ws,'import')(reg_pdf)
        #ws.factory("reg_Bin%d"%row+"_norm[%f]"%norm)
   	#datacard.write("shapes Bin%d *\t"%row + opts.output + "\t" + opts.workspace+":reg_Bin%d\n"%row)

importPdf(bkg,"model_Bkg")
datacard.write("shapes Bkg *\t" + opts.output + "\t" + opts.workspace+":model_Bkg\n")

## write out ws
ws.writeToFile(opts.output)

## write full datacard ## square matrix -- only 1 cat
datacard.write("-------------------------------------\n")
datacard.write("bin cat0\n")
datacard.write("observation -1\n")
datacard.write("-------------------------------------\n")

writeBkg=False # if true will be floated

datacard.write("bin")
for i in range(nbinsG):
   datacard.write("\tcat0")
if writeBkg:
	datacard.write("\tcat0") #bkg
datacard.write("\n")

datacard.write("process")
for i in range(nbinsG):
   datacard.write("\tBin%d"%i)
if writeBkg:
	datacard.write("\tBkg") #bkg
datacard.write("\n")

datacard.write("process")
for i in range(nbinsG):
   datacard.write("\t-%d"%i)
if writeBkg:
	datacard.write("\t%d"%nbinsG) #bkg
datacard.write("\n")

datacard.write("rate")
for i in range(nbinsG):
   datacard.write("\t1")
if writeBkg:
	datacard.write("\t1") #bkg
datacard.write("\n")

datacard.write("-------------------------------------\n")
# syst
datacard.write("# NO SYST\n")
for s in constrains:
	datacard.write("%s\n"%s)
datacard.write("-------------------------------------\n")

datacard.close()

import time
time.sleep(1)
cmd="text2workspace.py " + opts.datacard + " -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel "

for i in range(nbinsG):
   cmd+= " --PO map='.*Bin%d.*':r_Bin%d[1,-1,20]"%(i,i)
cmd+= " --PO map='.*Bkg.*':r_Bkg[1]"

cmd += " --X-allow-no-background "

print "Calling command:\n",cmd
st=call(cmd,shell=True)
if st==0: print "Datacard: Success"
else: print "Datacard: Error"

time.sleep(1)
import re
cmd= "combine -M MultiDimFit " + re.sub('.txt','.root',opts.datacard ) 
#cmd += " -S 0 "
cmd += " --floatOtherPOIs=1 "
cmd += " --saveWorkspace "

print "Run with:\n",cmd
st=call(cmd,shell=True)
if st==0: print "Fit: Success"
else: print "Fit: Error"

del ws

