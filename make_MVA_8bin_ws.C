// Config that uses 8 bins (Njets=7,...,14), and uses N7 for
//  the bases of calculating the other Njets bins.
// tt Njets shape systematic is included.
// Aron Soha
// May 7, 2018
// Sep 19, 2018 : This version now uses the new DeepESM descriminant (instead of the Fisher)
// Nov 2, 2018 : This version switches to using RooParametricHist for ttbar
// Nov 13, 2018 : Added helper function to construct RooArgList of RooFormulaVar's for tt bkg,
//                and feature to allow highest njets bin to be inclusive up to a set number.
// Nov 14, 2018 : Add multiplicative ttbar shape systematics with log normal nuisance parameters.
// Dec 12, 2018 : Switched to new fit function for tt background (solves correlation issue).

#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooFormulaVar.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "interface/RooParametricHist.h"

#include "TH1.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTreeFormulaManager.h"

#include <iostream>
#include <sstream>
#include <string>

class NuisanceParam
{
public:
    const RooAbsArg& r_name;
    const TH1D* h_r;
    const TH1D* h_rprime;
    RooArgList rprime_names;
    
    NuisanceParam(const RooAbsArg& r_name, const TH1D* h_r, const TH1D* h_rprime, RooArgList rprime_names) 
        : r_name(r_name)
        , h_r(h_r)
        , h_rprime(h_rprime)
        , rprime_names(rprime_names)
    {        
    }

    NuisanceParam(const RooAbsArg& r_name, const TH1D* h_r) 
        : r_name(r_name)
        , h_r(h_r)
        , h_rprime(nullptr)
    {        
    }
};

Double_t step(double_t x) {
  return 1;
}

void addNPs(stringstream& f, RooArgList& list, const double r, const RooAbsArg& NP)
{
    std::string paramNum = std::to_string(list.getSize());
    f << "*TMath::Power(" << r << ",@"+paramNum+")";
    list.add(NP);
}

void addNPs(stringstream& f, RooArgList& list, const double r, const RooAbsArg& NP, const double rprime, const RooAbsArg& NPg)
{
    std::string paramNum  = std::to_string(list.getSize());
    std::string paramNumG = std::to_string(list.getSize()+1);
    f << "*TMath::Power(" << r << ",@"+paramNum+")*TMath::Power(" << rprime << ",@"+paramNumG+"*@"+paramNum+")";
    list.add(NP);
    list.add(NPg);
}

template<typename T> void WriteHisto2WS(TFile* f, const std::string& histName, const std::string& wsName, const std::vector<std::string>& MVABins, const std::string& type = "")
{
    for(const auto& bin : MVABins)
    {
        T* h = (T*)f->Get( (bin+"_"+histName).c_str() );
        std::string name = bin+"_"+wsName;
        if( type != "") 
            name = bin+"_"+type+"_mcStat"+bin+type+"Bin"+wsName;
        h->SetName( name.c_str() );
        h->Write();
    }
}

void construct_formula(std::string procName, RooArgList& binlist, const RooArgList& paramlist, const std::vector<NuisanceParam>& NPs) 
{
  // Functional form:
  // f(x) = Njets bin x / Njets bin x-1 = a2 + [ (a1-a2)^(x-a0_val) / (a0-a2)^(x-a1_val) ]^(1/(a1_val-a0_val)) where a1 > a2,
  //   where x = 0 corresponds to 8.
  // a0 = Ratio of Njets=8 to Njets=7
  // a1 = Ratio of Njets=10 to Njets=9
  // a2 = Asymptotic value as Njets-->Inf

  // In terms of njets instead of ratio:
  // F(0) = N7
  // F(1) = N7*f(0)
  // F(2) = F(1)*f(1) = N7*f(0)*f(1)
  //     ...
  // N7 = Njets=7

  std::cout << "size of NPs: " << NPs.size() << std::endl;

  ROOT::v5::TFormula::SetMaxima(10000);

  int max_bin = 18; // 14 means just njets=14, 20 means last bin is inclusive up through njets=20

  // Will update this to go only to 12 jets, rather than 14. So 6 bins instead of 8
  for (int i=0; i<6; i++) 
  {
    std::stringstream form;
    RooArgList formArgList;

    form << "(@0";
    formArgList.add(paramlist[0]); // N7_tt for this MVA bin

    if (i>=1) { // for bin 1 and up
      for (int j=0; j<i; j++) {
	//form << "*(@3 + TMath::Power( TMath::Power( @2-@3 , "<<j<<"-0 ) / TMath::Power( @1-@3 , "<<j<<"-2 ) , 1/(2-0) ))";
	//form << "*(@3 + TMath::Power( TMath::Power( @2-@3 , "<<j<<" ) / TMath::Power( @1-@3 , "<<j-2<<" ) , 1/2 )* ( ( @1-@2 > 0) ? 1.0 : 0.0 ) * ((@2-@3 > 0) ? 1.0 : 0.0) )";
	//form << "*(@3 + TMath::Power( TMath::Power( @2-@3 , "<<j<<" ) / TMath::Power( @1-@3 , "<<j-2<<" ) , 1/2 )* ( ((@1-@2<0)(0.0)+(@1-@2>=0)(1.0))*((@2-@3<0)(0.0)+(@2-@3>=0)(1.0)) ) )";
	//form << "*(@3 + TMath::Power( TMath::Power( @2-@3 , "<<j<<" ) / TMath::Power( @1-@3 , "<<j-2<<" ) , 1/2 )* step(1.0))";
	//form << "*(@3 + (@1-@2>0)*(@2-@3>0)*TMath::Power( TMath::Power( @2-@3 , "<<j<<" ) / TMath::Power( @1-@3 , "<<j-2<<" ) , 1/2 ) )";
	//form << "*(@3 + TMath::Power( TMath::Power( @2-@3 , "<<j<<" ) / TMath::Power( @1-@3 , "<<j-2<<" ) , 1/2 ))";
	//form << "*(@3 + (@1-@2>0)*(@2-@3>0)*TMath::Power( TMath::Abs( TMath::Power( @2-@3 , "<<j<<" ) / TMath::Power( @1-@3 , "<<j-2<<" ) ) , 1/2 ))";
	//form << "*(@3 + TMath::Power( TMath::Abs( TMath::Power( @2-@3 , "<<j<<" ) / TMath::Power( @1-@3 , "<<j-2<<" ) ) , 1/2 ))";
	form << "*( @3>1 ? ( @2 - 1./@3 + TMath::Power(  TMath::Power( 1./@3 , "<<j<<" ) / TMath::Power( @1-@2+1./@3 , "<<j-2<<" ) , 1/2 ) ) : ( @2 - (2.-@3) + TMath::Power( TMath::Power( (2.-@3) , "<<j<<" ) / TMath::Power( @1-@2+(2.-@3) , "<<j-2<<" ) , 1/2 ) ) )";
      }
      formArgList.add(paramlist[1]); // a0_tt
      formArgList.add(paramlist[2]); // a1_tt
      formArgList.add(paramlist[3]); // d_tt
    } // end bin 1 and up

    // The last bin covers from njet=14 through njet=max_bin
    if (i==5) {
      for (int k=6; k<=max_bin-7; k++) {
	form << " + 1";
	for (int j=0; j<k; j++) {
	  form << "*( @3>1 ? ( @2 - 1./@3 + TMath::Power(  TMath::Power( 1./@3 , "<<j<<" ) / TMath::Power( @1-@2+1./@3 , "<<j-2<<" ) , 1/2 ) ) : ( @2 - (2.-@3) + TMath::Power( TMath::Power( (2.-@3) , "<<j<<" ) / TMath::Power( @1-@2+(2.-@3) , "<<j-2<<" ) , 1/2 ) ) )";
	}
      }
    }
    form << ")";

    // Add nuisance parameters
    for(unsigned int j = 0; j < NPs.size(); j++)
    {
        double r = NPs[j].h_r->GetBinContent(i+1);
        if(NPs[j].h_rprime)
        {
            double rprime = NPs[j].h_rprime->GetBinContent(i+1);
            addNPs(form, formArgList, r, NPs[j].r_name, rprime, NPs[j].rprime_names[i]);
        }
        else
        {
            addNPs(form, formArgList, r, NPs[j].r_name);
        }
    }

    // Create RooFormulaVar for this bin
    std::stringstream binName;
    binName << procName << "_b" << i;
    RooFormulaVar* binvar = new RooFormulaVar(binName.str().c_str(), "", form.str().c_str(), RooArgList(formArgList));
    binlist.add(*binvar);

    std::cout << "bin i = " << i << " , njets = " << i+7 << std::endl;
    std::cout << "process bin name : " << binName.str().c_str() << std::endl;
    std::cout << "Formula : " << form.str().c_str() << std::endl;
    formArgList.Print();
    std::cout << std::endl;
  }
}

void make_MVA_8bin_ws(const std::string year = "2016", const std::string infile_path = "Keras_V1.2.5_v4", const std::string model = "RPV", const std::string mass = "550", const std::string dataType = "pseudodata", const std::string syst = "", bool shared = true, bool TTonly = false) {
  using namespace RooFit;
  // Load the combine library to get access to the RooParametricHist
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
    
  // Output file and workspace 
  TFile *fOut = new TFile(("MVA_"+year+"_"+model+"_"+mass+"_ws.root").c_str(),"RECREATE");
  RooWorkspace *wspace = new RooWorkspace("wspace","wspace");

  // // njet is our variable, 6 bins, from 7 up through 12,
  // //   Note that njet=12 is inclusive as >=12
  // // D1, D2, D3, D4 are the MVA bins
  // wspace->factory("nj_D1[6.5,12.5]");
  // //wspace->factory("nj_D1[0,6]");
  // wspace->var("nj_D1")->setBins(6);
  // RooArgSet vars_D1(*wspace->var("nj_D1"));

  // wspace->factory("nj_D2[6.5,12.5]");
  // //wspace->factory("nj_D2[0,6]");
  // wspace->var("nj_D2")->setBins(6);
  // RooArgSet vars_D2(*wspace->var("nj_D2"));

  // wspace->factory("nj_D3[6.5,12.5]");
  // //wspace->factory("nj_D3[0,6]");
  // wspace->var("nj_D3")->setBins(6);
  // RooArgSet vars_D3(*wspace->var("nj_D3"));

  // wspace->factory("nj_D4[6.5,12.5]");
  // //wspace->factory("nj_D4[0,6]");
  // wspace->var("nj_D4")->setBins(6);
  // RooArgSet vars_D4(*wspace->var("nj_D4"));

  wspace->factory("CMS_th1x[0,6]");
  wspace->var("CMS_th1x")->setBins(6);
  RooArgSet vars_D1(*wspace->var("CMS_th1x"));
  RooArgSet vars_D2(*wspace->var("CMS_th1x"));
  RooArgSet vars_D3(*wspace->var("CMS_th1x"));
  RooArgSet vars_D4(*wspace->var("CMS_th1x"));
    
  // file for obtaining histograms
  TFile* file = TFile::Open((infile_path+"/njets_for_Aron.root").c_str());

  TH1D* data_th1_D1;
  if     (dataType == "data")        data_th1_D1 = (TH1D*)file->Get(("D1_data_h_njets_pt30_1l"+syst).c_str());  // Actual data -- be careful
  else if(dataType == "pseudodata" ) data_th1_D1 = (TH1D*)file->Get(("D1_pseudodata_h_njets_pt30_1l"+syst).c_str()); // without signal
  else if(dataType == "pseudodataS") data_th1_D1 = (TH1D*)file->Get(("D1_pseudodataS_"+model+"_"+mass+"_h_njets_pt30_1l"+syst).c_str());  // with signal
  else if(dataType == "pseudodataA" ) data_th1_D1 = (TH1D*)file->Get(("D1_pseudodataFunc28_24_236_h_njets_pt30_1l"+syst).c_str()); //
  else if(dataType == "pseudodataB" ) data_th1_D1 = (TH1D*)file->Get(("D1_pseudodataFunc28_24_18_h_njets_pt30_1l"+syst).c_str()); //
  else if(dataType == "pseudodataC" ) data_th1_D1 = (TH1D*)file->Get(("D1_pseudodataFunc28_24_-20_h_njets_pt30_1l"+syst).c_str()); //
  else                               data_th1_D1 = (TH1D*)file->Get(("D1_"+dataType+"_h_njets_pt30_1l"+syst).c_str());
  TH1D* otherMC_th1_D1 = (TH1D*)file->Get(("D1_OTHER_h_njets_pt30_1l"+syst).c_str());
  TH1D* qcdMC_th1_D1 = (TH1D*)file->Get(("D1_QCD_h_njets_pt30_1l"+syst).c_str());
  TH1D* ttxMC_th1_D1 = (TH1D*)file->Get(("D1_TTX_h_njets_pt30_1l"+syst).c_str());
  TH1D* sigMC_th1_D1 = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_h_njets_pt30_1l"+syst).c_str());

  TH1D* data_th1_D2;
  if     (dataType == "data")        data_th1_D2 = (TH1D*)file->Get(("D2_data_h_njets_pt30_1l"+syst).c_str());  // Actual data -- be careful
  else if(dataType == "pseudodata" ) data_th1_D2 = (TH1D*)file->Get(("D2_pseudodata_h_njets_pt30_1l"+syst).c_str()); // without signal
  else if(dataType == "pseudodataS") data_th1_D2 = (TH1D*)file->Get(("D2_pseudodataS_"+model+"_"+mass+"_h_njets_pt30_1l"+syst).c_str());  // with signal
  else if(dataType == "pseudodataA" ) data_th1_D2 = (TH1D*)file->Get(("D2_pseudodataFunc28_24_236_h_njets_pt30_1l"+syst).c_str()); //
  else if(dataType == "pseudodataB" ) data_th1_D2 = (TH1D*)file->Get(("D2_pseudodataFunc28_24_18_h_njets_pt30_1l"+syst).c_str()); //
  else if(dataType == "pseudodataC" ) data_th1_D2 = (TH1D*)file->Get(("D2_pseudodataFunc28_24_-20_h_njets_pt30_1l"+syst).c_str()); //
  else                               data_th1_D2 = (TH1D*)file->Get(("D2_"+dataType+"_h_njets_pt30_1l"+syst).c_str());
  TH1D* otherMC_th1_D2 = (TH1D*)file->Get(("D2_OTHER_h_njets_pt30_1l"+syst).c_str());
  TH1D* qcdMC_th1_D2 = (TH1D*)file->Get(("D2_QCD_h_njets_pt30_1l"+syst).c_str());
  TH1D* ttxMC_th1_D2 = (TH1D*)file->Get(("D2_TTX_h_njets_pt30_1l"+syst).c_str());
  TH1D* sigMC_th1_D2 = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_h_njets_pt30_1l"+syst).c_str());

  TH1D* data_th1_D3;
  if     (dataType == "data")        data_th1_D3 = (TH1D*)file->Get(("D3_data_h_njets_pt30_1l"+syst).c_str());  // Actual data -- be careful
  else if(dataType == "pseudodata" ) data_th1_D3 = (TH1D*)file->Get(("D3_pseudodata_h_njets_pt30_1l"+syst).c_str()); // without signal
  else if(dataType == "pseudodataS") data_th1_D3 = (TH1D*)file->Get(("D3_pseudodataS_"+model+"_"+mass+"_h_njets_pt30_1l"+syst).c_str());  // with signal
  else if(dataType == "pseudodataA" ) data_th1_D3 = (TH1D*)file->Get(("D3_pseudodataFunc28_24_236_h_njets_pt30_1l"+syst).c_str()); //
  else if(dataType == "pseudodataB" ) data_th1_D3 = (TH1D*)file->Get(("D3_pseudodataFunc28_24_18_h_njets_pt30_1l"+syst).c_str()); //
  else if(dataType == "pseudodataC" ) data_th1_D3 = (TH1D*)file->Get(("D3_pseudodataFunc28_24_-20_h_njets_pt30_1l"+syst).c_str()); //
  else                               data_th1_D3 = (TH1D*)file->Get(("D3_"+dataType+"_h_njets_pt30_1l"+syst).c_str());
  TH1D* otherMC_th1_D3 = (TH1D*)file->Get(("D3_OTHER_h_njets_pt30_1l"+syst).c_str());
  TH1D* qcdMC_th1_D3 = (TH1D*)file->Get(("D3_QCD_h_njets_pt30_1l"+syst).c_str());
  TH1D* ttxMC_th1_D3 = (TH1D*)file->Get(("D3_TTX_h_njets_pt30_1l"+syst).c_str());
  TH1D* sigMC_th1_D3 = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_h_njets_pt30_1l"+syst).c_str());

  TH1D* data_th1_D4;
  if     (dataType == "data")        data_th1_D4 = (TH1D*)file->Get(("D4_data_h_njets_pt30_1l"+syst).c_str());  // Actual data -- be careful
  else if(dataType == "pseudodata" ) data_th1_D4 = (TH1D*)file->Get(("D4_pseudodata_h_njets_pt30_1l"+syst).c_str()); // without signal
  else if(dataType == "pseudodataS") data_th1_D4 = (TH1D*)file->Get(("D4_pseudodataS_"+model+"_"+mass+"_h_njets_pt30_1l"+syst).c_str());  // with signal
  else if(dataType == "pseudodataA" ) data_th1_D4 = (TH1D*)file->Get(("D4_pseudodataFunc28_24_236_h_njets_pt30_1l"+syst).c_str()); //
  else if(dataType == "pseudodataB" ) data_th1_D4 = (TH1D*)file->Get(("D4_pseudodataFunc28_24_18_h_njets_pt30_1l"+syst).c_str()); //
  else if(dataType == "pseudodataC" ) data_th1_D4 = (TH1D*)file->Get(("D4_pseudodataFunc28_24_-20_h_njets_pt30_1l"+syst).c_str()); //
  else                               data_th1_D4 = (TH1D*)file->Get(("D4_"+dataType+"_h_njets_pt30_1l"+syst).c_str());
  TH1D* otherMC_th1_D4 = (TH1D*)file->Get(("D4_OTHER_h_njets_pt30_1l"+syst).c_str());
  TH1D* qcdMC_th1_D4 = (TH1D*)file->Get(("D4_QCD_h_njets_pt30_1l"+syst).c_str());
  TH1D* ttxMC_th1_D4 = (TH1D*)file->Get(("D4_TTX_h_njets_pt30_1l"+syst).c_str());
  TH1D* sigMC_th1_D4 = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_h_njets_pt30_1l"+syst).c_str());

  // tt bkg param setup
  //RooRealVar a0_tt("a0_tt","a0 of tt bkg shape",0.28,0.0,1.0);
  //RooRealVar a1_tt("a1_tt","a1 of tt bkg shape",0.24,0.0,1.0);
  ////RooRealVar a2_tt("a2_tt","a2 of tt bkg shape",0.10,-0.5,0.5);
  //RooRealVar a2_tt("a2_tt","a2 of tt bkg shape",0.10,-0.5,0.23);

  //RooRealVar a0_tt("a0_tt","a0 of tt bkg shape",0.28,0.26,0.30);
  //RooRealVar a1_tt("a1_tt","a1 of tt bkg shape",0.24,0.23,0.27);
  //RooRealVar a2_tt("a2_tt","a2 of tt bkg shape",0.10,-0.05,0.23);

  //RooRealVar a0_tt("a0_tt","a0 of tt bkg shape",0.28,0.26,0.50);
  //RooRealVar a1_tt("a1_tt","a1 of tt bkg shape",0.25,0.24,0.50);
  //RooRealVar a2_tt("a2_tt","a2 of tt bkg shape",0.10,-0.05,0.24);

  RooRealVar a0_tt(("a0_tt_"+year).c_str(),("a0 of tt bkg shape for "+year).c_str(),0.28,0.1,0.4);
  RooRealVar a1_tt(("a1_tt_"+year).c_str(),("a1 of tt bkg shape for "+year).c_str(),0.24,0.1,0.4);
  //RooRealVar a2_tt(("a2_tt_"+year).c_str(),("a2 of tt bkg shape for "+year).c_str(),0.20,0.0,0.4);
  RooRealVar d_tt(("d_tt_"+year).c_str(),("d of tt bkg shape for "+year).c_str(),20,-2500,2500);

  // also make the separate paramters per MVA bin for use in different fit setup 
  RooRealVar a0_tt_D1(("a0_tt_D1_"+year).c_str(),("a0 of tt bkg shape D1 for "+year).c_str(),0.28,0.0,1.0);
  RooRealVar a1_tt_D1(("a1_tt_D1_"+year).c_str(),("a1 of tt bkg shape D1 for "+year).c_str(),0.24,0.0,1.0);
  //RooRealVar a2_tt_D1(("a2_tt_D1_"+year).c_str(),("a2 of tt bkg shape D1 for "+year).c_str(),0.10,-0.5,0.5);
  RooRealVar  d_tt_D1(("d_tt_D1_" +year).c_str(),("d  of tt bkg shape D1 for "+year).c_str(),20,-2500,2500);

  RooRealVar a0_tt_D2(("a0_tt_D2_"+year).c_str(),("a0 of tt bkg shape D2 for "+year).c_str(),0.28,0.0,1.0);
  RooRealVar a1_tt_D2(("a1_tt_D2_"+year).c_str(),("a1 of tt bkg shape D2 for "+year).c_str(),0.24,0.0,1.0);
  //RooRealVar a2_tt_D2(("a2_tt_D2_"+year).c_str(),("a2 of tt bkg shape D2 for "+year).c_str(),0.10,-0.5,0.5);
  RooRealVar  d_tt_D2(("d_tt_D2_" +year).c_str(),("d  of tt bkg shape D2 for "+year).c_str(),20,-2500,2500);

  RooRealVar a0_tt_D3(("a0_tt_D3_"+year).c_str(),("a0 of tt bkg shape D3 for "+year).c_str(),0.28,0.0,1.0);
  RooRealVar a1_tt_D3(("a1_tt_D3_"+year).c_str(),("a1 of tt bkg shape D3 for "+year).c_str(),0.24,0.0,1.0);
  //RooRealVar a2_tt_D3(("a2_tt_D3_"+year).c_str(),("a2 of tt bkg shape D3 for "+year).c_str(),0.10,-0.5,0.5);
  RooRealVar  d_tt_D3(("d_tt_D3_" +year).c_str(),("d  of tt bkg shape D3 for "+year).c_str(),20,-2500,2500);

  RooRealVar a0_tt_D4(("a0_tt_D4_"+year).c_str(),("a0 of tt bkg shape D4 for "+year).c_str(),0.28,0.0,1.0);
  RooRealVar a1_tt_D4(("a1_tt_D4_"+year).c_str(),("a1 of tt bkg shape D4 for "+year).c_str(),0.24,0.0,1.0);
  //RooRealVar a2_tt_D4(("a2_tt_D4_"+year).c_str(),("a2 of tt bkg shape D4 for "+year).c_str(),0.10,-0.5,0.5);
  RooRealVar  d_tt_D4(("d_tt_D4_" +year).c_str(),("d  of tt bkg shape D4 for "+year).c_str(),20,-2500,2500);

  double_t n7_tt_portion_D1 = TTonly ? data_th1_D1->GetBinContent(1) : data_th1_D1->GetBinContent(1) - otherMC_th1_D1->GetBinContent(1) - qcdMC_th1_D1->GetBinContent(1) - ttxMC_th1_D1->GetBinContent(1);
  double_t n7_tt_portion_D1_low = n7_tt_portion_D1-6000;
  if (n7_tt_portion_D1_low<0) n7_tt_portion_D1_low = 0;
  RooRealVar N7_tt_D1(("N7_tt_D1_"+year).c_str(),"njets 7 for tt bkg in MVA D1",n7_tt_portion_D1,n7_tt_portion_D1_low,n7_tt_portion_D1+6000);

  double_t n7_tt_portion_D2 = TTonly ? data_th1_D2->GetBinContent(1) : data_th1_D2->GetBinContent(1) - otherMC_th1_D2->GetBinContent(1) - qcdMC_th1_D2->GetBinContent(1) - ttxMC_th1_D2->GetBinContent(1);
  double_t n7_tt_portion_D2_low = n7_tt_portion_D2-5000;
  if (n7_tt_portion_D2_low<0) n7_tt_portion_D2_low = 0;
  RooRealVar N7_tt_D2(("N7_tt_D2_"+year).c_str(),"njets 7 for tt bkg in MVA D2",n7_tt_portion_D2,n7_tt_portion_D2_low,n7_tt_portion_D2+5000);

  double_t n7_tt_portion_D3 = TTonly ? data_th1_D3->GetBinContent(1) : data_th1_D3->GetBinContent(1) - otherMC_th1_D3->GetBinContent(1) - qcdMC_th1_D3->GetBinContent(1) - ttxMC_th1_D3->GetBinContent(1);
  double_t n7_tt_portion_D3_low = n7_tt_portion_D3-4000;
  if (n7_tt_portion_D3_low<0) n7_tt_portion_D3_low = 0;
  RooRealVar N7_tt_D3(("N7_tt_D3_"+year).c_str(),"njets 7 for tt bkg in MVA D3",n7_tt_portion_D3,n7_tt_portion_D3_low,n7_tt_portion_D3+4000);

  double_t n7_tt_portion_D4 = TTonly ? data_th1_D4->GetBinContent(1) : data_th1_D4->GetBinContent(1) - otherMC_th1_D4->GetBinContent(1) - qcdMC_th1_D4->GetBinContent(1) - ttxMC_th1_D4->GetBinContent(1);
  double_t n7_tt_portion_D4_low = n7_tt_portion_D4-3000;
  if (n7_tt_portion_D4_low<0) n7_tt_portion_D4_low = 0;
  RooRealVar N7_tt_D4(("N7_tt_D4_"+year).c_str(),"njets 7 for tt bkg in MVA D4",n7_tt_portion_D4,n7_tt_portion_D4_low,n7_tt_portion_D4+3000);

  // tt shape systematic nuisance parameters
  wspace->factory(("np_tt_JECUp_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_JERUp_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_btg_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_lep_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_nom_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCR_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory("np_tt_pdf[0.0]"); // fully correlated
  wspace->factory("np_tt_FSR[0.0]"); // fully correlated
  wspace->factory("np_tt_ISR[0.0]"); // fully correlated
  wspace->factory("np_tt_scl[0.0]"); // fully correlated
  wspace->factory(("np_tt_ht_"+year+"[0.0]").c_str());// uncorrelated
  wspace->factory(("np_tt_httail_"+year+"[0.0]").c_str());// uncorrelated
  wspace->factory(("np_tt_htnjet_"+year+"[0.0]").c_str());// uncorrelated
  wspace->factory(("np_tt_pu_"+year+"[0.0]").c_str());// uncorrelated
  wspace->factory(("np_tt_JECDown_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_JERDown_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD1Bin1_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD1Bin2_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD1Bin3_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD1Bin4_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD1Bin5_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD1Bin6_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD2Bin1_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD2Bin2_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD2Bin3_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD2Bin4_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD2Bin5_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD2Bin6_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD3Bin1_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD3Bin2_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD3Bin3_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD3Bin4_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD3Bin5_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD3Bin6_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD4Bin1_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD4Bin2_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD4Bin3_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD4Bin4_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD4Bin5_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_qcdCRErrD4Bin6_"+year+"[0.0]").c_str()); // uncorrelated

  // Load in the histograms with the bin-by-bin ratios to be used in the ttbar shape systematics
  TFile* tt_syst_file = TFile::Open((infile_path+"/ttbar_systematics.root").c_str());

  TH1D* tt_syst_JECUp_D1 = (TH1D*)tt_syst_file->Get("D1_JECUp");
  TH1D* tt_syst_JECUp_D2 = (TH1D*)tt_syst_file->Get("D2_JECUp");
  TH1D* tt_syst_JECUp_D3 = (TH1D*)tt_syst_file->Get("D3_JECUp");
  TH1D* tt_syst_JECUp_D4 = (TH1D*)tt_syst_file->Get("D4_JECUp");

  TH1D* tt_syst_JECDown_D1 = (TH1D*)tt_syst_file->Get("D1_JECDown");
  TH1D* tt_syst_JECDown_D2 = (TH1D*)tt_syst_file->Get("D2_JECDown");
  TH1D* tt_syst_JECDown_D3 = (TH1D*)tt_syst_file->Get("D3_JECDown");
  TH1D* tt_syst_JECDown_D4 = (TH1D*)tt_syst_file->Get("D4_JECDown");

  TH1D* tt_syst_JERUp_D1 = (TH1D*)tt_syst_file->Get("D1_JERUp");
  TH1D* tt_syst_JERUp_D2 = (TH1D*)tt_syst_file->Get("D2_JERUp");
  TH1D* tt_syst_JERUp_D3 = (TH1D*)tt_syst_file->Get("D3_JERUp");
  TH1D* tt_syst_JERUp_D4 = (TH1D*)tt_syst_file->Get("D4_JERUp");

  TH1D* tt_syst_JERDown_D1 = (TH1D*)tt_syst_file->Get("D1_JERDown");
  TH1D* tt_syst_JERDown_D2 = (TH1D*)tt_syst_file->Get("D2_JERDown");
  TH1D* tt_syst_JERDown_D3 = (TH1D*)tt_syst_file->Get("D3_JERDown");
  TH1D* tt_syst_JERDown_D4 = (TH1D*)tt_syst_file->Get("D4_JERDown");

  TH1D* tt_syst_btg_D1 = (TH1D*)tt_syst_file->Get("D1_btg");
  TH1D* tt_syst_btg_D2 = (TH1D*)tt_syst_file->Get("D2_btg");
  TH1D* tt_syst_btg_D3 = (TH1D*)tt_syst_file->Get("D3_btg");
  TH1D* tt_syst_btg_D4 = (TH1D*)tt_syst_file->Get("D4_btg");

  TH1D* tt_syst_lep_D1 = (TH1D*)tt_syst_file->Get("D1_lep");
  TH1D* tt_syst_lep_D2 = (TH1D*)tt_syst_file->Get("D2_lep");
  TH1D* tt_syst_lep_D3 = (TH1D*)tt_syst_file->Get("D3_lep");
  TH1D* tt_syst_lep_D4 = (TH1D*)tt_syst_file->Get("D4_lep");

  TH1D* tt_syst_nom_D1 = (TH1D*)tt_syst_file->Get("D1_nom");
  TH1D* tt_syst_nom_D2 = (TH1D*)tt_syst_file->Get("D2_nom");
  TH1D* tt_syst_nom_D3 = (TH1D*)tt_syst_file->Get("D3_nom");
  TH1D* tt_syst_nom_D4 = (TH1D*)tt_syst_file->Get("D4_nom");

  TH1D* tt_syst_qcdCR_D1 = (TH1D*)tt_syst_file->Get("D1_qcdCR");
  TH1D* tt_syst_qcdCR_D2 = (TH1D*)tt_syst_file->Get("D2_qcdCR");
  TH1D* tt_syst_qcdCR_D3 = (TH1D*)tt_syst_file->Get("D3_qcdCR");
  TH1D* tt_syst_qcdCR_D4 = (TH1D*)tt_syst_file->Get("D4_qcdCR");

  TH1D* tt_syst_pdf_D1 = (TH1D*)tt_syst_file->Get("D1_pdf");
  TH1D* tt_syst_pdf_D2 = (TH1D*)tt_syst_file->Get("D2_pdf");
  TH1D* tt_syst_pdf_D3 = (TH1D*)tt_syst_file->Get("D3_pdf");
  TH1D* tt_syst_pdf_D4 = (TH1D*)tt_syst_file->Get("D4_pdf");

  TH1D* tt_syst_FSR_D1 = (TH1D*)tt_syst_file->Get("D1_FSR");
  TH1D* tt_syst_FSR_D2 = (TH1D*)tt_syst_file->Get("D2_FSR");
  TH1D* tt_syst_FSR_D3 = (TH1D*)tt_syst_file->Get("D3_FSR");
  TH1D* tt_syst_FSR_D4 = (TH1D*)tt_syst_file->Get("D4_FSR");

  TH1D* tt_syst_ISR_D1 = (TH1D*)tt_syst_file->Get("D1_ISR");
  TH1D* tt_syst_ISR_D2 = (TH1D*)tt_syst_file->Get("D2_ISR");
  TH1D* tt_syst_ISR_D3 = (TH1D*)tt_syst_file->Get("D3_ISR");
  TH1D* tt_syst_ISR_D4 = (TH1D*)tt_syst_file->Get("D4_ISR");

  TH1D* tt_syst_ht_D1 = (TH1D*)tt_syst_file->Get("D1_ht");
  TH1D* tt_syst_ht_D2 = (TH1D*)tt_syst_file->Get("D2_ht");
  TH1D* tt_syst_ht_D3 = (TH1D*)tt_syst_file->Get("D3_ht");
  TH1D* tt_syst_ht_D4 = (TH1D*)tt_syst_file->Get("D4_ht");

  TH1D* tt_syst_scl_D1 = (TH1D*)tt_syst_file->Get("D1_scl");
  TH1D* tt_syst_scl_D2 = (TH1D*)tt_syst_file->Get("D2_scl");
  TH1D* tt_syst_scl_D3 = (TH1D*)tt_syst_file->Get("D3_scl");
  TH1D* tt_syst_scl_D4 = (TH1D*)tt_syst_file->Get("D4_scl");
  
  TH1D* tt_syst_httail_D1 = (TH1D*)tt_syst_file->Get("D1_httail");
  TH1D* tt_syst_httail_D2 = (TH1D*)tt_syst_file->Get("D2_httail");
  TH1D* tt_syst_httail_D3 = (TH1D*)tt_syst_file->Get("D3_httail");
  TH1D* tt_syst_httail_D4 = (TH1D*)tt_syst_file->Get("D4_httail");

  TH1D* tt_syst_htnjet_D1 = (TH1D*)tt_syst_file->Get("D1_htnjet");
  TH1D* tt_syst_htnjet_D2 = (TH1D*)tt_syst_file->Get("D2_htnjet");
  TH1D* tt_syst_htnjet_D3 = (TH1D*)tt_syst_file->Get("D3_htnjet");
  TH1D* tt_syst_htnjet_D4 = (TH1D*)tt_syst_file->Get("D4_htnjet");
  
  TH1D* tt_syst_pu_D1 = (TH1D*)tt_syst_file->Get("D1_pu");
  TH1D* tt_syst_pu_D2 = (TH1D*)tt_syst_file->Get("D2_pu");
  TH1D* tt_syst_pu_D3 = (TH1D*)tt_syst_file->Get("D3_pu");
  TH1D* tt_syst_pu_D4 = (TH1D*)tt_syst_file->Get("D4_pu");

  TH1D* tt_qcdCRErr_D1 = (TH1D*)tt_syst_file->Get("D1_qcdCRErr");
  TH1D* tt_qcdCRErr_D2 = (TH1D*)tt_syst_file->Get("D2_qcdCRErr");
  TH1D* tt_qcdCRErr_D3 = (TH1D*)tt_syst_file->Get("D3_qcdCRErr");
  TH1D* tt_qcdCRErr_D4 = (TH1D*)tt_syst_file->Get("D4_qcdCRErr");

  // ----------------------  MVA bin 1  ------------------

  // Dataset with 8 bins -- now 6
  RooDataHist data_hist_D1("data_obs_D1","Data observed in MVA bin 1",vars_D1,data_th1_D1);
  wspace->import(data_hist_D1);

  //list of nuisance parameters for tt bkg D1
  const std::vector<NuisanceParam>& nuisanceParams_D1 = {
      {*wspace->var(("np_tt_JECUp_"+year).c_str()), tt_syst_JECUp_D1},
      {*wspace->var(("np_tt_JERUp_"+year).c_str()), tt_syst_JERUp_D1},
      {*wspace->var(("np_tt_btg_"+year).c_str()),   tt_syst_btg_D1},
      {*wspace->var(("np_tt_lep_"+year).c_str()),   tt_syst_lep_D1},
      {*wspace->var(("np_tt_nom_"+year).c_str()),   tt_syst_nom_D1},
      {
          *wspace->var(("np_tt_qcdCR_"+year).c_str()),   tt_syst_qcdCR_D1, tt_qcdCRErr_D1, 
          { 
              *wspace->var(("np_tt_qcdCRErrD1Bin1_"+year).c_str()),
              *wspace->var(("np_tt_qcdCRErrD1Bin2_"+year).c_str()),
              *wspace->var(("np_tt_qcdCRErrD1Bin3_"+year).c_str()),
              *wspace->var(("np_tt_qcdCRErrD1Bin4_"+year).c_str()),
              *wspace->var(("np_tt_qcdCRErrD1Bin5_"+year).c_str()),
              *wspace->var(("np_tt_qcdCRErrD1Bin6_"+year).c_str()),
          }
      },
      {*wspace->var("np_tt_pdf"),                     tt_syst_pdf_D1},
      {*wspace->var("np_tt_FSR"),                     tt_syst_FSR_D1},
      {*wspace->var("np_tt_ISR"),                     tt_syst_ISR_D1},
      {*wspace->var("np_tt_scl"),                     tt_syst_scl_D1},
      {*wspace->var(("np_tt_ht_"+year).c_str()),      tt_syst_ht_D1},
      {*wspace->var(("np_tt_httail_"+year).c_str()),  tt_syst_httail_D1},
      {*wspace->var(("np_tt_htnjet_"+year).c_str()),  tt_syst_htnjet_D1},
      {*wspace->var(("np_tt_pu_"+year).c_str()),      tt_syst_pu_D1},
      {*wspace->var(("np_tt_JECDown_"+year).c_str()), tt_syst_JECDown_D1},
      {*wspace->var(("np_tt_JERDown_"+year).c_str()), tt_syst_JERDown_D1}
  };
  
  //std::cout << "test" << std::endl;
  RooArgList *bkg_tt_bins_D1 = new RooArgList();
  std::string procName_D1 = "background_tt_D1_"+year;
  if (shared) 
  {
    //std::cout << "test shared" << std::endl;
    RooArgList parlist_D1(N7_tt_D1,a0_tt,a1_tt,d_tt);  // list of shape parameters for tt bkg
    //RooArgList parlist_D1(N7_tt_D1,a0_tt,a1_tt,a2_tt);  // list of shape parameters for tt bkg                                
    construct_formula(procName_D1,*bkg_tt_bins_D1,parlist_D1,nuisanceParams_D1);
    std::cout << "after constructing formula" << std::endl;
  } else 
  {
    RooArgList parlist_D1(N7_tt_D1,a0_tt_D1,a1_tt_D1,d_tt_D1);  // list of shape parameters for tt bkg                       
    //RooArgList parlist_D1(N7_tt_D1,a0_tt_D1,a1_tt_D1,a2_tt_D1);  // list of shape parameters for tt bkg
    construct_formula(procName_D1,*bkg_tt_bins_D1,parlist_D1,nuisanceParams_D1);
  }
  std::cout << "test" << std::endl;
  RooParametricHist background_tt_D1(procName_D1.c_str(),"",*wspace->var("CMS_th1x"),*bkg_tt_bins_D1,*data_th1_D1);
  wspace->import(background_tt_D1,RooFit::RecycleConflictNodes());
  std::stringstream procNameD1Norm;
  procNameD1Norm << procName_D1 << "_norm";
  RooAddition tt_norm_D1(procNameD1Norm.str().c_str(),"",*bkg_tt_bins_D1);
  wspace->import(tt_norm_D1,RooFit::RecycleConflictNodes());

  // ---------------------- MVA bin 2  ------------------
  
  // Dataset with 8 bins
  RooDataHist data_hist_D2("data_obs_D2","Data observed in MVA bin 2",vars_D2,data_th1_D2);
  wspace->import(data_hist_D2);
  
  //list of nuisance parameters for tt bkg D2
  const std::vector<NuisanceParam>& nuisanceParams_D2 = {
      {*wspace->var(("np_tt_JECUp_"+year).c_str()), tt_syst_JECUp_D2},
      {*wspace->var(("np_tt_JERUp_"+year).c_str()), tt_syst_JERUp_D2},
      {*wspace->var(("np_tt_btg_"+year).c_str()),   tt_syst_btg_D2},
      {*wspace->var(("np_tt_lep_"+year).c_str()),   tt_syst_lep_D2},
      {*wspace->var(("np_tt_nom_"+year).c_str()),   tt_syst_nom_D2},
      {
          *wspace->var(("np_tt_qcdCR_"+year).c_str()),   tt_syst_qcdCR_D2, tt_qcdCRErr_D2, 
          { 
              *wspace->var(("np_tt_qcdCRErrD2Bin1_"+year).c_str()),
              *wspace->var(("np_tt_qcdCRErrD2Bin2_"+year).c_str()),
              *wspace->var(("np_tt_qcdCRErrD2Bin3_"+year).c_str()),
              *wspace->var(("np_tt_qcdCRErrD2Bin4_"+year).c_str()),
              *wspace->var(("np_tt_qcdCRErrD2Bin5_"+year).c_str()),
              *wspace->var(("np_tt_qcdCRErrD2Bin6_"+year).c_str()),
          }
      },
      {*wspace->var("np_tt_pdf"),                     tt_syst_pdf_D2},
      {*wspace->var("np_tt_FSR"),                     tt_syst_FSR_D2},
      {*wspace->var("np_tt_ISR"),                     tt_syst_ISR_D2},
      {*wspace->var("np_tt_scl"),                     tt_syst_scl_D2},
      {*wspace->var(("np_tt_ht_"+year).c_str()),      tt_syst_ht_D2},
      {*wspace->var(("np_tt_httail_"+year).c_str()),  tt_syst_httail_D2},
      {*wspace->var(("np_tt_htnjet_"+year).c_str()),  tt_syst_htnjet_D2},
      {*wspace->var(("np_tt_pu_"+year).c_str()),      tt_syst_pu_D2},
      {*wspace->var(("np_tt_JECDown_"+year).c_str()), tt_syst_JECDown_D2},
      {*wspace->var(("np_tt_JERDown_"+year).c_str()), tt_syst_JERDown_D2}
  };
  
  RooArgList *bkg_tt_bins_D2 = new RooArgList();
  std::string procName_D2 = "background_tt_D2_"+year;
  if (shared) 
  {
    RooArgList parlist_D2(N7_tt_D2,a0_tt,a1_tt,d_tt);  // list of shape parameters for tt bkg
    //RooArgList parlist_D2(N7_tt_D2,a0_tt,a1_tt,a2_tt);  // list of shape parameters for tt bkg                                
    construct_formula(procName_D2,*bkg_tt_bins_D2,parlist_D2,nuisanceParams_D2);
  } else 
  {
    RooArgList parlist_D2(N7_tt_D2,a0_tt_D2,a1_tt_D2,d_tt_D2);  // list of shape parameters for tt bkg                       
    //RooArgList parlist_D2(N7_tt_D2,a0_tt_D2,a1_tt_D2,a2_tt_D2);  // list of shape parameters for tt bkg
    construct_formula(procName_D2,*bkg_tt_bins_D2,parlist_D2,nuisanceParams_D2);
  }
  RooParametricHist background_tt_D2(procName_D2.c_str(),"",*wspace->var("CMS_th1x"),*bkg_tt_bins_D2,*data_th1_D2);
  wspace->import(background_tt_D2,RooFit::RecycleConflictNodes());
  std::stringstream procNameD2Norm;
  procNameD2Norm << procName_D2 << "_norm";
  RooAddition tt_norm_D2(procNameD2Norm.str().c_str(),"",*bkg_tt_bins_D2);
  wspace->import(tt_norm_D2,RooFit::RecycleConflictNodes());
  
  // ---------------------- MVA bin 3  ------------------
  
  // Dataset with 8 bins
  RooDataHist data_hist_D3("data_obs_D3","Data observed in MVA bin 3",vars_D3,data_th1_D3);
  wspace->import(data_hist_D3);
  
  //list of nuisance parameters for tt bkg D3
  const std::vector<NuisanceParam>& nuisanceParams_D3 = {
      {*wspace->var(("np_tt_JECUp_"+year).c_str()), tt_syst_JECUp_D3},
      {*wspace->var(("np_tt_JERUp_"+year).c_str()), tt_syst_JERUp_D3},
      {*wspace->var(("np_tt_btg_"+year).c_str()),   tt_syst_btg_D3},
      {*wspace->var(("np_tt_lep_"+year).c_str()),   tt_syst_lep_D3},
      {*wspace->var(("np_tt_nom_"+year).c_str()),   tt_syst_nom_D3},
      {
          *wspace->var(("np_tt_qcdCR_"+year).c_str()),   tt_syst_qcdCR_D3, tt_qcdCRErr_D3, 
          { 
              *wspace->var(("np_tt_qcdCRErrD3Bin1_"+year).c_str()),
              *wspace->var(("np_tt_qcdCRErrD3Bin2_"+year).c_str()),
              *wspace->var(("np_tt_qcdCRErrD3Bin3_"+year).c_str()),
              *wspace->var(("np_tt_qcdCRErrD3Bin4_"+year).c_str()),
              *wspace->var(("np_tt_qcdCRErrD3Bin5_"+year).c_str()),
              *wspace->var(("np_tt_qcdCRErrD3Bin6_"+year).c_str()),
          }
      },
      {*wspace->var("np_tt_pdf"),                     tt_syst_pdf_D3},
      {*wspace->var("np_tt_FSR"),                     tt_syst_FSR_D3},
      {*wspace->var("np_tt_ISR"),                     tt_syst_ISR_D3},
      {*wspace->var("np_tt_scl"),                     tt_syst_scl_D3},
      {*wspace->var(("np_tt_ht_"+year).c_str()),      tt_syst_ht_D3},
      {*wspace->var(("np_tt_httail_"+year).c_str()),  tt_syst_httail_D3},
      {*wspace->var(("np_tt_htnjet_"+year).c_str()),  tt_syst_htnjet_D3},
      {*wspace->var(("np_tt_pu_"+year).c_str()),      tt_syst_pu_D3},
      {*wspace->var(("np_tt_JECDown_"+year).c_str()), tt_syst_JECDown_D3},
      {*wspace->var(("np_tt_JERDown_"+year).c_str()), tt_syst_JERDown_D3}
  };

  RooArgList *bkg_tt_bins_D3 = new RooArgList();
  std::string procName_D3 = "background_tt_D3_"+year;
  if (shared) 
  {
    RooArgList parlist_D3(N7_tt_D3,a0_tt,a1_tt,d_tt);  // list of shape parameters for tt bkg
    //RooArgList parlist_D3(N7_tt_D3,a0_tt,a1_tt,a2_tt);  // list of shape parameters for tt bkg                                
    construct_formula(procName_D3,*bkg_tt_bins_D3,parlist_D3,nuisanceParams_D3);
  } else 
  {
    RooArgList parlist_D3(N7_tt_D3,a0_tt_D3,a1_tt_D3,d_tt_D3);  // list of shape parameters for tt bkg                       
    //RooArgList parlist_D3(N7_tt_D3,a0_tt_D3,a1_tt_D3,a2_tt_D3);  // list of shape parameters for tt bkg
    construct_formula(procName_D3,*bkg_tt_bins_D3,parlist_D3,nuisanceParams_D3);
  }
  RooParametricHist background_tt_D3(procName_D3.c_str(),"",*wspace->var("CMS_th1x"),*bkg_tt_bins_D3,*data_th1_D3);
  wspace->import(background_tt_D3,RooFit::RecycleConflictNodes());
  std::stringstream procNameD3Norm;
  procNameD3Norm << procName_D3 << "_norm";
  RooAddition tt_norm_D3(procNameD3Norm.str().c_str(),"",*bkg_tt_bins_D3);
  wspace->import(tt_norm_D3,RooFit::RecycleConflictNodes());
  
  // ---------------------- MVA bin 4  ------------------
  
  // Dataset with 8 bins
  RooDataHist data_hist_D4("data_obs_D4","Data observed in MVA bin 4",vars_D4,data_th1_D4);
  wspace->import(data_hist_D4);
  
  //list of nuisance parameters for tt bkg D4
  const std::vector<NuisanceParam>& nuisanceParams_D4 = {
      {*wspace->var(("np_tt_JECUp_"+year).c_str()), tt_syst_JECUp_D4},
      {*wspace->var(("np_tt_JERUp_"+year).c_str()), tt_syst_JERUp_D4},
      {*wspace->var(("np_tt_btg_"+year).c_str()),   tt_syst_btg_D4},
      {*wspace->var(("np_tt_lep_"+year).c_str()),   tt_syst_lep_D4},
      {*wspace->var(("np_tt_nom_"+year).c_str()),   tt_syst_nom_D4},
      {
          *wspace->var(("np_tt_qcdCR_"+year).c_str()),   tt_syst_qcdCR_D4, tt_qcdCRErr_D4, 
          { 
              *wspace->var(("np_tt_qcdCRErrD4Bin1_"+year).c_str()),
              *wspace->var(("np_tt_qcdCRErrD4Bin2_"+year).c_str()),
              *wspace->var(("np_tt_qcdCRErrD4Bin3_"+year).c_str()),
              *wspace->var(("np_tt_qcdCRErrD4Bin4_"+year).c_str()),
              *wspace->var(("np_tt_qcdCRErrD4Bin5_"+year).c_str()),
              *wspace->var(("np_tt_qcdCRErrD4Bin6_"+year).c_str()),
          }
      },
      {*wspace->var("np_tt_pdf"),                     tt_syst_pdf_D4},
      {*wspace->var("np_tt_FSR"),                     tt_syst_FSR_D4},
      {*wspace->var("np_tt_ISR"),                     tt_syst_ISR_D4},
      {*wspace->var("np_tt_scl"),                     tt_syst_scl_D4},
      {*wspace->var(("np_tt_ht_"+year).c_str()),      tt_syst_ht_D4},
      {*wspace->var(("np_tt_httail_"+year).c_str()),  tt_syst_httail_D4},
      {*wspace->var(("np_tt_htnjet_"+year).c_str()),  tt_syst_htnjet_D4},
      {*wspace->var(("np_tt_pu_"+year).c_str()),      tt_syst_pu_D4},
      {*wspace->var(("np_tt_JECDown_"+year).c_str()), tt_syst_JECDown_D4},
      {*wspace->var(("np_tt_JERDown_"+year).c_str()), tt_syst_JERDown_D4}
  };
  
  RooArgList *bkg_tt_bins_D4 = new RooArgList();
  std::string procName_D4 = "background_tt_D4_"+year;
  if (shared) 
  {
    RooArgList parlist_D4(N7_tt_D4,a0_tt,a1_tt,d_tt);  // list of shape parameters for tt bkg
    //RooArgList parlist_D4(N7_tt_D4,a0_tt,a1_tt,a2_tt);  // list of shape parameters for tt bkg                                
    construct_formula(procName_D4,*bkg_tt_bins_D4,parlist_D4,nuisanceParams_D4);
  } else 
  {
    RooArgList parlist_D4(N7_tt_D4,a0_tt_D4,a1_tt_D4,d_tt_D4);  // list of shape parameters for tt bkg                       
    //RooArgList parlist_D4(N7_tt_D4,a0_tt_D4,a1_tt_D4,a2_tt_D4);  // list of shape parameters for tt bkg
    construct_formula(procName_D4,*bkg_tt_bins_D4,parlist_D4,nuisanceParams_D4);
  }
  RooParametricHist background_tt_D4(procName_D4.c_str(),"",*wspace->var("CMS_th1x"),*bkg_tt_bins_D4,*data_th1_D4);
  wspace->import(background_tt_D4,RooFit::RecycleConflictNodes());
  std::stringstream procNameD4Norm;
  procNameD4Norm << procName_D4 << "_norm";
  RooAddition tt_norm_D4(procNameD4Norm.str().c_str(),"",*bkg_tt_bins_D4);
  wspace->import(tt_norm_D4,RooFit::RecycleConflictNodes());
  
  // =================================================================================
     
  fOut->cd();
  
  // Shape histograms for signal
  sigMC_th1_D1->SetName("sigMC_th1_D1");
  sigMC_th1_D2->SetName("sigMC_th1_D2");
  sigMC_th1_D3->SetName("sigMC_th1_D3");
  sigMC_th1_D4->SetName("sigMC_th1_D4");
  sigMC_th1_D1->Write();
  sigMC_th1_D2->Write();
  sigMC_th1_D3->Write();
  sigMC_th1_D4->Write();
  
  // Shape histograms for other backgrounds
  otherMC_th1_D1->SetName("otherMC_th1_D1");
  otherMC_th1_D2->SetName("otherMC_th1_D2");
  otherMC_th1_D3->SetName("otherMC_th1_D3");
  otherMC_th1_D4->SetName("otherMC_th1_D4");
  otherMC_th1_D1->Write();
  otherMC_th1_D2->Write();
  otherMC_th1_D3->Write();
  otherMC_th1_D4->Write();

  // Shape histograms for QCD backgrounds
  qcdMC_th1_D1->SetName("qcdMC_th1_D1");
  qcdMC_th1_D2->SetName("qcdMC_th1_D2");
  qcdMC_th1_D3->SetName("qcdMC_th1_D3");
  qcdMC_th1_D4->SetName("qcdMC_th1_D4");
  qcdMC_th1_D1->Write();
  qcdMC_th1_D2->Write();
  qcdMC_th1_D3->Write();
  qcdMC_th1_D4->Write();

  // Shape histograms for TTX backgrounds
  ttxMC_th1_D1->SetName("ttxMC_th1_D1");
  ttxMC_th1_D2->SetName("ttxMC_th1_D2");
  ttxMC_th1_D3->SetName("ttxMC_th1_D3");
  ttxMC_th1_D4->SetName("ttxMC_th1_D4");
  ttxMC_th1_D1->Write();
  ttxMC_th1_D2->Write();
  ttxMC_th1_D3->Write();
  ttxMC_th1_D4->Write();
  
  // =================================================================================
  // Systematics
  
  // Signal systematics
  WriteHisto2WS<TH1D>(file, model+"_"+mass+"_h_njets_pt30_1l_JECUp",   "SIG_JEC_"+year+"Up",   {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, model+"_"+mass+"_h_njets_pt30_1l_JECDown", "SIG_JEC_"+year+"Down", {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, model+"_"+mass+"_h_njets_pt30_1l_JERUp",   "SIG_JER_"+year+"Up",   {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, model+"_"+mass+"_h_njets_pt30_1l_JERDown", "SIG_JER_"+year+"Down", {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, model+"_"+mass+"_h_njets_pt30_1l_btgUp",   "SIG_btg_"+year+"Up",   {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, model+"_"+mass+"_h_njets_pt30_1l_btgDown", "SIG_btg_"+year+"Down", {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, model+"_"+mass+"_h_njets_pt30_1l_lepUp",   "SIG_lep_"+year+"Up",   {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, model+"_"+mass+"_h_njets_pt30_1l_lepDown", "SIG_lep_"+year+"Down", {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, model+"_"+mass+"_h_njets_pt30_1l_pdfUp",   "SIG_pdfUp",            {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, model+"_"+mass+"_h_njets_pt30_1l_pdfDown", "SIG_pdfDown",          {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, model+"_"+mass+"_h_njets_pt30_1l_sclUp",   "SIG_sclUp",            {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, model+"_"+mass+"_h_njets_pt30_1l_sclDown", "SIG_sclDown",          {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, model+"_"+mass+"_h_njets_pt30_1l_puUp",    "SIG_pu_"+year+"Up",    {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, model+"_"+mass+"_h_njets_pt30_1l_puDown",  "SIG_pu_"+year+"Down",  {"D1","D2","D3","D4"});
  
  // "OTHER" background systematics
  WriteHisto2WS<TH1D>(file, "OTHER_h_njets_pt30_1l_JECUp",   "OTHER_JEC_"+year+"Up",   {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "OTHER_h_njets_pt30_1l_JECDown", "OTHER_JEC_"+year+"Down", {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "OTHER_h_njets_pt30_1l_JERUp",   "OTHER_JER_"+year+"Up",   {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "OTHER_h_njets_pt30_1l_JERDown", "OTHER_JER_"+year+"Down", {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "OTHER_h_njets_pt30_1l_btgUp",   "OTHER_btg_"+year+"Up",   {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "OTHER_h_njets_pt30_1l_btgDown", "OTHER_btg_"+year+"Down", {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "OTHER_h_njets_pt30_1l_lepUp",   "OTHER_lep_"+year+"Up",   {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "OTHER_h_njets_pt30_1l_lepDown", "OTHER_lep_"+year+"Down", {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "OTHER_h_njets_pt30_1l_pdfUp",   "OTHER_pdfUp",            {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "OTHER_h_njets_pt30_1l_pdfDown", "OTHER_pdfDown",          {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "OTHER_h_njets_pt30_1l_sclUp",   "OTHER_sclUp",            {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "OTHER_h_njets_pt30_1l_sclDown", "OTHER_sclDown",          {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "OTHER_h_njets_pt30_1l_htUp",    "OTHER_ht_"+year+"Up",    {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "OTHER_h_njets_pt30_1l_htDown",  "OTHER_ht_"+year+"Down",  {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "OTHER_h_njets_pt30_1l_puUp",    "OTHER_pu_"+year+"Up",    {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "OTHER_h_njets_pt30_1l_puDown",  "OTHER_pu_"+year+"Down",  {"D1","D2","D3","D4"});
    
  // "TTX" background systematics
  WriteHisto2WS<TH1D>(file, "TTX_h_njets_pt30_1l_JECUp",   "TTX_JEC_"+year+"Up",   {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "TTX_h_njets_pt30_1l_JECDown", "TTX_JEC_"+year+"Down", {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "TTX_h_njets_pt30_1l_JERUp",   "TTX_JER_"+year+"Up",   {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "TTX_h_njets_pt30_1l_JERDown", "TTX_JER_"+year+"Down", {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "TTX_h_njets_pt30_1l_btgUp",   "TTX_btg_"+year+"Up",   {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "TTX_h_njets_pt30_1l_btgDown", "TTX_btg_"+year+"Down", {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "TTX_h_njets_pt30_1l_lepUp",   "TTX_lep_"+year+"Up",   {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "TTX_h_njets_pt30_1l_lepDown", "TTX_lep_"+year+"Down", {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "TTX_h_njets_pt30_1l_pdfUp",   "TTX_pdfUp",            {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "TTX_h_njets_pt30_1l_pdfDown", "TTX_pdfDown",          {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "TTX_h_njets_pt30_1l_sclUp",   "TTX_sclUp",            {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "TTX_h_njets_pt30_1l_sclDown", "TTX_sclDown",          {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "TTX_h_njets_pt30_1l_htUp",    "TTX_ht_"+year+"Up",    {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "TTX_h_njets_pt30_1l_htDown",  "TTX_ht_"+year+"Down",  {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "TTX_h_njets_pt30_1l_puUp",    "TTX_pu_"+year+"Up",    {"D1","D2","D3","D4"});
  WriteHisto2WS<TH1D>(file, "TTX_h_njets_pt30_1l_puDown",  "TTX_pu_"+year+"Down",  {"D1","D2","D3","D4"});    
    
  // =================================================================================
  // Statistics-based Uncertainties
  
  // MC stat uncertainty histograms
  const std::vector<std::string>& nJetBin = {"1","2","3","4","5","6"};
  for(const auto& i : nJetBin)
  {
      // Particular signal model and mass point
      WriteHisto2WS<TH1D>(file, model+"_"+mass+"_mcStatBin"+i+"Up",   i+"_"+year+"Up",   {"D1","D2","D3","D4"}, "SIG");      
      WriteHisto2WS<TH1D>(file, model+"_"+mass+"_mcStatBin"+i+"Down", i+"_"+year+"Down", {"D1","D2","D3","D4"}, "SIG");      

      // OTHER backgrounds
      WriteHisto2WS<TH1D>(file, "OTHER_mcStatBin"+i+"Up",   i+"_"+year+"Up",   {"D1","D2","D3","D4"}, "OTHER");      
      WriteHisto2WS<TH1D>(file, "OTHER_mcStatBin"+i+"Down", i+"_"+year+"Down", {"D1","D2","D3","D4"}, "OTHER");      

      // QCD backgrounds
      WriteHisto2WS<TH1D>(file, "QCD_mcStatBin"+i+"Up",   i+"_"+year+"Up",   {"D1","D2","D3","D4"}, "QCD");      
      WriteHisto2WS<TH1D>(file, "QCD_mcStatBin"+i+"Down", i+"_"+year+"Down", {"D1","D2","D3","D4"}, "QCD");      

      // TTX backgrounds
      WriteHisto2WS<TH1D>(file, "TTX_mcStatBin"+i+"Up",   i+"_"+year+"Up",   {"D1","D2","D3","D4"}, "TTX");      
      WriteHisto2WS<TH1D>(file, "TTX_mcStatBin"+i+"Down", i+"_"+year+"Down", {"D1","D2","D3","D4"}, "TTX");      
  }

  //Finally write out the work space
  wspace->Write();
  
  // TCanvas *c1 = new TCanvas("c1","c1");
  // data_hist_D4.createHistogram("nj")->Draw("H");
  // TCanvas *c2 = new TCanvas("c2","c2");
  // ttMC_hist_D4.createHistogram("nj")->Draw("H");
  // TCanvas *c3 = new TCanvas("c3","c3");
  // otherMC_hist_D4.createHistogram("nj")->Draw("H");
  // TCanvas *c4 = new TCanvas("c4","c4");
  // sigMC_hist_D4.createHistogram("nj")->Draw("H");
  
  // TCanvas *c1 = new TCanvas("c1","c1");
  // data_hist_D4.createHistogram("nj")->Draw("H");
  // TCanvas *c2 = new TCanvas("c2","c2");
  // ttMC_hist_D4.createHistogram("nj")->Draw("H");
  // TCanvas *c3 = new TCanvas("c3","c3");
  // otherMC_hist_D4.createHistogram("nj")->Draw("H")
  // TCanvas *c4 = new TCanvas("c4","c4");
  // sigMC_hist_D4.createHistogram("nj")->Draw("H");
}

int main()
{
    make_MVA_8bin_ws();
}
