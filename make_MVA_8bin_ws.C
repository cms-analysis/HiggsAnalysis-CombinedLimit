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


Double_t step(double_t x) {
  return 1;
}

void construct_formula(string procName, RooArgList& binlist, const RooArgList& paramlist, const RooArgList& NPs, TH1D* h_syst[]) {

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

  ROOT::v5::TFormula::SetMaxima(10000);

  int max_bin = 18; // 14 means just njets=14, 20 means last bin is inclusive up through njets=20

  for (int i=0; i<8; i++) {

    stringstream form;
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
    if (i==7) {
      for (int k=8; k<=max_bin-7; k++) {
	form << " + 1";
	for (int j=0; j<k; j++) {
	  form << "*( @3>1 ? ( @2 - 1./@3 + TMath::Power(  TMath::Power( 1./@3 , "<<j<<" ) / TMath::Power( @1-@2+1./@3 , "<<j-2<<" ) , 1/2 ) ) : ( @2 - (2.-@3) + TMath::Power( TMath::Power( (2.-@3) , "<<j<<" ) / TMath::Power( @1-@2+(2.-@3) , "<<j-2<<" ) , 1/2 ) ) )";
	}
      }
    }
    form << ")";

    if (i==0) {
      form << "*TMath::Power(" << h_syst[0]->GetBinContent(i+1) << ",@1)";
      form << "*TMath::Power(" << h_syst[1]->GetBinContent(i+1) << ",@2)";
      form << "*TMath::Power(" << h_syst[2]->GetBinContent(i+1) << ",@3)";
      form << "*TMath::Power(" << h_syst[3]->GetBinContent(i+1) << ",@4)";
      form << "*TMath::Power(" << h_syst[4]->GetBinContent(i+1) << ",@5)";
      form << "*TMath::Power(" << h_syst[5]->GetBinContent(i+1) << ",@6)";
      form << "*TMath::Power(" << h_syst[6]->GetBinContent(i+1) << ",@7)";
      form << "*TMath::Power(" << h_syst[7]->GetBinContent(i+1) << ",@8)";
      form << "*TMath::Power(" << h_syst[8]->GetBinContent(i+1) << ",@9)";
    } else if (i>=1) {
      form << "*TMath::Power(" << h_syst[0]->GetBinContent(i+1) << ",@4)";
      form << "*TMath::Power(" << h_syst[1]->GetBinContent(i+1) << ",@5)";
      form << "*TMath::Power(" << h_syst[2]->GetBinContent(i+1) << ",@6)";
      form << "*TMath::Power(" << h_syst[3]->GetBinContent(i+1) << ",@7)";
      form << "*TMath::Power(" << h_syst[4]->GetBinContent(i+1) << ",@8)";
      form << "*TMath::Power(" << h_syst[5]->GetBinContent(i+1) << ",@9)";
      form << "*TMath::Power(" << h_syst[6]->GetBinContent(i+1) << ",@10)";
      form << "*TMath::Power(" << h_syst[7]->GetBinContent(i+1) << ",@11)";
      form << "*TMath::Power(" << h_syst[8]->GetBinContent(i+1) << ",@12)";
    }
    // nuisance parameters
    formArgList.add(NPs[0]);
    formArgList.add(NPs[1]);
    formArgList.add(NPs[2]);
    formArgList.add(NPs[3]);
    formArgList.add(NPs[4]);
    formArgList.add(NPs[5]);
    formArgList.add(NPs[6]);
    formArgList.add(NPs[7]);
    formArgList.add(NPs[8]);

    // Create RooFormulaVar for this bin
    stringstream binName;
    binName << procName << "_b" << i;
    RooFormulaVar* binvar = new RooFormulaVar(binName.str().c_str(), "", form.str().c_str(), RooArgList(formArgList));
    binlist.add(*binvar);

    cout << "bin i = " << i << " , njets = " << i+7 << endl;
    cout << "process bin name : " << binName.str().c_str() << endl;
    cout << "Formula : " << form.str().c_str() << endl;
    formArgList.Print();
    cout << endl;

  }

}


void make_MVA_8bin_ws(const string year = "2016", const string infile_path = "Keras_V1.2.5_v2", const string model = "RPV", const string mass = "550", const string dataType = "pseudodata", const string syst = "", bool shared = true, bool TTonly = false) {
  using namespace RooFit;
  // Load the combine library to get access to the RooParametricHist
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
    
  // Output file and workspace 
  TFile *fOut = new TFile(("MVA_"+year+"_"+model+"_"+mass+"_ws.root").c_str(),"RECREATE");
  RooWorkspace *wspace = new RooWorkspace("wspace","wspace");

  // // njet is our variable, 8 bins, from 7 up through 14,
  // //   Note that njet=14 is inclusive as >=14
  // // D1, D2, D3, D4 are the MVA bins
  // wspace->factory("nj_D1[6.5,14.5]");
  // //wspace->factory("nj_D1[0,8]");
  // wspace->var("nj_D1")->setBins(8);
  // RooArgSet vars_D1(*wspace->var("nj_D1"));

  // wspace->factory("nj_D2[6.5,14.5]");
  // //wspace->factory("nj_D2[0,8]");
  // wspace->var("nj_D2")->setBins(8);
  // RooArgSet vars_D2(*wspace->var("nj_D2"));

  // wspace->factory("nj_D3[6.5,14.5]");
  // //wspace->factory("nj_D3[0,8]");
  // wspace->var("nj_D3")->setBins(8);
  // RooArgSet vars_D3(*wspace->var("nj_D3"));

  // wspace->factory("nj_D4[6.5,14.5]");
  // //wspace->factory("nj_D4[0,8]");
  // wspace->var("nj_D4")->setBins(8);
  // RooArgSet vars_D4(*wspace->var("nj_D4"));

  wspace->factory("CMS_th1x[0,8]");
  wspace->var("CMS_th1x")->setBins(8);
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
  //TH1D* data_th1_D1 = (TH1D*)file->Get(("D1_pseudodataS_"+model+"_"+mass+"_h_njets_pt30_1l_JECUp").c_str()); // JEC UP with signal
  //TH1D* data_th1_D1 = (TH1D*)file->Get("D1_pseudodata_h_njets_pt30_1l_JECUp"); // JEC UP without signal
  TH1D* otherMC_th1_D1 = (TH1D*)file->Get(("D1_OTHER_h_njets_pt30_1l"+syst).c_str());
  TH1D* sigMC_th1_D1 = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_h_njets_pt30_1l"+syst).c_str());

  TH1D* data_th1_D2;
  if     (dataType == "data")        data_th1_D2 = (TH1D*)file->Get(("D2_data_h_njets_pt30_1l"+syst).c_str());  // Actual data -- be careful
  else if(dataType == "pseudodata" ) data_th1_D2 = (TH1D*)file->Get(("D2_pseudodata_h_njets_pt30_1l"+syst).c_str()); // without signal
  else if(dataType == "pseudodataS") data_th1_D2 = (TH1D*)file->Get(("D2_pseudodataS_"+model+"_"+mass+"_h_njets_pt30_1l"+syst).c_str());  // with signal
  else if(dataType == "pseudodataA" ) data_th1_D2 = (TH1D*)file->Get(("D2_pseudodataFunc28_24_236_h_njets_pt30_1l"+syst).c_str()); //
  else if(dataType == "pseudodataB" ) data_th1_D2 = (TH1D*)file->Get(("D2_pseudodataFunc28_24_18_h_njets_pt30_1l"+syst).c_str()); //
  else if(dataType == "pseudodataC" ) data_th1_D2 = (TH1D*)file->Get(("D2_pseudodataFunc28_24_-20_h_njets_pt30_1l"+syst).c_str()); //
  else                               data_th1_D2 = (TH1D*)file->Get(("D2_"+dataType+"_h_njets_pt30_1l"+syst).c_str());
  //TH1D* data_th1_D2 = (TH1D*)file->Get(("D2_pseudodataS_"+model+"_"+mass+"_h_njets_pt30_1l_JECUp").c_str()); // JEC UP with signal
  //TH1D* data_th1_D2 = (TH1D*)file->Get("D2_pseudodata_h_njets_pt30_1l_JECUp"); // JEC UP without signal
  TH1D* otherMC_th1_D2 = (TH1D*)file->Get(("D2_OTHER_h_njets_pt30_1l"+syst).c_str());
  TH1D* sigMC_th1_D2 = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_h_njets_pt30_1l"+syst).c_str());

  TH1D* data_th1_D3;
  if     (dataType == "data")        data_th1_D3 = (TH1D*)file->Get(("D3_data_h_njets_pt30_1l"+syst).c_str());  // Actual data -- be careful
  else if(dataType == "pseudodata" ) data_th1_D3 = (TH1D*)file->Get(("D3_pseudodata_h_njets_pt30_1l"+syst).c_str()); // without signal
  else if(dataType == "pseudodataS") data_th1_D3 = (TH1D*)file->Get(("D3_pseudodataS_"+model+"_"+mass+"_h_njets_pt30_1l"+syst).c_str());  // with signal
  else if(dataType == "pseudodataA" ) data_th1_D3 = (TH1D*)file->Get(("D3_pseudodataFunc28_24_236_h_njets_pt30_1l"+syst).c_str()); //
  else if(dataType == "pseudodataB" ) data_th1_D3 = (TH1D*)file->Get(("D3_pseudodataFunc28_24_18_h_njets_pt30_1l"+syst).c_str()); //
  else if(dataType == "pseudodataC" ) data_th1_D3 = (TH1D*)file->Get(("D3_pseudodataFunc28_24_-20_h_njets_pt30_1l"+syst).c_str()); //
  else                               data_th1_D3 = (TH1D*)file->Get(("D3_"+dataType+"_h_njets_pt30_1l"+syst).c_str());
  //TH1D* data_th1_D3 = (TH1D*)file->Get(("D3_pseudodataS_"+model+"_"+mass+"_h_njets_pt30_1l_JECUp").c_str()); // JEC UP with signal
  //TH1D* data_th1_D3 = (TH1D*)file->Get("D3_pseudodata_h_njets_pt30_1l_JECUp"); // JEC UP without signal
  TH1D* otherMC_th1_D3 = (TH1D*)file->Get(("D3_OTHER_h_njets_pt30_1l"+syst).c_str());
  TH1D* sigMC_th1_D3 = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_h_njets_pt30_1l"+syst).c_str());

  TH1D* data_th1_D4;
  if     (dataType == "data")        data_th1_D4 = (TH1D*)file->Get(("D4_data_h_njets_pt30_1l"+syst).c_str());  // Actual data -- be careful
  else if(dataType == "pseudodata" ) data_th1_D4 = (TH1D*)file->Get(("D4_pseudodata_h_njets_pt30_1l"+syst).c_str()); // without signal
  else if(dataType == "pseudodataS") data_th1_D4 = (TH1D*)file->Get(("D4_pseudodataS_"+model+"_"+mass+"_h_njets_pt30_1l"+syst).c_str());  // with signal
  else if(dataType == "pseudodataA" ) data_th1_D4 = (TH1D*)file->Get(("D4_pseudodataFunc28_24_236_h_njets_pt30_1l"+syst).c_str()); //
  else if(dataType == "pseudodataB" ) data_th1_D4 = (TH1D*)file->Get(("D4_pseudodataFunc28_24_18_h_njets_pt30_1l"+syst).c_str()); //
  else if(dataType == "pseudodataC" ) data_th1_D4 = (TH1D*)file->Get(("D4_pseudodataFunc28_24_-20_h_njets_pt30_1l"+syst).c_str()); //
  else                               data_th1_D4 = (TH1D*)file->Get(("D4_"+dataType+"_h_njets_pt30_1l"+syst).c_str());
  //TH1D* data_th1_D4 = (TH1D*)file->Get(("D4_pseudodataS_"+model+"_"+mass+"_h_njets_pt30_1l_JECUp").c_str()); // JEC UP with signal
  //TH1D* data_th1_D4 = (TH1D*)file->Get("D4_pseudodata_h_njets_pt30_1l_JECUp"); // JEC UP without signal
  TH1D* otherMC_th1_D4 = (TH1D*)file->Get(("D4_OTHER_h_njets_pt30_1l"+syst).c_str());
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

  // //double_t n7_tt_portion_D1 = data_th1_D1->GetBinContent(1) - otherMC_th1_D1->GetBinContent(1) - sigMC_th1_D1->GetBinContent(1);
  // double_t n7_tt_portion_D1 = data_th1_D1->GetBinContent(1) - otherMC_th1_D1->GetBinContent(1);
  // double_t n7_tt_portion_D1_low = n7_tt_portion_D1-2000;
  // if (n7_tt_portion_D1_low<0) n7_tt_portion_D1_low = 0;
  // RooRealVar N7_tt_D1("N7_tt_D1","njets 7 for tt bkg in MVA D1",n7_tt_portion_D1,n7_tt_portion_D1_low,n7_tt_portion_D1+1500);

  // //double_t n7_tt_portion_D2 = data_th1_D2->GetBinContent(1) - otherMC_th1_D2->GetBinContent(1) - sigMC_th1_D2->GetBinContent(1);
  // double_t n7_tt_portion_D2 = data_th1_D2->GetBinContent(1) - otherMC_th1_D2->GetBinContent(1);
  // double_t n7_tt_portion_D2_low = n7_tt_portion_D2-1500;
  // if (n7_tt_portion_D2_low<0) n7_tt_portion_D2_low = 0;
  // RooRealVar N7_tt_D2("N7_tt_D2","njets 7 for tt bkg in MVA D2",n7_tt_portion_D2,n7_tt_portion_D2_low,n7_tt_portion_D2+1500);

  // //double_t n7_tt_portion_D3 = data_th1_D3->GetBinContent(1) - otherMC_th1_D3->GetBinContent(1) - sigMC_th1_D3->GetBinContent(1);
  // double_t n7_tt_portion_D3 = data_th1_D3->GetBinContent(1) - otherMC_th1_D3->GetBinContent(1);
  // double_t n7_tt_portion_D3_low = n7_tt_portion_D3-1500;
  // if (n7_tt_portion_D3_low<0) n7_tt_portion_D3_low = 0;
  // RooRealVar N7_tt_D3("N7_tt_D3","njets 7 for tt bkg in MVA D3",n7_tt_portion_D3,n7_tt_portion_D3_low,n7_tt_portion_D3+1500);

  // //double_t n7_tt_portion_D4 = data_th1_D4->GetBinContent(1) - otherMC_th1_D4->GetBinContent(1) - sigMC_th1_D4->GetBinContent(1);
  // double_t n7_tt_portion_D4 = data_th1_D4->GetBinContent(1) - otherMC_th1_D4->GetBinContent(1);
  // double_t n7_tt_portion_D4_low = n7_tt_portion_D4-1500;
  // if (n7_tt_portion_D4_low<0) n7_tt_portion_D4_low = 0;
  // RooRealVar N7_tt_D4("N7_tt_D4","njets 7 for tt bkg in MVA D4",n7_tt_portion_D4,n7_tt_portion_D4_low,n7_tt_portion_D4+1500);



  //double_t n7_tt_portion_D1 = data_th1_D1->GetBinContent(1) - otherMC_th1_D1->GetBinContent(1) - sigMC_th1_D1->GetBinContent(1);
  double_t n7_tt_portion_D1 = TTonly ? data_th1_D1->GetBinContent(1) : data_th1_D1->GetBinContent(1) - otherMC_th1_D1->GetBinContent(1);
  double_t n7_tt_portion_D1_low = n7_tt_portion_D1-6000;
  if (n7_tt_portion_D1_low<0) n7_tt_portion_D1_low = 0;
  RooRealVar N7_tt_D1(("N7_tt_D1_"+year).c_str(),"njets 7 for tt bkg in MVA D1",n7_tt_portion_D1,n7_tt_portion_D1_low,n7_tt_portion_D1+6000);

  //double_t n7_tt_portion_D2 = data_th1_D2->GetBinContent(1) - otherMC_th1_D2->GetBinContent(1) - sigMC_th1_D2->GetBinContent(1);
  double_t n7_tt_portion_D2 = TTonly ? data_th1_D2->GetBinContent(1) : data_th1_D2->GetBinContent(1) - otherMC_th1_D2->GetBinContent(1);
  double_t n7_tt_portion_D2_low = n7_tt_portion_D2-5000;
  if (n7_tt_portion_D2_low<0) n7_tt_portion_D2_low = 0;
  RooRealVar N7_tt_D2(("N7_tt_D2_"+year).c_str(),"njets 7 for tt bkg in MVA D2",n7_tt_portion_D2,n7_tt_portion_D2_low,n7_tt_portion_D2+5000);

  //double_t n7_tt_portion_D3 = data_th1_D3->GetBinContent(1) - otherMC_th1_D3->GetBinContent(1) - sigMC_th1_D3->GetBinContent(1);
  double_t n7_tt_portion_D3 = TTonly ? data_th1_D3->GetBinContent(1) : data_th1_D3->GetBinContent(1) - otherMC_th1_D3->GetBinContent(1);
  double_t n7_tt_portion_D3_low = n7_tt_portion_D3-4000;
  if (n7_tt_portion_D3_low<0) n7_tt_portion_D3_low = 0;
  RooRealVar N7_tt_D3(("N7_tt_D3_"+year).c_str(),"njets 7 for tt bkg in MVA D3",n7_tt_portion_D3,n7_tt_portion_D3_low,n7_tt_portion_D3+4000);

  //double_t n7_tt_portion_D4 = data_th1_D4->GetBinContent(1) - otherMC_th1_D4->GetBinContent(1) - sigMC_th1_D4->GetBinContent(1);
  double_t n7_tt_portion_D4 = TTonly ? data_th1_D4->GetBinContent(1) : data_th1_D4->GetBinContent(1) - otherMC_th1_D4->GetBinContent(1);
  double_t n7_tt_portion_D4_low = n7_tt_portion_D4-3000;
  if (n7_tt_portion_D4_low<0) n7_tt_portion_D4_low = 0;
  RooRealVar N7_tt_D4(("N7_tt_D4_"+year).c_str(),"njets 7 for tt bkg in MVA D4",n7_tt_portion_D4,n7_tt_portion_D4_low,n7_tt_portion_D4+3000);



  // tt shape systematic nuisance parameters
  wspace->factory("np_tt_JEC[0.0]"); // fully correlated
  //wspace->factory("np_tt_JEC[0.0,0.0,0.0]");
  wspace->factory(("np_tt_JER_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_btg_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_lep_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory(("np_tt_nom_"+year+"[0.0]").c_str()); // uncorrelated
  wspace->factory("np_tt_pdf[0.0]"); // fully correlated
  wspace->factory("np_tt_FSR[0.0]"); // fully correlated
  if (year=="2016") wspace->factory("np_tt_ht_2016[0.0]");
  wspace->factory("np_tt_isr[0.0]"); // fully correlated
  if (year=="2017") wspace->factory("np_tt_scl_2017[0.0]");


  // Load in the histograms with the bin-by-bin ratios to be used in the ttbar shape systematics
  TFile* tt_syst_file = TFile::Open((infile_path+"/ttbar_systematics.root").c_str());

  TH1D* tt_syst_JEC_D1 = (TH1D*)tt_syst_file->Get("D1_JEC");
  TH1D* tt_syst_JEC_D2 = (TH1D*)tt_syst_file->Get("D2_JEC");
  TH1D* tt_syst_JEC_D3 = (TH1D*)tt_syst_file->Get("D3_JEC");
  TH1D* tt_syst_JEC_D4 = (TH1D*)tt_syst_file->Get("D4_JEC");

  TH1D* tt_syst_JER_D1 = (TH1D*)tt_syst_file->Get("D1_JER");
  TH1D* tt_syst_JER_D2 = (TH1D*)tt_syst_file->Get("D2_JER");
  TH1D* tt_syst_JER_D3 = (TH1D*)tt_syst_file->Get("D3_JER");
  TH1D* tt_syst_JER_D4 = (TH1D*)tt_syst_file->Get("D4_JER");

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

  TH1D* tt_syst_pdf_D1 = (TH1D*)tt_syst_file->Get("D1_pdf");
  TH1D* tt_syst_pdf_D2 = (TH1D*)tt_syst_file->Get("D2_pdf");
  TH1D* tt_syst_pdf_D3 = (TH1D*)tt_syst_file->Get("D3_pdf");
  TH1D* tt_syst_pdf_D4 = (TH1D*)tt_syst_file->Get("D4_pdf");

  TH1D* tt_syst_FSR_D1 = (TH1D*)tt_syst_file->Get("D1_FSR");
  TH1D* tt_syst_FSR_D2 = (TH1D*)tt_syst_file->Get("D2_FSR");
  TH1D* tt_syst_FSR_D3 = (TH1D*)tt_syst_file->Get("D3_FSR");
  TH1D* tt_syst_FSR_D4 = (TH1D*)tt_syst_file->Get("D4_FSR");

  TH1D* tt_syst_isr_D1 = (TH1D*)tt_syst_file->Get("D1_isr");
  TH1D* tt_syst_isr_D2 = (TH1D*)tt_syst_file->Get("D2_isr");
  TH1D* tt_syst_isr_D3 = (TH1D*)tt_syst_file->Get("D3_isr");
  TH1D* tt_syst_isr_D4 = (TH1D*)tt_syst_file->Get("D4_isr");

  TH1D* tt_syst_ht_D1;
  TH1D* tt_syst_ht_D2;
  TH1D* tt_syst_ht_D3;
  TH1D* tt_syst_ht_D4;
  if (year=="2016") {
    tt_syst_ht_D1 = (TH1D*)tt_syst_file->Get("D1_ht");
    tt_syst_ht_D2 = (TH1D*)tt_syst_file->Get("D2_ht");
    tt_syst_ht_D3 = (TH1D*)tt_syst_file->Get("D3_ht");
    tt_syst_ht_D4 = (TH1D*)tt_syst_file->Get("D4_ht");
  }

  TH1D* tt_syst_scl_D1;
  TH1D* tt_syst_scl_D2;
  TH1D* tt_syst_scl_D3;
  TH1D* tt_syst_scl_D4;
  if (year=="2017") {
    tt_syst_scl_D1 = (TH1D*)tt_syst_file->Get("D1_scl");
    tt_syst_scl_D2 = (TH1D*)tt_syst_file->Get("D2_scl");
    tt_syst_scl_D3 = (TH1D*)tt_syst_file->Get("D3_scl");
    tt_syst_scl_D4 = (TH1D*)tt_syst_file->Get("D4_scl");
  }


  // ----------------------  MVA bin 1  ------------------

  // Dataset with 8 bins
  RooDataHist data_hist_D1("data_obs_D1","Data observed in MVA bin 1",vars_D1,data_th1_D1);
  wspace->import(data_hist_D1);

  // ttbar bkg in D1
  RooArgList bkg_tt_syst_NP_D1(*wspace->var("np_tt_JEC"),
			       *wspace->var(("np_tt_JER_"+year).c_str()),
			       *wspace->var(("np_tt_btg_"+year).c_str()),
			       *wspace->var(("np_tt_lep_"+year).c_str()),
			       *wspace->var(("np_tt_nom_"+year).c_str()),
			       *wspace->var("np_tt_pdf"),
			       *wspace->var("np_tt_FSR"),
			       *wspace->var("np_tt_isr"));  // list of nuisance parameters for tt bkg
  if (year=="2016") {
    bkg_tt_syst_NP_D1.add(*wspace->var("np_tt_ht_2016"));
  }
  if (year=="2017") {
    bkg_tt_syst_NP_D1.add(*wspace->var("np_tt_scl_2017"));
  }
  TH1D* bkg_tt_syst_histos_D1[20];  // array of histograms containing each tt bkg shape uncertainty
  bkg_tt_syst_histos_D1[0] = tt_syst_JEC_D1;
  bkg_tt_syst_histos_D1[1] = tt_syst_JER_D1;
  bkg_tt_syst_histos_D1[2] = tt_syst_btg_D1;
  bkg_tt_syst_histos_D1[3] = tt_syst_lep_D1;
  bkg_tt_syst_histos_D1[4] = tt_syst_nom_D1;
  bkg_tt_syst_histos_D1[5] = tt_syst_pdf_D1;
  bkg_tt_syst_histos_D1[6] = tt_syst_FSR_D1;
  bkg_tt_syst_histos_D1[7] = tt_syst_isr_D1;
  if (year=="2016") {
    bkg_tt_syst_histos_D1[8] = tt_syst_ht_D1;
  }
  if (year=="2017") {
    bkg_tt_syst_histos_D1[8] = tt_syst_scl_D1;
  }
  RooArgList *bkg_tt_bins_D1 = new RooArgList();
  string procName_D1 = "background_tt_D1_"+year;
  if (shared) 
  {
    RooArgList parlist_D1(N7_tt_D1,a0_tt,a1_tt,d_tt);  // list of shape parameters for tt bkg
    //RooArgList parlist_D1(N7_tt_D1,a0_tt,a1_tt,a2_tt);  // list of shape parameters for tt bkg                                
    construct_formula(procName_D1,*bkg_tt_bins_D1,parlist_D1,bkg_tt_syst_NP_D1,bkg_tt_syst_histos_D1);
  } else 
  {
    RooArgList parlist_D1(N7_tt_D1,a0_tt_D1,a1_tt_D1,d_tt_D1);  // list of shape parameters for tt bkg                       
    //RooArgList parlist_D1(N7_tt_D1,a0_tt_D1,a1_tt_D1,a2_tt_D1);  // list of shape parameters for tt bkg
    construct_formula(procName_D1,*bkg_tt_bins_D1,parlist_D1,bkg_tt_syst_NP_D1,bkg_tt_syst_histos_D1);
  }
  RooParametricHist background_tt_D1(procName_D1.c_str(),"",*wspace->var("CMS_th1x"),*bkg_tt_bins_D1,*data_th1_D1);
  wspace->import(background_tt_D1,RooFit::RecycleConflictNodes());
  stringstream procNameD1Norm;
  procNameD1Norm << procName_D1 << "_norm";
  RooAddition tt_norm_D1(procNameD1Norm.str().c_str(),"",*bkg_tt_bins_D1);
  wspace->import(tt_norm_D1,RooFit::RecycleConflictNodes());

  // ---------------------- MVA bin 2  ------------------

  // Dataset with 8 bins
  RooDataHist data_hist_D2("data_obs_D2","Data observed in MVA bin 2",vars_D2,data_th1_D2);
  wspace->import(data_hist_D2);

  // ttbar bkg in D2
  RooArgList bkg_tt_syst_NP_D2(*wspace->var("np_tt_JEC"),
			       *wspace->var(("np_tt_JER_"+year).c_str()),
			       *wspace->var(("np_tt_btg_"+year).c_str()),
			       *wspace->var(("np_tt_lep_"+year).c_str()),
			       *wspace->var(("np_tt_nom_"+year).c_str()),
			       *wspace->var("np_tt_pdf"),
			       *wspace->var("np_tt_FSR"),
			       *wspace->var("np_tt_isr"));  // list of nuisance parameters for tt bkg
  if (year=="2016") {
    bkg_tt_syst_NP_D2.add(*wspace->var("np_tt_ht_2016"));
  }
  if (year=="2017") {
    bkg_tt_syst_NP_D2.add(*wspace->var("np_tt_scl_2017"));
  }
  TH1D* bkg_tt_syst_histos_D2[20];  // array of histograms containing each tt bkg shape uncertainty
  bkg_tt_syst_histos_D2[0] = tt_syst_JEC_D2;
  bkg_tt_syst_histos_D2[1] = tt_syst_JER_D2;
  bkg_tt_syst_histos_D2[2] = tt_syst_btg_D2;
  bkg_tt_syst_histos_D2[3] = tt_syst_lep_D2;
  bkg_tt_syst_histos_D2[4] = tt_syst_nom_D2;
  bkg_tt_syst_histos_D2[5] = tt_syst_pdf_D2;
  bkg_tt_syst_histos_D2[6] = tt_syst_FSR_D2;
  bkg_tt_syst_histos_D2[7] = tt_syst_isr_D2;
  if (year=="2016") {
    bkg_tt_syst_histos_D2[8] = tt_syst_ht_D2;
  }
  if (year=="2017") {
    bkg_tt_syst_histos_D2[8] = tt_syst_scl_D2;
  }
  RooArgList *bkg_tt_bins_D2 = new RooArgList();
  string procName_D2 = "background_tt_D2_"+year;
  if (shared) 
  {
    RooArgList parlist_D2(N7_tt_D2,a0_tt,a1_tt,d_tt);  // list of shape parameters for tt bkg
    //RooArgList parlist_D2(N7_tt_D2,a0_tt,a1_tt,a2_tt);  // list of shape parameters for tt bkg                                
    construct_formula(procName_D2,*bkg_tt_bins_D2,parlist_D2,bkg_tt_syst_NP_D2,bkg_tt_syst_histos_D2);
  } else 
  {
    RooArgList parlist_D2(N7_tt_D2,a0_tt_D2,a1_tt_D2,d_tt_D2);  // list of shape parameters for tt bkg                       
    //RooArgList parlist_D2(N7_tt_D2,a0_tt_D2,a1_tt_D2,a2_tt_D2);  // list of shape parameters for tt bkg
    construct_formula(procName_D2,*bkg_tt_bins_D2,parlist_D2,bkg_tt_syst_NP_D2,bkg_tt_syst_histos_D2);
  }
  RooParametricHist background_tt_D2(procName_D2.c_str(),"",*wspace->var("CMS_th1x"),*bkg_tt_bins_D2,*data_th1_D2);
  wspace->import(background_tt_D2,RooFit::RecycleConflictNodes());
  stringstream procNameD2Norm;
  procNameD2Norm << procName_D2 << "_norm";
  RooAddition tt_norm_D2(procNameD2Norm.str().c_str(),"",*bkg_tt_bins_D2);
  wspace->import(tt_norm_D2,RooFit::RecycleConflictNodes());

  // ---------------------- MVA bin 3  ------------------

  // Dataset with 8 bins
  RooDataHist data_hist_D3("data_obs_D3","Data observed in MVA bin 3",vars_D3,data_th1_D3);
  wspace->import(data_hist_D3);

  // ttbar bkg in D3
  RooArgList bkg_tt_syst_NP_D3(*wspace->var("np_tt_JEC"),
			       *wspace->var(("np_tt_JER_"+year).c_str()),
			       *wspace->var(("np_tt_btg_"+year).c_str()),
			       *wspace->var(("np_tt_lep_"+year).c_str()),
			       *wspace->var(("np_tt_nom_"+year).c_str()),
			       *wspace->var("np_tt_pdf"),
			       *wspace->var("np_tt_FSR"),
			       *wspace->var("np_tt_isr"));  // list of nuisance parameters for tt bkg
  if (year=="2016") {
    bkg_tt_syst_NP_D3.add(*wspace->var("np_tt_ht_2016"));
  }
  if (year=="2017") {
    bkg_tt_syst_NP_D3.add(*wspace->var("np_tt_scl_2017"));
  }
  TH1D* bkg_tt_syst_histos_D3[20];  // array of histograms containing each tt bkg shape uncertainty
  bkg_tt_syst_histos_D3[0] = tt_syst_JEC_D3;
  bkg_tt_syst_histos_D3[1] = tt_syst_JER_D3;
  bkg_tt_syst_histos_D3[2] = tt_syst_btg_D3;
  bkg_tt_syst_histos_D3[3] = tt_syst_lep_D3;
  bkg_tt_syst_histos_D3[4] = tt_syst_nom_D3;
  bkg_tt_syst_histos_D3[5] = tt_syst_pdf_D3;
  bkg_tt_syst_histos_D3[6] = tt_syst_FSR_D3;
  bkg_tt_syst_histos_D3[7] = tt_syst_isr_D3;
  if (year=="2016") {
    bkg_tt_syst_histos_D3[8] = tt_syst_ht_D3;
  }
  if (year=="2017") {
    bkg_tt_syst_histos_D3[8] = tt_syst_scl_D3;
  }
  RooArgList *bkg_tt_bins_D3 = new RooArgList();
  string procName_D3 = "background_tt_D3_"+year;
  if (shared) 
  {
    RooArgList parlist_D3(N7_tt_D3,a0_tt,a1_tt,d_tt);  // list of shape parameters for tt bkg
    //RooArgList parlist_D3(N7_tt_D3,a0_tt,a1_tt,a2_tt);  // list of shape parameters for tt bkg                                
    construct_formula(procName_D3,*bkg_tt_bins_D3,parlist_D3,bkg_tt_syst_NP_D3,bkg_tt_syst_histos_D3);
  } else 
  {
    RooArgList parlist_D3(N7_tt_D3,a0_tt_D3,a1_tt_D3,d_tt_D3);  // list of shape parameters for tt bkg                       
    //RooArgList parlist_D3(N7_tt_D3,a0_tt_D3,a1_tt_D3,a2_tt_D3);  // list of shape parameters for tt bkg
    construct_formula(procName_D3,*bkg_tt_bins_D3,parlist_D3,bkg_tt_syst_NP_D3,bkg_tt_syst_histos_D3);
  }
  RooParametricHist background_tt_D3(procName_D3.c_str(),"",*wspace->var("CMS_th1x"),*bkg_tt_bins_D3,*data_th1_D3);
  wspace->import(background_tt_D3,RooFit::RecycleConflictNodes());
  stringstream procNameD3Norm;
  procNameD3Norm << procName_D3 << "_norm";
  RooAddition tt_norm_D3(procNameD3Norm.str().c_str(),"",*bkg_tt_bins_D3);
  wspace->import(tt_norm_D3,RooFit::RecycleConflictNodes());

  // ---------------------- MVA bin 4  ------------------

  // Dataset with 8 bins
  RooDataHist data_hist_D4("data_obs_D4","Data observed in MVA bin 4",vars_D4,data_th1_D4);
  wspace->import(data_hist_D4);

  // ttbar bkg in D4
  RooArgList bkg_tt_syst_NP_D4(*wspace->var("np_tt_JEC"),
			       *wspace->var(("np_tt_JER_"+year).c_str()),
			       *wspace->var(("np_tt_btg_"+year).c_str()),
			       *wspace->var(("np_tt_lep_"+year).c_str()),
			       *wspace->var(("np_tt_nom_"+year).c_str()),
			       *wspace->var("np_tt_pdf"),
			       *wspace->var("np_tt_FSR"),
			       *wspace->var("np_tt_isr"));  // list of nuisance parameters for tt bkg
  if (year=="2016") {
    bkg_tt_syst_NP_D4.add(*wspace->var("np_tt_ht_2016"));
  }
  if (year=="2017") {
    bkg_tt_syst_NP_D4.add(*wspace->var("np_tt_scl_2017"));
  }
  TH1D* bkg_tt_syst_histos_D4[20];  // array of histograms containing each tt bkg shape uncertainty
  bkg_tt_syst_histos_D4[0] = tt_syst_JEC_D4;
  bkg_tt_syst_histos_D4[1] = tt_syst_JER_D4;
  bkg_tt_syst_histos_D4[2] = tt_syst_btg_D4;
  bkg_tt_syst_histos_D4[3] = tt_syst_lep_D4;
  bkg_tt_syst_histos_D4[4] = tt_syst_nom_D4;
  bkg_tt_syst_histos_D4[5] = tt_syst_pdf_D4;
  bkg_tt_syst_histos_D4[6] = tt_syst_FSR_D4;
  bkg_tt_syst_histos_D4[7] = tt_syst_isr_D4;
  if (year=="2016") {
    bkg_tt_syst_histos_D4[8] = tt_syst_ht_D4;
  }
  if (year=="2017") {
    bkg_tt_syst_histos_D4[8] = tt_syst_scl_D4;
  }
  RooArgList *bkg_tt_bins_D4 = new RooArgList();
  string procName_D4 = "background_tt_D4_"+year;
  if (shared) 
  {
    RooArgList parlist_D4(N7_tt_D4,a0_tt,a1_tt,d_tt);  // list of shape parameters for tt bkg
    //RooArgList parlist_D4(N7_tt_D4,a0_tt,a1_tt,a2_tt);  // list of shape parameters for tt bkg                                
    construct_formula(procName_D4,*bkg_tt_bins_D4,parlist_D4,bkg_tt_syst_NP_D4,bkg_tt_syst_histos_D4);
  } else 
  {
    RooArgList parlist_D4(N7_tt_D4,a0_tt_D4,a1_tt_D4,d_tt_D4);  // list of shape parameters for tt bkg                       
    //RooArgList parlist_D4(N7_tt_D4,a0_tt_D4,a1_tt_D4,a2_tt_D4);  // list of shape parameters for tt bkg
    construct_formula(procName_D4,*bkg_tt_bins_D4,parlist_D4,bkg_tt_syst_NP_D4,bkg_tt_syst_histos_D4);
  }
  RooParametricHist background_tt_D4(procName_D4.c_str(),"",*wspace->var("CMS_th1x"),*bkg_tt_bins_D4,*data_th1_D4);
  wspace->import(background_tt_D4,RooFit::RecycleConflictNodes());
  stringstream procNameD4Norm;
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

  // =================================================================================
  // Systematics

  // Signal systematics

  TH1D* D1_SIG_JECUp = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_h_njets_pt30_1l_JECUp").c_str());
  TH1D* D2_SIG_JECUp = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_h_njets_pt30_1l_JECUp").c_str());
  TH1D* D3_SIG_JECUp = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_h_njets_pt30_1l_JECUp").c_str());
  TH1D* D4_SIG_JECUp = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_h_njets_pt30_1l_JECUp").c_str());
  D1_SIG_JECUp->SetName("D1_SIG_JECUp");
  D2_SIG_JECUp->SetName("D2_SIG_JECUp");
  D3_SIG_JECUp->SetName("D3_SIG_JECUp");
  D4_SIG_JECUp->SetName("D4_SIG_JECUp");
  D1_SIG_JECUp->Write();
  D2_SIG_JECUp->Write();
  D3_SIG_JECUp->Write();
  D4_SIG_JECUp->Write();
  TH1D* D1_SIG_JECDown = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_h_njets_pt30_1l_JECDown").c_str());
  TH1D* D2_SIG_JECDown = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_h_njets_pt30_1l_JECDown").c_str());
  TH1D* D3_SIG_JECDown = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_h_njets_pt30_1l_JECDown").c_str());
  TH1D* D4_SIG_JECDown = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_h_njets_pt30_1l_JECDown").c_str());
  D1_SIG_JECDown->SetName("D1_SIG_JECDown");
  D2_SIG_JECDown->SetName("D2_SIG_JECDown");
  D3_SIG_JECDown->SetName("D3_SIG_JECDown");
  D4_SIG_JECDown->SetName("D4_SIG_JECDown");
  D1_SIG_JECDown->Write();
  D2_SIG_JECDown->Write();
  D3_SIG_JECDown->Write();
  D4_SIG_JECDown->Write();

  TH1D* D1_SIG_JERUp = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_h_njets_pt30_1l_JERUp").c_str());
  TH1D* D2_SIG_JERUp = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_h_njets_pt30_1l_JERUp").c_str());
  TH1D* D3_SIG_JERUp = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_h_njets_pt30_1l_JERUp").c_str());
  TH1D* D4_SIG_JERUp = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_h_njets_pt30_1l_JERUp").c_str());
  D1_SIG_JERUp->SetName(("D1_SIG_JER_"+year+"Up").c_str());
  D2_SIG_JERUp->SetName(("D2_SIG_JER_"+year+"Up").c_str());
  D3_SIG_JERUp->SetName(("D3_SIG_JER_"+year+"Up").c_str());
  D4_SIG_JERUp->SetName(("D4_SIG_JER_"+year+"Up").c_str());
  D1_SIG_JERUp->Write();
  D2_SIG_JERUp->Write();
  D3_SIG_JERUp->Write();
  D4_SIG_JERUp->Write();
  TH1D* D1_SIG_JERDown = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_h_njets_pt30_1l_JERDown").c_str());
  TH1D* D2_SIG_JERDown = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_h_njets_pt30_1l_JERDown").c_str());
  TH1D* D3_SIG_JERDown = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_h_njets_pt30_1l_JERDown").c_str());
  TH1D* D4_SIG_JERDown = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_h_njets_pt30_1l_JERDown").c_str());
  D1_SIG_JERDown->SetName(("D1_SIG_JER_"+year+"Down").c_str());
  D2_SIG_JERDown->SetName(("D2_SIG_JER_"+year+"Down").c_str());
  D3_SIG_JERDown->SetName(("D3_SIG_JER_"+year+"Down").c_str());
  D4_SIG_JERDown->SetName(("D4_SIG_JER_"+year+"Down").c_str());
  D1_SIG_JERDown->Write();
  D2_SIG_JERDown->Write();
  D3_SIG_JERDown->Write();
  D4_SIG_JERDown->Write();

  TH1D* D1_SIG_btgUp = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_h_njets_pt30_1l_btgUp").c_str());
  TH1D* D2_SIG_btgUp = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_h_njets_pt30_1l_btgUp").c_str());
  TH1D* D3_SIG_btgUp = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_h_njets_pt30_1l_btgUp").c_str());
  TH1D* D4_SIG_btgUp = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_h_njets_pt30_1l_btgUp").c_str());
  D1_SIG_btgUp->SetName(("D1_SIG_btg_"+year+"Up").c_str());
  D2_SIG_btgUp->SetName(("D2_SIG_btg_"+year+"Up").c_str());
  D3_SIG_btgUp->SetName(("D3_SIG_btg_"+year+"Up").c_str());
  D4_SIG_btgUp->SetName(("D4_SIG_btg_"+year+"Up").c_str());
  D1_SIG_btgUp->Write();
  D2_SIG_btgUp->Write();
  D3_SIG_btgUp->Write();
  D4_SIG_btgUp->Write();
  TH1D* D1_SIG_btgDown = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_h_njets_pt30_1l_btgDown").c_str());
  TH1D* D2_SIG_btgDown = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_h_njets_pt30_1l_btgDown").c_str());
  TH1D* D3_SIG_btgDown = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_h_njets_pt30_1l_btgDown").c_str());
  TH1D* D4_SIG_btgDown = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_h_njets_pt30_1l_btgDown").c_str());
  D1_SIG_btgDown->SetName(("D1_SIG_btg_"+year+"Down").c_str());
  D2_SIG_btgDown->SetName(("D2_SIG_btg_"+year+"Down").c_str());
  D3_SIG_btgDown->SetName(("D3_SIG_btg_"+year+"Down").c_str());
  D4_SIG_btgDown->SetName(("D4_SIG_btg_"+year+"Down").c_str());
  D1_SIG_btgDown->Write();
  D2_SIG_btgDown->Write();
  D3_SIG_btgDown->Write();
  D4_SIG_btgDown->Write();

  TH1D* D1_SIG_lepUp = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_h_njets_pt30_1l_lepUp").c_str());
  TH1D* D2_SIG_lepUp = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_h_njets_pt30_1l_lepUp").c_str());
  TH1D* D3_SIG_lepUp = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_h_njets_pt30_1l_lepUp").c_str());
  TH1D* D4_SIG_lepUp = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_h_njets_pt30_1l_lepUp").c_str());
  D1_SIG_lepUp->SetName(("D1_SIG_lep_"+year+"Up").c_str());
  D2_SIG_lepUp->SetName(("D2_SIG_lep_"+year+"Up").c_str());
  D3_SIG_lepUp->SetName(("D3_SIG_lep_"+year+"Up").c_str());
  D4_SIG_lepUp->SetName(("D4_SIG_lep_"+year+"Up").c_str());
  D1_SIG_lepUp->Write();
  D2_SIG_lepUp->Write();
  D3_SIG_lepUp->Write();
  D4_SIG_lepUp->Write();
  TH1D* D1_SIG_lepDown = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_h_njets_pt30_1l_lepDown").c_str());
  TH1D* D2_SIG_lepDown = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_h_njets_pt30_1l_lepDown").c_str());
  TH1D* D3_SIG_lepDown = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_h_njets_pt30_1l_lepDown").c_str());
  TH1D* D4_SIG_lepDown = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_h_njets_pt30_1l_lepDown").c_str());
  D1_SIG_lepDown->SetName(("D1_SIG_lep_"+year+"Down").c_str());
  D2_SIG_lepDown->SetName(("D2_SIG_lep_"+year+"Down").c_str());
  D3_SIG_lepDown->SetName(("D3_SIG_lep_"+year+"Down").c_str());
  D4_SIG_lepDown->SetName(("D4_SIG_lep_"+year+"Down").c_str());
  D1_SIG_lepDown->Write();
  D2_SIG_lepDown->Write();
  D3_SIG_lepDown->Write();
  D4_SIG_lepDown->Write();

  TH1D* D1_SIG_pdfUp = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_h_njets_pt30_1l_pdfUp").c_str());
  TH1D* D2_SIG_pdfUp = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_h_njets_pt30_1l_pdfUp").c_str());
  TH1D* D3_SIG_pdfUp = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_h_njets_pt30_1l_pdfUp").c_str());
  TH1D* D4_SIG_pdfUp = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_h_njets_pt30_1l_pdfUp").c_str());
  D1_SIG_pdfUp->SetName("D1_SIG_pdfUp");
  D2_SIG_pdfUp->SetName("D2_SIG_pdfUp");
  D3_SIG_pdfUp->SetName("D3_SIG_pdfUp");
  D4_SIG_pdfUp->SetName("D4_SIG_pdfUp");
  D1_SIG_pdfUp->Write();
  D2_SIG_pdfUp->Write();
  D3_SIG_pdfUp->Write();
  D4_SIG_pdfUp->Write();
  TH1D* D1_SIG_pdfDown = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_h_njets_pt30_1l_pdfDown").c_str());
  TH1D* D2_SIG_pdfDown = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_h_njets_pt30_1l_pdfDown").c_str());
  TH1D* D3_SIG_pdfDown = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_h_njets_pt30_1l_pdfDown").c_str());
  TH1D* D4_SIG_pdfDown = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_h_njets_pt30_1l_pdfDown").c_str());
  D1_SIG_pdfDown->SetName("D1_SIG_pdfDown");
  D2_SIG_pdfDown->SetName("D2_SIG_pdfDown");
  D3_SIG_pdfDown->SetName("D3_SIG_pdfDown");
  D4_SIG_pdfDown->SetName("D4_SIG_pdfDown");
  D1_SIG_pdfDown->Write();
  D2_SIG_pdfDown->Write();
  D3_SIG_pdfDown->Write();
  D4_SIG_pdfDown->Write();

  TH1D* D1_SIG_sclUp;
  TH1D* D2_SIG_sclUp;
  TH1D* D3_SIG_sclUp;
  TH1D* D4_SIG_sclUp;
  TH1D* D1_SIG_sclDown;
  TH1D* D2_SIG_sclDown;
  TH1D* D3_SIG_sclDown;
  TH1D* D4_SIG_sclDown;
  if (year=="2017") {
    D1_SIG_sclUp = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_h_njets_pt30_1l_sclUp").c_str());
    D2_SIG_sclUp = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_h_njets_pt30_1l_sclUp").c_str());
    D3_SIG_sclUp = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_h_njets_pt30_1l_sclUp").c_str());
    D4_SIG_sclUp = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_h_njets_pt30_1l_sclUp").c_str());
    D1_SIG_sclUp->SetName("D1_SIG_scl_2017Up");
    D2_SIG_sclUp->SetName("D2_SIG_scl_2017Up");
    D3_SIG_sclUp->SetName("D3_SIG_scl_2017Up");
    D4_SIG_sclUp->SetName("D4_SIG_scl_2017Up");
    D1_SIG_sclUp->Write();
    D2_SIG_sclUp->Write();
    D3_SIG_sclUp->Write();
    D4_SIG_sclUp->Write();
    D1_SIG_sclDown = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_h_njets_pt30_1l_sclDown").c_str());
    D2_SIG_sclDown = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_h_njets_pt30_1l_sclDown").c_str());
    D3_SIG_sclDown = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_h_njets_pt30_1l_sclDown").c_str());
    D4_SIG_sclDown = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_h_njets_pt30_1l_sclDown").c_str());
    D1_SIG_sclDown->SetName("D1_SIG_scl_2017Down");
    D2_SIG_sclDown->SetName("D2_SIG_scl_2017Down");
    D3_SIG_sclDown->SetName("D3_SIG_scl_2017Down");
    D4_SIG_sclDown->SetName("D4_SIG_scl_2017Down");
    D1_SIG_sclDown->Write();
    D2_SIG_sclDown->Write();
    D3_SIG_sclDown->Write();
    D4_SIG_sclDown->Write();
  }


  // "OTHER" background systematics

  TH1D* D1_OTHER_JECUp = (TH1D*)file->Get("D1_OTHER_h_njets_pt30_1l_JECUp");
  TH1D* D2_OTHER_JECUp = (TH1D*)file->Get("D2_OTHER_h_njets_pt30_1l_JECUp");
  TH1D* D3_OTHER_JECUp = (TH1D*)file->Get("D3_OTHER_h_njets_pt30_1l_JECUp");
  TH1D* D4_OTHER_JECUp = (TH1D*)file->Get("D4_OTHER_h_njets_pt30_1l_JECUp");
  D1_OTHER_JECUp->SetName("D1_OTHER_JECUp");
  D2_OTHER_JECUp->SetName("D2_OTHER_JECUp");
  D3_OTHER_JECUp->SetName("D3_OTHER_JECUp");
  D4_OTHER_JECUp->SetName("D4_OTHER_JECUp");
  D1_OTHER_JECUp->Write();
  D2_OTHER_JECUp->Write();
  D3_OTHER_JECUp->Write();
  D4_OTHER_JECUp->Write();
  TH1D* D1_OTHER_JECDown = (TH1D*)file->Get("D1_OTHER_h_njets_pt30_1l_JECDown");
  TH1D* D2_OTHER_JECDown = (TH1D*)file->Get("D2_OTHER_h_njets_pt30_1l_JECDown");
  TH1D* D3_OTHER_JECDown = (TH1D*)file->Get("D3_OTHER_h_njets_pt30_1l_JECDown");
  TH1D* D4_OTHER_JECDown = (TH1D*)file->Get("D4_OTHER_h_njets_pt30_1l_JECDown");
  D1_OTHER_JECDown->SetName("D1_OTHER_JECDown");
  D2_OTHER_JECDown->SetName("D2_OTHER_JECDown");
  D3_OTHER_JECDown->SetName("D3_OTHER_JECDown");
  D4_OTHER_JECDown->SetName("D4_OTHER_JECDown");
  D1_OTHER_JECDown->Write();
  D2_OTHER_JECDown->Write();
  D3_OTHER_JECDown->Write();
  D4_OTHER_JECDown->Write();

  TH1D* D1_OTHER_JERUp = (TH1D*)file->Get("D1_OTHER_h_njets_pt30_1l_JERUp");
  TH1D* D2_OTHER_JERUp = (TH1D*)file->Get("D2_OTHER_h_njets_pt30_1l_JERUp");
  TH1D* D3_OTHER_JERUp = (TH1D*)file->Get("D3_OTHER_h_njets_pt30_1l_JERUp");
  TH1D* D4_OTHER_JERUp = (TH1D*)file->Get("D4_OTHER_h_njets_pt30_1l_JERUp");
  D1_OTHER_JERUp->SetName(("D1_OTHER_JER_"+year+"Up").c_str());
  D2_OTHER_JERUp->SetName(("D2_OTHER_JER_"+year+"Up").c_str());
  D3_OTHER_JERUp->SetName(("D3_OTHER_JER_"+year+"Up").c_str());
  D4_OTHER_JERUp->SetName(("D4_OTHER_JER_"+year+"Up").c_str());
  D1_OTHER_JERUp->Write();
  D2_OTHER_JERUp->Write();
  D3_OTHER_JERUp->Write();
  D4_OTHER_JERUp->Write();
  TH1D* D1_OTHER_JERDown = (TH1D*)file->Get("D1_OTHER_h_njets_pt30_1l_JERDown");
  TH1D* D2_OTHER_JERDown = (TH1D*)file->Get("D2_OTHER_h_njets_pt30_1l_JERDown");
  TH1D* D3_OTHER_JERDown = (TH1D*)file->Get("D3_OTHER_h_njets_pt30_1l_JERDown");
  TH1D* D4_OTHER_JERDown = (TH1D*)file->Get("D4_OTHER_h_njets_pt30_1l_JERDown");
  D1_OTHER_JERDown->SetName(("D1_OTHER_JER_"+year+"Down").c_str());
  D2_OTHER_JERDown->SetName(("D2_OTHER_JER_"+year+"Down").c_str());
  D3_OTHER_JERDown->SetName(("D3_OTHER_JER_"+year+"Down").c_str());
  D4_OTHER_JERDown->SetName(("D4_OTHER_JER_"+year+"Down").c_str());
  D1_OTHER_JERDown->Write();
  D2_OTHER_JERDown->Write();
  D3_OTHER_JERDown->Write();
  D4_OTHER_JERDown->Write();

  TH1D* D1_OTHER_btgUp = (TH1D*)file->Get("D1_OTHER_h_njets_pt30_1l_btgUp");
  TH1D* D2_OTHER_btgUp = (TH1D*)file->Get("D2_OTHER_h_njets_pt30_1l_btgUp");
  TH1D* D3_OTHER_btgUp = (TH1D*)file->Get("D3_OTHER_h_njets_pt30_1l_btgUp");
  TH1D* D4_OTHER_btgUp = (TH1D*)file->Get("D4_OTHER_h_njets_pt30_1l_btgUp");
  D1_OTHER_btgUp->SetName(("D1_OTHER_btg_"+year+"Up").c_str());
  D2_OTHER_btgUp->SetName(("D2_OTHER_btg_"+year+"Up").c_str());
  D3_OTHER_btgUp->SetName(("D3_OTHER_btg_"+year+"Up").c_str());
  D4_OTHER_btgUp->SetName(("D4_OTHER_btg_"+year+"Up").c_str());
  D1_OTHER_btgUp->Write();
  D2_OTHER_btgUp->Write();
  D3_OTHER_btgUp->Write();
  D4_OTHER_btgUp->Write();
  TH1D* D1_OTHER_btgDown = (TH1D*)file->Get("D1_OTHER_h_njets_pt30_1l_btgDown");
  TH1D* D2_OTHER_btgDown = (TH1D*)file->Get("D2_OTHER_h_njets_pt30_1l_btgDown");
  TH1D* D3_OTHER_btgDown = (TH1D*)file->Get("D3_OTHER_h_njets_pt30_1l_btgDown");
  TH1D* D4_OTHER_btgDown = (TH1D*)file->Get("D4_OTHER_h_njets_pt30_1l_btgDown");
  D1_OTHER_btgDown->SetName(("D1_OTHER_btg_"+year+"Down").c_str());
  D2_OTHER_btgDown->SetName(("D2_OTHER_btg_"+year+"Down").c_str());
  D3_OTHER_btgDown->SetName(("D3_OTHER_btg_"+year+"Down").c_str());
  D4_OTHER_btgDown->SetName(("D4_OTHER_btg_"+year+"Down").c_str());
  D1_OTHER_btgDown->Write();
  D2_OTHER_btgDown->Write();
  D3_OTHER_btgDown->Write();
  D4_OTHER_btgDown->Write();

  TH1D* D1_OTHER_lepUp = (TH1D*)file->Get("D1_OTHER_h_njets_pt30_1l_lepUp");
  TH1D* D2_OTHER_lepUp = (TH1D*)file->Get("D2_OTHER_h_njets_pt30_1l_lepUp");
  TH1D* D3_OTHER_lepUp = (TH1D*)file->Get("D3_OTHER_h_njets_pt30_1l_lepUp");
  TH1D* D4_OTHER_lepUp = (TH1D*)file->Get("D4_OTHER_h_njets_pt30_1l_lepUp");
  D1_OTHER_lepUp->SetName(("D1_OTHER_lep_"+year+"Up").c_str());
  D2_OTHER_lepUp->SetName(("D2_OTHER_lep_"+year+"Up").c_str());
  D3_OTHER_lepUp->SetName(("D3_OTHER_lep_"+year+"Up").c_str());
  D4_OTHER_lepUp->SetName(("D4_OTHER_lep_"+year+"Up").c_str());
  D1_OTHER_lepUp->Write();
  D2_OTHER_lepUp->Write();
  D3_OTHER_lepUp->Write();
  D4_OTHER_lepUp->Write();
  TH1D* D1_OTHER_lepDown = (TH1D*)file->Get("D1_OTHER_h_njets_pt30_1l_lepDown");
  TH1D* D2_OTHER_lepDown = (TH1D*)file->Get("D2_OTHER_h_njets_pt30_1l_lepDown");
  TH1D* D3_OTHER_lepDown = (TH1D*)file->Get("D3_OTHER_h_njets_pt30_1l_lepDown");
  TH1D* D4_OTHER_lepDown = (TH1D*)file->Get("D4_OTHER_h_njets_pt30_1l_lepDown");
  D1_OTHER_lepDown->SetName(("D1_OTHER_lep_"+year+"Down").c_str());
  D2_OTHER_lepDown->SetName(("D2_OTHER_lep_"+year+"Down").c_str());
  D3_OTHER_lepDown->SetName(("D3_OTHER_lep_"+year+"Down").c_str());
  D4_OTHER_lepDown->SetName(("D4_OTHER_lep_"+year+"Down").c_str());
  D1_OTHER_lepDown->Write();
  D2_OTHER_lepDown->Write();
  D3_OTHER_lepDown->Write();
  D4_OTHER_lepDown->Write();

  TH1D* D1_OTHER_htUp;
  TH1D* D2_OTHER_htUp;
  TH1D* D3_OTHER_htUp;
  TH1D* D4_OTHER_htUp;
  TH1D* D1_OTHER_htDown;
  TH1D* D2_OTHER_htDown;
  TH1D* D3_OTHER_htDown;
  TH1D* D4_OTHER_htDown;
  if (year=="2016") {
    D1_OTHER_htUp = (TH1D*)file->Get("D1_OTHER_h_njets_pt30_1l_htUp");
    D2_OTHER_htUp = (TH1D*)file->Get("D2_OTHER_h_njets_pt30_1l_htUp");
    D3_OTHER_htUp = (TH1D*)file->Get("D3_OTHER_h_njets_pt30_1l_htUp");
    D4_OTHER_htUp = (TH1D*)file->Get("D4_OTHER_h_njets_pt30_1l_htUp");
    D1_OTHER_htUp->SetName("D1_OTHER_ht_2016Up");
    D2_OTHER_htUp->SetName("D2_OTHER_ht_2016Up");
    D3_OTHER_htUp->SetName("D3_OTHER_ht_2016Up");
    D4_OTHER_htUp->SetName("D4_OTHER_ht_2016Up");
    D1_OTHER_htUp->Write();
    D2_OTHER_htUp->Write();
    D3_OTHER_htUp->Write();
    D4_OTHER_htUp->Write();
    D1_OTHER_htDown = (TH1D*)file->Get("D1_OTHER_h_njets_pt30_1l_htDown");
    D2_OTHER_htDown = (TH1D*)file->Get("D2_OTHER_h_njets_pt30_1l_htDown");
    D3_OTHER_htDown = (TH1D*)file->Get("D3_OTHER_h_njets_pt30_1l_htDown");
    D4_OTHER_htDown = (TH1D*)file->Get("D4_OTHER_h_njets_pt30_1l_htDown");
    D1_OTHER_htDown->SetName("D1_OTHER_ht_2016Down");
    D2_OTHER_htDown->SetName("D2_OTHER_ht_2016Down");
    D3_OTHER_htDown->SetName("D3_OTHER_ht_2016Down");
    D4_OTHER_htDown->SetName("D4_OTHER_ht_2016Down");
    D1_OTHER_htDown->Write();
    D2_OTHER_htDown->Write();
    D3_OTHER_htDown->Write();
    D4_OTHER_htDown->Write();
  }

  TH1D* D1_OTHER_sclUp;
  TH1D* D2_OTHER_sclUp;
  TH1D* D3_OTHER_sclUp;
  TH1D* D4_OTHER_sclUp;
  TH1D* D1_OTHER_sclDown;
  TH1D* D2_OTHER_sclDown;
  TH1D* D3_OTHER_sclDown;
  TH1D* D4_OTHER_sclDown;
  if (year=="2017") {
    D1_OTHER_sclUp = (TH1D*)file->Get("D1_OTHER_h_njets_pt30_1l_sclUp");
    D2_OTHER_sclUp = (TH1D*)file->Get("D2_OTHER_h_njets_pt30_1l_sclUp");
    D3_OTHER_sclUp = (TH1D*)file->Get("D3_OTHER_h_njets_pt30_1l_sclUp");
    D4_OTHER_sclUp = (TH1D*)file->Get("D4_OTHER_h_njets_pt30_1l_sclUp");
    D1_OTHER_sclUp->SetName("D1_OTHER_scl_2017Up");
    D2_OTHER_sclUp->SetName("D2_OTHER_scl_2017Up");
    D3_OTHER_sclUp->SetName("D3_OTHER_scl_2017Up");
    D4_OTHER_sclUp->SetName("D4_OTHER_scl_2017Up");
    D1_OTHER_sclUp->Write();
    D2_OTHER_sclUp->Write();
    D3_OTHER_sclUp->Write();
    D4_OTHER_sclUp->Write();
    D1_OTHER_sclDown = (TH1D*)file->Get("D1_OTHER_h_njets_pt30_1l_sclDown");
    D2_OTHER_sclDown = (TH1D*)file->Get("D2_OTHER_h_njets_pt30_1l_sclDown");
    D3_OTHER_sclDown = (TH1D*)file->Get("D3_OTHER_h_njets_pt30_1l_sclDown");
    D4_OTHER_sclDown = (TH1D*)file->Get("D4_OTHER_h_njets_pt30_1l_sclDown");
    D1_OTHER_sclDown->SetName("D1_OTHER_scl_2017Down");
    D2_OTHER_sclDown->SetName("D2_OTHER_scl_2017Down");
    D3_OTHER_sclDown->SetName("D3_OTHER_scl_2017Down");
    D4_OTHER_sclDown->SetName("D4_OTHER_scl_2017Down");
    D1_OTHER_sclDown->Write();
    D2_OTHER_sclDown->Write();
    D3_OTHER_sclDown->Write();
    D4_OTHER_sclDown->Write();
  }



  // =================================================================================
  // Statistics-based Uncertainties

  // MC stat uncertainty histograms for the particular signal model and mass point
  TH1D* D1_SIG_mcStatBin1Up = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_mcStatBin1Up").c_str());
  TH1D* D1_SIG_mcStatBin2Up = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_mcStatBin2Up").c_str());
  TH1D* D1_SIG_mcStatBin3Up = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_mcStatBin3Up").c_str());
  TH1D* D1_SIG_mcStatBin4Up = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_mcStatBin4Up").c_str());
  TH1D* D1_SIG_mcStatBin5Up = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_mcStatBin5Up").c_str());
  TH1D* D1_SIG_mcStatBin6Up = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_mcStatBin6Up").c_str());
  TH1D* D1_SIG_mcStatBin7Up = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_mcStatBin7Up").c_str());
  TH1D* D1_SIG_mcStatBin8Up = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_mcStatBin8Up").c_str());
  D1_SIG_mcStatBin1Up->SetName(("D1_SIG_mcStatD1SIGBin1_"+year+"Up").c_str());
  D1_SIG_mcStatBin2Up->SetName(("D1_SIG_mcStatD1SIGBin2_"+year+"Up").c_str());
  D1_SIG_mcStatBin3Up->SetName(("D1_SIG_mcStatD1SIGBin3_"+year+"Up").c_str());
  D1_SIG_mcStatBin4Up->SetName(("D1_SIG_mcStatD1SIGBin4_"+year+"Up").c_str());
  D1_SIG_mcStatBin5Up->SetName(("D1_SIG_mcStatD1SIGBin5_"+year+"Up").c_str());
  D1_SIG_mcStatBin6Up->SetName(("D1_SIG_mcStatD1SIGBin6_"+year+"Up").c_str());
  D1_SIG_mcStatBin7Up->SetName(("D1_SIG_mcStatD1SIGBin7_"+year+"Up").c_str());
  D1_SIG_mcStatBin8Up->SetName(("D1_SIG_mcStatD1SIGBin8_"+year+"Up").c_str());
  D1_SIG_mcStatBin1Up->Write();
  D1_SIG_mcStatBin2Up->Write();
  D1_SIG_mcStatBin3Up->Write();
  D1_SIG_mcStatBin4Up->Write();
  D1_SIG_mcStatBin5Up->Write();
  D1_SIG_mcStatBin6Up->Write();
  D1_SIG_mcStatBin7Up->Write();
  D1_SIG_mcStatBin8Up->Write();
  TH1D* D1_SIG_mcStatBin1Down = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_mcStatBin1Down").c_str());
  TH1D* D1_SIG_mcStatBin2Down = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_mcStatBin2Down").c_str());
  TH1D* D1_SIG_mcStatBin3Down = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_mcStatBin3Down").c_str());
  TH1D* D1_SIG_mcStatBin4Down = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_mcStatBin4Down").c_str());
  TH1D* D1_SIG_mcStatBin5Down = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_mcStatBin5Down").c_str());
  TH1D* D1_SIG_mcStatBin6Down = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_mcStatBin6Down").c_str());
  TH1D* D1_SIG_mcStatBin7Down = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_mcStatBin7Down").c_str());
  TH1D* D1_SIG_mcStatBin8Down = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_mcStatBin8Down").c_str());
  D1_SIG_mcStatBin1Down->SetName(("D1_SIG_mcStatD1SIGBin1_"+year+"Down").c_str());
  D1_SIG_mcStatBin2Down->SetName(("D1_SIG_mcStatD1SIGBin2_"+year+"Down").c_str());
  D1_SIG_mcStatBin3Down->SetName(("D1_SIG_mcStatD1SIGBin3_"+year+"Down").c_str());
  D1_SIG_mcStatBin4Down->SetName(("D1_SIG_mcStatD1SIGBin4_"+year+"Down").c_str());
  D1_SIG_mcStatBin5Down->SetName(("D1_SIG_mcStatD1SIGBin5_"+year+"Down").c_str());
  D1_SIG_mcStatBin6Down->SetName(("D1_SIG_mcStatD1SIGBin6_"+year+"Down").c_str());
  D1_SIG_mcStatBin7Down->SetName(("D1_SIG_mcStatD1SIGBin7_"+year+"Down").c_str());
  D1_SIG_mcStatBin8Down->SetName(("D1_SIG_mcStatD1SIGBin8_"+year+"Down").c_str());
  D1_SIG_mcStatBin1Down->Write();
  D1_SIG_mcStatBin2Down->Write();
  D1_SIG_mcStatBin3Down->Write();
  D1_SIG_mcStatBin4Down->Write();
  D1_SIG_mcStatBin5Down->Write();
  D1_SIG_mcStatBin6Down->Write();
  D1_SIG_mcStatBin7Down->Write();
  D1_SIG_mcStatBin8Down->Write();

  TH1D* D2_SIG_mcStatBin1Up = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_mcStatBin1Up").c_str());
  TH1D* D2_SIG_mcStatBin2Up = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_mcStatBin2Up").c_str());
  TH1D* D2_SIG_mcStatBin3Up = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_mcStatBin3Up").c_str());
  TH1D* D2_SIG_mcStatBin4Up = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_mcStatBin4Up").c_str());
  TH1D* D2_SIG_mcStatBin5Up = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_mcStatBin5Up").c_str());
  TH1D* D2_SIG_mcStatBin6Up = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_mcStatBin6Up").c_str());
  TH1D* D2_SIG_mcStatBin7Up = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_mcStatBin7Up").c_str());
  TH1D* D2_SIG_mcStatBin8Up = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_mcStatBin8Up").c_str());
  D2_SIG_mcStatBin1Up->SetName(("D2_SIG_mcStatD2SIGBin1_"+year+"Up").c_str());
  D2_SIG_mcStatBin2Up->SetName(("D2_SIG_mcStatD2SIGBin2_"+year+"Up").c_str());
  D2_SIG_mcStatBin3Up->SetName(("D2_SIG_mcStatD2SIGBin3_"+year+"Up").c_str());
  D2_SIG_mcStatBin4Up->SetName(("D2_SIG_mcStatD2SIGBin4_"+year+"Up").c_str());
  D2_SIG_mcStatBin5Up->SetName(("D2_SIG_mcStatD2SIGBin5_"+year+"Up").c_str());
  D2_SIG_mcStatBin6Up->SetName(("D2_SIG_mcStatD2SIGBin6_"+year+"Up").c_str());
  D2_SIG_mcStatBin7Up->SetName(("D2_SIG_mcStatD2SIGBin7_"+year+"Up").c_str());
  D2_SIG_mcStatBin8Up->SetName(("D2_SIG_mcStatD2SIGBin8_"+year+"Up").c_str());
  D2_SIG_mcStatBin1Up->Write();
  D2_SIG_mcStatBin2Up->Write();
  D2_SIG_mcStatBin3Up->Write();
  D2_SIG_mcStatBin4Up->Write();
  D2_SIG_mcStatBin5Up->Write();
  D2_SIG_mcStatBin6Up->Write();
  D2_SIG_mcStatBin7Up->Write();
  D2_SIG_mcStatBin8Up->Write();
  TH1D* D2_SIG_mcStatBin1Down = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_mcStatBin1Down").c_str());
  TH1D* D2_SIG_mcStatBin2Down = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_mcStatBin2Down").c_str());
  TH1D* D2_SIG_mcStatBin3Down = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_mcStatBin3Down").c_str());
  TH1D* D2_SIG_mcStatBin4Down = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_mcStatBin4Down").c_str());
  TH1D* D2_SIG_mcStatBin5Down = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_mcStatBin5Down").c_str());
  TH1D* D2_SIG_mcStatBin6Down = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_mcStatBin6Down").c_str());
  TH1D* D2_SIG_mcStatBin7Down = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_mcStatBin7Down").c_str());
  TH1D* D2_SIG_mcStatBin8Down = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_mcStatBin8Down").c_str());
  D2_SIG_mcStatBin1Down->SetName(("D2_SIG_mcStatD2SIGBin1_"+year+"Down").c_str());
  D2_SIG_mcStatBin2Down->SetName(("D2_SIG_mcStatD2SIGBin2_"+year+"Down").c_str());
  D2_SIG_mcStatBin3Down->SetName(("D2_SIG_mcStatD2SIGBin3_"+year+"Down").c_str());
  D2_SIG_mcStatBin4Down->SetName(("D2_SIG_mcStatD2SIGBin4_"+year+"Down").c_str());
  D2_SIG_mcStatBin5Down->SetName(("D2_SIG_mcStatD2SIGBin5_"+year+"Down").c_str());
  D2_SIG_mcStatBin6Down->SetName(("D2_SIG_mcStatD2SIGBin6_"+year+"Down").c_str());
  D2_SIG_mcStatBin7Down->SetName(("D2_SIG_mcStatD2SIGBin7_"+year+"Down").c_str());
  D2_SIG_mcStatBin8Down->SetName(("D2_SIG_mcStatD2SIGBin8_"+year+"Down").c_str());
  D2_SIG_mcStatBin1Down->Write();
  D2_SIG_mcStatBin2Down->Write();
  D2_SIG_mcStatBin3Down->Write();
  D2_SIG_mcStatBin4Down->Write();
  D2_SIG_mcStatBin5Down->Write();
  D2_SIG_mcStatBin6Down->Write();
  D2_SIG_mcStatBin7Down->Write();
  D2_SIG_mcStatBin8Down->Write();

  TH1D* D3_SIG_mcStatBin1Up = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_mcStatBin1Up").c_str());
  TH1D* D3_SIG_mcStatBin2Up = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_mcStatBin2Up").c_str());
  TH1D* D3_SIG_mcStatBin3Up = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_mcStatBin3Up").c_str());
  TH1D* D3_SIG_mcStatBin4Up = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_mcStatBin4Up").c_str());
  TH1D* D3_SIG_mcStatBin5Up = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_mcStatBin5Up").c_str());
  TH1D* D3_SIG_mcStatBin6Up = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_mcStatBin6Up").c_str());
  TH1D* D3_SIG_mcStatBin7Up = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_mcStatBin7Up").c_str());
  TH1D* D3_SIG_mcStatBin8Up = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_mcStatBin8Up").c_str());
  D3_SIG_mcStatBin1Up->SetName(("D3_SIG_mcStatD3SIGBin1_"+year+"Up").c_str());
  D3_SIG_mcStatBin2Up->SetName(("D3_SIG_mcStatD3SIGBin2_"+year+"Up").c_str());
  D3_SIG_mcStatBin3Up->SetName(("D3_SIG_mcStatD3SIGBin3_"+year+"Up").c_str());
  D3_SIG_mcStatBin4Up->SetName(("D3_SIG_mcStatD3SIGBin4_"+year+"Up").c_str());
  D3_SIG_mcStatBin5Up->SetName(("D3_SIG_mcStatD3SIGBin5_"+year+"Up").c_str());
  D3_SIG_mcStatBin6Up->SetName(("D3_SIG_mcStatD3SIGBin6_"+year+"Up").c_str());
  D3_SIG_mcStatBin7Up->SetName(("D3_SIG_mcStatD3SIGBin7_"+year+"Up").c_str());
  D3_SIG_mcStatBin8Up->SetName(("D3_SIG_mcStatD3SIGBin8_"+year+"Up").c_str());
  D3_SIG_mcStatBin1Up->Write();
  D3_SIG_mcStatBin2Up->Write();
  D3_SIG_mcStatBin3Up->Write();
  D3_SIG_mcStatBin4Up->Write();
  D3_SIG_mcStatBin5Up->Write();
  D3_SIG_mcStatBin6Up->Write();
  D3_SIG_mcStatBin7Up->Write();
  D3_SIG_mcStatBin8Up->Write();
  TH1D* D3_SIG_mcStatBin1Down = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_mcStatBin1Down").c_str());
  TH1D* D3_SIG_mcStatBin2Down = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_mcStatBin2Down").c_str());
  TH1D* D3_SIG_mcStatBin3Down = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_mcStatBin3Down").c_str());
  TH1D* D3_SIG_mcStatBin4Down = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_mcStatBin4Down").c_str());
  TH1D* D3_SIG_mcStatBin5Down = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_mcStatBin5Down").c_str());
  TH1D* D3_SIG_mcStatBin6Down = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_mcStatBin6Down").c_str());
  TH1D* D3_SIG_mcStatBin7Down = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_mcStatBin7Down").c_str());
  TH1D* D3_SIG_mcStatBin8Down = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_mcStatBin8Down").c_str());
  D3_SIG_mcStatBin1Down->SetName(("D3_SIG_mcStatD3SIGBin1_"+year+"Down").c_str());
  D3_SIG_mcStatBin2Down->SetName(("D3_SIG_mcStatD3SIGBin2_"+year+"Down").c_str());
  D3_SIG_mcStatBin3Down->SetName(("D3_SIG_mcStatD3SIGBin3_"+year+"Down").c_str());
  D3_SIG_mcStatBin4Down->SetName(("D3_SIG_mcStatD3SIGBin4_"+year+"Down").c_str());
  D3_SIG_mcStatBin5Down->SetName(("D3_SIG_mcStatD3SIGBin5_"+year+"Down").c_str());
  D3_SIG_mcStatBin6Down->SetName(("D3_SIG_mcStatD3SIGBin6_"+year+"Down").c_str());
  D3_SIG_mcStatBin7Down->SetName(("D3_SIG_mcStatD3SIGBin7_"+year+"Down").c_str());
  D3_SIG_mcStatBin8Down->SetName(("D3_SIG_mcStatD3SIGBin8_"+year+"Down").c_str());
  D3_SIG_mcStatBin1Down->Write();
  D3_SIG_mcStatBin2Down->Write();
  D3_SIG_mcStatBin3Down->Write();
  D3_SIG_mcStatBin4Down->Write();
  D3_SIG_mcStatBin5Down->Write();
  D3_SIG_mcStatBin6Down->Write();
  D3_SIG_mcStatBin7Down->Write();
  D3_SIG_mcStatBin8Down->Write();

  TH1D* D4_SIG_mcStatBin1Up = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_mcStatBin1Up").c_str());
  TH1D* D4_SIG_mcStatBin2Up = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_mcStatBin2Up").c_str());
  TH1D* D4_SIG_mcStatBin3Up = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_mcStatBin3Up").c_str());
  TH1D* D4_SIG_mcStatBin4Up = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_mcStatBin4Up").c_str());
  TH1D* D4_SIG_mcStatBin5Up = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_mcStatBin5Up").c_str());
  TH1D* D4_SIG_mcStatBin6Up = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_mcStatBin6Up").c_str());
  TH1D* D4_SIG_mcStatBin7Up = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_mcStatBin7Up").c_str());
  TH1D* D4_SIG_mcStatBin8Up = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_mcStatBin8Up").c_str());
  D4_SIG_mcStatBin1Up->SetName(("D4_SIG_mcStatD4SIGBin1_"+year+"Up").c_str());
  D4_SIG_mcStatBin2Up->SetName(("D4_SIG_mcStatD4SIGBin2_"+year+"Up").c_str());
  D4_SIG_mcStatBin3Up->SetName(("D4_SIG_mcStatD4SIGBin3_"+year+"Up").c_str());
  D4_SIG_mcStatBin4Up->SetName(("D4_SIG_mcStatD4SIGBin4_"+year+"Up").c_str());
  D4_SIG_mcStatBin5Up->SetName(("D4_SIG_mcStatD4SIGBin5_"+year+"Up").c_str());
  D4_SIG_mcStatBin6Up->SetName(("D4_SIG_mcStatD4SIGBin6_"+year+"Up").c_str());
  D4_SIG_mcStatBin7Up->SetName(("D4_SIG_mcStatD4SIGBin7_"+year+"Up").c_str());
  D4_SIG_mcStatBin8Up->SetName(("D4_SIG_mcStatD4SIGBin8_"+year+"Up").c_str());
  D4_SIG_mcStatBin1Up->Write();
  D4_SIG_mcStatBin2Up->Write();
  D4_SIG_mcStatBin3Up->Write();
  D4_SIG_mcStatBin4Up->Write();
  D4_SIG_mcStatBin5Up->Write();
  D4_SIG_mcStatBin6Up->Write();
  D4_SIG_mcStatBin7Up->Write();
  D4_SIG_mcStatBin8Up->Write();
  TH1D* D4_SIG_mcStatBin1Down = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_mcStatBin1Down").c_str());
  TH1D* D4_SIG_mcStatBin2Down = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_mcStatBin2Down").c_str());
  TH1D* D4_SIG_mcStatBin3Down = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_mcStatBin3Down").c_str());
  TH1D* D4_SIG_mcStatBin4Down = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_mcStatBin4Down").c_str());
  TH1D* D4_SIG_mcStatBin5Down = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_mcStatBin5Down").c_str());
  TH1D* D4_SIG_mcStatBin6Down = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_mcStatBin6Down").c_str());
  TH1D* D4_SIG_mcStatBin7Down = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_mcStatBin7Down").c_str());
  TH1D* D4_SIG_mcStatBin8Down = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_mcStatBin8Down").c_str());
  D4_SIG_mcStatBin1Down->SetName(("D4_SIG_mcStatD4SIGBin1_"+year+"Down").c_str());
  D4_SIG_mcStatBin2Down->SetName(("D4_SIG_mcStatD4SIGBin2_"+year+"Down").c_str());
  D4_SIG_mcStatBin3Down->SetName(("D4_SIG_mcStatD4SIGBin3_"+year+"Down").c_str());
  D4_SIG_mcStatBin4Down->SetName(("D4_SIG_mcStatD4SIGBin4_"+year+"Down").c_str());
  D4_SIG_mcStatBin5Down->SetName(("D4_SIG_mcStatD4SIGBin5_"+year+"Down").c_str());
  D4_SIG_mcStatBin6Down->SetName(("D4_SIG_mcStatD4SIGBin6_"+year+"Down").c_str());
  D4_SIG_mcStatBin7Down->SetName(("D4_SIG_mcStatD4SIGBin7_"+year+"Down").c_str());
  D4_SIG_mcStatBin8Down->SetName(("D4_SIG_mcStatD4SIGBin8_"+year+"Down").c_str());
  D4_SIG_mcStatBin1Down->Write();
  D4_SIG_mcStatBin2Down->Write();
  D4_SIG_mcStatBin3Down->Write();
  D4_SIG_mcStatBin4Down->Write();
  D4_SIG_mcStatBin5Down->Write();
  D4_SIG_mcStatBin6Down->Write();
  D4_SIG_mcStatBin7Down->Write();
  D4_SIG_mcStatBin8Down->Write();


  // MC stat uncertainty histograms for OTHER backgrounds
  TH1D* D1_OTHER_mcStatBin1Up = (TH1D*)file->Get("D1_OTHER_mcStatBin1Up");
  TH1D* D1_OTHER_mcStatBin2Up = (TH1D*)file->Get("D1_OTHER_mcStatBin2Up");
  TH1D* D1_OTHER_mcStatBin3Up = (TH1D*)file->Get("D1_OTHER_mcStatBin3Up");
  TH1D* D1_OTHER_mcStatBin4Up = (TH1D*)file->Get("D1_OTHER_mcStatBin4Up");
  TH1D* D1_OTHER_mcStatBin5Up = (TH1D*)file->Get("D1_OTHER_mcStatBin5Up");
  TH1D* D1_OTHER_mcStatBin6Up = (TH1D*)file->Get("D1_OTHER_mcStatBin6Up");
  TH1D* D1_OTHER_mcStatBin7Up = (TH1D*)file->Get("D1_OTHER_mcStatBin7Up");
  TH1D* D1_OTHER_mcStatBin8Up = (TH1D*)file->Get("D1_OTHER_mcStatBin8Up");
  D1_OTHER_mcStatBin1Up->SetName(("D1_OTHER_mcStatD1OTHERBin1_"+year+"Up").c_str());
  D1_OTHER_mcStatBin2Up->SetName(("D1_OTHER_mcStatD1OTHERBin2_"+year+"Up").c_str());
  D1_OTHER_mcStatBin3Up->SetName(("D1_OTHER_mcStatD1OTHERBin3_"+year+"Up").c_str());
  D1_OTHER_mcStatBin4Up->SetName(("D1_OTHER_mcStatD1OTHERBin4_"+year+"Up").c_str());
  D1_OTHER_mcStatBin5Up->SetName(("D1_OTHER_mcStatD1OTHERBin5_"+year+"Up").c_str());
  D1_OTHER_mcStatBin6Up->SetName(("D1_OTHER_mcStatD1OTHERBin6_"+year+"Up").c_str());
  D1_OTHER_mcStatBin7Up->SetName(("D1_OTHER_mcStatD1OTHERBin7_"+year+"Up").c_str());
  D1_OTHER_mcStatBin8Up->SetName(("D1_OTHER_mcStatD1OTHERBin8_"+year+"Up").c_str());
  D1_OTHER_mcStatBin1Up->Write();
  D1_OTHER_mcStatBin2Up->Write();
  D1_OTHER_mcStatBin3Up->Write();
  D1_OTHER_mcStatBin4Up->Write();
  D1_OTHER_mcStatBin5Up->Write();
  D1_OTHER_mcStatBin6Up->Write();
  D1_OTHER_mcStatBin7Up->Write();
  D1_OTHER_mcStatBin8Up->Write();
  TH1D* D1_OTHER_mcStatBin1Down = (TH1D*)file->Get("D1_OTHER_mcStatBin1Down");
  TH1D* D1_OTHER_mcStatBin2Down = (TH1D*)file->Get("D1_OTHER_mcStatBin2Down");
  TH1D* D1_OTHER_mcStatBin3Down = (TH1D*)file->Get("D1_OTHER_mcStatBin3Down");
  TH1D* D1_OTHER_mcStatBin4Down = (TH1D*)file->Get("D1_OTHER_mcStatBin4Down");
  TH1D* D1_OTHER_mcStatBin5Down = (TH1D*)file->Get("D1_OTHER_mcStatBin5Down");
  TH1D* D1_OTHER_mcStatBin6Down = (TH1D*)file->Get("D1_OTHER_mcStatBin6Down");
  TH1D* D1_OTHER_mcStatBin7Down = (TH1D*)file->Get("D1_OTHER_mcStatBin7Down");
  TH1D* D1_OTHER_mcStatBin8Down = (TH1D*)file->Get("D1_OTHER_mcStatBin8Down");
  D1_OTHER_mcStatBin1Down->SetName(("D1_OTHER_mcStatD1OTHERBin1_"+year+"Down").c_str());
  D1_OTHER_mcStatBin2Down->SetName(("D1_OTHER_mcStatD1OTHERBin2_"+year+"Down").c_str());
  D1_OTHER_mcStatBin3Down->SetName(("D1_OTHER_mcStatD1OTHERBin3_"+year+"Down").c_str());
  D1_OTHER_mcStatBin4Down->SetName(("D1_OTHER_mcStatD1OTHERBin4_"+year+"Down").c_str());
  D1_OTHER_mcStatBin5Down->SetName(("D1_OTHER_mcStatD1OTHERBin5_"+year+"Down").c_str());
  D1_OTHER_mcStatBin6Down->SetName(("D1_OTHER_mcStatD1OTHERBin6_"+year+"Down").c_str());
  D1_OTHER_mcStatBin7Down->SetName(("D1_OTHER_mcStatD1OTHERBin7_"+year+"Down").c_str());
  D1_OTHER_mcStatBin8Down->SetName(("D1_OTHER_mcStatD1OTHERBin8_"+year+"Down").c_str());
  D1_OTHER_mcStatBin1Down->Write();
  D1_OTHER_mcStatBin2Down->Write();
  D1_OTHER_mcStatBin3Down->Write();
  D1_OTHER_mcStatBin4Down->Write();
  D1_OTHER_mcStatBin5Down->Write();
  D1_OTHER_mcStatBin6Down->Write();
  D1_OTHER_mcStatBin7Down->Write();
  D1_OTHER_mcStatBin8Down->Write();

  TH1D* D2_OTHER_mcStatBin1Up = (TH1D*)file->Get("D2_OTHER_mcStatBin1Up");
  TH1D* D2_OTHER_mcStatBin2Up = (TH1D*)file->Get("D2_OTHER_mcStatBin2Up");
  TH1D* D2_OTHER_mcStatBin3Up = (TH1D*)file->Get("D2_OTHER_mcStatBin3Up");
  TH1D* D2_OTHER_mcStatBin4Up = (TH1D*)file->Get("D2_OTHER_mcStatBin4Up");
  TH1D* D2_OTHER_mcStatBin5Up = (TH1D*)file->Get("D2_OTHER_mcStatBin5Up");
  TH1D* D2_OTHER_mcStatBin6Up = (TH1D*)file->Get("D2_OTHER_mcStatBin6Up");
  TH1D* D2_OTHER_mcStatBin7Up = (TH1D*)file->Get("D2_OTHER_mcStatBin7Up");
  TH1D* D2_OTHER_mcStatBin8Up = (TH1D*)file->Get("D2_OTHER_mcStatBin8Up");
  D2_OTHER_mcStatBin1Up->SetName(("D2_OTHER_mcStatD2OTHERBin1_"+year+"Up").c_str());
  D2_OTHER_mcStatBin2Up->SetName(("D2_OTHER_mcStatD2OTHERBin2_"+year+"Up").c_str());
  D2_OTHER_mcStatBin3Up->SetName(("D2_OTHER_mcStatD2OTHERBin3_"+year+"Up").c_str());
  D2_OTHER_mcStatBin4Up->SetName(("D2_OTHER_mcStatD2OTHERBin4_"+year+"Up").c_str());
  D2_OTHER_mcStatBin5Up->SetName(("D2_OTHER_mcStatD2OTHERBin5_"+year+"Up").c_str());
  D2_OTHER_mcStatBin6Up->SetName(("D2_OTHER_mcStatD2OTHERBin6_"+year+"Up").c_str());
  D2_OTHER_mcStatBin7Up->SetName(("D2_OTHER_mcStatD2OTHERBin7_"+year+"Up").c_str());
  D2_OTHER_mcStatBin8Up->SetName(("D2_OTHER_mcStatD2OTHERBin8_"+year+"Up").c_str());
  D2_OTHER_mcStatBin1Up->Write();
  D2_OTHER_mcStatBin2Up->Write();
  D2_OTHER_mcStatBin3Up->Write();
  D2_OTHER_mcStatBin4Up->Write();
  D2_OTHER_mcStatBin5Up->Write();
  D2_OTHER_mcStatBin6Up->Write();
  D2_OTHER_mcStatBin7Up->Write();
  D2_OTHER_mcStatBin8Up->Write();
  TH1D* D2_OTHER_mcStatBin1Down = (TH1D*)file->Get("D2_OTHER_mcStatBin1Down");
  TH1D* D2_OTHER_mcStatBin2Down = (TH1D*)file->Get("D2_OTHER_mcStatBin2Down");
  TH1D* D2_OTHER_mcStatBin3Down = (TH1D*)file->Get("D2_OTHER_mcStatBin3Down");
  TH1D* D2_OTHER_mcStatBin4Down = (TH1D*)file->Get("D2_OTHER_mcStatBin4Down");
  TH1D* D2_OTHER_mcStatBin5Down = (TH1D*)file->Get("D2_OTHER_mcStatBin5Down");
  TH1D* D2_OTHER_mcStatBin6Down = (TH1D*)file->Get("D2_OTHER_mcStatBin6Down");
  TH1D* D2_OTHER_mcStatBin7Down = (TH1D*)file->Get("D2_OTHER_mcStatBin7Down");
  TH1D* D2_OTHER_mcStatBin8Down = (TH1D*)file->Get("D2_OTHER_mcStatBin8Down");
  D2_OTHER_mcStatBin1Down->SetName(("D2_OTHER_mcStatD2OTHERBin1_"+year+"Down").c_str());
  D2_OTHER_mcStatBin2Down->SetName(("D2_OTHER_mcStatD2OTHERBin2_"+year+"Down").c_str());
  D2_OTHER_mcStatBin3Down->SetName(("D2_OTHER_mcStatD2OTHERBin3_"+year+"Down").c_str());
  D2_OTHER_mcStatBin4Down->SetName(("D2_OTHER_mcStatD2OTHERBin4_"+year+"Down").c_str());
  D2_OTHER_mcStatBin5Down->SetName(("D2_OTHER_mcStatD2OTHERBin5_"+year+"Down").c_str());
  D2_OTHER_mcStatBin6Down->SetName(("D2_OTHER_mcStatD2OTHERBin6_"+year+"Down").c_str());
  D2_OTHER_mcStatBin7Down->SetName(("D2_OTHER_mcStatD2OTHERBin7_"+year+"Down").c_str());
  D2_OTHER_mcStatBin8Down->SetName(("D2_OTHER_mcStatD2OTHERBin8_"+year+"Down").c_str());
  D2_OTHER_mcStatBin1Down->Write();
  D2_OTHER_mcStatBin2Down->Write();
  D2_OTHER_mcStatBin3Down->Write();
  D2_OTHER_mcStatBin4Down->Write();
  D2_OTHER_mcStatBin5Down->Write();
  D2_OTHER_mcStatBin6Down->Write();
  D2_OTHER_mcStatBin7Down->Write();
  D2_OTHER_mcStatBin8Down->Write();

  TH1D* D3_OTHER_mcStatBin1Up = (TH1D*)file->Get("D3_OTHER_mcStatBin1Up");
  TH1D* D3_OTHER_mcStatBin2Up = (TH1D*)file->Get("D3_OTHER_mcStatBin2Up");
  TH1D* D3_OTHER_mcStatBin3Up = (TH1D*)file->Get("D3_OTHER_mcStatBin3Up");
  TH1D* D3_OTHER_mcStatBin4Up = (TH1D*)file->Get("D3_OTHER_mcStatBin4Up");
  TH1D* D3_OTHER_mcStatBin5Up = (TH1D*)file->Get("D3_OTHER_mcStatBin5Up");
  TH1D* D3_OTHER_mcStatBin6Up = (TH1D*)file->Get("D3_OTHER_mcStatBin6Up");
  TH1D* D3_OTHER_mcStatBin7Up = (TH1D*)file->Get("D3_OTHER_mcStatBin7Up");
  TH1D* D3_OTHER_mcStatBin8Up = (TH1D*)file->Get("D3_OTHER_mcStatBin8Up");
  D3_OTHER_mcStatBin1Up->SetName(("D3_OTHER_mcStatD3OTHERBin1_"+year+"Up").c_str());
  D3_OTHER_mcStatBin2Up->SetName(("D3_OTHER_mcStatD3OTHERBin2_"+year+"Up").c_str());
  D3_OTHER_mcStatBin3Up->SetName(("D3_OTHER_mcStatD3OTHERBin3_"+year+"Up").c_str());
  D3_OTHER_mcStatBin4Up->SetName(("D3_OTHER_mcStatD3OTHERBin4_"+year+"Up").c_str());
  D3_OTHER_mcStatBin5Up->SetName(("D3_OTHER_mcStatD3OTHERBin5_"+year+"Up").c_str());
  D3_OTHER_mcStatBin6Up->SetName(("D3_OTHER_mcStatD3OTHERBin6_"+year+"Up").c_str());
  D3_OTHER_mcStatBin7Up->SetName(("D3_OTHER_mcStatD3OTHERBin7_"+year+"Up").c_str());
  D3_OTHER_mcStatBin8Up->SetName(("D3_OTHER_mcStatD3OTHERBin8_"+year+"Up").c_str());
  D3_OTHER_mcStatBin1Up->Write();
  D3_OTHER_mcStatBin2Up->Write();
  D3_OTHER_mcStatBin3Up->Write();
  D3_OTHER_mcStatBin4Up->Write();
  D3_OTHER_mcStatBin5Up->Write();
  D3_OTHER_mcStatBin6Up->Write();
  D3_OTHER_mcStatBin7Up->Write();
  D3_OTHER_mcStatBin8Up->Write();
  TH1D* D3_OTHER_mcStatBin1Down = (TH1D*)file->Get("D3_OTHER_mcStatBin1Down");
  TH1D* D3_OTHER_mcStatBin2Down = (TH1D*)file->Get("D3_OTHER_mcStatBin2Down");
  TH1D* D3_OTHER_mcStatBin3Down = (TH1D*)file->Get("D3_OTHER_mcStatBin3Down");
  TH1D* D3_OTHER_mcStatBin4Down = (TH1D*)file->Get("D3_OTHER_mcStatBin4Down");
  TH1D* D3_OTHER_mcStatBin5Down = (TH1D*)file->Get("D3_OTHER_mcStatBin5Down");
  TH1D* D3_OTHER_mcStatBin6Down = (TH1D*)file->Get("D3_OTHER_mcStatBin6Down");
  TH1D* D3_OTHER_mcStatBin7Down = (TH1D*)file->Get("D3_OTHER_mcStatBin7Down");
  TH1D* D3_OTHER_mcStatBin8Down = (TH1D*)file->Get("D3_OTHER_mcStatBin8Down");
  D3_OTHER_mcStatBin1Down->SetName(("D3_OTHER_mcStatD3OTHERBin1_"+year+"Down").c_str());
  D3_OTHER_mcStatBin2Down->SetName(("D3_OTHER_mcStatD3OTHERBin2_"+year+"Down").c_str());
  D3_OTHER_mcStatBin3Down->SetName(("D3_OTHER_mcStatD3OTHERBin3_"+year+"Down").c_str());
  D3_OTHER_mcStatBin4Down->SetName(("D3_OTHER_mcStatD3OTHERBin4_"+year+"Down").c_str());
  D3_OTHER_mcStatBin5Down->SetName(("D3_OTHER_mcStatD3OTHERBin5_"+year+"Down").c_str());
  D3_OTHER_mcStatBin6Down->SetName(("D3_OTHER_mcStatD3OTHERBin6_"+year+"Down").c_str());
  D3_OTHER_mcStatBin7Down->SetName(("D3_OTHER_mcStatD3OTHERBin7_"+year+"Down").c_str());
  D3_OTHER_mcStatBin8Down->SetName(("D3_OTHER_mcStatD3OTHERBin8_"+year+"Down").c_str());
  D3_OTHER_mcStatBin1Down->Write();
  D3_OTHER_mcStatBin2Down->Write();
  D3_OTHER_mcStatBin3Down->Write();
  D3_OTHER_mcStatBin4Down->Write();
  D3_OTHER_mcStatBin5Down->Write();
  D3_OTHER_mcStatBin6Down->Write();
  D3_OTHER_mcStatBin7Down->Write();
  D3_OTHER_mcStatBin8Down->Write();

  TH1D* D4_OTHER_mcStatBin1Up = (TH1D*)file->Get("D4_OTHER_mcStatBin1Up");
  TH1D* D4_OTHER_mcStatBin2Up = (TH1D*)file->Get("D4_OTHER_mcStatBin2Up");
  TH1D* D4_OTHER_mcStatBin3Up = (TH1D*)file->Get("D4_OTHER_mcStatBin3Up");
  TH1D* D4_OTHER_mcStatBin4Up = (TH1D*)file->Get("D4_OTHER_mcStatBin4Up");
  TH1D* D4_OTHER_mcStatBin5Up = (TH1D*)file->Get("D4_OTHER_mcStatBin5Up");
  TH1D* D4_OTHER_mcStatBin6Up = (TH1D*)file->Get("D4_OTHER_mcStatBin6Up");
  TH1D* D4_OTHER_mcStatBin7Up = (TH1D*)file->Get("D4_OTHER_mcStatBin7Up");
  TH1D* D4_OTHER_mcStatBin8Up = (TH1D*)file->Get("D4_OTHER_mcStatBin8Up");
  D4_OTHER_mcStatBin1Up->SetName(("D4_OTHER_mcStatD4OTHERBin1_"+year+"Up").c_str());
  D4_OTHER_mcStatBin2Up->SetName(("D4_OTHER_mcStatD4OTHERBin2_"+year+"Up").c_str());
  D4_OTHER_mcStatBin3Up->SetName(("D4_OTHER_mcStatD4OTHERBin3_"+year+"Up").c_str());
  D4_OTHER_mcStatBin4Up->SetName(("D4_OTHER_mcStatD4OTHERBin4_"+year+"Up").c_str());
  D4_OTHER_mcStatBin5Up->SetName(("D4_OTHER_mcStatD4OTHERBin5_"+year+"Up").c_str());
  D4_OTHER_mcStatBin6Up->SetName(("D4_OTHER_mcStatD4OTHERBin6_"+year+"Up").c_str());
  D4_OTHER_mcStatBin7Up->SetName(("D4_OTHER_mcStatD4OTHERBin7_"+year+"Up").c_str());
  D4_OTHER_mcStatBin8Up->SetName(("D4_OTHER_mcStatD4OTHERBin8_"+year+"Up").c_str());
  D4_OTHER_mcStatBin1Up->Write();
  D4_OTHER_mcStatBin2Up->Write();
  D4_OTHER_mcStatBin3Up->Write();
  D4_OTHER_mcStatBin4Up->Write();
  D4_OTHER_mcStatBin5Up->Write();
  D4_OTHER_mcStatBin6Up->Write();
  D4_OTHER_mcStatBin7Up->Write();
  D4_OTHER_mcStatBin8Up->Write();
  TH1D* D4_OTHER_mcStatBin1Down = (TH1D*)file->Get("D4_OTHER_mcStatBin1Down");
  TH1D* D4_OTHER_mcStatBin2Down = (TH1D*)file->Get("D4_OTHER_mcStatBin2Down");
  TH1D* D4_OTHER_mcStatBin3Down = (TH1D*)file->Get("D4_OTHER_mcStatBin3Down");
  TH1D* D4_OTHER_mcStatBin4Down = (TH1D*)file->Get("D4_OTHER_mcStatBin4Down");
  TH1D* D4_OTHER_mcStatBin5Down = (TH1D*)file->Get("D4_OTHER_mcStatBin5Down");
  TH1D* D4_OTHER_mcStatBin6Down = (TH1D*)file->Get("D4_OTHER_mcStatBin6Down");
  TH1D* D4_OTHER_mcStatBin7Down = (TH1D*)file->Get("D4_OTHER_mcStatBin7Down");
  TH1D* D4_OTHER_mcStatBin8Down = (TH1D*)file->Get("D4_OTHER_mcStatBin8Down");
  D4_OTHER_mcStatBin1Down->SetName(("D4_OTHER_mcStatD4OTHERBin1_"+year+"Down").c_str());
  D4_OTHER_mcStatBin2Down->SetName(("D4_OTHER_mcStatD4OTHERBin2_"+year+"Down").c_str());
  D4_OTHER_mcStatBin3Down->SetName(("D4_OTHER_mcStatD4OTHERBin3_"+year+"Down").c_str());
  D4_OTHER_mcStatBin4Down->SetName(("D4_OTHER_mcStatD4OTHERBin4_"+year+"Down").c_str());
  D4_OTHER_mcStatBin5Down->SetName(("D4_OTHER_mcStatD4OTHERBin5_"+year+"Down").c_str());
  D4_OTHER_mcStatBin6Down->SetName(("D4_OTHER_mcStatD4OTHERBin6_"+year+"Down").c_str());
  D4_OTHER_mcStatBin7Down->SetName(("D4_OTHER_mcStatD4OTHERBin7_"+year+"Down").c_str());
  D4_OTHER_mcStatBin8Down->SetName(("D4_OTHER_mcStatD4OTHERBin8_"+year+"Down").c_str());
  D4_OTHER_mcStatBin1Down->Write();
  D4_OTHER_mcStatBin2Down->Write();
  D4_OTHER_mcStatBin3Down->Write();
  D4_OTHER_mcStatBin4Down->Write();
  D4_OTHER_mcStatBin5Down->Write();
  D4_OTHER_mcStatBin6Down->Write();
  D4_OTHER_mcStatBin7Down->Write();
  D4_OTHER_mcStatBin8Down->Write();


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
  // otherMC_hist_D4.createHistogram("nj")->Draw("H");

  // TCanvas *c4 = new TCanvas("c4","c4");
  // sigMC_hist_D4.createHistogram("nj")->Draw("H");


}
