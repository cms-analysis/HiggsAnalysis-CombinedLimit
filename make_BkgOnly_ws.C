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
// Dec 20, 2018 : This branch does is intended for "background only" fits to the tt bar background (instead of pseudodata)

void construct_formula(string procName, RooArgList& binlist, const RooArgList& paramlist) {

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



  int max_bin = 18; // 14 means just njets=14, 20 means last bin is inclusive up through njets=20

  for (int i=0; i<8; i++) {

    stringstream form;
    RooArgList formArgList;

    form << "(@0";
    formArgList.add(paramlist[0]); // N7_tt for this MVA bin

    if (i>=1) { // for bin 1 and up
      for (int j=0; j<i; j++) {
	//form << "*(@3 + TMath::Power( TMath::Power( @2-@3 , "<<j<<"-0 ) / TMath::Power( @1-@3 , "<<j<<"-2 ) , 1/(2-0) ))";
	form << "*(@3 + TMath::Power( TMath::Power( @2-@3 , "<<j<<" ) / TMath::Power( @1-@3 , "<<j-2<<" ) , 1/2 ))";
      }
      formArgList.add(paramlist[1]); // a0_tt
      formArgList.add(paramlist[2]); // a1_tt
      formArgList.add(paramlist[3]); // a2_tt
    } // end bin 1 and up

    // The last bin covers from njet=14 through njet=max_bin
    if (i==7) {
      for (int k=8; k<=max_bin-7; k++) {
	form << " + 1";
	for (int j=0; j<k; j++) {
	  //form << "*(@3 + TMath::Power( TMath::Power( @2-@3 , "<<j<<"-0 ) / TMath::Power( @1-@3 , "<<j<<"-2 ) , 1/(2-0) ))";
	  form << "*(@3 + TMath::Power( TMath::Power( @2-@3 , "<<j<<" ) / TMath::Power( @1-@3 , "<<j-2<<" ) , 1/2 ))";
	}
      }
    }
    form << ")";

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


void make_BkgOnly_ws(string syst = "") {
  // syst should be left empty for no systematic:
  //   .x make_BkgOnly_ws.C();
  // Or set to "_btgUp", "_btgDown", "_lepUp", "_lepDown", for example:
  //   .x make_BkgOnly_ws.C("_btgUp");



  using namespace RooFit;
  // Load the combine library to get access to the RooParametricHist
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
    
  // Output file and workspace 
  TFile *fOut = new TFile("MVA_BkgOnly_ws.root","RECREATE");
  RooWorkspace *wspace = new RooWorkspace("wspace","wspace");

  // // njet is our variable, 8 bins, from 7 up through 14,
  // //   Note that njet=14 is inclusive as >=14
  // // D1, D2, D3, D4 are the MVA bins
  wspace->factory("CMS_th1x[0,8]");
  wspace->var("CMS_th1x")->setBins(8);
  RooArgSet vars_D1(*wspace->var("CMS_th1x"));
  RooArgSet vars_D2(*wspace->var("CMS_th1x"));
  RooArgSet vars_D3(*wspace->var("CMS_th1x"));
  RooArgSet vars_D4(*wspace->var("CMS_th1x"));

    
  // file for obtaining histograms
  TFile* file = TFile::Open("Keras_V1.2.4/njets_rebin_for_Aron.root");

  string tt_file = "---------------------------------------------";
  
  tt_file = "D1_TT_h_njets_pt30_1l"+syst;
  TH1D* data_th1_D1 = (TH1D*)file->Get(tt_file.c_str()); // treat tt MC as the data
  TH1D* sigMC_th1_D1 = (TH1D*)file->Get("D1_RPV_550_h_njets_pt30_1l");

  tt_file = "D2_TT_h_njets_pt30_1l"+syst;
  TH1D* data_th1_D2 = (TH1D*)file->Get(tt_file.c_str()); // treat tt MC as the data
  TH1D* sigMC_th1_D2 = (TH1D*)file->Get("D2_RPV_550_h_njets_pt30_1l");

  tt_file = "D3_TT_h_njets_pt30_1l"+syst;
  TH1D* data_th1_D3 = (TH1D*)file->Get(tt_file.c_str()); // treat tt MC as the data
  TH1D* sigMC_th1_D3 = (TH1D*)file->Get("D3_RPV_550_h_njets_pt30_1l");

  tt_file = "D4_TT_h_njets_pt30_1l"+syst;
  TH1D* data_th1_D4 = (TH1D*)file->Get(tt_file.c_str()); // treat tt MC as the data
  TH1D* sigMC_th1_D4 = (TH1D*)file->Get("D4_RPV_550_h_njets_pt30_1l");


  // tt bkg param setup
  RooRealVar a0_tt("a0_tt","a0 of tt bkg shape",0.28,0.0,1.0);
  RooRealVar a1_tt("a1_tt","a1 of tt bkg shape",0.24,0.0,1.0);
  RooRealVar a2_tt("a2_tt","a2 of tt bkg shape",0.10,-0.5,0.5);

  // RooRealVar a0_tt_D1("a0_tt_D1","a0 of tt bkg shape D1",0.28,0.0,1.0);
  // RooRealVar a1_tt_D1("a1_tt_D1","a1 of tt bkg shape D1",0.24,0.0,1.0);
  // RooRealVar a2_tt_D1("a2_tt_D1","a2 of tt bkg shape D1",0.10,-0.5,0.5);

  // RooRealVar a0_tt_D2("a0_tt_D2","a0 of tt bkg shape D2",0.28,0.0,1.0);
  // RooRealVar a1_tt_D2("a1_tt_D2","a1 of tt bkg shape D2",0.24,0.0,1.0);
  // RooRealVar a2_tt_D2("a2_tt_D2","a2 of tt bkg shape D2",0.10,-0.5,0.5);

  // RooRealVar a0_tt_D3("a0_tt_D3","a0 of tt bkg shape D3",0.28,0.0,1.0);
  // RooRealVar a1_tt_D3("a1_tt_D3","a1 of tt bkg shape D3",0.24,0.0,1.0);
  // RooRealVar a2_tt_D3("a2_tt_D3","a2 of tt bkg shape D3",0.10,-0.5,0.5);

  // RooRealVar a0_tt_D4("a0_tt_D4","a0 of tt bkg shape D4",0.28,0.0,1.0);
  // RooRealVar a1_tt_D4("a1_tt_D4","a1 of tt bkg shape D4",0.24,0.0,1.0);
  // RooRealVar a2_tt_D4("a2_tt_D4","a2 of tt bkg shape D4",0.10,-0.5,0.5);


  double_t n7_tt_portion_D1 = data_th1_D1->GetBinContent(1) - sigMC_th1_D1->GetBinContent(1);
  RooRealVar N7_tt_D1("N7_tt_D1","njets 7 for tt bkg in MVA D1",n7_tt_portion_D1,n7_tt_portion_D1-1500,n7_tt_portion_D1+1500);

  double_t n7_tt_portion_D2 = data_th1_D2->GetBinContent(1) - sigMC_th1_D2->GetBinContent(1);
  RooRealVar N7_tt_D2("N7_tt_D2","njets 7 for tt bkg in MVA D2",n7_tt_portion_D2,n7_tt_portion_D2-1500,n7_tt_portion_D2+1500);

  double_t n7_tt_portion_D3 = data_th1_D3->GetBinContent(1) - sigMC_th1_D3->GetBinContent(1);
  RooRealVar N7_tt_D3("N7_tt_D3","njets 7 for tt bkg in MVA D3",n7_tt_portion_D3,n7_tt_portion_D3-1500,n7_tt_portion_D3+1500);

  double_t n7_tt_portion_D4 = data_th1_D4->GetBinContent(1) - sigMC_th1_D4->GetBinContent(1);
  RooRealVar N7_tt_D4("N7_tt_D4","njets 7 for tt bkg in MVA D4",n7_tt_portion_D4,n7_tt_portion_D4-1500,n7_tt_portion_D4+1500);

  //RooRealVar N7_tt_D1("N7_tt_D1","njets 7 for tt bkg in MVA D1",1.0,1.0,1.0);
  //RooRealVar N7_tt_D2("N7_tt_D2","njets 7 for tt bkg in MVA D2",1.0,1.0,1.0);
  //RooRealVar N7_tt_D3("N7_tt_D3","njets 7 for tt bkg in MVA D3",1.0,1.0,1.0);
  //RooRealVar N7_tt_D4("N7_tt_D4","njets 7 for tt bkg in MVA D4",1.0,1.0,1.0);



  // ----------------------  MVA bin 1  ------------------

  // Dataset with 8 bins
  RooDataHist data_hist_D1("data_obs_D1","Data observed in MVA bin 1",vars_D1,data_th1_D1);
  wspace->import(data_hist_D1);

  // ttbar bkg in D1
  RooArgList parlist_D1(N7_tt_D1,a0_tt,a1_tt,a2_tt);  // list of shape parameters for tt bkg
  //RooArgList parlist_D1(N7_tt_D1,a0_tt_D1,a1_tt_D1,a2_tt_D1);  // list of shape parameters for tt bkg
  RooArgList *bkg_tt_bins_D1 = new RooArgList();
  string procName_D1 = "background_tt_D1";
  construct_formula(procName_D1,*bkg_tt_bins_D1,parlist_D1);
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
  RooArgList parlist_D2(N7_tt_D2,a0_tt,a1_tt,a2_tt);  // list of shape parameters for tt bkg
  //RooArgList parlist_D2(N7_tt_D2,a0_tt_D2,a1_tt_D2,a2_tt_D2);  // list of shape parameters for tt bkg
  RooArgList *bkg_tt_bins_D2 = new RooArgList();
  string procName_D2 = "background_tt_D2";
  construct_formula(procName_D2,*bkg_tt_bins_D2,parlist_D2);
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
  RooArgList parlist_D3(N7_tt_D3,a0_tt,a1_tt,a2_tt);  // list of shape parameters for tt bkg
  //RooArgList parlist_D3(N7_tt_D3,a0_tt_D3,a1_tt_D3,a2_tt_D3);  // list of shape parameters for tt bkg
  RooArgList *bkg_tt_bins_D3 = new RooArgList();
  string procName_D3 = "background_tt_D3";
  construct_formula(procName_D3,*bkg_tt_bins_D3,parlist_D3);
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
  RooArgList parlist_D4(N7_tt_D4,a0_tt,a1_tt,a2_tt);  // list of shape parameters for tt bkg
  //RooArgList parlist_D4(N7_tt_D4,a0_tt_D4,a1_tt_D4,a2_tt_D4);  // list of shape parameters for tt bkg
  RooArgList *bkg_tt_bins_D4 = new RooArgList();
  string procName_D4 = "background_tt_D4";
  construct_formula(procName_D4,*bkg_tt_bins_D4,parlist_D4);
  RooParametricHist background_tt_D4(procName_D4.c_str(),"",*wspace->var("CMS_th1x"),*bkg_tt_bins_D4,*data_th1_D4);
  wspace->import(background_tt_D4,RooFit::RecycleConflictNodes());
  stringstream procNameD4Norm;
  procNameD4Norm << procName_D4 << "_norm";
  RooAddition tt_norm_D4(procNameD4Norm.str().c_str(),"",*bkg_tt_bins_D4);
  wspace->import(tt_norm_D4,RooFit::RecycleConflictNodes());

  // =================================================================================

   
  fOut->cd();

  sigMC_th1_D1->SetName("sigMC_th1_D1");
  sigMC_th1_D2->SetName("sigMC_th1_D2");
  sigMC_th1_D3->SetName("sigMC_th1_D3");
  sigMC_th1_D4->SetName("sigMC_th1_D4");

  sigMC_th1_D1->Write();
  sigMC_th1_D2->Write();
  sigMC_th1_D3->Write();
  sigMC_th1_D4->Write();


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
