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

    if (i==0) {
      form << "*TMath::Power(" << h_syst[0]->GetBinContent(i+1) << ",@1)";
      form << "*TMath::Power(" << h_syst[1]->GetBinContent(i+1) << ",@2)";
      form << "*TMath::Power(" << h_syst[2]->GetBinContent(i+1) << ",@3)";
      form << "*TMath::Power(" << h_syst[3]->GetBinContent(i+1) << ",@4)";
    } else if (i>=1) {
      form << "*TMath::Power(" << h_syst[0]->GetBinContent(i+1) << ",@4)";
      form << "*TMath::Power(" << h_syst[1]->GetBinContent(i+1) << ",@5)";
      form << "*TMath::Power(" << h_syst[2]->GetBinContent(i+1) << ",@6)";
      form << "*TMath::Power(" << h_syst[3]->GetBinContent(i+1) << ",@7)";
    }
    // nuisance parameters
    formArgList.add(NPs[0]);
    formArgList.add(NPs[1]);
    formArgList.add(NPs[2]);
    formArgList.add(NPs[3]);

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


void make_MVA_8bin_ws(const string infile_path = "Keras_V1.2.4", const string model = "RPV", const string mass = "550", const string syst = "") {
  using namespace RooFit;
  // Load the combine library to get access to the RooParametricHist
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
    
  // Output file and workspace 
  TFile *fOut = new TFile(("MVA_"+model+"_"+mass+"_ws.root").c_str(),"RECREATE");
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
  TFile* file = TFile::Open((infile_path+"/njets_rebin_for_Aron.root").c_str());

  TH1D* data_th1_D1 = (TH1D*)file->Get(("D1_pseudodataS_"+model+"_"+mass+"_h_njets_pt30_1l"+syst).c_str());
  //TH1D* data_th1_D1 = (TH1D*)file->Get(("D1_pseudodata_h_njets_pt30_1l"+syst).c_str());
  TH1D* otherMC_th1_D1 = (TH1D*)file->Get(("D1_OTHER_h_njets_pt30_1l"+syst).c_str());
  TH1D* sigMC_th1_D1 = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_h_njets_pt30_1l"+syst).c_str());

  TH1D* data_th1_D2 = (TH1D*)file->Get(("D2_pseudodataS_"+model+"_"+mass+"_h_njets_pt30_1l"+syst).c_str());
  //TH1D* data_th1_D2 = (TH1D*)file->Get(("D2_pseudodata_h_njets_pt30_1l"+syst).c_str());
  TH1D* otherMC_th1_D2 = (TH1D*)file->Get(("D2_OTHER_h_njets_pt30_1l"+syst).c_str());
  TH1D* sigMC_th1_D2 = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_h_njets_pt30_1l"+syst).c_str());

  TH1D* data_th1_D3 = (TH1D*)file->Get(("D3_pseudodataS_"+model+"_"+mass+"_h_njets_pt30_1l"+syst).c_str());
  //TH1D* data_th1_D3 = (TH1D*)file->Get(("D3_pseudodata_h_njets_pt30_1l"+syst).c_str());
  TH1D* otherMC_th1_D3 = (TH1D*)file->Get(("D3_OTHER_h_njets_pt30_1l"+syst).c_str());
  TH1D* sigMC_th1_D3 = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_h_njets_pt30_1l"+syst).c_str());

  TH1D* data_th1_D4 = (TH1D*)file->Get(("D4_pseudodataS_"+model+"_"+mass+"_h_njets_pt30_1l"+syst).c_str());
  //TH1D* data_th1_D4 = (TH1D*)file->Get(("D4_pseudodata_h_njets_pt30_1l"+syst).c_str());
  TH1D* otherMC_th1_D4 = (TH1D*)file->Get(("D4_OTHER_h_njets_pt30_1l"+syst).c_str());
  TH1D* sigMC_th1_D4 = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_h_njets_pt30_1l"+syst).c_str());


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


  double_t n7_tt_portion_D1 = data_th1_D1->GetBinContent(1) - otherMC_th1_D1->GetBinContent(1) - sigMC_th1_D1->GetBinContent(1);
  RooRealVar N7_tt_D1("N7_tt_D1","njets 7 for tt bkg in MVA D1",n7_tt_portion_D1,n7_tt_portion_D1-1500,n7_tt_portion_D1+1500);

  double_t n7_tt_portion_D2 = data_th1_D2->GetBinContent(1) - otherMC_th1_D2->GetBinContent(1) - sigMC_th1_D2->GetBinContent(1);
  RooRealVar N7_tt_D2("N7_tt_D2","njets 7 for tt bkg in MVA D2",n7_tt_portion_D2,n7_tt_portion_D2-1500,n7_tt_portion_D2+1500);

  double_t n7_tt_portion_D3 = data_th1_D3->GetBinContent(1) - otherMC_th1_D3->GetBinContent(1) - sigMC_th1_D3->GetBinContent(1);
  RooRealVar N7_tt_D3("N7_tt_D3","njets 7 for tt bkg in MVA D3",n7_tt_portion_D3,n7_tt_portion_D3-1500,n7_tt_portion_D3+1500);

  double_t n7_tt_portion_D4 = data_th1_D4->GetBinContent(1) - otherMC_th1_D4->GetBinContent(1) - sigMC_th1_D4->GetBinContent(1);
  RooRealVar N7_tt_D4("N7_tt_D4","njets 7 for tt bkg in MVA D4",n7_tt_portion_D4,n7_tt_portion_D4-1500,n7_tt_portion_D4+1500);

  //RooRealVar N7_tt_D1("N7_tt_D1","njets 7 for tt bkg in MVA D1",1.0,1.0,1.0);
  //RooRealVar N7_tt_D2("N7_tt_D2","njets 7 for tt bkg in MVA D2",1.0,1.0,1.0);
  //RooRealVar N7_tt_D3("N7_tt_D3","njets 7 for tt bkg in MVA D3",1.0,1.0,1.0);
  //RooRealVar N7_tt_D4("N7_tt_D4","njets 7 for tt bkg in MVA D4",1.0,1.0,1.0);



  // tt shape systematic nuisance parameters
  //wspace->factory("np_tt_JECup[0.0]");
  //wspace->factory("np_tt_JECdown[0.0]");
  //wspace->factory("np_tt_JERup[0.0]");
  //wspace->factory("np_tt_JERdown[0.0]");
  wspace->factory("np_tt_JECup[0.0,0,0]");
  wspace->factory("np_tt_JECdown[0.0,0,0]");
  wspace->factory("np_tt_JERup[0.0,0,0]");
  wspace->factory("np_tt_JERdown[0.0,0,0]");

  //wspace->factory("np_tt_JECup_D1[0.0]");
  //wspace->factory("np_tt_JECup_D2[0.0]");
  //wspace->factory("np_tt_JECup_D3[0.0]");
  //wspace->factory("np_tt_JECup_D4[0.0]");
  //wspace->factory("np_tt_JECdown_D1[0.0]");
  //wspace->factory("np_tt_JECdown_D2[0.0]");
  //wspace->factory("np_tt_JECdown_D3[0.0]");
  //wspace->factory("np_tt_JECdown_D4[0.0]");
  //wspace->factory("np_tt_JERup_D1[0.0]");
  //wspace->factory("np_tt_JERup_D2[0.0]");
  //wspace->factory("np_tt_JERup_D3[0.0]");
  //wspace->factory("np_tt_JERup_D4[0.0]");
  //wspace->factory("np_tt_JERdown_D1[0.0]");
  //wspace->factory("np_tt_JERdown_D2[0.0]");
  //wspace->factory("np_tt_JERdown_D3[0.0]");
  //wspace->factory("np_tt_JERdown_D4[0.0]");

  // Assemble the tt systematic histograms, containing the ratio
  //   "actual after njet variation / average over the MVA bins after njet variation"
  
  // JEC up

  TH1D* ttMC_syst_JECup_D1 = (TH1D*)file->Get("D1_TT_h_njets_pt30_1l_JECUp");
  TH1D* ttMC_syst_JECup_D2 = (TH1D*)file->Get("D2_TT_h_njets_pt30_1l_JECUp");
  TH1D* ttMC_syst_JECup_D3 = (TH1D*)file->Get("D3_TT_h_njets_pt30_1l_JECUp");
  TH1D* ttMC_syst_JECup_D4 = (TH1D*)file->Get("D4_TT_h_njets_pt30_1l_JECUp");

  TH1D* ttMC_syst_JECup_AVG = (TH1D*)ttMC_syst_JECup_D1->Clone("ttMC_syst_JECup_AVG");
  ttMC_syst_JECup_AVG->Add(ttMC_syst_JECup_D2);
  ttMC_syst_JECup_AVG->Add(ttMC_syst_JECup_D3);
  ttMC_syst_JECup_AVG->Add(ttMC_syst_JECup_D4);
  ttMC_syst_JECup_AVG->Scale(1/ttMC_syst_JECup_AVG->Integral());

  ttMC_syst_JECup_D1->Scale(1/ttMC_syst_JECup_D1->Integral());
  ttMC_syst_JECup_D2->Scale(1/ttMC_syst_JECup_D2->Integral());
  ttMC_syst_JECup_D3->Scale(1/ttMC_syst_JECup_D3->Integral());
  ttMC_syst_JECup_D4->Scale(1/ttMC_syst_JECup_D4->Integral());

  ttMC_syst_JECup_D1->Divide(ttMC_syst_JECup_AVG);
  ttMC_syst_JECup_D2->Divide(ttMC_syst_JECup_AVG);
  ttMC_syst_JECup_D3->Divide(ttMC_syst_JECup_AVG);
  ttMC_syst_JECup_D4->Divide(ttMC_syst_JECup_AVG);

  // TCanvas *cJECup = new TCanvas("cJECup","cJECup",800,800);
  // cJECup->Divide(2,2);
  // cJECup->cd(1);
  // ttMC_syst_JECup_D1->Draw();
  // cJECup->cd(2);
  // ttMC_syst_JECup_D2->Draw();
  // cJECup->cd(3);
  // ttMC_syst_JECup_D3->Draw();
  // cJECup->cd(4);
  // ttMC_syst_JECup_D4->Draw();
  // cJECup->SaveAs("ratios_JECup.png");

  // JEC down

  TH1D* ttMC_syst_JECdown_D1 = (TH1D*)file->Get("D1_TT_h_njets_pt30_1l_JECDown");
  TH1D* ttMC_syst_JECdown_D2 = (TH1D*)file->Get("D2_TT_h_njets_pt30_1l_JECDown");
  TH1D* ttMC_syst_JECdown_D3 = (TH1D*)file->Get("D3_TT_h_njets_pt30_1l_JECDown");
  TH1D* ttMC_syst_JECdown_D4 = (TH1D*)file->Get("D4_TT_h_njets_pt30_1l_JECDown");

  TH1D* ttMC_syst_JECdown_AVG = (TH1D*)ttMC_syst_JECdown_D1->Clone("ttMC_syst_JECdown_AVG");
  ttMC_syst_JECdown_AVG->Add(ttMC_syst_JECdown_D2);
  ttMC_syst_JECdown_AVG->Add(ttMC_syst_JECdown_D3);
  ttMC_syst_JECdown_AVG->Add(ttMC_syst_JECdown_D4);
  ttMC_syst_JECdown_AVG->Scale(1/ttMC_syst_JECdown_AVG->Integral());

  ttMC_syst_JECdown_D1->Scale(1/ttMC_syst_JECdown_D1->Integral());
  ttMC_syst_JECdown_D2->Scale(1/ttMC_syst_JECdown_D2->Integral());
  ttMC_syst_JECdown_D3->Scale(1/ttMC_syst_JECdown_D3->Integral());
  ttMC_syst_JECdown_D4->Scale(1/ttMC_syst_JECdown_D4->Integral());

  ttMC_syst_JECdown_D1->Divide(ttMC_syst_JECdown_AVG);
  ttMC_syst_JECdown_D2->Divide(ttMC_syst_JECdown_AVG);
  ttMC_syst_JECdown_D3->Divide(ttMC_syst_JECdown_AVG);
  ttMC_syst_JECdown_D4->Divide(ttMC_syst_JECdown_AVG);

  // TCanvas *cJECdown = new TCanvas("cJECdown","cJECdown",800,800);
  // cJECdown->Divide(2,2);
  // cJECdown->cd(1);
  // ttMC_syst_JECdown_D1->Draw();
  // cJECdown->cd(2);
  // ttMC_syst_JECdown_D2->Draw();
  // cJECdown->cd(3);
  // ttMC_syst_JECdown_D3->Draw();
  // cJECdown->cd(4);
  // ttMC_syst_JECdown_D4->Draw();
  // cJECdown->SaveAs("ratios_JECdown.png");

  // JER up

  TH1D* ttMC_syst_JERup_D1 = (TH1D*)file->Get("D1_TT_h_njets_pt30_1l_JERUp");
  TH1D* ttMC_syst_JERup_D2 = (TH1D*)file->Get("D2_TT_h_njets_pt30_1l_JERUp");
  TH1D* ttMC_syst_JERup_D3 = (TH1D*)file->Get("D3_TT_h_njets_pt30_1l_JERUp");
  TH1D* ttMC_syst_JERup_D4 = (TH1D*)file->Get("D4_TT_h_njets_pt30_1l_JERUp");

  TH1D* ttMC_syst_JERup_AVG = (TH1D*)ttMC_syst_JERup_D1->Clone("ttMC_syst_JERup_AVG");
  ttMC_syst_JERup_AVG->Add(ttMC_syst_JERup_D2);
  ttMC_syst_JERup_AVG->Add(ttMC_syst_JERup_D3);
  ttMC_syst_JERup_AVG->Add(ttMC_syst_JERup_D4);
  ttMC_syst_JERup_AVG->Scale(1/ttMC_syst_JERup_AVG->Integral());

  ttMC_syst_JERup_D1->Scale(1/ttMC_syst_JERup_D1->Integral());
  ttMC_syst_JERup_D2->Scale(1/ttMC_syst_JERup_D2->Integral());
  ttMC_syst_JERup_D3->Scale(1/ttMC_syst_JERup_D3->Integral());
  ttMC_syst_JERup_D4->Scale(1/ttMC_syst_JERup_D4->Integral());

  ttMC_syst_JERup_D1->Divide(ttMC_syst_JERup_AVG);
  ttMC_syst_JERup_D2->Divide(ttMC_syst_JERup_AVG);
  ttMC_syst_JERup_D3->Divide(ttMC_syst_JERup_AVG);
  ttMC_syst_JERup_D4->Divide(ttMC_syst_JERup_AVG);

  // TCanvas *cJERup = new TCanvas("cJERup","cJERup",800,800);
  // cJERup->Divide(2,2);
  // cJERup->cd(1);
  // ttMC_syst_JERup_D1->Draw();
  // cJERup->cd(2);
  // ttMC_syst_JERup_D2->Draw();
  // cJERup->cd(3);
  // ttMC_syst_JERup_D3->Draw();
  // cJERup->cd(4);
  // ttMC_syst_JERup_D4->Draw();
  // cJERup->SaveAs("ratios_JERup.png");

  // JER down

  TH1D* ttMC_syst_JERdown_D1 = (TH1D*)file->Get("D1_TT_h_njets_pt30_1l_JERDown");
  TH1D* ttMC_syst_JERdown_D2 = (TH1D*)file->Get("D2_TT_h_njets_pt30_1l_JERDown");
  TH1D* ttMC_syst_JERdown_D3 = (TH1D*)file->Get("D3_TT_h_njets_pt30_1l_JERDown");
  TH1D* ttMC_syst_JERdown_D4 = (TH1D*)file->Get("D4_TT_h_njets_pt30_1l_JERDown");

  TH1D* ttMC_syst_JERdown_AVG = (TH1D*)ttMC_syst_JERdown_D1->Clone("ttMC_syst_JERdown_AVG");
  ttMC_syst_JERdown_AVG->Add(ttMC_syst_JERdown_D2);
  ttMC_syst_JERdown_AVG->Add(ttMC_syst_JERdown_D3);
  ttMC_syst_JERdown_AVG->Add(ttMC_syst_JERdown_D4);
  ttMC_syst_JERdown_AVG->Scale(1/ttMC_syst_JERdown_AVG->Integral());

  ttMC_syst_JERdown_D1->Scale(1/ttMC_syst_JERdown_D1->Integral());
  ttMC_syst_JERdown_D2->Scale(1/ttMC_syst_JERdown_D2->Integral());
  ttMC_syst_JERdown_D3->Scale(1/ttMC_syst_JERdown_D3->Integral());
  ttMC_syst_JERdown_D4->Scale(1/ttMC_syst_JERdown_D4->Integral());

  ttMC_syst_JERdown_D1->Divide(ttMC_syst_JERdown_AVG);
  ttMC_syst_JERdown_D2->Divide(ttMC_syst_JERdown_AVG);
  ttMC_syst_JERdown_D3->Divide(ttMC_syst_JERdown_AVG);
  ttMC_syst_JERdown_D4->Divide(ttMC_syst_JERdown_AVG);

  // TCanvas *cJERdown = new TCanvas("cJERdown","cJERdown",800,800);
  // cJERdown->Divide(2,2);
  // cJERdown->cd(1);
  // ttMC_syst_JERdown_D1->Draw();
  // cJERdown->cd(2);
  // ttMC_syst_JERdown_D2->Draw();
  // cJERdown->cd(3);
  // ttMC_syst_JERdown_D3->Draw();
  // cJERdown->cd(4);
  // ttMC_syst_JERdown_D4->Draw();
  // cJERdown->SaveAs("ratios_JERdown.png");

  // ----------------------  MVA bin 1  ------------------

  // Dataset with 8 bins
  RooDataHist data_hist_D1("data_obs_D1","Data observed in MVA bin 1",vars_D1,data_th1_D1);
  wspace->import(data_hist_D1);

  // ttbar bkg in D1
  RooArgList parlist_D1(N7_tt_D1,a0_tt,a1_tt,a2_tt);  // list of shape parameters for tt bkg
  //RooArgList parlist_D1(N7_tt_D1,a0_tt_D1,a1_tt_D1,a2_tt_D1);  // list of shape parameters for tt bkg
  RooArgList bkg_tt_syst_NP_D1(*wspace->var("np_tt_JECup"),*wspace->var("np_tt_JECdown"),
			       *wspace->var("np_tt_JERup"),*wspace->var("np_tt_JERdown"));  // list of nuisance parameters for tt bkg
  TH1D* bkg_tt_syst_histos_D1[20];  // array of histograms containing each tt bkg shape uncertainty
  bkg_tt_syst_histos_D1[0] = ttMC_syst_JECup_D1;
  bkg_tt_syst_histos_D1[1] = ttMC_syst_JECdown_D1;
  bkg_tt_syst_histos_D1[2] = ttMC_syst_JERup_D1;
  bkg_tt_syst_histos_D1[3] = ttMC_syst_JERdown_D1;
  RooArgList *bkg_tt_bins_D1 = new RooArgList();
  string procName_D1 = "background_tt_D1";
  construct_formula(procName_D1,*bkg_tt_bins_D1,parlist_D1,bkg_tt_syst_NP_D1,bkg_tt_syst_histos_D1);
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
  RooArgList bkg_tt_syst_NP_D2(*wspace->var("np_tt_JECup"),*wspace->var("np_tt_JECdown"),
			       *wspace->var("np_tt_JERup"),*wspace->var("np_tt_JERdown"));  // list of nuisance parameters for tt bkg
  TH1D* bkg_tt_syst_histos_D2[20];  // array of histograms containing each tt bkg shape uncertainty
  bkg_tt_syst_histos_D2[0] = ttMC_syst_JECup_D2;
  bkg_tt_syst_histos_D2[1] = ttMC_syst_JECdown_D2;
  bkg_tt_syst_histos_D2[2] = ttMC_syst_JERup_D2;
  bkg_tt_syst_histos_D2[3] = ttMC_syst_JERdown_D2;
  RooArgList *bkg_tt_bins_D2 = new RooArgList();
  string procName_D2 = "background_tt_D2";
  construct_formula(procName_D2,*bkg_tt_bins_D2,parlist_D2,bkg_tt_syst_NP_D2,bkg_tt_syst_histos_D2);
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
  RooArgList bkg_tt_syst_NP_D3(*wspace->var("np_tt_JECup"),*wspace->var("np_tt_JECdown"),
			       *wspace->var("np_tt_JERup"),*wspace->var("np_tt_JERdown"));  // list of nuisance parameters for tt bkg
  TH1D* bkg_tt_syst_histos_D3[20];  // array of histograms containing each tt bkg shape uncertainty
  bkg_tt_syst_histos_D3[0] = ttMC_syst_JECup_D3;
  bkg_tt_syst_histos_D3[1] = ttMC_syst_JECdown_D3;
  bkg_tt_syst_histos_D3[2] = ttMC_syst_JERup_D3;
  bkg_tt_syst_histos_D3[3] = ttMC_syst_JERdown_D3;
  RooArgList *bkg_tt_bins_D3 = new RooArgList();
  string procName_D3 = "background_tt_D3";
  construct_formula(procName_D3,*bkg_tt_bins_D3,parlist_D3,bkg_tt_syst_NP_D3,bkg_tt_syst_histos_D3);
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
  RooArgList bkg_tt_syst_NP_D4(*wspace->var("np_tt_JECup"),*wspace->var("np_tt_JECdown"),
			       *wspace->var("np_tt_JERup"),*wspace->var("np_tt_JERdown"));  // list of nuisance parameters for tt bkg
  TH1D* bkg_tt_syst_histos_D4[20];  // array of histograms containing each tt bkg shape uncertainty
  bkg_tt_syst_histos_D4[0] = ttMC_syst_JECup_D4;
  bkg_tt_syst_histos_D4[1] = ttMC_syst_JECdown_D4;
  bkg_tt_syst_histos_D4[2] = ttMC_syst_JERup_D4;
  bkg_tt_syst_histos_D4[3] = ttMC_syst_JERdown_D4;
  RooArgList *bkg_tt_bins_D4 = new RooArgList();
  string procName_D4 = "background_tt_D4";
  construct_formula(procName_D4,*bkg_tt_bins_D4,parlist_D4,bkg_tt_syst_NP_D4,bkg_tt_syst_histos_D4);
  RooParametricHist background_tt_D4(procName_D4.c_str(),"",*wspace->var("CMS_th1x"),*bkg_tt_bins_D4,*data_th1_D4);
  wspace->import(background_tt_D4,RooFit::RecycleConflictNodes());
  stringstream procNameD4Norm;
  procNameD4Norm << procName_D4 << "_norm";
  RooAddition tt_norm_D4(procNameD4Norm.str().c_str(),"",*bkg_tt_bins_D4);
  wspace->import(tt_norm_D4,RooFit::RecycleConflictNodes());

  // =================================================================================

   
  fOut->cd();

  otherMC_th1_D1->SetName("otherMC_th1_D1");
  sigMC_th1_D1->SetName("sigMC_th1_D1");
  otherMC_th1_D2->SetName("otherMC_th1_D2");
  sigMC_th1_D2->SetName("sigMC_th1_D2");
  otherMC_th1_D3->SetName("otherMC_th1_D3");
  sigMC_th1_D3->SetName("sigMC_th1_D3");
  otherMC_th1_D4->SetName("otherMC_th1_D4");
  sigMC_th1_D4->SetName("sigMC_th1_D4");

  otherMC_th1_D1->Write();
  sigMC_th1_D1->Write();
  otherMC_th1_D2->Write();
  sigMC_th1_D2->Write();
  otherMC_th1_D3->Write();
  sigMC_th1_D3->Write();
  otherMC_th1_D4->Write();
  sigMC_th1_D4->Write();

  TH1D* D1_SIG_JERUp = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_h_njets_pt30_1l_JERUp").c_str());
  TH1D* D1_SIG_JERDown = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_h_njets_pt30_1l_JERDown").c_str());
  TH1D* D2_SIG_JERUp = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_h_njets_pt30_1l_JERUp").c_str());
  TH1D* D2_SIG_JERDown = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_h_njets_pt30_1l_JERDown").c_str());
  TH1D* D3_SIG_JERUp = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_h_njets_pt30_1l_JERUp").c_str());
  TH1D* D3_SIG_JERDown = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_h_njets_pt30_1l_JERDown").c_str());
  TH1D* D4_SIG_JERUp = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_h_njets_pt30_1l_JERUp").c_str());
  TH1D* D4_SIG_JERDown = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_h_njets_pt30_1l_JERDown").c_str());

  D1_SIG_JERUp->SetName("D1_SIG_JERUp");
  D1_SIG_JERDown->SetName("D1_SIG_JERDown");
  D2_SIG_JERUp->SetName("D2_SIG_JERUp");
  D2_SIG_JERDown->SetName("D2_SIG_JERDown");
  D3_SIG_JERUp->SetName("D3_SIG_JERUp");
  D3_SIG_JERDown->SetName("D3_SIG_JERDown");
  D4_SIG_JERUp->SetName("D4_SIG_JERUp");
  D4_SIG_JERDown->SetName("D4_SIG_JERDown");

  D1_SIG_JERUp->Write();
  D1_SIG_JERDown->Write();
  D2_SIG_JERUp->Write();
  D2_SIG_JERDown->Write();
  D3_SIG_JERUp->Write();
  D3_SIG_JERDown->Write();
  D4_SIG_JERUp->Write();
  D4_SIG_JERDown->Write();

  TH1D* D1_SIG_JECUp = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_h_njets_pt30_1l_JECUp").c_str());
  TH1D* D1_SIG_JECDown = (TH1D*)file->Get(("D1_"+model+"_"+mass+"_h_njets_pt30_1l_JECDown").c_str());
  TH1D* D2_SIG_JECUp = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_h_njets_pt30_1l_JECUp").c_str());
  TH1D* D2_SIG_JECDown = (TH1D*)file->Get(("D2_"+model+"_"+mass+"_h_njets_pt30_1l_JECDown").c_str());
  TH1D* D3_SIG_JECUp = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_h_njets_pt30_1l_JECUp").c_str());
  TH1D* D3_SIG_JECDown = (TH1D*)file->Get(("D3_"+model+"_"+mass+"_h_njets_pt30_1l_JECDown").c_str());
  TH1D* D4_SIG_JECUp = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_h_njets_pt30_1l_JECUp").c_str());
  TH1D* D4_SIG_JECDown = (TH1D*)file->Get(("D4_"+model+"_"+mass+"_h_njets_pt30_1l_JECDown").c_str());

  D1_SIG_JECUp->SetName("D1_SIG_JECUp");
  D1_SIG_JECDown->SetName("D1_SIG_JECDown");
  D2_SIG_JECUp->SetName("D2_SIG_JECUp");
  D2_SIG_JECDown->SetName("D2_SIG_JECDown");
  D3_SIG_JECUp->SetName("D3_SIG_JECUp");
  D3_SIG_JECDown->SetName("D3_SIG_JECDown");
  D4_SIG_JECUp->SetName("D4_SIG_JECUp");
  D4_SIG_JECDown->SetName("D4_SIG_JECDown");

  D1_SIG_JECUp->Write();
  D1_SIG_JECDown->Write();
  D2_SIG_JECUp->Write();
  D2_SIG_JECDown->Write();
  D3_SIG_JECUp->Write();
  D3_SIG_JECDown->Write();
  D4_SIG_JECUp->Write();
  D4_SIG_JECDown->Write();

  TH1D* D1_OTHER_JERUp = (TH1D*)file->Get("D1_OTHER_h_njets_pt30_1l_JERUp");
  TH1D* D1_OTHER_JERDown = (TH1D*)file->Get("D1_OTHER_h_njets_pt30_1l_JERDown");
  TH1D* D2_OTHER_JERUp = (TH1D*)file->Get("D2_OTHER_h_njets_pt30_1l_JERUp");
  TH1D* D2_OTHER_JERDown = (TH1D*)file->Get("D2_OTHER_h_njets_pt30_1l_JERDown");
  TH1D* D3_OTHER_JERUp = (TH1D*)file->Get("D3_OTHER_h_njets_pt30_1l_JERUp");
  TH1D* D3_OTHER_JERDown = (TH1D*)file->Get("D3_OTHER_h_njets_pt30_1l_JERDown");
  TH1D* D4_OTHER_JERUp = (TH1D*)file->Get("D4_OTHER_h_njets_pt30_1l_JERUp");
  TH1D* D4_OTHER_JERDown = (TH1D*)file->Get("D4_OTHER_h_njets_pt30_1l_JERDown");

  D1_OTHER_JERUp->SetName("D1_OTHER_JERUp");
  D1_OTHER_JERDown->SetName("D1_OTHER_JERDown");
  D2_OTHER_JERUp->SetName("D2_OTHER_JERUp");
  D2_OTHER_JERDown->SetName("D2_OTHER_JERDown");
  D3_OTHER_JERUp->SetName("D3_OTHER_JERUp");
  D3_OTHER_JERDown->SetName("D3_OTHER_JERDown");
  D4_OTHER_JERUp->SetName("D4_OTHER_JERUp");
  D4_OTHER_JERDown->SetName("D4_OTHER_JERDown");

  D1_OTHER_JERUp->Write();
  D1_OTHER_JERDown->Write();
  D2_OTHER_JERUp->Write();
  D2_OTHER_JERDown->Write();
  D3_OTHER_JERUp->Write();
  D3_OTHER_JERDown->Write();
  D4_OTHER_JERUp->Write();
  D4_OTHER_JERDown->Write();

  TH1D* D1_OTHER_JECUp = (TH1D*)file->Get("D1_OTHER_h_njets_pt30_1l_JECUp");
  TH1D* D1_OTHER_JECDown = (TH1D*)file->Get("D1_OTHER_h_njets_pt30_1l_JECDown");
  TH1D* D2_OTHER_JECUp = (TH1D*)file->Get("D2_OTHER_h_njets_pt30_1l_JECUp");
  TH1D* D2_OTHER_JECDown = (TH1D*)file->Get("D2_OTHER_h_njets_pt30_1l_JECDown");
  TH1D* D3_OTHER_JECUp = (TH1D*)file->Get("D3_OTHER_h_njets_pt30_1l_JECUp");
  TH1D* D3_OTHER_JECDown = (TH1D*)file->Get("D3_OTHER_h_njets_pt30_1l_JECDown");
  TH1D* D4_OTHER_JECUp = (TH1D*)file->Get("D4_OTHER_h_njets_pt30_1l_JECUp");
  TH1D* D4_OTHER_JECDown = (TH1D*)file->Get("D4_OTHER_h_njets_pt30_1l_JECDown");

  D1_OTHER_JECUp->SetName("D1_OTHER_JECUp");
  D1_OTHER_JECDown->SetName("D1_OTHER_JECDown");
  D2_OTHER_JECUp->SetName("D2_OTHER_JECUp");
  D2_OTHER_JECDown->SetName("D2_OTHER_JECDown");
  D3_OTHER_JECUp->SetName("D3_OTHER_JECUp");
  D3_OTHER_JECDown->SetName("D3_OTHER_JECDown");
  D4_OTHER_JECUp->SetName("D4_OTHER_JECUp");
  D4_OTHER_JECDown->SetName("D4_OTHER_JECDown");

  D1_OTHER_JECUp->Write();
  D1_OTHER_JECDown->Write();
  D2_OTHER_JECUp->Write();
  D2_OTHER_JECDown->Write();
  D3_OTHER_JECUp->Write();
  D3_OTHER_JECDown->Write();
  D4_OTHER_JECUp->Write();
  D4_OTHER_JECDown->Write();


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
