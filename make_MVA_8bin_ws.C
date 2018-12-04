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


void construct_formula(string procName, RooArgList& binlist, const RooArgList& paramlist, const RooArgList& NPs, TH1D* h_syst[]) {

  int max_bin = 20; // 14 means just njets=14, 20 means last bin is inclusive up through njets=20

  for (int i=1; i<=8; i++) {

    stringstream form;
    RooArgList formArgList;

    form << "(@0";
    formArgList.add(paramlist[0]); // N7_tt for this MVA bin

    if (i>=2) {
      form << "*@1";
      formArgList.add(paramlist[1]); // p0_tt
    }

    if (i>=3) { // for bin 3 and up
      for (int j=3; j<=i; j++) {
	form << "*(@2+(@1-@2)*exp((" << j << "-2)*@3))";
      }
      formArgList.add(paramlist[2]); // p1_tt
      formArgList.add(paramlist[3]); // p2_tt
    } // end bin 3 and up

    // The last bin covers from njet=14 through njet=max_bin
    if (i==8) {
      for (int k=9; k<=max_bin-6; k++) {
	form << " + @0*@1";
	for (int j=3; j<=k; j++) {
	  form << "*(@2+(@1-@2)*exp((" << j << "-2)*@3))";
	}
      }
    }
    form << ")";

    if (i==1) {
      form << "*TMath::Power(" << h_syst[0]->GetBinContent(i) << ",@1)";
      form << "*TMath::Power(" << h_syst[1]->GetBinContent(i) << ",@2)";
      form << "*TMath::Power(" << h_syst[2]->GetBinContent(i) << ",@3)";
      form << "*TMath::Power(" << h_syst[3]->GetBinContent(i) << ",@4)";
    } else if (i==2) {
      form << "*TMath::Power(" << h_syst[0]->GetBinContent(i) << ",@2)";
      form << "*TMath::Power(" << h_syst[1]->GetBinContent(i) << ",@3)";
      form << "*TMath::Power(" << h_syst[2]->GetBinContent(i) << ",@4)";
      form << "*TMath::Power(" << h_syst[3]->GetBinContent(i) << ",@5)";
    } else if (i>=3) {
      form << "*TMath::Power(" << h_syst[0]->GetBinContent(i) << ",@4)";
      form << "*TMath::Power(" << h_syst[1]->GetBinContent(i) << ",@5)";
      form << "*TMath::Power(" << h_syst[2]->GetBinContent(i) << ",@6)";
      form << "*TMath::Power(" << h_syst[3]->GetBinContent(i) << ",@7)";
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

    cout << "bin i = " << i << " , njets = " << i+6 << endl;
    cout << "process bin name : " << binName.str().c_str() << endl;
    cout << "Formula : " << form.str().c_str() << endl;
    formArgList.Print();
    cout << endl;

  }

}


void make_MVA_8bin_ws() {
  using namespace RooFit;
  // Load the combine library to get access to the RooParametricHist
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
    
  // Output file and workspace 
  TFile *fOut = new TFile("MVA_ws.root","RECREATE");
  RooWorkspace *wspace = new RooWorkspace("wspace","wspace");

  // njet is our variable, 8 bins, from 7 up through 14,
  //   Note that njet=14 is inclusive as >=14
  // D1, D2, D3, D4 are the MVA bins
  wspace->factory("nj_D1[6.5,14.5]");
  //wspace->factory("nj_D1[0,8]");
  wspace->var("nj_D1")->setBins(8);
  RooArgSet vars_D1(*wspace->var("nj_D1"));

  wspace->factory("nj_D2[6.5,14.5]");
  //wspace->factory("nj_D2[0,8]");
  wspace->var("nj_D2")->setBins(8);
  RooArgSet vars_D2(*wspace->var("nj_D2"));

  wspace->factory("nj_D3[6.5,14.5]");
  //wspace->factory("nj_D3[0,8]");
  wspace->var("nj_D3")->setBins(8);
  RooArgSet vars_D3(*wspace->var("nj_D3"));

  wspace->factory("nj_D4[6.5,14.5]");
  //wspace->factory("nj_D4[0,8]");
  wspace->var("nj_D4")->setBins(8);
  RooArgSet vars_D4(*wspace->var("nj_D4"));


  //wspace->factory("CMS_th1x[0,8]");
  //wspace->var("CMS_th1x")->setBins(8);
  //RooArgSet vars_D1(*wspace->var("CMS_th1x"));
  //RooArgSet vars_D2(*wspace->var("CMS_th1x"));
  //RooArgSet vars_D3(*wspace->var("CMS_th1x"));
  //RooArgSet vars_D4(*wspace->var("CMS_th1x"));

    
  // file for obtaining histograms
  TFile* file = TFile::Open("njets_for_Aron_V1.2.3_JEC_Nov28.root");

  TH1* data_th1_D1 = 0;
  file->GetObject("h_njets_pt30_1l_deepESMbin1_pseudodataS_RPV_550",data_th1_D1);
  //file->GetObject("h_njets_pt30_1l_deepESMbin1_pseudodata",data_th1_D1);
  //file->GetObject("h_njets_pt30_1l_deepESMbin1_TT",data_th1_D1);
  TH1* ttMC_th1_D1 = 0;
  file->GetObject("h_njets_pt30_1l_deepESMbin1_TT",ttMC_th1_D1);
  TH1D* ttMC_syst_th1_D1 = (TH1D*)ttMC_th1_D1->Clone("ttMC_syst_th1_D1");
  //file->GetObject("tt_D1_syst",ttMC_syst_th1_D1);
  TH1* otherMC_th1_D1 = 0;
  file->GetObject("h_njets_pt30_1l_deepESMbin1_other",otherMC_th1_D1);
  TH1* sigMC_th1_D1 = 0;
  file->GetObject("h_njets_pt30_1l_deepESMbin1_RPV_550",sigMC_th1_D1);

  TH1* data_th1_D2 = 0;
  file->GetObject("h_njets_pt30_1l_deepESMbin2_pseudodataS_RPV_550",data_th1_D2);
  //file->GetObject("h_njets_pt30_1l_deepESMbin2_pseudodata",data_th1_D2);
  //file->GetObject("h_njets_pt30_1l_deepESMbin2_TT",data_th1_D2);
  TH1* ttMC_th1_D2 = 0;
  file->GetObject("h_njets_pt30_1l_deepESMbin2_TT",ttMC_th1_D2);
  TH1D* ttMC_syst_th1_D2 = (TH1D*)ttMC_th1_D2->Clone("ttMC_syst_th1_D2");
  //file->GetObject("tt_D2_syst",ttMC_syst_th1_D2);
  TH1* otherMC_th1_D2 = 0;
  file->GetObject("h_njets_pt30_1l_deepESMbin2_other",otherMC_th1_D2);
  TH1* sigMC_th1_D2 = 0;
  file->GetObject("h_njets_pt30_1l_deepESMbin2_RPV_550",sigMC_th1_D2);

  TH1* data_th1_D3 = 0;
  file->GetObject("h_njets_pt30_1l_deepESMbin3_pseudodataS_RPV_550",data_th1_D3);
  //file->GetObject("h_njets_pt30_1l_deepESMbin3_pseudodata",data_th1_D3);
  //file->GetObject("h_njets_pt30_1l_deepESMbin3_TT",data_th1_D3);
  TH1* ttMC_th1_D3 = 0;
  file->GetObject("h_njets_pt30_1l_deepESMbin3_TT",ttMC_th1_D3);
  TH1D* ttMC_syst_th1_D3 = (TH1D*)ttMC_th1_D3->Clone("ttMC_syst_th1_D3");
  //file->GetObject("tt_D3_syst",ttMC_syst_th1_D3);
  TH1* otherMC_th1_D3 = 0;
  file->GetObject("h_njets_pt30_1l_deepESMbin3_other",otherMC_th1_D3);
  TH1* sigMC_th1_D3 = 0;
  file->GetObject("h_njets_pt30_1l_deepESMbin3_RPV_550",sigMC_th1_D3);

  TH1* data_th1_D4 = 0;
  file->GetObject("h_njets_pt30_1l_deepESMbin4_pseudodataS_RPV_550",data_th1_D4);
  //file->GetObject("h_njets_pt30_1l_deepESMbin4_pseudodata",data_th1_D4);
  //file->GetObject("h_njets_pt30_1l_deepESMbin4_TT",data_th1_D4);
  TH1* ttMC_th1_D4 = 0;
  file->GetObject("h_njets_pt30_1l_deepESMbin4_TT",ttMC_th1_D4);
  TH1D* ttMC_syst_th1_D4 = (TH1D*)ttMC_th1_D4->Clone("ttMC_syst_th1_D4");
  //file->GetObject("tt_D4_syst",ttMC_syst_th1_D4);
  TH1* otherMC_th1_D4 = 0;
  file->GetObject("h_njets_pt30_1l_deepESMbin4_other",otherMC_th1_D4);
  TH1* sigMC_th1_D4 = 0;
  file->GetObject("h_njets_pt30_1l_deepESMbin4_RPV_550",sigMC_th1_D4);


  // tt bkg param setup
  RooRealVar p0_tt("p0_tt","p0 of tt bkg shape",0.35,0.0,1.0);
  RooRealVar p1_tt("p1_tt","p1 of tt bkg shape",0.1);
  RooRealVar p2_tt("p2_tt","p2 of tt bkg shape",-0.25,-1.0,0.0);


  double_t n7_tt_portion_D1 = data_th1_D1->GetBinContent(1) - otherMC_th1_D1->GetBinContent(1) - sigMC_th1_D1->GetBinContent(1);
  RooRealVar N7_tt_D1("N7_tt_D1","njets 7 for tt bkg in MVA D1",n7_tt_portion_D1,n7_tt_portion_D1-1500,n7_tt_portion_D1+1500);

  double_t n7_tt_portion_D2 = data_th1_D2->GetBinContent(1) - otherMC_th1_D2->GetBinContent(1) - sigMC_th1_D2->GetBinContent(1);
  RooRealVar N7_tt_D2("N7_tt_D2","njets 7 for tt bkg in MVA D2",n7_tt_portion_D2,n7_tt_portion_D2-1500,n7_tt_portion_D2+1500);

  double_t n7_tt_portion_D3 = data_th1_D3->GetBinContent(1) - otherMC_th1_D3->GetBinContent(1) - sigMC_th1_D3->GetBinContent(1);
  RooRealVar N7_tt_D3("N7_tt_D3","njets 7 for tt bkg in MVA D3",n7_tt_portion_D3,n7_tt_portion_D3-1500,n7_tt_portion_D3+1500);

  double_t n7_tt_portion_D4 = data_th1_D4->GetBinContent(1) - otherMC_th1_D4->GetBinContent(1) - sigMC_th1_D4->GetBinContent(1);
  RooRealVar N7_tt_D4("N7_tt_D4","njets 7 for tt bkg in MVA D4",n7_tt_portion_D4,n7_tt_portion_D4-1500,n7_tt_portion_D4+1500);

  // tt shape systematic nuisance parameters
  wspace->factory("np_tt_JECup_D1[0.0]");
  wspace->factory("np_tt_JECup_D2[0.0]");
  wspace->factory("np_tt_JECup_D3[0.0]");
  wspace->factory("np_tt_JECup_D4[0.0]");
  wspace->factory("np_tt_JECdown_D1[0.0]");
  wspace->factory("np_tt_JECdown_D2[0.0]");
  wspace->factory("np_tt_JECdown_D3[0.0]");
  wspace->factory("np_tt_JECdown_D4[0.0]");
  wspace->factory("np_tt_JERup_D1[0.0]");
  wspace->factory("np_tt_JERup_D2[0.0]");
  wspace->factory("np_tt_JERup_D3[0.0]");
  wspace->factory("np_tt_JERup_D4[0.0]");
  wspace->factory("np_tt_JERdown_D1[0.0]");
  wspace->factory("np_tt_JERdown_D2[0.0]");
  wspace->factory("np_tt_JERdown_D3[0.0]");
  wspace->factory("np_tt_JERdown_D4[0.0]");

  // Assemble the tt systematic histograms, containing the ratio
  //   "actual after njet variation / average over the MVA bins after njet variation"

  TH1D* ttMC_syst_JECup_D1 = (TH1D*)ttMC_th1_D1->Clone("ttMC_syst_JECup_D1");
  TH1D* ttMC_syst_JECdown_D1 = (TH1D*)ttMC_th1_D1->Clone("ttMC_syst_JECdown_D1");
  TH1D* ttMC_syst_JERup_D1 = (TH1D*)ttMC_th1_D1->Clone("ttMC_syst_JERup_D1");
  TH1D* ttMC_syst_JERdown_D1 = (TH1D*)ttMC_th1_D1->Clone("ttMC_syst_JERdown_D1");

  TH1D* ttMC_syst_JECup_D2 = (TH1D*)ttMC_th1_D2->Clone("ttMC_syst_JECup_D2");
  TH1D* ttMC_syst_JECdown_D2 = (TH1D*)ttMC_th1_D2->Clone("ttMC_syst_JECdown_D2");
  TH1D* ttMC_syst_JERup_D2 = (TH1D*)ttMC_th1_D2->Clone("ttMC_syst_JERup_D2");
  TH1D* ttMC_syst_JERdown_D2 = (TH1D*)ttMC_th1_D2->Clone("ttMC_syst_JERdown_D2");

  TH1D* ttMC_syst_JECup_D3 = (TH1D*)ttMC_th1_D3->Clone("ttMC_syst_JECup_D3");
  TH1D* ttMC_syst_JECdown_D3 = (TH1D*)ttMC_th1_D3->Clone("ttMC_syst_JECdown_D3");
  TH1D* ttMC_syst_JERup_D3 = (TH1D*)ttMC_th1_D3->Clone("ttMC_syst_JERup_D3");
  TH1D* ttMC_syst_JERdown_D3 = (TH1D*)ttMC_th1_D3->Clone("ttMC_syst_JERdown_D3");

  TH1D* ttMC_syst_JECup_D4 = (TH1D*)ttMC_th1_D4->Clone("ttMC_syst_JECup_D4");
  TH1D* ttMC_syst_JECdown_D4 = (TH1D*)ttMC_th1_D4->Clone("ttMC_syst_JECdown_D4");
  TH1D* ttMC_syst_JERup_D4 = (TH1D*)ttMC_th1_D4->Clone("ttMC_syst_JERup_D4");
  TH1D* ttMC_syst_JERdown_D4 = (TH1D*)ttMC_th1_D4->Clone("ttMC_syst_JERdown_D4");
  
  // JEC up

  file->GetObject("h_njets_pt30_1l_deepESMbin1JECup_TT",ttMC_syst_JECup_D1);
  file->GetObject("h_njets_pt30_1l_deepESMbin2JECup_TT",ttMC_syst_JECup_D2);
  file->GetObject("h_njets_pt30_1l_deepESMbin3JECup_TT",ttMC_syst_JECup_D3);
  file->GetObject("h_njets_pt30_1l_deepESMbin4JECup_TT",ttMC_syst_JECup_D4);

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

  ttMC_syst_JECup_D3->SetBinContent(8,ttMC_syst_JECup_D3->GetBinContent(7));
  ttMC_syst_JECup_D4->SetBinContent(8,ttMC_syst_JECup_D4->GetBinContent(7));
  
  TCanvas *cJECup = new TCanvas("cJECup","cJECup",1200,1200);
  cJECup->Divide(2,2);
  cJECup->cd(1);
  ttMC_syst_JECup_D1->Draw();
  cJECup->cd(2);
  ttMC_syst_JECup_D2->Draw();
  cJECup->cd(3);
  ttMC_syst_JECup_D3->Draw();
  cJECup->cd(4);
  ttMC_syst_JECup_D4->Draw();

  // JEC down

  file->GetObject("h_njets_pt30_1l_deepESMbin1JECdown_TT",ttMC_syst_JECdown_D1);
  file->GetObject("h_njets_pt30_1l_deepESMbin2JECdown_TT",ttMC_syst_JECdown_D2);
  file->GetObject("h_njets_pt30_1l_deepESMbin3JECdown_TT",ttMC_syst_JECdown_D3);
  file->GetObject("h_njets_pt30_1l_deepESMbin4JECdown_TT",ttMC_syst_JECdown_D4);

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

  //ttMC_syst_JECdown_D3->SetBinContent(8,ttMC_syst_JECdown_D3->GetBinContent(7));
  ttMC_syst_JECdown_D4->SetBinContent(8,ttMC_syst_JECdown_D4->GetBinContent(7));
  
  TCanvas *cJECdown = new TCanvas("cJECdown","cJECdown",1200,1200);
  cJECdown->Divide(2,2);
  cJECdown->cd(1);
  ttMC_syst_JECdown_D1->Draw();
  cJECdown->cd(2);
  ttMC_syst_JECdown_D2->Draw();
  cJECdown->cd(3);
  ttMC_syst_JECdown_D3->Draw();
  cJECdown->cd(4);
  ttMC_syst_JECdown_D4->Draw();

  // JER up

  file->GetObject("h_njets_pt30_1l_deepESMbin1JERup_TT",ttMC_syst_JERup_D1);
  file->GetObject("h_njets_pt30_1l_deepESMbin2JERup_TT",ttMC_syst_JERup_D2);
  file->GetObject("h_njets_pt30_1l_deepESMbin3JERup_TT",ttMC_syst_JERup_D3);
  file->GetObject("h_njets_pt30_1l_deepESMbin4JERup_TT",ttMC_syst_JERup_D4);

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

  //ttMC_syst_JERup_D3->SetBinContent(8,ttMC_syst_JERup_D3->GetBinContent(7));
  ttMC_syst_JERup_D4->SetBinContent(8,ttMC_syst_JERup_D4->GetBinContent(7));
  
  TCanvas *cJERup = new TCanvas("cJERup","cJERup",1200,1200);
  cJERup->Divide(2,2);
  cJERup->cd(1);
  ttMC_syst_JERup_D1->Draw();
  cJERup->cd(2);
  ttMC_syst_JERup_D2->Draw();
  cJERup->cd(3);
  ttMC_syst_JERup_D3->Draw();
  cJERup->cd(4);
  ttMC_syst_JERup_D4->Draw();

  // JER down

  file->GetObject("h_njets_pt30_1l_deepESMbin1JERdown_TT",ttMC_syst_JERdown_D1);
  file->GetObject("h_njets_pt30_1l_deepESMbin2JERdown_TT",ttMC_syst_JERdown_D2);
  file->GetObject("h_njets_pt30_1l_deepESMbin3JERdown_TT",ttMC_syst_JERdown_D3);
  file->GetObject("h_njets_pt30_1l_deepESMbin4JERdown_TT",ttMC_syst_JERdown_D4);

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

  //ttMC_syst_JERdown_D3->SetBinContent(8,ttMC_syst_JERdown_D3->GetBinContent(7));
  ttMC_syst_JERdown_D4->SetBinContent(8,ttMC_syst_JERdown_D4->GetBinContent(7));
  
  TCanvas *cJERdown = new TCanvas("cJERdown","cJERdown",1200,1200);
  cJERdown->Divide(2,2);
  cJERdown->cd(1);
  ttMC_syst_JERdown_D1->Draw();
  cJERdown->cd(2);
  ttMC_syst_JERdown_D2->Draw();
  cJERdown->cd(3);
  ttMC_syst_JERdown_D3->Draw();
  cJERdown->cd(4);
  ttMC_syst_JERdown_D4->Draw();

  // ----------------------  MVA bin 1  ------------------

  // Dataset with 8 bins
  RooDataHist data_hist_D1("data_obs_D1","Data observed in MVA bin 1",vars_D1,data_th1_D1);
  wspace->import(data_hist_D1);

  // ttbar bkg in D1
  RooArgList parlist_D1(N7_tt_D1,p0_tt,p1_tt,p2_tt);  // list of shape parameters for tt bkg
  RooArgList bkg_tt_syst_NP_D1(*wspace->var("np_tt_JECup_D1"),*wspace->var("np_tt_JECdown_D1"),
			       *wspace->var("np_tt_JERup_D1"),*wspace->var("np_tt_JERdown_D1"));  // list of nuisance parameters for tt bkg
  TH1D* bkg_tt_syst_histos_D1[20];  // array of histograms containing each tt bkg shape uncertainty
  bkg_tt_syst_histos_D1[0] = ttMC_syst_JECup_D1;
  bkg_tt_syst_histos_D1[1] = ttMC_syst_JECdown_D1;
  bkg_tt_syst_histos_D1[2] = ttMC_syst_JERup_D1;
  bkg_tt_syst_histos_D1[3] = ttMC_syst_JERdown_D1;
  RooArgList *bkg_tt_bins_D1 = new RooArgList();
  string procName_D1 = "background_tt_D1";
  construct_formula(procName_D1,*bkg_tt_bins_D1,parlist_D1,bkg_tt_syst_NP_D1,bkg_tt_syst_histos_D1);
  RooParametricHist background_tt_D1(procName_D1.c_str(),"",*wspace->var("nj_D1"),*bkg_tt_bins_D1,*data_th1_D1);
  wspace->import(background_tt_D1,RooFit::RecycleConflictNodes());
  stringstream procNameD1Norm;
  procNameD1Norm << procName_D1 << "_norm";
  RooAddition tt_norm_D1(procNameD1Norm.str().c_str(),"",*bkg_tt_bins_D1);
  wspace->import(tt_norm_D1,RooFit::RecycleConflictNodes());

  // other bkg in D1
  // RooArgList *bkg_other_bins_D1 = new RooArgList();
  // RooRealVar *otherb1D1 = new RooRealVar("bkg_other_b1D1","",otherMC_th1_D1->GetBinContent(1));
  // RooRealVar *otherb2D1 = new RooRealVar("bkg_other_b2D1","",otherMC_th1_D1->GetBinContent(2));
  // RooRealVar *otherb3D1 = new RooRealVar("bkg_other_b3D1","",otherMC_th1_D1->GetBinContent(3));
  // RooRealVar *otherb4D1 = new RooRealVar("bkg_other_b4D1","",otherMC_th1_D1->GetBinContent(4));
  // RooRealVar *otherb5D1 = new RooRealVar("bkg_other_b5D1","",otherMC_th1_D1->GetBinContent(5));
  // RooRealVar *otherb6D1 = new RooRealVar("bkg_other_b6D1","",otherMC_th1_D1->GetBinContent(6));
  // RooRealVar *otherb7D1 = new RooRealVar("bkg_other_b7D1","",otherMC_th1_D1->GetBinContent(7));
  // RooRealVar *otherb8D1 = new RooRealVar("bkg_other_b8D1","",otherMC_th1_D1->GetBinContent(8));
  // bkg_other_bins_D1->add(*otherb1D1);
  // bkg_other_bins_D1->add(*otherb2D1);
  // bkg_other_bins_D1->add(*otherb3D1);
  // bkg_other_bins_D1->add(*otherb4D1);
  // bkg_other_bins_D1->add(*otherb5D1);
  // bkg_other_bins_D1->add(*otherb6D1);
  // bkg_other_bins_D1->add(*otherb7D1);
  // bkg_other_bins_D1->add(*otherb8D1);
  // RooParametricHist background_other_D1("background_other_D1","Other background PDF in MVA bin 1",*wspace->var("nj_D1"),*bkg_other_bins_D1,*data_th1_D1);
  // wspace->import(background_other_D1,RooFit::RecycleConflictNodes());
  // RooAddition other_norm_D1("background_other_D1_norm","",*bkg_other_bins_D1);
  // wspace->import(other_norm_D1,RooFit::RecycleConflictNodes());

  RooDataHist otherMC_hist_D1("otherMC_hist_D1","other MC observed in MVA bin 1",vars_D1,otherMC_th1_D1);
  wspace->import(*(new RooHistPdf("background_other_D1","",vars_D1,otherMC_hist_D1)));

  // signal in D1
  // RooArgList *sig_bins_D1 = new RooArgList();
  // RooRealVar *sigb1D1 = new RooRealVar("sig_b1D1","",sigMC_th1_D1->GetBinContent(1));
  // RooRealVar *sigb2D1 = new RooRealVar("sig_b2D1","",sigMC_th1_D1->GetBinContent(2));
  // RooRealVar *sigb3D1 = new RooRealVar("sig_b3D1","",sigMC_th1_D1->GetBinContent(3));
  // RooRealVar *sigb4D1 = new RooRealVar("sig_b4D1","",sigMC_th1_D1->GetBinContent(4));
  // RooRealVar *sigb5D1 = new RooRealVar("sig_b5D1","",sigMC_th1_D1->GetBinContent(5));
  // RooRealVar *sigb6D1 = new RooRealVar("sig_b6D1","",sigMC_th1_D1->GetBinContent(6));
  // RooRealVar *sigb7D1 = new RooRealVar("sig_b7D1","",sigMC_th1_D1->GetBinContent(7));
  // RooRealVar *sigb8D1 = new RooRealVar("sig_b8D1","",sigMC_th1_D1->GetBinContent(8));
  // sig_bins_D1->add(*sigb1D1);
  // sig_bins_D1->add(*sigb2D1);
  // sig_bins_D1->add(*sigb3D1);
  // sig_bins_D1->add(*sigb4D1);
  // sig_bins_D1->add(*sigb5D1);
  // sig_bins_D1->add(*sigb6D1);
  // sig_bins_D1->add(*sigb7D1);
  // sig_bins_D1->add(*sigb8D1);
  // RooParametricHist signal_D1("signal_D1","Signal PDF in MVA bin 1",*wspace->var("nj_D1"),*sig_bins_D1,*data_th1_D1);
  // wspace->import(signal_D1,RooFit::RecycleConflictNodes());
  // RooAddition signal_norm_D1("signal_D1_norm","",*sig_bins_D1);
  // wspace->import(signal_norm_D1,RooFit::RecycleConflictNodes());

  RooDataHist sigMC_hist_D1("sigMC_hist_D1","signal MC in MVA bin 1",vars_D1,sigMC_th1_D1);
  wspace->import(*(new RooHistPdf("signal_D1","",vars_D1,sigMC_hist_D1)));

  // ---------------------- MVA bin 2  ------------------

  // Dataset with 8 bins
  RooDataHist data_hist_D2("data_obs_D2","Data observed in MVA bin 2",vars_D2,data_th1_D2);
  wspace->import(data_hist_D2);

  // ttbar bkg in D2
  RooArgList parlist_D2(N7_tt_D2,p0_tt,p1_tt,p2_tt);  // list of shape parameters for tt bkg
  RooArgList bkg_tt_syst_NP_D2(*wspace->var("np_tt_JECup_D2"),*wspace->var("np_tt_JECdown_D2"),
			       *wspace->var("np_tt_JERup_D2"),*wspace->var("np_tt_JERdown_D2"));  // list of nuisance parameters for tt bkg
  TH1D* bkg_tt_syst_histos_D2[20];  // array of histograms containing each tt bkg shape uncertainty
  bkg_tt_syst_histos_D2[0] = ttMC_syst_JECup_D2;
  bkg_tt_syst_histos_D2[1] = ttMC_syst_JECdown_D2;
  bkg_tt_syst_histos_D2[2] = ttMC_syst_JERup_D2;
  bkg_tt_syst_histos_D2[3] = ttMC_syst_JERdown_D2;
  RooArgList *bkg_tt_bins_D2 = new RooArgList();
  string procName_D2 = "background_tt_D2";
  construct_formula(procName_D2,*bkg_tt_bins_D2,parlist_D2,bkg_tt_syst_NP_D2,bkg_tt_syst_histos_D2);
  RooParametricHist background_tt_D2(procName_D2.c_str(),"",*wspace->var("nj_D2"),*bkg_tt_bins_D2,*data_th1_D2);
  wspace->import(background_tt_D2,RooFit::RecycleConflictNodes());
  stringstream procNameD2Norm;
  procNameD2Norm << procName_D2 << "_norm";
  RooAddition tt_norm_D2(procNameD2Norm.str().c_str(),"",*bkg_tt_bins_D2);
  wspace->import(tt_norm_D2,RooFit::RecycleConflictNodes());

  // other bkg in D2
  // RooArgList *bkg_other_bins_D2 = new RooArgList();
  // RooRealVar *otherb1D2 = new RooRealVar("bkg_other_b1D2","",otherMC_th1_D2->GetBinContent(1));
  // RooRealVar *otherb2D2 = new RooRealVar("bkg_other_b2D2","",otherMC_th1_D2->GetBinContent(2));
  // RooRealVar *otherb3D2 = new RooRealVar("bkg_other_b3D2","",otherMC_th1_D2->GetBinContent(3));
  // RooRealVar *otherb4D2 = new RooRealVar("bkg_other_b4D2","",otherMC_th1_D2->GetBinContent(4));
  // RooRealVar *otherb5D2 = new RooRealVar("bkg_other_b5D2","",otherMC_th1_D2->GetBinContent(5));
  // RooRealVar *otherb6D2 = new RooRealVar("bkg_other_b6D2","",otherMC_th1_D2->GetBinContent(6));
  // RooRealVar *otherb7D2 = new RooRealVar("bkg_other_b7D2","",otherMC_th1_D2->GetBinContent(7));
  // RooRealVar *otherb8D2 = new RooRealVar("bkg_other_b8D2","",otherMC_th1_D2->GetBinContent(8));
  // bkg_other_bins_D2->add(*otherb1D2);
  // bkg_other_bins_D2->add(*otherb2D2);
  // bkg_other_bins_D2->add(*otherb3D2);
  // bkg_other_bins_D2->add(*otherb4D2);
  // bkg_other_bins_D2->add(*otherb5D2);
  // bkg_other_bins_D2->add(*otherb6D2);
  // bkg_other_bins_D2->add(*otherb7D2);
  // bkg_other_bins_D2->add(*otherb8D2);
  // RooParametricHist background_other_D2("background_other_D2","Other background PDF in MVA bin 2",*wspace->var("nj_D2"),*bkg_other_bins_D2,*data_th1_D2);
  // wspace->import(background_other_D2,RooFit::RecycleConflictNodes());
  // RooAddition other_norm_D2("background_other_D2_norm","",*bkg_other_bins_D2);
  // wspace->import(other_norm_D2,RooFit::RecycleConflictNodes());

  RooDataHist otherMC_hist_D2("otherMC_hist_D2","other MC observed in MVA bin 2",vars_D2,otherMC_th1_D2);
  wspace->import(*(new RooHistPdf("background_other_D2","",vars_D2,otherMC_hist_D2)));

  // signal in D2
  // RooArgList *sig_bins_D2 = new RooArgList();
  // RooRealVar *sigb1D2 = new RooRealVar("sig_b1D2","",sigMC_th1_D2->GetBinContent(1));
  // RooRealVar *sigb2D2 = new RooRealVar("sig_b2D2","",sigMC_th1_D2->GetBinContent(2));
  // RooRealVar *sigb3D2 = new RooRealVar("sig_b3D2","",sigMC_th1_D2->GetBinContent(3));
  // RooRealVar *sigb4D2 = new RooRealVar("sig_b4D2","",sigMC_th1_D2->GetBinContent(4));
  // RooRealVar *sigb5D2 = new RooRealVar("sig_b5D2","",sigMC_th1_D2->GetBinContent(5));
  // RooRealVar *sigb6D2 = new RooRealVar("sig_b6D2","",sigMC_th1_D2->GetBinContent(6));
  // RooRealVar *sigb7D2 = new RooRealVar("sig_b7D2","",sigMC_th1_D2->GetBinContent(7));
  // RooRealVar *sigb8D2 = new RooRealVar("sig_b8D2","",sigMC_th1_D2->GetBinContent(8));
  // sig_bins_D2->add(*sigb1D2);
  // sig_bins_D2->add(*sigb2D2);
  // sig_bins_D2->add(*sigb3D2);
  // sig_bins_D2->add(*sigb4D2);
  // sig_bins_D2->add(*sigb5D2);
  // sig_bins_D2->add(*sigb6D2);
  // sig_bins_D2->add(*sigb7D2);
  // sig_bins_D2->add(*sigb8D2);
  // RooParametricHist signal_D2("signal_D2","Signal PDF in MVA bin 2",*wspace->var("nj_D2"),*sig_bins_D2,*data_th1_D2);
  // wspace->import(signal_D2,RooFit::RecycleConflictNodes());
  // RooAddition signal_norm_D2("signal_D2_norm","",*sig_bins_D2);
  // wspace->import(signal_norm_D2,RooFit::RecycleConflictNodes());

  RooDataHist sigMC_hist_D2("sigMC_hist_D2","signal MC in MVA bin 2",vars_D2,sigMC_th1_D2);
  wspace->import(*(new RooHistPdf("signal_D2","",vars_D2,sigMC_hist_D2)));

  // ---------------------- MVA bin 3  ------------------

  // Dataset with 8 bins
  RooDataHist data_hist_D3("data_obs_D3","Data observed in MVA bin 3",vars_D3,data_th1_D3);
  wspace->import(data_hist_D3);

  // ttbar bkg in D3
  RooArgList parlist_D3(N7_tt_D3,p0_tt,p1_tt,p2_tt);  // list of shape parameters for tt bkg
  RooArgList bkg_tt_syst_NP_D3(*wspace->var("np_tt_JECup_D3"),*wspace->var("np_tt_JECdown_D3"),
			       *wspace->var("np_tt_JERup_D3"),*wspace->var("np_tt_JERdown_D3"));  // list of nuisance parameters for tt bkg
  TH1D* bkg_tt_syst_histos_D3[20];  // array of histograms containing each tt bkg shape uncertainty
  bkg_tt_syst_histos_D3[0] = ttMC_syst_JECup_D3;
  bkg_tt_syst_histos_D3[1] = ttMC_syst_JECdown_D3;
  bkg_tt_syst_histos_D3[2] = ttMC_syst_JERup_D3;
  bkg_tt_syst_histos_D3[3] = ttMC_syst_JERdown_D3;
  RooArgList *bkg_tt_bins_D3 = new RooArgList();
  string procName_D3 = "background_tt_D3";
  construct_formula(procName_D3,*bkg_tt_bins_D3,parlist_D3,bkg_tt_syst_NP_D3,bkg_tt_syst_histos_D3);
  RooParametricHist background_tt_D3(procName_D3.c_str(),"",*wspace->var("nj_D3"),*bkg_tt_bins_D3,*data_th1_D3);
  wspace->import(background_tt_D3,RooFit::RecycleConflictNodes());
  stringstream procNameD3Norm;
  procNameD3Norm << procName_D3 << "_norm";
  RooAddition tt_norm_D3(procNameD3Norm.str().c_str(),"",*bkg_tt_bins_D3);
  wspace->import(tt_norm_D3,RooFit::RecycleConflictNodes());

  // other bkg in D3
  // RooArgList *bkg_other_bins_D3 = new RooArgList();
  // RooRealVar *otherb1D3 = new RooRealVar("bkg_other_b1D3","",otherMC_th1_D3->GetBinContent(1));
  // RooRealVar *otherb2D3 = new RooRealVar("bkg_other_b2D3","",otherMC_th1_D3->GetBinContent(2));
  // RooRealVar *otherb3D3 = new RooRealVar("bkg_other_b3D3","",otherMC_th1_D3->GetBinContent(3));
  // RooRealVar *otherb4D3 = new RooRealVar("bkg_other_b4D3","",otherMC_th1_D3->GetBinContent(4));
  // RooRealVar *otherb5D3 = new RooRealVar("bkg_other_b5D3","",otherMC_th1_D3->GetBinContent(5));
  // RooRealVar *otherb6D3 = new RooRealVar("bkg_other_b6D3","",otherMC_th1_D3->GetBinContent(6));
  // RooRealVar *otherb7D3 = new RooRealVar("bkg_other_b7D3","",otherMC_th1_D3->GetBinContent(7));
  // RooRealVar *otherb8D3 = new RooRealVar("bkg_other_b8D3","",otherMC_th1_D3->GetBinContent(8));
  // bkg_other_bins_D3->add(*otherb1D3);
  // bkg_other_bins_D3->add(*otherb2D3);
  // bkg_other_bins_D3->add(*otherb3D3);
  // bkg_other_bins_D3->add(*otherb4D3);
  // bkg_other_bins_D3->add(*otherb5D3);
  // bkg_other_bins_D3->add(*otherb6D3);
  // bkg_other_bins_D3->add(*otherb7D3);
  // bkg_other_bins_D3->add(*otherb8D3);
  // RooParametricHist background_other_D3("background_other_D3","Other background PDF in MVA bin 3",*wspace->var("nj_D3"),*bkg_other_bins_D3,*data_th1_D3);
  // wspace->import(background_other_D3,RooFit::RecycleConflictNodes());
  // RooAddition other_norm_D3("background_other_D3_norm","",*bkg_other_bins_D3);
  // wspace->import(other_norm_D3,RooFit::RecycleConflictNodes());

  RooDataHist otherMC_hist_D3("otherMC_hist_D3","other MC observed in MVA bin 3",vars_D3,otherMC_th1_D3);
  wspace->import(*(new RooHistPdf("background_other_D3","",vars_D3,otherMC_hist_D3)));

  // signal in D3
  // RooArgList *sig_bins_D3 = new RooArgList();
  // RooRealVar *sigb1D3 = new RooRealVar("sig_b1D3","",sigMC_th1_D3->GetBinContent(1));
  // RooRealVar *sigb2D3 = new RooRealVar("sig_b2D3","",sigMC_th1_D3->GetBinContent(2));
  // RooRealVar *sigb3D3 = new RooRealVar("sig_b3D3","",sigMC_th1_D3->GetBinContent(3));
  // RooRealVar *sigb4D3 = new RooRealVar("sig_b4D3","",sigMC_th1_D3->GetBinContent(4));
  // RooRealVar *sigb5D3 = new RooRealVar("sig_b5D3","",sigMC_th1_D3->GetBinContent(5));
  // RooRealVar *sigb6D3 = new RooRealVar("sig_b6D3","",sigMC_th1_D3->GetBinContent(6));
  // RooRealVar *sigb7D3 = new RooRealVar("sig_b7D3","",sigMC_th1_D3->GetBinContent(7));
  // RooRealVar *sigb8D3 = new RooRealVar("sig_b8D3","",sigMC_th1_D3->GetBinContent(8));
  // sig_bins_D3->add(*sigb1D3);
  // sig_bins_D3->add(*sigb2D3);
  // sig_bins_D3->add(*sigb3D3);
  // sig_bins_D3->add(*sigb4D3);
  // sig_bins_D3->add(*sigb5D3);
  // sig_bins_D3->add(*sigb6D3);
  // sig_bins_D3->add(*sigb7D3);
  // sig_bins_D3->add(*sigb8D3);
  // RooParametricHist signal_D3("signal_D3","Signal PDF in MVA bin 3",*wspace->var("nj_D3"),*sig_bins_D3,*data_th1_D3);
  // wspace->import(signal_D3,RooFit::RecycleConflictNodes());
  // RooAddition signal_norm_D3("signal_D3_norm","",*sig_bins_D3);
  // wspace->import(signal_norm_D3,RooFit::RecycleConflictNodes());

  RooDataHist sigMC_hist_D3("sigMC_hist_D3","signa3 MC in MVA bin 3",vars_D3,sigMC_th1_D3);
  wspace->import(*(new RooHistPdf("signal_D3","",vars_D3,sigMC_hist_D3)));

  // ---------------------- MVA bin 4  ------------------

  // Dataset with 8 bins
  RooDataHist data_hist_D4("data_obs_D4","Data observed in MVA bin 4",vars_D4,data_th1_D4);
  wspace->import(data_hist_D4);

  // ttbar bkg in D4
  RooArgList parlist_D4(N7_tt_D4,p0_tt,p1_tt,p2_tt);  // list of shape parameters for tt bkg
  RooArgList bkg_tt_syst_NP_D4(*wspace->var("np_tt_JECup_D4"),*wspace->var("np_tt_JECdown_D4"),
			       *wspace->var("np_tt_JERup_D4"),*wspace->var("np_tt_JERdown_D4"));  // list of nuisance parameters for tt bkg
  TH1D* bkg_tt_syst_histos_D4[20];  // array of histograms containing each tt bkg shape uncertainty
  bkg_tt_syst_histos_D4[0] = ttMC_syst_JECup_D4;
  bkg_tt_syst_histos_D4[1] = ttMC_syst_JECdown_D4;
  bkg_tt_syst_histos_D4[2] = ttMC_syst_JERup_D4;
  bkg_tt_syst_histos_D4[3] = ttMC_syst_JERdown_D4;
  RooArgList *bkg_tt_bins_D4 = new RooArgList();
  string procName_D4 = "background_tt_D4";
  construct_formula(procName_D4,*bkg_tt_bins_D4,parlist_D4,bkg_tt_syst_NP_D4,bkg_tt_syst_histos_D4);
  RooParametricHist background_tt_D4(procName_D4.c_str(),"",*wspace->var("nj_D4"),*bkg_tt_bins_D4,*data_th1_D4);
  wspace->import(background_tt_D4,RooFit::RecycleConflictNodes());
  stringstream procNameD4Norm;
  procNameD4Norm << procName_D4 << "_norm";
  RooAddition tt_norm_D4(procNameD4Norm.str().c_str(),"",*bkg_tt_bins_D4);
  wspace->import(tt_norm_D4,RooFit::RecycleConflictNodes());

  // other bkg in D4
  // RooArgList *bkg_other_bins_D4 = new RooArgList();
  // RooRealVar *otherb1D4 = new RooRealVar("bkg_other_b1D4","",otherMC_th1_D4->GetBinContent(1));
  // RooRealVar *otherb2D4 = new RooRealVar("bkg_other_b2D4","",otherMC_th1_D4->GetBinContent(2));
  // RooRealVar *otherb3D4 = new RooRealVar("bkg_other_b3D4","",otherMC_th1_D4->GetBinContent(3));
  // RooRealVar *otherb4D4 = new RooRealVar("bkg_other_b4D4","",otherMC_th1_D4->GetBinContent(4));
  // RooRealVar *otherb5D4 = new RooRealVar("bkg_other_b5D4","",otherMC_th1_D4->GetBinContent(5));
  // RooRealVar *otherb6D4 = new RooRealVar("bkg_other_b6D4","",otherMC_th1_D4->GetBinContent(6));
  // RooRealVar *otherb7D4 = new RooRealVar("bkg_other_b7D4","",otherMC_th1_D4->GetBinContent(7));
  // RooRealVar *otherb8D4 = new RooRealVar("bkg_other_b8D4","",otherMC_th1_D4->GetBinContent(8));
  // bkg_other_bins_D4->add(*otherb1D4);
  // bkg_other_bins_D4->add(*otherb2D4);
  // bkg_other_bins_D4->add(*otherb3D4);
  // bkg_other_bins_D4->add(*otherb4D4);
  // bkg_other_bins_D4->add(*otherb5D4);
  // bkg_other_bins_D4->add(*otherb6D4);
  // bkg_other_bins_D4->add(*otherb7D4);
  // bkg_other_bins_D4->add(*otherb8D4);
  // RooParametricHist background_other_D4("background_other_D4","Other background PDF in MVA bin 4",*wspace->var("nj_D4"),*bkg_other_bins_D4,*data_th1_D4);
  // wspace->import(background_other_D4,RooFit::RecycleConflictNodes());
  // RooAddition other_norm_D4("background_other_D4_norm","",*bkg_other_bins_D4);
  // wspace->import(other_norm_D4,RooFit::RecycleConflictNodes());

  RooDataHist otherMC_hist_D4("otherMC_hist_D4","other MC observed in MVA bin 4",vars_D4,otherMC_th1_D4);
  wspace->import(*(new RooHistPdf("background_other_D4","",vars_D4,otherMC_hist_D4)));

  // signal in D4
  // RooArgList *sig_bins_D4 = new RooArgList();
  // RooRealVar *sigb1D4 = new RooRealVar("sig_b1D4","",sigMC_th1_D4->GetBinContent(1));
  // RooRealVar *sigb2D4 = new RooRealVar("sig_b2D4","",sigMC_th1_D4->GetBinContent(2));
  // RooRealVar *sigb3D4 = new RooRealVar("sig_b3D4","",sigMC_th1_D4->GetBinContent(3));
  // RooRealVar *sigb4D4 = new RooRealVar("sig_b4D4","",sigMC_th1_D4->GetBinContent(4));
  // RooRealVar *sigb5D4 = new RooRealVar("sig_b5D4","",sigMC_th1_D4->GetBinContent(5));
  // RooRealVar *sigb6D4 = new RooRealVar("sig_b6D4","",sigMC_th1_D4->GetBinContent(6));
  // RooRealVar *sigb7D4 = new RooRealVar("sig_b7D4","",sigMC_th1_D4->GetBinContent(7));
  // RooRealVar *sigb8D4 = new RooRealVar("sig_b8D4","",sigMC_th1_D4->GetBinContent(8));
  // sig_bins_D4->add(*sigb1D4);
  // sig_bins_D4->add(*sigb2D4);
  // sig_bins_D4->add(*sigb3D4);
  // sig_bins_D4->add(*sigb4D4);
  // sig_bins_D4->add(*sigb5D4);
  // sig_bins_D4->add(*sigb6D4);
  // sig_bins_D4->add(*sigb7D4);
  // sig_bins_D4->add(*sigb8D4);
  // RooParametricHist signal_D4("signal_D4","Signal PDF in MVA bin 4",*wspace->var("nj_D4"),*sig_bins_D4,*data_th1_D4);
  // wspace->import(signal_D4,RooFit::RecycleConflictNodes());
  // RooAddition signal_norm_D4("signal_D4_norm","",*sig_bins_D4);
  // wspace->import(signal_norm_D4,RooFit::RecycleConflictNodes());

  RooDataHist sigMC_hist_D4("sigMC_hist_D4","signal MC in MVA bin 4",vars_D4,sigMC_th1_D4);
  wspace->import(*(new RooHistPdf("signal_D4","",vars_D4,sigMC_hist_D4)));

  // =================================================================================

   
  fOut->cd();

  // otherMC_th1_D1->SetName("otherMC_th1_D1");
  // sigMC_th1_D1->SetName("sigMC_th1_D1");
  // otherMC_th1_D2->SetName("otherMC_th1_D2");
  // sigMC_th1_D2->SetName("sigMC_th1_D2");
  // otherMC_th1_D3->SetName("otherMC_th1_D3");
  // sigMC_th1_D3->SetName("sigMC_th1_D3");
  // otherMC_th1_D4->SetName("otherMC_th1_D4");
  // sigMC_th1_D4->SetName("sigMC_th1_D4");

  // otherMC_th1_D1->Write();
  // sigMC_th1_D1->Write();
  // otherMC_th1_D2->Write();
  // sigMC_th1_D2->Write();
  // otherMC_th1_D3->Write();
  // sigMC_th1_D3->Write();
  // otherMC_th1_D4->Write();
  // sigMC_th1_D4->Write();

  // TH1* otherUpSystD1 = 0;
  // file->GetObject("OtherD1Up",otherUpSystD1);
  // otherUpSystD1->Write();
  // TH1* otherDownSystD1 = 0;
  // file->GetObject("OtherD1Down",otherDownSystD1);
  // otherDownSystD1->Write();

  // TH1* otherUpSystD2 = 0;
  // file->GetObject("OtherD2Up",otherUpSystD2);
  // otherUpSystD2->Write();
  // TH1* otherDownSystD2 = 0;
  // file->GetObject("OtherD2Down",otherDownSystD2);
  // otherDownSystD2->Write();

  // TH1* otherUpSystD3 = 0;
  // file->GetObject("OtherD3Up",otherUpSystD3);
  // otherUpSystD3->Write();
  // TH1* otherDownSystD3 = 0;
  // file->GetObject("OtherD3Down",otherDownSystD3);
  // otherDownSystD3->Write();

  // TH1* otherUpSystD4 = 0;
  // file->GetObject("OtherD4Up",otherUpSystD4);
  // otherUpSystD4->Write();
  // TH1* otherDownSystD4 = 0;
  // file->GetObject("OtherD4Down",otherDownSystD4);
  // otherDownSystD4->Write();



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
