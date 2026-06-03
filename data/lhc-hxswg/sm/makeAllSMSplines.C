
const std::string location = "./";  //"$CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/data/lhc-hxswg/sm/";
const int nxs = 9;
const int nbr = 11;
const int nbr4f = 14;

void makeBRSplines(){

	RooRealVar MH("MH","MH",125.09,120,130); MH.setConstant();
        TFile *fout = new TFile("sm_br_yr4.root","RECREATE");
	
	RooWorkspace ws("br","br");

	std::string br[nbr]   =  {"hbb","htt","hmm","hcc","hss","htoptop","hgluglu","hgg","hzg","hww","hzz"};

	for (int bri=0;bri<nbr;bri++){
		std::string name = br[bri];
		RooSpline1D spline(Form("%s",name.c_str()),Form("file br/BR4.txt, x=%d, y=%d",0,1),MH,Form("%s/br/BR4.txt",location.c_str()),0,bri+1,1,"CSPLINE");
		ws.import(spline);
	}

	// 4 fermions too
	std::string br4f[nbr4f] = {"llll_e_m_t","llll_e_m","4e","2e2m","llvv_e_m_t","llvv_e_m","evev","evmv","llqq_e_m_t","llqq_e_m","l+vqq_e_m","vvqq","qqqq","ffff"};

	for (int bri=0;bri<nbr4f;bri++){
		std::string name = br4f[bri];
		RooSpline1D spline(Form("%s",name.c_str()),Form("file br/BR4_4f.txt, x=%d, y=%d",0,1),MH,Form("%s/br/BR4_4f.txt",location.c_str()),0,bri+1,1,"CSPLINE");
		ws.import(spline);
	}


	ws.Print();
	fout->cd();
	ws.Write();
	fout->Close();
}

void makeXSSplines(std::string sqrts="13"){

        TFile *fout = new TFile(Form("sm_yr5_%sTeV.root",sqrts.c_str()),"RECREATE");
	RooWorkspace ws(Form("xs_%sTeV",sqrts.c_str()),Form("xs_%sTeV",sqrts.c_str()));

	RooRealVar MH("MH","MH",125.09,120,130); MH.setConstant();

	std::string xs[nxs] = {"WH","ZH","bbH","ggH","ttH","vbfH","ggZH","tHW","tHq"};

	for (int xsi=0;xsi<nxs;xsi++){
		std::string name = xs[xsi];
		if (name=="WH"){
		  RooSpline1D splineP(Form("WplusH_%sTeV",sqrts.c_str()),Form("file %sTeV/%sTeV-%s.txt, x=%d, y=%d",sqrts.c_str(),sqrts.c_str(),name.c_str(),0,7),MH,Form("%s/xs/%sTeV/%sTeV-%s.txt",location.c_str(),sqrts.c_str(),sqrts.c_str(),name.c_str()),0,7,1,"CSPLINE");		
		  RooSpline1D splineM(Form("WminusH_%sTeV",sqrts.c_str()),Form("file %sTeV/%sTeV-%s.txt, x=%d, y=%d",sqrts.c_str(),sqrts.c_str(),name.c_str(),0,8),MH,Form("%s/xs/%sTeV/%sTeV-%s.txt",location.c_str(),sqrts.c_str(),sqrts.c_str(),name.c_str()),0,8,1,"CSPLINE");
		  ws.import(splineM);
		  ws.import(splineP);

		 }
		RooSpline1D spline(Form("%s_%sTeV",name.c_str(),sqrts.c_str()),Form("file %sTeV/%sTeV-%s.txt, x=%d, y=%d",sqrts.c_str(),sqrts.c_str(),name.c_str(),0,1),MH,Form("%s/xs/%sTeV/%sTeV-%s.txt",location.c_str(),sqrts.c_str(),sqrts.c_str(),name.c_str()),0,1,1,"CSPLINE");		
		ws.import(spline);
	}

	// make the spline for qqZH 
	RooFormulaVar spline_qqZH(Form("qqZH_%sTeV",sqrts.c_str()),Form("qqZH (ZH-ggZH) - %sTeV",sqrts.c_str()),"@0-@1",RooArgList(*ws.function(Form("ZH_%sTeV",sqrts.c_str())),*ws.function(Form("ggZH_%sTeV",sqrts.c_str()))));
	ws.import(spline_qqZH);
	
	ws.Print();
	fout->cd();
	ws.Write();
	fout->Close();
}

void makeAllSMSplines(){

	gROOT->ProcessLine(".L $CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so");

	// Make XS 
	makeXSSplines("14");	
	makeXSSplines("13.6");	
	makeXSSplines("13");	
	makeXSSplines("8");	
	makeXSSplines("7");	

	// And BR
	makeBRSplines();
}

