TLine *getMedian(TH1F *b){
	double medx[1] ={0.5};
	double medy[1] ={0.};

	b->GetQuantiles(1,medy,medx);
	std::cout << " Median - bonly " << medy[0] << std::endl;
	int binc = b->FindBin(medy[0]);
	TLine *l = new TLine(medy[0],b->GetMinimum(),medy[0],b->GetBinContent(binc));
	l->SetLineColor(1);
	l->SetLineStyle(2);
	return l;
}

double getMedianVal(TH1F *b){
	double medx[1] ={0.5};
	double medy[1] ={0.};

	b->GetQuantiles(1,medy,medx);
	return medy[0];
}

TH1F *getQuantHist(TH1F *b, double qMin, double qMax){
	
	TH1F *hist68 = (TH1F*) b->Clone();
	hist68->SetName(Form("%s_%g_%g",b->GetName(),qMin,qMax));

	double medx[2] ={qMin,qMax};
	double medy[2] ={0.,0.};

	b->GetQuantiles(2,medy,medx);
	std::cout << Form(" Quantiles, %g,%g - bonly ",medx[0],medx[1]) << medy[0] << ", " << medy[1] << std::endl;

	for (int bin=1;bin<=b->GetNbinsX();bin++){
		double binc = b->GetBinCenter(bin);
		if (binc < medy[0] || binc > medy[1] ) {
			hist68->SetBinContent(bin,0) ;
		}
	}
	
	hist68->SetFillStyle(3001);
	return hist68;
}


TPaveText *cmsprel;
void spam(const char *text=0, double x1=0.17, double y1=0.89, double x2=0.58, double y2=0.94, int textAlign=12, bool fill=true, float fontSize=0) {
   if (TString(text).Contains("#splitline")) {
       if (y1 > 0.5) y1 -= 0.06; else y2 += 0.06;
   }
   TPaveText *cmsprel = new TPaveText(x1,y1,x2,y2,"brtlNDC");
   if (fontSize == 0) fontSize = 0.04;
   cmsprel->SetTextSize(fontSize);
   cmsprel->SetFillColor(0);
   cmsprel->SetFillStyle((fill || text == 0) ? 1001 : 0);
   cmsprel->SetLineStyle(0);
   cmsprel->SetLineColor(0);
   cmsprel->SetLineWidth(0);
   cmsprel->SetTextAlign(textAlign);
   cmsprel->SetTextFont(42);
   cmsprel->SetTextColor( 1 );
   cmsprel->AddText(text);
   cmsprel->SetBorderSize(0);
   cmsprel->Draw("same");
}


TH1F *tail(TH1F *dist, double cut,int mode) {
    TH1F *ret = (TH1F*) dist->Clone();
    if (mode==0){
      for (int b = 1, bmax = dist->GetXaxis()->FindBin(cut); b <= bmax; ++b) {
        ret->SetBinContent(b, 0);
        ret->SetBinError(b, 0);
      }
    } else {
      for (int b = dist->GetXaxis()->FindBin(cut); b <= dist->GetNbinsX(); ++b) {
        ret->SetBinContent(b, 0);
        ret->SetBinError(b, 0);
      }
    }
    return ret;
}

double tailReal(TTree *t, std::string br, double cv, int mode){

    int tpass;
    int ttotal;
    if (br=="qB") {
    	if (mode==0) tpass = t->Draw("max(2*q,0)>>qTMP",Form("weight*(type==-1 && max(2*q,0)> %g)",cv));
        else tpass = t->Draw("2*q>>qTMP",Form("weight*(type==-1 && 2*q < %g)",cv));
	ttotal = t->Draw("2*q>>whatever","weight*(type==-1)");
    } else if (br=="qS"){
        if (mode==0) tpass = t->Draw("max(2*q,0)>>qTMP",Form("weight*(type==+1 && max(2*q,0)> %g)",cv),"");
        else tpass  =  t->Draw("2*q>>qTMP",Form("weight*(type==+1 && 2*q < %g)",cv),"");
	ttotal = t->Draw("2*q>>whatever","weight*(type==+1)");
    }
    std::cout << "tpass = " << tpass <<", ttotal = " << ttotal<< std::endl;
    std::cout << "frac pass = " <<  ((double)tpass)/ttotal << std::endl;
    return ((double)tpass)/ttotal;
}



TCanvas *q0Plot(float mass, std::string poinam , int rebin=0) {

    if (gFile == 0) { std::cerr << "You must have a file open " << std::endl; return 0; }
    TTree *t = (TTree *) gFile->Get("q");
    if (t == 0) { std::cerr << "File " << gFile->GetName() << " does not contain a tree called 'q'" << std::endl; return 0; }

    TCanvas *c1 = new TCanvas("c1","c1");
    c1->SetBottomMargin(0.15);
    
    TH1F *qB;
    TH1F *qS;

    t->Draw("max(2*q,0)>>qB","weight*(type==-1)");
    qB = (TH1F*) gROOT->FindObject("qB")->Clone();
    qB->SetName("NullHyp");
    qB->Print();

    double yMin = 4.9/qB->Integral();

    t->Draw("max(2*q,0)>>qObs","weight*(type==0)");
    double qObs = ((TH1F*) gROOT->FindObject("qObs"))->GetMean();

    TArrow *qO = new TArrow(qObs, 0.2, qObs, yMin*1.05, 0.01, "---|>");
    
    if (rebin>0){
	qB->Rebin(rebin);
    }

    double  nB = qB->Integral();
    qB->Scale(1.0/qB->Integral());

    gStyle->SetOptStat(0); c1->SetLogy(1);
    qB->SetLineColor(kRed);
    qB->SetLineWidth(2);
    qO->SetLineColor(kBlack);
    qO->SetLineWidth(3);

    double clB;
    clB  = tailReal(t,"qB",qObs,0);

    double clBerr  = sqrt(clB*(1-clB)/nB);
    double sig = RooStats::PValueToSignificance(clB);
    double sigerr = 0.5*( TMath::Abs(RooStats::PValueToSignificance(clB+clBerr)-sig) + TMath::Abs(RooStats::PValueToSignificance(clB-clBerr)-sig));

    printf("P-val (CLb)  = %.4f +/- %.4f\n", clB , clBerr);
    printf("Signif  = %.1f +/- %.2f sigma\n",  sig, sigerr);

    // Worst way to calculate !
    TH1F *qB1 = tail(qB, qObs,0);
    qB1->SetFillColor(qB1->GetLineColor()); qB1->SetLineWidth(0); qB1->SetFillStyle(3345);

    TLegend *leg1 = new TLegend(.45,.68,.93,.85);
    leg1->SetFillColor(0);
    leg1->SetShadowColor(0);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.04);
    leg1->SetLineColor(1);

    leg1->AddEntry(qB, "expected for Null", "L");
    leg1->AddEntry(qO, "observed value", "L");

    TLegend *leg2 = new TLegend(.58,.67,.93,.54);
    leg2->SetFillColor(0);
    leg2->SetShadowColor(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.04);
    leg2->SetLineColor(1);
    leg2->AddEntry(qB1, Form("CL_{b}   = %.4f (%.1f#sigma)", clB,sig), "F");
    qB->Draw();
    qB1->Draw("HIST SAME"); 
    qO->Draw(); 
    qB->Draw("AXIS SAME");
    qB->GetYaxis()->SetRangeUser(yMin, 2.0);
    leg1->Draw();
    leg2->Draw();
    qB->SetTitle("");
    qB->GetYaxis()->SetTitle("");
    qB->GetXaxis()->SetTitle(Form("q_{0}(m_{H} = %g GeV)",mass));
    qB->GetXaxis()->SetTitleOffset(1.05);
    
    c1->SetName(Form("q0_%.1f",mass));

    return c1;
}



TCanvas *qmuPlot(float mass, std::string poinam, double poival, int mode=0, int invert=0,int rebin=0, int runExpected_=0, double quantileExpected_=0.5) {
    if (gFile == 0) { std::cerr << "You must have a file open " << std::endl; return 0; }
    TTree *t = (TTree *) gFile->Get("q");
    if (t == 0) { std::cerr << "File " << gFile->GetName() << " does not contain a tree called 'q'" << std::endl; return 0; }

    TCanvas *c1 = new TCanvas("c1","c1");
    c1->SetBottomMargin(0.15);
    
    TH1F *qB;
    TH1F *qS;

    if (mode==0) t->Draw("max(2*q,0)>>qB","weight*(type==-1)");
    else t->Draw("2*q>>qB","weight*(type==-1)");

    if (invert) { // swap null for alt!
       qS = (TH1F*) gROOT->FindObject("qB")->Clone();
       qS->SetName("NullHyp");
       qS->Print();  // in root6, without this Print, the integral is 0?!?!?!
    } else { 
        if (gROOT->FindObject("qB") == nullptr)
            qB=new TH1F("qB","empty",10,0,10);
        else
            qB = (TH1F*) gROOT->FindObject("qB")->Clone();
        qB->SetName("NullHyp");
        qB->Print();
    }

    if (mode==0) t->Draw("max(2*q,0)>>qS","weight*(type==+1)","SAME");
    else t->Draw("2*q>>qS","weight*(type==+1)","SAME");

    if (invert) { // swap null for alt!
       qB = (TH1F*) gROOT->FindObject("qS")->Clone(); 
       qB->SetName("AltHyp");
       qB->Print();
    } else { 
        if (gROOT->FindObject("qS") == nullptr)
            qS=new TH1F("qS","empty",10,0,10);
        else
            qS = (TH1F*) gROOT->FindObject("qS")->Clone();
        qS->SetName("AltHyp");
        qS->Print();
    }
    double yMin = 4.9/qB->Integral();

    if (mode==0)t->Draw("max(2*q,0)>>qObs","weight*(type==0)");
    else t->Draw("2*q>>qObs","weight*(type==0)");
    double qObs = (gROOT->FindObject("qObs")==nullptr) ? -1: ((TH1F*) gROOT->FindObject("qObs"))->GetMean();
    if (runExpected_) {
      
      double medx[1] ={1.-quantileExpected_};
      double medy[1] ={0.};
      qB->GetQuantiles(1,medy,medx);
      qObs = medy[0];
    }

    TArrow *qO = new TArrow(qObs, 0.2, qObs, yMin*1.05, 0.01, "---|>");
    
    if (rebin>0){
	qB->Rebin(rebin);
	qS->Rebin(rebin);
    }

    double nS = qS->Integral(), nB = qB->Integral();
    qB->Scale(1.0/qB->Integral());
    qS->Scale(1.0/qS->Integral());

    gStyle->SetOptStat(0); c1->SetLogy(1);
    qS->SetLineColor(kBlue);
    qS->SetLineWidth(2);
    qB->SetLineColor(kRed);
    qB->SetLineWidth(2);
    qO->SetLineColor(kBlack);
    qO->SetLineWidth(3);

    double clSB;
    double clB;
    if (invert){ 
    	clSB = tailReal(t,"qB",qObs,mode);
    	clB  = tailReal(t,"qS",qObs,mode);
    } else {
    	clSB = tailReal(t,"qS",qObs,mode);
    	clB  = tailReal(t,"qB",qObs,mode);
    }

    //double clSB = qS1->Integral(), clB = qB1->Integral(), 
    double clS = clSB/clB;
    double clSBerr = sqrt(clSB*(1-clSB)/nS);
    double clBerr  = sqrt(clB*(1-clB)/nB);
    double clSerr  = clS * TMath::Hypot(clBerr/clB, clSBerr/clSB);
    printf("CLs+b = %.4f +/- %.4f\n", clSB, clSBerr);
    printf("CLb   = %.4f +/- %.4f\n", clB , clBerr);
    printf("CLs   = %.4f +/- %.4f\n", clS , clSerr);

    // Worst way to calculate !
    TH1F *qS1 = tail(qS, qObs,mode); 
    TH1F *qB1 = tail(qB, qObs,mode);
    qS1->SetFillColor(64); qS1->SetLineWidth(0); qS1->SetFillStyle(1001);
    qB1->SetFillColor(qB1->GetLineColor()); qB1->SetLineWidth(0); qB1->SetFillStyle(3345);

    TLegend *leg1 = new TLegend(.45,.68,.93,.85);
    leg1->SetFillColor(0);
    leg1->SetShadowColor(0);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.04);
    leg1->SetLineColor(1);
    if (mode==0){
     leg1->AddEntry(qS, "expected for sig+bkg",  "L");
     leg1->AddEntry(qB, "expected for bkg-only", "L");
    } else {
     if (invert){ 
         leg1->AddEntry(qS, "expected for Alt",  "L");
    	 leg1->AddEntry(qB, "expected for Null", "L");
     } else {
         leg1->AddEntry(qS, "expected for Null",  "L");
    	 leg1->AddEntry(qB, "expected for Alt", "L");
     }
    }
    leg1->AddEntry(qO, "observed value", "L");

    TLegend *leg2 = new TLegend(.63,.67,.93,.48);
    leg2->SetFillColor(0);
    leg2->SetShadowColor(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.04);
    leg2->SetLineColor(1);
    //if (mode==0) 
    leg2->AddEntry(qS1, Form("CL_{s+b} = %.4f", clSB), "F"); 
    //if (mode==0) 
    leg2->AddEntry(qB1, Form("CL_{b}   = %.4f", clB), "F");
    leg2->AddEntry("",  Form("CL_{s}   = %.4f", clS), "");

    qB->Draw();

    // Draw bands on the Background distribution 
    TLine *median_b = getMedian(qB); 
    TH1F *oneSig    = getQuantHist(qB,0.16,0.84);   oneSig->SetFillColor(kGreen+1);
    TH1F *twoSig    = getQuantHist(qB,0.025,0.975); twoSig->SetFillColor(kYellow);
    twoSig->Draw("histFsame");
    oneSig->Draw("histFsame");
    //qB->Draw("same");

    qS->Draw("SAME"); 
    qS1->Draw("HIST SAME"); 
    qB1->Draw("HIST SAME"); 
    qO->Draw(); 
    qB->Draw("AXIS SAME");
    qB->GetYaxis()->SetRangeUser(yMin, 2.0);
    median_b->Draw();
    leg1->Draw();
    leg2->Draw();
    qB->SetTitle("");
    qB->GetYaxis()->SetTitle("");
    qB->GetXaxis()->SetTitle(Form("q_{%s}(%s = %g, m_{H} = %g GeV)",poinam.c_str(),poinam.c_str(),poival,mass));
    qB->GetXaxis()->SetTitleOffset(1.05);


    //c1->Print(Form("qmu_example_%.1f.pdf",mass));
    //c1->Print(Form("qmu_example_%.1f.png",mass));
    //c1->Print(Form("qmu_example_%.1f.eps",mass));

    c1->SetName(Form("qmu_%.1f_%g",mass,poival));

    return c1;

}
