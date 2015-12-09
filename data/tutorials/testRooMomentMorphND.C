#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooHistPdf.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMomentMorphND.h"

#include "TCanvas.h"
#include "TSystem.h"

//
void testRooMomentMorphND(RooMomentMorphND::Setting morphSetting,Bool_t useHorizontalMorphing)
{
  using namespace RooFit;

  gSystem->Load("libHiggsAnalysisCombinedLimit");

  RooWorkspace w("w", 1);

  w.factory("obs[0,400]");

  //smeared gaussian distributions
  w.factory("RooGaussian::gaussian00(obs,156.8,15)");
  w.factory("RooGaussian::gaussian01(obs,160.0,15)");
  w.factory("RooGaussian::gaussian02(obs,163.3,15)");

  w.factory("RooGaussian::gaussian10(obs,169.05,16)");
  w.factory("RooGaussian::gaussian11(obs,172.50,16)");
  w.factory("RooGaussian::gaussian12(obs,175.95,16)");

  w.factory("RooGaussian::gaussian20(obs,181.3,17)");
  w.factory("RooGaussian::gaussian21(obs,185.0,17)");
  w.factory("RooGaussian::gaussian22(obs,188.7,17)");


  RooBinning firstDim(2,0,2);   //#bins = #nodes -1
  RooBinning secondDim(2,0,2);
  RooMomentMorphND::Grid refGrid(firstDim,secondDim);
  refGrid.addPdf( *(w.pdf("gaussian00")),0,0);
  refGrid.addPdf( *(w.pdf("gaussian01")),0,1);
  refGrid.addPdf( *(w.pdf("gaussian02")),0,2);
  refGrid.addPdf( *(w.pdf("gaussian10")),1,0);
  refGrid.addPdf( *(w.pdf("gaussian11")),1,1);
  refGrid.addPdf( *(w.pdf("gaussian12")),1,2);
  refGrid.addPdf( *(w.pdf("gaussian20")),2,0);
  refGrid.addPdf( *(w.pdf("gaussian21")),2,1);
  refGrid.addPdf( *(w.pdf("gaussian22")),2,2);

  w.factory("alpha[0,2]");
  w.factory("beta[0,2]");
  RooMomentMorphND pdf("morphndpdf","morphndpdf",
		       RooArgList( *(w.var("alpha")), *(w.var("beta")) ),
		       RooArgList( *(w.var("obs")) ),
		       refGrid,
		       morphSetting);
  pdf.useHorizontalMorphing(useHorizontalMorphing);

  //morph the template
  w.import(pdf);

  //for testing purposes
  w.factory(Form("RooGaussian::gaussiantest0v5_0(obs,%f,%f)",(156.8+169.05)/2., (15+16)/2.));
  w.factory(Form("RooGaussian::gaussiantest0v5_1(obs,%f,%f)",(160.0+172.50)/2., (15+16)/2.));

  //test the pdf
  TCanvas *c= new TCanvas("c","c",1000,1000);
  c->Divide(2,2);

  c->cd(1);
  RooPlot* frame1 = w.var("obs")->frame();
  w.var("alpha")->setVal(0);
  w.var("beta")->setVal(0);
  w.pdf("gaussian00")->plotOn(frame1, LineColor(kGray), LineStyle(kSolid));
  w.pdf("morphndpdf")->plotOn(frame1, LineColor(kRed), LineStyle(kDashed));
  frame1->SetTitle("Closure on (#alpha,#beta)=(0,0) (grid point)");
  frame1->Draw();

  c->cd(2);
  RooPlot* frame2 = w.var("obs")->frame();
  w.var("alpha")->setVal(1);
  w.var("beta")->setVal(2.0);
  w.pdf("gaussian12")->plotOn(frame2, LineColor(kGray), LineStyle(kSolid));
  w.pdf("morphndpdf")->plotOn(frame2, LineColor(kRed), LineStyle(kDashed));
  frame2->SetTitle("Closure on (#alpha,#beta)=(1,2) (grid point)");
  frame2->Draw();

  c->cd(3);
  RooPlot* frame3 = w.var("obs")->frame();
  w.var("alpha")->setVal(0.5);
  w.var("beta")->setVal(1.0);
  w.pdf("gaussiantest0v5_1")->plotOn(frame3, LineColor(kGray), LineStyle(kSolid));
  w.pdf("morphndpdf")->plotOn(frame3, LineColor(kRed), LineStyle(kDashed));
  frame3->SetTitle("Interpolation to (#alpha,#beta)=(0.5,1.0)");
  frame3->Draw();

  c->cd(4);
  RooPlot* frame4 = w.var("obs")->frame();
  w.var("alpha")->setVal(0.5);
  w.var("beta")->setVal(0.0);
  w.pdf("gaussiantest0v5_0")->plotOn(frame4, LineColor(kGray), LineStyle(kSolid));
  w.pdf("morphndpdf")->plotOn(frame4, LineColor(kRed), LineStyle(kDashed));
  frame4->SetTitle("Interpolation to (#alpha,#beta)=(0.5,0.0)");
  frame4->Draw();
}
