{
   
   // library containing the RooSplineND 
   gSystem->Load("libHiggsAnalysisCombinedLimit.so");

   TTree *tree = new TTree("tree_vals","tree_vals");
   float xb,yb,fb;

   tree->Branch("f",&fb,"f/Float_t");
   tree->Branch("x",&xb,"x/Float_t");
   tree->Branch("y",&yb,"y/Float_t");
   
   TRandom3 *r = new TRandom3();
   int nentries = 20; // just use a regular grid of 20x20

   double xmin = -3.2;
   double xmax = 3.2;
   double ymin = -3.2;
   double ymax = 3.2;

   for (int n=0;n<nentries;n++){
    for (int k=0;k<nentries;k++){

      xb=xmin+n*((xmax-xmin)/nentries);
      yb=ymin+k*((ymax-ymin)/nentries);
      // Gaussian * cosine function radial in "F(x^2+y^2)"
      double R = (xb*xb)+(yb*yb);
      fb = 0.1*TMath::Exp(-1*(R)/9)*TMath::Cos(2.5*TMath::Sqrt(R));
      tree->Fill();
     }
   }
   
   // 2D graph of points in tree
   TGraph2D *p0 = new TGraph2D();
   p0->SetMarkerSize(0.8);
   p0->SetMarkerStyle(20);
   
   int c0=0;
   for (int p=0;p<tree->GetEntries();p++){
        tree->GetEntry(p);
        p0->SetPoint(c0,xb,yb,fb);
        c0++;
        }


   // ------------------------------ THIS IS WHERE WE BUILD THE SPLINE ------------------------ //
   // Create 2 Real-vars, one for each of the parameters of the spline 
   // The variables MUST be named the same as the corresponding branches in the tree
   RooRealVar x("x","x",0.1,xmin,xmax); 
   RooRealVar y("y","y",0.1,ymin,ymax);

   
   // And the spline - arguments are 
   // Required ->   name, title, arglist of dependants, input tree, 
   // Optional ->  function branch name, interpolation width (tunable parameter), rescale Axis bool, cutstring 
   // The tunable parameter gives the radial basis a "width", over which the interpolation will be effectively taken 
   
   // the reascale Axis bool (if true) will first try to rescale the points so that they are of order 1 in range
   // This can be helpful if for example one dimension is in much larger units than another.
   
   // The cutstring is just a ROOT string which can be used to apply cuts to the tree in case only a sub-set of the points should be used 
   
   RooArgList args(x,y);
   RooSplineND *spline = new RooSplineND("spline","spline",args,tree,"f",1,true);
      // ----------------------------------------------------------------------------------------- //
   

   //TGraph *gr = spline->getGraph("x",0.1); // Return 1D graph. Will be a slice of the spline for fixed y generated at steps of 0.1
   
   // Plot the 2D spline 
   TGraph2D *gr = new TGraph2D();
   int pt = 0;
   for (double xx=xmin;xx<xmax;xx+=0.1){
     for (double yy=xmin;yy<ymax;yy+=0.1){
        x.setVal(xx);
        y.setVal(yy);
      gr->SetPoint(pt,xx,yy,spline.getVal());
      pt++;
     }
   }

   gr->SetTitle("");

   gr->SetLineColor(1);
   //p0->SetTitle("0.1 exp(-(x{^2}+y{^2})/9) #times Cos(2.5#sqrt{x^{2}+y^{2}})");
   gr->Draw("surf");
   gr->GetXaxis()->SetTitle("x");
   gr->GetYaxis()->SetTitle("y");
   p0->Draw("Pcolsame");

   //p0->Draw("surfsame");
   TLegend *leg = new TLegend(0.2,0.82,0.82,0.98);
   leg->SetFillColor(0);
   leg->AddEntry(p0,"0.1 exp(-(x{^2}+y{^2})/9) #times Cos(2.5#sqrt{x^{2}+y^{2}})","p");
   leg->AddEntry(gr,"RooSplineND (N=2) interpolation","L");
   leg->Draw();
}