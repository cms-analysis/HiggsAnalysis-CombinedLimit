void makeRooMultiPdfWorkspace(){

   // Load the combine Library 
   gSystem->Load("libHiggsAnalysisCombinedLimit.so");

   // Open the dummy H->gg workspace 
   TFile *f_hgg = TFile::Open("toyhgg_in.root");
   RooWorkspace *w_hgg = (RooWorkspace*)f_hgg->Get("multipdf");

   // The observable (CMS_hgg_mass in the workspace)
   RooRealVar *mass =  w_hgg->var("CMS_hgg_mass");

   // Get three of the functions inside, exponential, linear polynomial, power law
   RooAbsPdf *pdf_exp = w_hgg->pdf("env_pdf_1_8TeV_exp1");
   RooAbsPdf *pdf_pol = w_hgg->pdf("env_pdf_1_8TeV_bern2");
   RooAbsPdf *pdf_pow = w_hgg->pdf("env_pdf_1_8TeV_pow1");


   // Fit the functions to the data to set the "prefit" state (note this can and should be redone with combine when doing 
   // bias studies as one typically throws toys from the "best-fit"
   RooDataSet *data = (RooDataSet*)w_hgg->data("roohist_data_mass_cat1_toy1_cutrange__CMS_hgg_mass");
   pdf_exp->fitTo(*data);  // index 0
   pdf_pow->fitTo(*data); // index 1 
   pdf_pol->fitTo(*data);   // index 2

   // Make a plot (data is a toy dataset)
   RooPlot *plot = mass->frame();   data->plotOn(plot);
   pdf_exp->plotOn(plot,RooFit::LineColor(kGreen));
   pdf_pol->plotOn(plot,RooFit::LineColor(kBlue));
   pdf_pow->plotOn(plot,RooFit::LineColor(kRed));
   plot->SetTitle("PDF fits to toy data");
   plot->Draw();

   // Make a RooCategory object. This will control which of the pdfs is "active"
   RooCategory cat("pdf_index","Index of Pdf which is active");

   // Make a RooMultiPdf object. The order of the pdfs will be the order of their index, ie for below 
   // 0 == exponential
   // 1 == linear function
   // 2 == powerlaw
   RooArgList mypdfs;
   mypdfs.add(*pdf_exp);
   mypdfs.add(*pdf_pol);
   mypdfs.add(*pdf_pow);
   
   RooMultiPdf multipdf("roomultipdf","All Pdfs",cat,mypdfs);
   
   // As usual make an extended term for the background with _norm for freely floating yield
   RooRealVar norm("roomultipdf_norm","Number of background events",0,10000);
   
   // Save to a new workspace
   TFile *fout = new TFile("background_pdfs.root","RECREATE");
   RooWorkspace wout("backgrounds","backgrounds");
   wout.import(cat);
   wout.import(norm);
   wout.import(multipdf);
   wout.Print();
   wout.Write();

}