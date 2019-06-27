
#include "histio.c"


   void toy_validation( const char* input_file = "output-files/qcdcr-syst-parameters.root" ) {

      char outfile[1000] ;
      gStyle -> SetOptStat( "mr" ) ;

      TFile infile( input_file, "read" ) ;
      if ( !infile.IsOpen() ) { printf("\n\n *** bad input file: %s\n\n", input_file ) ; gSystem -> Exit(-1) ; }

      TTree* tt = (TTree*) infile.Get( "qcdCR_syst_parameters" ) ;
      if ( tt == 0x0 ) { printf("\n\n *** can't find qcdCR_syst_parameters TTree in %s\n\n", input_file ) ; gSystem -> Exit(-1) ; }

      TCanvas* can = new TCanvas( "can", "", 900, 700 ) ;
      TCanvas* can2 = new TCanvas( "can2", "", 900, 700 ) ;

      int Njets ;
      int MVAbin ;
      int year ;
      double coef0 ;
      double coef1 ;
      double coef2 ;
      double coef3 ;
      tt -> SetBranchAddress( "Njets", &Njets ) ;
      tt -> SetBranchAddress( "MVAbin", &MVAbin ) ;
      tt -> SetBranchAddress( "year", &year ) ;
      tt -> SetBranchAddress( "coef0", &coef0 ) ;
      tt -> SetBranchAddress( "coef1", &coef1 ) ;
      tt -> SetBranchAddress( "coef2", &coef2 ) ;
      tt -> SetBranchAddress( "coef3", &coef3 ) ;

      TRandom tran(12345) ;

      double c0[2][4][6] ;
      double c1[2][4][6] ;
      double c2[2][4][6] ;
      double c3[2][4][6] ;
      for ( int ei = 0; ei < tt -> GetEntries() ; ei++ ) {

         tt -> GetEntry( ei ) ;
         printf("  %3d :  %d, D%d, Njets = %d :  coef0 = %7.4f  coef1 = %7.4f  coef2 = %7.4f  coef3 = %7.4f\n", ei, year, MVAbin, Njets, coef0, coef1, coef2, coef3 ) ;
         int yi=year-2016 ;
         int di=MVAbin-1 ;
         int nji = Njets-1 ;
         c0[yi][di][nji] = coef0 ;
         c1[yi][di][nji] = coef1 ;
         c2[yi][di][nji] = coef2 ;
         c3[yi][di][nji] = coef3 ;

      } // ei

      for ( int yi=0; yi<2; yi++ ) {
         for ( int di=0; di<4; di++ ) {
            char hname[100] ;
            int nbx(6) ;
            int nby(200) ;
            sprintf( hname, "h_toy_%d_D%d", yi+2016, di+1 ) ;
            TH2F* hp = new TH2F( hname, hname, nbx, 0., 6., nby, 0., 2. ) ;
            int ntoy = 10000 ;
            for ( int ti=0; ti<ntoy; ti++ ) {
               double theta1 = tran.Gaus(0.,1.) ;
               double theta2 = tran.Gaus(0.,1.) ;
               double theta3 = tran.Gaus(0.,1.) ;
               for ( int xi=0; xi<6; xi++ ) {
                  double x = xi+0.5 ;
                  double fval = c0[yi][di][xi] + theta1 * c1[yi][di][xi] + theta2 * c2[yi][di][xi] + theta3 * c3[yi][di][xi] ;
                  hp -> Fill( x, fval ) ;
               } // xi
            } // ti

            char tgename[1000] ;
            sprintf( tgename, "tge_%d_D%d_fit", yi+2016, di+1 ) ;
            TGraphErrors* tge = (TGraphErrors*) infile.Get( tgename ) ;
            if ( tge == 0x0 ) { printf("\n\n *** can't find %s\n\n", tgename ) ; gSystem -> Exit(-1) ; }

            can -> cd() ;
            hp -> Draw("colz") ;
            tge -> SetLineColor(1) ;
            tge -> SetLineWidth(2) ;
            tge -> Draw( "same" ) ;
            tge -> Draw( "|| same" ) ;
            sprintf( outfile, "output-files/validation-plot1-%d-D%d.pdf", yi+2016, di+1 ) ;
            can -> SaveAs( outfile ) ;

            double* fvals = tge -> GetY() ;
            double* ferrs = tge -> GetEY() ;

            can2 -> cd() ;
            can2 -> Clear() ;
            can2 -> Divide(3,2) ;

            TText* ttlabel = new TText() ;

            for ( int nji=0; nji<6; nji++ ) {
               can2 -> cd( nji+1 ) ;
               char pname[1000] ;
               sprintf( pname, "h_%d_D%d_Njets%d", yi+2016, di+1, nji+1 ) ;
               TH1* hproj = hp -> ProjectionY( pname, nji+1, nji+1 ) ;
               hproj -> Draw() ;
               char label[100] ;
               sprintf( label, "Fit = %7.4f +/- %7.4f", fvals[nji], ferrs[nji] ) ;
               ttlabel -> DrawTextNDC( 0.15, 0.80, label ) ;
            } // nji

            sprintf( outfile, "output-files/validation-plot2-%d-D%d.pdf", yi+2016, di+1 ) ;
            can2 -> SaveAs( outfile ) ;



         } // di
      } // yi

      saveHist( "output-files/toy-validation.root", "h*" ) ;


   } // toy_validation


