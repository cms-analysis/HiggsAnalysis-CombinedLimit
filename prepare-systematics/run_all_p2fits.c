

#include "histio.c"

   TH1* get_hist( const char* hname ) ;

  //--------

   void run_all_p2fits( bool add_full_dev = true,
                        const char* input_file_16 = "../fit-input-root-files/2016/ttbar_systematics.root",
                        const char* input_file_17 = "../fit-input-root-files/2017/ttbar_systematics.root" ) {

      gDirectory -> Delete( "h*" ) ;

      gSystem -> Exec( "mkdir -p output-files" ) ;

      loadHist( input_file_16, "2016" ) ;
      loadHist( input_file_17, "2017" ) ;

      gDirectory -> ls( "*qcdCR*" ) ;

      TCanvas* can = new TCanvas( "can", "", 900, 700 ) ;


      TTree* tt_out = new TTree( "qcdCR_syst_parameters", "constants for the QCD CR systematic" ) ;
      int Njets ;
      int MVAbin ;
      int year ;
      double coef0 ;
      double coef1 ;
      double coef2 ;
      double coef3 ;
      tt_out -> Branch( "year", &year, "year/I" ) ;
      tt_out -> Branch( "Njets", &Njets, "Njets/I" ) ;
      tt_out -> Branch( "MVAbin", &MVAbin, "MVAbin/I" ) ;
      tt_out -> Branch( "coef0", &coef0, "coef0/D" ) ;
      tt_out -> Branch( "coef1", &coef1, "coef1/D" ) ;
      tt_out -> Branch( "coef2", &coef2, "coef2/D" ) ;
      tt_out -> Branch( "coef3", &coef3, "coef3/D" ) ;

      vector<TGraphErrors*> tge_vec ;

      for ( year=2016; year<=2017; year++ ) {
         for ( int di=1; di<=4; di++ ) {

            printf("\n\n ====== %d, D%d\n\n", year, di ) ;

            char hname[1000] ;

            sprintf( hname, "D%d_qcdCR_%d", di, year ) ;
            TH1* hp = get_hist( hname ) ;

            if ( add_full_dev ) {
               sprintf( hname, "D%d_qcdCRErr_%d", di, year ) ;
               TH1* hp_err = get_hist( hname ) ;
               printf("\n\n Adding full deviation errors:\n") ;
               for ( int bi=1; bi<=(hp_err->GetNbinsX()); bi++ ) {
                  printf("   bin %2d :   %7.4f +/- %7.4f , total error %7.4f\n", bi, hp->GetBinContent( bi ) , hp -> GetBinError( bi ) , hp_err -> GetBinContent( bi ) - 1. ) ;
                  hp -> SetBinError( bi, (hp_err -> GetBinContent( bi ) - 1.) ) ;
               } // bi
            } // add_full_dev?

            printf("\n Running fit.\n") ;

            TF1* tf1 = new TF1( "tf1", "[0]+[1]*x+[2]*x*x", 0., 6. ) ;
            tf1 -> SetParameter( 0, 1. ) ;
            tf1 -> SetParameter( 1, 0. ) ;
            tf1 -> SetParameter( 2, 0. ) ;

            TFitResultPtr tfr = hp -> Fit( tf1, "S" ) ;

            double p0 = tfr -> Parameter(0) ;
            double p1 = tfr -> Parameter(1) ;
            double p2 = tfr -> Parameter(2) ;
            printf("\n\n p0 = %7.3f , p1 = %9.6f , p2 = %10.8f\n\n", p0, p1, p2 ) ;

            TMatrixDSym cov = tfr -> GetCovarianceMatrix() ;
            printf(" Fit covariance matrix:\n") ;
            cov.Print() ;

           //--- Save the fit values and errors at each Njet point in the output file.
            double x_val[6] ;
            double f_val[6] ;
            double x_err[6] ;
            double f_err[6] ;

            int n_hist_bins = hp -> GetNbinsX() ;

            for ( int bi = 1; bi <= n_hist_bins ; bi++ ) {

               double x = hp -> GetXaxis() -> GetBinCenter( bi ) ;

               TMatrixD partial_derivative_row_vec(1,3) ;
               TMatrixD partial_derivative_col_vec(3,1) ;

               partial_derivative_row_vec[0][0] = 1 ;
               partial_derivative_row_vec[0][1] = x ;
               partial_derivative_row_vec[0][2] = x*x ;

               partial_derivative_col_vec[0][0] = 1 ;
               partial_derivative_col_vec[1][0] = x ;
               partial_derivative_col_vec[2][0] = x*x ;

               TMatrixD f_err2 = partial_derivative_row_vec * cov * partial_derivative_col_vec ;

               f_err2.Print() ;

               x_val[bi-1] = x ;
               x_err[bi-1] = 0.5 ;
               f_val[bi-1] = p0 + p1 * x + p2 * x * x ;
               f_err[bi-1] = sqrt( f_err2[0][0] ) ;

               printf( "  %2d :  x = %7.1f ,  f = %7.3f +/- %7.3f\n", bi, x, f_val[bi-1], f_err[bi-1] ) ;

            } // bi

            TGraphErrors* tge = new TGraphErrors( n_hist_bins, x_val, f_val, x_err, f_err ) ;
            char tgname[100] ;
            sprintf( tgname, "tge_%d_D%d_fit", year, di ) ;
            tge -> SetName( tgname ) ;
            tge -> SetFillColor(29) ;

            tge_vec.emplace_back( tge ) ;

            hp -> Draw() ;
            tge -> Draw( "3 same" ) ;
            hp -> Draw( "same" ) ;

            char fname[1000] ;
            sprintf( fname, "output-files/fit-plot-%d-D%d.pdf", year, di ) ;
            can -> SaveAs( fname ) ;




            TVectorD eigenVals ;
            TMatrixD eigenVecs = cov.EigenVectors( eigenVals ) ;
            printf("  Eigen Vals: %.8f, %.8f, %.8f\n", eigenVals(0), eigenVals(1), eigenVals(2) ) ;

            printf("  Eigen Vec matrix:\n" ) ;
            printf("     columns are eigen vectors.\n") ;
            for ( int i=0; i<3; i++ ) {
               for ( int j=0; j<3; j++ ) {
                  printf("   (%d,%d) = %7.3f | ", i, j, eigenVecs(i,j) ) ;
               } // j
               printf("\n") ;
            } // i

            TMatrixD eigenVecsT = eigenVecs ;
            eigenVecsT.T() ; // Note!  calling T() changes the contents!  Don't do this... TMatrixD eigenVecsT = eigenVecs.T()

            TMatrixD par_vec_col(3,1) ;
            par_vec_col[0][0] = p0 ;
            par_vec_col[1][0] = p1 ;
            par_vec_col[2][0] = p2 ;

            TMatrixD par_prime_vec_col = eigenVecsT * par_vec_col ;
            printf("\n\n par_prime_vec_col:\n") ;
            par_prime_vec_col.Print() ;

            for ( int nji=1; nji <= (hp -> GetNbinsX()) ; nji++ ) {

               double x = hp -> GetXaxis() -> GetBinCenter( nji ) ;

               Njets = nji ;
               MVAbin = di ;

               coef0 = 0. ;

               double c1 = eigenVecs(0,0) + eigenVecs(1,0) * x + eigenVecs(2,0) * x * x ;
               double c2 = eigenVecs(0,1) + eigenVecs(1,1) * x + eigenVecs(2,1) * x * x ;
               double c3 = eigenVecs(0,2) + eigenVecs(1,2) * x + eigenVecs(2,2) * x * x ;

               coef0 = c1 * par_prime_vec_col[0][0]  + c2 * par_prime_vec_col[1][0]  + c3 * par_prime_vec_col[2][0] ;

               coef1 = c1 * sqrt( eigenVals(0) ) ;
               coef2 = c2 * sqrt( eigenVals(1) ) ;
               coef3 = c3 * sqrt( eigenVals(2) ) ;

               printf("   %d , D%d , Nj%d , x = %.1f ,  coef0 = %10.6f, coef1 = %10.6f, coef2 = %10.6f, coef3 = %10.6f\n",
                   year, di, nji, x, coef0, coef1, coef2, coef3 ) ;

               tt_out  -> Fill() ;

            } // nji


         } // di
      } // year

      saveHist( "output-files/qcdcr-syst-parameters.root", "*" ) ;

      TFile of( "output-files/qcdcr-syst-parameters.root", "update" ) ;
      for ( int i=0; i<tge_vec.size(); i++ ) {
         tge_vec.at(i)->Write() ;
      }
      of.Close() ;


   } // run_all_p2fits

  //==============================================================================================================

   TH1* get_hist( const char* hname ) {
      TH1* hp = (TH1*) gDirectory -> FindObject( hname ) ;
      if ( hp == 0x0 ) {
         printf("\n\n *** can't find hist with name  %s\n\n", hname ) ;
         gSystem -> Exit(-1) ;
      }
      return hp ;
   }

  //==============================================================================================================


