
#include "histio.c"

   void p2fit1( int dbin_ind = 1, bool add_full_dev = false, const char* infile = "../fit-input-root-files/2016/ttbar_systematics.root" ) {

      loadHist( infile ) ;

      char hname[100] ;
      sprintf( hname, "D%d_qcdCR", dbin_ind ) ;
      TH1* hp = (TH1*) gDirectory -> FindObject( hname ) ;
      if ( hp == 0x0 ) { printf("\n\n *** can't find hist %s in %s\n\n", hname, infile ) ; gSystem -> Exit(-1) ; }


      if ( add_full_dev ) {
         sprintf( hname, "D%d_qcdCRErr", dbin_ind ) ;
         TH1* hp_err = (TH1*) gDirectory -> FindObject( hname ) ;
         if ( hp_err == 0x0 ) { printf("\n\n *** can't find hist %s in %s\n\n", hname, infile ) ; gSystem -> Exit(-1) ; }
         printf("\n\n Adding full deviation errors:\n") ;
         for ( int bi=1; bi<=(hp_err->GetNbinsX()); bi++ ) {
            printf("   bin %2d :   %7.4f +/- %7.4f , full error %7.4f\n", bi, hp->GetBinContent( bi ) , hp -> GetBinError( bi ) , hp_err -> GetBinContent( bi ) - 1. ) ;
            hp -> SetBinError( bi, (hp_err -> GetBinContent( bi ) - 1.) ) ;
         }
      }


      TCanvas* can = new TCanvas( "can", "", 700, 900 ) ;
      can -> Divide(1,2) ;

      can -> cd(1) ;

      hp -> DrawCopy() ;
      gPad -> SetGridy(1) ;

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
      cov.Print() ;

      int n_hist_bins = hp -> GetNbinsX() ;
      double xmin = hp -> GetXaxis() -> GetBinLowEdge(1) ;
      double xmax = hp -> GetXaxis() -> GetBinUpEdge( n_hist_bins ) ;

      double func_err[6] ;

      double x_val[6] ;
      double f_val[6] ;
      double x_err[6] ;
      double f_err[6] ;

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
      tge -> SetFillColor(11) ;


      tge -> Draw("3 same" ) ;
      tge -> Draw("same" ) ;
      hp -> Draw("same") ;

     //-----
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

      TMatrixD eigenVec1_col(3,1) ;
      eigenVec1_col[0][0] = eigenVecs[0][0] ;
      eigenVec1_col[1][0] = eigenVecs[1][0] ;
      eigenVec1_col[2][0] = eigenVecs[2][0] ;

      TMatrixD check_cov_ev1 = cov * eigenVec1_col ;
      printf("  Check:  cov * ev1 = \n" ) ;
      check_cov_ev1.Print() ;

      TMatrixD eigenVec1_lambda1_col(3,1) ;
      eigenVec1_lambda1_col = eigenVals(0) * eigenVec1_col ;
      printf("  Check:  lambda1 * ev1 = \n" ) ;
      eigenVec1_lambda1_col.Print() ;

      printf("\n\n eigenVecs\n") ;
      eigenVecs.Print() ;

      TMatrixD eigenVecsT = eigenVecs ;
      eigenVecsT.T() ;

      TMatrixD check_unit_from_eigenvecs5 = eigenVecsT*eigenVecs ;
      printf("\n\n eigenVecsT*eigenVecs\n") ;
      check_unit_from_eigenvecs5.Print() ;

      TMatrixD diagonalCov = eigenVecsT * cov * eigenVecs ;
      printf("\n\n Diagonal cov:\n") ;
      diagonalCov.Print() ;

     //---------------

      TMatrixD par_vec_col(3,1) ;
      par_vec_col[0][0] = p0 ;
      par_vec_col[1][0] = p1 ;
      par_vec_col[2][0] = p2 ;

      TMatrixD par_prime_vec_col = eigenVecsT * par_vec_col ;
      printf("\n\n par_prime_vec_col:\n") ;
      par_prime_vec_col.Print() ;

      TMatrixD par_vec_from_par_prime_check(3,1) ;
      par_vec_from_par_prime_check = eigenVecs * par_prime_vec_col ;
      printf("\n\n par_vec_from_par_prime_check\n") ;
      par_vec_from_par_prime_check.Print() ;

      TRandom tran(12345) ;
      int npx(200) ;
      int npy = npx ;

      TH2F* h_toycheck = new TH2F( "h_toycheck", "", npx, xmin, xmax, npy, 0., 2. ) ;
      int ntoy = 10000 ;

      for ( int ti=0; ti<ntoy; ti++ ) {

         TMatrixD random_par_prim_vec_col(3,1) ;
         random_par_prim_vec_col[0][0] = tran.Gaus( par_prime_vec_col(0,0), sqrt( eigenVals(0) ) ) ;
         random_par_prim_vec_col[1][0] = tran.Gaus( par_prime_vec_col(1,0), sqrt( eigenVals(1) ) ) ;
         random_par_prim_vec_col[2][0] = tran.Gaus( par_prime_vec_col(2,0), sqrt( eigenVals(2) ) ) ;

         TMatrixD random_par_vec_col(3,1) ;
         random_par_vec_col = eigenVecs * random_par_prim_vec_col ;
         double this_p0 = random_par_vec_col(0,0) ;
         double this_p1 = random_par_vec_col(1,0) ;
         double this_p2 = random_par_vec_col(2,0) ;
         if (ti < 10 ) printf("  toy %4d :  p0 = %10.7f ,  p1 = %10.7f  ,  p2 = %10.7f\n", ti, this_p0, this_p1, this_p2 ) ;

         for ( int xi=0; xi<npx; xi++ ) {
            double x = xmin + (xi+0.5)*(xmax-xmin)/(1.*npx) ;
            double f = this_p0 + this_p1 * x + this_p2 * x * x ;
            h_toycheck -> Fill( x, f ) ;
         } // xi

      } // ti

      can -> cd(2) ;
      h_toycheck -> Draw("colz") ;
      gPad -> SetGridx(1) ;
      gPad -> SetGridy(1) ;


     //------
      TCanvas* can2 = new TCanvas( "can2", "", 900, 700 ) ;
      can2 -> Divide(3,2) ;
      gStyle -> SetOptStat("mr") ;
      gStyle -> SetStatW( 0.35 ) ;
      gStyle -> SetStatH( 0.20 ) ;

      TText* tt_label = new TText() ;
      for ( int bi=0; bi<n_hist_bins; bi++ ) {
         int h2bin = h_toycheck -> GetXaxis() -> FindBin( x_val[bi] ) ;
         sprintf( hname, "h_toy_xbin%d", bi ) ;
         TH1* hp = h_toycheck -> ProjectionY( hname, h2bin, h2bin ) ;
         hp -> SetFillColor(11) ;
         can2 -> cd( bi+1 ) ;
         hp -> Draw() ;
         gPad -> SetGridx(1) ;
         hp -> Draw( "axis same" ) ;
         hp -> Draw( "axig same" ) ;
         printf("  %2d :  x = %5.1f ,  function %9.5f +/- %9.5f\n", bi, x_val[bi], f_val[bi], f_err[bi] ) ;
         char label[1000] ;
         sprintf( label, "Func. = %7.4f +/- %7.4f", f_val[bi], f_err[bi] ) ;
         tt_label -> DrawTextNDC( 0.15, 0.80, label ) ;
      } // bi



   } // p2fit1



