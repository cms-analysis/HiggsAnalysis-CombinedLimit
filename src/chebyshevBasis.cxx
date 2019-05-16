#include "HiggsAnalysis/CombinedLimit/interface/chebyshevBasis.h"
#include <map>

ClassImp(chebyshevBasis)

chebyshevBasis::chebyshevBasis( const char *name,
                                const char *title,
                                const Int_t _orderX,
                                const Int_t _orderY,
                                const TH2 &_binning ) :
    orderX(_orderX), orderY(_orderY)

{
    std::string thisName;
    binning = static_cast<TH2F*>(_binning.Clone());
    // Build list of RooRealVar Coefficients
    for (Int_t ix=0; ix <= orderX; ix++){
        for (Int_t iy=0; iy <= orderY; iy++){
            thisName = "ChebCoeff_x"+std::to_string(ix)+"y"+std::to_string(iy);
            RooRealVar* thisRRV = new RooRealVar(thisName.c_str(),thisName.c_str(),0.01,-10.0,10.0);
            coefList.add(*thisRRV);
        }
    }
    

    // Make the var map variables
    xmin = _binning.GetXaxis()->GetXmin();
    xmax = _binning.GetXaxis()->GetXmax();
    ymin = _binning.GetYaxis()->GetXmin();
    ymax = _binning.GetYaxis()->GetXmax();
    slope_x = 2/(xmax-xmin);
    slope_y = 2/(ymax-ymin);

}

RooAddition chebyshevBasis::getBinVal(const float xCenter, const float yCenter) const {
    // Map "real" axis values to [-1,1] range
    std::pair<float, float> mappedvals = mapToChebyshev(xCenter,yCenter);

    Int_t xbin = binning->GetXaxis()->FindBin(xCenter);
    Int_t ybin = binning->GetYaxis()->FindBin(yCenter);

    std::string basisName;
    std::string coeffName;
    std::string formulaName;

    RooArgList toSum("toSum");
    for (Int_t ix=0; ix<=orderX; ix++){
        for (Int_t iy=0; iy<=orderY; iy++){
            // std::cout << "DEBUG: " << ix << ", " << iy << std::endl;
            // Bin value
            basisName = "ChebBasisVal_x"+std::to_string(ix)+"y"+std::to_string(iy)+"_bin"+std::to_string(xbin)+"-"+std::to_string(ybin);
            RooConstVar* basisConst = new RooConstVar(basisName.c_str(),basisName.c_str(),Eval2DChebyshev(mappedvals.first,mappedvals.second,ix,iy));
            // Coefficient
            coeffName = "ChebCoeff_x"+std::to_string(ix)+"y"+std::to_string(iy);
            RooAbsArg* coeff = coefList.find(coeffName.c_str());

            RooArgList* formulaInput = new RooArgList(*coeff,*basisConst);

            formulaName = "ChebFormula_x"+std::to_string(ix)+"y"+std::to_string(iy)+"_bin"+std::to_string(xbin)+"-"+std::to_string(ybin);
            RooFormulaVar* thisCheb = new RooFormulaVar(formulaName.c_str(),formulaName.c_str(),
                                    "@0*@1+abs(@0)",
                                    *formulaInput);    // a_ij * T_i(x)*T_j(y) + abs(a_ij)

            // std::cout << "DEBUG: " << formulaName << std::endl;

            toSum.add(*thisCheb);
        }
    }

    std::string thisName = "ChebSum_"+std::to_string(xbin)+"-"+std::to_string(ybin);
    RooAddition binVal(thisName.c_str(),thisName.c_str(),toSum);

    return binVal;
}

float chebyshevBasis::Eval2DChebyshev(float x, float y, int thisOrderX, int thisOrderY) const{
    // Initialize polynomials
    float Tx_n2 = 0; // T_n-2
    float Tx_n1 = 0; // T_n-1
    float Tx_n = 0;  // T_n
    float Ty_n2 = 0; // T_n-2
    float Ty_n1 = 0; // T_n-1
    float Ty_n = 0;  // T_n
    float Tx = 0;
    float Ty = 0;

    // Construct polynomial in x
    for(Int_t i=0; i<=thisOrderX; ++i){
        // T_0 == 1
        if (i == 0) {
            Tx_n = 1;
        }
        // T_1 = x
        else if (i == 1) {
            Tx_n1 = 1;
            Tx_n = x;
        }
        // T_n = 2*x*T_n-1 - T_n-2
        else {
            Tx_n2 = Tx_n1;
            Tx_n1 = Tx_n;
            Tx_n = 2*x*Tx_n1-Tx_n2;
        }
    }

    // Construct polynomial in y
    for (Int_t j=0; j<=thisOrderY; ++j){
        // T_0 == 1
        if (j == 0) {
            Ty_n = 1;
        }
        // T_1 = x
        else if (j == 1) {
            Ty_n1 = 1;
            Ty_n = y;
        }
        // T_n = 2*x*T_n-1 - T_n-2
        else {
            Ty_n2 = Ty_n1;
            Ty_n1 = Ty_n;
            Ty_n = 2*y*Ty_n1-Ty_n2;
        }
    }
    // Only want to shift if order is > 0
    Tx = Tx_n;
    Ty = Ty_n;

    return Tx*Ty;
}


std::pair<float, float> chebyshevBasis::mapToChebyshev(float ix,float iy) const {
    
    float xp = slope_x*(ix-xmax)+1;  // Map x
    float yp = slope_y*(iy-ymax)+1;  // Map y
    
    return std::make_pair(xp,yp);
}

void chebyshevBasis::drawBasis(std::string file_name) {

    TFile outfile(file_name.c_str(),"RECREATE");

    TCanvas basisCan("basisCan","basisCan",800,700);
    std::string canName;
    outfile.cd();
    
    std::string polyKey;
    std::string cheb_name;

    for (Int_t oX=0; oX<=orderX; ++oX){
        for (Int_t oY=0; oY<=orderY; ++oY) {
            polyKey = "x"+std::to_string(oX)+"y"+std::to_string(oY);
            
            cheb_name = "cheb_Tx"+std::to_string(oX)+"_Ty"+std::to_string(oY);
            const char *cheb_name_cstr = cheb_name.c_str();
            TH2F* cheb2D = new TH2F(cheb_name_cstr, cheb_name_cstr, 100, -1.0, 1.0, 100, -1.0, 1.0);

            // Loop over TH2F bins and evaluate
            for (Int_t xbin=1; xbin<=cheb2D->GetNbinsX(); xbin++){
                for (Int_t ybin=1; ybin<=cheb2D->GetNbinsY(); ybin++){
                    float xcenter = cheb2D->GetXaxis()->GetBinCenter(xbin);
                    float ycenter = cheb2D->GetYaxis()->GetBinCenter(ybin);
                    float val = Eval2DChebyshev(xcenter,ycenter,oX,oY);
                   
                    cheb2D->SetBinContent(xbin,ybin,val);
                }
            }

            cheb2D->Draw("surf");
            cheb2D->Write();

            canName = "basis_plots/x"+std::to_string(oX)+"y"+std::to_string(oY)+".png";
            basisCan.Print(canName.c_str(),"png");
        }
    }

    outfile.Close();
}