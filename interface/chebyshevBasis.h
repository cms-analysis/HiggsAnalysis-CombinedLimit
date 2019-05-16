#ifndef CHEBYSHEVBASIS
#define CHEBYSHEVBASIS

#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "RooFit.h"
#include "RooFormulaVar.h"
#include "RooConstVar.h"
#include "RooArgList.h"
#include "RooAddition.h"
#include "RooAbsArg.h"
#include "Riostream.h" 
#include "RooRealVar.h"
#include <map>

class chebyshevBasis : public TObject {
public:
    chebyshevBasis() {};
    chebyshevBasis( const char *name,
                    const char *title,
                    const Int_t _orderX,
                    const Int_t _orderY,
                    const TH2 &_binning );
    inline virtual ~chebyshevBasis (){};
    // float getChebyshevVal(float xval, float yval, float thisOrderX, float thisOrderY) const; 
    RooAddition getBinVal(const float xCenter, const float yCenter) const;
    void drawBasis(std::string file_name);

protected:
    // variables
    TH2F* binning;
    Int_t orderX;
    Int_t orderY;
    RooArgList coefList;
    Double_t xmin;
    Double_t xmax;
    Double_t ymin;
    Double_t ymax;
    float slope_x;
    float slope_y;
    // std::map<std::string, TH2F> polyHists;  

    // functions
    std::pair<float, float> mapToChebyshev(float ix,float iy) const;
    float Eval2DChebyshev(float x, float y, int thisOrderX, int thisOrderY) const;
    // TH2F Make2DChebyshev(const int& thisOrderX, const int& thisOrderY) const;
    

private:
    ClassDef(chebyshevBasis, 1)
};

#endif