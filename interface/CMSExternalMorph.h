#ifndef CMSExternalMorph_h
#define CMSExternalMorph_h
#include <vector>

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRealProxy.h"

class CMSExternalMorph : public RooAbsReal {
  public:
    CMSExternalMorph();
    /*
     * All subclasses need to provide an edges array of length nbins+1
     * of the observable (x)
     * TODO: CMSHistFunc and CMSHistSum do not check the binning is compatible
     * with their binning other than having the correct length
     */
    CMSExternalMorph(
        const char* name,
        const char* title,
        RooRealVar& x,
        const std::vector<double>& edges
    );
    CMSExternalMorph(CMSExternalMorph const& other, const char* name = 0);
    ~CMSExternalMorph() override;

    /* Batch accessor for CMSHistFunc / CMSHistSum, to be overriden by concrete
     * implementations. hasChanged() should indicate whether or not
     * batchGetBinValues() would return a new vector, given the state of
     * any dependent variables. 
     */
    virtual bool hasChanged() const = 0;
    virtual const std::vector<double>& batchGetBinValues() const = 0;

  protected:
    RooRealProxy x_;
    std::vector<double> edges_;

    double evaluate() const override;

  private:
    ClassDefOverride(CMSExternalMorph, 1)
};

#endif // CMSExternalMorph_h
