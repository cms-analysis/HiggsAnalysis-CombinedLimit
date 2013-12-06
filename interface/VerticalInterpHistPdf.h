#ifndef ROO_VERTICA_INTERP_HIST
#define ROO_VERTICA_INTERP_HIST

/** Vertical interpolation between multiple histograms (or binned non-parametric pdfs, which eventually means histograms...)
    Does integration internally */

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "TH1.h"
#include "../interface/SimpleCacheSentry.h"
#include "../interface/FastTemplate.h"
#include <cmath>

class FastVerticalInterpHistPdf;
class FastVerticalInterpHistPdf2Base;
class FastVerticalInterpHistPdf2;
class FastVerticalInterpHistPdf2D2;

class VerticalInterpHistPdf : public RooAbsPdf {
public:

  VerticalInterpHistPdf() ;
  VerticalInterpHistPdf(const char *name, const char *title, const RooRealVar &x, const RooArgList& funcList, const RooArgList& coefList, Double_t smoothRegion=1., Int_t smoothAlgo=1) ;
  VerticalInterpHistPdf(const VerticalInterpHistPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new VerticalInterpHistPdf(*this,newname) ; }
  virtual ~VerticalInterpHistPdf() ;

  Bool_t selfNormalized() const { return kTRUE; }

  Double_t evaluate() const ;

  const RooArgList& funcList() const { return _funcList ; }
  const RooArgList& coefList() const { return _coefList ; }

  const RooRealProxy &x() const { return _x; }
  friend class FastVerticalInterpHistPdf;
protected:
  RooRealProxy   _x;
  RooListProxy _funcList ;   //  List of component FUNCs
  RooListProxy _coefList ;  //  List of coefficients
  Double_t     _smoothRegion;
  Int_t        _smoothAlgo;
  TIterator* _funcIter ;     //! Iterator over FUNC list
  TIterator* _coefIter ;    //! Iterator over coefficient list

  // TH1 containing the histogram of this pdf
  mutable SimpleCacheSentry _sentry; // !not to be serialized
  mutable TH1  *_cacheTotal;     //! not to be serialized
  // For additive morphing, histograms of fNominal, fUp and fDown
  // For multiplicative morphing, histograms of fNominal, log(fUp/fNominal), -log(fDown/fNominal);
  mutable TH1  **_cacheSingle; //! not to be serialized
  mutable int  *_cacheSingleGood; //! not to be serialized
private:

  ClassDef(VerticalInterpHistPdf,1) // 

protected:
  void setupCaches() const ;
  void syncTotal() const ;
  void syncComponent(int which) const ;
  // return a smooth function that is equal to +/-1 for |x| >= smoothRegion_ and it's null in zero
  inline double smoothStepFunc(double x) const { 
    if (fabs(x) >= _smoothRegion) return x > 0 ? +1 : -1;
    double xnorm = x/_smoothRegion, xnorm2 = xnorm*xnorm;
    return 0.125 * xnorm * (xnorm2 * (3.*xnorm2 - 10.) + 15);
  }
};

class FastVerticalInterpHistPdfBase : public RooAbsPdf {
public:

  FastVerticalInterpHistPdfBase() ;
  FastVerticalInterpHistPdfBase(const char *name, const char *title, const RooArgSet &obs, const RooArgList& funcList, const RooArgList& coefList, Double_t smoothRegion=1., Int_t smoothAlgo=1) ;
  FastVerticalInterpHistPdfBase(const FastVerticalInterpHistPdfBase& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const = 0; 
  virtual ~FastVerticalInterpHistPdfBase() ;

  Bool_t selfNormalized() const { return kTRUE; }

  Double_t evaluate() const = 0;

  const RooArgList& funcList() const { return _funcList ; }
  const RooArgList& coefList() const { return _coefList ; }

  /// Must be public, for serialization
  struct Morph { FastTemplate sum; FastTemplate diff; };

  friend class FastVerticalInterpHistPdf2Base;
protected:
  RooRealProxy   _x;
  RooListProxy _funcList ;   //  List of component FUNCs
  RooListProxy _coefList ;  //  List of coefficients
  Double_t     _smoothRegion;
  Int_t        _smoothAlgo;
  TIterator* _funcIter ;     //! Iterator over FUNC list
  TIterator* _coefIter ;    //! Iterator over coefficient list

  // TH1 containing the histogram of this pdf
  mutable bool              _init; //! not to be serialized
  mutable SimpleCacheSentry _sentry; //! not to be serialized

  // For additive morphing, histograms of (fUp-f0)+(fDown-f0) and (fUp-f0)-(fDown-f0)
  // For multiplicative morphing, log(fUp/f0)+log(fDown/f0),  log(fUp/f0)-log(fDown/f0)
  mutable std::vector<Morph> _morphs;  //! not to be serialized
  mutable std::vector<RooAbsReal *> _morphParams; //! not to be serialized

  // Prepare morphing data for a triplet of templates
  void syncMorph(Morph &out, const FastTemplate &nominal, FastTemplate &lo, FastTemplate &hi) const;

  // Do the vertical morphing from nominal value and morphs into cache. 
  // Do not normalize yet, as that depends on the dimension of the template
  void syncTotal(FastTemplate &cache, const FastTemplate &cacheNominal, const FastTemplate &cacheNominalLog) const ;

  // return a smooth function that is equal to +/-1 for |x| >= smoothRegion_ and it's null in zero
  inline double smoothStepFunc(double x) const { 
    if (fabs(x) >= _smoothRegion) return x > 0 ? +1 : -1;
    double xnorm = x/_smoothRegion, xnorm2 = xnorm*xnorm;
    return 0.125 * xnorm * (xnorm2 * (3.*xnorm2 - 10.) + 15);
  }

private:
  ClassDef(FastVerticalInterpHistPdfBase,2) // 
};


struct FastVerticalInterpHistPdfV;
class FastVerticalInterpHistPdf : public FastVerticalInterpHistPdfBase {
public:

  FastVerticalInterpHistPdf() : FastVerticalInterpHistPdfBase() {}
  FastVerticalInterpHistPdf(const char *name, const char *title, const RooRealVar &x, const RooArgList& funcList, const RooArgList& coefList, Double_t smoothRegion=1., Int_t smoothAlgo=1) :
    FastVerticalInterpHistPdfBase(name,title,RooArgSet(x),funcList,coefList,smoothRegion,smoothAlgo),
    _x("x","Independent variable",this,const_cast<RooRealVar&>(x)),
    _cache(), _cacheNominal(), _cacheNominalLog()  {}
  FastVerticalInterpHistPdf(const FastVerticalInterpHistPdf& other, const char* name=0) :
    FastVerticalInterpHistPdfBase(other, name),
    _x("x",this,other._x),
    _cache(other._cache), _cacheNominal(other._cacheNominal), _cacheNominalLog(other._cacheNominalLog)  {}
  virtual TObject* clone(const char* newname) const { return new FastVerticalInterpHistPdf(*this,newname) ; }
  virtual ~FastVerticalInterpHistPdf() {}

  Double_t evaluate() const ;

  Bool_t hasCache()     const { return _cache.size() > 0; }
  Bool_t isCacheReady() const { return _cache.size() > 0 && _init; }
  friend class FastVerticalInterpHistPdfV;
  friend class FastVerticalInterpHistPdf2;
protected:
  RooRealProxy   _x;

  /// Cache of the result
  mutable FastHisto _cache; //! not to be serialized
  /// Cache of nominal pdf (additive morphing) and its bin-by-bin logarithm (multiplicative)
  mutable FastHisto _cacheNominal; //! not to be serialized
  mutable FastHisto _cacheNominalLog; //! not to be serialized

  void setupCaches() const ;
  void syncTotal() const ;
  void syncComponent(int which) const ;
  void syncNominal() const ;
  void syncComponents(int dimension) const ;

private:
  ClassDef(FastVerticalInterpHistPdf,1) // 
};

class FastVerticalInterpHistPdfV {
    public: 
        FastVerticalInterpHistPdfV(const FastVerticalInterpHistPdf &, const RooAbsData &data) ;
        void fill(std::vector<Double_t> &out) const ;
    private:
        const FastVerticalInterpHistPdf & hpdf_;
        int begin_, end_, nbins_;
        struct Block { 
            int index, begin, end; 
            Block(int i, int begin_, int end_) : index(i), begin(begin_), end(end_) {}
        };
        std::vector<Block> blocks_;
        std::vector<int> bins_;
};

class FastVerticalInterpHistPdf2D : public FastVerticalInterpHistPdfBase {
public:

  FastVerticalInterpHistPdf2D() : FastVerticalInterpHistPdfBase() {}
  /// If conditional = false, the pdf is normalized integrating on (x,y). 
  /// If conditional = true, the pdf is separately normalized integrating on (y) for each specific (x) bin
  FastVerticalInterpHistPdf2D(const char *name, const char *title, const RooRealVar &x, const RooRealVar &y, bool conditional, const RooArgList& funcList, const RooArgList& coefList, Double_t smoothRegion=1., Int_t smoothAlgo=1) :
    FastVerticalInterpHistPdfBase(name,title,RooArgSet(x,y),funcList,coefList,smoothRegion,smoothAlgo),
    _x("x","Independent variable",this,const_cast<RooRealVar&>(x)),
    _y("y","Independent variable",this,const_cast<RooRealVar&>(y)),
    _conditional(conditional),
    _cache(), _cacheNominal(), _cacheNominalLog()  {}
  FastVerticalInterpHistPdf2D(const FastVerticalInterpHistPdf2D& other, const char* name=0) :
    FastVerticalInterpHistPdfBase(other, name),
    _x("x",this,other._x), _y("y",this,other._y), _conditional(other._conditional),
    _cache(other._cache), _cacheNominal(other._cacheNominal), _cacheNominalLog(other._cacheNominalLog)  {}
  virtual TObject* clone(const char* newname) const { return new FastVerticalInterpHistPdf2D(*this,newname) ; }
  virtual ~FastVerticalInterpHistPdf2D() {}

  const RooRealVar & x() const { return dynamic_cast<const RooRealVar &>(_x.arg()); }
  const RooRealVar & y() const { return dynamic_cast<const RooRealVar &>(_y.arg()); }
  Bool_t conditional() const { return _conditional; }

  Double_t evaluate() const ;

  Bool_t hasCache()     const { return _cache.size() > 0; }
  Bool_t isCacheReady() const { return _cache.size() > 0 && _init; }

  friend class FastVerticalInterpHistPdf2D2;
protected:
  RooRealProxy _x, _y;
  bool _conditional;

  /// Cache of the result
  mutable FastHisto2D _cache; //! not to be serialized
  /// Cache of nominal pdf (additive morphing) and its bin-by-bin logarithm (multiplicative)
  mutable FastHisto2D _cacheNominal; //! not to be serialized
  mutable FastHisto2D _cacheNominalLog; //! not to be serialized

  void setupCaches() const ;
  void syncTotal() const ;
  void syncComponent(int which) const ;
  void syncNominal() const ;
  void syncComponents(int dimension) const ;

private:
  ClassDef(FastVerticalInterpHistPdf2D,1) // 
};


class FastVerticalInterpHistPdf2Base : public RooAbsPdf {
public:

  FastVerticalInterpHistPdf2Base() ;
  FastVerticalInterpHistPdf2Base(const char *name, const char *title, const RooArgSet &obs, const TList & funcList, const RooArgList& coefList, Double_t smoothRegion=1., Int_t smoothAlgo=1) ;
  FastVerticalInterpHistPdf2Base(const FastVerticalInterpHistPdf2Base& other, const char* name=0) ;
  explicit FastVerticalInterpHistPdf2Base(const FastVerticalInterpHistPdfBase& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const = 0; 
  virtual ~FastVerticalInterpHistPdf2Base() ;

  Bool_t selfNormalized() const { return kTRUE; }

  Double_t evaluate() const = 0;

  const RooArgList& coefList() const { return _coefList ; }

  virtual void setActiveBins(unsigned int bins) {}

  /// Must be public, for serialization
  typedef FastVerticalInterpHistPdfBase::Morph Morph;
protected:
  RooListProxy _coefList ;  //  List of coefficients
  Double_t     _smoothRegion;
  Int_t        _smoothAlgo;

  // set to true after _morphParams and _sentry are initialized 
  mutable bool _initBase;   //! not to be serialized

  // to check if parameters change
  mutable SimpleCacheSentry _sentry; //! not to be serialized

  // For additive morphing, histograms of (fUp-f0)+(fDown-f0) and (fUp-f0)-(fDown-f0)
  // For multiplicative morphing, log(fUp/f0)+log(fDown/f0),  log(fUp/f0)-log(fDown/f0)
  // NOTE: it's the responsibility of the daughter to make sure these are initialized!!!
  std::vector<Morph> _morphs;  

  // Coefficients of the list in _coefList, already dynamic_cast'ed and in a vector
  mutable std::vector<const RooAbsReal *> _morphParams; //! not to be serialized

  // Prepare morphing data for a triplet of templates
  void initMorph(Morph &out, const FastTemplate &nominal, FastTemplate &lo, FastTemplate &hi) const;

  // Do the vertical morphing from nominal value and morphs into cache. 
  // Do not normalize yet, as that depends on the dimension of the template
  void syncTotal(FastTemplate &cache, const FastTemplate &cacheNominal, const FastTemplate &cacheNominalLog) const ;

  // return a smooth function that is equal to +/-1 for |x| >= smoothRegion_ and it's null in zero
  inline double smoothStepFunc(double x) const { 
    if (fabs(x) >= _smoothRegion) return x > 0 ? +1 : -1;
    double xnorm = x/_smoothRegion, xnorm2 = xnorm*xnorm;
    return 0.125 * xnorm * (xnorm2 * (3.*xnorm2 - 10.) + 15);
  }

  // initialize the morphParams and the sentry. to be called by the daughter class, sets also _initBase to true
  void initBase() const ; 

private:
  ClassDef(FastVerticalInterpHistPdf2Base,1) // 
};


struct FastVerticalInterpHistPdf2V;
class FastVerticalInterpHistPdf2 : public FastVerticalInterpHistPdf2Base {
public:

  FastVerticalInterpHistPdf2() : FastVerticalInterpHistPdf2Base() {}
  FastVerticalInterpHistPdf2(const char *name, const char *title, const RooRealVar &x, const TList & funcList, const RooArgList& coefList, Double_t smoothRegion=1., Int_t smoothAlgo=1) ;

  FastVerticalInterpHistPdf2(const FastVerticalInterpHistPdf2& other, const char* name=0) :
    FastVerticalInterpHistPdf2Base(other, name),
    _x("x",this,other._x),
    _cache(other._cache), _cacheNominal(other._cacheNominal), _cacheNominalLog(other._cacheNominalLog)  {}
  explicit FastVerticalInterpHistPdf2(const FastVerticalInterpHistPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new FastVerticalInterpHistPdf2(*this,newname) ; }
  virtual ~FastVerticalInterpHistPdf2() {}

  virtual void setActiveBins(unsigned int bins) ;
  Double_t evaluate() const ;

  friend class FastVerticalInterpHistPdf2V;
protected:
  RooRealProxy   _x;

  /// Cache of the result
  mutable FastHisto _cache; //! not to be serialized

  /// Cache of nominal pdf (additive morphing) and its bin-by-bin logarithm (multiplicative)
  FastHisto _cacheNominal; 
  FastHisto _cacheNominalLog; 

  void syncTotal() const ;
  void initNominal(TObject *nominal) ;
  void initComponent(int which, TObject *hi, TObject *lo) ;

private:
  ClassDef(FastVerticalInterpHistPdf2,1) // 
};
class FastVerticalInterpHistPdf2V {
    public: 
        FastVerticalInterpHistPdf2V(const FastVerticalInterpHistPdf2 &, const RooAbsData &data) ;
        void fill(std::vector<Double_t> &out) const ;
    private:
        const FastVerticalInterpHistPdf2 & hpdf_;
        int begin_, end_, nbins_;
        struct Block { 
            int index, begin, end; 
            Block(int i, int begin_, int end_) : index(i), begin(begin_), end(end_) {}
        };
        std::vector<Block> blocks_;
        std::vector<int> bins_;
};



class FastVerticalInterpHistPdf2D2 : public FastVerticalInterpHistPdf2Base {
public:

  FastVerticalInterpHistPdf2D2() : FastVerticalInterpHistPdf2Base() {}
  FastVerticalInterpHistPdf2D2(const char *name, const char *title, const RooRealVar &x, const RooRealVar &y, bool conditional, const TList & funcList, const RooArgList& coefList, Double_t smoothRegion=1., Int_t smoothAlgo=1) ;

  FastVerticalInterpHistPdf2D2(const FastVerticalInterpHistPdf2D2& other, const char* name=0) :
    FastVerticalInterpHistPdf2Base(other, name),
    _x("x",this,other._x), _y("y",this,other._y), _conditional(other._conditional),
    _cache(other._cache), _cacheNominal(other._cacheNominal), _cacheNominalLog(other._cacheNominalLog)  {}
  explicit FastVerticalInterpHistPdf2D2(const FastVerticalInterpHistPdf2D& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new FastVerticalInterpHistPdf2D2(*this,newname) ; }
  virtual ~FastVerticalInterpHistPdf2D2() {}

  const RooRealVar & x() const { return dynamic_cast<const RooRealVar &>(_x.arg()); }
  const RooRealVar & y() const { return dynamic_cast<const RooRealVar &>(_y.arg()); }
  Bool_t conditional() const { return _conditional; }

  Double_t evaluate() const ;
protected:
  RooRealProxy _x, _y;
  bool _conditional;

  /// Cache of the result
  mutable FastHisto2D _cache; //! not to be serialized

  /// Cache of nominal pdf (additive morphing) and its bin-by-bin logarithm (multiplicative)
  FastHisto2D _cacheNominal; 
  FastHisto2D _cacheNominalLog; 

  void syncTotal() const ;
  void initNominal(TObject *nominal) ;
  void initComponent(int which, TObject *hi, TObject *lo) ;

private:
  ClassDef(FastVerticalInterpHistPdf2D2,1) // 
};




#endif
