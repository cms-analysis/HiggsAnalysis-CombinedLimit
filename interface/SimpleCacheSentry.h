#ifndef ROO_SIMPLE_CACHE_SENTRY
#define ROO_SIMPLE_CACHE_SENTRY

#include "RooRealVar.h"
#include "RooSetProxy.h"

class SimpleCacheSentry : public RooAbsArg {
    public:
        SimpleCacheSentry() ;
        SimpleCacheSentry(const char *name, const char *title) ;
        SimpleCacheSentry(const RooRealVar &var) ;
        SimpleCacheSentry(const RooAbsCollection &vars) ;
        SimpleCacheSentry(const RooAbsArg &func, const RooArgSet *obs=0) ;
        SimpleCacheSentry(const SimpleCacheSentry &other, const char *newname = 0) ;
        RooSetProxy & deps() { return _deps; }
        const RooArgSet & deps() const { return _deps; }
        void addArg(const RooAbsArg &arg) { _deps.add(arg); }
        void addVar(const RooRealVar &var) { _deps.add(var); } 
        void addVars(const RooAbsCollection &vars) ; 
        void addFunc(const RooAbsArg &func, const RooArgSet *obs=0) ;
        bool good() const { return !isValueDirty(); } 
        bool empty() const { return _deps.getSize() == 0; }
        void reset() { clearValueDirty(); } 
        // base class methods to be implemented
        TObject* clone(const char* newname) const override { return new SimpleCacheSentry(*this, newname); }
        RooAbsArg *createFundamental(const char* newname=0) const override { return 0; }
        Bool_t readFromStream(std::istream& is, Bool_t compact, Bool_t verbose=kFALSE) override { return false; }
        void writeToStream(std::ostream& os, Bool_t compact) const override { }
        Bool_t operator==(const RooAbsArg& other) const override { return this == &other; }
        void syncCache(const RooArgSet* nset=0) override {}
        void copyCache(const RooAbsArg* source, Bool_t valueOnly=kFALSE, Bool_t setValDirty=kTRUE) override {}
        void attachToTree(TTree& t, Int_t bufSize=32000) override {}
        void attachToVStore(RooVectorDataStore& vstore) override {}
        void setTreeBranchStatus(TTree& t, Bool_t active) override {}
        void fillTreeBranch(TTree& t) override {}
        Bool_t isIdentical(const RooAbsArg& other, Bool_t assumeSameType=kFALSE) const override;
    private:
        RooSetProxy _deps;
        ClassDefOverride(SimpleCacheSentry,1) 
};

#endif
