#ifndef TH1Keys_h
#define TH1Keys_h

#include <TH1.h>
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooNDKeysPdf.h>

class TH1Keys : public TH1 {
    public:
       TH1Keys();
       TH1Keys(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup, TString options = "a", Double_t rho = 1.5);
       TH1Keys(const char *name,const char *title,Int_t nbinsx,const Float_t  *xbins, TString options = "a", Double_t rho = 1.5) ;
       TH1Keys(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins, TString options = "a", Double_t rho = 1.5) ;
       TH1Keys(const TH1Keys &other);
       ~TH1Keys() override;

       TH1 * GetHisto() { if (!isCacheGood_) FillH1(); return cache_; }
       const TH1 * GetHisto() const { if (!isCacheGood_) FillH1(); return cache_; }

       Int_t    Fill(Double_t x) override { return Fill(x,1.0); }
       Int_t    Fill(Double_t x, Double_t w) override;
       void     FillN(Int_t ntimes, const Double_t *x, const Double_t *w, Int_t stride=1) override;
#if ROOT_VERSION_CODE <  ROOT_VERSION(5,34,00)
       virtual void     Add(const TH1 *h1, Double_t c1=1);
       virtual void     Add(const TH1 *h, const TH1 *h2, Double_t c1=1, Double_t c2=1) { dont("Add with two arguments"); } 
#else
       Bool_t   Add(const TH1 *h1, Double_t c1=1) override;
       Bool_t   Add(const TH1 *h, const TH1 *h2, Double_t c1=1, Double_t c2=1) override { dont("Add with two arguments"); return false;} 
#endif
       void     AddBinContent(Int_t bin) override { AddBinContent(bin, 1.0); }
       void     AddBinContent(Int_t bin, Double_t w) override { dont("AddBinContent"); }
       void     Copy(TObject &hnew) const override { dont("Copy"); }

       virtual TH1     *DrawCopy(Option_t *option="") const { dont("DrawCopy"); return 0; }

       Double_t GetBinContent(Int_t bin) const override { return GetHisto()->GetBinContent(bin); }
       Double_t GetBinContent(Int_t bin, Int_t) const override { return GetHisto()->GetBinContent(bin); }
       Double_t GetBinContent(Int_t bin, Int_t, Int_t) const override {return GetHisto()->GetBinContent(bin); }

       Double_t GetEntries() const override { return dataset_->numEntries(); }

       void     Reset(Option_t *option="") override ; 
       void     SetBinContent(Int_t bin, Double_t content) override { dont("SetBinContent"); }
       void     SetBinContent(Int_t bin, Int_t, Double_t content) override        { SetBinContent(bin,content); }
       void     SetBinContent(Int_t bin, Int_t, Int_t, Double_t content) override { SetBinContent(bin,content); }
       void     SetBinsLength(Int_t n=-1) override { dont("SetBinLength"); }
       void     Scale(Double_t c1=1, Option_t *option="") override;

       ClassDefOverride(TH1Keys,1)  //

    private:
        Double_t    min_, max_;
        RooRealVar *x_, *w_;
        RooArgSet   point_;
        RooDataSet *dataset_;
        Double_t    underflow_, overflow_;
        Double_t    globalScale_;

	TString          options_;
        Double_t           rho_;

        mutable TH1 *cache_;
        mutable bool isCacheGood_;

        void FillH1() const;

        void dont(const char *) const ;
}; // class

#endif
