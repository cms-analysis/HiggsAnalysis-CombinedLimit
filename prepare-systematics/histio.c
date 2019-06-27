#ifndef histio_c
#define histio_c

#include "TClass.h"
#include "TList.h"
#include "TFile.h"
#include "TIterator.h"
#include "TRegexp.h"
#include "TDirectory.h"
#include "TObject.h"
#include "TSystem.h"
#include "TKey.h"
#include "TH1.h"
//#include "histio.h"

#include <iostream>
#include <iomanip>

void histio()
{
}

void saveHist(const char* filename, const char* pat, bool delete_hists = false )
{
  printf("  Saving histograms in %s\n", filename ) ; fflush( stdout ) ;
  TList* list = gDirectory->GetList() ;
  TIterator* iter = list->MakeIterator();

  TRegexp re(pat,kTRUE) ;

  TFile outf(filename,"RECREATE") ;
  TObject* obj ;
  while((obj=iter->Next())) {    
    if (TString(obj->GetName()).Index(re)>=0) {
      obj->Write() ;
      std::cout << "." ;
    //  cout << setw(9) << counter++ << " : " << obj->GetName() << std::endl ;
    }
  }
  std::cout << std::endl ;
  outf.Close() ;

  delete iter ;
  if ( delete_hists ) {
     printf("\n\n Deleting histograms.\n\n") ;
     gDirectory->DeleteAll();//This line will remove histrograms from memory. Without this line, when we create a histogram, I will be saved everytime we call saveHist function. Important when we use "run_all.C"
  } else {
     printf("\n\n Not deleting histograms.\n\n") ;
  }

}


//void loadHist(const char* filename, const char* pfx, const char* pat, Bool_t doAdd, Double_t scaleFactor)
void loadHist(const char* filename="in.root", const char* pfx=0, const char* pat="*", Bool_t doAdd=kFALSE, Double_t scaleFactor=-1.0)
{
  ///////////cout << " Reading histograms from file: " << filename << endl ;
  TFile inf(filename) ;
  //inf.ReadAll() ;
  if ( ! inf.IsOpen() ) { gSystem -> Exit(-1) ; }
  TList* list = inf.GetListOfKeys() ;
  TIterator* iter = list->MakeIterator();

  TRegexp re(pat,kTRUE) ;
  //////std::cout << "pat = " << pat << std::endl ;

  gDirectory->cd("Rint:") ;

  TObject* obj ;
  TKey* key ;
  ///////////std::cout << "doAdd = " << (doAdd?"T":"F") << std::endl ;
  ///////////std::cout << "loadHist: reading." ;
  while((key=(TKey*)iter->Next())) {
   
    Int_t ridx = TString(key->GetName()).Index(re) ;    
    if (ridx==-1) {
      continue ;
    }

    obj = inf.Get(key->GetName()) ;
    TObject* clone ;
    if (pfx) {

      // Find existing TH1-derived objects
      TObject* oldObj = 0 ;
      if (doAdd){
	//////oldObj = gDirectory->Get(Form("%s_%s",pfx,obj->GetName())) ;
	oldObj = gDirectory->Get(Form("%s_%s",obj->GetName(),pfx)) ;
	if (oldObj && !oldObj->IsA()->InheritsFrom(TH1::Class())) {
	  oldObj = 0 ;
	}
      }
      if (oldObj) {
	clone = oldObj ;
   //// ((TH1*)clone)->Add((TH1*)obj) ;
        if ( scaleFactor > 0 ) {
           ((TH1*)clone)->Sumw2() ;
           ((TH1*)clone)->Add((TH1*)obj, scaleFactor) ;
        } else {
           ((TH1*)clone)->Add((TH1*)obj) ;
        }
      } else {
	/////clone = obj->Clone(Form("%s_%s",pfx,obj->GetName())) ;
	clone = obj->Clone(Form("%s_%s",obj->GetName(),pfx)) ;
      }


    } else {

      // Find existing TH1-derived objects
      TObject* oldObj = 0 ;
      if (doAdd){
	oldObj = gDirectory->Get(key->GetName()) ;
	if (oldObj && !oldObj->IsA()->InheritsFrom(TH1::Class())) {
	  oldObj = 0 ;
	}
      }

      if (oldObj) {
	clone = oldObj ;
 /////  ((TH1*)clone)->Add((TH1*)obj) ;
        if ( scaleFactor > 0 ) {
           ((TH1*)clone)->Sumw2() ;
           ((TH1*)clone)->Add((TH1*)obj, scaleFactor) ;
        } else {
           ((TH1*)clone)->Add((TH1*)obj) ;
        }
      } else {
	clone = obj->Clone() ;
      }
    }
    if ( scaleFactor > 0 && !doAdd ) {
       ((TH1*) clone)->Sumw2() ;
       ((TH1*) clone)->Scale(scaleFactor) ;
    }
    if (!gDirectory->GetList()->FindObject(clone)) {
      gDirectory->Append(clone) ;
    }
    //////std::cout << "." ;
    //////std::cout.flush() ;
  }
  ////////std::cout << std::endl;
  inf.Close() ;
  delete iter ;
}
#endif
