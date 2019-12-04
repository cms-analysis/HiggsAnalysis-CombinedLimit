#include "CombineHarvester/CombineTools/interface/ValidationTools.h"
#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <fstream>
#include <map>
#include "boost/format.hpp"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooAbsReal.h"
#include "RooAbsData.h"
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"

namespace ch {
using json = nlohmann::json;

void PrintSystematic(ch::Systematic *syst){
  std::cout<<"Systematic "<<syst->name()<<" for process "<<syst->process()<<" in region "<<syst->bin()<<" :";
}

void PrintProc(ch::Process *proc){
  std::cout<<"Process "<<proc->process()<<" in region "<<proc->bin()<<" :";
}



void ValidateShapeUncertaintyDirection(CombineHarvester& cb, json& jsobj){
  cb.ForEachSyst([&](ch::Systematic *sys){
    if(sys->type()=="shape" && ( (sys->value_u() > 1. && sys->value_d() > 1.) || (sys->value_u() < 1. && sys->value_d() < 1.))){
      jsobj["uncertVarySameDirect"][sys->name()][sys->bin()][sys->process()]={{"value_u",sys->value_u()},{"value_d",sys->value_d()}};
    }
  });
}

void ValidateShapeUncertaintyDirection(CombineHarvester& cb){
  cb.ForEachSyst([&](ch::Systematic *sys){
    if(sys->type()=="shape" && ( (sys->value_u() > 1. && sys->value_d() > 1.) || (sys->value_u() < 1. && sys->value_d() < 1.))){
      PrintSystematic(sys);
      std::cout<<" Up/Down normalisations go in the same direction: up variation: "<<sys->value_u()<<", down variation: "<<sys->value_d()<<std::endl;
    }
  });
}

void ValidateShapeTemplates(CombineHarvester& cb, json& jsobj){
  cb.ForEachSyst([&](ch::Systematic *sys){
  const TH1* hist_u;
  const TH1* hist_d;
  if(sys->type()=="shape" && ( fabs(sys->value_u() - sys->value_d()) < 0.0000001)){
      hist_u = sys->shape_u();
      hist_d = sys->shape_d();
      bool is_same=1;
      for(int i=1;i<=hist_u->GetNbinsX();i++){
        if(fabs(hist_u->GetBinContent(i))+fabs(hist_d->GetBinContent(i))>0){
          if(2*double(fabs(hist_u->GetBinContent(i)-hist_d->GetBinContent(i)))/(fabs(hist_u->GetBinContent(i))+fabs(hist_d->GetBinContent(i)))>0.001) is_same = 0;
        }
      }
      if(is_same){
        jsobj["uncertTemplSame"][sys->name()][sys->bin()][sys->process()]={{"value_u",sys->value_u()},{"value_d",sys->value_d()}};
      }
    }
  });
}

void ValidateShapeTemplates(CombineHarvester& cb){
  cb.ForEachSyst([&](ch::Systematic *sys){
    const TH1* hist_u;
    const TH1* hist_d;
    if(sys->type()=="shape" && ( fabs(sys->value_u() - sys->value_d()) < 0.0000001)){
      hist_u = sys->shape_u();
      hist_d = sys->shape_d();
      bool is_same=1;
      for(int i=1;i<=hist_u->GetNbinsX();i++){
        if(fabs(hist_u->GetBinContent(i))+fabs(hist_d->GetBinContent(i))>0){
          if(2*double(fabs(hist_u->GetBinContent(i)-hist_d->GetBinContent(i)))/(fabs(hist_u->GetBinContent(i))+fabs(hist_d->GetBinContent(i)))>0.001) is_same = 0;
        }
      }
      if(is_same){
        PrintSystematic(sys);
        std::cout<<"Up/Down templates are identical: up variation: "<<sys->value_u()<<", down variation: "<<sys->value_d()<<std::endl;
      }
    }
  });
}



void CheckEmptyShapes(CombineHarvester& cb, json& jsobj){
  std::vector<ch::Process*> empty_procs;
  auto bins = cb.bin_set();
  cb.ForEachProc([&](ch::Process *proc){
    if(proc->rate()==0.){
      empty_procs.push_back(proc); 
      if (jsobj["emptyProcessShape"][proc->bin()] !=NULL){
        jsobj["emptyProcessShape"][proc->bin()].push_back(proc->process());
      } else {
        jsobj["emptyProcessShape"][proc->bin()] = {proc->process()};
      }
   }
  });
  cb.ForEachSyst([&](ch::Systematic *sys){
    bool no_check=0;
    for( unsigned int i=0; i< empty_procs.size(); i++){
      if ( MatchingProcess(*sys,*empty_procs.at(i)) ) no_check=1;
    }
    if(!no_check){
      if(sys->type()=="shape" &&  (sys->value_u()==0. || sys->value_d()==0.)){
        jsobj["emptySystematicShape"][sys->name()][sys->bin()][sys->process()]={{"value_u",sys->value_u()},{"value_d",sys->value_d()}};
      }
    }
  });
}

void CheckEmptyShapes(CombineHarvester& cb){
  std::vector<ch::Process*> empty_procs;
  cb.ForEachProc([&](ch::Process *proc){
    if(proc->rate()==0.){
      empty_procs.push_back(proc); 
      PrintProc(proc);
      std::cout<<" has 0 yield"<<std::endl;
   }
  });
  cb.ForEachSyst([&](ch::Systematic *sys){
    bool no_check=0;
    for( unsigned int i=0; i< empty_procs.size(); i++){
      if ( MatchingProcess(*sys,*empty_procs.at(i)) ) no_check=1;
    }
    if(!no_check){
      if(sys->type()=="shape" &&  (sys->value_u()==0. || sys->value_d()==0.)){
        PrintSystematic(sys);
        std::cout<<" At least one empty histogram: up variation: "<<sys->value_u()<<" Down variation: "<<sys->value_d()<<std::endl;
      }
    }
  });
}

void CheckNormEff(CombineHarvester& cb, double maxNormEff){
  std::vector<ch::Process*> empty_procs;
  cb.ForEachProc([&](ch::Process *proc){
    if(proc->rate()==0.){
      empty_procs.push_back(proc); 
   }
  });
  cb.ForEachSyst([&](ch::Systematic *sys){
    bool no_check=0;
    for( unsigned int i=0; i< empty_procs.size(); i++){
      if ( MatchingProcess(*sys,*empty_procs.at(i)) ) no_check=1;
    }
    if(!no_check && ((sys->type()=="shape" &&  (sys->value_u()-1 > maxNormEff || sys->value_u()-1 < -maxNormEff  || sys->value_d()-1>maxNormEff || sys->value_d()-1< -maxNormEff)) || (sys->type()=="lnN" && (sys->value_u()-1 > maxNormEff || sys->value_u()-1 < - maxNormEff) ))){
      PrintSystematic(sys);
      std::cout<<"Uncertainty has a large normalisation effect: up variation: "<<sys->value_u()<<" Down variation: "<<sys->value_d()<<std::endl;
    }
  });
}

void CheckNormEff(CombineHarvester& cb, double maxNormEff, json& jsobj){
  std::vector<ch::Process*> empty_procs;
  cb.ForEachProc([&](ch::Process *proc){
    if(proc->rate()==0.){
      empty_procs.push_back(proc); 
   }
  });
  cb.ForEachSyst([&](ch::Systematic *sys){
    bool no_check=0;
    for( unsigned int i=0; i< empty_procs.size(); i++){
      if ( MatchingProcess(*sys,*empty_procs.at(i)) ) no_check=1;
    }
    if(!no_check && ((sys->type()=="shape" &&  (std::abs(sys->value_u()-1) > maxNormEff || std::abs(sys->value_d()-1)>maxNormEff)) || (sys->type()=="lnN" && (std::abs(sys->value_u()-1) > maxNormEff) ))){
      jsobj["largeNormEff"][sys->name()][sys->bin()][sys->process()]={{"value_u",sys->value_u()},{"value_d",sys->value_d()}};
    }
  });
}

void CheckEmptyBins(CombineHarvester& cb){
  TH1F tothist;
  auto bins = cb.bin_set();
  for(auto b : bins){
    auto cb_bin_backgrounds = cb.cp().bin({b}).backgrounds();
    tothist = cb_bin_backgrounds.GetShape();
    for(int i=1;i<=tothist.GetNbinsX();i++){
      if(tothist.GetBinContent(i)<=0){ 
        std::cout<<"Channel "<<b<<" bin "<<i<<" of the templates is empty in background"<<std::endl;
      }
    }
  }
}

void CheckEmptyBins(CombineHarvester& cb, json& jsobj){
  TH1F tothist;
  auto bins = cb.bin_set();
  for(auto b : bins){
    auto cb_bin_backgrounds = cb.cp().bin({b}).backgrounds();
    tothist = cb_bin_backgrounds.GetShape();
    for(int i=1;i<=tothist.GetNbinsX();i++){
      if(tothist.GetBinContent(i)<=0){ 
        if (jsobj["emptyBkgBin"][b] !=NULL){
          jsobj["emptyBkgBin"][b].push_back(i);
        } else {
          jsobj["emptyBkgBin"][b] = {i};
        }
      }
    }
  }
}


void CheckSizeOfShapeEffect(CombineHarvester& cb){
  double diff_lim=0.001;
  cb.ForEachSyst([&](ch::Systematic *sys){
    const TH1* hist_u;
    const TH1* hist_d;
    TH1F hist_nom;
    if(sys->type()=="shape"){
      hist_u = sys->shape_u();
      hist_d = sys->shape_d();
      hist_nom=cb.cp().bin({sys->bin()}).process({sys->process()}).GetShape();
      hist_nom.Scale(1./hist_nom.Integral());
      double up_diff=0;
      double down_diff=0;
      for(int i=1;i<=hist_u->GetNbinsX();i++){
        if(fabs(hist_u->GetBinContent(i))+fabs(hist_nom.GetBinContent(i))>0){
          up_diff+=2*double(fabs(hist_u->GetBinContent(i)-hist_nom.GetBinContent(i)))/(fabs(hist_u->GetBinContent(i))+fabs(hist_nom.GetBinContent(i)));
        }
        if(fabs(hist_d->GetBinContent(i))+fabs(hist_nom.GetBinContent(i))>0){
          down_diff+=2*double(fabs(hist_d->GetBinContent(i)-hist_nom.GetBinContent(i)))/(fabs(hist_d->GetBinContent(i))+fabs(hist_nom.GetBinContent(i)));
        }
      }
      if(up_diff<diff_lim && down_diff<diff_lim){
        PrintSystematic(sys);
        std::cout<<"Uncertainty probably has no genuine shape effect. Summed relative difference per bin between normalised nominal and up shape: "<<up_diff<<" between normalised nominal and down shape: "<<down_diff<<std::endl;
      }
    }
  });
}


void CheckSizeOfShapeEffect(CombineHarvester& cb, json& jsobj){
  double diff_lim=0.001;
  cb.ForEachSyst([&](ch::Systematic *sys){
    const TH1* hist_u;
    const TH1* hist_d;
    TH1F hist_nom;
    if(sys->type()=="shape"){
      hist_u = sys->shape_u();
      hist_d = sys->shape_d();
      hist_nom=cb.cp().bin({sys->bin()}).process({sys->process()}).GetShape();
      hist_nom.Scale(1./hist_nom.Integral());
      double up_diff=0;
      double down_diff=0;
      for(int i=1;i<=hist_u->GetNbinsX();i++){
        if(fabs(hist_u->GetBinContent(i))+fabs(hist_nom.GetBinContent(i))>0){
          up_diff+=2*double(fabs(hist_u->GetBinContent(i)-hist_nom.GetBinContent(i)))/(fabs(hist_u->GetBinContent(i))+fabs(hist_nom.GetBinContent(i)));
        }
        if(fabs(hist_d->GetBinContent(i))+fabs(hist_nom.GetBinContent(i))>0){
          down_diff+=2*double(fabs(hist_d->GetBinContent(i)-hist_nom.GetBinContent(i)))/(fabs(hist_d->GetBinContent(i))+fabs(hist_nom.GetBinContent(i)));
        }
      }
      if(up_diff<diff_lim && down_diff<diff_lim) jsobj["smallShapeEff"][sys->name()][sys->bin()][sys->process()]={{"diff_u",up_diff},{"diff_d",down_diff}};
    } 
  });
}

void CheckSmallSignals(CombineHarvester& cb,double minSigFrac){
  auto bins = cb.bin_set();
  for(auto b : bins){
    auto cb_bin_signals = cb.cp().bin({b}).signals();
    auto cb_bin_backgrounds = cb.cp().bin({b}).backgrounds();
    auto cb_bin = cb.cp().bin({b}); 
    double sigrate = cb_bin_signals.GetRate();
    for(auto p : cb_bin_signals.process_set()){
      if(cb_bin_signals.cp().process({p}).GetRate() < minSigFrac*sigrate){
        std::cout<<"Very small signal process. In bin "<<b<<" signal process "<<p<<" has yield "<<cb_bin_signals.cp().process({p}).GetRate()<<". Total signal rate in this bin is "<<sigrate<<std::endl;
      }
    }
  }
}


void CheckSmallSignals(CombineHarvester& cb, double minSigFrac, json& jsobj){
  auto bins = cb.bin_set();
  for(auto b : bins){
    auto cb_bin_signals = cb.cp().bin({b}).signals();
    auto cb_bin_backgrounds = cb.cp().bin({b}).backgrounds();
    auto cb_bin = cb.cp().bin({b}); 
    double sigrate = cb_bin_signals.GetRate();
    for(auto p : cb_bin_signals.process_set()){
      if(cb_bin_signals.cp().process({p}).GetRate() < minSigFrac*sigrate){
        jsobj["smallSignalProc"][b][p]={{"sigrate_tot",sigrate},{"procrate",cb_bin_signals.cp().process({p}).GetRate()}};
      }
    }
  }
}
  

void ValidateCards(CombineHarvester& cb, std::string const& filename, double maxNormEff, double minSigFrac){
 json output_js; 
 bool is_shape_card=1;
 cb.ForEachProc([&](ch::Process *proc){
   if(proc->pdf()||!(proc->shape())){
     is_shape_card=0;
    }
 });
 if(is_shape_card){     
   ValidateShapeUncertaintyDirection(cb, output_js);      
   CheckSizeOfShapeEffect(cb, output_js);
   ValidateShapeTemplates(cb,output_js);
   CheckEmptyBins(cb,output_js);
 } else {
   std::cout<<"Not a shape-based datacard / shape-based datacard using RooDataHist. Skipping checks on systematic shapes."<<std::endl;
 }
 CheckEmptyShapes(cb, output_js);      
 CheckNormEff(cb, maxNormEff, output_js);
 CheckSmallSignals(cb,minSigFrac, output_js);
 std::ofstream outfile(filename);
 outfile <<std::setw(4)<<output_js<<std::endl;
}

}
