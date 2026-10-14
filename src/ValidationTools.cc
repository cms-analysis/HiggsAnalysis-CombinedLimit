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
        std::cout<<" Up/Down templates are identical: up variation: "<<sys->value_u()<<", down variation: "<<sys->value_d()<<std::endl;
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
      jsobj["emptyProcessShape"][proc->bin()].push_back(proc->process());
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
      std::cout<<" Uncertainty has a large normalisation effect: up variation: "<<sys->value_u()<<" Down variation: "<<sys->value_d()<<std::endl;
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
  auto bins = cb.bin_set();
  for(auto b : bins){
    auto cb_bin_backgrounds = cb.cp().bin({b}).backgrounds();
    const TH1F& tothist = cb_bin_backgrounds.GetShape();
    for(int i=1;i<=tothist.GetNbinsX();i++){
      if(tothist.GetBinContent(i)<=0){
        std::cout<<"Channel "<<b<<" bin "<<i<<" of the templates is empty in background"<<std::endl;
      }
    }
  }
}

void CheckEmptyBins(CombineHarvester& cb, json& jsobj){
  auto bins = cb.bin_set();
  for(auto b : bins){
    auto cb_bin_backgrounds = cb.cp().bin({b}).backgrounds();
    const TH1F& tothist = cb_bin_backgrounds.GetShape();
    for(int i=1;i<=tothist.GetNbinsX();i++){
      if(tothist.GetBinContent(i)<=0){
        jsobj["emptyBkgBin"][b].push_back(i);
      }
    }
  }
}

bool HistErrorsAreSqrtN(const TH1* hist_norm, double rate=1.0) {
  TH1* hist = (TH1*) hist_norm->Clone();
  hist->Scale(rate);
  for(int i=1; i<hist->GetNcells(); ++i) { // exclude under-/overflow
    double binContent = std::abs(hist->GetBinContent(i));
    if(binContent!=0 && !hist->IsBinUnderflow(i) && !hist->IsBinOverflow(i)){
      double binError   = hist->GetBinError(i);
      double errSquared = binError*binError;
      double threshold  = 1e-6 * std::max(binContent,errSquared);
      if(std::abs(binContent-errSquared)>threshold){
        return false; // at least one bin error is not sqrt(N)
      }
    }
  }
  return true; // all bin errors are sqrt(N) within 1e-6 tolerance
}

/**
 * Check if bin errors are set correctly.
 * 
 * To compute the bin-by-bin uncertainties correctly with autoMCStats,
 * the bin errors are expected to be the square root of the sum of weights, sqrt(sumw2).
 * This method checks the bin errors are not a Poisson error or sqrt(N),
 * and if the relative uncertainty is not too large (indicating too large weights, or a bug).
 */
void CheckBinErrors(CombineHarvester& cb, double maxRelBinErr){
  const std::vector<std::string> bins_bbb = Set2Vec(cb.GetAutoMCStatsBins());
  auto cb_bin = cb.cp().bin(bins_bbb);
  cb_bin.ForEachProc([&](ch::Process *proc){
    const TH1* shape = proc->shape();
    if(!shape){
      return;
    }
    const int errOpt = shape->GetBinErrorOption();
    if(errOpt!=TH1F::kNormal){
      std::string errEnum = (errOpt==TH1F::kPoisson ? "TH1::kPoisson" : errOpt==TH1F::kPoisson2 ? "TH1::kPoisson2" : std::to_string(errOpt));
      PrintProc(proc);
      std::cout<<" Bin error option is "<<errEnum<<"="<<errOpt<<"!=TH1::kNormal, but sqrt(sumw2) is needed for autoMCStats"<<std::endl;
      return;
    }
    if(shape->GetNcells()>=3 && shape->GetSumw2N()==0){
      PrintProc(proc);
      std::cout<<" Sumw2 is not defined, which is needed for autoMCStats"<<std::endl;
      return;
    }
    if(HistErrorsAreSqrtN(shape,proc->rate())){
      PrintProc(proc);
      std::cout<<" Bin errors are sqrt(N) instead of sqrt(sumw2), needed for autoMCStats"<<std::endl;
      return;
    }
  });
  if(maxRelBinErr<=0){
    return;
  }
  for(auto b : cb_bin.bin_set()){
    const TH1F& tothist = cb.cp().bin({b}).backgrounds().GetShape();
    for(int i=1;i<=tothist.GetNbinsX();i++){
      if(tothist.GetBinContent(i)!=0){
        double relErr = tothist.GetBinError(i) / tothist.GetBinContent(i);
        if(relErr>maxRelBinErr){
          std::cout<<"Channel "<<b<<" bin "<<i<<" of the templates is has large relative bin error: "<<relErr<<" > "<<maxRelBinErr<<std::endl;
        }
      }
    }
  }
}

void CheckBinErrors(CombineHarvester& cb, double maxRelBinErr, json& jsobj){
  std::vector<std::string> bins = Set2Vec(cb.GetAutoMCStatsBins());
  auto cb_bin = cb.cp().bin(bins);
  cb_bin.ForEachProc([&](ch::Process *proc){
    const TH1* shape = proc->shape();
    if(!shape){
      return;
    }
    const int errOpt = shape->GetBinErrorOption();
    if(errOpt!=TH1F::kNormal){
      std::string errEnum = (errOpt==TH1F::kPoisson ? "TH1::kPoisson" : errOpt==TH1F::kPoisson2 ? "TH1::kPoisson2" : std::to_string(errOpt));
      jsobj["binErrorIsNotSumw2"][proc->bin()][proc->process()] = "BinErrorOption="+errEnum;
      return;
    }
    if(shape->GetNcells()>=3 && shape->GetSumw2N()==0){
      jsobj["binErrorIsNotSumw2"][proc->bin()][proc->process()] = "NoSumw2Defined";
      return;
    }
    if(HistErrorsAreSqrtN(shape,proc->rate())){
      jsobj["binErrorIsNotSumw2"][proc->bin()][proc->process()] = "BinErrorIsSqrtN";
      return;
    }
  });
  for(auto b : cb_bin.bin_set()){
    const TH1F& tothist = cb.cp().bin({b}).backgrounds().GetShape();
    for(int i=1;i<=tothist.GetNbinsX();i++){
      if(tothist.GetBinContent(i)!=0){
        double relErr = tothist.GetBinError(i) / tothist.GetBinContent(i);
        if(relErr>maxRelBinErr){
          jsobj["largeBinError"][b].push_back(i);
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
        std::cout<<" Uncertainty probably has no genuine shape effect. Summed relative difference per bin between normalised nominal and up shape: "<<up_diff<<" between normalised nominal and down shape: "<<down_diff<<" . If you are using 1-bin shapes you can ignore this warning, but you can consider using lnN instead of a shape uncertainty"<<std::endl;
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
      if (hist_u->GetNbinsX() == 1) jsobj["smallShapeEff1bin"][sys->name()][sys->bin()][sys->process()]={{"diff_u",up_diff},{"diff_d",down_diff}};
      else {if(up_diff<diff_lim && down_diff<diff_lim) jsobj["smallShapeEff"][sys->name()][sys->bin()][sys->process()]={{"diff_u",up_diff},{"diff_d",down_diff}}; }
    }
  });
}

void CheckSmallSignals(CombineHarvester& cb, double minSigFrac){
  auto bins = cb.bin_set();
  for(auto b : bins){
    auto cb_bin_signals = cb.cp().bin({b}).signals();
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
    auto cb_bin = cb.cp().bin({b});
    double sigrate = cb_bin_signals.GetRate();
    for(auto p : cb_bin_signals.process_set()){
      if(cb_bin_signals.cp().process({p}).GetRate() < minSigFrac*sigrate){
        jsobj["smallSignalProc"][b][p]={{"sigrate_tot",sigrate},{"procrate",cb_bin_signals.cp().process({p}).GetRate()}};
      }
    }
  }
}

void ValidateCards(CombineHarvester& cb, std::string const& filename, double maxNormEff, double minSigFrac, double maxRelBinErr=2.){
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
    ValidateShapeTemplates(cb, output_js);
    CheckEmptyBins(cb, output_js);
    CheckBinErrors(cb, maxRelBinErr, output_js); // for accurate autoMCStats
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
