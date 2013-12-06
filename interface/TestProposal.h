#ifndef HiggsAnalysis_CombinedLimit_TestProposal_h
#define HiggsAnalysis_CombinedLimit_TestProposal_h

#include <Rtypes.h>

struct RooRealVar;

#include <RooArgSet.h>
#include <vector>

#include <RooStats/ProposalFunction.h>

class TestProposal : public RooStats::ProposalFunction {

   public:
      TestProposal() : RooStats::ProposalFunction() {}
      TestProposal(double divisor, const RooRealVar *alwaysStepMe=0) ;
      TestProposal(double divisor, const RooArgList &alwaysStepMe) ;

      void setDiscreteModels(RooRealVar *indicatorVar, RooArgSet vars, const std::vector<RooArgSet> &points) {
        discreteModelIndicator_ = indicatorVar;
        discreteModelVars_.removeAll();
        discreteModelVars_.add(vars);
        discreteModelPoints_.resize(points.size());
        for (unsigned int i = 0, n = points.size(); i < n; ++i) {
            discreteModelPoints_[i].removeAll();
            discreteModelPoints_[i].addClone(points[i]);        
        }
      } 

      // Populate xPrime with a new proposed point
      virtual void Propose(RooArgSet& xPrime, RooArgSet& x);

      // Determine whether or not the proposal density is symmetric for
      // points x1 and x2 - that is, whether the probabilty of reaching x2
      // from x1 is equal to the probability of reaching x1 from x2
      virtual Bool_t IsSymmetric(RooArgSet& x1, RooArgSet& x2) ;

      // Return the probability of proposing the point x1 given the starting
      // point x2
      virtual Double_t GetProposalDensity(RooArgSet& x1, RooArgSet& x2);

      virtual ~TestProposal() {}

      ClassDef(TestProposal,1) // A concrete implementation of ProposalFunction, that uniformly samples the parameter space.
    
    private:
      double divisor_, poiDivisor_;
      RooArgList alwaysStepMe_;

      RooRealVar *discreteModelIndicator_;
      RooArgSet   discreteModelVars_;
      std::vector<RooArgSet> discreteModelPoints_;
};

#endif
