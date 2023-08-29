#include "../interface/CMSInterferenceFunc.h"
#include "TBuffer.h"
#include <Eigen/Dense>

class _InterferenceEval {
  public:
    _InterferenceEval(std::vector<std::vector<double>> scaling_in, size_t ncoef) {
        binscaling_.reserve(scaling_in.size());
        for(const auto& bin : scaling_in) {
            Eigen::MatrixXd mat(ncoef, ncoef);
            size_t k=0;
            for(size_t i=0; i<ncoef; i++) {
                for(size_t j=0; j<=i; j++) {
                    mat(i, j) = mat(j, i) = bin[k++];
                }
            }
            binscaling_.push_back(std::move(mat));
        }
        coefficients_.resize(ncoef);
        values_.resize(binscaling_.size());
    };
    inline void setCoefficient(size_t i, double val) { coefficients_[i] = val; };
    void computeValues() {
        for (size_t i=0; i < binscaling_.size(); ++i) {
            values_[i] = coefficients_.dot(binscaling_[i]*coefficients_);
        }
    };
    const std::vector<double>& getValues() const { return values_; };

  private:
    std::vector<Eigen::MatrixXd> binscaling_;
    Eigen::VectorXd coefficients_;
    std::vector<double> values_;
};


CMSInterferenceFunc::CMSInterferenceFunc() {};

CMSInterferenceFunc::CMSInterferenceFunc(
    CMSInterferenceFunc const& other, const char* name
  ) :
    RooAbsReal(other, name), x_("x", this, other.x_),
    coefficients_("coefficients", this, other.coefficients_),
    edges_(other.edges_), binscaling_(other.binscaling_),
    sentry_(name ? TString(name) + "_sentry" : TString(other.GetName())+"_sentry", "")
{
}

CMSInterferenceFunc::CMSInterferenceFunc(
    const char* name, const char* title, RooRealVar& x,
    RooArgList const& coefficients, const std::vector<double>& edges,
    const std::vector<std::vector<double>> binscaling
  ) :
    RooAbsReal(name, title), x_("x", "", this, x),
    coefficients_("coefficients", "", this),
    edges_(edges), binscaling_(binscaling),
    sentry_(TString(name) + "_sentry", "")
{
    coefficients_.add(coefficients);
}

CMSInterferenceFunc::~CMSInterferenceFunc() = default;

void CMSInterferenceFunc::printMultiline(
        std::ostream& os, Int_t contents, Bool_t verbose, TString indent
        ) const {
    RooAbsReal::printMultiline(os, contents, verbose, indent);
    os << ">> Sentry: " << (sentry_.good() ? "clean" : "dirty") << "\n";
    sentry_.Print("v");
}

void CMSInterferenceFunc::initialize() const {
    // take the opportunity to validate persistent data
    size_t nbins = edges_.size() - 1;
    size_t ncoef = coefficients_.getSize();

    for (size_t i=0; i < ncoef; ++i) {
      if ( coefficients_.at(i) == nullptr ) {
        throw std::invalid_argument("Lost coefficient " + std::to_string(i));
      }
      if ( not coefficients_.at(i)->InheritsFrom("RooAbsReal") ) {
        throw std::invalid_argument(
            "Coefficient " + std::to_string(i) + " is not a RooAbsReal"
            );
      }
    }
    if ( binscaling_.size() != nbins ) {
        throw std::invalid_argument(
                "Number of bins as determined from bin edges ("
                + std::to_string(nbins) + ") does not match bin"
                " scaling array (" + std::to_string(binscaling_.size()) + ")"
                );
    }
    for (size_t i=0; i < nbins; ++i) {
        if ( binscaling_[i].size() != ncoef*(ncoef+1)/2 ) {
            throw std::invalid_argument(
                    "Length of bin scaling matrix lower triangle for bin " + std::to_string(i)
                    + " (" + std::to_string(binscaling_[i].size() ) + ") is not consistent"
                    + " with the length of the coefficient array (" + std::to_string(ncoef) + ")"
                    );
        }
    }

    sentry_.SetName(TString(GetName()) + "_sentry");
    sentry_.addVars(coefficients_);
    sentry_.setValueDirty();
    evaluator_ = std::make_unique<_InterferenceEval>(binscaling_, ncoef);
}

void CMSInterferenceFunc::updateCache() const {
    for (int i=0; i < coefficients_.getSize(); ++i) {
        auto* coef = static_cast<RooAbsReal const*>(coefficients_.at(i));
        if ( coef == nullptr ) throw std::runtime_error("Lost coef!");
        evaluator_->setCoefficient(i, coef->getVal());
    }
    evaluator_->computeValues();
    sentry_.reset();
}

double CMSInterferenceFunc::evaluate() const {
    if ( not evaluator_ ) initialize();
    if ( not sentry_.good() ) updateCache();

    auto it = std::upper_bound(std::begin(edges_), std::end(edges_), x_->getVal());
    if ( (it == std::begin(edges_)) or (it == std::end(edges_)) ) {
        return 0.0;
    }
    size_t idx = std::distance(std::begin(edges_), it) - 1;
    return evaluator_->getValues()[idx];
}

const std::vector<double>& CMSInterferenceFunc::batchGetBinValues() const {
    // we don't really expect the cache to be valid, as upstream callers are
    // managing their own and calling this only when dirty, but let's check anyway
    if ( not evaluator_ ) initialize();
    if ( not sentry_.good() ) updateCache();
    return evaluator_->getValues();
}

