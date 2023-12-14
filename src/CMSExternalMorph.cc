#include "../interface/CMSExternalMorph.h"

CMSExternalMorph::CMSExternalMorph() {}

CMSExternalMorph::CMSExternalMorph(
        const char* name,
        const char* title,
        RooRealVar& x,
        const std::vector<double>& edges
        ) :
    RooAbsReal(name, title),
    x_("x", "", this, x),
    edges_(edges)
{
}

CMSExternalMorph::CMSExternalMorph(CMSExternalMorph const& other, const char* name) :
    RooAbsReal(other, name),
    x_("x", this, other.x_),
    edges_(other.edges_)
{
}

CMSExternalMorph::~CMSExternalMorph() = default;

double CMSExternalMorph::evaluate() const {
    auto it = std::upper_bound(std::begin(edges_), std::end(edges_), x_->getVal());
    if ( (it == std::begin(edges_)) or (it == std::end(edges_)) ) {
        return 0.0;
    }
    size_t idx = std::distance(std::begin(edges_), it) - 1;
    return batchGetBinValues()[idx];
}
