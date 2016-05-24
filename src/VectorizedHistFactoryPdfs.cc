#include "HiggsAnalysis/CombinedLimit/interface/VectorizedHistFactoryPdfs.h"
#include <memory>
#include <RooRealVar.h>

cacheutils::VectorizedHistFunc::VectorizedHistFunc(const RooHistFunc &pdf, bool includeZeroWeights) :
    pdf_(&pdf), data_(0), includeZeroWeights_(includeZeroWeights)
{
}

void
cacheutils::VectorizedHistFunc::fill() 
{
    RooArgSet obs(*data_->get());
    std::auto_ptr<RooArgSet> params(pdf_->getObservables(obs));

    yvals_.reserve(data_->numEntries());
    for (unsigned int i = 0, n = data_->numEntries(); i < n; ++i) {
        params->assignValueOnly(*data_->get(i), true);
        if (data_->weight() || includeZeroWeights_) {
            yvals_.push_back(pdf_->getVal());        
        }
    }
}

const std::vector<Double_t> & 
cacheutils::VectorizedHistFunc::eval(const RooAbsData &data) 
{
    if (data_ != &data) {
        data_ = &data;
        fill();
    }
    return yvals_;
}

#if 0
cacheutils::VectorizedHistFunc::VectorizedHistFunc(const RooHistFunc &pdf, const RooAbsData &data, bool includeZeroWeights) 
    
{
    RooArgSet obs(*data.get());
    std::auto_ptr<RooArgSet> params(pdf.getParameters(data));
    if (params->getSize() != 0) throw std::invalid_argument("HistFunc with parameters !?");

    //const RooRealVar * x_ = dynamic_cast<const RooRealVar*>(obs.first());

    yvals_.reserve(data.numEntries());
    for (unsigned int i = 0, n = data.numEntries(); i < n; ++i) {
        obs.assignValueOnly(*data.get(i), true);
        if (data.weight() || includeZeroWeights) {
            //xvals_.push_back(x_->getVal());        
            yvals_.push_back(pdf.getVal());        
        }
    }
}

void
cacheutils::VectorizedHistFunc::fill(std::vector<Double_t> &out) const 
{
    out = yvals_;
}
#endif

cacheutils::VectorizedParamHistFunc::VectorizedParamHistFunc(const ParamHistFunc &pdf, const RooAbsData &data, bool includeZeroWeights) 
    
{
    RooArgSet obs(*data.get());
    yvars_.reserve(data.numEntries());
    for (unsigned int i = 0, n = data.numEntries(); i < n; ++i) {
        obs.assignValueOnly(*data.get(i), true);
        if (data.weight() || includeZeroWeights) {
            RooRealVar * rrv = & pdf.getParameter();
            yvars_.push_back(rrv);        
        }
    }
}

void 
cacheutils::VectorizedParamHistFunc::fill(std::vector<Double_t> &out) const 
{
    out.resize(yvars_.size());
    for (unsigned int i = 0, n = yvars_.size(); i < n; ++i) {
        out[i] = yvars_[i]->getVal();
    }
}

namespace {
    class PiecewiseInterpolationWithAccessor : public PiecewiseInterpolation {
        public:
            PiecewiseInterpolationWithAccessor(const PiecewiseInterpolation &p) :
                PiecewiseInterpolation(p) {}
            const RooAbsReal & nominal() const { return _nominal.arg(); }
            int getInterpCode(int i) const { return _interpCode[i]; }
            bool positiveDefinite() const { return _positiveDefinite; }
    };
}

cacheutils::CachingPiecewiseInterpolation::CachingPiecewiseInterpolation(const PiecewiseInterpolation &pdf, const RooArgSet &obs) :
    pdf_(&pdf)
{
    PiecewiseInterpolationWithAccessor fixme(pdf);
    const RooArgList & highList  = pdf.highList();
    const RooArgList & lowList   = pdf.lowList();
    const RooAbsReal & nominal   = fixme.nominal(); 
    const RooArgList & coeffs = pdf.paramList();
    //std::cout << "Making a CachingPiecewiseInterpolation for " << pdf.GetName() << " with " << highList.getSize() << " pdfs, " << coeffs.getSize() << " coeffs." << std::endl;
    positiveDefinite_ = fixme.positiveDefinite();
    cachingPdfNominal_.reset(makeCachingPdf(const_cast<RooAbsReal*>(&nominal),&obs));
    for (int i = 0, n = highList.getSize(); i < n; ++i) {
        RooAbsReal *pdfiHi = (RooAbsReal*) highList.at(i);
        cachingPdfsHi_.push_back(makeCachingPdf(pdfiHi, &obs));
        RooAbsReal *pdfiLo = (RooAbsReal*) lowList.at(i);
        cachingPdfsLow_.push_back(makeCachingPdf(pdfiLo, &obs));
        coeffs_.push_back((RooAbsReal*) coeffs.at(i));
        codes_.push_back(fixme.getInterpCode(i));
        //std::cout << "      PiecewiseInterpolation Adding " << pdf.GetName() << "[" << i << "] Hi    : " << pdfiHi->ClassName() << " " << pdfiHi->GetName() << " using " << typeid(cachingPdfsHi_.back()).name() << std::endl;
        //std::cout << "      PiecewiseInterpolation Adding " << pdf.GetName() << "[" << i << "] Hi    : " << pdfiHi->ClassName() << " " << pdfiHi->GetName() << " using " << typeid(cachingPdfsHi_.back()).name() << std::endl;
        //std::cout << "      PiecewiseInterpolation Adding " << pdf.GetName() << "[" << i << "] Coeff : " << coeffs_.back()->ClassName() << " " << coeffs_.back()->GetName()  << std::endl;
        //std::cout << "      PiecewiseInterpolation Adding " << pdf.GetName() << "[" << i << "] Code  : " << codes_.back()  << std::endl;
    }
}

cacheutils::CachingPiecewiseInterpolation::~CachingPiecewiseInterpolation()
{
}


const std::vector<Double_t> & cacheutils::CachingPiecewiseInterpolation::eval(const RooAbsData &data)
{
    const std::vector<Double_t> & nominal = cachingPdfNominal_->eval(data);
    unsigned int size = nominal.size();
    work_.resize(size);
    std::copy(nominal.begin(), nominal.end(), work_.begin());
    for (int i = 0, n =  coeffs_.size(); i < n; ++i) {
        double param = coeffs_[i]->getVal();
        int    code  = codes_[i];
        switch(code){
            case 0: 
                {
                    if (param > 0) {
                        const std::vector<Double_t> & hi = cachingPdfsHi_[i].eval(data);
                        for (unsigned int j = 0; j < size; ++j) {
                            work_[j] += param * (hi[j] - nominal[j]);
                        }
                    } else {
                        const std::vector<Double_t> & lo = cachingPdfsLow_[i].eval(data);
                        for (unsigned int j = 0; j < size; ++j) {
                            work_[j] += param * (nominal[j] - lo[j]);
                        }
                    }
                } break;
            case 2:
                {
                    if (param > 0) {
                        const std::vector<Double_t> & hi = cachingPdfsHi_[i].eval(data);
                        for (unsigned int j = 0; j < size; ++j) {
                            work_[j] *= std::pow(hi[j]/nominal[j], param);
                        }
                    } else {
                        const std::vector<Double_t> & lo = cachingPdfsLow_[i].eval(data);
                        for (unsigned int j = 0; j < size; ++j) {
                            work_[j] *= std::pow(lo[j]/nominal[j], -param);
                        }
                    }
                } break;
            case 4:
                {
                    if (param > 1.) {
                        const std::vector<Double_t> & hi = cachingPdfsHi_[i].eval(data);
                        for (unsigned int j = 0; j < size; ++j) {
                            work_[j] += param * (hi[j] - nominal[j]);
                        }
                    } else if (param < -1.){
                        const std::vector<Double_t> & lo = cachingPdfsLow_[i].eval(data);
                        for (unsigned int j = 0; j < size; ++j) {
                            work_[j] += param * (nominal[j] - lo[j]);
                        }
                    } else {
                        const std::vector<Double_t> & hi = cachingPdfsHi_[i].eval(data);
                        const std::vector<Double_t> & lo = cachingPdfsLow_[i].eval(data);
                        for (unsigned int j = 0; j < size; ++j) {
                            double eps_plus  = (hi[j]-nominal[j]);
                            double eps_minus = (nominal[j] - lo[j]);
                            double S = 0.5*(eps_plus+eps_minus);
                            double A = 0.0625*(eps_plus-eps_minus);
                            double val = nominal[j] + param * (S + param * A * ( 15 + param * param * (-10 + param * param * 3  ) ) );
                            if (val < 0) val = 0;
                            work_[j] += (val - nominal[j]);
                        }
                    }
                } break;
            default:
                std::cout << "Interpolation code " << code << " not implemented. Sorry" << std::endl;
                throw std::invalid_argument("Bad interpolation code in CachingPiecewiseInterpolation");
        }
    }
    if (positiveDefinite_) {
        for (unsigned int j = 0; j < size; ++j) {
            if (work_[j] < 0) work_[j] = 0;
        }
    }
    return work_;
}

void cacheutils::CachingPiecewiseInterpolation::setDataDirty()
{
    cachingPdfNominal_->setDataDirty();
    for (CachingPdfBase &pdf : cachingPdfsHi_) pdf.setDataDirty();
    for (CachingPdfBase &pdf : cachingPdfsLow_) pdf.setDataDirty();
}

void cacheutils::CachingPiecewiseInterpolation::setIncludeZeroWeights(bool includeZeroWeights) 
{
    cachingPdfNominal_->setIncludeZeroWeights(includeZeroWeights);
    for (CachingPdfBase &pdf : cachingPdfsHi_) pdf.setIncludeZeroWeights(includeZeroWeights);
    for (CachingPdfBase &pdf : cachingPdfsLow_) pdf.setIncludeZeroWeights(includeZeroWeights);
}



