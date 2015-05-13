#include "HiggsAnalysis/CombinedLimit/interface/RooSpline1D.h"

#include <stdexcept>

#include <fstream>
#include <sstream>

RooSpline1D::RooSpline1D(const char *name, const char *title, RooAbsReal &xvar, const char *path, const unsigned short xcol, const unsigned short ycol, const unsigned short skipLines, const char *algo) :
        RooAbsReal(name,title),
        xvar_("xvar","Variable", this, xvar),
        x_(), y_(), type_(algo),
        interp_(0)
{
        std::ifstream file( path, std::ios::in);
        std::string line;

        for(int lineno=0; std::getline(file, line); lineno++){
        	if(lineno<skipLines) continue;
            std::istringstream ss(line);
            std::istream_iterator<std::string > begin(ss), end;
            std::vector<std::string> tokens(begin, end);

            x_.push_back(atof(tokens[xcol].c_str()));
            y_.push_back(atof(tokens[ycol].c_str()));

//            std::cout << lineno << ": " << line << std::endl;
        }

        file.close();
}


RooSpline1D::RooSpline1D(const char *name, const char *title, RooAbsReal &xvar, unsigned int npoints, const double *xvals, const double *yvals, const char *algo) :
        RooAbsReal(name,title),
        xvar_("xvar","Variable", this, xvar), 
        x_(npoints), y_(npoints), type_(algo),
        interp_(0)
{ 
    for (unsigned int i = 0; i < npoints; ++i) {
        x_[i] = xvals[i];
        y_[i] = yvals[i];
    }
}

RooSpline1D::RooSpline1D(const char *name, const char *title, RooAbsReal &xvar, unsigned int npoints, const float *xvals, const float *yvals, const char *algo) :
        RooAbsReal(name,title),
        xvar_("xvar","Variable", this, xvar), 
        x_(npoints), y_(npoints), type_(algo),
        interp_(0)
{ 
    for (unsigned int i = 0; i < npoints; ++i) {
        x_[i] = xvals[i];
        y_[i] = yvals[i];
    }
}

RooSpline1D::RooSpline1D(const RooSpline1D &other, const char *newname) :
    RooAbsReal(other,newname),
    xvar_("xvar",this,other.xvar_),
    x_(other.x_), y_(other.y_), type_(other.type_),
    interp_(0)
{
}

RooSpline1D::~RooSpline1D() 
{
    delete interp_;
}


TObject *RooSpline1D::clone(const char *newname) const 
{
    return new RooSpline1D(*this, newname);
}

void RooSpline1D::init() const {
    delete interp_;
    if      (type_ == "CSPLINE") interp_ = new ROOT::Math::Interpolator(x_, y_, ROOT::Math::Interpolation::kCSPLINE);
    else if (type_ == "LINEAR") interp_ = new ROOT::Math::Interpolator(x_, y_, ROOT::Math::Interpolation::kLINEAR);
    else if (type_ == "POLYNOMIAL") interp_ = new ROOT::Math::Interpolator(x_, y_, ROOT::Math::Interpolation::kPOLYNOMIAL);
    else if (type_ == "CSPLINE_PERIODIC") interp_ = new ROOT::Math::Interpolator(x_, y_, ROOT::Math::Interpolation::kCSPLINE_PERIODIC);
    else if (type_ == "AKIMA") interp_ = new ROOT::Math::Interpolator(x_, y_, ROOT::Math::Interpolation::kAKIMA);
    else if (type_ == "AKIMA_PERIODIC") interp_ = new ROOT::Math::Interpolator(x_, y_, ROOT::Math::Interpolation::kAKIMA_PERIODIC);
    else throw std::invalid_argument("Unknown interpolation type '"+type_+"'");
}

Double_t RooSpline1D::evaluate() const {
    if (interp_ == 0) init();
    return interp_->Eval(xvar_);
}


ClassImp(RooSpline1D)
