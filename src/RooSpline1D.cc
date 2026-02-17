#include "../interface/RooSpline1D.h"

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


#include <RooFitHS3/RooJSONFactoryWSTool.h>
#include <RooFit/Detail/JSONInterface.h>
#include <RooFitHS3/JSONIO.h>

using RooFit::Detail::JSONNode;

class RooSpline1DStreamer : public RooFit::JSONIO::Exporter {
public:
   std::string const &key() const override {
     const static std::string keystring = "spline";
     return keystring;                          
   }    

   bool exportObject(RooJSONFactoryWSTool *, const RooAbsArg *func, RooFit::Detail::JSONNode &elem) const override
   {
      auto const *sp = static_cast<::RooSpline1D const *>(func);

      elem["type"] << key();

      // Independent variable
      elem["x"] << sp->x().GetName();

      // Unified schema field: interpolation
      TString interpolation = sp->type();
      if (interpolation == "CSPLINE") interpolation = "poly3";
      if (interpolation == "POLYNOMIAL") interpolation = "polyN";
      if (interpolation == "LINEAR") interpolation = "lin";
      interpolation.ToLower();
      elem["interpolation"] << interpolation;

      // Knots as primitive arrays
      auto const &xv = sp->xVals();
      auto const &yv = sp->yVals();

      auto &x0 = elem["x0"].set_seq();
      auto &y0 = elem["y0"].set_seq();

      // Be strict: these should always be same length, but don’t silently truncate.
      if (xv.size() != yv.size()) {
         RooJSONFactoryWSTool::error("RooSpline1D '" + std::string(sp->GetName()) +
                                    "' has x/y size mismatch (" + std::to_string(xv.size()) +
                                    " vs " + std::to_string(yv.size()) + ")");
      }

      for (std::size_t i = 0; i < xv.size(); ++i) {
         x0.append_child() << xv[i];
         y0.append_child() << yv[i];
      }

      return true;
   }
};


ClassImp(RooSpline1D)

class RooSpline1DFactory : public RooFit::JSONIO::Importer {
public:
   bool importArg(RooJSONFactoryWSTool *tool, const JSONNode &p) const override
   {
      const std::string name(RooJSONFactoryWSTool::name(p));

      // Required fields
      if (!p.has_child("x")) {
         RooJSONFactoryWSTool::error("no x given in '" + name + "'");
      }
      if (!p.has_child("x0") || !p.has_child("y0")) {
         RooJSONFactoryWSTool::error("no x0/y0 given in '" + name + "'");
      }

      // Unified schema: reject log flags unless false/missing
      const bool logx = p.has_child("logx") ? p["logx"].val_bool() : false;
      const bool logy = p.has_child("logy") ? p["logy"].val_bool() : false;
      if (logx || logy) {
         RooJSONFactoryWSTool::error("RooSpline1D cannot represent logx/logy in '" + name + "'");
      }

      RooAbsReal *x = tool->requestArg<RooAbsReal>(p, "x");

      // Unified schema field: interpolation
      std::string interpolation = p.has_child("interpolation") ? p["interpolation"].val() : "poly3";

      // Reject poly5 (no quintic mode in RooSpline1D)
      if (interpolation == "poly5") {
         RooJSONFactoryWSTool::error("RooSpline1D cannot represent interpolation 'poly5' in '" + name + "'");
      }

      // Map unified poly3 -> RooSpline1D type
      TString algo = interpolation;
      if (algo == "poly3") algo = "CSPLINE";
      if (algo == "polyN") algo = "POLYNOMIAL";      
      if (algo == "lin") algo = "LINEAR";      
      algo.ToUpper();

      // Validate against the RooSpline1D::init() supported set
      const bool supported =
         (algo == "CSPLINE") ||
         (algo == "LINEAR") ||
         (algo == "POLYNOMIAL") ||
         (algo == "CSPLINE_PERIODIC") ||
         (algo == "AKIMA") ||
         (algo == "AKIMA_PERIODIC");

      if (!supported) {
         RooJSONFactoryWSTool::error("unsupported interpolation '" + interpolation + "' (mapped to algo '" + algo +
                                    "') for RooSpline1D in '" + name + "'");
      }

      // Read knots
      std::vector<double> x0;
      std::vector<double> y0;
      x0.reserve(p["x0"].num_children());
      y0.reserve(p["y0"].num_children());

      for (const auto &v : p["x0"].children()) x0.push_back(v.val_double());
      for (const auto &v : p["y0"].children()) y0.push_back(v.val_double());

      if (x0.size() != y0.size()) {
         RooJSONFactoryWSTool::error("x0/y0 size mismatch in '" + name + "': x0 has " + std::to_string(x0.size()) +
                                    ", y0 has " + std::to_string(y0.size()));
      }
      if (x0.size() < 2) {
         RooJSONFactoryWSTool::error("need at least 2 knots in '" + name + "'");
      }

      // Construct:
      // RooSpline1D(name, title, xvar, npoints, xvals, yvals, algo)
      tool->wsEmplace<::RooSpline1D>(name.c_str(), *x,
                                    static_cast<unsigned int>(x0.size()),
                                    x0.data(), y0.data(),
                                    algo.Data());
     
      return true;
   }
};


#include "static_execute.h"

STATIC_EXECUTE([]() {
  using namespace RooFit::JSONIO;
  
  registerImporter<RooSpline1DFactory>("spline", false);
  registerExporter<RooSpline1DStreamer>(RooSpline1D::Class(), false);
 });
  

