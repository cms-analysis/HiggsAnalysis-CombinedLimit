#ifndef HiggsAnalysis_CombinedLimit_CombineUtils_h
#define HiggsAnalysis_CombinedLimit_CombineUtils_h

#include <memory>

class RooAbsReal;
class RooAbsPdf;
class RooAbsData;
class RooArgSet;

std::unique_ptr<RooAbsReal> combineCreateNLL(RooAbsPdf& pdf,
                                             RooAbsData& data,
                                             RooArgSet const* constraint = nullptr,
                                             bool offset = false,
                                             bool warnAboutDifferentBackend = true);

// Class that makes the utility functions also available in Python (adding classes to the dictionaries is less brittle than adding functions)
class CombineUtils {
public:
  static std::unique_ptr<RooAbsReal> combineCreateNLL(RooAbsPdf& pdf,
                                                      RooAbsData& data,
                                                      RooArgSet const* constraint = nullptr,
                                                      bool offset = false,
                                                      bool warnAboutDifferentBackend = true) {
    return ::combineCreateNLL(pdf, data, constraint, offset, warnAboutDifferentBackend);
  }
};

#endif
