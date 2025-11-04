#ifndef HiggsAnalysis_CombinedLimit_CombineCodegenImpl_h
#define HiggsAnalysis_CombinedLimit_CombineCodegenImpl_h

#include <string>

class AsymPow;
class FastVerticalInterpHistPdf2;
class FastVerticalInterpHistPdf2D2;
class ProcessNormalization;
class VerticalInterpPdf;

namespace RooFit::Experimental {

  class CodegenContext;

  void codegenImpl(AsymPow& arg, CodegenContext& ctx);
  void codegenImpl(FastVerticalInterpHistPdf2& arg, CodegenContext& ctx);
  void codegenImpl(FastVerticalInterpHistPdf2D2& arg, CodegenContext& ctx);
  void codegenImpl(ProcessNormalization& arg, CodegenContext& ctx);
  void codegenImpl(VerticalInterpPdf& arg, CodegenContext& ctx);

  std::string codegenIntegralImpl(VerticalInterpPdf& arg, int code, const char* rangeName, CodegenContext& ctx);

}  // namespace RooFit::Experimental

#endif
