#ifndef HiggsAnalysis_CombinedLimit_CombineCodegenImpl_h
#define HiggsAnalysis_CombinedLimit_CombineCodegenImpl_h

#include <ROOT/RConfig.hxx> // for ROOT_VERSION

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,35,0)
# define COMBINE_DECLARE_CODEGEN_IMPL(CLASS_NAME) \
    namespace RooFit { namespace Experimental { void codegenImpl(CLASS_NAME &arg, CodegenContext &ctx); }}
# define COMBINE_DECLARE_CODEGEN_INTEGRAL_IMPL(CLASS_NAME) \
    namespace RooFit { namespace Experimental { std::string codegenIntegralImpl(CLASS_NAME &arg, int code, const char *rangeName, CodegenContext &ctx); }}
# define COMBINE_DECLARE_TRANSLATE
# define COMBINE_DECLARE_ANALYTICAL_INTEGRAL
#elif ROOT_VERSION_CODE >= ROOT_VERSION(6,32,0)
# define COMBINE_DECLARE_CODEGEN_IMPL(CLASS_NAME)
# define COMBINE_DECLARE_CODEGEN_INTEGRAL_IMPL(CLASS_NAME)
# define COMBINE_DECLARE_TRANSLATE \
    void translate(RooFit::Detail::CodeSquashContext &ctx) const override;
# define COMBINE_DECLARE_ANALYTICAL_INTEGRAL \
    std::string buildCallToAnalyticIntegral(Int_t code, const char *rangeName, RooFit::Detail::CodeSquashContext &ctx) const override;
#else
# define COMBINE_DECLARE_CODEGEN_IMPL(_)
# define COMBINE_DECLARE_CODEGEN_INTEGRAL_IMPL(_)
# define COMBINE_DECLARE_TRANSLATE
# define COMBINE_DECLARE_ANALYTICAL_INTEGRAL
#endif

#endif
