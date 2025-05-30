#ifndef HiggsAnalysis_CombinedLimit_CombineCodegenImpl_h
#define HiggsAnalysis_CombinedLimit_CombineCodegenImpl_h

#include <ROOT/RConfig.hxx> // for ROOT_VERSION

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,35,0)
# define COMBINE_DECLARE_CODEGEN_IMPL(CLASS_NAME) \
    namespace RooFit { namespace Experimental { void codegenImpl(CLASS_NAME &arg, CodegenContext &ctx); }}
# define COMBINE_DECLARE_TRANSLATE
#elif ROOT_VERSION_CODE >= ROOT_VERSION(6,32,0)
# define COMBINE_DECLARE_CODEGEN_IMPL(CLASS_NAME)
# define COMBINE_DECLARE_TRANSLATE void translate(RooFit::Detail::CodeSquashContext &ctx) const override;
#else
# define COMBINE_DECLARE_CODEGEN_IMPL(_)
# define COMBINE_DECLARE_TRANSLATE
#endif

#endif
