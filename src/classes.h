#include "HiggsAnalysis/CombinedLimit/interface/TestProposal.h"
#include "HiggsAnalysis/CombinedLimit/interface/DebugProposal.h"
#include "HiggsAnalysis/CombinedLimit/interface/VerticalInterpPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/VerticalInterpHistPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/AsymPow.h"
#include "HiggsAnalysis/CombinedLimit/interface/CombDataSetFactory.h"
#include "HiggsAnalysis/CombinedLimit/interface/TH1Keys.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooSimultaneousOpt.h"
#include "HiggsAnalysis/CombinedLimit/interface/SimpleCacheSentry.h"
#include "HiggsAnalysis/CombinedLimit/interface/th1fmorph.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4LRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/HWWLVJRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/HGGRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZGRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/SequentialMinimizer.h"
#include "HiggsAnalysis/CombinedLimit/interface/ProcessNormalization.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooSpline1D.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooSplineND.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooScaleLOSM.h"
#include "HiggsAnalysis/CombinedLimit/interface/rVrFLikelihood.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooBernsteinFast.h"
#include "HiggsAnalysis/CombinedLimit/interface/SimpleGaussianConstraint.h"
#include "HiggsAnalysis/CombinedLimit/interface/AtlasPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4L_RooSpinZeroPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4L_RooSpinZeroPdf_phase.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4L_RooSpinZeroPdf_2D.h"
#include "HiggsAnalysis/CombinedLimit/interface/HWWLVJJRooPdfs.h"

namespace {
    struct dictionary {
	RooBernsteinFast<1> my_RooBernsteinFast_1;
	RooBernsteinFast<2> my_RooBernsteinFast_2;
	RooBernsteinFast<3> my_RooBernsteinFast_3;
	RooBernsteinFast<4> my_RooBernsteinFast_4;
	RooBernsteinFast<5> my_RooBernsteinFast_5;
	RooBernsteinFast<6> my_RooBernsteinFast_6;
	RooBernsteinFast<7> my_RooBernsteinFast_7;
    };
}
