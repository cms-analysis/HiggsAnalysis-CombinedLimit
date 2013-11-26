#include "../interface/TestProposal.h"
#include "../interface/DebugProposal.h"
#include "../interface/VerticalInterpPdf.h"
#include "../interface/VerticalInterpHistPdf.h"
#include "../interface/AsymPow.h"
#include "../interface/CombDataSetFactory.h"
#include "../interface/TH1Keys.h"
#include "../interface/RooSimultaneousOpt.h"
#include "../interface/SimpleCacheSentry.h"
#include "../interface/th1fmorph.h"
#include "../interface/HZZ4LRooPdfs.h"
#include "../interface/HWWLVJRooPdfs.h"
#include "../interface/HZZ2L2QRooPdfs.h"
#include "../interface/HGGRooPdfs.h"
#include "../interface/HZGRooPdfs.h"
#include "../interface/SequentialMinimizer.h"
#include "../interface/ProcessNormalization.h"
#include "../interface/RooSpline1D.h"
#include "../interface/RooScaleLOSM.h"
#include "../interface/rVrFLikelihood.h"
#include "../interface/RooMultiPdf.h"
#include "../interface/RooBernsteinFast.h"
#include "../interface/SimpleGaussianConstraint.h"

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ class TestProposal+;
#pragma link C++ class DebugProposal+;
#pragma link C++ class VerticalInterpPdf+;
#pragma link C++ class VerticalInterpHistPdf+;
#pragma link C++ class FastTemplate+;
#pragma link C++ class FastHisto+;
#pragma link C++ class FastHisto2D+;
#pragma link C++ class FastVerticalInterpHistPdfBase+;
#pragma link C++ class FastVerticalInterpHistPdfBase::Morph+;
#pragma link C++ class FastVerticalInterpHistPdf+;
#pragma link C++ class FastVerticalInterpHistPdf2D+;
#pragma link C++ class FastVerticalInterpHistPdf2Base+;
#pragma link C++ class FastVerticalInterpHistPdf2+;
#pragma link C++ class FastVerticalInterpHistPdf2D2+;
#pragma link C++ class AsymPow+;
#pragma link C++ class CombDataSetFactory+;
#pragma link C++ class TH1Keys+;
#pragma link C++ class RooSimultaneousOpt+;
#pragma link C++ class SimpleGaussianConstraint+;
#pragma link C++ class SimpleCacheSentry+;
#pragma link C++ function th1fmorph;
#pragma link C++ class RooqqZZPdf+;
#pragma link C++ class RooggZZPdf+;
#pragma link C++ class RooqqZZPdf_v2+;
#pragma link C++ class RooggZZPdf_v2+;
#pragma link C++ class RooBetaFunc_v2+;
#pragma link C++ class Roo4lMasses2D+;
#pragma link C++ class Roo4lMasses2D_Bkg+;
#pragma link C++ class Roo4lMasses2D_BkgGGZZ+;
#pragma link C++ class RooFourMuMassShapePdf2+;
#pragma link C++ class RooFourEMassShapePdf2+;
#pragma link C++ class RooTwoETwoMuMassShapePdf2+;
#pragma link C++ class RooFourMuMassRes+;
#pragma link C++ class RooFourEMassRes+;
#pragma link C++ class RooTwoETwoMuMassRes+;
#pragma link C++ class RooRelBW1+;
#pragma link C++ class RooRelBWUF+;
#pragma link C++ class RooRelBWUF_SM4+;
#pragma link C++ class RooRelBWUFParam+;
#pragma link C++ class RooRelBWUFParamWidth+;
#pragma link C++ class RooRelBW+;
#pragma link C++ class RooRelBWHighMass+;
#pragma link C++ class RooDoubleCB+;
#pragma link C++ class RooaDoubleCBxBW+;
#pragma link C++ class RooCB+;
#pragma link C++ class RooFermi+;
#pragma link C++ class Triangle+;
#pragma link C++ class RooLevelledExp+;
#pragma link C++ class RooPower+;
#pragma link C++ class cmsmath::SequentialMinimizer+;
#pragma link C++ class ProcessNormalization+;
#pragma link C++ class RooSpline1D+;
#pragma link C++ class RooScaleLOSM+;
#pragma link C++ class RooScaleHGamGamLOSM+;
#pragma link C++ class RooScaleHGluGluLOSM+;
#pragma link C++ class RooScaleHGamGamLOSMPlusX+;
#pragma link C++ class RooScaleHGluGluLOSMPlusX+;
#pragma link C++ class RooTsallis;
#pragma link C++ class RooErfExpPdf+;
#pragma link C++ class RooAlpha+;
#pragma link C++ class RooAlphaExp+;
#pragma link C++ class RooBWRunPdf+;
#pragma link C++ class RooErfPow2Pdf+;
#pragma link C++ class RooAlpha4ErfPow2Pdf+;
#pragma link C++ class RooErfPowPdf+;
#pragma link C++ class RooAlpha4ErfPowPdf+;
#pragma link C++ class RooPow2Pdf+;
#pragma link C++ class RooPowPdf+;
#pragma link C++ class RooQCDPdf+;
#pragma link C++ class RooUser1Pdf+;
#pragma link C++ class RooAlpha4ErfPowExpPdf+;
#pragma link C++ class RooErfPowExpPdf+;
#pragma link C++ class RooStepBernstein+;
#pragma link C++ class RooGaussStepBernstein+;
#pragma link C++ class rVrFLikelihood+;
#pragma link C++ class RooMultiPdf+;
#pragma link C++ class RooBernsteinFast<1>+;
#pragma link C++ class RooBernsteinFast<2>+;
#pragma link C++ class RooBernsteinFast<3>+;
#pragma link C++ class RooBernsteinFast<4>+;
#pragma link C++ class RooBernsteinFast<5>+;
#pragma link C++ class RooBernsteinFast<6>+;
#pragma link C++ class RooBernsteinFast<7>+;
#pragma link C++ class RooAnaExpNPdf+;
#pragma link C++ class RooAlpha42ExpPdf+;
#pragma link C++ class Roo2ExpPdf+;
#pragma link C++ class RooExpNPdf+;
#pragma link C++ class RooAlpha4ExpNPdf+;
#pragma link C++ class RooExpTailPdf+;
#pragma link C++ class RooAlpha4ExpTailPdf+;
#pragma link C++ class RooAlpha4GausExpPdf+;
#pragma link C++ class RooGausExpPdf+;
#pragma link C++ class RooDoubleCrystalBall+;

#endif
