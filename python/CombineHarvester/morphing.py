from __future__ import absolute_import
from os import environ

# Prevent cppyy's check for the PCH
environ['CLING_STANDARD_PCH'] = 'none'
import cppyy

cppyy.load_reflection_info("libCombineHarvesterCombinePdfs")

def BuildRooMorphing(ws, cb, bin, process, mass_var, norm_postfix='norm', allow_morph=True, verbose=False, force_template_limit=False, file=None):
    return cppyy.gbl.ch.BuildRooMorphing(ws, cb, bin, process, mass_var, norm_postfix, allow_morph, verbose, force_template_limit, file);

def BuildCMSHistFuncFactory(ws, cb, mass_vars, process):
  f = cppyy.gbl.ch.CMSHistFuncFactory()
  morphVars = cppyy.gbl.map['std::string','RooAbsReal*']()
  morphVars[process] = mass_vars
  f.SetHorizontalMorphingVariable(morphVars)
  f.Run(cb, ws)

