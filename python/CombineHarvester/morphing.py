from __future__ import absolute_import
import cppyy

cppyy.include("CombineHarvester/CombinePdfs/interface/MorphFunctions.h")


def BuildRooMorphing(ws, cb, bin, process, mass_var, norm_postfix='norm', allow_morph=True, verbose=False, force_template_limit=False, file=None):
    return cppyy.gbl.BuildRooMorphing(ws, cb, bin, process, mass_var, norm_postfix, allow_morph, verbose, force_template_limit, file);
