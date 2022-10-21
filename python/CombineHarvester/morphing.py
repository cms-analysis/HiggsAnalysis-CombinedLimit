try:
    from libCombineHarvesterCombinePdfs import *
except ImportError:
    raise ImportError(
        'Module is missing: you need to compile the C++ -> python shared library by running make')


def BuildRooMorphing(ws, cb, bin, process, mass_var, norm_postfix='norm', allow_morph=True, verbose=False, force_template_limit=False, file=None):
    return BuildRooMorphingX(ws, cb, bin, process, mass_var, norm_postfix, allow_morph, verbose, force_template_limit, file)

def BuildCMSHistFuncFactory(ws, cb, mass_vars, process):
        return BuildCMSHistFuncFactoryX(ws, cb, mass_vars, process)
