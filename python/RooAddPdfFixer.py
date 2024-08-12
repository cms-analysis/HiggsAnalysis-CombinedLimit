"""
Python wrapper around RooAddPdfFixer in utils.
Deprecated since 112x because ROOT bug has been fixed
https://sft.its.cern.ch/jira/browse/ROOT-6008
"""


def FixAll(workspace):
    raise RuntimeError("Since CMSSW_11_2_X, this utility is unnecessary. It will be removed in a future release.")
