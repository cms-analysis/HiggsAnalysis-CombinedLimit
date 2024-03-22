from __future__ import absolute_import
import os
from WMCore.Configuration import Configuration


config = Configuration()

config.section_("General")
config.General.requestName = ""
# if (args.workArea != ''):
#   config.General.workArea = args.workArea

config.section_("JobType")
config.JobType.pluginName = "PrivateMC"
config.JobType.psetName = os.environ["CMSSW_BASE"] + "/src/HiggsAnalysis.CombinedLimit/scripts/do_nothing_cfg.py"
config.JobType.scriptExe = ""
config.JobType.inputFiles = [
    os.environ["CMSSW_BASE"] + "/src/HiggsAnalysis.CombinedLimit/scripts/FrameworkJobReport.xml",
    os.environ["CMSSW_BASE"] + "/src/HiggsAnalysis.CombinedLimit/scripts/copyRemoteWorkspace.sh",
    os.environ["CMSSW_BASE"] + "/bin/" + os.environ["SCRAM_ARCH"] + "/combine",
]
config.JobType.outputFiles = ["combine_output.tar"]
# config.JobType.maxMemoryMB = args.maxMemory

config.section_("Data")
config.Data.outputPrimaryDataset = "Combine"
config.Data.splitting = "EventBased"
config.Data.unitsPerJob = 1
config.Data.totalUnits = 1
config.Data.publication = False
config.Data.outputDatasetTag = ""

config.section_("User")

config.section_("Site")
config.Site.blacklist = ["T3_IT_Bologna", "T3_US_UMiss"]
config.Site.storageSite = "T2_CH_CERN"
