# For STXStoEFT models

# Python script to hold dictionary of processes in STXS Stages

#using convention in https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsTemplateCrossSection

#stage0_procs = {}

stage1_procs = { 
  'ggH':['ggH_VBFTOPO_JET3VETO', 'ggH_VBFTOPO_JET3',
            'ggH_0J',
            'ggH_1J_PTH_0_60', 'ggH_1J_PTH_60_120', 'ggH_1J_PTH_120_200', 'ggH_1J_PTH_GT200',
            'ggH_GE2J_PTH_0_60', 'ggH_GE2J_PTH_60_120', 'ggH_GE2J_PTH_120_200', 'ggH_GE2J_PTH_GT200'],
  'qqH':['qqH_VBFTOPO_JET3VETO' , 'qqH_VBFTOPO_JET3',
            'qqH_VH2JET', 'qqH_REST', 'qqH_PTJET1_GT200'],
  'WH_lep':['WH_lep_PTV_0_150', 'WH_lep_PTV_150_250_0J', 'WH_lep_PTV_150_250_GE1J', 'WH_lep_PTV_GT250'],
  'ZH_lep':['ZH_lep_PTV_0_150', 'ZH_lep_PTV_150_250_0J', 'ZH_lep_PTV_150_250_GE1J', 'ZH_lep_PTV_GT250'],
  #'ggZH_lep':['ggZH_lep_PTV_0_150', 'ggZH_lep_PTV_GT150_0J', 'ggZH_lep_PTV_GT150_GE1J'],
  'VH_had':['WH_had_VBFTOPO_JET3VETO', 'WH_had_VBFTOPO_JET3',
            'WH_had_VH2JET', 'WH_had_REST', 'WH_had_PTJET1_GT200',
            'ZH_had_VBFTOPO_JET3VETO', 'ZH_had_VBFTOPO_JET3',
            'ZH_had_VH2JET', 'ZH_had_REST', 'ZH_had_PTJET1_GT200'],
  'ttH':['ttH'],
  'other':['tHq','tHW','bbH']
}

# stage1p1_procs = {}
