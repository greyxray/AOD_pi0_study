from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from WMCore.DataStructs.LumiList import LumiList

config = config()

config.General.requestName = 'QCD_Pt_20toInf_MuEnrichedPt15_HPS_matched_byReference_full' #
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True # whether or not to transfer the output files to the storage site.
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis' # Specifies if this task is running an analysis ('Analysis') on an existing dataset or is running MC event generation ('PrivateMC').
config.JobType.psetName = 'AOD_pi0_study/AOD_pi0/python/ConfFileWithHPSTracks_cfg.py' #ConfFileWithHPSTracks_cfg_CRAB3 'pset_tutorial_analysis.py'

config.Data.inputDataset = '/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'# Analysis job type only supports the following splitting algorithms (plus 'Automatic' as of CMSSW_7_2_X): ['FileBased', 'LumiBased', 'EventAwareLumiBased', 'Automatic'].
config.Data.unitsPerJob = 1
# No config.Data.lumiMask for MC
# No config.Data.runRange for MC
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'CRAB3_tutorial_May2015_Data_analysis_QCD_Pt_20toInf_MuEnrichedPt15_HPS_matched_byReference_full'

config.Site.storageSite = 'T2_DE_RWTH'