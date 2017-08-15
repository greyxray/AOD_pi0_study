from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from WMCore.DataStructs.LumiList import LumiList

config = config()

config.General.requestName = 'JetHTdata_HPS_matched_byReference_full' #
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True # whether or not to transfer the output files to the storage site.
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis' # Specifies if this task is running an analysis ('Analysis') on an existing dataset or is running MC event generation ('PrivateMC').
config.JobType.psetName = 'AOD_pi0_study/AOD_pi0/python/ConfFileWithHPSTracks_cfg.py' #ConfFileWithHPSTracks_cfg_CRAB3 'pset_tutorial_analysis.py'

config.Data.inputDataset = '/JetHT/Run2016B-PromptReco-v2/AOD'#'/SingleMu/Run2012B-13Jul2012-v1/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 20
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt'#'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt'
# config.Data.runRange = '273158-276811'# from lumimask 273158-276811 actually lumimask is enough
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'CRAB3_tutorial_May2015_Data_analysis_JetHTdata_HPS_matched_byReference_full'

config.Site.storageSite = 'T2_DE_RWTH'