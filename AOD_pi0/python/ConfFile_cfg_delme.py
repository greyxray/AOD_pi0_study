# Use the tracks_and_vertices.root file as input.
import FWCore.ParameterSet.Config as cms
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag #from Configuration.AlCa.GlobalTag import GlobalTag
#from RecoVertex.V0Producer.generalV0Candidates_cfi import *

process = cms.Process("KSHORTS")

# Use the tracks_and_vertices.root file as input.
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring("file:/.automount/home/home__home2/institut_3b/hlushchenko/Work/CMSSW_8_0_26_patch1/src/AOD_pi0_study/0E2AE912-1C0E-E611-86FB-002590743042.root"))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
studyroot = {
	'SUSYMC': 
		{
			'isData':False, 
			'fileName': 'file:/.automount/home/home__home2/institut_3b/hlushchenko/Work/CMSSW_8_0_26_patch1/src/AOD_pi0_study/0E2AE912-1C0E-E611-86FB-002590743042.root',
			'output_rootfile_name': "out_simAOD_SUSYMC.root"
		},
	'JetHTdata': # JetHT AOD dataset=/JetHT*/*PromptReco-v2/AOD*  => dataset=/JetHT/Run2016B-PromptReco-v2/AOD
		{
			'isData': True, 
			'fileName': 'file:/.automount/home/home__home2/institut_3b/hlushchenko/Work/CMSSW_8_0_20/src/Kappa/Skimming/higgsTauTau/16DA718F-DA19-E611-BCEE-02163E01376E.root',
			'output_rootfile_name': "out_AOD_JetHTdata.root"
		},
	'kappaminidata':
		{
			'isData': True, 
			'fileName': 'file:/.automount/home/home__home2/institut_3b/hlushchenko/Work/CMSSW_8_0_20/src/Kappa/Skimming/higgsTauTau/miniAOD-prod_PAT.root',
			'output_rootfile_name': "out_miniAOD_kappaminidata.root"
		},
	'JetHTdata2': #ONLY ON NFS # JetHT AOD dataset=/JetHT*/*PromptReco-v2/AOD*  => dataset=/JetHT/Run2016B-PromptReco-v2/AOD
		{
			'isData': True, 
			'fileName': 'file:/nfs/dust/cms/user/glusheno/CMSSW_8_0_20/src/Kappa/Skimming/higgsTauTau/16DA718F-DA19-E611-BCEE-02163E01376E.root',
			'output_rootfile_name': "out_AOD_JetHTdata2.root"
		},
	'JetHTdataFull':
		{
			'isData': True,
			'fileName': ['/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/150/00000/FC972EB3-D819-E611-94F9-02163E0134F4.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/02C276F7-DA19-E611-8189-02163E0141DD.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/065F7CCE-DA19-E611-90AD-02163E0125FC.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/086C3FED-DA19-E611-BC4E-02163E012611.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/087EA43F-DA19-E611-8610-02163E011DF8.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/0C69B99A-DA19-E611-9EAB-02163E011808.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/0E75A2A8-DA19-E611-8148-02163E011D55.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/10BEC6A7-DA19-E611-A373-02163E01456A.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/12BCFF67-EB19-E611-B2C4-02163E01417D.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/16DA718F-DA19-E611-BCEE-02163E01376E.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/18788F8B-DA19-E611-BA71-02163E01264D.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/1A8F7BFA-E819-E611-910F-02163E014683.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/1C0996EA-DA19-E611-827F-02163E014129.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/1E09B0CE-DA19-E611-8411-02163E0125FC.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/20E6C9D5-DA19-E611-888D-02163E012B1F.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/26D91FF0-DF19-E611-93F0-02163E013895.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/28177EED-DA19-E611-998B-02163E011FDE.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/28EE6368-DA19-E611-9F3F-02163E01456A.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/2A6A2E06-DB19-E611-9127-02163E014501.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/2ACD0DBF-DA19-E611-98F9-02163E011FAB.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/2C3C32B1-DA19-E611-BE23-02163E014151.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/2CCD4200-DB19-E611-8910-02163E01441D.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/32B30915-DC19-E611-9685-02163E0143CF.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/34A253B0-DD19-E611-834E-02163E014614.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/34D1220F-DB19-E611-8212-02163E01456A.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/36E4FFBD-DA19-E611-858B-02163E01376E.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/38406C52-DB19-E611-90DC-02163E0139C3.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/3A56B173-DA19-E611-8616-02163E01456A.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/3A8A29DA-DB19-E611-A5CF-02163E01348D.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/3EA5C90C-DB19-E611-8809-02163E0144F0.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/48F10CC8-DB19-E611-BA67-02163E014501.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/4AB035DA-DA19-E611-8C38-02163E013917.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/4CA6F7C0-DA19-E611-BA60-02163E0138B2.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/4E24C919-DB19-E611-A7C4-02163E014540.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/4E945EE7-DA19-E611-9084-02163E0139C8.root'],
	        'output_rootfile_name': "out_AOD_JetHTdata_FULL.root"
		}
}
filekey  = 'JetHTdata'
isData = studyroot[filekey]['isData']
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(studyroot[filekey]['fileName']))
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)# -1 for all events
)


# Suppress messages that are less important than ERRORs.
#process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring("cout"),
    cout = cms.untracked.PSet(threshold = cms.untracked.string("INFO")))

# Load part of the CMSSW reconstruction sequence to make vertexing possible.
# We'll need the CMS geometry and magnetic field to follow the true, non-helical shapes of tracks through the detector.
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')# 80X_dataRun2_Prompt_v16 - doesn't work
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

# Copy most of the vertex producer's parameters, but accept tracks with progressively more strict quality.
process.load("RecoVertex.V0Producer.generalV0Candidates_cfi")
process.load("HPStracks.HPStracksProducer.HPSTracks_cfi") # gives hpsTracks
