import FWCore.ParameterSet.Config as cms
'''
to run in CRAB3 I am not sure if parameters created in crabConfig_ step are relevant 
in range of this CMSSW parameter-set file, so I will have to use a copy-past code here
with no flexibility
'''
pName = "KSHORTS"
process = cms.Process(pName)

process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')# 80X_dataRun2_Prompt_v16 - doesn't work
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

studyroot = {
	'SUSYGluGluToHToTauTau_Full':
		{
			'isData': False,
			'fileName': 
			[
				'/store/mc/RunIIFall15DR76/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/0A56906C-73A7-E511-ADB2-00259073E512.root',
				'/store/mc/RunIIFall15DR76/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/0EB2F543-96A6-E511-8155-44A84225CDA4.root',
				'/store/mc/RunIIFall15DR76/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/38C656EE-7EA6-E511-8F28-44A84225CABC.root',
				'/store/mc/RunIIFall15DR76/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/644D005E-FCA7-E511-B160-00259073E52A.root',
				'/store/mc/RunIIFall15DR76/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/AE5DE371-93A6-E511-88F2-44A84225D0B7.root',
				'/store/mc/RunIIFall15DR76/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/C65AE49E-82A6-E511-91FF-20CF3027A561.root',
				'/store/mc/RunIIFall15DR76/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/CE0370C6-8CA6-E511-829E-00259073BB58.root',
				'/store/mc/RunIIFall15DR76/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/D4C0DDDF-D2A7-E511-9AFC-0CC47A4D767C.root',
				'/store/mc/RunIIFall15DR76/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/D6DA8FA0-A9A6-E511-83A0-44A84223FF3C.root'#,
				#'/store/mc/RunIIFall15DR76/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/F83430DF-9EA6-E511-8DE0-44A84225CABC.root' #migh be corrupted
			],
			'output_rootfile_name': "out_simAOD_SUSYGluGluToHToTauTau_Full.root"
		},
	'SUSYGluGluToHToTauTau':
		{
			'isData': False,
			'fileName': 'file:/.automount/home/home__home2/institut_3b/hlushchenko/Work/CMSSW_8_0_26_patch1/src/AOD_pi0_study/F83430DF-9EA6-E511-8DE0-44A84225CABC.root',
			'output_rootfile_name': "out_simAOD_SUSYGluGluToHToTauTau_siggleRemote.root"
		},
	'QCD_Pt_20toInf_MuEnrichedPt15': 
		{
			'isData':False, 
			'fileName': '/store/mc/RunIISpring16DR80/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/0C88808A-930D-E611-8DE7-B083FED0FFCF.root',#'file:/.automount/home/home__home2/institut_3b/hlushchenko/Work/CMSSW_8_0_26_patch1/src/AOD_pi0_study/0C88808A-930D-E611-8DE7-B083FED0FFCF.root',#'fileName': cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/0C88808A-930D-E611-8DE7-B083FED0FFCF.root'),
			'output_rootfile_name': "out_simAOD_QCD_Pt_20toInf_MuEnrichedPt15.root"
		},
	'QCD_Pt_20toInf_MuEnrichedPt15_Full':
		{
			'isData':False, 
			'fileName': 
			[
				'/store/mc/RunIISpring16DR80/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/0053436E-790D-E611-BA31-B499BAAC03BA.root', 
				'/store/mc/RunIISpring16DR80/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/00FEF7A6-3B0D-E611-8C73-00221983E092.root', 
				'/store/mc/RunIISpring16DR80/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/02B4D52C-980D-E611-8D65-001E4F1BC725.root', 
				'/store/mc/RunIISpring16DR80/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/02C16DFE-8F0D-E611-8CA1-44A84225C4EB.root', 
				'/store/mc/RunIISpring16DR80/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/04FB8C8C-C40D-E611-AF6C-00259057490C.root', 
				'/store/mc/RunIISpring16DR80/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/087A984B-AD0E-E611-9DA9-44A842CFC9F3.root', 
				'/store/mc/RunIISpring16DR80/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/08A931F3-570D-E611-ABF5-02163E012767.root', 
				'/store/mc/RunIISpring16DR80/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/0C5D0DD3-570D-E611-B4CB-02163E0137B9.root', 
				'/store/mc/RunIISpring16DR80/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/0C88808A-930D-E611-8DE7-B083FED0FFCF.root', 
				'file:/.automount/home/home__home2/institut_3b/hlushchenko/Work/CMSSW_8_0_26_patch1/src/AOD_pi0_study/0E2AE912-1C0E-E611-86FB-002590743042.root'
			],
			'output_rootfile_name': "out_simAOD_QCD_Pt_20toInf_MuEnrichedPt15_Full.root"
		},
	'JetHTdata': # JetHT AOD dataset=/JetHT*/*PromptReco-v2/AOD*  => dataset=/JetHT/Run2016B-PromptReco-v2/AOD
		{
			'isData': True, 
			'fileName': '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/150/00000/FC972EB3-D819-E611-94F9-02163E0134F4.root',#'file:/.automount/home/home__home2/institut_3b/hlushchenko/Work/CMSSW_8_0_20/src/Kappa/Skimming/higgsTauTau/16DA718F-DA19-E611-BCEE-02163E01376E.root',
			'output_rootfile_name': "out_AOD_JetHTdata_With_HPS_temp.root"
		}, # 10 out of 1315
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
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/086C3FED-DA19-E611-BC4E-02163E012611.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/087EA43F-DA19-E611-8610-02163E011DF8.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/0C69B99A-DA19-E611-9EAB-02163E011808.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/0E75A2A8-DA19-E611-8148-02163E011D55.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/10BEC6A7-DA19-E611-A373-02163E01456A.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/12BCFF67-EB19-E611-B2C4-02163E01417D.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/16DA718F-DA19-E611-BCEE-02163E01376E.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/18788F8B-DA19-E611-BA71-02163E01264D.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/1A8F7BFA-E819-E611-910F-02163E014683.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/1C0996EA-DA19-E611-827F-02163E014129.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/1E09B0CE-DA19-E611-8411-02163E0125FC.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/20E6C9D5-DA19-E611-888D-02163E012B1F.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/26D91FF0-DF19-E611-93F0-02163E013895.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/28177EED-DA19-E611-998B-02163E011FDE.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/28EE6368-DA19-E611-9F3F-02163E01456A.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/2A6A2E06-DB19-E611-9127-02163E014501.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/2ACD0DBF-DA19-E611-98F9-02163E011FAB.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/2C3C32B1-DA19-E611-BE23-02163E014151.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/2CCD4200-DB19-E611-8910-02163E01441D.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/32B30915-DC19-E611-9685-02163E0143CF.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/34A253B0-DD19-E611-834E-02163E014614.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/34D1220F-DB19-E611-8212-02163E01456A.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/36E4FFBD-DA19-E611-858B-02163E01376E.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/38406C52-DB19-E611-90DC-02163E0139C3.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/3A56B173-DA19-E611-8616-02163E01456A.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/3A8A29DA-DB19-E611-A5CF-02163E01348D.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/3EA5C90C-DB19-E611-8809-02163E0144F0.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/48F10CC8-DB19-E611-BA67-02163E014501.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/4AB035DA-DA19-E611-8C38-02163E013917.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/4CA6F7C0-DA19-E611-BA60-02163E0138B2.root',
				        # '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/4E24C919-DB19-E611-A7C4-02163E014540.root',
				        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/4E945EE7-DA19-E611-9084-02163E0139C8.root'
				        ],
	        'output_rootfile_name': "out_AOD_JetHTdata_FULL_With_HPS.root"
		}
}
filekey  = 'QCD_Pt_20toInf_MuEnrichedPt15_Full'
isData = studyroot[filekey]['isData']
#HERE
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(studyroot[filekey]['fileName']))
process.source = cms.Source("PoolSource",  fileNames = cms.untracked.vstring())
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))#HERE
process.MessageLogger = cms.Service("MessageLogger", destinations = cms.untracked.vstring("cout"), cout = cms.untracked.PSet(threshold = cms.untracked.string("ERROR")))
print("filekey: ", filekey)

print("HPSTracksLable is now hpsTracks")
process.load("HPStracks.HPStracksProducer.HPSTracks_cfi") # gives hpsTracks

print("SecondaryVerticesFromNewV0")
process.load("RecoVertex.V0Producer.generalV0Candidates_cfi")
process.SecondaryVerticesFromNewV0 = process.generalV0Candidates.clone( 
	# which beamSpot to reference
	beamSpot = cms.InputTag('offlineBeamSpot'),
	# reference primary vertex instead of beamSpot
	useVertex = cms.bool(True), # By def False
	# which vertex collection to use
	vertices = cms.InputTag('offlinePrimaryVertices'),
	# which TrackCollection to use for vertexing "generalV0Candidates","Lambda","RECO" || vector<reco::VertexCompositeCandidate>    "generalV0Candidates"       "Lambda"          "RECO"   
	trackRecoAlgorithm = cms.InputTag("hpsTracks", "HPSTracks", pName),#The standard track collection (label "generalTracks") is not saved in the MiniAOD event content.
	# which V0s to reconstruct
	doKShorts = cms.bool(True),
	doLambdas = cms.bool(False),

	# which vertex fitting algorithm to use
	# True -> KalmanVertexFitter (recommended)
	# False -> AdaptiveVertexFitter (not recommended)
	#vertexFitter = cms.bool(True),

	# use the refitted tracks returned from the KVF for V0Candidate kinematics
	# this is automatically set to False if using the AdaptiveVertexFitter
	#useRefTracks = cms.bool(True),

	#trackQualities = cms.vstring("loose"), # This is non-existing parameter which is prob will be stored in the output.root file
	# -- cuts on initial track collection --
	# Track normalized Chi2 <
	#tkChi2Cut = cms.double(10.),
	# Number of valid hits on track >=
	#tkNHitsCut = cms.int32(7),
	#tkNhitsCut = cms.int32(-9999999), #patch, variable broken in 72X when reading old file
	# Pt of track >
	#tkPtCut = cms.double(0.35),
	# Track impact parameter significance >
	####### tkIPSigXYCut = cms.double(-1),# was 2
	#tkIPSigZCut = cms.double(-1.),

	# -- cuts on the vertex --
	# Vertex chi2 <
	#vtxChi2Cut = cms.double(15.),
	# XY decay distance significance >
	vtxDecaySigXYCut = cms.double(10),#10 - DPG
	# XYZ decay distance significance >
	vtxDecaySigXYZCut = cms.double(-1.),

	# -- miscellaneous cuts --
	# POCA distance between tracks <
	#tkDCACut = cms.double(2.),
	# invariant mass of track pair - assuming both tracks are charged pions <
	#mPiPiCut = cms.double(0.6),
	# check if either track has a hit radially inside the vertex position minus this number times the sigma of the vertex fit
	# note: Set this to -1 to disable this cut, which MUST be done if you want to run V0Producer on the AOD track collection!
	innerHitPosCut = cms.double(-1.)
	# cos(angleXY) between x and p of V0 candidate >
	#cosThetaXYCut = cms.double(0.9998),
	# cos(angleXYZ) between x and p of V0 candidate >
	#cosThetaXYZCut = cms.double(-2.),

	# -- cuts on the V0 candidate mass --
	# V0 mass window +- pdg value
	#kShortMassCut = cms.double(0.07),
	#lambdaMassCut = cms.double(0.05)
)

available_v0 = {
					'new': 
					{
						"process_link": process.SecondaryVerticesFromNewV0,
						"collectionName": "SecondaryVerticesFromNewV0",
						"lable": "Kshort",
						"Process": pName,
						"newv0": True
					},
					'old': 
					{
						"process_link": 1,
						"collectionName": "generalV0Candidates",
						"lable": "Kshort",
						"Process": "RECO",
						"newv0": False
					}
				}
which_v0 = available_v0['new']

print("demo AOD_pi0")
process.demo = cms.EDAnalyzer('AOD_pi0',
	# # data, year, period, skim
	 IsData = cms.untracked.bool(isData),
	 OutFileName = cms.untracked.string(studyroot[filekey]['output_rootfile_name']),
	# Year = cms.untracked.uint32(2016),
	# Period = cms.untracked.string('Run2016B'),
	# Skim = cms.untracked.uint32(0),
	# # switches of collections
	GenParticles = cms.untracked.bool(False),
	RecPrimVertex = cms.untracked.bool(True),
	RecBeamSpot = cms.untracked.bool(True),
	RecTrack = cms.untracked.bool(False),
	RecPhoton = cms.untracked.bool(False),
	RecPiZero = cms.untracked.bool(True),
	RecSecVertex = cms.untracked.bool(True),
	RecJet = cms.untracked.bool(False),
	RecV0 = cms.untracked.bool(True),
	# # JEC
	# # 
	# # Collections
	# JetCollectionTag = cms.InputTag(""),
	 PFCandidateCollectionTag = cms.InputTag("particleFlow"),
	 GenParticleCollectionTag = cms.InputTag("genParticles"),
	 BeamSpotCollectionTag =  cms.InputTag("offlineBeamSpot"),
	 PVCollectionTag = cms.InputTag("offlinePrimaryVertices"),
	 KshortCollectionTag_stand = cms.InputTag("generalV0Candidates","Kshort","RECO"),#the lable RECO is only for the original collection
	 KshortCollectionTag = cms.InputTag(which_v0['collectionName'], which_v0['lable'], which_v0['Process']),#the lable RECO is only for the original collection
	 LambdaCollectionTag = cms.InputTag("generalV0Candidates","Lambda","RECO"),
	 TauPiZeroCollectionTag = cms.InputTag("hpsPFTauProducer","pizeros","RECO"),
	# SecVertexCollectionTag = cms.InputTag("inclusiveSecondaryVertices"),
	# tracks
	RecTrackPtMin = cms.untracked.double(1.0),
	RecTrackEtaMax = cms.untracked.double(2.4),
	RecTrackDxyMax = cms.untracked.double(1.0),
	RecTrackDzMax = cms.untracked.double(1.0),
	RecTrackNum = cms.untracked.int32(0),
	# pi0s
	RecPiZeroPtMin = cms.untracked.double(1.0),
	RecPiZeroEtaMax = cms.untracked.double(2.5),
	RecPizeroNum = cms.untracked.int32(0),
	# photons
	RecPhotonPtMin = cms.untracked.double(1.),
	RecPhotonEtaMax = cms.untracked.double(2.5),
	RecPhotonNum = cms.untracked.int32(0),
	# electrons
	RecElectronPtMin = cms.untracked.double(8.),
	RecElectronEtaMax = cms.untracked.double(2.6),
	RecElectronNum = cms.untracked.int32(0),
	# taus
	RecTauPtMin = cms.untracked.double(18),
	RecTauEtaMax = cms.untracked.double(2.4),                                      
	RecTauNum = cms.untracked.int32(0),
	# jets
	RecJetPtMin = cms.untracked.double(18.),
	RecJetEtaMax = cms.untracked.double(5.2),
	RecJetNum = cms.untracked.int32(0),
	# other
	Debug = cms.untracked.bool(False),
	Mute = cms.untracked.bool(True),
	Match_KsV0_to_HPS = cms.untracked.bool(True),
	HPSTrackTag = cms.InputTag("HPSTrackLable", "HPSTracks", "HPSTRACKS")#,
	#tkIPSigXYCut = cms.double(-1),# was 2
	#vtxDecaySigXYCut = cms.double(10)#10
)


#HERE
#this should probably go only together with grid
'''
process.output = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("path")),
    outputCommands = cms.untracked.vstring("drop *",
                                           "keep *_*_*_KSHORTS",
                                           "keep *_offlineBeamSpot_*_*",
                                           "keep *_offlinePrimaryVertices_*_*",
                                           "keep *_offlinePrimaryVerticesWithBS_*_*",
    ),
    fileName = cms.untracked.string("output.root"))
'''

# Run all three versions of the algorithm.
if which_v0['newv0']: 
	process.path = cms.Path( process.hpsTracks * which_v0["process_link"] * process.demo)#
else: 
	process.path = cms.Path(process.demo)
# Writer to a new file called output.root.  Save only the new K-shorts and the primary vertices (for later exercises).


process.TFileService = cms.Service("TFileService", fileName = cms.string(studyroot[filekey]['output_rootfile_name']) )

process.endpath = cms.EndPath(process.output)

