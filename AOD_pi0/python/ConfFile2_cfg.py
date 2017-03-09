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

# loose
process.SecondaryVerticesFromLooseTracks = process.generalV0Candidates.clone( 
     # which beamSpot to reference
   beamSpot = cms.InputTag('offlineBeamSpot'),
   # reference primary vertex instead of beamSpot
   useVertex = cms.bool(False),
   # which vertex collection to use
   vertices = cms.InputTag('offlinePrimaryVertices'),
   # which TrackCollection to use for vertexing
    trackRecoAlgorithm = cms.InputTag("generalTracks"),#The standard track collection (label "generalTracks") is not saved in the MiniAOD event content.
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
   tkIPSigXYCut = cms.double(-1),#2
   #tkIPSigZCut = cms.double(-1.),

   # -- cuts on the vertex --
   # Vertex chi2 <
   #vtxChi2Cut = cms.double(15.),
   # XY decay distance significance >
   vtxDecaySigXYCut = cms.double(-1),#10
   # XYZ decay distance significance >
   #vtxDecaySigXYZCut = cms.double(-1.),

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

# tight
process.SecondaryVerticesFromTightTracks = process.SecondaryVerticesFromLooseTracks.clone(
    trackQualities = cms.vstring("tight"),
    )

# highPurity
process.SecondaryVerticesFromHighPurityTracks = process.SecondaryVerticesFromLooseTracks.clone(
    trackQualities = cms.vstring("highPurity"),
    )

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
	# BeamSpotCollectionTag =  cms.InputTag("offlineBeamSpot"),
	 PVCollectionTag = cms.InputTag("offlinePrimaryVertices"),
	 KshortCollectionTag = cms.InputTag("SecondaryVerticesFromHighPurityTracks","Kshort","KSHORTS"),#the lable RECO is only for the original collection
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
	Debug = cms.untracked.bool(True),
	Mute = cms.untracked.bool(False),
	Match_KsV0_to_HPS = cms.untracked.bool(True)
)

# Run all three versions of the algorithm.
process.path = cms.Path(process.SecondaryVerticesFromHighPurityTracks * process.demo)

# Writer to a new file called output.root.  Save only the new K-shorts and the primary vertices (for later exercises).
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
process.endpath = cms.EndPath(process.output)
'''
