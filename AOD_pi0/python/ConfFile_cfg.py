import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
          'file:/nfs/dust/cms/user/glusheno/CMSSW_8_0_20/src/Kappa/Skimming/higgsTauTau/16DA718F-DA19-E611-BCEE-02163E01376E.root'        
    )
)

process.demo = cms.EDAnalyzer('AOD_pi0',
	# # data, year, period, skim
	# IsData = cms.untracked.bool(True),
	# Year = cms.untracked.uint32(2016),
	# Period = cms.untracked.string('Run2016B'),
	# Skim = cms.untracked.uint32(0),
	# # switches of collections
	# GenParticles = cms.untracked.bool(False),
	# RecPrimVertex = cms.untracked.bool(True),
	# RecBeamSpot = cms.untracked.bool(True),
	# RecTrack = cms.untracked.bool(False),
	# RecPhoton = cms.untracked.bool(False),
	# RecPiZero = cms.untracked.bool(True),
	# RecSecVertex = cms.untracked.bool(True),
	# RecJet = cms.untracked.bool(False),
	# RecV0 = cms.untracked.bool(True),
	# # JEC
	# # 
	# # Collections
	# JetCollectionTag = cms.InputTag(""),
	# PFCandidateCollectionTag = cms.InputTag("particleFlow"),
	# GenParticleCollectionTag = cms.InputTag("genParticles"),
	# BeamSpotCollectionTag =  cms.InputTag("offlineBeamSpot"),
	# PVCollectionTag = cms.InputTag("offlinePrimaryVertices"),
	 KshortCollectionTag = cms.InputTag("generalV0Candidates","Kshort","RECO"),
	 LambdaCollectionTag = cms.InputTag("generalV0Candidates","Lambda","RECO"),
	 TauPiZeroCollectionTag = cms.InputTag("hpsPFTauProducer","pizeros","RECO")
	# SecVertexCollectionTag = cms.InputTag("inclusiveSecondaryVertices"),
	# # tracks
	# RecTrackPtMin = cms.untracked.double(1.0),
	# RecTrackEtaMax = cms.untracked.double(2.4),
	# RecTrackDxyMax = cms.untracked.double(1.0),
	# RecTrackDzMax = cms.untracked.double(1.0),
	# RecTrackNum = cms.untracked.int32(0),
	# # pi0s
	# RecPiZeroPtMin = cms.untracked.double(1.0),
	# RecPiZeroEtaMax = cms.untracked.double(2.5),
	# RecPizeroNum = cms.untracked.int32(0),
	# # photons
	# RecPhotonPtMin = cms.untracked.double(1.),
	# RecPhotonEtaMax = cms.untracked.double(2.5),
	# RecPhotonNum = cms.untracked.int32(0),
	# # electrons
	# RecElectronPtMin = cms.untracked.double(8.),
	# RecElectronEtaMax = cms.untracked.double(2.6),
	# RecElectronNum = cms.untracked.int32(0),
	# # taus
	# RecTauPtMin = cms.untracked.double(18),
	# RecTauEtaMax = cms.untracked.double(2.4),                                      
	# RecTauNum = cms.untracked.int32(0),
	# # jets
	# RecJetPtMin = cms.untracked.double(18.),
	# RecJetEtaMax = cms.untracked.double(5.2),
	# RecJetNum = cms.untracked.int32(0)
)


process.p = cms.Path(process.demo)
