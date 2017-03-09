# Use the tracks_and_vertices.root file as input.
import FWCore.ParameterSet.Config as cms

process = cms.Process("KSHORTS")

# Use the tracks_and_vertices.root file as input.
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring("file:/.automount/home/home__home2/institut_3b/hlushchenko/Work/CMSSW_8_0_26_patch1/src/AOD_pi0_study/0E2AE912-1C0E-E611-86FB-002590743042.root"))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

# Suppress messages that are less important than ERRORs.
process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring("cout"),
    cout = cms.untracked.PSet(threshold = cms.untracked.string("INFO")))

# Load part of the CMSSW reconstruction sequence to make vertexing possible.
# We'll need the CMS geometry and magnetic field to follow the true, non-helical shapes of tracks through the detector.
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

# Copy most of the vertex producer's parameters, but accept tracks with progressively more strict quality.
process.load("RecoVertex.V0Producer.generalV0Candidates_cfi")

# loose
process.SecondaryVerticesFromLooseTracks = process.generalV0Candidates.clone(
    trackRecoAlgorithm = cms.InputTag("generalTracks"),
    selectKshorts = cms.bool(True),
    selectLambdas = cms.bool(False),
    trackQualities = cms.vstring("loose"),
    innerHitPosCut = cms.double(-1.),
    tkNhitsCut = cms.int32(-9999999) #patch, variable broken in 72X when reading old file
    )

# tight
process.SecondaryVerticesFromTightTracks = process.SecondaryVerticesFromLooseTracks.clone(
    trackQualities = cms.vstring("tight"),
    )

# highPurity
process.SecondaryVerticesFromHighPurityTracks = process.SecondaryVerticesFromLooseTracks.clone(
    trackQualities = cms.vstring("highPurity"),
    )

# Run all three versions of the algorithm.
process.path = cms.Path(process.SecondaryVerticesFromLooseTracks *
                        process.SecondaryVerticesFromTightTracks *
                        process.SecondaryVerticesFromHighPurityTracks)

# Writer to a new file called output.root.  Save only the new K-shorts and the primary vertices (for later exercises).
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
