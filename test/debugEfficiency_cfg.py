
import FWCore.ParameterSet.Config as cms

process = cms.Process("debugEfficiencyLoss")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100000)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/home/veelken/Phase2HLT/CMSSW_11_1_0_pre6/src/HLTTrigger/Phase2HLTPFTaus/test/selEvents_debugEfficiencyLoss_RECO_2020Jul03.root'
    ),
    ##eventsToProcess = cms.untracked.VEventRange(
    ##    '1:11:1454'
    ##) 
)

inputFilePath = '/hdfs/cms/store/user/rdewanje/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5/HLTConfig_w_offlineVtxCollection_woHGCal_VBFHTT_Phase2HLTTDRWinter20_PU200_CMSSW_11_1_0_pre8_isoFix2/200707_101538/' # without HGCal energy regression
outputFileName = "selEvents_debugEfficiencyLoss_absEta1p6to2p4_isoFix2_2020Jul07.txt"

#--------------------------------------------------------------------------------
# set input files
from HLTTrigger.TallinnHLTPFTauAnalyzer.tools.jobTools import getInputFileNames
print("Searching for input files in path = '%s'" % inputFilePath)
inputFileNames = getInputFileNames(inputFilePath)
print("Found %i input files." % len(inputFileNames))
process.source.fileNames = cms.untracked.vstring(inputFileNames)
#--------------------------------------------------------------------------------

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.analysisSequence = cms.Sequence()

process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.analysisSequence += process.genParticles

process.load("PhysicsTools.JetMCAlgos.TauGenJets_cfi")
process.tauGenJets.GenParticles = cms.InputTag('genParticles')
process.analysisSequence += process.tauGenJets

process.load("PhysicsTools.JetMCAlgos.TauGenJetsDecayModeSelectorAllHadrons_cfi")
process.analysisSequence += process.tauGenJetsSelectorAllHadrons

process.selectedGenHadTaus = cms.EDFilter("GenJetSelector",
  src = cms.InputTag('tauGenJetsSelectorAllHadrons'),
  cut = cms.string('pt > 20. & abs(eta) > 2.0 & abs(eta) < 2.4'),
  filter = cms.bool(False)
)
process.analysisSequence += process.selectedGenHadTaus

process.hltPFTausByROI = cms.EDFilter("PFTauAntiOverlapSelector",
  src = cms.InputTag('hltSelectedHpsPFTaus8HitsMaxDeltaZWithOfflineVertices'),
  srcNotToBeFiltered = cms.VInputTag('selectedGenHadTaus'),
  dRmin = cms.double(0.3),
  invert = cms.bool(True),
  filter = cms.bool(True)                                                          
)
process.analysisSequence += process.hltPFTausByROI

from RecoTauTag.RecoTau.PFRecoTauQualityCuts_cfi import PFTauQualityCuts
hltQualityCuts = PFTauQualityCuts.clone()
hltQualityCuts.signalQualityCuts.minTrackPt = cms.double(0.9)
hltQualityCuts.isolationQualityCuts.minTrackPt = cms.double(0.9)
hltQualityCuts.isolationQualityCuts.maxDeltaZ = cms.double(0.2)
hltQualityCuts.isolationQualityCuts.maxDeltaZToLeadTrack = cms.double(-1.)
hltQualityCuts.isolationQualityCuts.minTrackHits = cms.uint32(8)
hltQualityCuts.primaryVertexSrc = cms.InputTag('offlinePrimaryVertices')

from RecoTauTag.RecoTau.TauDiscriminatorTools import noPrediscriminants
process.hltPFTauChargedIsoPtSumByROI = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
  PFTauProducer = cms.InputTag('hltPFTausByROI'),
  particleFlowSrc = cms.InputTag('particleFlowTmp'),
  vertexSrc = cms.InputTag('offlinePrimaryVertices'),
  qualityCuts = hltQualityCuts,
  Prediscriminants = noPrediscriminants,
  ApplyDiscriminationByTrackerIsolation = cms.bool(True),
  ApplyDiscriminationByECALIsolation = cms.bool(False),
  ApplyDiscriminationByWeightedECALIsolation = cms.bool(False),
  enableHGCalWorkaround = cms.bool(True),
  WeightECALIsolation = cms.double(0.),
  minTauPtForNoIso = cms.double(-99.),
  applyOccupancyCut = cms.bool(False),
  maximumOccupancy = cms.uint32(0),
  applySumPtCut = cms.bool(False),
  maximumSumPtCut = cms.double(-1.),
  applyRelativeSumPtCut = cms.bool(False),
  relativeSumPtCut = cms.double(-1.),
  relativeSumPtOffset = cms.double(0.),
  storeRawOccupancy = cms.bool(False),
  storeRawSumPt = cms.bool(True),
  storeRawPUsumPt = cms.bool(False),
  storeRawFootprintCorrection = cms.bool(False),
  storeRawPhotonSumPt_outsideSignalCone = cms.bool(False),
  customOuterCone = cms.double(-1.),
  applyPhotonPtSumOutsideSignalConeCut = cms.bool(False),
  maxAbsPhotonSumPt_outsideSignalCone = cms.double(1.e+9),
  maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
  applyFootprintCorrection = cms.bool(False),
  footprintCorrections = cms.VPSet(),
  applyDeltaBetaCorrection = cms.bool(False),
  deltaBetaPUTrackPtCutOverride = cms.bool(False),
  deltaBetaPUTrackPtCutOverride_val = cms.double(0.5),
  isoConeSizeForDeltaBeta = cms.double(0.5),
  deltaBetaFactor = cms.string("0.38"),
  applyRhoCorrection = cms.bool(False),
  rhoProducer = cms.InputTag("NotUsed"),
  rhoConeSize = cms.double(0.357),
  rhoUEOffsetCorrection = cms.double(0.),
  UseAllPFCandsForWeights = cms.bool(False),
  verbosity = cms.int32(0)
)
process.analysisSequence += process.hltPFTauChargedIsoPtSumByROI

process.hltPFTausFailingTrigger = cms.EDProducer("MyPFTauSelector",
  src = cms.InputTag('hltPFTausByROI'),
  src_sumChargedIso = cms.InputTag('hltPFTauChargedIsoPtSumByROI'),
  min_pt = cms.double(20.),
  max_pt = cms.double(-1.),
  min_absEta = cms.double(-1.),
  max_absEta = cms.double(2.4),
  min_leadTrackPt = cms.double(5.),
  max_leadTrackPt = cms.double(-1.),
  min_relChargedIso = cms.double(-1.),
  max_relChargedIso = cms.double(0.05),
  min_absChargedIso = cms.double(-1.),
  max_absChargedIso = cms.double(-1.),
  invert = cms.bool(True)
  ##invert = cms.bool(False)
)
process.analysisSequence += process.hltPFTausFailingTrigger

process.hltPFTausFailingTriggerFilter = cms.EDFilter("CandViewCountFilter",
  src = cms.InputTag('hltPFTausFailingTrigger'),
  minNumber = cms.uint32(1)
)
process.analysisSequence += process.hltPFTausFailingTriggerFilter

process.genTaus = cms.EDFilter("GenParticleSelector",
  src = cms.InputTag('genParticles'),
  cut = cms.string('abs(pdgId) = 15 & status = 2 & pt > 20. & abs(eta) > 2.0 & abs(eta) < 2.4'),
  stableOnly = cms.bool(True),
  filter = cms.bool(False)
)
process.analysisSequence += process.genTaus

process.genTausByROI = cms.EDFilter("GenParticleAntiOverlapSelector",
  src = cms.InputTag('genTaus'),
  srcNotToBeFiltered = cms.VInputTag('hltPFTausFailingTrigger'),
  dRmin = cms.double(0.3),
  invert = cms.bool(True),
  filter = cms.bool(False)                                                          
)
process.analysisSequence += process.genTausByROI

process.dumpGenTaus = cms.EDAnalyzer("DumpGenParticles",
  src = cms.InputTag('genTausByROI')
)
process.analysisSequence += process.dumpGenTaus

process.genHadTausByROI = cms.EDFilter("GenJetAntiOverlapSelector",
  src = cms.InputTag('selectedGenHadTaus'),
  srcNotToBeFiltered = cms.VInputTag('hltPFTausFailingTrigger'),
  dRmin = cms.double(0.3),
  invert = cms.bool(True),
  filter = cms.bool(False)                                                          
)
process.analysisSequence += process.genHadTausByROI

process.dumpGenHadTaus = cms.EDAnalyzer("DumpGenTaus",
  src = cms.InputTag('genHadTausByROI')
)
process.analysisSequence += process.dumpGenHadTaus

process.dumpOnlinePFTaus = cms.EDAnalyzer("DumpRecoPFTaus",
  src = cms.InputTag('hltPFTausByROI'),
  src_sumChargedIso = cms.InputTag('hltPFTauChargedIsoPtSumByROI'),
  src_discriminators = cms.VInputTag()
)
process.analysisSequence += process.dumpOnlinePFTaus

##process.dumpOnlinePFTaus2 = cms.EDAnalyzer("DumpRecoPFTaus",
##  src = cms.InputTag('hltSelectedHpsPFTaus8HitsMaxDeltaZWithOfflineVertices'),
##  src_sumChargedIso = cms.InputTag('hltPFTauChargedIsoPtSum8HitsMaxDeltaZWithOfflineVertices'),
##  src_discriminators = cms.VInputTag()
##)
##process.analysisSequence += process.dumpOnlinePFTaus2

process.stableGenParticles = cms.EDFilter("GenParticleSelector",
  src = cms.InputTag('genParticles'),
  cut = cms.string('status = 1 & abs(pdgId) != 12 & abs(pdgId) != 14 & abs(pdgId) != 16'),
  stableOnly = cms.bool(True),
  filter = cms.bool(False)
)
process.analysisSequence += process.stableGenParticles

process.stableGenParticlesByROI = cms.EDFilter("GenParticleAntiOverlapSelector",
  src = cms.InputTag('stableGenParticles'),
  srcNotToBeFiltered = cms.VInputTag('hltPFTausFailingTrigger'),
  dRmin = cms.double(0.3),
  invert = cms.bool(True),
  filter = cms.bool(False)                                                          
)
process.analysisSequence += process.stableGenParticlesByROI

process.dumpGenParticles = cms.EDAnalyzer("DumpGenParticles",
  src = cms.InputTag('stableGenParticlesByROI')
)
process.analysisSequence += process.dumpGenParticles

# CV: select subset of reco::Track and reco::PFCandidate objects within dR < 0.8 cones around generator-level hadronic tau decays
#     in order to reduce computing time
process.hltPFCandidatesByROI = cms.EDFilter("PFCandidateAntiOverlapSelector",
  src = cms.InputTag('particleFlowTmp'),
  srcNotToBeFiltered = cms.VInputTag('hltPFTausFailingTrigger'),
  dRmin = cms.double(0.3),
  invert = cms.bool(True),
  filter = cms.bool(False)
)
process.analysisSequence += process.hltPFCandidatesByROI

process.dumpOnlinePFCandidates = cms.EDAnalyzer("DumpRecoPFCandidates",
  src = cms.InputTag('hltPFCandidatesByROI'),
  min_pt = cms.double(-1.),
  max_absEta = cms.double(3.0)
)
process.analysisSequence += process.dumpOnlinePFCandidates

process.hltTracksByROI = cms.EDProducer("RecoTrackAntiOverlapSelector",
  src = cms.InputTag('generalTracks'),
  srcNotToBeFiltered = cms.VInputTag('hltPFTausFailingTrigger'),
  dRmin = cms.double(0.3),
  invert = cms.bool(True)
)
process.analysisSequence += process.hltTracksByROI

process.dumpOnlineTracks = cms.EDAnalyzer("DumpRecoTracks",
  src = cms.InputTag('hltTracksByROI')
)
process.analysisSequence += process.dumpOnlineTracks

process.dumpRho = cms.EDAnalyzer("DumpDouble",
  src = cms.InputTag('hltKT6PFJets:rho')
)
process.analysisSequence += process.dumpRho

process.runLumiSectionEventNumberAnalyzer = cms.EDAnalyzer("RunLumiSectionEventNumberAnalyzer",
  output = cms.string(outputFileName),
  separator = cms.string(":")
)
process.analysisSequence += process.runLumiSectionEventNumberAnalyzer

process.p = cms.Path(process.analysisSequence)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)