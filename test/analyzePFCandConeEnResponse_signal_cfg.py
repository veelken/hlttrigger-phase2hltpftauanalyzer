
import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzePFCandConeEnResponse")

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
    )
)

##inputFilePath = '/hdfs/cms/store/user/rdewanje/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5/HLTConfig_w_offlineVtxCollection_HGCalFix_VBFHTT_Phase2HLTTDRWinter20_PU200_CMSSW_11_1_0_pre8/200627_142633/' # with HGCal energy regression
inputFilePath = '/hdfs/cms/store/user/rdewanje/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5/HLTConfig_w_offlineVtxCollection_woHGCal_VBFHTT_Phase2HLTTDRWinter20_PU200_CMSSW_11_1_0_pre8_isoFix2/200707_101538/' # without HGCal energy regression
outputFileName = "analyzePFCandConeEnResponse_signal_2020Jul07.root"

#--------------------------------------------------------------------------------
# set input files
from HLTrigger.TallinnHLTPFTauAnalyzer.tools.jobTools import getInputFileNames
print("Searching for input files in path = '%s'" % inputFilePath)
inputFileNames = getInputFileNames(inputFilePath)
print("Found %i input files." % len(inputFileNames))
process.source.fileNames = cms.untracked.vstring(inputFileNames)
#--------------------------------------------------------------------------------

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.analysisSequence = cms.Sequence()

#--------------------------------------------------------------------------------
process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.analysisSequence += process.genParticles

process.load("PhysicsTools.JetMCAlgos.TauGenJets_cfi")
process.tauGenJets.GenParticles = cms.InputTag('genParticles')
process.analysisSequence += process.tauGenJets

process.load("PhysicsTools.JetMCAlgos.TauGenJetsDecayModeSelectorAllHadrons_cfi")
process.analysisSequence += process.tauGenJetsSelectorAllHadrons

process.selectedGenHadTaus = cms.EDFilter("GenJetSelector",
  src = cms.InputTag('tauGenJetsSelectorAllHadrons'),
  cut = cms.string('pt > 20. & abs(eta) < 2.4'),
  filter = cms.bool(False)
)
process.analysisSequence += process.selectedGenHadTaus

# CV: select subset of reco::Track and reco::PFCandidate objects within dR < 0.8 cones around generator-level hadronic tau decays
#     in order to reduce computing time
process.generalTracksByROI = cms.EDProducer("RecoTrackAntiOverlapSelector",
  src = cms.InputTag('generalTracks'),
  srcNotToBeFiltered = cms.VInputTag('selectedGenHadTaus'),
  dRmin = cms.double(0.5),
  invert = cms.bool(True)
)
process.analysisSequence += process.generalTracksByROI

process.particleFlowTmpByROI = cms.EDFilter("PFCandidateAntiOverlapSelector",
  src = cms.InputTag('particleFlowTmp'),
  srcNotToBeFiltered = cms.VInputTag('selectedGenHadTaus'),
  dRmin = cms.double(0.5),
  invert = cms.bool(True),
  filter = cms.bool(False)
)
process.analysisSequence += process.particleFlowTmpByROI

process.load("RecoMET.Configuration.GenMETParticles_cff")
process.analysisSequence += process.genParticlesForMETAllVisible

process.load("RecoMET.Configuration.RecoGenMET_cff")
process.analysisSequence += process.genMetTrue

process.genParticlesForMETAllVisibleByROI = cms.EDFilter("CandidateViewAntiOverlapSelector",
  src = cms.InputTag('genParticlesForMETAllVisible'),
  srcNotToBeFiltered = cms.VInputTag('selectedGenHadTaus'),
  dRmin = cms.double(0.5),
  invert = cms.bool(True),
  filter = cms.bool(False)                                                          
)
process.analysisSequence += process.genParticlesForMETAllVisibleByROI

process.analyzePFCandConeEnResponse = cms.EDAnalyzer("RecoPFCandidateConeEnResponseAnalyzer",
  srcPFCandidates = cms.InputTag('particleFlowTmpByROI'),
  srcTracks = cms.InputTag('generalTracksByROI'),
  srcGenParticles = cms.InputTag('genParticlesForMETAllVisibleByROI'),
  srcGenHadTaus = cms.InputTag('selectedGenHadTaus'),
  dRcones = cms.vdouble(0.15, 0.25, 0.35),
  dqmDirectory = cms.string("RecoPFCandidateConeEnResponseAnalyzer")
)
process.analysisSequence += process.analyzePFCandConeEnResponse

process.analyzePFMEt = cms.EDAnalyzer("RecoMEtResolutionAnalyzer",
  srcRecMEt = cms.InputTag('hltPFMET'),
  srcGenMEt = cms.InputTag('genMetTrue'),
  dqmDirectory = cms.string("RecoMEtResolutionAnalyzer/hltPFMET"),
)
process.analysisSequence += process.analyzePFMEt

process.analyzePFMEtTypeOne = cms.EDAnalyzer("RecoMEtResolutionAnalyzer",
  srcRecMEt = cms.InputTag('hltPFMETTypeOne'),
  srcGenMEt = cms.InputTag('genMetTrue'),
  dqmDirectory = cms.string("RecoMEtResolutionAnalyzer/hltPFMETTypeOne"),
)
process.analysisSequence += process.analyzePFMEtTypeOne

process.analyzePuppiMEt = cms.EDAnalyzer("RecoMEtResolutionAnalyzer",
  srcRecMEt = cms.InputTag('hltPuppiMET'),
  srcGenMEt = cms.InputTag('genMetTrue'),
  dqmDirectory = cms.string("RecoMEtResolutionAnalyzer/hltPuppiMET"),
)
process.analysisSequence += process.analyzePuppiMEt

process.analyzePuppiMEtTypeOne = cms.EDAnalyzer("RecoMEtResolutionAnalyzer",
  srcRecMEt = cms.InputTag('hltPuppiMETTypeOne'),
  srcGenMEt = cms.InputTag('genMetTrue'),
  dqmDirectory = cms.string("RecoMEtResolutionAnalyzer/hltPuppiMETTypeOne"),
)
process.analysisSequence += process.analyzePFMEtTypeOne
#--------------------------------------------------------------------------------

process.load("DQMServices.Core.DQMStore_cfi")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string(outputFileName)
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
