
import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzePFTausSignal")

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
        'file:/home/veelken/Phase2HLT/CMSSW_11_1_0_pre6/src/HLTTrigger/Phase2HLTPFTaus/test/step3_RAW2DIGI_RECO.root'
    )
)

inputFilePath = '/hdfs/cms/store/user/rdewanje/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5/HLTConfig_w_offlineVtxCollection_HGCalFix_VBFHTT_Phase2HLTTDRWinter20_PU200_CMSSW_11_1_0_pre8/200627_142633/'
processName = "qqH_htt"
srcVertices = 'offlinePrimaryVertices'
#srcVertices = 'hltPhase2PixelVertices'
outputFileName = "analyzePFTauResponse_signal_2020Jul01v2.root"

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

#--------------------------------------------------------------------------------
process.load("PhysicsTools.JetMCAlgos.TauGenJets_cfi")
process.tauGenJets.GenParticles = cms.InputTag('prunedGenParticles')
process.analysisSequence += process.tauGenJets

process.load("PhysicsTools.JetMCAlgos.TauGenJetsDecayModeSelectorAllHadrons_cfi")
process.analysisSequence += process.tauGenJetsSelectorAllHadrons

process.selectedGenHadTaus = cms.EDFilter("GenJetSelector",
  src = cms.InputTag('tauGenJetsSelectorAllHadrons'),
  cut = cms.string('pt > 20. & abs(eta) < 2.4'),
  filter = cms.bool(False)
)
process.analysisSequence += process.selectedGenHadTaus

process.genMatchedOfflinePFTaus = cms.EDFilter("PATTauAntiOverlapSelector",
  src = cms.InputTag('slimmedTaus'),
  srcNotToBeFiltered = cms.VInputTag('selectedGenHadTaus'),
  dRmin = cms.double(0.3),
  invert = cms.bool(True),
  filter = cms.bool(True)                                                          
)
process.analysisSequence += process.genMatchedOfflinePFTaus

process.selectedOfflinePFTaus = cms.EDFilter("PATTauSelector",
  src = cms.InputTag('genMatchedOfflinePFTaus'),
  cut = cms.string("pt > 20. & abs(eta) < 2.4 & tauID('decayModeFinding') > 0.5 & tauID('byLooseCombinedIsolationDeltaBetaCorr3Hits') > 0.5")
)
process.analysisSequence += process.selectedOfflinePFTaus

process.selectedOfflinePFTauFilter = cms.EDFilter("CandViewCountFilter",
  src = cms.InputTag('selectedOfflinePFTaus'),
  minNumber = cms.uint32(1)
)
process.analysisSequence += process.selectedOfflinePFTauFilter

process.offlineMatchedGenHadTaus = cms.EDFilter("GenJetAntiOverlapSelector",
  src = cms.InputTag('selectedGenHadTaus'),
  srcNotToBeFiltered = cms.VInputTag('selectedOfflinePFTaus'),
  dRmin = cms.double(0.3),
  invert = cms.bool(True),
  filter = cms.bool(True)                                                          
)
process.analysisSequence += process.offlineMatchedGenHadTaus

process.genVertex = cms.EDProducer("GenVertexProducer",
  src = cms.InputTag('prunedGenParticles'),
  pdgIds = cms.vint32(-15, +15) # CV: use { -15, +15 } for signal, empty list for background                                    
) 
process.analysisSequence += process.genVertex

for algorithm in [ "hps", "shrinking-cone" ]:

  pfTauLabel = None
  if algorithm == "shrinking-cone":
    pfTauLabel = "PFTau"
  elif algorithm == "hps":
    pfTauLabel = "HpsPFTau"
  else:
    raise ValueError("Invalid parameter algorithm = '%s' !!" % algorithm)

  for isolation_maxDeltaZOption in [ "primaryVertex", "leadTrack" ]:
    ##for isolation_minTrackHits in [ 3, 5, 8 ]:  
    for isolation_minTrackHits in [ 8 ]: 

      suffix = "%iHits" % isolation_minTrackHits
      if isolation_maxDeltaZOption == "primaryVertex":
        suffix += "MaxDeltaZ"
      elif isolation_maxDeltaZOption == "leadTrack":
        suffix += "MaxDeltaZToLeadTrack"
      else:
        raise ValueError("Invalid parameter isolation_maxDeltaZOption = '%s' !!" % isolation_maxDeltaZOption)
      if srcVertices == "offlinePrimaryVertices":
        suffix += "WithOfflineVertices"
      elif srcVertices == "hltPhase2PixelVertices":
        suffix += "WithOnlineVertices"
      elif srcVertices == "hltPhase2TrimmedPixelVertices":
        suffix += "WithOnlineVerticesTrimmed"
      else:
        raise ValueError("Invalid parameter srcVertices = '%s' !!" % srcVertices)  

      moduleName_wrtGenHadTaus = "analyze%sResponse%sWrtGenHadTaus" % (pfTauLabel, suffix)
      module_wrtGenHadTaus = cms.EDAnalyzer("RecoPFTauResponseAnalyzer",
        srcPFTaus = cms.InputTag('hltSelected%ss%s' % (pfTauLabel, suffix)),
        srcPFTauSumChargedIso = cms.InputTag('hlt%sChargedIsoPtSum%s' % (pfTauLabel, suffix)),
        srcRefTaus = cms.InputTag('offlineMatchedGenHadTaus'),
        typeRefTaus = cms.string("gen"),                                                                            
        dqmDirectory = cms.string("%s/%sResponseAnalyzer%s_wrtGenHadTaus" % (srcVertices, pfTauLabel, suffix))
      )
      setattr(process, moduleName_wrtGenHadTaus, module_wrtGenHadTaus)
      process.analysisSequence += module_wrtGenHadTaus

      moduleName_wrtOfflineTaus = "analyze%sResponse%sWrtOfflineTaus" % (pfTauLabel, suffix)
      module_wrtOfflineTaus = cms.EDAnalyzer("RecoPFTauResponseAnalyzer",
        srcPFTaus = cms.InputTag('hltSelected%ss%s' % (pfTauLabel, suffix)),
        srcPFTauSumChargedIso = cms.InputTag('hlt%sChargedIsoPtSum%s' % (pfTauLabel, suffix)),
        srcRefTaus = cms.InputTag('selectedOfflinePFTaus'),
        typeRefTaus = cms.string("offline"),                                                                   
        dqmDirectory = cms.string("%s/%sResponseAnalyzer%s_wrtOfflineTaus" % (srcVertices, pfTauLabel, suffix))
      )
      setattr(process, moduleName_wrtOfflineTaus, module_wrtOfflineTaus)
      process.analysisSequence += module_wrtOfflineTaus
#--------------------------------------------------------------------------------

process.load("DQMServices.Core.DQMStore_cfi")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string(outputFileName)
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
