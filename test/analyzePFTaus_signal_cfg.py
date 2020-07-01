
import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzePFTausSignal")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
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
outputFileName = "analyzePFTaus_signal_%s_%s_DEBUG.root" % (processName, srcVertices)

##inputFilePath = None
##inputFileNames = $inputFileNames
##processName = "$processName"
##srcVertices = '$srcVertices'
##outputFileName = "$outputFileName"

#--------------------------------------------------------------------------------
# set input files
if inputFilePath:
    from HLTTrigger.TallinnHLTPFTauAnalyzer.tools import getInputFileNames
    print("Searching for input files in path = '%s'" % inputFilePath)
    inputFileNames = getInputFileNames.getInputFileNames(inputFilePath)
    print("Found %i input files." % len(inputFileNames))
else:
    print("Processing %i input files: %s" % (len(inputFileNames), inputFileNames))
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

process.analyzeVertices1 = cms.EDAnalyzer("RecoVertexAnalyzer",
  srcGenVertex_z = cms.InputTag('genVertex:z0'),
  srcRecVertices = cms.InputTag('offlinePrimaryVertices'), 
  dqmDirectory = cms.string("recoVertexAnalyzer/offlinePrimaryVertices")
)
process.analysisSequence += process.analyzeVertices1

process.analyzeVertices2 = process.analyzeVertices1.clone(
  srcRecVertices = cms.InputTag('hltPhase2PixelVertices'), 
  dqmDirectory = cms.string("recoVertexAnalyzer/hltPhase2PixelVertices")
)
process.analysisSequence += process.analyzeVertices2

process.analyzeVertices3 = process.analyzeVertices1.clone(
  srcRecVertices = cms.InputTag('hltPhase2TrimmedPixelVertices'), 
  dqmDirectory = cms.string("recoVertexAnalyzer/hltPhase2TrimmedPixelVertices")
)
process.analysisSequence += process.analyzeVertices3

# CV: select subset of reco::Track and reco::PFCandidate objects within dR < 0.8 cones around generator-level hadronic tau decays
#     in order to reduce computing time
process.generalTracksGenHadTauMatched = cms.EDProducer("RecoTrackAntiOverlapSelector",
  src = cms.InputTag('generalTracks'),
  srcNotToBeFiltered = cms.VInputTag('selectedGenHadTaus'),
  dRmin = cms.double(0.8),
  invert = cms.bool(True)
)
process.analysisSequence += process.generalTracksGenHadTauMatched

process.packedPFCandidatesGenHadTauMatched = cms.EDFilter("PackedCandidateAntiOverlapSelector",
  src = cms.InputTag('packedPFCandidates'),
  srcNotToBeFiltered = cms.VInputTag('selectedGenHadTaus'),
  dRmin = cms.double(0.8),
  invert = cms.bool(True),
  filter = cms.bool(False)
)
process.analysisSequence += process.packedPFCandidatesGenHadTauMatched

process.particleFlowTmpGenHadTauMatched = cms.EDFilter("PFCandidateAntiOverlapSelector",
  src = cms.InputTag('particleFlowTmp'),
  srcNotToBeFiltered = cms.VInputTag('selectedGenHadTaus'),
  dRmin = cms.double(0.8),
  invert = cms.bool(True),
  filter = cms.bool(False)
)
process.analysisSequence += process.particleFlowTmpGenHadTauMatched

# CV: reco::Track collections reconstructed offline and at HLT level not yet separately stored in ROOT file,
#     which is why the same collections are used for "Offline" and "HLT" inputs
process.analyzeTracksWrtRecVertex = cms.EDAnalyzer("RecoTrackAnalyzer",
  srcGenTaus = cms.InputTag('selectedGenHadTaus'),
  vtxMode = cms.string("recVtx"),
  srcOfflineVertices = cms.InputTag('offlineSlimmedPrimaryVertices'),                                       
  srcOfflineTracks = cms.InputTag('generalTracksGenHadTauMatched'),                                                     
  srcOfflinePFCands = cms.InputTag('packedPFCandidatesGenHadTauMatched'),
  srcOnlineVertices = cms.InputTag('hltPhase2PixelVertices'),
  #srcOnlineVertices = cms.InputTag('offlinePrimaryVertices'),                                                  
  srcOnlineTracks = cms.InputTag('generalTracksGenHadTauMatched'),
  srcOnlinePFCands = cms.InputTag('particleFlowTmpGenHadTauMatched'),
  dqmDirectory = cms.string("recoTrackAnalyzerWrtRecVertex"),
  debug = cms.bool(False)                                     
)
process.analysisSequence += process.analyzeTracksWrtRecVertex

process.genVertex = cms.EDProducer("GenVertexProducer",
  src = cms.InputTag('prunedGenParticles'),
  pdgIds = cms.vint32(-15, +15) # CV: use { -15, +15 } for signal, empty list for background                                    
) 
process.analysisSequence += process.genVertex

process.analyzeTracksWrtGenVertex = cms.EDAnalyzer("RecoTrackAnalyzer",
  srcGenTaus = cms.InputTag('selectedGenHadTaus'),
  vtxMode = cms.string("genVtx"),
  srcGenVertex_position = cms.InputTag('genVertex:position'),                                                   
  srcOfflineTracks = cms.InputTag('generalTracksGenHadTauMatched'),
  srcOfflinePFCands = cms.InputTag('packedPFCandidatesGenHadTauMatched'),
  srcOnlineVertices = cms.InputTag('hltPhase2PixelVertices'),
  #srcOnlineVertices = cms.InputTag('offlinePrimaryVertices'),  
  srcOnlineTracks = cms.InputTag('generalTracksGenHadTauMatched'),
  srcOnlinePFCands = cms.InputTag('particleFlowTmpGenHadTauMatched'),
  dqmDirectory = cms.string("recoTrackAnalyzerWrtGenVertex"),
  debug = cms.bool(False)                                     
)
process.analysisSequence += process.analyzeTracksWrtGenVertex

for algorithm in [ "hps", "shrinking-cone" ]:

  pfTauLabel = None
  if algorithm == "shrinking-cone":
    pfTauLabel = "PFTau"
  elif algorithm == "hps":
    pfTauLabel = "HpsPFTau"
  else:
    raise ValueError("Invalid parameter algorithm = '%s' !!" % algorithm)

  for isolation_maxDeltaZOption in [ "primaryVertex", "leadTrack" ]:
    for isolation_minTrackHits in [ 3, 5, 8 ]:  

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

      moduleName_PFTauAnalyzerSignal_wrtGenHadTaus = "analyze%ss%sWrtGenHadTaus" % (pfTauLabel, suffix)
      modulePF_PFTauAnalyzerSignal_wrtGenHadTaus = cms.EDAnalyzer("RecoPFTauAnalyzerSignal",
        srcPFTaus = cms.InputTag('hltSelected%ss%s' % (pfTauLabel, suffix)),
        srcPFTauSumChargedIso = cms.InputTag('hlt%sChargedIsoPtSum%s' % (pfTauLabel, suffix)),
        srcDenominator = cms.InputTag('offlineMatchedGenHadTaus'),
        typeDenominator = cms.string("gen"),                                                                            
        lumiScale = cms.double(1.),
        dqmDirectory = cms.string("%s/%sAnalyzerSignal%s_wrtGenHadTaus" % (srcVertices, pfTauLabel, suffix))
      )
      setattr(process, moduleName_PFTauAnalyzerSignal_wrtGenHadTaus, modulePF_PFTauAnalyzerSignal_wrtGenHadTaus)
      process.analysisSequence += modulePF_PFTauAnalyzerSignal_wrtGenHadTaus

      from HLTTrigger.Phase2HLTPFTaus.PFTauPairProducer_cfi import PFTauPairs
      moduleName_PFTauPairProducer = "hlt%sPairs%s" % (pfTauLabel, suffix)
      module_PFTauPairProducer = PFTauPairs.clone(
        srcPFTaus = cms.InputTag('hltSelected%ss%s' % (pfTauLabel, suffix)),
        srcPFTauSumChargedIso = cms.InputTag('hlt%sChargedIsoPtSum%s' % (pfTauLabel, suffix))
      )
      setattr(process, moduleName_PFTauPairProducer, module_PFTauPairProducer)
      process.analysisSequence += module_PFTauPairProducer

      moduleName_PFTauPairAnalyzer_wrtGenHadTaus = "analyze%sPairs%sWrtGenHadTaus" % (pfTauLabel, suffix)
      module_PFTauPairAnalyzer_wrtGenHadTaus = cms.EDAnalyzer("RecoPFTauPairAnalyzer",
        srcPFTauPairs = cms.InputTag(moduleName_PFTauPairProducer),
        srcRefTaus = cms.InputTag('offlineMatchedGenHadTaus'),
        min_refTau_pt = cms.double(20.),
        max_refTau_pt = cms.double(1.e+3),                                                                
        min_refTau_absEta = cms.double(-1.),
        max_refTau_absEta = cms.double(2.4),     
        lumiScale = cms.double(1.),                                               
        dqmDirectory = cms.string("%s/%sPairAnalyzer%s_wrtGenHadTaus" % (srcVertices, pfTauLabel, suffix))
      )
      setattr(process, moduleName_PFTauPairAnalyzer_wrtGenHadTaus, module_PFTauPairAnalyzer_wrtGenHadTaus)
      process.analysisSequence += module_PFTauPairAnalyzer_wrtGenHadTaus

      moduleName_PFTauAnalyzerSignal_wrtOfflineTaus = "analyze%ss%sWrtOfflineTaus" % (pfTauLabel, suffix)
      modulePF_PFTauAnalyzerSignal_wrtOfflineTaus = cms.EDAnalyzer("RecoPFTauAnalyzerSignal",
        srcPFTaus = cms.InputTag('hltSelected%ss%s' % (pfTauLabel, suffix)),
        srcPFTauSumChargedIso = cms.InputTag('hlt%sChargedIsoPtSum%s' % (pfTauLabel, suffix)),
        srcDenominator = cms.InputTag('selectedOfflinePFTaus'),
        typeDenominator = cms.string("offline"),  
        lumiScale = cms.double(1.),   
        dqmDirectory = cms.string("%s/%sAnalyzerSignal%s_wrtOfflineTaus" % (srcVertices, pfTauLabel, suffix))
      )
      setattr(process, moduleName_PFTauAnalyzerSignal_wrtOfflineTaus, modulePF_PFTauAnalyzerSignal_wrtOfflineTaus)
      process.analysisSequence += modulePF_PFTauAnalyzerSignal_wrtOfflineTaus

      moduleName_PFTauPairAnalyzer_wrtOfflineTaus = "analyze%sPairs%sWrtOfflineTaus" % (pfTauLabel, suffix)
      module_PFTauPairAnalyzer_wrtOfflineTaus = cms.EDAnalyzer("RecoPFTauPairAnalyzer",
        srcPFTauPairs = cms.InputTag(moduleName_PFTauPairProducer),
        srcRefTaus = cms.InputTag('selectedOfflinePFTaus'),
        min_refTau_pt = cms.double(20.),
        max_refTau_pt = cms.double(1.e+3),                                                                
        min_refTau_absEta = cms.double(-1.),
        max_refTau_absEta = cms.double(2.4),  
        lumiScale = cms.double(1.),                                                              
        dqmDirectory = cms.string("%s/%sPairAnalyzer%s_wrtOfflineTaus" % (srcVertices, pfTauLabel, suffix))
      )
      setattr(process, moduleName_PFTauPairAnalyzer_wrtOfflineTaus, module_PFTauPairAnalyzer_wrtOfflineTaus)
      process.analysisSequence += module_PFTauPairAnalyzer_wrtOfflineTaus

      moduleName_PFTauIsolationAnalyzer = "analyze%sIsolation%s" % (pfTauLabel, suffix)
      module_PFTauIsolationAnalyzer = cms.EDAnalyzer("RecoPFTauIsolationAnalyzer",
        srcPFTaus = cms.InputTag('hltSelected%ss%s' % (pfTauLabel, suffix)),
        srcPFTauSumChargedIso = cms.InputTag('hlt%sChargedIsoPtSum%s' % (pfTauLabel, suffix)),
        srcPFTauSumNeutralIso = cms.InputTag('hlt%sNeutralIsoPtSum%s' % (pfTauLabel, suffix)),
        srcGenTaus = cms.InputTag(''),
        dRmatch = cms.double(0.3),                                                            
        srcRho = cms.InputTag('hltKT6PFJets:rho'),
        #inputFileName_rhoCorr = cms.string("HLTTrigger/TallinnHLTPFTauAnalyzer/data/rhoCorr.root"),
        #histogramName_rhoCorr = cms.string("DQMData/RhoCorrAnalyzer/neutralPFCandPt_vs_absEta"),
        inputFileName_rhoCorr = cms.string(""),
        histogramName_rhoCorr = cms.string(""),     
        lumiScale = cms.double(1.),                                              
        dqmDirectory = cms.string("%s/%sIsolationAnalyzer%s" % (srcVertices, pfTauLabel, suffix))
      )
      setattr(process, moduleName_PFTauIsolationAnalyzer, module_PFTauIsolationAnalyzer)
      ##process.analysisSequence += module_PFTauIsolationAnalyzer
#--------------------------------------------------------------------------------

process.load("DQMServices.Core.DQMStore_cfi")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string(outputFileName)
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
