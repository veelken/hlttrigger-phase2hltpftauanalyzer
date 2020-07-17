
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
        'file:/home/veelken/Phase2HLT/CMSSW_11_1_0/src/HLTrigger/Phase2HLTPFTaus/test/step3_RAW2DIGI_RECO.root'
    )
)

inputFilePath = '/hdfs/cms/store/user/rdewanje/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5/HLTConfig_w_offlineVtxCollection_HGCalFix_VBFHTT_Phase2HLTTDRWinter20_PU200_CMSSW_11_1_0_pre8/200627_142633/'
processName = "qqH_htt"
hlt_srcVertices = 'offlinePrimaryVertices'
#hlt_srcVertices = 'hltPhase2PixelVertices'
#hlt_algorithms = [ "hps", "shrinking-cone" ]
hlt_algorithms = [ "hps" ]
hlt_isolation_maxDeltaZOptions = [ "primaryVertex", "leadTrack" ]
hlt_isolation_minTrackHits = 8
outputFileName = "analyzePFTaus_signal_%s_%s_DEBUG.root" % (processName, hlt_srcVertices)

##inputFilePath = None
##inputFileNames = $inputFileNames
##processName = "$processName"
##hlt_srcVertices = '$hlt_srcVertices'
##hlt_algorithms = [ "$hlt_algorithm" ]
##hlt_isolation_maxDeltaZOptions = [ "$hlt_isolation_maxDeltaZOption" ]
##hlt_isolation_minTrackHits = $hlt_isolation_minTrackHits
##outputFileName = "$outputFileName"

#--------------------------------------------------------------------------------
# set input files
if inputFilePath:
    from HLTrigger.TallinnHLTPFTauAnalyzer.tools.jobTools import getInputFileNames
    print("Searching for input files in path = '%s'" % inputFilePath)
    inputFileNames = getInputFileNames(inputFilePath)
    print("Found %i input files." % len(inputFileNames))
    #process.source.fileNames = cms.untracked.vstring(inputFileNames)
else:
    print("Processing %i input files: %s" % (len(inputFileNames), inputFileNames))
    process.source.fileNames = cms.untracked.vstring(inputFileNames)
#--------------------------------------------------------------------------------

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.analysisSequence = cms.Sequence()

#--------------------------------------------------------------------------------
# CV: match offline recontructed and generator-level taus
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
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# CV: fill plots related to tracking and vertexing
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
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# CV: fill HLT tau plots
for hlt_algorithm in hlt_algorithms:

  hlt_pfTauLabel = None
  if hlt_algorithm == "shrinking-cone":
    hlt_pfTauLabel = "PFTau"
  elif hlt_algorithm == "hps":
    hlt_pfTauLabel = "HpsPFTau"
  else:
    raise ValueError("Invalid parameter hlt_algorithm = '%s' !!" % hlt_algorithm)

  for hlt_isolation_maxDeltaZOption in hlt_isolation_maxDeltaZOptions:

    suffix = "%iHits" % hlt_isolation_minTrackHits
    if hlt_isolation_maxDeltaZOption == "primaryVertex":
      suffix += "MaxDeltaZ"
    elif hlt_isolation_maxDeltaZOption == "leadTrack":
      suffix += "MaxDeltaZToLeadTrack"
    else:
      raise ValueError("Invalid parameter hlt_isolation_maxDeltaZOption = '%s' !!" % hlt_isolation_maxDeltaZOption)
    if hlt_srcVertices == "offlinePrimaryVertices":
      suffix += "WithOfflineVertices"
    elif hlt_srcVertices == "hltPhase2PixelVertices":
      suffix += "WithOnlineVertices"
    elif hlt_srcVertices == "hltPhase2TrimmedPixelVertices":
      suffix += "WithOnlineVerticesTrimmed"
    else:
      raise ValueError("Invalid parameter hlt_srcVertices = '%s' !!" % hlt_srcVertices)

    moduleName_PFTauAnalyzerSignal_wrtGenHadTaus = "analyze%ss%sWrtGenHadTaus" % (hlt_pfTauLabel, suffix)
    module_PFTauAnalyzerSignal_wrtGenHadTaus = cms.EDAnalyzer("RecoPFTauAnalyzerSignal",
      srcPFTaus = cms.InputTag('hltSelected%ss%s' % (hlt_pfTauLabel, suffix)),
      srcPFTauSumChargedIso = cms.InputTag('hlt%sChargedIsoPtSum%s' % (hlt_pfTauLabel, suffix)),
      srcDenominator = cms.InputTag('offlineMatchedGenHadTaus'),
      typeDenominator = cms.string("gen"),                                                                            
      lumiScale = cms.double(1.),
      dqmDirectory = cms.string("%s/%sAnalyzerSignal%s_wrtGenHadTaus" % (hlt_srcVertices, hlt_pfTauLabel, suffix))
    )
    setattr(process, moduleName_PFTauAnalyzerSignal_wrtGenHadTaus, module_PFTauAnalyzerSignal_wrtGenHadTaus)
    process.analysisSequence += module_PFTauAnalyzerSignal_wrtGenHadTaus

    from HLTrigger.Phase2HLTPFTaus.PFTauPairProducer_cfi import PFTauPairs
    moduleName_PFTauPairProducer = "hlt%sPairs%s" % (hlt_pfTauLabel, suffix)
    module_PFTauPairProducer = PFTauPairs.clone(
      srcPFTaus = cms.InputTag('hltSelected%ss%s' % (hlt_pfTauLabel, suffix)),
      srcPFTauSumChargedIso = cms.InputTag('hlt%sChargedIsoPtSum%s' % (hlt_pfTauLabel, suffix))
    )
    setattr(process, moduleName_PFTauPairProducer, module_PFTauPairProducer)
    process.analysisSequence += module_PFTauPairProducer

    moduleName_PFTauPairAnalyzer_wrtGenHadTaus = "analyze%sPairs%sWrtGenHadTaus" % (hlt_pfTauLabel, suffix)
    module_PFTauPairAnalyzer_wrtGenHadTaus = cms.EDAnalyzer("RecoPFTauPairAnalyzer",
      srcPFTauPairs = cms.InputTag(moduleName_PFTauPairProducer),
      srcRefTaus = cms.InputTag('offlineMatchedGenHadTaus'),
      min_refTau_pt = cms.double(20.),
      max_refTau_pt = cms.double(1.e+3),                                                                
      min_refTau_absEta = cms.double(-1.),
      max_refTau_absEta = cms.double(2.4),     
      lumiScale = cms.double(1.),                                               
      dqmDirectory = cms.string("%s/%sPairAnalyzer%s_wrtGenHadTaus" % (hlt_srcVertices, hlt_pfTauLabel, suffix))
    )
    setattr(process, moduleName_PFTauPairAnalyzer_wrtGenHadTaus, module_PFTauPairAnalyzer_wrtGenHadTaus)
    process.analysisSequence += module_PFTauPairAnalyzer_wrtGenHadTaus

    moduleName_PFTauAnalyzerSignal_wrtOfflineTaus = "analyze%ss%sWrtOfflineTaus" % (hlt_pfTauLabel, suffix)
    modulePF_PFTauAnalyzerSignal_wrtOfflineTaus = cms.EDAnalyzer("RecoPFTauAnalyzerSignal",
    srcPFTaus = cms.InputTag('hltSelected%ss%s' % (hlt_pfTauLabel, suffix)),
      srcPFTauSumChargedIso = cms.InputTag('hlt%sChargedIsoPtSum%s' % (hlt_pfTauLabel, suffix)),
      srcDenominator = cms.InputTag('selectedOfflinePFTaus'),
      typeDenominator = cms.string("offline"),  
      lumiScale = cms.double(1.),   
      dqmDirectory = cms.string("%s/%sAnalyzerSignal%s_wrtOfflineTaus" % (hlt_srcVertices, hlt_pfTauLabel, suffix))
    )
    setattr(process, moduleName_PFTauAnalyzerSignal_wrtOfflineTaus, modulePF_PFTauAnalyzerSignal_wrtOfflineTaus)
    process.analysisSequence += modulePF_PFTauAnalyzerSignal_wrtOfflineTaus

    moduleName_PFTauPairAnalyzer_wrtOfflineTaus = "analyze%sPairs%sWrtOfflineTaus" % (hlt_pfTauLabel, suffix)
    module_PFTauPairAnalyzer_wrtOfflineTaus = cms.EDAnalyzer("RecoPFTauPairAnalyzer",
      srcPFTauPairs = cms.InputTag(moduleName_PFTauPairProducer),
      srcRefTaus = cms.InputTag('selectedOfflinePFTaus'),
      min_refTau_pt = cms.double(20.),
      max_refTau_pt = cms.double(1.e+3),                                                                
      min_refTau_absEta = cms.double(-1.),
      max_refTau_absEta = cms.double(2.4),  
      lumiScale = cms.double(1.),                                                              
      dqmDirectory = cms.string("%s/%sPairAnalyzer%s_wrtOfflineTaus" % (hlt_srcVertices, hlt_pfTauLabel, suffix))
    )
    setattr(process, moduleName_PFTauPairAnalyzer_wrtOfflineTaus, module_PFTauPairAnalyzer_wrtOfflineTaus)
    process.analysisSequence += module_PFTauPairAnalyzer_wrtOfflineTaus

    moduleName_PFTauIsolationAnalyzer = "analyze%sIsolation%s" % (hlt_pfTauLabel, suffix)
    module_PFTauIsolationAnalyzer = cms.EDAnalyzer("RecoPFTauIsolationAnalyzer",
      srcPFTaus = cms.InputTag('hltSelected%ss%s' % (hlt_pfTauLabel, suffix)),
      srcPFTauSumChargedIso = cms.InputTag('hlt%sChargedIsoPtSum%s' % (hlt_pfTauLabel, suffix)),
      srcPFTauSumNeutralIso = cms.InputTag('hlt%sNeutralIsoPtSum%s' % (hlt_pfTauLabel, suffix)),
      srcGenTaus = cms.InputTag(''),
      dRmatch = cms.double(0.3),                                                            
      srcRho = cms.InputTag('hltKT6PFJets:rho'),
      #inputFileName_rhoCorr = cms.string("HLTrigger/TallinnHLTPFTauAnalyzer/data/rhoCorr.root"),
      #histogramName_rhoCorr = cms.string("DQMData/RhoCorrAnalyzer/neutralPFCandPt_vs_absEta"),
      inputFileName_rhoCorr = cms.string(""),
      histogramName_rhoCorr = cms.string(""),     
      lumiScale = cms.double(1.),                                              
      dqmDirectory = cms.string("%s/%sIsolationAnalyzer%s" % (hlt_srcVertices, hlt_pfTauLabel, suffix))
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
