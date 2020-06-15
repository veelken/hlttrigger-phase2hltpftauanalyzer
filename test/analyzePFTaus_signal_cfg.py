
import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzeTallinnL1PFTausSignal")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50000)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/home/veelken/Phase2HLT/CMSSW_11_1_0_pre6/src/HLTTrigger/Phase2HLTPFTaus/test/step3_RAW2DIGI_RECO.root'
    )
)

#--------------------------------------------------------------------------------
# set input files

import os
import re

inputFilePaths = [
    '/hdfs/cms/store/user/rdewanje/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5/HLTConfig_Christian_VBFHTT_Phase2HLTTDRWinter20_PU200_CMSSW_11_1_0_pre6_numCores8_maxMem16kMB_T2_EE_Estonia_blacklist/200612_212847/0000/'
]

inputFile_regex = r"[a-zA-Z0-9_/:.-]*step3_RAW2DIGI_RECO_[a-zA-Z0-9-_]+.root"

# check if name of inputFile matches regular expression
inputFileNames = []
for inputFilePath in inputFilePaths:
    files = [ "".join([ "file:", inputFilePath, file ]) for file in os.listdir(inputFilePath) ]
    for file in files:
        inputFile_matcher = re.compile(inputFile_regex)
        if inputFile_matcher.match(file):
            inputFileNames.append(file)
print "inputFileNames = %s" % inputFileNames 

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

process.analyzeVertices = cms.EDAnalyzer("RecoVertexAnalyzer",
  srcGenVertex_z = cms.InputTag('genVertex:z0'),
  #srcHLTVertices = cms.InputTag('hltPixelVertices'),
  srcHLTVertices = cms.InputTag('offlinePrimaryVertices'),                        
  srcOfflineVertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
  dqmDirectory = cms.string("recoVertexAnalyzer")
)
process.analysisSequence += process.analyzeVertices

# CV: reco::Track and reco::PFCandidate collections reconstructed offline and at HLT level not yet separately stored in ROOT file,
#     which is why the same collections are used for "Offline" and "HLT" inputs

process.analyzeTracksWrtRecVertex = cms.EDAnalyzer("RecoTrackAnalyzer",
  srcGenTaus = cms.InputTag('selectedGenHadTaus'),
  vtxMode = cms.string("recVtx"),
  srcOfflineVertices = cms.InputTag('offlineSlimmedPrimaryVertices'),                                       
  srcOfflineTracks = cms.InputTag('generalTracks'),                                                     
  srcOfflinePFCands = cms.InputTag('packedPFCandidates'),
  #srcHLTVertices = cms.InputTag('hltPixelVertices'),
  srcHLTVertices = cms.InputTag('offlinePrimaryVertices'),                                                  
  srcHLTTracks = cms.InputTag('generalTracks'),
  srcHLTPFCands = cms.InputTag('particleFlowTmp'),
  dqmDirectory = cms.string("recoTrackAnalyzerWrtRecVertex"),
  debug = cms.bool(False)                                     
)
##process.analysisSequence += process.analyzeTracksWrtRecVertex

process.genVertex = cms.EDProducer("GenVertexProducer",
  src = cms.InputTag('prunedGenParticles'),
  pdgIds = cms.vint32(-15, +15) # CV: use { -15, +15 } for signal, empty list for background                                    
) 
process.analysisSequence += process.genVertex

process.analyzeTracksWrtGenVertex = cms.EDAnalyzer("RecoTrackAnalyzer",
  srcGenTaus = cms.InputTag('selectedGenHadTaus'),
  vtxMode = cms.string("genVtx"),
  srcGenVertex_position = cms.InputTag('genVertex:position'),                                                   
  srcOfflineTracks = cms.InputTag('generalTracks'),
  srcOfflinePFCands = cms.InputTag('packedPFCandidates'),
  srcHLTTracks = cms.InputTag('generalTracks'),
  srcHLTPFCands = cms.InputTag('particleFlowTmp'),
  dqmDirectory = cms.string("recoTrackAnalyzerWrtGenVertex"),
  debug = cms.bool(False)                                     
)
##process.analysisSequence += process.analyzeTracksWrtGenVertex

for algorithm in [ "hps", "shrinking-cone" ]:
  pfTauLabel = None
  if algorithm == "shrinking-cone":
    pfTauLabel = "PFTau"
  elif algorithm == "hps":
    pfTauLabel = "HpsPFTau"
  else:
    raise ValueError("Invalid parameter algorithm = '%s' !!" % algorithm)

  for srcVertices in [ "offlinePrimaryVertices", "hltPhase2PixelVertices" ]:
    suffix = None
    if srcVertices == "offlinePrimaryVertices":
      suffix = "WithOfflineVertices"
    elif srcVertices == "hltPhase2PixelVertices":
      suffix = "WithOnlineVertices"
    else:
      raise ValueError("Invalid parameter srcVertices = '%s' !!" % srcVertices)

    moduleName_PFTauAnalyzerSignal_wrtGenHadTaus = "analyze%ss%sWrtGenHadTaus" % (pfTauLabel, suffix)
    modulePF_PFTauAnalyzerSignal_wrtGenHadTaus = cms.EDAnalyzer("RecoPFTauAnalyzerSignal",
      srcPFTaus = cms.InputTag('hltSelected%ss%s' % (pfTauLabel, suffix)),
      srcPFTauSumChargedIso = cms.InputTag('hlt%sChargedIsoPtSum%s' % (pfTauLabel, suffix)),
      srcDenominator = cms.InputTag('offlineMatchedGenHadTaus'),
      typeDenominator = cms.string("gen"),                                                                            
      dqmDirectory = cms.string("%sAnalyzerSignal%s_wrtGenHadTaus" % (pfTauLabel, suffix))
    )
    setattr(process, moduleName_PFTauAnalyzerSignal_wrtGenHadTaus, modulePF_PFTauAnalyzerSignal_wrtGenHadTaus)
    process.analysisSequence += modulePF_PFTauAnalyzerSignal_wrtGenHadTaus

    moduleName_PFTauPairAnalyzer_wrtGenHadTaus = "analyze%sPairs%sWrtGenHadTaus" % (pfTauLabel, suffix)
    module_PFTauPairAnalyzer_wrtGenHadTaus = cms.EDAnalyzer("RecoPFTauPairAnalyzer",
      srcPFTaus = cms.InputTag('hltSelected%ss%s' % (pfTauLabel, suffix)),
      srcPFTauSumChargedIso = cms.InputTag('hlt%sChargedIsoPtSum%s' % (pfTauLabel, suffix)),
      srcRefTaus = cms.InputTag('offlineMatchedGenHadTaus'),
      min_refTau_pt = cms.double(20.),
      max_refTau_pt = cms.double(1.e+3),                                                                
      min_refTau_absEta = cms.double(-1.),
      max_refTau_absEta = cms.double(2.4),                                                                
      dqmDirectory = cms.string("%sPairAnalyzer%s_wrtGenHadTaus" % (pfTauLabel, suffix))
    )
    setattr(process, moduleName_PFTauPairAnalyzer_wrtGenHadTaus, module_PFTauPairAnalyzer_wrtGenHadTaus)
    process.analysisSequence += module_PFTauPairAnalyzer_wrtGenHadTaus

    moduleName_PFTauAnalyzerSignal_wrtOfflineTaus = "analyze%ss%sWrtOfflineTaus" % (pfTauLabel, suffix)
    modulePF_PFTauAnalyzerSignal_wrtOfflineTaus = cms.EDAnalyzer("RecoPFTauAnalyzerSignal",
      srcPFTaus = cms.InputTag('hltSelected%ss%s' % (pfTauLabel, suffix)),
      srcPFTauSumChargedIso = cms.InputTag('hlt%sChargedIsoPtSum%s' % (pfTauLabel, suffix)),
      srcDenominator = cms.InputTag('selectedOfflinePFTaus'),
      typeDenominator = cms.string("offline"),  
      dqmDirectory = cms.string("%sAnalyzerSignal%s_wrtOfflineTaus" % (pfTauLabel, suffix))
    )
    setattr(process, moduleName_PFTauAnalyzerSignal_wrtOfflineTaus, modulePF_PFTauAnalyzerSignal_wrtOfflineTaus)
    process.analysisSequence += modulePF_PFTauAnalyzerSignal_wrtOfflineTaus

    moduleName_PFTauPairAnalyzer_wrtOfflineTaus = "analyze%sPairs%sWrtOfflineTaus" % (pfTauLabel, suffix)
    module_PFTauPairAnalyzer_wrtOfflineTaus = cms.EDAnalyzer("RecoPFTauPairAnalyzer",
      srcPFTaus = cms.InputTag('hltSelected%ss%s' % (pfTauLabel, suffix)),
      srcPFTauSumChargedIso = cms.InputTag('hlt%sChargedIsoPtSum%s' % (pfTauLabel, suffix)),
      srcRefTaus = cms.InputTag('selectedOfflinePFTaus'),
      min_refTau_pt = cms.double(20.),
      max_refTau_pt = cms.double(1.e+3),                                                                
      min_refTau_absEta = cms.double(-1.),
      max_refTau_absEta = cms.double(2.4),                                                                
      dqmDirectory = cms.string("%sPairAnalyzer%s_wrtOfflineTaus" % (pfTauLabel, suffix))
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
      dqmDirectory = cms.string("%sIsolationAnalyzer%s" % (pfTauLabel, suffix))
    )
    setattr(process, moduleName_PFTauIsolationAnalyzer, module_PFTauIsolationAnalyzer)
    process.analysisSequence += module_PFTauIsolationAnalyzer
#--------------------------------------------------------------------------------

process.load("DQMServices.Core.DQMStore_cfi")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('analyzePFTaus_signal_2020Jun03.root')
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
