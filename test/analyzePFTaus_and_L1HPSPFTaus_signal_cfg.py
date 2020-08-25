
import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzePFTausSignal")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/home/veelken/Phase2HLT/CMSSW_11_1_0/src/HLTrigger/Phase2HLTPFTaus/test/step3_RAW2DIGI_RECO.root'
    )
)

inputFilePath = '/hdfs/cms/store/user/rdewanje/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5/HLTConfig_VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5_wOfflineVtx_wDeepTau3/'
processName = "qqH_htt"
hlt_srcVertices = 'offlinePrimaryVertices'
#hlt_srcVertices = 'hltPhase2PixelVertices'
#hlt_algorithms = [ "hps", "shrinking-cone" ]
hlt_algorithms = [ "hps" ]
hlt_isolation_maxDeltaZOptions = [ "primaryVertex", "leadTrack" ]
hlt_isolation_minTrackHits = 8
l1_useStrips = True
outputFileName = "analyzePFTaus_signal_%s_%s_DEBUG.root" % (processName, hlt_srcVertices)

##inputFilePath = None
##inputFileNames = $inputFileNames
##processName = "$processName"
##hlt_srcVertices = '$hlt_srcVertices'
##hlt_algorithms = [ "$hlt_algorithm" ]
##hlt_isolation_maxDeltaZOptions = [ "$hlt_isolation_maxDeltaZOption" ]
##hlt_isolation_minTrackHits = $hlt_isolation_minTrackHits
##l1_useStrips = '$l1_useStrips'
##outputFileName = "$outputFileName"

#hasMiniAOD = False
hasMiniAOD = True

#--------------------------------------------------------------------------------
# set input files
if inputFilePath:
    from HLTrigger.TallinnHLTPFTauAnalyzer.tools.jobTools import getInputFileNames
    print("Searching for input files in path = '%s'" % inputFilePath)
    inputFileNames = getInputFileNames(inputFilePath)
    print("Found %i input files." % len(inputFileNames))
    process.source.fileNames = cms.untracked.vstring(inputFileNames)
else:
    print("Processing %i input files: %s" % (len(inputFileNames), inputFileNames))
    process.source.fileNames = cms.untracked.vstring(inputFileNames)
#--------------------------------------------------------------------------------

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.analysisSequence = cms.Sequence()

#--------------------------------------------------------------------------------
# CV: compute event weights
process.lumiScale = cms.EDProducer("EvtWeightProducerLumiScale",
  lumiScale = cms.double(1.) # CV: no event weights needed for signal
)
process.analysisSequence += process.lumiScale
#--------------------------------------------------------------------------------

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

if hasMiniAOD:
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

if hasMiniAOD:
  srcGenHadTaus = 'offlineMatchedGenHadTaus'
else:
  srcGenHadTaus = 'selectedGenHadTaus'
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# CV: fill plots related to tracking and vertexing
process.genVertex = cms.EDProducer("GenVertexProducer",
  src = cms.InputTag('genParticles'),
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

if hasMiniAOD:
  # CV: select subset of reco::Track and reco::PFCandidate objects within dR < 0.8 cones around generator-level hadronic tau decays
  #     in order to reduce computing time
  process.generalTracksGenHadTauMatched = cms.EDProducer("RecoTrackAntiOverlapSelector",
    src = cms.InputTag('generalTracks'),
    srcNotToBeFiltered = cms.VInputTag(srcGenHadTaus),
    dRmin = cms.double(0.8),
    invert = cms.bool(True)
  )
  process.analysisSequence += process.generalTracksGenHadTauMatched

  process.packedPFCandidatesGenHadTauMatched = cms.EDFilter("PackedCandidateAntiOverlapSelector",
    src = cms.InputTag('packedPFCandidates'),
    srcNotToBeFiltered = cms.VInputTag(srcGenHadTaus),
    dRmin = cms.double(0.8),
    invert = cms.bool(True),
    filter = cms.bool(False)
  )
  process.analysisSequence += process.packedPFCandidatesGenHadTauMatched

  process.particleFlowTmpGenHadTauMatched = cms.EDFilter("PFCandidateAntiOverlapSelector",
    src = cms.InputTag('particleFlowTmp'),
    srcNotToBeFiltered = cms.VInputTag(srcGenHadTaus),
    dRmin = cms.double(0.8),
    invert = cms.bool(True),
    filter = cms.bool(False)
  )
  process.analysisSequence += process.particleFlowTmpGenHadTauMatched

  # CV: reco::Track collections reconstructed offline and at HLT level not yet separately stored in ROOT file,
  #     which is why the same collections are used for "Offline" and "HLT" inputs
  process.analyzeTracksWrtRecVertex = cms.EDAnalyzer("RecoTrackAnalyzer",
    srcGenTaus = cms.InputTag(srcGenHadTaus),
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
    src = cms.InputTag('genParticles'),
    pdgIds = cms.vint32(-15, +15) # CV: use { -15, +15 } for signal, empty list for background                                    
  ) 
  process.analysisSequence += process.genVertex

  process.analyzeTracksWrtGenVertex = cms.EDAnalyzer("RecoTrackAnalyzer",
    srcGenTaus = cms.InputTag(srcGenHadTaus),
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
# CV: fill L1 tau plots
moduleNameBase = "L1HPSPFTauProducer"
moduleLabel = None
if l1_useStrips:
    moduleLabel = "WithStrips"
else:
    moduleLabel = "WithoutStrips"
l1_srcPFTaus = moduleNameBase + moduleLabel + "PF"

moduleName_L1HPSPFTauAnalyzerSignalWrtGenHadTaus = "analyzeL1HPSPFTaus" + moduleLabel + "PFWrtGenHadTaus"
module_L1HPSPFTauAnalyzerSignalWrtGenHadTaus = cms.EDAnalyzer("L1HPSPFTauAnalyzerSignal",
  srcNumerator = cms.InputTag(l1_srcPFTaus),
  srcDenominator = cms.InputTag(srcGenHadTaus),
  typeDenominator = cms.string("gen"),                                                                            
  dqmDirectory = cms.string("L1HPSPFTauAnalyzerSignal" + moduleLabel + "PF_wrtGenHadTaus")
)
setattr(process, moduleName_L1HPSPFTauAnalyzerSignalWrtGenHadTaus, module_L1HPSPFTauAnalyzerSignalWrtGenHadTaus)
process.analysisSequence += getattr(process, moduleName_L1HPSPFTauAnalyzerSignalWrtGenHadTaus)

moduleName_L1HPSPFTauPairAnalyzerWrtGenHadTaus = "analyzeL1HPSPFTauPairs" + moduleLabel + "PFWrtGenHadTaus"
module_L1HPSPFTauPairAnalyzerWrtGenHadTaus = cms.EDAnalyzer("L1HPSPFTauPairAnalyzer",
  srcL1PFTaus = cms.InputTag(l1_srcPFTaus),
  srcRefTaus = cms.InputTag(srcGenHadTaus),
  min_refTau_pt = cms.double(20.),
  max_refTau_pt = cms.double(1.e+3),                                                                
  min_refTau_absEta = cms.double(-1.),
  max_refTau_absEta = cms.double(2.4),                                                                
  dqmDirectory = cms.string("L1HPSPFTauPairAnalyzer" + moduleLabel + "PF_wrtGenHadTaus")
)
setattr(process, moduleName_L1HPSPFTauPairAnalyzerWrtGenHadTaus, module_L1HPSPFTauPairAnalyzerWrtGenHadTaus)
process.analysisSequence += getattr(process, moduleName_L1HPSPFTauPairAnalyzerWrtGenHadTaus)

if hasMiniAOD:
  moduleName_L1HPSPFTauAnalyzerSignalWrtOfflineTaus = "analyzeL1HPSPFTaus" + moduleLabel + "PFWrtOfflineTaus"
  module_L1HPSPFTauAnalyzerSignalWrtOfflineTaus = cms.EDAnalyzer("L1HPSPFTauAnalyzerSignal",
    srcNumerator = cms.InputTag(l1_srcPFTaus),
    srcDenominator = cms.InputTag('selectedOfflinePFTaus'),
    typeDenominator = cms.string("offline"),                                                                   
    dqmDirectory = cms.string("L1HPSPFTauAnalyzerSignal" + moduleLabel + "PF_wrtOfflineTaus")
  )
  setattr(process, moduleName_L1HPSPFTauAnalyzerSignalWrtOfflineTaus, module_L1HPSPFTauAnalyzerSignalWrtOfflineTaus)
  process.analysisSequence += getattr(process, moduleName_L1HPSPFTauAnalyzerSignalWrtOfflineTaus)

  moduleName_L1HPSPFTauPairAnalyzerWrtOfflineTaus = "analyzeL1HPSPFTauPairs" + moduleLabel + "PFWrtOfflineTaus"
  module_L1HPSPFTauPairAnalyzerWrtOfflineTaus = cms.EDAnalyzer("L1HPSPFTauPairAnalyzer",
    srcL1PFTaus = cms.InputTag(l1_srcPFTaus),
    srcRefTaus = cms.InputTag('selectedOfflinePFTaus'),
    min_refTau_pt = cms.double(20.),
    max_refTau_pt = cms.double(1.e+3),                                                                
    min_refTau_absEta = cms.double(-1.),
    max_refTau_absEta = cms.double(2.4),
    dqmDirectory = cms.string("L1HPSPFTauPairAnalyzer" + moduleLabel + "PF_wrtOfflineTaus")
  )
  setattr(process, moduleName_L1HPSPFTauPairAnalyzerWrtOfflineTaus, module_L1HPSPFTauPairAnalyzerWrtOfflineTaus)
  process.analysisSequence += getattr(process, moduleName_L1HPSPFTauPairAnalyzerWrtOfflineTaus)

moduleName_L1HPSPFTauIsolationAnalyzer = "analyzeTallinL1PFTauIsolation" + moduleLabel + "PF"
module_L1HPSPFTauIsolationAnalyzer = cms.EDAnalyzer("L1HPSPFTauIsolationAnalyzer",
  srcL1PFTaus  = cms.InputTag(l1_srcPFTaus),
  srcGenTaus = cms.InputTag(''),
  dRmatch = cms.double(0.3),                                                            
  srcRho = cms.InputTag('kt6L1PFJetsPF:rho'),
  #inputFileName_rhoCorr = cms.string("L1Trigger/L1HPSPFTauAnalyzer/data/rhoCorr.root"),
  #histogramName_rhoCorr = cms.string("DQMData/RhoCorrAnalyzerPF/neutralPFCandPt_vs_absEta"),
  inputFileName_rhoCorr = cms.string(""),
  histogramName_rhoCorr = cms.string(""),                                                      
  dqmDirectory = cms.string("L1HPSPFTauIsolationAnalyzer" + moduleLabel + "PF")
)
setattr(process, moduleName_L1HPSPFTauIsolationAnalyzer, module_L1HPSPFTauIsolationAnalyzer)
process.analysisSequence += getattr(process, moduleName_L1HPSPFTauIsolationAnalyzer)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# CV: fill HLT tau plots (with and without matching to L1 taus)
for hlt_algorithm in hlt_algorithms:
  for hlt_matchToL1 in [ True, False ]:

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

      #----------------------------------------------------------------------------
      # CV: add DeepTau tau ID discriminator
      from HLTrigger.TallinnHLTPFTauAnalyzer.tools.addDeepTauDiscriminator import addDeepTauDiscriminator
      hlt_srcPFTaus = 'hltSelected%ss%s' % (hlt_pfTauLabel, suffix)
      hlt_srcPFJets = 'hlt%sAK4PFJets%s' % (hlt_pfTauLabel, suffix)
      deepTauSequenceName = "hltDeep%sSequence%s" % (hlt_pfTauLabel, suffix)
      deepTauSequence = addDeepTauDiscriminator(process, hlt_srcPFTaus, hlt_srcPFJets, hlt_srcVertices, hlt_pfTauLabel, suffix, deepTauSequenceName)
      process.analysisSequence += deepTauSequence
      #----------------------------------------------------------------------------

      hlt_srcPFTaus = 'hltSelected%ss%s' % (hlt_pfTauLabel, suffix)
      hlt_srcPFTauSumChargedIso = 'hltSelected%sChargedIsoPtSum%s' % (hlt_pfTauLabel, suffix)
      hlt_srcPFTauSumNeutralIso = 'hltSelected%sNeutralIsoPtSum%s' % (hlt_pfTauLabel, suffix)
      hlt_srcUpdatedPatTaus = 'hltUpdatedPat%ss%s' % (hlt_pfTauLabel, suffix)
      if hlt_matchToL1:
        suffix += "MatchedToL1"

        if not hasattr(process, "L1HPSPFTausPassingTrigger"):
          process.L1HPSPFTausPassingTrigger = cms.EDProducer("L1HPSPFTauSelector",
            src = cms.InputTag(l1_srcPFTaus),
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
            invert = cms.bool(False)
          )
          process.analysisSequence += process.L1HPSPFTausPassingTrigger

        moduleName_hltPFTausMatchedToL1 = "%sMatchedToL1" % hlt_srcPFTaus
        module_hltPFTausMatchedToL1 = cms.EDFilter("PFTauAntiOverlapSelector",
          src = cms.InputTag(hlt_srcPFTaus),
          srcNotToBeFiltered = cms.VInputTag('L1HPSPFTausPassingTrigger'),
          dRmin = cms.double(0.3),
          invert = cms.bool(True),
          filter = cms.bool(False)
        )
        setattr(process, moduleName_hltPFTausMatchedToL1, module_hltPFTausMatchedToL1)
        process.analysisSequence += module_hltPFTausMatchedToL1
        hlt_srcPFTaus = moduleName_hltPFTausMatchedToL1

        from HLTrigger.Phase2HLTPFTaus.PFRecoTauChargedIsoPtSum_cfi import hltPFTauChargedIsoPtSum
        moduleName_hltPFTauChargedIsoPtSum = "%sMatchedToL1" % hlt_srcPFTauSumChargedIso
        module_hltPFTauChargedIsoPtSum = hltPFTauChargedIsoPtSum.clone()
        module_hltPFTauChargedIsoPtSum.PFTauProducer = cms.InputTag(hlt_srcPFTaus)
        module_hltPFTauChargedIsoPtSum.particleFlowSrc = cms.InputTag('particleFlowTmp')
        module_hltPFTauChargedIsoPtSum.vertexSrc = cms.InputTag(hlt_srcVertices)
        hlt_isolation_maxDeltaZ            = None
        hlt_isolation_maxDeltaZToLeadTrack = None
        if hlt_isolation_maxDeltaZOption == "primaryVertex":
          hlt_isolation_maxDeltaZ            =  0.15 # value optimized for offline tau reconstruction at higher pileup expected during LHC Phase-2
          hlt_isolation_maxDeltaZToLeadTrack = -1.   # disabled
        elif hlt_isolation_maxDeltaZOption == "leadTrack":
          hlt_isolation_maxDeltaZ            = -1.   # disabled
          hlt_isolation_maxDeltaZToLeadTrack =  0.15 # value optimized for offline tau reconstruction at higher pileup expected during LHC Phase-2
        else:
          raise ValueError("Invalid parameter hlt_isolation_maxDeltaZOption = '%s' !!" % hlt_isolation_maxDeltaZOption)
        module_hltPFTauChargedIsoPtSum.qualityCuts.isolationQualityCuts.maxDeltaZ = cms.double(hlt_isolation_maxDeltaZ)
        module_hltPFTauChargedIsoPtSum.qualityCuts.isolationQualityCuts.maxDeltaZToLeadTrack = cms.double(hlt_isolation_maxDeltaZToLeadTrack)
        module_hltPFTauChargedIsoPtSum.qualityCuts.primaryVertexSrc = cms.InputTag(hlt_srcVertices)
        module_hltPFTauChargedIsoPtSum.qualityCuts.isolationQualityCuts.minTrackHits = cms.uint32(hlt_isolation_minTrackHits)
        setattr(process, moduleName_hltPFTauChargedIsoPtSum, module_hltPFTauChargedIsoPtSum)
        process.analysisSequence += module_hltPFTauChargedIsoPtSum
        hlt_srcPFTauSumChargedIso = moduleName_hltPFTauChargedIsoPtSum

        moduleName_hltPFTauNeutralIsoPtSum = "%sMatchedToL1" % hlt_srcPFTauSumNeutralIso
        module_hltPFTauNeutralIsoPtSum = module_hltPFTauChargedIsoPtSum.clone(
          ApplyDiscriminationByTrackerIsolation = cms.bool(False),
          ApplyDiscriminationByECALIsolation = cms.bool(True),
          WeightECALIsolation = cms.double(1.)
        )
        setattr(process, moduleName_hltPFTauNeutralIsoPtSum, module_hltPFTauNeutralIsoPtSum)
        process.analysisSequence += module_hltPFTauNeutralIsoPtSum
        hlt_srcPFTauSumNeutralIso = moduleName_hltPFTauNeutralIsoPtSum

        moduleName_hltUpdatedPatTausMatchedToL1 = "%sMatchedToL1" % hlt_srcUpdatedPatTaus
        module_hltUpdatedPatTausMatchedToL1 = cms.EDFilter("PATTauAntiOverlapSelector",
          src = cms.InputTag(hlt_srcUpdatedPatTaus),
          srcNotToBeFiltered = cms.VInputTag('L1HPSPFTausPassingTrigger'),
          dRmin = cms.double(0.3),
          invert = cms.bool(True),
          filter = cms.bool(False)
        )
        setattr(process, moduleName_hltUpdatedPatTausMatchedToL1, module_hltUpdatedPatTausMatchedToL1)
        process.analysisSequence += module_hltUpdatedPatTausMatchedToL1
        hlt_srcUpdatedPatTaus = moduleName_hltUpdatedPatTausMatchedToL1

      moduleName_PFTauAnalyzerSignal_recoSumChargedIso_wrtGenHadTaus = "analyze%ss%sRecoSumChargedIsoWrtGenHadTaus" % (hlt_pfTauLabel, suffix)
      module_PFTauAnalyzerSignal_recoSumChargedIso_wrtGenHadTaus = cms.EDAnalyzer("RecoPFTauAnalyzerSignal",
        srcPFTaus = cms.InputTag(hlt_srcPFTaus),
        srcPFTauDiscriminator = cms.InputTag(hlt_srcPFTauSumChargedIso),
        srcDenominator = cms.InputTag(srcGenHadTaus),
        typeDenominator = cms.string("gen"),                                                                            
        min_pt_denominator = cms.double(45.),
        max_pt_denominator = cms.double(-1.),
        min_pt_numerator = cms.vdouble( 20., 25., 30., 35., 40., 45. ),
        max_pt_numerator = cms.vdouble( -1., -1., -1., -1., -1., -1. ),
        min_absEta = cms.vdouble( -1.,   1.4,   1.4, -1.,    -1.  ),
        max_absEta = cms.vdouble(  1.4,  2.172, 2.4,  2.172,  2.4 ),
        decayModes = cms.vstring("oneProng0Pi0", "oneProng1Pi0", "oneProng2Pi0", "threeProng0Pi0", "threeProng1Pi0", "all"),
        min_leadTrackPt = cms.vdouble(  1.,  2.,  5. ),
        max_leadTrackPt = cms.vdouble( -1., -1., -1. ),
        min_relDiscriminator = cms.vdouble( -1.,   -1.,   -1.,   -1.,   -1.,   -1.,   -1.   ),
        max_relDiscriminator = cms.vdouble( -1.,    0.40,  0.20,  0.10,  0.05,  0.02,  0.01 ),
        min_absDiscriminator = cms.vdouble(),
        max_absDiscriminator = cms.vdouble(),                                                                                
        src_evtWeight = cms.InputTag('lumiScale'),
        dqmDirectory = cms.string("%s/%sAnalyzerSignal%s_recoSumChargedIso_wrtGenHadTaus" % (hlt_srcVertices, hlt_pfTauLabel, suffix))
      )
      setattr(process, moduleName_PFTauAnalyzerSignal_recoSumChargedIso_wrtGenHadTaus, module_PFTauAnalyzerSignal_recoSumChargedIso_wrtGenHadTaus)
      process.analysisSequence += module_PFTauAnalyzerSignal_recoSumChargedIso_wrtGenHadTaus

      moduleName_PFTauAnalyzerSignal_patSumChargedIso_wrtGenHadTaus = "analyze%ss%sPatSumChargedIsoWrtGenHadTaus" % (hlt_pfTauLabel, suffix)
      module_PFTauAnalyzerSignal_patSumChargedIso_wrtGenHadTaus = cms.EDAnalyzer("PATTauAnalyzerSignal",
        srcPFTaus = cms.InputTag(hlt_srcUpdatedPatTaus),
        pfTauDiscriminator = cms.string('chargedIsoPtSum'),
        srcDenominator = cms.InputTag(srcGenHadTaus),
        typeDenominator = cms.string("gen"),         
        min_pt_denominator = cms.double(45.),
        max_pt_denominator = cms.double(-1.),
        min_pt_numerator = cms.vdouble( 20., 25., 30., 35., 40., 45. ),
        max_pt_numerator = cms.vdouble( -1., -1., -1., -1., -1., -1. ),
        min_absEta = cms.vdouble( -1.,   1.4,   1.4, -1.,    -1.  ),
        max_absEta = cms.vdouble(  1.4,  2.172, 2.4,  2.172,  2.4 ),
        decayModes = cms.vstring("oneProng0Pi0", "oneProng1Pi0", "oneProng2Pi0", "threeProng0Pi0", "threeProng1Pi0", "all"),
        min_leadTrackPt = cms.vdouble(  1.,  2.,  5. ),
        max_leadTrackPt = cms.vdouble( -1., -1., -1. ),
        min_relDiscriminator = cms.vdouble( -1.,   -1.,   -1.,   -1.,   -1.,   -1.,   -1.   ),
        max_relDiscriminator = cms.vdouble( -1.,    0.40,  0.20,  0.10,  0.05,  0.02,  0.01 ),
        min_absDiscriminator = cms.vdouble(),
        max_absDiscriminator = cms.vdouble(),                                                                                
        src_evtWeight = cms.InputTag('lumiScale'),
        dqmDirectory = cms.string("%s/%sAnalyzerSignal%s_patSumChargedIso_wrtGenHadTaus" % (hlt_srcVertices, hlt_pfTauLabel, suffix))
      )
      setattr(process, moduleName_PFTauAnalyzerSignal_patSumChargedIso_wrtGenHadTaus, module_PFTauAnalyzerSignal_patSumChargedIso_wrtGenHadTaus)
      process.analysisSequence += module_PFTauAnalyzerSignal_patSumChargedIso_wrtGenHadTaus

      moduleName_PFTauAnalyzerSignal_patDeepTau_wrtGenHadTaus = "analyze%ss%sPatDeepTauWrtGenHadTaus" % (hlt_pfTauLabel, suffix)
      module_PFTauAnalyzerSignal_patDeepTau_wrtGenHadTaus = module_PFTauAnalyzerSignal_patSumChargedIso_wrtGenHadTaus.clone(
        pfTauDiscriminator = cms.string('byDeepTau2017v2VSjetraw'),
        min_relDiscriminator = cms.vdouble(),
        max_relDiscriminator = cms.vdouble(), 
        min_absDiscriminator = cms.vdouble(  0.2599605,  0.4249705,  0.5983682,  0.7848675,  0.8834768,  0.9308689,  0.9573137,  0.9733927 ),
        max_absDiscriminator = cms.vdouble( -1.,        -1.,        -1.,        -1.,        -1.,        -1.,        -1.,        -1.        ),
        dqmDirectory = cms.string("%s/%sAnalyzerSignal%s_patDeepTau_wrtGenHadTaus" % (hlt_srcVertices, hlt_pfTauLabel, suffix))
      )
      setattr(process, moduleName_PFTauAnalyzerSignal_patDeepTau_wrtGenHadTaus, module_PFTauAnalyzerSignal_patDeepTau_wrtGenHadTaus)
      process.analysisSequence += module_PFTauAnalyzerSignal_patDeepTau_wrtGenHadTaus

      from HLTrigger.Phase2HLTPFTaus.PFTauPairProducer_cfi import PFTauPairs
      moduleName_PFTauPairProducer = "hlt%sPairs%s" % (hlt_pfTauLabel, suffix)
      module_PFTauPairProducer = PFTauPairs.clone(
        srcPFTaus = cms.InputTag(hlt_srcPFTaus),
        srcPFTauSumChargedIso = cms.InputTag(hlt_srcPFTauSumChargedIso)
      )
      setattr(process, moduleName_PFTauPairProducer, module_PFTauPairProducer)
      process.analysisSequence += module_PFTauPairProducer

      moduleName_PFTauPairAnalyzer_wrtGenHadTaus = "analyze%sPairs%sWrtGenHadTaus" % (hlt_pfTauLabel, suffix)
      module_PFTauPairAnalyzer_wrtGenHadTaus = cms.EDAnalyzer("RecoPFTauPairAnalyzer",
        srcPFTauPairs = cms.InputTag(moduleName_PFTauPairProducer),
        srcRefTaus = cms.InputTag(srcGenHadTaus),
        min_refTau_pt = cms.double(20.),
        max_refTau_pt = cms.double(1.e+3),                                                                
        min_refTau_absEta = cms.double(-1.),
        max_refTau_absEta = cms.double(2.4),     
        lumiScale = cms.double(1.),                                               
        dqmDirectory = cms.string("%s/%sPairAnalyzer%s_wrtGenHadTaus" % (hlt_srcVertices, hlt_pfTauLabel, suffix))
      )
      setattr(process, moduleName_PFTauPairAnalyzer_wrtGenHadTaus, module_PFTauPairAnalyzer_wrtGenHadTaus)
      process.analysisSequence += module_PFTauPairAnalyzer_wrtGenHadTaus

      if hasMiniAOD:
        moduleName_PFTauAnalyzerSignal_recoSumChargedIso_wrtOfflineTaus = "analyze%ss%sRecoSumChargedIsoWrtOfflineTaus" % (hlt_pfTauLabel, suffix)
        modulePF_PFTauAnalyzerSignal_recoSumChargedIso_wrtOfflineTaus = module_PFTauAnalyzerSignal_recoSumChargedIso_wrtGenHadTaus.clone(
          srcDenominator = cms.InputTag('selectedOfflinePFTaus'),
          typeDenominator = cms.string("offline"),  
          dqmDirectory = cms.string("%s/%sAnalyzerSignal%s_recoSumChargedIso_wrtOfflineTaus" % (hlt_srcVertices, hlt_pfTauLabel, suffix))
        )
        setattr(process, moduleName_PFTauAnalyzerSignal_recoSumChargedIso_wrtOfflineTaus, modulePF_PFTauAnalyzerSignal_recoSumChargedIso_wrtOfflineTaus)
        process.analysisSequence += modulePF_PFTauAnalyzerSignal_recoSumChargedIso_wrtOfflineTaus

        moduleName_PFTauAnalyzerSignal_patSumChargedIso_wrtOfflineTaus = "analyze%ss%sPatSumChargedIsoWrtOfflineTaus" % (hlt_pfTauLabel, suffix)
        modulePF_PFTauAnalyzerSignal_patSumChargedIso_wrtOfflineTaus = module_PFTauAnalyzerSignal_patSumChargedIso_wrtGenHadTaus.clone(
          srcDenominator = cms.InputTag('selectedOfflinePFTaus'),
          typeDenominator = cms.string("offline"),  
          dqmDirectory = cms.string("%s/%sAnalyzerSignal%s_patSumChargedIso_wrtOfflineTaus" % (hlt_srcVertices, hlt_pfTauLabel, suffix))
        )
        setattr(process, moduleName_PFTauAnalyzerSignal_patSumChargedIso_wrtOfflineTaus, modulePF_PFTauAnalyzerSignal_patSumChargedIso_wrtOfflineTaus)
        process.analysisSequence += modulePF_PFTauAnalyzerSignal_patSumChargedIso_wrtOfflineTaus

        moduleName_PFTauAnalyzerSignal_patDeepTau_wrtOfflineTaus = "analyze%ss%sPatDeepTauWrtOfflineTaus" % (hlt_pfTauLabel, suffix)
        modulePF_PFTauAnalyzerSignal_patDeepTau_wrtOfflineTaus = module_PFTauAnalyzerSignal_patDeepTau_wrtGenHadTaus.clone(
          srcDenominator = cms.InputTag('selectedOfflinePFTaus'),
          typeDenominator = cms.string("offline"),  
          dqmDirectory = cms.string("%s/%sAnalyzerSignal%s_patDeepTau_wrtOfflineTaus" % (hlt_srcVertices, hlt_pfTauLabel, suffix))
        )
        setattr(process, moduleName_PFTauAnalyzerSignal_patDeepTau_wrtOfflineTaus, modulePF_PFTauAnalyzerSignal_patDeepTau_wrtOfflineTaus)
        process.analysisSequence += modulePF_PFTauAnalyzerSignal_patDeepTau_wrtOfflineTaus

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
        srcPFTaus = cms.InputTag(hlt_srcPFTaus),
        srcPFTauSumChargedIso = cms.InputTag(hlt_srcPFTauSumChargedIso),
        srcPFTauSumNeutralIso = cms.InputTag(hlt_srcPFTauSumNeutralIso),
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
