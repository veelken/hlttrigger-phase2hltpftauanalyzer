
import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzePFTausBackground")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000000)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/home/veelken/Phase2HLT/CMSSW_11_1_0/src/HLTrigger/Phase2HLTPFTaus/test/step3_RAW2DIGI_RECO.root'
    ),
##    eventsToProcess = cms.untracked.VEventRange(
##        '1:262:90464'
##    ) 
)

inputFilePath = '/hdfs/cms/store/user/rdewanje/MinBias_TuneCP5_14TeV-pythia8/HLTConfig_MinBias_TuneCP5_14TeV-pythia8_wOfflineVtx_wL1/'
#inputFilePath = None
inputFileNames = []
#processName = "minbias"
processName = "QCD"
lumiScale = 2.8e+7 # 28 MHz
hlt_srcVertices = 'offlinePrimaryVertices'
#hlt_srcVertices = 'hltPhase2PixelVertices'
#hlt_algorithms = [ "hps", "shrinking-cone" ]
hlt_algorithms = [ "hps" ]
hlt_isolation_maxDeltaZOptions = [ "primaryVertex", "leadTrack" ]
hlt_isolation_minTrackHits = 8
l1_useStrips = True
outputFileName = "analyzePFTaus_background_%s_%s_DEBUG.root" % (processName, hlt_srcVertices)

##inputFilePath = None
##inputFileNames = $inputFileNames
##processName = "$processName"
##lumiScale = $lumiScale
##hlt_srcVertices = '$hlt_srcVertices'
##hlt_algorithms = [ "$hlt_algorithm" ]
##hlt_isolation_maxDeltaZOptions = [ "$hlt_isolation_maxDeltaZOption" ]
##hlt_isolation_minTrackHits = $hlt_isolation_minTrackHits
##l1_useStrips = '$l1_useStrips'
##outputFileName = "$outputFileName"

#--------------------------------------------------------------------------------
# set input files
if inputFilePath:
    from HLTrigger.TallinnHLTPFTauAnalyzer.tools.jobTools import getInputFileNames
    print("Searching for input files in path = '%s'" % inputFilePath)
    inputFileNames = getInputFileNames(inputFilePath)
    print("Found %i input files." % len(inputFileNames))
    process.source.fileNames = cms.untracked.vstring(inputFileNames)
elif len(inputFileNames) > 0:
    print("Processing %i input files: %s" % (len(inputFileNames), inputFileNames))
    process.source.fileNames = cms.untracked.vstring(inputFileNames)
#--------------------------------------------------------------------------------

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.analysisSequence = cms.Sequence()

#--------------------------------------------------------------------------------
# CV: fill L1 tau plots
moduleNameBase = "L1HPSPFTauProducer"
moduleLabel = None
if l1_useStrips:
    moduleLabel = "WithStrips"
else:
    moduleLabel = "WithoutStrips"
l1_srcPFTaus = moduleNameBase + moduleLabel + "PF"

moduleName_L1HPSPFTauAnalyzerBackground = "analyzeL1HPSPFTaus" + moduleLabel + "PF"
module_L1HPSPFTauAnalyzerBackground = cms.EDAnalyzer("L1HPSPFTauAnalyzerBackground",
  srcL1PFTaus = cms.InputTag(l1_srcPFTaus),
  dqmDirectory = cms.string("L1HPSPFTauAnalyzerBackground" + moduleLabel + "PF")
)
setattr(process, moduleName_L1HPSPFTauAnalyzerBackground, module_L1HPSPFTauAnalyzerBackground)
process.analysisSequence += getattr(process, moduleName_L1HPSPFTauAnalyzerBackground)

moduleName_L1HPSPFTauPairAnalyzer = "analyzeL1HPSPFTauPairs" + moduleLabel + "PF"
module_L1HPSPFTauPairAnalyzer = cms.EDAnalyzer("L1HPSPFTauPairAnalyzer",
  srcL1PFTaus = cms.InputTag(moduleNameBase + moduleLabel + "PF"),
  srcRefTaus = cms.InputTag(''),
  dqmDirectory = cms.string("L1HPSPFTauPairAnalyzer" + moduleLabel + "PF")
)
setattr(process, moduleName_L1HPSPFTauPairAnalyzer, module_L1HPSPFTauPairAnalyzer)
process.analysisSequence += getattr(process, moduleName_L1HPSPFTauPairAnalyzer)

moduleName_L1HPSPFTauIsolationAnalyzer = "analyzeTallinL1PFTauIsolation" + moduleLabel + "PF"
module_L1HPSPFTauIsolationAnalyzer = cms.EDAnalyzer("L1HPSPFTauIsolationAnalyzer",
  srcL1PFTaus  = cms.InputTag(moduleNameBase + moduleLabel + "PF"),
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

def get_suffix(hlt_srcVertices, hlt_isolation_maxDeltaZOption, hlt_isolation_minTrackHits):
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
    return suffix

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

      suffix = get_suffix(hlt_srcVertices, hlt_isolation_maxDeltaZOption, hlt_isolation_minTrackHits)

      hlt_srcPFTaus = 'hltSelected%ss%s' % (hlt_pfTauLabel, suffix)
      hlt_srcPFTauSumChargedIso = 'hlt%sChargedIsoPtSum%s' % (hlt_pfTauLabel, suffix)
      hlt_srcPFTauSumNeutralIso = 'hlt%sNeutralIsoPtSum%s' % (hlt_pfTauLabel, suffix)
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

      moduleName_PFTauAnalyzerBackground = "analyze%ss%s" % (hlt_pfTauLabel, suffix)
      module_PFTauAnalyzerBackground = cms.EDAnalyzer("RecoPFTauAnalyzerBackground",
        srcPFTaus = cms.InputTag(hlt_srcPFTaus),
        srcPFTauSumChargedIso = cms.InputTag(hlt_srcPFTauSumChargedIso),
        min_pt = cms.double(20.),
        max_absEta = cms.double(2.4), 
        lumiScale = cms.double(lumiScale),
        dqmDirectory = cms.string("%s/%s/%sAnalyzerBackground%s" % (processName, hlt_srcVertices, hlt_pfTauLabel, suffix))
      )
      setattr(process, moduleName_PFTauAnalyzerBackground, module_PFTauAnalyzerBackground)
      process.analysisSequence += module_PFTauAnalyzerBackground

      from HLTrigger.Phase2HLTPFTaus.PFTauPairProducer_cfi import PFTauPairs
      moduleName_PFTauPairProducer = "hlt%sPairs%s" % (hlt_pfTauLabel, suffix)
      module_PFTauPairProducer = PFTauPairs.clone(
        srcPFTaus = cms.InputTag(hlt_srcPFTaus),
        srcPFTauSumChargedIso = cms.InputTag(hlt_srcPFTauSumChargedIso)
      )
      setattr(process, moduleName_PFTauPairProducer, module_PFTauPairProducer)
      process.analysisSequence += module_PFTauPairProducer

      moduleName_PFTauPairAnalyzer = "analyze%sPairs%s" % (hlt_pfTauLabel, suffix)
      module_PFTauPairAnalyzer = cms.EDAnalyzer("RecoPFTauPairAnalyzer",
        srcPFTauPairs = cms.InputTag(moduleName_PFTauPairProducer),
        srcRefTaus = cms.InputTag(''),
        lumiScale = cms.double(lumiScale),
        dqmDirectory = cms.string("%s/%s/%sPairAnalyzer%s" % (processName, hlt_srcVertices, hlt_pfTauLabel, suffix))
      )
      setattr(process, moduleName_PFTauPairAnalyzer, module_PFTauPairAnalyzer)
      process.analysisSequence += module_PFTauPairAnalyzer

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
        lumiScale = cms.double(lumiScale),                    
        dqmDirectory = cms.string("%s/%s/%sIsolationAnalyzer%s" % (processName, hlt_srcVertices, hlt_pfTauLabel, suffix))
      )
      setattr(process, moduleName_PFTauIsolationAnalyzer, module_PFTauIsolationAnalyzer)
      ##process.analysisSequence += module_PFTauIsolationAnalyzer
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# CV: fill histogram of generator PtHat information when running on QCD multijet MC samples
#    (generator PtHat information will be valid for QCD multijet MC samples only)
if processName == "QCD":

    print("Adding GenPtHatAnalyzer.")
    process.genPtHatSequence = cms.Sequence()

    hlt_isolation_maxDeltaZOption = None
    if hlt_srcVertices == 'hltPhase2PixelVertices' and "leadTrack" in hlt_isolation_maxDeltaZOptions:
      hlt_isolation_maxDeltaZOption = "leadTrack"
    else:
      hlt_isolation_maxDeltaZOption = "primaryVertex"

    suffix = get_suffix(hlt_srcVertices, hlt_isolation_maxDeltaZOption, hlt_isolation_minTrackHits)

    process.hltHpsPFTausPassingTrigger = cms.EDProducer("MyPFTauSelector",
      src = cms.InputTag('hltSelectedHpsPFTaus%sMatchedToL1' % suffix),
      src_sumChargedIso = cms.InputTag('hltHpsPFTauChargedIsoPtSum%sMatchedToL1' % suffix),
      min_pt = cms.double(30.),
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
    process.genPtHatSequence += process.hltHpsPFTausPassingTrigger

    process.hltHpsPFTauPairsPassingTrigger = PFTauPairs.clone(
      srcPFTaus = cms.InputTag('hltHpsPFTausPassingTrigger'),
      srcPFTauSumChargedIso = cms.InputTag('')
    )
    process.genPtHatSequence += process.hltHpsPFTauPairsPassingTrigger

    process.hltHpsPFTauPairFilter = cms.EDFilter("CandViewCountFilter",
      src = cms.InputTag('hltHpsPFTauPairsPassingTrigger'),
      minNumber = cms.uint32(1)
    )
    process.genPtHatSequence += process.hltHpsPFTauPairFilter

    process.genPtHatAnalzer = cms.EDAnalyzer("GenPtHatAnalyzer",
      src = cms.InputTag('generator'),
      lumiScale = cms.double(lumiScale_genPtHat),
      dqmDirectory = cms.string("GenPtHatAnalyzer")
    )
    process.genPtHatSequence += process.genPtHatAnalzer

    process.genPtHatPath = cms.Path(process.genPtHatSequence)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# CV: check events that pass HLT ditau trigger with L1 trigger matching,
#     but fail the same HLT ditau trigger in case no L1 trigger matching is applied

process.debugSequence = cms.Sequence()

hlt_isolation_maxDeltaZOption = None
if hlt_srcVertices == 'hltPhase2PixelVertices' and "leadTrack" in hlt_isolation_maxDeltaZOptions:
  hlt_isolation_maxDeltaZOption = "leadTrack"
else:
  hlt_isolation_maxDeltaZOption = "primaryVertex"

suffix = get_suffix(hlt_srcVertices, hlt_isolation_maxDeltaZOption, hlt_isolation_minTrackHits)

process.hltHpsPFTausPassingPtGt125 = cms.EDProducer("MyPFTauSelector",
  src = cms.InputTag('hltSelectedHpsPFTaus%s' % suffix),
  src_sumChargedIso = cms.InputTag('hltHpsPFTauChargedIsoPtSum%s' % suffix),
  min_pt = cms.double(125.),
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
process.debugSequence += process.hltHpsPFTausPassingPtGt125

process.hltHpsPFTauPairsPassingPtGt125 = PFTauPairs.clone(
  srcPFTaus = cms.InputTag('hltHpsPFTausPassingPtGt125'),
  srcPFTauSumChargedIso = cms.InputTag('')
)
process.debugSequence += process.hltHpsPFTauPairsPassingPtGt125

process.hltHpsPFTauPairsPassingPtGt125Veto = cms.EDFilter("MyCandViewCountFilter",
  src = cms.InputTag('hltHpsPFTauPairsPassingPtGt125'),
  minNumber = cms.int32(-1),
  maxNumber = cms.int32(0)
)
process.debugSequence += process.hltHpsPFTauPairsPassingPtGt125Veto

process.hltHpsPFTausMatchedToL1PassingPtGt125 = process.hltHpsPFTausPassingPtGt125.clone(
  src = cms.InputTag('hltSelectedHpsPFTaus%sMatchedToL1' % suffix),
  src_sumChargedIso = cms.InputTag('hltHpsPFTauChargedIsoPtSum%sMatchedToL1' % suffix)
)
process.debugSequence += process.hltHpsPFTausMatchedToL1PassingPtGt125

process.hltHpsPFTauPairsMatchedToL1PassingPtGt125 = PFTauPairs.clone(
  srcPFTaus = cms.InputTag('hltHpsPFTausMatchedToL1PassingPtGt125'),
  srcPFTauSumChargedIso = cms.InputTag('')
)
process.debugSequence += process.hltHpsPFTauPairsMatchedToL1PassingPtGt125

process.hltHpsPFTauPairsMatchedToL1PassingPtGt125Filter = cms.EDFilter("MyCandViewCountFilter",
  src = cms.InputTag('hltHpsPFTauPairsMatchedToL1PassingPtGt125'),
  minNumber = cms.int32(1),
  maxNumber = cms.int32(-1)
)
process.debugSequence += process.hltHpsPFTauPairsMatchedToL1PassingPtGt125Filter

process.dumpOnlineHpsPFTaus = cms.EDAnalyzer("DumpRecoPFTaus",
  src = cms.InputTag('hltSelectedHpsPFTaus%s' % suffix),
  src_sumChargedIso = cms.InputTag('hltHpsPFTauChargedIsoPtSum%s' % suffix),
  src_discriminators = cms.VInputTag()
)
process.debugSequence += process.dumpOnlineHpsPFTaus

process.dumpOnlineHpsPFTausMatchedToL1 = cms.EDAnalyzer("DumpRecoPFTaus",
  src = cms.InputTag('hltSelectedHpsPFTaus%sMatchedToL1' % suffix),
  src_sumChargedIso = cms.InputTag('hltHpsPFTauChargedIsoPtSum%sMatchedToL1' % suffix),
  src_discriminators = cms.VInputTag()
)
process.debugSequence += process.dumpOnlineHpsPFTausMatchedToL1

import os
process.runLumiSectionEventNumberAnalyzer = cms.EDAnalyzer("RunLumiSectionEventNumberAnalyzer",
  output = cms.string(os.path.join("/home/veelken/Phase2HLT/CMSSW_11_1_0/src/HLTrigger/TallinnHLTPFTauAnalyzer/test/DEBUG/", outputFileName.replace(".root", ".txt"))),
  separator = cms.string(":")
)
process.debugSequence += process.runLumiSectionEventNumberAnalyzer

process.debugPath = cms.Path(process.debugSequence)
#--------------------------------------------------------------------------------

process.load("DQMServices.Core.DQMStore_cfi")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string(outputFileName)
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
