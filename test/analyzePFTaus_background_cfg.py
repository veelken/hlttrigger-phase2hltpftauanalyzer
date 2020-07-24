
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
outputFileName = "analyzePFTaus_background_%s_%s_DEBUG.root" % (processName, hlt_srcVertices)

##inputFilePath = None
##inputFileNames = $inputFileNames
##processName = "$processName"
##lumiScale = $lumiScale
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
    process.source.fileNames = cms.untracked.vstring(inputFileNames)
elif len(inputFileNames) > 0:
    print("Processing %i input files: %s" % (len(inputFileNames), inputFileNames))
    process.source.fileNames = cms.untracked.vstring(inputFileNames)
#--------------------------------------------------------------------------------

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.analysisSequence = cms.Sequence()

#--------------------------------------------------------------------------------
# CV: fill HLT tau plots

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

  hlt_pfTauLabel = None
  if hlt_algorithm == "shrinking-cone":
    hlt_pfTauLabel = "PFTau"
  elif hlt_algorithm == "hps":
    hlt_pfTauLabel = "HpsPFTau"
  else:
    raise ValueError("Invalid parameter hlt_algorithm = '%s' !!" % hlt_algorithm)

  for hlt_isolation_maxDeltaZOption in hlt_isolation_maxDeltaZOptions:

    suffix = get_suffix(hlt_srcVertices, hlt_isolation_maxDeltaZOption, hlt_isolation_minTrackHits)

    moduleName_PFTauAnalyzerBackground = "analyze%ss%s" % (hlt_pfTauLabel, suffix)
    module_PFTauAnalyzerBackground = cms.EDAnalyzer("RecoPFTauAnalyzerBackground",
      srcPFTaus = cms.InputTag('hltSelected%ss%s' % (hlt_pfTauLabel, suffix)),
      srcPFTauSumChargedIso = cms.InputTag('hlt%sChargedIsoPtSum%s' % (hlt_pfTauLabel, suffix)),
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
      srcPFTaus = cms.InputTag('hltSelected%ss%s' % (hlt_pfTauLabel, suffix)),
      srcPFTauSumChargedIso = cms.InputTag('hlt%sChargedIsoPtSum%s' % (hlt_pfTauLabel, suffix))
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
    elif "primaryVertex" in hlt_isolation_maxDeltaZOptions:
      hlt_isolation_maxDeltaZOption = "primaryVertex"
    elif "leadTrack" in hlt_isolation_maxDeltaZOptions:
      hlt_isolation_maxDeltaZOption = "leadTrack"
    else:
      raise ValueError("Invalid parameter hlt_isolation_maxDeltaZOptions = %s !!" % hlt_isolation_maxDeltaZOptions)
    suffix = get_suffix(hlt_srcVertices, hlt_isolation_maxDeltaZOption, hlt_isolation_minTrackHits)

    module_hltPFTausMatchedToL1 = getattr(process, "hltSelected%ss%sMatchedToL1" % (hlt_pfTauLabel, suffix))
    process.analysisSequence += module_hltPFTausMatchedToL1

    module_hltPFTauChargedIsoPtSum = getattr(process, "hlt%sChargedIsoPtSum%sMatchedToL1" % (hlt_pfTauLabel, suffix))
    process.analysisSequence += module_hltPFTauChargedIsoPtSum

    process.hltHpsPFTausPassingTrigger = cms.EDProducer("MyPFTauSelector",
      src = cms.InputTag('hltSelectedHpsPFTaus%s' % suffix),
      src_sumChargedIso = cms.InputTag('hltHpsPFTauChargedIsoPtSum%s' % suffix),
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
      lumiScale = cms.double(lumiScale),
      dqmDirectory = cms.string("GenPtHatAnalyzer")
    )
    process.genPtHatSequence += process.genPtHatAnalzer
  
    process.genPtHatPath = cms.Path(process.genPtHatSequence)
#--------------------------------------------------------------------------------

process.load("DQMServices.Core.DQMStore_cfi")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string(outputFileName)
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
