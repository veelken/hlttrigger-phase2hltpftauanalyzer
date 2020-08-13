
import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzePFTausBackground")

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
    ),
##    eventsToProcess = cms.untracked.VEventRange(
##        '1:262:90464'
##    ) 
)

inputFilePath = '/hdfs/cms/store/user/rdewanje/MinBias_TuneCP5_14TeV-pythia8/HLTConfig_MinBias_TuneCP5_14TeV-pythia8_wOfflineVtx_wDeepTau2/'
#inputFilePath = None
inputFileNames = []
#processName = "minbias"
processName = "QCD"
#sampleName = "minbias"
sampleName = "qcd_pt30to50"
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
##sampleName = "$sampleName"
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
# CV: remove minbias MC events that overlap with QCD MC samples

if processName == "minbias":
    process.load("HLTrigger.TallinnHLTPFTauAnalyzer.minbiasMCFilter_cff")
    process.analysisSequence += process.minbiasMCFilterSequence
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# CV: fill HLT tau plots

from HLTrigger.TallinnHLTPFTauAnalyzer.tools.get_suffix import get_suffix

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

    #----------------------------------------------------------------------------
    # CV: add DeepTau tau ID discriminator
    from HLTrigger.TallinnHLTPFTauAnalyzer.tools.addDeepTauDiscriminator import addDeepTauDiscriminator
    hlt_srcPFTaus = 'hltSelected%ss%s' % (hlt_pfTauLabel, suffix)
    hlt_srcPFJets = 'hlt%sAK4PFJets%s' % (hlt_pfTauLabel, suffix)
    deepTauSequenceName = "hltDeep%sSequence%s" % (hlt_pfTauLabel, suffix)
    deepTauSequence = addDeepTauDiscriminator(process, hlt_srcPFTaus, hlt_srcPFJets, hlt_srcVertices, hlt_pfTauLabel, suffix, deepTauSequenceName)
    process.analysisSequence += deepTauSequence
    #----------------------------------------------------------------------------

    moduleName_PFTauAnalyzerBackground_sumChargedIso = "analyze%ss%sSumChargedIso" % (hlt_pfTauLabel, suffix)
    module_PFTauAnalyzerBackground_sumChargedIso = cms.EDAnalyzer("RecoPFTauAnalyzerBackground",
      srcPFTaus = cms.InputTag('hltSelected%ss%s' % (hlt_pfTauLabel, suffix)),
      srcPFTauDiscriminator = cms.InputTag('hltSelected%sChargedIsoPtSum%s' % (hlt_pfTauLabel, suffix)),
      min_pt = cms.double(20.),
      max_pt = cms.double(-1.),
      min_absEta = cms.vdouble( -1.,   1.4,   1.4, -1.,    -1.  ),
      max_absEta = cms.vdouble(  1.4,  2.172, 2.4,  2.172,  2.4 ),
      min_leadTrackPt = cms.vdouble(  1.,  2.,  5. ),
      max_leadTrackPt = cms.vdouble( -1., -1., -1. ),
      min_relDiscriminator = cms.vdouble( -1.,   -1.,   -1.,   -1.,   -1.,   -1.,   -1.   ),
      max_relDiscriminator = cms.vdouble( -1.,    0.40,  0.20,  0.10,  0.05,  0.02,  0.01 ),
      min_absDiscriminator = cms.vdouble(),
      max_absDiscriminator = cms.vdouble(),
      min_dzValues = cms.vdouble( -1.  ),
      max_dzValues = cms.vdouble(  0.2 ),
      lumiScale = cms.double(lumiScale),
      dqmDirectory = cms.string("%s/%s/%sAnalyzerBackground%s_sumChargedIso" % (processName, hlt_srcVertices, hlt_pfTauLabel, suffix))
    )
    setattr(process, moduleName_PFTauAnalyzerBackground_sumChargedIso, module_PFTauAnalyzerBackground_sumChargedIso)
    process.analysisSequence += module_PFTauAnalyzerBackground_sumChargedIso

    moduleName_PFTauAnalyzerBackground_deepTau = "analyze%ss%sDeepTau" % (hlt_pfTauLabel, suffix)
    module_PFTauAnalyzerBackground_deepTau = cms.EDAnalyzer("PATTauAnalyzerBackground",
      srcPFTaus = cms.InputTag('hltUpdatedPat%ss%s' % (hlt_pfTauLabel, suffix)),
      pfTauDiscriminator = cms.string('byDeepTau2017v2VSjetraw'),
      min_pt = cms.double(20.),
      max_pt = cms.double(-1.),
      min_absEta = cms.vdouble( -1.,   1.4,   1.4, -1.,    -1.  ),
      max_absEta = cms.vdouble(  1.4,  2.172, 2.4,  2.172,  2.4 ),
      min_leadTrackPt = cms.vdouble(  1.,  2.,  5. ),
      max_leadTrackPt = cms.vdouble( -1., -1., -1. ),
      min_relDiscriminator = cms.vdouble(),
      max_relDiscriminator = cms.vdouble(), 
      min_absDiscriminator = cms.vdouble( -1.,        -1.,        -1.,        -1.,        -1.,        -1.,        -1.,        -1.        ),
      max_absDiscriminator = cms.vdouble(  0.2599605,  0.4249705,  0.5983682,  0.7848675,  0.8834768,  0.9308689,  0.9573137,  0.9733927 ),
      min_dzValues = cms.vdouble( -1.  ),
      max_dzValues = cms.vdouble(  0.2 ),
      lumiScale = cms.double(lumiScale),
      dqmDirectory = cms.string("%s/%s/%sAnalyzerBackground%s_deepTau" % (processName, hlt_srcVertices, hlt_pfTauLabel, suffix))
    )
    setattr(process, moduleName_PFTauAnalyzerBackground_deepTau, module_PFTauAnalyzerBackground_deepTau)
    process.analysisSequence += module_PFTauAnalyzerBackground_deepTau

    from HLTrigger.Phase2HLTPFTaus.PFTauPairProducer_cfi import PFTauPairs
    moduleName_PFTauPairProducer = "hlt%sPairs%s" % (hlt_pfTauLabel, suffix)
    module_PFTauPairProducer = PFTauPairs.clone(
      srcPFTaus = cms.InputTag('hltSelected%ss%s' % (hlt_pfTauLabel, suffix)),
      srcPFTauSumChargedIso = cms.InputTag('hltSelected%sChargedIsoPtSum%s' % (hlt_pfTauLabel, suffix))
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
      srcPFTauSumChargedIso = cms.InputTag('hltSelected%sChargedIsoPtSum%s' % (hlt_pfTauLabel, suffix)),
      srcPFTauSumNeutralIso = cms.InputTag('hltSelected%sNeutralIsoPtSum%s' % (hlt_pfTauLabel, suffix)),
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
    from HLTrigger.TallinnHLTPFTauAnalyzer.tools.addGenPtHatSequence import addGenPtHatSequence
    addGenPtHatSequence(process, hlt_srcVertices, hlt_isolation_maxDeltaZOptions, hlt_isolation_minTrackHits, lumiScale, False)
    process.genPtHatPath = cms.Path(process.genPtHatSequence)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# CV: fill histogram of generator-level jet pT distribution 
#     to check "stitching" of QCD samples in bins of generator PtHat
process.genJetAnalyzer = cms.EDAnalyzer("GenJetAnalyzer",
  src = cms.InputTag('ak4GenJetsNoNu'),
  lumiScale = cms.double(lumiScale),                    
  dqmDirectory = cms.string("%s/GenJetAnalyzer" % sampleName)
)
process.analysisSequence += process.genJetAnalyzer
#--------------------------------------------------------------------------------

process.load("DQMServices.Core.DQMStore_cfi")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string(outputFileName)
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
