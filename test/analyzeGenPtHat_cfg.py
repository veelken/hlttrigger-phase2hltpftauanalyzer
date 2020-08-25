
import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzeGenPtHat")

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

inputFilePath = '/hdfs/cms/store/user/rdewanje/MinBias_TuneCP5_14TeV-pythia8/HLTConfig_MinBias_TuneCP5_14TeV-pythia8_wOfflineVtx_wDeepTau3/'
#inputFilePath = '/hdfs/cms/store/user/rdewanje/QCD_Pt_30to50_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_30to50_TuneCP5_14TeV_pythia8_wOfflineVtx_wDeepTau3/'
#inputFilePath = None
inputFileNames = []
sampleName = "minbias"
#sampleName = "qcd_pt30to50"
outputFileName = "analyzeGenPtHat_%s_DEBUG.root" % sampleName

##inputFilePath = None
##inputFileNames = $inputFileNames
##sampleName = "$sampleName"
##lumiScale = $lumiScale
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

process.genPtHatAnalzer = cms.EDAnalyzer("GenPtHatAnalyzer",
   src_genEventInfo = cms.InputTag('generator'),
   src_genJets = cms.InputTag('ak4GenJetsNoNu'),
   src_pileupSummaryInfo = cms.InputTag('addPileupInfo'),
   lumiScale = cms.double(1.),
   dqmDirectory = cms.string("GenPtHatAnalyzer/%s" % sampleName)
)
process.analysisSequence += process.genPtHatAnalzer

process.dumpGenPtHat = cms.EDAnalyzer("DumpGenPtHat",
   src_genEventInfo = cms.InputTag('generator'),
   src_genJets = cms.InputTag('ak4GenJetsNoNu'),
   src_pileupSummaryInfo = cms.InputTag('addPileupInfo')
)
##process.analysisSequence += process.dumpGenPtHat

process.load("DQMServices.Core.DQMStore_cfi")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string(outputFileName)
)

process.p = cms.Path(process.analysisSequence + process.savePlots)
