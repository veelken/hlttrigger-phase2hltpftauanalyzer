
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

inputFilePath = '/hdfs/cms/store/user/rdewanje/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5/HLTConfig_VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5_wOfflineVtx_wDeepTau2/'
processName = "qqH_htt"
hlt_srcVertices = 'offlinePrimaryVertices'
#hlt_srcVertices = 'hltPhase2PixelVertices'
#hlt_algorithms = [ "hps", "shrinking-cone" ]
hlt_algorithms = [ "hps" ]
hlt_isolation_maxDeltaZOptions = [ "primaryVertex", "leadTrack" ]
hlt_isolation_minTrackHits = 8
outputFileName = "analyzePFTaus_signal_%s_%s_DEBUG.root" % (processName, hlt_srcVertices)

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

hlt_algorithm = "hps"
hlt_pfTauLabel = "HpsPFTau"

hlt_isolation_maxDeltaZOption = "primaryVertex"

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

process.p = cms.Path(process.analysisSequence)
