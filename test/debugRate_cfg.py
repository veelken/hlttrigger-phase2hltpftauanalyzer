
import FWCore.ParameterSet.Config as cms

process = cms.Process("debugRate")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000000)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/home/veelken/Phase2HLT/CMSSW_11_1_0_pre6/src/HLTTrigger/Phase2HLTPFTaus/test/step3_RAW2DIGI_RECO.root'
    ),
##    eventsToProcess = cms.untracked.VEventRange(
##        '1:262:90464'
##    ) 
)

#--------------------------------------------------------------------------------
# set input files

import os
import re

inputFilePaths = [
    '/hdfs/cms/store/user/rdewanje/MinBias_TuneCP5_14TeV-pythia8/HLTConfig_Christian_MINBIAS_Phase2HLTTDRWinter20_PU200_CMSSW_11_1_0_pre6_numCores8_maxMem16kMB_T2_EE_Estonia_blacklist/200612_212910/0000/',
    '/hdfs/cms/store/user/rdewanje/MinBias_TuneCP5_14TeV-pythia8/HLTConfig_Christian_MINBIAS_Phase2HLTTDRWinter20_PU200_CMSSW_11_1_0_pre6_numCores8_maxMem16kMB_T2_EE_Estonia_blacklist/200612_212910/0001/'
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

from HLTTrigger.Phase2HLTPFTaus.PFTauPairProducer_cfi import PFTauPairs
process.hltHpsPFTauPairsWithOfflineVertices = PFTauPairs.clone(
    srcPFTaus = cms.InputTag('hltSelectedHpsPFTausWithOfflineVertices'),
    srcPFTauSumChargedIso = cms.InputTag('hltHpsPFTauChargedIsoPtSumWithOfflineVertices')
)
process.analysisSequence += process.hltHpsPFTauPairsWithOfflineVertices

process.hltSelectedHpsPFTauPairsWithOfflineVertices = cms.EDFilter("CandViewSelector",
    src = cms.InputTag('hltHpsPFTauPairsWithOfflineVertices'),
    cut = cms.string("leadPFTau.pt > 30.0 & abs(leadPFTau.eta) < 2.4 & leadPFTau_sumChargedIso < (0.10*leadPFTau.pt) & subleadPFTau.pt > 30.0 & abs(subleadPFTau.eta) < 2.4 & subleadPFTau_sumChargedIso < (0.10*subleadPFTau.pt)"),
    filter = cms.bool(False)
)
process.analysisSequence += process.hltSelectedHpsPFTauPairsWithOfflineVertices

process.zeroHpsPFTauPairsWithOfflineVertices = cms.EDFilter("PATCandViewCountFilter",
    src = cms.InputTag('hltSelectedHpsPFTauPairsWithOfflineVertices'),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(0)
)
process.analysisSequence += process.zeroHpsPFTauPairsWithOfflineVertices

process.hltHpsPFTauPairsWithOnlineVertices = process.hltHpsPFTauPairsWithOfflineVertices.clone(
    srcPFTaus = cms.InputTag('hltSelectedHpsPFTausWithOnlineVertices'),
    srcPFTauSumChargedIso = cms.InputTag('hltHpsPFTauChargedIsoPtSumWithOnlineVertices')
)
process.analysisSequence += process.hltHpsPFTauPairsWithOnlineVertices

process.hltSelectedHpsPFTauPairsWithOnlineVertices = process.hltSelectedHpsPFTauPairsWithOfflineVertices.clone(
    src = cms.InputTag('hltHpsPFTauPairsWithOnlineVertices')
)
process.analysisSequence += process.hltSelectedHpsPFTauPairsWithOnlineVertices

process.oneHpsPFTauPairWithOnlineVertices = cms.EDFilter("PATCandViewCountFilter",
    src = cms.InputTag('hltSelectedHpsPFTauPairsWithOnlineVertices'),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1000)
)
process.analysisSequence += process.oneHpsPFTauPairWithOnlineVertices

process.dumpOfflineVertices = cms.EDAnalyzer("DumpRecoVertices",
    src = cms.InputTag('offlinePrimaryVertices')
)
process.analysisSequence += process.dumpOfflineVertices

process.dumpHpsPFTauPairsWithOfflineVertices = cms.EDAnalyzer("DumpRecoPFTauPairs",
    src = cms.InputTag('hltHpsPFTauPairsWithOfflineVertices')
)
process.analysisSequence += process.dumpHpsPFTauPairsWithOfflineVertices

process.dumpOnlineVertices = process.dumpOfflineVertices.clone(
    src = cms.InputTag('hltPhase2PixelVertices')
)
process.analysisSequence += process.dumpOnlineVertices

process.dumpHpsPFTauPairsWithOnlineVertices = process.dumpHpsPFTauPairsWithOfflineVertices.clone(
    src = cms.InputTag('hltHpsPFTauPairsWithOnlineVertices')
)
process.analysisSequence += process.dumpHpsPFTauPairsWithOnlineVertices

process.debugRate_path = cms.Path(process.analysisSequence)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

# Output definition
process.RECOoutput = cms.OutputModule("PoolOutputModule",
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('debugRate_RECO.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('debugRate_path')
    ),
    outputCommands = cms.untracked.vstring(
        'keep *',
    )
)

process.RECOoutput_step = cms.EndPath(process.RECOoutput)
process.endjob_step = cms.EndPath(process.endOfProcess)
