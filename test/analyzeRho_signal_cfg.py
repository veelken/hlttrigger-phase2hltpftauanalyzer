
import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzeRho")

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

#--------------------------------------------------------------------------------
# set input files
##
##import os
##import re
##
##inputFilePath = '/hdfs/cms/store/user/sbhowmik/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack/PhaseIIMTDTDRAutumn18MiniAOD_20190524/190524_111901/0000/'
##
##inputFile_regex = r"[a-zA-Z0-9_/:.-]*NTuple_TallinnL1PFTauProducer_[a-zA-Z0-9-_]+.root"
##
# check if name of inputFile matches regular expression
##inputFileNames = []
##files = [ "".join([ "file:", inputFilePath, file ]) for file in os.listdir(inputFilePath) ]
##for file in files:
##    inputFile_matcher = re.compile(inputFile_regex)
##    if inputFile_matcher.match(file):
##        inputFileNames.append(file)
##print "inputFileNames = %s" % inputFileNames 
##
##process.source.fileNames = cms.untracked.vstring(inputFileNames)
#--------------------------------------------------------------------------------

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.analysisSequence = cms.Sequence()

# CV: rho and rhoNeutral not yet separately reconstructed and stored in ROOT file,
#     which is why the same collections are used for "rho" and "rhoNeutral" inputs

process.analyzeRho = cms.EDAnalyzer("RhoAnalyzer",
  srcRho = cms.InputTag('hltKT6PFJets:rho'),
  srcRhoNeutral = cms.InputTag('hltKT6PFJets:rho'),
  dqmDirectory = cms.string("RhoAnalyzer"),
)
process.analysisSequence += process.analyzeRho

process.load("DQMServices.Core.DQMStore_cfi")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('analyzeRho_signal_2020Jun03.root')
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
