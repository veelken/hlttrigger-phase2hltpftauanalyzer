
import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzeTallinnL1PFTausBackground")

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
##inputFilePath = ''
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

#--------------------------------------------------------------------------------
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

    moduleName_PFTauAnalyzerBackground = "analyze%ss%s" % (pfTauLabel, suffix)
    module_PFTauAnalyzerBackground = cms.EDAnalyzer("RecoPFTauAnalyzerBackground",
      srcPFTaus = cms.InputTag('hltSelected%ss%s' % (pfTauLabel, suffix)),
      srcPFTauSumChargedIso = cms.InputTag('hlt%sChargedIsoPtSum%s' % (pfTauLabel, suffix)),
      dqmDirectory = cms.string("%sAnalyzerBackground%s" % (pfTauLabel, suffix))
    )
    setattr(process, moduleName_PFTauAnalyzerBackground, module_PFTauAnalyzerBackground)
    process.analysisSequence += module_PFTauAnalyzerBackground

    moduleName_PFTauPairAnalyzer = "analyze%sPairs%s" % (pfTauLabel, suffix)
    module_PFTauPairAnalyzer = cms.EDAnalyzer("RecoPFTauPairAnalyzer",
      srcPFTaus = cms.InputTag('hltSelected%ss%s' % (pfTauLabel, suffix)),
      srcPFTauSumChargedIso = cms.InputTag('hlt%sChargedIsoPtSum%s' % (pfTauLabel, suffix)),
      srcRefTaus = cms.InputTag(''),
      dqmDirectory = cms.string("%sPairAnalyzer%s" % (pfTauLabel, suffix))
    )
    setattr(process, moduleName_PFTauPairAnalyzer, module_PFTauPairAnalyzer)
    process.analysisSequence += module_PFTauPairAnalyzer

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
    outputFileName = cms.string('analyzePFTaus_background_2020Jun03.root')
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

#process.options = cms.untracked.PSet(
#    wantSummary = cms.untracked.bool(True)
#)
