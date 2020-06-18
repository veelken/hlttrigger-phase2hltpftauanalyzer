
import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzePFTausSignal")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100000)
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

for algorithm in [ "hps", "shrinking-cone" ]:
  pfTauLabel = None
  if algorithm == "shrinking-cone":
    pfTauLabel = "PFTau"
  elif algorithm == "hps":
    pfTauLabel = "HpsPFTau"
  else:
    raise ValueError("Invalid parameter algorithm = '%s' !!" % algorithm)

  for srcVertices in [ "offlinePrimaryVertices", "hltPhase2PixelVertices", "hltPhase2TrimmedPixelVertices" ]:
    suffix = None
    if srcVertices == "offlinePrimaryVertices":
      suffix = "WithOfflineVertices"
    elif srcVertices == "hltPhase2PixelVertices":
      suffix = "WithOnlineVertices"
    elif srcVertices == "hltPhase2TrimmedPixelVertices":
      suffix = "WithOnlineVerticesTrimmed"
    else:
      raise ValueError("Invalid parameter srcVertices = '%s' !!" % srcVertices)

    moduleName_wrtGenHadTaus = "analyze%sResponse%sWrtGenHadTaus" % (pfTauLabel, suffix)
    module_wrtGenHadTaus = cms.EDAnalyzer("RecoPFTauResponseAnalyzer",
      srcPFTaus = cms.InputTag('hltSelected%ss%s' % (pfTauLabel, suffix)),
      srcPFTauSumChargedIso = cms.InputTag('hlt%sChargedIsoPtSum%s' % (pfTauLabel, suffix)),
      srcRefTaus = cms.InputTag('offlineMatchedGenHadTaus'),
      typeRefTaus = cms.string("gen"),                                                                            
      dqmDirectory = cms.string("%sResponseAnalyzer%s_wrtGenHadTaus" % (pfTauLabel, suffix))
    )
    setattr(process, moduleName_wrtGenHadTaus, module_wrtGenHadTaus)
    process.analysisSequence += module_wrtGenHadTaus

    moduleName_wrtOfflineTaus = "analyze%sResponse%sWrtOfflineTaus" % (pfTauLabel, suffix)
    module_wrtOfflineTaus = cms.EDAnalyzer("RecoPFTauResponseAnalyzer",
      srcPFTaus = cms.InputTag('hltSelected%ss%s' % (pfTauLabel, suffix)),
      srcPFTauSumChargedIso = cms.InputTag('hlt%sChargedIsoPtSum%s' % (pfTauLabel, suffix)),
      srcRefTaus = cms.InputTag('selectedOfflinePFTaus'),
      typeRefTaus = cms.string("offline"),                                                                   
      dqmDirectory = cms.string("%sResponseAnalyzer%s_wrtOfflineTaus" % (pfTauLabel, suffix))
    )
    setattr(process, moduleName_wrtOfflineTaus, module_wrtOfflineTaus)
    process.analysisSequence += module_wrtOfflineTaus
#--------------------------------------------------------------------------------

process.load("DQMServices.Core.DQMStore_cfi")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('analyzePFTauResponse_signal_2020Jun18.root')
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
