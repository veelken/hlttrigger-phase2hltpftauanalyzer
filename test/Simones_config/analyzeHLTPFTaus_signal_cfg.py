
import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzeHLTPFTausSignal")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/home/ram/HLT_TAU_PHASE2/CMSSW_11_1_0_pre3_HLTTauConfigAnnCatherine/src/VBFHTAUTAU/hltoutput_hlt_phase2_TDR_Winter20_VBF_HTauTau_PU200.root'
    )
)

sample = "qqH" # SM VBF Higgs->tautau
#sample = "ggH" #  SM Higgs->tautau produced via gluon fusion 

#--------------------------------------------------------------------------------
# set input files
##
##import os
##import re
##
##inputFilePath = None
##if sample == "qqH":
##    inputFilePath = '/hdfs/cms/store/user/sbhowmik/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack/PhaseIIMTDTDRAutumn18MiniAOD_20190617/190618_084235/0000/'
##elif sample == "ggH":
##    inputFilePath = '/hdfs/cms/store/user/sbhowmik/GluGluHToTauTau_M125_14TeV_powheg_pythia8/GluGluHToTauTau_PhaseIIMTDTDRAutumn18MiniAOD_20190617/190618_090007/0000/'
##else:
##    raise ValueError("Invalid sample = '%s' !!" % sample)
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
#process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2023_realistic_v1', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.analysisSequence = cms.Sequence()

#--------------------------------------------------------------------------------
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

process.analyzeHLTPFTausWrtGenHadTaus = cms.EDAnalyzer("TallinnL1PFTauAnalyzerSignal",
  srcNumerator = cms.InputTag('hltHpsSelectedPFTausTrackPt1MediumChargedIsolationReg'),
  srcDenominator = cms.InputTag('offlineMatchedGenHadTaus'),
  typeDenominator = cms.string("gen"),                                                                            
  dqmDirectory = cms.string("hltPFTauAnalyzerSignalWrtGenHadTaus")
)
process.analysisSequence += process.analyzeHLTPFTausWrtGenHadTaus

process.analyzeHLTPFTauPairsWrtGenHadTaus = cms.EDAnalyzer("TallinnL1PFTauPairAnalyzer",
  srcL1PFTaus = cms.InputTag('hltHpsSelectedPFTausTrackPt1MediumChargedIsolationReg'),
  srcRefTaus = cms.InputTag('offlineMatchedGenHadTaus'),
  min_refTau_pt = cms.double(20.),
  max_refTau_pt = cms.double(1.e+3),                                                                
  min_refTau_absEta = cms.double(-1.),
  max_refTau_absEta = cms.double(2.4),                                                                
  dqmDirectory = cms.string("hltPFTauPairAnalyzerWrtGenHadTaus")
)
process.analysisSequence += process.analyzeHLTPFTauPairsWrtGenHadTaus

process.analyzeHLTPFTausWrtOfflineTaus = cms.EDAnalyzer("TallinnL1PFTauAnalyzerSignal",
  srcNumerator = cms.InputTag('hltHpsSelectedPFTausTrackPt1MediumChargedIsolationReg'),
  srcDenominator = cms.InputTag('selectedOfflinePFTaus'),
  typeDenominator = cms.string("gen"),                                                                            
  dqmDirectory = cms.string("hltPFTauAnalyzerSignalWrtOfflineTaus")
)
process.analysisSequence += process.analyzeHLTPFTausWrtOfflineTaus

process.analyzeHLTPFTauPairsWrtOfflineTaus = cms.EDAnalyzer("TallinnL1PFTauPairAnalyzer",
  srcL1PFTaus = cms.InputTag('hltHpsSelectedPFTausTrackPt1MediumChargedIsolationReg'),
  srcRefTaus = cms.InputTag('selectedOfflinePFTaus'),
  min_refTau_pt = cms.double(20.),
  max_refTau_pt = cms.double(1.e+3),                                                                
  min_refTau_absEta = cms.double(-1.),
  max_refTau_absEta = cms.double(2.4),                                                                
  dqmDirectory = cms.string("hltPFTauPairAnalyzerWrtOfflineTaus")
)
process.analysisSequence += process.analyzeHLTPFTauPairsWrtOfflineTaus
#--------------------------------------------------------------------------------

process.DQMStore = cms.Service("DQMStore")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('hltPFTauAnalyzer_signal_%s_2020Jun02.root' % sample)
)

process.p = cms.Path(process.analysisSequence + process.savePlots)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
