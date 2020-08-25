import FWCore.ParameterSet.Config as cms

min_pT_hat = 30.0 # CV: minimum threshold on pT_hat for QCD MC samples

minbiasMCFilterSequence = cms.Sequence()

ak4SelGenJets = cms.EDFilter("CandViewSelector",
  src = cms.InputTag('ak4GenJetsNoNu'),
  cut = cms.string("abs(eta) < 2.4 & pt > %2.1f" % min_pT_hat),
  filter = cms.bool(False)
)
minbiasMCFilterSequence += ak4SelGenJets

minbiasMCFilter = cms.EDFilter("MyCandViewCountFilter",
  src = cms.InputTag('ak4SelGenJets'),
  maxNumber = cms.int32(2)
)
minbiasMCFilterSequence += minbiasMCFilter
