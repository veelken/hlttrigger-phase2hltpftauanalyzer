import FWCore.ParameterSet.Config as cms

import os

from HLTrigger.TallinnHLTPFTauAnalyzer.tools.get_suffix import get_suffix
from HLTrigger.Phase2HLTPFTaus.PFTauPairProducer_cfi import PFTauPairs

def addGenPtHatSequence(process, hlt_srcVertices, hlt_isolation_maxDeltaZOptions, hlt_isolation_minTrackHits, lumiScale, matchToL1):

    process.genPtHatSequence = cms.Sequence()

    process.genPtHatAnalzerBeforeCuts = cms.EDAnalyzer("GenPtHatAnalyzer",
      src_genEventInfo = cms.InputTag('generator'),
      src_genJets = cms.InputTag('ak4GenJetsNoNu'),
      src_pileupSummaryInfo = cms.InputTag('addPileupInfo'),
      lumiScale = cms.double(lumiScale),
      dqmDirectory = cms.string("GenPtHatAnalyzer/beforeCuts")
    )
    process.genPtHatSequence += process.genPtHatAnalzerBeforeCuts

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
    if matchToL1:
      suffix += "MatchedToL1"

    moduleName_hltPFTaus = "hltSelectedHpsPFTaus%s" % suffix
    if hasattr(process, moduleName_hltPFTaus):
        module_hltPFTaus = getattr(process, moduleName_hltPFTaus)
        process.analysisSequence += module_hltPFTaus

    moduleName_hltPFTauChargedIsoPtSum = "hltSelectedHpsPFTauChargedIsoPtSum%s" % suffix
    if hasattr(process, moduleName_hltPFTauChargedIsoPtSum):
        module_hltPFTauChargedIsoPtSum = getattr(process, moduleName_hltPFTauChargedIsoPtSum)
        process.analysisSequence += module_hltPFTauChargedIsoPtSum

    process.hltHpsPFTausPassingTrigger = cms.EDProducer("MyPFTauSelector",
      src = cms.InputTag(moduleName_hltPFTaus),
      src_sumChargedIso = cms.InputTag(moduleName_hltPFTauChargedIsoPtSum),
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

    process.genPtHatAnalzerAfterCuts = process.genPtHatAnalzerBeforeCuts.clone(
      dqmDirectory = cms.string("GenPtHatAnalyzer/afterCuts")
    )
    process.genPtHatSequence += process.genPtHatAnalzerAfterCuts

    process.genPtHatPath = cms.Path(process.genPtHatSequence)
