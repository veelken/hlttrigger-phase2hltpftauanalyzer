import FWCore.ParameterSet.Config as cms

from HLTrigger.TallinnHLTPFTauAnalyzer.tools.get_suffix import get_suffix
from HLTrigger.Phase2HLTPFTaus.PFTauPairProducer_cfi import PFTauPairs

def addDebugSequenceBackground(process, hlt_srcVertices, hlt_isolation_maxDeltaZOptions, hlt_isolation_minTrackHits, outputFileName):

    process.debugSequence = cms.Sequence()

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

    src_hltPFTaus = "hltSelectedHpsPFTaus%s" % suffix
    src_hltPFTauChargedIsoPtSum = "hltSelectedHpsPFTauChargedIsoPtSum%s" % suffix
 
    src_hltPFTaus_matchedToL1 = "hltSelectedHpsPFTaus%sMatchedToL1" % suffix
    src_hltPFTauChargedIsoPtSum_matchedToL1 = "hltSelectedHpsPFTauChargedIsoPtSum%sMatchedToL1" % suffix

    process.hltHpsPFTausPassingPtGt125 = cms.EDProducer("MyPFTauSelector",
      src = cms.InputTag(src_hltPFTaus),
      src_sumChargedIso = cms.InputTag(src_hltPFTauChargedIsoPtSum),
      min_pt = cms.double(125.),
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
    process.debugSequence += process.hltHpsPFTausPassingPtGt125

    process.hltHpsPFTauPairsPassingPtGt125 = PFTauPairs.clone(
      srcPFTaus = cms.InputTag('hltHpsPFTausPassingPtGt125'),
      srcPFTauSumChargedIso = cms.InputTag('')
    )
    process.debugSequence += process.hltHpsPFTauPairsPassingPtGt125

    process.hltHpsPFTauPairsPassingPtGt125Veto = cms.EDFilter("MyCandViewCountFilter",
      src = cms.InputTag('hltHpsPFTauPairsPassingPtGt125'),
      minNumber = cms.int32(-1),
      maxNumber = cms.int32(0)
    )
    process.debugSequence += process.hltHpsPFTauPairsPassingPtGt125Veto

    process.hltHpsPFTausMatchedToL1PassingPtGt125 = process.hltHpsPFTausPassingPtGt125.clone(
      src = cms.InputTag(src_hltPFTaus_matchedToL1),
      src_sumChargedIso = cms.InputTag(src_hltPFTauChargedIsoPtSum_matchedToL1)
    )
    process.debugSequence += process.hltHpsPFTausMatchedToL1PassingPtGt125

    process.hltHpsPFTauPairsMatchedToL1PassingPtGt125 = PFTauPairs.clone(
      srcPFTaus = cms.InputTag('hltHpsPFTausMatchedToL1PassingPtGt125'),
      srcPFTauSumChargedIso = cms.InputTag('')
    )
    process.debugSequence += process.hltHpsPFTauPairsMatchedToL1PassingPtGt125

    process.hltHpsPFTauPairsMatchedToL1PassingPtGt125Filter = cms.EDFilter("MyCandViewCountFilter",
      src = cms.InputTag('hltHpsPFTauPairsMatchedToL1PassingPtGt125'),
      minNumber = cms.int32(1),
      maxNumber = cms.int32(-1)
    )
    process.debugSequence += process.hltHpsPFTauPairsMatchedToL1PassingPtGt125Filter

    process.dumpOnlineHpsPFTaus = cms.EDAnalyzer("DumpRecoPFTaus",
      src = cms.InputTag(src_hltPFTaus),
      src_sumChargedIso = cms.InputTag(src_hltPFTauChargedIsoPtSum),
      src_discriminators = cms.VInputTag()
    )
    process.debugSequence += process.dumpOnlineHpsPFTaus

    process.dumpOnlineHpsPFTausMatchedToL1 = cms.EDAnalyzer("DumpRecoPFTaus",
      src = cms.InputTag(src_hltPFTaus_matchedToL1),
      src_sumChargedIso = cms.InputTag(src_hltPFTauChargedIsoPtSum_matchedToL1),
      src_discriminators = cms.VInputTag()
    )
    process.debugSequence += process.dumpOnlineHpsPFTausMatchedToL1

    process.runLumiSectionEventNumberAnalyzer = cms.EDAnalyzer("RunLumiSectionEventNumberAnalyzer",
      output = cms.string(outputFileName),
      separator = cms.string(":")
    )
    process.debugSequence += process.runLumiSectionEventNumberAnalyzer
