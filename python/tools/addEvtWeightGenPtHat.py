import FWCore.ParameterSet.Config as cms

from HLTrigger.TallinnHLTPFTauAnalyzer.samples_cfi import background_samples

def addEvtWeightGenPtHat(process, hlt_srcVertices):

    process.stitchingWeight = cms.EDProducer("EvtWeightProducerGenPtHatStitching",
      src_genEventInfo = cms.InputTag('generator'),
      src_pileupSummaryInfo = cms.InputTag('addPileupInfo'),
      samples = cms.PSet(),
      pT_hat_bins = cms.vdouble(0., 30., 50., 80., 120., 170., 300.),
      bxFrequency = cms.double(2.8e+7) # 28 MHz bunch-crossing frequency
    )

    for sampleName, sample in background_samples.items(): 
        if sampleName in [ "minbias", "qcd_pt30to50", "qcd_pt50to80", "qcd_pt80to120", "qcd_pt120to170", "qcd_pt170to300" ]:
            setattr(process.stitchingWeight.samples, sampleName, cms.PSet(
              crossSection = cms.double(sample['crossSection']),
              numEvents = cms.int32(sample['samples'][hlt_srcVertices]['numEvents'])
            ))
