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
            kFactor_pT_hat = None
            pT_hat_bin = None            
            if sampleName == "minbias":
                kFactor_pT_hat = 1.
                pT_hat_bin = -1
            elif sampleName == "qcd_pt30to50":
                kFactor_pT_hat = sample['kFactor_pT_hat']
                pT_hat_bin =  1
            elif sampleName == "qcd_pt50to80":
                kFactor_pT_hat = sample['kFactor_pT_hat']
                pT_hat_bin =  2
            elif sampleName == "qcd_pt80to120":
                kFactor_pT_hat = sample['kFactor_pT_hat']
                pT_hat_bin =  3
            elif sampleName == "qcd_pt120to170":
                kFactor_pT_hat = sample['kFactor_pT_hat']
                pT_hat_bin =  4
            elif sampleName == "qcd_pt170to300":
                kFactor_pT_hat = sample['kFactor_pT_hat']
                pT_hat_bin =  5
            else:
                raise ValueError("Invalid sample = '%s' !!" % sampleName)
            ##print("Using pT_hat k-factor of %1.2f for sample = '%s'." % (kFactor_pT_hat, sampleName))
            setattr(process.stitchingWeight.samples, sampleName, cms.PSet(
              crossSection = cms.double(sample['crossSection']*kFactor_pT_hat),
              numEvents = cms.uint32(sample['samples'][hlt_srcVertices]['numEvents']),
              pT_hat_bin = cms.int32(pT_hat_bin)
            ))
