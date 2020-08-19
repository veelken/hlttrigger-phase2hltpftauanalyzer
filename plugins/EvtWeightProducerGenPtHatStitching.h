#ifndef HLTrigger_TallinnHLTPFTauAnalyzer_EvtWeightProducerGenPtHatStitching_h
#define HLTrigger_TallinnHLTPFTauAnalyzer_EvtWeightProducerGenPtHatStitching_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include <TH1.h>

class EvtWeightProducerGenPtHatStitching : public edm::EDProducer 
{
 public:
  explicit EvtWeightProducerGenPtHatStitching(const edm::ParameterSet& cfg);
  ~EvtWeightProducerGenPtHatStitching();

 private:
  void produce(edm::Event& evt, const edm::EventSetup& es);

  std::string moduleLabel_;

  edm::InputTag src_genEventInfo_;
  edm::EDGetTokenT<GenEventInfoProduct> token_genEventInfo_;

  edm::InputTag src_pileupSummaryInfo_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> token_pileupSummaryInfo_;

  struct sampleEntryType
  {
    std::string name_;
    double crossSection_;
    int numEvents_;
    int pT_hat_bin_;
  };
  std::vector<sampleEntryType> samples_;
  sampleEntryType sample_minbias_;
  std::vector<sampleEntryType> samples_qcd_;

  std::vector<double> p_k_;

  std::vector<double> binning_pT_hat_;
  int numBins_pT_hat_;

  double bxFrequency_;

  TH1* histogram_X_;
};

#endif
