#include "HLTrigger/TallinnHLTPFTauAnalyzer/plugins/EvtWeightProducerLumiScale.h"

EvtWeightProducerLumiScale::EvtWeightProducerLumiScale(const edm::ParameterSet& cfg) 
{
  lumiScale_ = cfg.getParameter<double>("lumiScale");

  produces<double>();
}

EvtWeightProducerLumiScale::~EvtWeightProducerLumiScale()
{}

void EvtWeightProducerLumiScale::produce(edm::Event& evt, const edm::EventSetup& es)
{
  std::unique_ptr<double> evtWeight(new double());

  *evtWeight = lumiScale_;

  evt.put(std::move(evtWeight));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(EvtWeightProducerLumiScale);
