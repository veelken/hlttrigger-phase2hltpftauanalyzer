#ifndef HLTrigger_TallinnHLTPFTauAnalyzer_EvtWeightProducerLumiScale_h
#define HLTrigger_TallinnHLTPFTauAnalyzer_EvtWeightProducerLumiScale_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

class EvtWeightProducerLumiScale : public edm::EDProducer 
{
 public:
  explicit EvtWeightProducerLumiScale(const edm::ParameterSet& cfg);
  ~EvtWeightProducerLumiScale();

 private:
  void produce(edm::Event& evt, const edm::EventSetup& es);

  double lumiScale_;
};

#endif
