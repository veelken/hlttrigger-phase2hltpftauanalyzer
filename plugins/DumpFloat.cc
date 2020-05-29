#include "HLTTrigger/TallinnHLTPFTauAnalyzer/plugins/DumpFloat.h"

#include "DataFormats/Common/interface/Handle.h"

#include <iostream>
#include <iomanip>

DumpFloat::DumpFloat(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<float>(src_);
}

DumpFloat::~DumpFloat()
{}

void DumpFloat::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  std::cout << "<DumpFloat::analyze (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
  std::cout << " src = " << src_ << std::endl;

  edm::Handle<float> value;
  evt.getByToken(token_, value);
  
  std::cout << "value = " << (*value) << std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DumpFloat);





