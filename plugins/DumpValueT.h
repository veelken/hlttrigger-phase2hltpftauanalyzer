#ifndef HLTrigger_TallinnHLTPFTauAnalyzer_DumpValueT_h
#define HLTrigger_TallinnHLTPFTauAnalyzer_DumpValueT_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"

#include <string>
#include <iostream>
#include <iomanip>

template <typename T>
class DumpValueT : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit DumpValueT(const edm::ParameterSet& cfg)
    : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
  {
    src_ = cfg.getParameter<edm::InputTag>("src");
    token_ = consumes<T>(src_);
  }

  // destructor
  ~DumpValueT() {}
    
 private:
  void analyze(const edm::Event& evt, const edm::EventSetup& es)
  {
    std::cout << "<DumpValueT::analyze (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
    std::cout << " src = " << src_ << std::endl;

    edm::Handle<T> value;
    evt.getByToken(token_, value);
  
    std::cout << "value = " << (*value) << std::endl;
  }

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<T> token_;
};

#endif   
