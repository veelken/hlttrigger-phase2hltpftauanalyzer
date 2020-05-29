#ifndef HLTTrigger_TallinnHLTPFTauAnalyzer_DumpFloat_h
#define HLTTrigger_TallinnHLTPFTauAnalyzer_DumpFloat_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <string>

class DumpFloat : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit DumpFloat(const edm::ParameterSet&);
    
  // destructor
  ~DumpFloat();
    
 private:
  void analyze(const edm::Event&, const edm::EventSetup&);

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<float> token_;
};

#endif   
