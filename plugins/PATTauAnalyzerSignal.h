#ifndef HLTrigger_TallinnHLTPFTauAnalyzer_PATTauAnalyzerSignal_h
#define HLTrigger_TallinnHLTPFTauAnalyzer_PATTauAnalyzerSignal_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/Tau.h"                           // pat::Tau, pat::TauCollection
#include "HLTrigger/TallinnHLTPFTauAnalyzer/interface/BaseTauAnalyzerSignal.h" // BaseTauAnalyzerSignal

#include <vector>    // std::vector
#include <string>    // std::string

using namespace dqm::implementation;

class PATTauAnalyzerSignal : public BaseTauAnalyzerSignal
{
 public:
  // constructor 
  explicit PATTauAnalyzerSignal(const edm::ParameterSet&);
    
  // destructor
  ~PATTauAnalyzerSignal();
    
 private:
  std::vector<BaseTau> buildBaseTaus(const edm::Event&, const edm::EventSetup&);

  edm::InputTag srcPFTaus_;
  edm::EDGetTokenT<pat::TauCollection> tokenPFTaus_;
  std::string pfTauDiscriminator_;
};

#endif   
