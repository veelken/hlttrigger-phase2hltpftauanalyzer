#ifndef HLTrigger_TallinnHLTPFTauAnalyzer_PATTauAnalyzerBackground_h
#define HLTrigger_TallinnHLTPFTauAnalyzer_PATTauAnalyzerBackground_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/PatCandidates/interface/Tau.h"                               // pat::Tau, pat::TauCollection
#include "HLTrigger/TallinnHLTPFTauAnalyzer/interface/BaseTauAnalyzerBackground.h" // BaseTauAnalyzerBackground
#include "HLTrigger/TallinnHLTPFTauAnalyzer/interface/histogramAuxFunctions.h"     // fillWithOverFlow()

class PATTauAnalyzerBackground : public BaseTauAnalyzerBackground
{
 public:
  // constructor 
  explicit PATTauAnalyzerBackground(const edm::ParameterSet&);
    
  // destructor
  ~PATTauAnalyzerBackground();
    
 private:
  std::vector<BaseTau> buildBaseTaus(const edm::Event&, const edm::EventSetup&);

  edm::InputTag srcPFTaus_;
  edm::EDGetTokenT<pat::TauCollection> tokenPFTaus_;
  std::string pfTauDiscriminator_;
};

#endif   
