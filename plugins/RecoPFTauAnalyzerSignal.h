#ifndef HLTrigger_TallinnHLTPFTauAnalyzer_RecoPFTauAnalyzerSignal_h
#define HLTrigger_TallinnHLTPFTauAnalyzer_RecoPFTauAnalyzerSignal_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TauReco/interface/PFTau.h"                               // reco::PFTau
#include "DataFormats/TauReco/interface/PFTauFwd.h"                            // reco::PFTauCollection
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"                  // reco::PFTauDiscriminator
#include "HLTrigger/TallinnHLTPFTauAnalyzer/interface/BaseTauAnalyzerSignal.h" // BaseTauAnalyzerSignal

#include <vector>    // std::vector
#include <string>    // std::string

using namespace dqm::implementation;

class RecoPFTauAnalyzerSignal : public BaseTauAnalyzerSignal
{
 public:
  // constructor 
  explicit RecoPFTauAnalyzerSignal(const edm::ParameterSet&);
    
  // destructor
  ~RecoPFTauAnalyzerSignal();
    
 private:
  std::vector<BaseTau> buildBaseTaus(const edm::Event&, const edm::EventSetup&);

  edm::InputTag srcPFTaus_;
  edm::EDGetTokenT<reco::PFTauCollection> tokenPFTaus_;
  edm::InputTag srcPFTauDiscriminator_;
  edm::EDGetTokenT<reco::PFTauDiscriminator> tokenPFTauDiscriminator_;
};

#endif   
