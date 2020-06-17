#ifndef HLTTrigger_TallinnHLTPFTauAnalyzer_DumpRecoPFTauPairs_h
#define HLTTrigger_TallinnHLTPFTauAnalyzer_DumpRecoPFTauPairs_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Phase2HLTPFTaus/interface/PFTauPair.h"    // reco::PFTauPair
#include "DataFormats/Phase2HLTPFTaus/interface/PFTauPairFwd.h" // reco::PFTauPairCollection

#include <vector>
#include <string>

class DumpRecoPFTauPairs : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit DumpRecoPFTauPairs(const edm::ParameterSet&);
    
  // destructor
  ~DumpRecoPFTauPairs();
    
 private:
  void analyze(const edm::Event&, const edm::EventSetup&);

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<reco::PFTauPairCollection> token_;
};

#endif   
