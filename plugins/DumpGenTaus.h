#ifndef HLTTrigger_TallinnHLTPFTauAnalyzer_DumpGenTaus_h
#define HLTTrigger_TallinnHLTPFTauAnalyzer_DumpGenTaus_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include <vector>
#include <string>

class DumpGenTaus : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit DumpGenTaus(const edm::ParameterSet&);
    
  // destructor
  ~DumpGenTaus();
    
 private:
  void analyze(const edm::Event&, const edm::EventSetup&);

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<reco::GenJetCollection> token_;
};

#endif   
