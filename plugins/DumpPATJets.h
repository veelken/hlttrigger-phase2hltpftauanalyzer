#ifndef HLTTrigger_TallinnHLTPFTauAnalyzer_DumpPATJets_h
#define HLTTrigger_TallinnHLTPFTauAnalyzer_DumpPATJets_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include <vector>
#include <string>

class DumpPATJets : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit DumpPATJets(const edm::ParameterSet&);
    
  // destructor
  ~DumpPATJets();
    
 private:
  void analyze(const edm::Event&, const edm::EventSetup&);

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<pat::JetCollection> token_;
};

#endif   
