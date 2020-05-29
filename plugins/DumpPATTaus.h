#ifndef HLTTrigger_TallinnHLTPFTauAnalyzer_DumpPATTaus_h
#define HLTTrigger_TallinnHLTPFTauAnalyzer_DumpPATTaus_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include <vector>
#include <string>

class DumpPATTaus : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit DumpPATTaus(const edm::ParameterSet&);
    
  // destructor
  ~DumpPATTaus();
    
 private:
  void analyze(const edm::Event&, const edm::EventSetup&);

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<pat::TauCollection> token_;
};

#endif   
