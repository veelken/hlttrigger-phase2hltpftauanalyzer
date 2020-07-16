#ifndef HLTrigger_TallinnHLTPFTauAnalyzer_DumpPackedCandidates_h
#define HLTrigger_TallinnHLTPFTauAnalyzer_DumpPackedCandidates_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // pat::PackedCandidate, pat::PackedCandidateCollection

#include <vector>
#include <string>

class DumpPackedCandidates : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit DumpPackedCandidates(const edm::ParameterSet&);
    
  // destructor
  ~DumpPackedCandidates();
    
 private:
  void analyze(const edm::Event&, const edm::EventSetup&);

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> token_;

  double minPt_;
};

#endif   
