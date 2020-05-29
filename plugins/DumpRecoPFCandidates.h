#ifndef HLTTrigger_TallinnHLTPFTauAnalyzer_DumpRecoPFCandidates_h
#define HLTTrigger_TallinnHLTPFTauAnalyzer_DumpRecoPFCandidates_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"    // reco::PFCandidate
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h" // reco::PFCandidateCollection

#include <vector>
#include <string>

class DumpRecoPFCandidates : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit DumpRecoPFCandidates(const edm::ParameterSet&);
    
  // destructor
  ~DumpRecoPFCandidates();
    
 private:
  void analyze(const edm::Event&, const edm::EventSetup&);

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<reco::PFCandidateCollection> token_;

  double minPt_;
  double maxAbsEta_;
};

#endif   
