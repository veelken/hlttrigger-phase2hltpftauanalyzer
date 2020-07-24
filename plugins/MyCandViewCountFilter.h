#ifndef HLTrigger_TallinnHLTPFTauAnalyzer_MyCandViewCountFilter_h
#define HLTrigger_TallinnHLTPFTauAnalyzer_MyCandViewCountFilter_h

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"    // reco::Candidate
#include "DataFormats/Candidate/interface/CandidateFwd.h" // reco::CandidateView

class MyCandViewCountFilter : public edm::EDFilter
{
 public:
  explicit MyCandViewCountFilter(const edm::ParameterSet&);
  virtual ~MyCandViewCountFilter();
    
 private:
  bool filter(edm::Event& evt, const edm::EventSetup& es);

  edm::InputTag src_;
  edm::EDGetTokenT<reco::CandidateView> token_;

  int minNumber_;
  int maxNumber_;
};

#endif   
