#ifndef HLTrigger_TallinnHLTPFTauAnalyzer_CandidateViewAntiOverlapSelector_h
#define HLTrigger_TallinnHLTPFTauAnalyzer_CandidateViewAntiOverlapSelector_h

#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h" 

#include "DataFormats/Math/interface/deltaR.h"

#include <vector>

class CandidateViewAntiOverlapSelector : public edm::stream::EDFilter<>
{
 public:
  explicit CandidateViewAntiOverlapSelector(const edm::ParameterSet& cfg);
  ~CandidateViewAntiOverlapSelector();

  bool filter(edm::Event& evt, const edm::EventSetup& es);

  static void fillDescriptions(edm::ConfigurationDescriptions& desc);

 private:
  edm::InputTag src_;
  edm::EDGetTokenT<edm::View<reco::Candidate>> token_;
  
  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag srcNotToBeFiltered_;
  typedef std::vector<edm::EDGetTokenT<reco::CandidateView>> vCandidateViewToken;
  vCandidateViewToken tokensNotToBeFiltered_;

  //when invert is TRUE the selector looks for overlapping objects
  bool invert_;

  double dRmin_;

  bool filter_;
};

#endif
