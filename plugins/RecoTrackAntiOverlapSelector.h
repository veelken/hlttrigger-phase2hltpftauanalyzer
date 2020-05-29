#ifndef HLTTrigger_TallinnHLTPFTauAnalyzer_RecoTrackAntiOverlapSelector_h
#define HLTTrigger_TallinnHLTPFTauAnalyzer_RecoTrackAntiOverlapSelector_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h" // reco::CandidateView
#include "DataFormats/Candidate/interface/Candidate.h"    // reco::Candidate
#include "DataFormats/TrackReco/interface/Track.h"        // reco::Track
#include "DataFormats/TrackReco/interface/TrackFwd.h"     // reco::TrackCollection

#include "DataFormats/Math/interface/deltaR.h"

#include <vector>

class RecoTrackAntiOverlapSelector : public edm::EDProducer 
{
 public:
  explicit RecoTrackAntiOverlapSelector(const edm::ParameterSet& cfg);
  ~RecoTrackAntiOverlapSelector();

 private:
  void produce(edm::Event& evt, const edm::EventSetup& es);

  edm::InputTag src_;
  edm::EDGetTokenT<reco::TrackCollection> token_;

  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag srcNotToBeFiltered_;
  typedef std::vector<edm::EDGetTokenT<reco::CandidateView>> vCandidateViewToken;
  vCandidateViewToken tokensNotToBeFiltered_;

  //when invert is TRUE the selector looks for overlapping objects
  bool invert_;

  double dRmin_;
};

#endif
