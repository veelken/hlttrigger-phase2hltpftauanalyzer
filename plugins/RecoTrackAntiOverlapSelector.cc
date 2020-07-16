#include "HLTrigger/TallinnHLTPFTauAnalyzer/plugins/RecoTrackAntiOverlapSelector.h"

RecoTrackAntiOverlapSelector::RecoTrackAntiOverlapSelector(const edm::ParameterSet& cfg) 
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<reco::TrackCollection>(src_);
  srcNotToBeFiltered_ = cfg.getParameter<vInputTag>("srcNotToBeFiltered");
  for ( auto it : srcNotToBeFiltered_ )
  {
    tokensNotToBeFiltered_.push_back(consumes<reco::CandidateView>(it));
  }
  dRmin_ = cfg.getParameter<double>("dRmin");
  invert_ = ( cfg.exists("invert") ) ?
    cfg.getParameter<bool>("invert") : false;

  produces<reco::TrackCollection>();
}

RecoTrackAntiOverlapSelector::~RecoTrackAntiOverlapSelector()
{}

void RecoTrackAntiOverlapSelector::produce(edm::Event& evt, const edm::EventSetup& es)
{
  std::unique_ptr<reco::TrackCollection> tracks_selected(new reco::TrackCollection());

  edm::Handle<reco::TrackCollection> tracks;
  evt.getByToken(token_, tracks);

  std::vector<bool> isOverlap(tracks->size());
    
  for ( auto it : tokensNotToBeFiltered_ )
  {
    edm::Handle<reco::CandidateView> particlesNotToBeFiltered;
    evt.getByToken(it, particlesNotToBeFiltered);
      
    for ( reco::CandidateView::const_iterator particleNotToBeFiltered = particlesNotToBeFiltered->begin(); 
	  particleNotToBeFiltered != particlesNotToBeFiltered->end();  ++particleNotToBeFiltered ) 
    {
      size_t numTracks = tracks->size();
      for ( size_t idxTrack = 0; idxTrack < numTracks; ++idxTrack )
      {
	const reco::Track& track = tracks->at(idxTrack);
	double dR = reco::deltaR(track.eta(), track.phi(), particleNotToBeFiltered->eta(), particleNotToBeFiltered->phi());	  
	if ( dR < dRmin_ ) isOverlap[idxTrack] = true;
      }
    }
  }
    
  size_t numTracks = tracks->size();
  for ( size_t idxTrack = 0; idxTrack < numTracks; ++idxTrack )
  {
    const reco::Track& track = tracks->at(idxTrack);
    if ( (invert_ == false && !isOverlap[idxTrack]) ||
	 (invert_ == true  &&  isOverlap[idxTrack]) ) 
    {
      tracks_selected->push_back(track); 
    }
  }

  evt.put(std::move(tracks_selected));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RecoTrackAntiOverlapSelector);
