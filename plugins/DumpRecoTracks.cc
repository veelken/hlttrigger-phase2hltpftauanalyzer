#include "HLTrigger/TallinnHLTPFTauAnalyzer/plugins/DumpRecoTracks.h"

#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

DumpRecoTracks::DumpRecoTracks(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<reco::TrackCollection>(src_);
}

DumpRecoTracks::~DumpRecoTracks()
{}

void DumpRecoTracks::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  std::cout << "<DumpRecoTracks::analyze (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
  std::cout << " src = " << src_ << std::endl;

  edm::Handle<reco::TrackCollection> tracks;
  evt.getByToken(token_, tracks);
  
  size_t numTracks = tracks->size();
  for ( size_t idxTrack = 0; idxTrack < numTracks; ++idxTrack ) 
  {
    const reco::Track& track = tracks->at(idxTrack);
    std::cout << "Track #" << idxTrack << ":" 
	      << " pT = " << track.pt() << ", eta = " << track.eta() << ", phi = " << track.phi() << std::endl;
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DumpRecoTracks);





