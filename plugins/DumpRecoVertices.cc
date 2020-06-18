#include "HLTTrigger/TallinnHLTPFTauAnalyzer/plugins/DumpRecoVertices.h"

#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

DumpRecoVertices::DumpRecoVertices(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<reco::VertexCollection>(src_);
}

DumpRecoVertices::~DumpRecoVertices()
{}

namespace
{
  double compSumTrackPt2(const reco::Vertex& vertex)
  {
    double sumTrackPt2 = 0.;
    for ( reco::Vertex::trackRef_iterator track = vertex.tracks_begin(); track != vertex.tracks_end(); ++track )
    {
      double track_pt = (*track)->pt();
      sumTrackPt2 += (track_pt*track_pt);
    }
    return sumTrackPt2;
  }
}

void DumpRecoVertices::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  std::cout << "<DumpRecoVertices::analyze (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
  std::cout << " src = " << src_ << std::endl;

  edm::Handle<reco::VertexCollection> vertices;
  evt.getByToken(token_, vertices);
  
  size_t numVertices = vertices->size();
  for ( size_t idxVertex = 0; idxVertex < numVertices; ++idxVertex ) 
  {
    const reco::Vertex& vertex = vertices->at(idxVertex);
/*
    std::cout << "offlineVertex #" << idxVertex << ":" 
	      << " z0 = " << vertex.position().z() << " (sum track pT^2 = " << compSumTrackPt2(vertex) << ")" << std::endl;
    size_t idxTrack = 0;
    for ( reco::Vertex::trackRef_iterator track = vertex.tracks_begin(); track != vertex.tracks_end(); ++track )
    {
      std::cout << " offlineTrack #" << idxTrack << ":" 
		<< " pT = " << (*track)->pt() << ", eta = " << (*track)->eta() << ", phi = " << (*track)->phi() << std::endl;
      ++idxTrack;
    }
 */
    std::cout << "offlineVertex #" << idxVertex << ":" 
	      << " z0 = " << vertex.position().z() << std::endl;
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DumpRecoVertices);





