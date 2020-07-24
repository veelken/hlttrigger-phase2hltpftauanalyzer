#include "HLTrigger/TallinnHLTPFTauAnalyzer/plugins/MyCandViewCountFilter.h"

MyCandViewCountFilter::MyCandViewCountFilter(const edm::ParameterSet& cfg)
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<reco::CandidateView>(src_);

  minNumber_ = ( cfg.exists("minNumber") ) ? cfg.getParameter<int>("minNumber") : -1;
  maxNumber_ = ( cfg.exists("maxNumber") ) ? cfg.getParameter<int>("maxNumber") : -1;
}

MyCandViewCountFilter::~MyCandViewCountFilter()
{}

bool MyCandViewCountFilter::filter(edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<reco::CandidateView> particles;
  evt.getByToken(token_, particles);

  int numParticles = particles->size();
  if ( (minNumber_ < 0 || numParticles >= minNumber_) &&
       (maxNumber_ < 0 || numParticles <= maxNumber_) )
  {
    return true;
  }
  else
  {
    return false;
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(MyCandViewCountFilter);
