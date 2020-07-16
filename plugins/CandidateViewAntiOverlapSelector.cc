#include "HLTrigger/TallinnHLTPFTauAnalyzer/plugins/CandidateViewAntiOverlapSelector.h"

#include <FWCore/ParameterSet/interface/ConfigurationDescriptions.h>
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/Common/interface/Handle.h"

CandidateViewAntiOverlapSelector::CandidateViewAntiOverlapSelector(const edm::ParameterSet& cfg)
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<edm::View<reco::Candidate>>(src_);
  srcNotToBeFiltered_ = cfg.getParameter<vInputTag>("srcNotToBeFiltered");
  for ( auto it : srcNotToBeFiltered_ )
  {
    tokensNotToBeFiltered_.push_back(consumes<reco::CandidateView>(it));
  }
  dRmin_ = cfg.getParameter<double>("dRmin");
  invert_ = ( cfg.exists("invert") ) ?
    cfg.getParameter<bool>("invert") : false;
  filter_ = cfg.getParameter<bool>("filter");
  produces<reco::CandidatePtrVector>(); 
}

CandidateViewAntiOverlapSelector::~CandidateViewAntiOverlapSelector()
{}

bool 
CandidateViewAntiOverlapSelector::filter(edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<edm::View<reco::Candidate>> particlesToBeFiltered;
  evt.getByToken(token_, particlesToBeFiltered);

  auto selected = std::make_unique<reco::CandidatePtrVector>();

  std::vector<bool> isOverlap(particlesToBeFiltered->size());

  for ( auto it : tokensNotToBeFiltered_ )
  {
    edm::Handle<reco::CandidateView> particlesNotToBeFiltered;
    evt.getByToken(it, particlesNotToBeFiltered);
    //std::cout << "#particlesNotToBeFiltered = " << particlesNotToBeFiltered->size() << std::endl;
      
    for ( reco::CandidateView::const_iterator particleNotToBeFiltered = particlesNotToBeFiltered->begin(); 
	  particleNotToBeFiltered != particlesNotToBeFiltered->end();  ++particleNotToBeFiltered ) 
    {
      size_t numParticlesToBeFiltered = particlesToBeFiltered->size();
      for ( size_t idxParticleToBeFiltered = 0; idxParticleToBeFiltered < numParticlesToBeFiltered; ++idxParticleToBeFiltered )
      {
        const reco::Candidate& particleToBeFiltered = particlesToBeFiltered->at(idxParticleToBeFiltered);
        double dR = reco::deltaR(particleToBeFiltered.eta(), particleToBeFiltered.phi(), particleNotToBeFiltered->eta(), particleNotToBeFiltered->phi());	  
        //std::cout << "dR = " << dR << std::endl;	
        if ( dR < dRmin_ ) isOverlap[idxParticleToBeFiltered] = true;
      }
    }
  }

  //for ( size_t idx = 0; idx < isOverlap.size(); ++idx )
  //{
  //  std::cout << "isOverlap[" << idx << "] = " << isOverlap[idx] << std::endl;
  //}
    
  size_t numParticlesToBeFiltered = particlesToBeFiltered->size();
  for ( size_t idxParticleToBeFiltered = 0; idxParticleToBeFiltered < numParticlesToBeFiltered; ++idxParticleToBeFiltered )
  {
    if ( (invert_ == false && !isOverlap[idxParticleToBeFiltered]) ||
         (invert_ == true  &&  isOverlap[idxParticleToBeFiltered]) ) 
    {
      selected->push_back(particlesToBeFiltered->ptrAt(idxParticleToBeFiltered));
    }
  }

  bool retVal = ( filter_ && selected->size() == 0 ) ? false : true;

  evt.put(std::move(selected));

  return retVal;
}

void 
CandidateViewAntiOverlapSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src")->setComment("collection of particles to be filtered");
  desc.add<vInputTag>("srcNotToBeFiltered")->setComment("collections of particles NOT to be filtered");
  desc.add<double>("dRmin")->setComment("size of dR cone used to match particles to be filtered to particles NOT to be filtered");
  desc.add<bool>("invert", false)->setComment("to be set to TRUE (FALSE) when you want to select overlapping (non-overlapping) objects");
  desc.add<bool>("filter")->setComment("to be set to TRUE (FALSE) when you want to abort the event processing in case no particle from the collection specified by the 'src' parameter passes the selection");
  descriptions.add("CandidateViewAntiOverlapSelector", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(CandidateViewAntiOverlapSelector);
