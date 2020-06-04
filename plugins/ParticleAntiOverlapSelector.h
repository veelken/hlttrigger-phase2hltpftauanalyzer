#ifndef HLTTrigger_TallinnHLTPFTauAnalyzer_ParticleAntiOverlapSelector_h
#define HLTTrigger_TallinnHLTPFTauAnalyzer_ParticleAntiOverlapSelector_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h" 
#include "DataFormats/Candidate/interface/Candidate.h" 

#include "DataFormats/Math/interface/deltaR.h"

#include <vector>

template <class T, class TCollection = std::vector<T>>
class ParticleAntiOverlapSelector
{
 public:
  typedef TCollection collection;

  explicit ParticleAntiOverlapSelector(const edm::ParameterSet& cfg, edm::ConsumesCollector&& iC)
  {
    srcNotToBeFiltered_ = cfg.getParameter<vInputTag>("srcNotToBeFiltered");
    for ( auto it : srcNotToBeFiltered_ )
    {
      tokensNotToBeFiltered_.push_back(iC.consumes<reco::CandidateView>(it));
    }
    dRmin_ = cfg.getParameter<double>("dRmin");
    invert_ = ( cfg.exists("invert") ) ?
      cfg.getParameter<bool>("invert") : false;
  }

  typename std::vector<const T*>::const_iterator begin() const { return selected_.begin(); }
  typename std::vector<const T*>::const_iterator end() const { return selected_.end(); }

  void select(const edm::Handle<TCollection>& particlesToBeFiltered, const edm::Event& evt, const edm::EventSetup& es)
  {
    selected_.clear();
    //std::cout << "#particlesToBeFiltered = " << particlesToBeFiltered->size() << std::endl;

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
	  const T& particleToBeFiltered = particlesToBeFiltered->at(idxParticleToBeFiltered);
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
      const T& particleToBeFiltered = particlesToBeFiltered->at(idxParticleToBeFiltered);
      if ( (invert_ == false && !isOverlap[idxParticleToBeFiltered]) ||
           (invert_ == true  &&  isOverlap[idxParticleToBeFiltered]) ) 
      {
	selected_.push_back(&particleToBeFiltered); 
      }
    }
  }

  size_t size() const { return selected_.size(); }

 private:
  std::vector<const T*> selected_;

  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag srcNotToBeFiltered_;
  typedef std::vector<edm::EDGetTokenT<reco::CandidateView>> vCandidateViewToken;
  vCandidateViewToken tokensNotToBeFiltered_;

  //when invert is TRUE the selector looks for overlapping objects
  bool invert_;

  double dRmin_;
};

#endif
