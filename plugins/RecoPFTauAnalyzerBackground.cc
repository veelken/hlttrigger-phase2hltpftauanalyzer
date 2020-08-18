#include "HLTrigger/TallinnHLTPFTauAnalyzer/plugins/RecoPFTauAnalyzerBackground.h"

#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

RecoPFTauAnalyzerBackground::RecoPFTauAnalyzerBackground(const edm::ParameterSet& cfg)
  : BaseTauAnalyzerBackground(cfg)
{
  srcPFTaus_ = cfg.getParameter<edm::InputTag>("srcPFTaus");
  tokenPFTaus_ = consumes<reco::PFTauCollection>(srcPFTaus_);
  srcPFTauDiscriminator_ = cfg.getParameter<edm::InputTag>("srcPFTauDiscriminator");
  tokenPFTauDiscriminator_ = consumes<reco::PFTauDiscriminator>(srcPFTauDiscriminator_);
}

RecoPFTauAnalyzerBackground::~RecoPFTauAnalyzerBackground()
{}

std::vector<BaseTau> 
RecoPFTauAnalyzerBackground::buildBaseTaus(const edm::Event& evt, const edm::EventSetup& es)
{
  std::vector<BaseTau> baseTaus;

  edm::Handle<reco::PFTauCollection> pfTaus;
  evt.getByToken(tokenPFTaus_, pfTaus);
  edm::Handle<reco::PFTauDiscriminator> pfTauDiscriminator;
  evt.getByToken(tokenPFTauDiscriminator_, pfTauDiscriminator);

  size_t numPFTaus = pfTaus->size();
  for ( size_t idxPFTau = 0; idxPFTau < numPFTaus; ++idxPFTau ) 
  { 
    reco::PFTauRef pfTauRef(pfTaus, idxPFTau);
    if ( (min_pt_ < 0. || pfTauRef->pt()                                 >= min_pt_ ) &&
         (max_pt_ < 0. || pfTauRef->pt()                                 <= max_pt_ ) &&
         (                pfTauRef->leadPFChargedHadrCand().isNonnull()             && 
                          pfTauRef->leadPFChargedHadrCand()->bestTrack()            ) )
    {
      double discriminator = (*pfTauDiscriminator)[pfTauRef];
      double zVtx = pfTauRef->leadPFChargedHadrCand()->bestTrack()->vertex().z();
      baseTaus.push_back(BaseTau(pfTauRef->p4(), pfTauRef->leadPFChargedHadrCand()->p4(), discriminator, zVtx));
    }
  }
  
  return baseTaus;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(RecoPFTauAnalyzerBackground);
