#include "HLTrigger/TallinnHLTPFTauAnalyzer/plugins/RecoPFTauAnalyzerSignal.h"

#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

RecoPFTauAnalyzerSignal::RecoPFTauAnalyzerSignal(const edm::ParameterSet& cfg)
  : BaseTauAnalyzerSignal(cfg)
{
  srcPFTaus_ = cfg.getParameter<edm::InputTag>("srcPFTaus");
  tokenPFTaus_ = consumes<reco::PFTauCollection>(srcPFTaus_);
  srcPFTauDiscriminator_ = cfg.getParameter<edm::InputTag>("srcPFTauDiscriminator");
  tokenPFTauDiscriminator_ = consumes<reco::PFTauDiscriminator>(srcPFTauDiscriminator_);
}

RecoPFTauAnalyzerSignal::~RecoPFTauAnalyzerSignal()
{}

std::vector<BaseTau> 
RecoPFTauAnalyzerSignal::buildBaseTaus(const edm::Event& evt, const edm::EventSetup& es)
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
    if ( pfTauRef->leadPFChargedHadrCand().isNonnull() && pfTauRef->leadPFChargedHadrCand()->bestTrack() )
    {
      double discriminator = (*pfTauDiscriminator)[pfTauRef];
      double zVtx = pfTauRef->leadPFChargedHadrCand()->bestTrack()->vertex().z();
      baseTaus.push_back(BaseTau(pfTauRef->p4(), pfTauRef->leadPFChargedHadrCand()->p4(), discriminator, zVtx));
    }
  }
  
  return baseTaus;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(RecoPFTauAnalyzerSignal);
