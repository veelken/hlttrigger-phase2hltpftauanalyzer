#include "HLTrigger/TallinnHLTPFTauAnalyzer/plugins/PATTauAnalyzerSignal.h"

#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

PATTauAnalyzerSignal::PATTauAnalyzerSignal(const edm::ParameterSet& cfg)
  : BaseTauAnalyzerSignal(cfg)
{
  srcPFTaus_ = cfg.getParameter<edm::InputTag>("srcPFTaus");
  tokenPFTaus_ = consumes<pat::TauCollection>(srcPFTaus_);
  pfTauDiscriminator_ = cfg.getParameter<std::string>("pfTauDiscriminator");
}

PATTauAnalyzerSignal::~PATTauAnalyzerSignal()
{}

std::vector<BaseTau> 
PATTauAnalyzerSignal::buildBaseTaus(const edm::Event& evt, const edm::EventSetup& es)
{
  std::vector<BaseTau> baseTaus;

  edm::Handle<pat::TauCollection> pfTaus;
  evt.getByToken(tokenPFTaus_, pfTaus);
  
  size_t numPFTaus = pfTaus->size();
  for ( size_t idxPFTau = 0; idxPFTau < numPFTaus; ++idxPFTau ) 
  { 
    const pat::Tau& pfTau = pfTaus->at(idxPFTau);
    if ( pfTau.leadChargedHadrCand().isNonnull() )
    {
      double discriminator = pfTau.tauID(pfTauDiscriminator_);
      double zVtx = pfTau.leadChargedHadrCand()->vertex().z();
      baseTaus.push_back(BaseTau(pfTau.p4(), pfTau.leadChargedHadrCand()->p4(), discriminator, zVtx));
    }
  }

  return baseTaus;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(PATTauAnalyzerSignal);
