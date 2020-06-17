#include "HLTTrigger/TallinnHLTPFTauAnalyzer/plugins/DumpRecoPFTauPairs.h"

#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

DumpRecoPFTauPairs::DumpRecoPFTauPairs(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<reco::PFTauPairCollection>(src_);
}

DumpRecoPFTauPairs::~DumpRecoPFTauPairs()
{}

void DumpRecoPFTauPairs::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  std::cout << "<DumpRecoPFTauPairs::analyze (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
  std::cout << " src = " << src_ << std::endl;

  edm::Handle<reco::PFTauPairCollection> tauPairs;
  evt.getByToken(token_, tauPairs);
  
  size_t numTauPairs = tauPairs->size();
  for ( size_t idxTauPair = 0; idxTauPair < numTauPairs; ++idxTauPair ) 
  {
    const reco::PFTauPair& tauPair = tauPairs->at(idxTauPair);
    std::cout << "PFTauPair #" << idxTauPair << ": " << tauPair << std::endl;
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DumpRecoPFTauPairs);





