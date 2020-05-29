#include "HLTTrigger/TallinnHLTPFTauAnalyzer/plugins/DumpGenTaus.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h" // JetMCTagUtils::genTauDecayMode()

#include <TMath.h>

#include <iostream>
#include <iomanip>
#include <algorithm> // std::sort

DumpGenTaus::DumpGenTaus(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<reco::GenJetCollection>(src_);
}

DumpGenTaus::~DumpGenTaus()
{}

namespace
{
  bool
  isHigherPt(const reco::GenJet& genTau1,
	     const reco::GenJet& genTau2)
  {
    return genTau1.pt() > genTau2.pt();
  }
}

void DumpGenTaus::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  std::cout << "<DumpGenTaus::analyze> (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
  std::cout << " src = " << src_ << std::endl;

  edm::Handle<reco::GenJetCollection> genTaus;
  evt.getByToken(token_, genTaus);

  // sort genTau collection by decreasing pT
  reco::GenJetCollection genTaus_sorted = *genTaus;
  std::sort(genTaus_sorted.begin(), genTaus_sorted.end(), isHigherPt);
  
  size_t numTaus = genTaus_sorted.size();
  for ( size_t idxTau = 0; idxTau < numTaus; ++idxTau ) 
  {
    const reco::GenJet& genTau = genTaus_sorted.at(idxTau);
    std::cout << "genTau #" << idxTau  << ": pT = " << genTau.pt() << ", eta = " << genTau.eta() << ", phi = " << genTau.phi() 
	      << " (decayMode = " << JetMCTagUtils::genTauDecayMode(genTau) << ")" << std::endl;
    std::vector<const reco::GenParticle*> daughters = genTau.getGenConstituents();
    std::cout << "daughters:" << std::endl;
    for ( auto daughter : daughters ) 
    {
      std::cout << " pdgId = " << daughter->pdgId() << ": pT = " << daughter->pt() << ", eta = " << daughter->eta() << ", phi = " << daughter->phi() << std::endl;
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DumpGenTaus);





