#include "HLTrigger/TallinnHLTPFTauAnalyzer/plugins/DumpPATJets.h"

#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

DumpPATJets::DumpPATJets(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<pat::JetCollection>(src_);
}

DumpPATJets::~DumpPATJets()
{}

void DumpPATJets::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  std::cout << "<DumpPATJets::analyze (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
  std::cout << " src = " << src_ << std::endl;

  edm::Handle<pat::JetCollection> jets;
  evt.getByToken(token_, jets);
  
  for ( size_t idxJet = 0; idxJet < jets->size(); ++idxJet ) {
    const pat::Jet& jet = jets->at(idxJet);
    std::cout << "jet #" << idxJet << ": pT = " << jet.pt() << ", eta = " << jet.eta() << ", phi = " << jet.phi() << "," 
	      << " mass = " << jet.mass() << std::endl;
    for ( size_t idxDaughter = 0; idxDaughter < jet.numberOfDaughters(); ++idxDaughter ) {
      const reco::Candidate* daughter = dynamic_cast<const reco::Candidate*>(jet.daughter(idxDaughter));
      std::cout << " daughter #" << idxDaughter << ": pT = " << daughter->pt() << ", eta = " << daughter->eta() << ", phi = " << daughter->phi() << "," 
		<< " pdgId = " << daughter->pdgId() << std::endl;
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DumpPATJets);





