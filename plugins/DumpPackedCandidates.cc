#include "HLTrigger/TallinnHLTPFTauAnalyzer/plugins/DumpPackedCandidates.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

DumpPackedCandidates::DumpPackedCandidates(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<pat::PackedCandidateCollection>(src_);
  minPt_ = ( cfg.exists("min_pt") ) ? cfg.getParameter<double>("min_pt") : -1.;
}

DumpPackedCandidates::~DumpPackedCandidates()
{}

void DumpPackedCandidates::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  std::cout << "<DumpPackedCandidates::analyze (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
  std::cout << " src = " << src_ << std::endl;

  edm::Handle<pat::PackedCandidateCollection> packedCandidates;
  evt.getByToken(token_, packedCandidates);
  
  size_t numPackedCandidates = packedCandidates->size();
  for ( size_t idxPackedCandidate = 0; idxPackedCandidate < numPackedCandidates; ++idxPackedCandidate ) {
    const pat::PackedCandidate& packedCandidate = packedCandidates->at(idxPackedCandidate);
    if ( !(packedCandidate.pt() > minPt_) ) continue;
    std::cout << "PFCandidate #" << idxPackedCandidate << ":" 
	      << " pT = " << packedCandidate.pt() << ", eta = " << packedCandidate.eta() << ", phi = " << packedCandidate.phi() << ","
	      << " pdgId = " << packedCandidate.pdgId() << std::endl;
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DumpPackedCandidates);





