#include "HLTTrigger/TallinnHLTPFTauAnalyzer/plugins/DumpRecoPFCandidates.h"

#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

DumpRecoPFCandidates::DumpRecoPFCandidates(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<reco::PFCandidateCollection>(src_);
  minPt_ = ( cfg.exists("min_pt") ) ? cfg.getParameter<double>("min_pt") : -1.;
  maxAbsEta_ = ( cfg.exists("max_absEta") ) ? cfg.getParameter<double>("max_absEta") : 1.e+3;
}

DumpRecoPFCandidates::~DumpRecoPFCandidates()
{}

void DumpRecoPFCandidates::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  std::cout << "<DumpRecoPFCandidates::analyze (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
  std::cout << " src = " << src_ << std::endl;

  edm::Handle<reco::PFCandidateCollection> pfCands;
  evt.getByToken(token_, pfCands);
  
  reco::Candidate::LorentzVector sumP4;

  size_t numPFCands = pfCands->size();
  for ( size_t idxPFCand = 0; idxPFCand < numPFCands; ++idxPFCand ) 
  {
    const reco::PFCandidate& pfCand = pfCands->at(idxPFCand);
    if ( (minPt_ < 0. || pfCand.pt() > minPt_) && (maxAbsEta_ < 0. || std::fabs(pfCand.eta()) < maxAbsEta_) )
    {
      std::string type_string;
      if      ( pfCand.particleId() == reco::PFCandidate::h     ) type_string = "PFChargedHadron";
      else if ( pfCand.particleId() == reco::PFCandidate::e     ) type_string = "PFElectron";
      else if ( pfCand.particleId() == reco::PFCandidate::h0    ) type_string = "PFNeutralHadron";
      else if ( pfCand.particleId() == reco::PFCandidate::gamma ) type_string = "PFPhoton";
      else if ( pfCand.particleId() == reco::PFCandidate::mu    ) type_string = "PFMuon";
      else                                                        type_string = "N/A";
      std::cout << "PFCandidate #" << idxPFCand << " (type = " << type_string << "):" 
	        << " pT = " << pfCand.pt() << ", eta = " << pfCand.eta() << ", phi = " << pfCand.phi() << std::endl;
      std::cout << " vertex: x = " << pfCand.vertex().x() << ", y = " << pfCand.vertex().y() << ", z = " << pfCand.vertex().z() << std::endl;
      if ( pfCand.bestTrack() )
      {
        const reco::Track* bestTrack = pfCand.bestTrack();
        std::cout << " bestTrack: pT = " << bestTrack->pt() << ", eta = " << bestTrack->eta() << ", phi = " << bestTrack->phi() << std::endl;
        std::cout << " vertex: x = " << bestTrack->vertex().x() << ", y = " << bestTrack->vertex().y() << ", z = " << bestTrack->vertex().z() << std::endl;
      }
      std::cout << " calorimeter energy: ECAL = " << pfCand.ecalEnergy() << ", HCAL = " << pfCand.hcalEnergy() << std::endl;
    }
    sumP4 += pfCand.p4();
  }
  std::cout << "(sum: E = " << sumP4.energy() << ", pT = " << sumP4.pt() << ", eta = " << sumP4.eta() << ", phi = " << sumP4.phi() << ")" << std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DumpRecoPFCandidates);





