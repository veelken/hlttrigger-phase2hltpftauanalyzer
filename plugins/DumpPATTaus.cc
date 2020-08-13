#include "HLTrigger/TallinnHLTPFTauAnalyzer/plugins/DumpPATTaus.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h" // reco::PFCandidate::h, ::e, ::mu, ::gamma, ::h0

#include <TMath.h>

#include <iostream>
#include <iomanip>

DumpPATTaus::DumpPATTaus(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<pat::TauCollection>(src_);
}

DumpPATTaus::~DumpPATTaus()
{}

namespace
{
  //void printPFCand(const reco::PFCandidate& pfCand, const reco::Candidate::Point& primaryVertexPos)
  //{
  //  int particleId_int = pfCand.particleId();
  //  std::string particleId_string;
  //  if      ( particleId_int == reco::PFCandidate::h     ) particleId_string = "PFChargedHadron";
  //  else if ( particleId_int == reco::PFCandidate::e     ) particleId_string = "PFElectron";
  //  else if ( particleId_int == reco::PFCandidate::mu    ) particleId_string = "PFMuon";
  //  else if ( particleId_int == reco::PFCandidate::gamma ) particleId_string = "PFGamma";
  //  else if ( particleId_int == reco::PFCandidate::h0    ) particleId_string = "PFNeutralHadron";
  //  else                                                   particleId_string = "N/A";
  //  std::cout << " " << particleId_string << " with "
  //	        << " pT = " << pfCand.pt() << ", eta = " << pfCand.eta() << ", phi = " << pfCand.phi();
  //  if ( pfCand.charge() != 0 )
  //  {
  //    std::cout << " (dz = " << std::fabs(pfCand.vertex().z() - primaryVertexPos.z()) << ")";
  //  }
  //  std::cout << std::endl;
  //}
  //
  //void printPFCands(const std::vector<reco::PFCandidatePtr>& pfCands, const reco::Candidate::Point& primaryVertexPos)
  //{
  //  for ( auto pfCand : pfCands )
  //  {
  //    printPFCand(*pfCand, primaryVertexPos);
  //  }
  //}

  void printPackedPFCand(const reco::Candidate& pfCand, const reco::Candidate::Point& primaryVertexPos)
  {
    std::cout << " pdgId = " << pfCand.pdgId() << ":"
	      << " pT = " << pfCand.pt() << ", eta = " << pfCand.eta() << ", phi = " << pfCand.phi();
    if ( pfCand.charge() != 0 )
    {
      std::cout << " (dz = " << std::fabs(pfCand.vertex().z() - primaryVertexPos.z()) << ")";
    }
    std::cout << std::endl;
  }

  void printPackedPFCands(const reco::CandidatePtrVector& pfCands, const reco::Candidate::Point& primaryVertexPos)
  {
    for ( auto pfCand : pfCands )
    {
      printPackedPFCand(*pfCand, primaryVertexPos);
    }
  }
}

void DumpPATTaus::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  std::cout << "<DumpPATTaus::analyze (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
  std::cout << " src = " << src_ << std::endl;

  edm::Handle<pat::TauCollection> taus;
  evt.getByToken(token_, taus);

  std::vector<std::string> discriminatorsToCheck;
  discriminatorsToCheck.push_back(std::string("decayModeFinding"));
  discriminatorsToCheck.push_back(std::string("decayModeFindingNewDMs"));
  discriminatorsToCheck.push_back(std::string("byLooseIsolationMVArun2v1DBoldDMwLT2017v2"));
  discriminatorsToCheck.push_back(std::string("byMediumIsolationMVArun2v1DBoldDMwLT2017v2"));
  discriminatorsToCheck.push_back(std::string("byTightIsolationMVArun2v1DBoldDMwLT2017v2"));
  discriminatorsToCheck.push_back(std::string("chargedIsoPtSum"));
  discriminatorsToCheck.push_back(std::string("neutralIsoPtSum"));
  discriminatorsToCheck.push_back(std::string("byDeepTau2017v2VSjetraw"));
  
  size_t numTaus = taus->size();
  for ( size_t idxTau = 0; idxTau < numTaus; ++idxTau ) 
  {
    const pat::Tau& tau = taus->at(idxTau);

    std::cout << "offlinePFTau #" << idxTau << ":" 
	      << " pT = " << tau.pt() << ", eta = " << tau.eta() << ", phi = " << tau.phi() << "," 
	      << " mass = " << tau.mass() << " (charge = " << tau.charge() << ")" << std::endl;
    int tauDecayMode_int = tau.decayMode();
    std::string tauDecayMode_string;
    if      ( tauDecayMode_int == reco::PFTau::kOneProng0PiZero   ) tauDecayMode_string = "OneProng0PiZero";
    else if ( tauDecayMode_int == reco::PFTau::kOneProng1PiZero   ) tauDecayMode_string = "OneProng1PiZero";
    else if ( tauDecayMode_int == reco::PFTau::kOneProng2PiZero   ) tauDecayMode_string = "OneProng2PiZero";
    else if ( tauDecayMode_int == reco::PFTau::kTwoProng0PiZero   ) tauDecayMode_string = "TwoProng0PiZero";
    else if ( tauDecayMode_int == reco::PFTau::kTwoProng1PiZero   ) tauDecayMode_string = "TwoProng1PiZero";
    else if ( tauDecayMode_int == reco::PFTau::kThreeProng0PiZero ) tauDecayMode_string = "ThreeProng0PiZero";
    else if ( tauDecayMode_int == reco::PFTau::kThreeProng0PiZero ) tauDecayMode_string = "ThreeProng1PiZero";
    else tauDecayMode_string = "Rare";
    std::cout << "decay mode = " << tauDecayMode_string << std::endl;
    //std::cout << "associated jet:" 
    //	        << " pT = " << tau.p4Jet().pt() << ", eta = " << tau.p4Jet().eta() << ", phi = " << tau.p4Jet().phi() << "," 
    //	        << " mass = " << tau.p4Jet().mass() << std::endl;
    //std::cout << "vertex: x = " << tau.vertex().x() << ", y = " << tau.vertex().y() << ", z = " << tau.vertex().z() << std::endl;
    //if ( tau.leadTauChargedHadronCandidate().isNonnull() ) 
    //{
    //  std::cout << "lead. ChargedPFCand:" << std::endl;
    //  tau.leadTauChargedHadronCandidate()->print(std::cout);
    //}
    std::cout << "lead. ChargedPFCand:" << std::endl;
    if ( tau.leadChargedHadrCand().isNonnull() ) 
    {
      printPackedPFCand(*tau.leadChargedHadrCand(), tau.vertex());
    }
    else 
    {
      std::cout << " N/A" << std::endl;
    }
    std::cout << "signalPFCands:" << std::endl;
    //printPFCandidates(tau.signalPFCands(), tau.vertex());
    printPackedPFCands(tau.signalCands(), tau.vertex());
    std::cout << "isolationPFCandidates:" << std::endl;
    //printPFCandidates(tau.isolationPFCands(), tau.vertex());
    printPackedPFCands(tau.isolationCands(), tau.vertex());
    std::cout << "discriminators:" << std::endl;
    for ( auto discriminatorToCheck : discriminatorsToCheck ) 
    {
      if ( tau.isTauIDAvailable(discriminatorToCheck) )
      {
    	std::cout << "  " << discriminatorToCheck << " = " << tau.tauID(discriminatorToCheck) << std::endl;
      }
      else
      {
    	std::cout << "  " << discriminatorToCheck << ": N/A" << std::endl;
      }
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DumpPATTaus);





