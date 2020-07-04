#include "HLTTrigger/TallinnHLTPFTauAnalyzer/plugins/RecoPFCandidateConeEnResponseAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

RecoPFCandidateConeEnResponseAnalyzer::RecoPFCandidateConeEnResponseAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  srcPFCandidates_ = cfg.getParameter<edm::InputTag>("srcPFCandidates");
  tokenPFCandidates_ = consumes<reco::PFCandidateCollection>(srcPFCandidates_);
  srcTracks_ = cfg.getParameter<edm::InputTag>("srcTracks");
  tokenTracks_ = consumes<reco::TrackCollection>(srcTracks_);
  srcGenParticles_ = cfg.getParameter<edm::InputTag>("srcGenParticles");
  tokenGenParticles_ = consumes<edm::View<reco::Candidate>>(srcGenParticles_);
  srcGenHadTaus_ = cfg.getParameter<edm::InputTag>("srcGenHadTaus");
  tokenGenHadTaus_ = consumes<reco::GenJetCollection>(srcGenHadTaus_);

  dRcones_ = cfg.getParameter<std::vector<double>>("dRcones");

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

RecoPFCandidateConeEnResponseAnalyzer::~RecoPFCandidateConeEnResponseAnalyzer()
{
  for ( auto responsePlot : responsePlots_ ) 
  {
    delete responsePlot;
  }
}

void RecoPFCandidateConeEnResponseAnalyzer::beginJob()
{
  if ( !edm::Service<dqm::legacy::DQMStore>().isAvailable() ) {
    throw cms::Exception("RecoPFCandidateConeEnResponseAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<dqm::legacy::DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory_.data());

  std::vector<std::string> decayModes = { "oneProng0Pi0", "oneProng1Pi0", "oneProng2Pi0", "threeProng0Pi0", "threeProng1Pi0", "all" };
  for ( auto decayMode : decayModes )
  {
    TString dqmDirectory = dqmDirectory_.data();
    std::string decayMode_capitalized = decayMode;
    decayMode_capitalized[0] = toupper(decayMode_capitalized[0]);	
    dqmDirectory.Append(Form("/gen%sTau", decayMode_capitalized.data()));

    dqmStore.setCurrentFolder(dqmDirectory.Data());
    for ( auto dRcone : dRcones_ )
    {
      responsePlotEntryType* responsePlots = new responsePlotEntryType(decayMode, dRcone);
      responsePlots->bookHistograms(dqmStore);
      responsePlots_.push_back(responsePlots);
    }
  }
}

void RecoPFCandidateConeEnResponseAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<reco::PFCandidateCollection> pfCands;
  evt.getByToken(tokenPFCandidates_, pfCands);

  edm::Handle<reco::TrackCollection> tracks;
  evt.getByToken(tokenTracks_, tracks);

  edm::Handle<edm::View<reco::Candidate>> genParticles;
  evt.getByToken(tokenGenParticles_, genParticles);

  edm::Handle<reco::GenJetCollection> genHadTaus;
  evt.getByToken(tokenGenHadTaus_, genHadTaus);

  const double evtWeight = 1.;

  for ( auto responsePlot : responsePlots_ ) 
  {    
    responsePlot->fillHistograms(*pfCands, *tracks, *genParticles, *genHadTaus, evtWeight);
  }
}

void RecoPFCandidateConeEnResponseAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(RecoPFCandidateConeEnResponseAnalyzer);
