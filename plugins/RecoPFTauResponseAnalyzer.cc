#include "HLTTrigger/TallinnHLTPFTauAnalyzer/plugins/RecoPFTauResponseAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

RecoPFTauResponseAnalyzer::RecoPFTauResponseAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  srcPFTaus_ = cfg.getParameter<edm::InputTag>("srcPFTaus");
  tokenPFTaus_ = consumes<reco::PFTauCollection>(srcPFTaus_);
  srcPFTauSumChargedIso_ = cfg.getParameter<edm::InputTag>("srcPFTauSumChargedIso");
  tokenPFTauSumChargedIso_ = consumes<reco::PFTauDiscriminator>(srcPFTauSumChargedIso_);
  srcRefTaus_ = cfg.getParameter<edm::InputTag>("srcRefTaus");
  std::string typeRefTaus_string = cfg.getParameter<std::string>("typeRefTaus");
  if ( typeRefTaus_string == "gen"     ) 
  {
    typeRefTaus_ = kGen;
    tokenRefTaus_gen_ = consumes<reco::GenJetCollection>(srcRefTaus_);
  }
  else if ( typeRefTaus_string == "offline" ) 
  {
    typeRefTaus_ = kOffline;
    tokenRefTaus_offline_ = consumes<pat::TauCollection>(srcRefTaus_);
  }
  else
  {
    throw cms::Exception("TallinnL1PFTauResponseAnalyzer") 
      << " Invalid Configuration parameter 'typeRefTaus' = " << typeRefTaus_string << " !!\n";;
  }

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

RecoPFTauResponseAnalyzer::~RecoPFTauResponseAnalyzer()
{
  for ( auto responsePlot : responsePlots_ ) 
  {
    delete responsePlot;
  }
}

void RecoPFTauResponseAnalyzer::beginJob()
{
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    throw cms::Exception("RecoPFTauResponseAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory_.data());

  std::vector<std::string> decayModes = { "oneProng0Pi0", "oneProng1Pi0", "oneProng2Pi0", "threeProng0Pi0", "threeProng1Pi0", "all" };
  std::vector<double> min_pts = { 20., 25., 30., 35., 40. };
  for ( auto decayMode : decayModes )
  {
    for ( auto min_pt : min_pts )
    {
      responsePlots_.push_back(new responsePlotEntryType(min_pt, 1.0, decayMode, 0.40, -1.)); // vLoose
      responsePlots_.push_back(new responsePlotEntryType(min_pt, 1.0, decayMode, 0.20, -1.)); // Loose
      responsePlots_.push_back(new responsePlotEntryType(min_pt, 1.0, decayMode, 0.10, -1.)); // Medium
      responsePlots_.push_back(new responsePlotEntryType(min_pt, 1.0, decayMode, 0.05, -1.)); // Tight
    
      responsePlots_.push_back(new responsePlotEntryType(min_pt, 1.4, decayMode, 0.40, -1.)); // vLoose
      responsePlots_.push_back(new responsePlotEntryType(min_pt, 1.4, decayMode, 0.20, -1.)); // Loose
      responsePlots_.push_back(new responsePlotEntryType(min_pt, 1.4, decayMode, 0.10, -1.)); // Medium
      responsePlots_.push_back(new responsePlotEntryType(min_pt, 1.4, decayMode, 0.05, -1.)); // Tight
    }
  }

  for ( auto responsePlot : responsePlots_ ) 
  {
    responsePlot->bookHistograms(dqmStore);
  }
}

void RecoPFTauResponseAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<reco::PFTauCollection> pfTaus;
  evt.getByToken(tokenPFTaus_, pfTaus);
  
  edm::Handle<reco::PFTauDiscriminator> pfTauSumChargedIso;
  evt.getByToken(tokenPFTauSumChargedIso_, pfTauSumChargedIso);

  const double evtWeight = 1.;

  if ( typeRefTaus_ == kGen )
  {
    edm::Handle<reco::GenJetCollection> refTaus_gen;
    evt.getByToken(tokenRefTaus_gen_, refTaus_gen);

    for ( auto responsePlot : responsePlots_ ) 
    {    
      responsePlot->fillHistograms(pfTaus, *pfTauSumChargedIso, *refTaus_gen, evtWeight);
    }
  }

  if ( typeRefTaus_ == kOffline )
  {
    edm::Handle<pat::TauCollection> refTaus_offline;
    evt.getByToken(tokenRefTaus_offline_, refTaus_offline);

    for ( auto responsePlot : responsePlots_ ) 
    {    
      responsePlot->fillHistograms(pfTaus, *pfTauSumChargedIso, *refTaus_offline, evtWeight);
    }
  }
}

void RecoPFTauResponseAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(RecoPFTauResponseAnalyzer);
