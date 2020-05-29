#include "HLTTrigger/TallinnHLTPFTauAnalyzer/plugins/RecoPFTauAnalyzerBackground.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

RecoPFTauAnalyzerBackground::RecoPFTauAnalyzerBackground(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  srcPFTaus_ = cfg.getParameter<edm::InputTag>("srcPFTaus");
  tokenPFTaus_ = consumes<reco::PFTauCollection>(srcPFTaus_);
  srcPFTauSumChargedIso_ = cfg.getParameter<edm::InputTag>("srcPFTauSumChargedIso");
  tokenPFTauSumChargedIso_ = consumes<reco::PFTauDiscriminator>(srcPFTauSumChargedIso_);

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

RecoPFTauAnalyzerBackground::~RecoPFTauAnalyzerBackground()
{
  for ( auto ratePlot : ratePlots_ ) 
  {
    delete ratePlot;
  }
}

void RecoPFTauAnalyzerBackground::beginJob()
{
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    throw cms::Exception("RecoPFTauAnalyzerBackground") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory_.data());

  ratePlots_.push_back(new ratePlotEntryType(-1.,  1.4,   0.40, -1., 0.4)); // vLoose
  ratePlots_.push_back(new ratePlotEntryType(-1.,  1.4,   0.20, -1., 0.4)); // Loose
  ratePlots_.push_back(new ratePlotEntryType(-1.,  1.4,   0.10, -1., 0.4)); // Medium
  ratePlots_.push_back(new ratePlotEntryType(-1.,  1.4,   0.05, -1., 0.4)); // Tight
  ratePlots_.push_back(new ratePlotEntryType(-1.,  1.4,   0.02, -1., 0.4)); // vTight
  ratePlots_.push_back(new ratePlotEntryType(-1.,  1.4,   0.01, -1., 0.4)); // vvTight
  
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.172, 0.40, -1., 0.4)); // vLoose
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.172, 0.20, -1., 0.4)); // Loose
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.172, 0.10, -1., 0.4)); // Medium
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.172, 0.05, -1., 0.4)); // Tight
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.172, 0.02, -1., 0.4)); // vTight
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.172, 0.01, -1., 0.4)); // vvTight

  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.4,   0.40, -1., 0.4)); // vLoose
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.4,   0.20, -1., 0.4)); // Loose
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.4,   0.10, -1., 0.4)); // Medium
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.4,   0.05, -1., 0.4)); // Tight
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.4,   0.02, -1., 0.4)); // vTight
  ratePlots_.push_back(new ratePlotEntryType( 1.4, 2.4,   0.01, -1., 0.4)); // vvTight

  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.172, 0.40, -1., 0.4)); // vLoose
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.172, 0.20, -1., 0.4)); // Loose
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.172, 0.10, -1., 0.4)); // Medium
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.172, 0.05, -1., 0.4)); // Tight
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.172, 0.02, -1., 0.4)); // vTight
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.172, 0.01, -1., 0.4)); // vvTight

  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.4,   0.40, -1., 0.4)); // vLoose
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.4,   0.20, -1., 0.4)); // Loose
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.4,   0.10, -1., 0.4)); // Medium
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.4,   0.05, -1., 0.4)); // Tight
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.4,   0.02, -1., 0.4)); // vTight
  ratePlots_.push_back(new ratePlotEntryType(-1.,  2.4,   0.01, -1., 0.4)); // vvTight

  for ( auto ratePlot : ratePlots_ ) 
  {
    ratePlot->bookHistograms(dqmStore);
  }
}

void RecoPFTauAnalyzerBackground::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<reco::PFTauCollection> pfTaus;
  evt.getByToken(tokenPFTaus_, pfTaus);
  
  edm::Handle<reco::PFTauDiscriminator> pfTauSumChargedIso;
  evt.getByToken(tokenPFTauSumChargedIso_, pfTauSumChargedIso);
  
  const double evtWeight = 1.;

  for ( auto ratePlot : ratePlots_ ) 
  {
    ratePlot->fillHistograms(pfTaus, *pfTauSumChargedIso, evtWeight);
  }
}

void RecoPFTauAnalyzerBackground::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(RecoPFTauAnalyzerBackground);
