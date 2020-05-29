#include "HLTTrigger/TallinnHLTPFTauAnalyzer/plugins/PATTauAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

PATTauAnalyzer::PATTauAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<pat::TauCollection>(src_);

  min_pt_ = cfg.getParameter<double>("min_pt");
  max_pt_ = cfg.getParameter<double>("max_pt");
  min_absEta_ = cfg.getParameter<double>("min_absEta");
  max_absEta_ = cfg.getParameter<double>("max_absEta");

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

PATTauAnalyzer::~PATTauAnalyzer()
{
  for ( auto plot : plots_ ) 
  {
    delete plot;
  }
}

void PATTauAnalyzer::beginJob()
{
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    throw cms::Exception("PATTauAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory_.data());

  std::vector<std::string> decayModes = { "oneProng0Pi0", "oneProng1Pi0", "oneProng2Pi0", "threeProng0Pi0", "threeProng1Pi0", "all" };
  for ( auto decayMode : decayModes )
  {    
    plots_.push_back(new plotEntryType(min_pt_, max_pt_, min_absEta_, max_absEta_, decayMode, 0)); // leadingTau
    plots_.push_back(new plotEntryType(min_pt_, max_pt_, min_absEta_, max_absEta_, decayMode, 1)); // subleadingTau
  }

  for ( auto plot : plots_ ) 
  {
    plot->bookHistograms(dqmStore);
  }

  TString histogramName_deltaR = "deltaR";
  me_deltaR_ = dqmStore.book1D(histogramName_deltaR.Data(), histogramName_deltaR.Data(), 50, 0., 5.); 
  histogram_deltaR_ = me_deltaR_->getTH1();
  assert(histogram_deltaR_);
}

namespace
{
  bool
  isHigherPt(const pat::Tau& offlineTau1,
	     const pat::Tau& offlineTau2)
  {
    return offlineTau1.pt() > offlineTau2.pt();
  }
}

void PATTauAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<pat::TauCollection> offlineTaus;
  evt.getByToken(token_, offlineTaus);
  if ( !(offlineTaus->size() >= 2) ) return;

  // sort offlineTau collection by decreasing pT
  pat::TauCollection offlineTaus_sorted = *offlineTaus;
  std::sort(offlineTaus_sorted.begin(), offlineTaus_sorted.end(), isHigherPt);

  const double evtWeight = 1.;

  for ( auto plot : plots_ ) 
  {    
    plot->fillHistograms(offlineTaus_sorted, evtWeight);
  }

  const pat::Tau& offlineTau_lead = offlineTaus_sorted[0];
  const pat::Tau& offlineTau_sublead = offlineTaus_sorted[1];
  double dR = reco::deltaR(offlineTau_lead.eta(), offlineTau_lead.phi(), offlineTau_sublead.eta(), offlineTau_sublead.phi());
  histogram_deltaR_->Fill(dR, evtWeight);
}

void PATTauAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(PATTauAnalyzer);
