#include "HLTrigger/TallinnHLTPFTauAnalyzer/plugins/GenHadTauAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Math/interface/deltaR.h" // reco::deltaR

#include <TMath.h>

#include <iostream>
#include <iomanip>

GenHadTauAnalyzer::GenHadTauAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<reco::GenJetCollection>(src_);

  min_pt_ = cfg.getParameter<double>("min_pt");
  max_pt_ = cfg.getParameter<double>("max_pt");
  min_absEta_ = cfg.getParameter<double>("min_absEta");
  max_absEta_ = cfg.getParameter<double>("max_absEta");

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

GenHadTauAnalyzer::~GenHadTauAnalyzer()
{
  for ( auto plot : plots_ ) 
  {
    delete plot;
  }
}

void GenHadTauAnalyzer::beginJob()
{
  if ( !edm::Service<dqm::legacy::DQMStore>().isAvailable() ) {
    throw cms::Exception("GenHadTauAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<dqm::legacy::DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory_.data());

  std::vector<std::string> decayModes = { "oneProng0Pi0", "oneProng1Pi0", "oneProng2Pi0", "threeProng0Pi0", "threeProng1Pi0", "all" };
  for ( auto decayMode : decayModes )
  {    
    plots_.push_back(new plotEntryType(min_pt_, max_pt_, min_absEta_, max_absEta_, decayMode, 0)); // leadingHadTau
    plots_.push_back(new plotEntryType(min_pt_, max_pt_, min_absEta_, max_absEta_, decayMode, 1)); // subleadingHadTau
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
  isHigherPt(const reco::GenJet& genHadTau1,
	     const reco::GenJet& genHadTau2)
  {
    return genHadTau1.pt() > genHadTau2.pt();
  }
}

void GenHadTauAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<reco::GenJetCollection> genHadTaus;
  evt.getByToken(token_, genHadTaus);
  if ( !(genHadTaus->size() >= 2) ) return;
  
  // sort genTau collection by decreasing pT
  reco::GenJetCollection genHadTaus_sorted = *genHadTaus;
  std::sort(genHadTaus_sorted.begin(), genHadTaus_sorted.end(), isHigherPt);

  const double evtWeight = 1.;

  for ( auto plot : plots_ ) 
  {    
    plot->fillHistograms(genHadTaus_sorted, evtWeight);
  }

  const reco::GenJet& genHadTau_lead = genHadTaus_sorted[0];
  const reco::GenJet& genHadTau_sublead = genHadTaus_sorted[1];
  double dR = reco::deltaR(genHadTau_lead.eta(), genHadTau_lead.phi(), genHadTau_sublead.eta(), genHadTau_sublead.phi());
  histogram_deltaR_->Fill(dR, evtWeight);
}

void GenHadTauAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(GenHadTauAnalyzer);
