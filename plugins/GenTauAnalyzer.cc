#include "HLTTrigger/TallinnHLTPFTauAnalyzer/plugins/GenTauAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Math/interface/deltaR.h" // reco::deltaR

#include <TMath.h>

#include <iostream>
#include <iomanip>

GenTauAnalyzer::GenTauAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<reco::GenParticleCollection>(src_);

  min_pt_ = cfg.getParameter<double>("min_pt");
  max_pt_ = cfg.getParameter<double>("max_pt");
  min_absEta_ = cfg.getParameter<double>("min_absEta");
  max_absEta_ = cfg.getParameter<double>("max_absEta");

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

GenTauAnalyzer::~GenTauAnalyzer()
{
  for ( auto plot : plots_ ) 
  {
    delete plot;
  }
}

void GenTauAnalyzer::beginJob()
{
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    throw cms::Exception("GenTauAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory_.data());

  plots_.push_back(new plotEntryType(min_pt_, max_pt_, min_absEta_, max_absEta_, 0)); // leadingTau
  plots_.push_back(new plotEntryType(min_pt_, max_pt_, min_absEta_, max_absEta_, 1)); // subleadingTau

  for ( auto plot : plots_ ) 
  {
    plot->bookHistograms(dqmStore);
  }
  
  TString histogramName_diTau_pt = "diTau_pt";
  me_diTau_pt_ = dqmStore.book1D(histogramName_diTau_pt.Data(), histogramName_diTau_pt.Data(), 50, 0., 500.);
  histogram_diTau_pt_ = me_diTau_pt_->getTH1();
  assert(histogram_diTau_pt_);
  TString histogramName_diTau_eta = "diTau_eta";
  me_diTau_eta_ = dqmStore.book1D(histogramName_diTau_eta.Data(), histogramName_diTau_eta.Data(), 30, -3., +3.);
  histogram_diTau_eta_ = me_diTau_eta_->getTH1();
  assert(histogram_diTau_eta_);
  TString histogramName_diTau_phi = "diTau_phi";
  me_diTau_phi_ = dqmStore.book1D(histogramName_diTau_phi.Data(), histogramName_diTau_phi.Data(), 18, -TMath::Pi(), +TMath::Pi());
  histogram_diTau_phi_ = me_diTau_phi_->getTH1();
  assert(histogram_diTau_phi_);
  TString histogramName_diTau_mass = "diTau_mass";
  me_diTau_mass_ = dqmStore.book1D(histogramName_diTau_mass.Data(), histogramName_diTau_mass.Data(), 250, 0., 250.);
  histogram_diTau_mass_ = me_diTau_mass_->getTH1();
  assert(histogram_diTau_mass_);
  TString histogramName_deltaR = "deltaR";
  me_deltaR_ = dqmStore.book1D(histogramName_deltaR.Data(), histogramName_deltaR.Data(), 50, 0., 5.); 
  histogram_deltaR_ = me_deltaR_->getTH1();
  assert(histogram_deltaR_);
}

namespace
{
  bool
  isHigherPt(const reco::GenParticle& genTau1,
	     const reco::GenParticle& genTau2)
  {
    return genTau1.pt() > genTau2.pt();
  }
}

void GenTauAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<reco::GenParticleCollection> genTaus;
  evt.getByToken(token_, genTaus);
  if ( !(genTaus->size() >= 2) ) return;
  
  // sort genTau collection by decreasing pT
  reco::GenParticleCollection genTaus_sorted = *genTaus;
  std::sort(genTaus_sorted.begin(), genTaus_sorted.end(), isHigherPt);

  const double evtWeight = 1.;

  for ( auto plot : plots_ ) 
  {    
    plot->fillHistograms(genTaus_sorted, evtWeight);
  }

  const reco::GenParticle& genTau_lead = genTaus_sorted[0];
  const reco::GenParticle& genTau_sublead = genTaus_sorted[1];
  reco::Candidate::LorentzVector genDiTauP4 = genTau_lead.p4() + genTau_sublead.p4();
  histogram_diTau_pt_->Fill(genDiTauP4.pt(), evtWeight);
  histogram_diTau_eta_->Fill(genDiTauP4.eta(), evtWeight);
  histogram_diTau_phi_->Fill(genDiTauP4.phi(), evtWeight);
  histogram_diTau_mass_->Fill(genDiTauP4.mass(), evtWeight);
  double dR = reco::deltaR(genTau_lead.eta(), genTau_lead.phi(), genTau_sublead.eta(), genTau_sublead.phi());
  histogram_deltaR_->Fill(dR, evtWeight);
}

void GenTauAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(GenTauAnalyzer);
