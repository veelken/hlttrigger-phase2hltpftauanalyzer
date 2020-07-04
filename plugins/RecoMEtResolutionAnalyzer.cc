#include "HLTTrigger/TallinnHLTPFTauAnalyzer/plugins/RecoMEtResolutionAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Common/interface/Handle.h"

#include "HLTTrigger/TallinnHLTPFTauAnalyzer/interface/histogramAuxFunctions.h" // fillWithOverFlow

#include <TMath.h>

#include <iostream>
#include <iomanip>

RecoMEtResolutionAnalyzer::RecoMEtResolutionAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  srcRecMEt_ = cfg.getParameter<edm::InputTag>("srcRecMEt");
  tokenRecMEt_ = consumes<reco::PFMETCollection>(srcRecMEt_);
  srcGenMEt_ = cfg.getParameter<edm::InputTag>("srcGenMEt");
  tokenGenMEt_ = consumes<reco::GenMETCollection>(srcGenMEt_);
  
  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

RecoMEtResolutionAnalyzer::~RecoMEtResolutionAnalyzer()
{}

void RecoMEtResolutionAnalyzer::beginJob()
{
  if ( !edm::Service<dqm::legacy::DQMStore>().isAvailable() ) {
    throw cms::Exception("RecoMEtResolutionAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<dqm::legacy::DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory_.data());

  me_recMEt_Pt_ = dqmStore.book1D("recMEt_Pt", "recMEt_Pt", 100, 0., 500.);
  histogram_recMEt_Pt_ = me_recMEt_Pt_->getTH1();
  assert(histogram_recMEt_Pt_);
  me_genMEt_Pt_ = dqmStore.book1D("genMEt_Pt", "genMEt_Pt", 100, 0., 500.);
  histogram_genMEt_Pt_ = me_genMEt_Pt_->getTH1();
  assert(histogram_genMEt_Pt_);
  me_deltaMEt_Px_ = dqmStore.book1D("deltaMEt_Px", "deltaMEt_Px", 100, -250., +250.);
  histogram_deltaMEt_Px_ = me_deltaMEt_Px_->getTH1(); 
  me_deltaMEt_Py_ = dqmStore.book1D("deltaMEt_Py", "deltaMEt_Py", 100, -250., +250.);
  histogram_deltaMEt_Py_ = me_deltaMEt_Py_->getTH1(); 
}

void RecoMEtResolutionAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<reco::PFMETCollection> recMEt;
  evt.getByToken(tokenRecMEt_, recMEt);
  const reco::Candidate::LorentzVector& recMEtP4 = recMEt->front().p4();

  edm::Handle<reco::GenMETCollection> genMEt;
  evt.getByToken(tokenGenMEt_, genMEt);
  const reco::Candidate::LorentzVector& genMEtP4 = genMEt->front().p4();

  const double evtWeight = 1.;

  fillWithOverFlow(histogram_recMEt_Pt_, recMEtP4.pt(), evtWeight);
  fillWithOverFlow(histogram_genMEt_Pt_, genMEtP4.pt(), evtWeight);
  fillWithOverFlow(histogram_deltaMEt_Px_, recMEtP4.px() - genMEtP4.px(), evtWeight);
  fillWithOverFlow(histogram_deltaMEt_Py_, recMEtP4.py() - genMEtP4.py(), evtWeight);
}

void RecoMEtResolutionAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(RecoMEtResolutionAnalyzer);
