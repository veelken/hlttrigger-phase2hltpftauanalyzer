#include "HLTTrigger/TallinnHLTPFTauAnalyzer/plugins/RhoAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

RhoAnalyzer::RhoAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_rho_ = cfg.getParameter<edm::InputTag>("srcRho");
  token_rho_ = consumes<double>(src_rho_);

  src_rhoNeutral_ = cfg.getParameter<edm::InputTag>("srcRhoNeutral");
  token_rhoNeutral_ = consumes<double>(src_rhoNeutral_);

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

RhoAnalyzer::~RhoAnalyzer()
{}

void RhoAnalyzer::beginJob()
{
  if ( !edm::Service<dqm::legacy::DQMStore>().isAvailable() ) {
    throw cms::Exception("RhoAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<dqm::legacy::DQMStore>());

  dqmStore.setCurrentFolder(dqmDirectory_);

  me_rho_ = dqmStore.book1D("rho", "rho", 200, 0., 100.);
  histogram_rho_ = me_rho_->getTH1();
  assert(histogram_rho_);

  me_rhoNeutral_ = dqmStore.book1D("rhoNeutral", "rhoNeutral", 200, 0., 100.);
  histogram_rhoNeutral_ = me_rhoNeutral_->getTH1();
  assert(histogram_rhoNeutral_);
}
    
void RhoAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<double> rho;
  evt.getByToken(token_rho_, rho);
  
  edm::Handle<double> rhoNeutral;
  evt.getByToken(token_rhoNeutral_, rhoNeutral);

  const double evtWeight = 1.;

  histogram_rho_->Fill(*rho, evtWeight);
  
  histogram_rhoNeutral_->Fill(*rhoNeutral, evtWeight);
}

void RhoAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(RhoAnalyzer);


