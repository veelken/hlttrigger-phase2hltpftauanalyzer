#include "HLTrigger/TallinnHLTPFTauAnalyzer/plugins/EvtWeightAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include "HLTrigger/TallinnHLTPFTauAnalyzer/interface/histogramAuxFunctions.h" // fillWithOverFlow

#include <TMath.h> // TMath::Log10, TMath::Max

EvtWeightAnalyzer::EvtWeightAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<double>(src_);

  numBins_evtWeight_ = cfg.getParameter<int>("numBins_evtWeight");
  xMin_evtWeight_ = cfg.getParameter<double>("xMin_evtWeight");
  xMax_evtWeight_ = cfg.getParameter<double>("xMax_evtWeight");

  numBins_log10evtWeight_ = cfg.getParameter<int>("numBins_log10evtWeight");
  xMin_log10evtWeight_ = cfg.getParameter<double>("xMin_log10evtWeight");
  xMax_log10evtWeight_ = cfg.getParameter<double>("xMax_log10evtWeight");

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

EvtWeightAnalyzer::~EvtWeightAnalyzer()
{}

void EvtWeightAnalyzer::beginJob()
{
  if ( !edm::Service<dqm::legacy::DQMStore>().isAvailable() ) {
    throw cms::Exception("EvtWeightAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<dqm::legacy::DQMStore>());

  dqmStore.setCurrentFolder(dqmDirectory_);

  me_evtWeight_ = dqmStore.book1D("evtWeight", "evtWeight", numBins_evtWeight_, xMin_evtWeight_, xMax_evtWeight_);
  histogram_evtWeight_ = me_evtWeight_->getTH1();
  assert(histogram_evtWeight_);

  me_log10evtWeight_ = dqmStore.book1D("log10evtWeight", "log10evtWeight", numBins_log10evtWeight_, xMin_log10evtWeight_, xMax_log10evtWeight_);
  histogram_log10evtWeight_ = me_log10evtWeight_->getTH1();
  assert(histogram_log10evtWeight_);
}
    
void EvtWeightAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<double> evtWeight;
  evt.getByToken(token_, evtWeight);

  fillWithOverFlow(histogram_evtWeight_, *evtWeight, 1.);

  double log10evtWeight = TMath::Log10(TMath::Max(1.e-37, *evtWeight));
  fillWithOverFlow(histogram_log10evtWeight_, log10evtWeight, 1.);
}

void EvtWeightAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(EvtWeightAnalyzer);


