#include "HLTTrigger/TallinnHLTPFTauAnalyzer/plugins/GenPtHatAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include "HLTTrigger/TallinnHLTPFTauAnalyzer/interface/histogramAuxFunctions.h" // fillWithOverFlow

GenPtHatAnalyzer::GenPtHatAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<GenEventInfoProduct>(src_);

  lumiScale_ = cfg.getParameter<double>("lumiScale");

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

GenPtHatAnalyzer::~GenPtHatAnalyzer()
{}

void GenPtHatAnalyzer::beginJob()
{
  if ( !edm::Service<dqm::legacy::DQMStore>().isAvailable() ) {
    throw cms::Exception("GenPtHatAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<dqm::legacy::DQMStore>());

  dqmStore.setCurrentFolder(dqmDirectory_);

  me_genPtHat_ = dqmStore.book1D("genPtHat", "genPtHat", 1002, -1.5, 1000.5);
  histogram_genPtHat_ = me_genPtHat_->getTH1();
  assert(histogram_genPtHat_);
}
    
void GenPtHatAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<GenEventInfoProduct> genEventInfoProduct;
  evt.getByToken(token_, genEventInfoProduct);

  const double evtWeight = lumiScale_;

  // CV: code to access generator PtHat information taken from
  //       https://cmssdt.cern.ch/lxr/source/Calibration/HcalCalibAlgos/test/GammaJetAnalysis.cc
  //    (generator PtHat information will be valid for QCD multijet MC samples only)
  double genPtHat = -1.;
  if ( genEventInfoProduct.isValid() && genEventInfoProduct->hasBinningValues() )
  {
    genPtHat = genEventInfoProduct->binningValues()[0];
  }
  
  fillWithOverFlow(histogram_genPtHat_, genPtHat, evtWeight);
}

void GenPtHatAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(GenPtHatAnalyzer);


