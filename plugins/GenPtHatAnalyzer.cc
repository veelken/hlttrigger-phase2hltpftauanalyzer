#include "HLTrigger/TallinnHLTPFTauAnalyzer/plugins/GenPtHatAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include "HLTrigger/TallinnHLTPFTauAnalyzer/interface/histogramAuxFunctions.h" // fillWithOverFlow

GenPtHatAnalyzer::GenPtHatAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_genEventInfo_ = cfg.getParameter<edm::InputTag>("src_genEventInfo");
  token_genEventInfo_ = consumes<GenEventInfoProduct>(src_genEventInfo_);

  src_genJets_ = cfg.getParameter<edm::InputTag>("src_genJets");
  token_genJets_ = consumes<reco::GenJetCollection>(src_genJets_);

  src_pileupSummaryInfo_ = cfg.getParameter<edm::InputTag>("src_pileupSummaryInfo");
  token_pileupSummaryInfo_ = consumes<std::vector<PileupSummaryInfo>>(src_pileupSummaryInfo_);

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

  me_genPtHat_hardscatter_ = dqmStore.book1D("genPtHat_hardscatter", "genPtHat_hardscatter", 1002, -2., 1000.);
  histogram_genPtHat_hardscatter_ = me_genPtHat_hardscatter_->getTH1();
  assert(histogram_genPtHat_hardscatter_);

  me_leadGenJetPt_vs_genPtHat_ = dqmStore.book2D("leadGenJetPt_vs_genPtHat", "leadGenJetPt_vs_genPtHat", 60, 0., 300., 60, 0., 300.);
  histogram_leadGenJetPt_vs_genPtHat_ = dynamic_cast<TH2*>(me_leadGenJetPt_vs_genPtHat_->getTH1());
  assert(histogram_leadGenJetPt_vs_genPtHat_);

  me_subleadGenJetPt_vs_genPtHat_ = dqmStore.book2D("subleadGenJetPt_vs_genPtHat", "subleadGenJetPt_vs_genPtHat", 60, 0., 300., 60, 0., 300.);
  histogram_subleadGenJetPt_vs_genPtHat_ = dynamic_cast<TH2*>(me_subleadGenJetPt_vs_genPtHat_->getTH1());
  assert(histogram_subleadGenJetPt_vs_genPtHat_);

  me_genPtHat_pileup_ = dqmStore.book1D("genPtHat_pileup", "genPtHat_pileup", 1002, -2., 1000.);
  histogram_genPtHat_pileup_ = me_genPtHat_pileup_->getTH1();
  assert(histogram_genPtHat_pileup_);
}
    
namespace
{
  bool
  isHigherPt(const reco::GenJet* genJet1,
             const reco::GenJet* genJet2)
  {
    return genJet1->pt() > genJet2->pt();
  }
}

void GenPtHatAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<GenEventInfoProduct> genEventInfoProduct;
  evt.getByToken(token_genEventInfo_, genEventInfoProduct);

  edm::Handle<reco::GenJetCollection> genJets;
  evt.getByToken(token_genJets_, genJets);

  edm::Handle<std::vector<PileupSummaryInfo>> pileupSummaryInfos;
  evt.getByToken(token_pileupSummaryInfo_, pileupSummaryInfos);

  const double evtWeight = lumiScale_;

  double genPtHat_hardscatter = -1.e+3;
  bool genPtHat_hardscatter_isValid = false;
  if ( genEventInfoProduct.isValid() )
  {
    genPtHat_hardscatter = genEventInfoProduct->qScale();
    genPtHat_hardscatter_isValid = true;
  }
  
  fillWithOverFlow(histogram_genPtHat_hardscatter_, genPtHat_hardscatter, evtWeight);

  std::vector<const reco::GenJet*> genJets_sorted;
  for ( reco::GenJetCollection::const_iterator genJet = genJets->begin();
        genJet != genJets->end(); ++genJet ) {
    genJets_sorted.push_back(&(*genJet));
  }
  std::sort(genJets_sorted.begin(), genJets_sorted.end(), isHigherPt);

  const reco::GenJet* genJet_lead    = ( genJets_sorted.size() >= 1 ) ? genJets_sorted[0] : nullptr;
  const reco::GenJet* genJet_sublead = ( genJets_sorted.size() >= 2 ) ? genJets_sorted[1] : nullptr;
  if ( genPtHat_hardscatter_isValid && genJet_lead )
  { 
    fillWithOverFlow2d(histogram_leadGenJetPt_vs_genPtHat_, genPtHat_hardscatter, genJet_lead->pt(), evtWeight);
  }
  if ( genPtHat_hardscatter_isValid && genJet_sublead )
  {
    fillWithOverFlow2d(histogram_subleadGenJetPt_vs_genPtHat_, genPtHat_hardscatter, genJet_sublead->pt(), evtWeight);
  }
  
  for ( std::vector<PileupSummaryInfo>::const_iterator pileupSummaryInfo = pileupSummaryInfos->begin();
        pileupSummaryInfo != pileupSummaryInfos->end(); ++pileupSummaryInfo ) {
    if ( pileupSummaryInfo->getBunchCrossing() == 0 ) // CV: in-time pileup
    {
      const std::vector<float>& pileupSummaryInfo_genPtHat = pileupSummaryInfo->getPU_pT_hats();
      for ( std::vector<float>::const_iterator genPtHat_pileup = pileupSummaryInfo_genPtHat.begin();
            genPtHat_pileup != pileupSummaryInfo_genPtHat.end(); ++genPtHat_pileup ) {
        fillWithOverFlow(histogram_genPtHat_pileup_, *genPtHat_pileup, evtWeight);
      }
    }
  }
}

void GenPtHatAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(GenPtHatAnalyzer);


