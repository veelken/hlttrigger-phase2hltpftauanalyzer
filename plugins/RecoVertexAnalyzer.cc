#include "HLTTrigger/TallinnHLTPFTauAnalyzer/plugins/RecoVertexAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

RecoVertexAnalyzer::RecoVertexAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  src_genVertex_z_ = cfg.getParameter<edm::InputTag>("srcGenVertex_z");
  token_genVertex_z_ = consumes<float>(src_genVertex_z_);

  src_hltVertices_ = cfg.getParameter<edm::InputTag>("srcHLTVertices");
  token_hltVertices_ = consumes<reco::VertexCollection>(src_hltVertices_);

  src_offlineVertices_ = cfg.getParameter<edm::InputTag>("srcOfflineVertices");
  token_offlineVertices_ = consumes<reco::VertexCollection>(src_offlineVertices_);

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

RecoVertexAnalyzer::~RecoVertexAnalyzer()
{
  delete hltVertexPlots_;

  delete offlineVertexPlots_;
}

void RecoVertexAnalyzer::beginJob()
{
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    throw cms::Exception("JetToTauFakeRateAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  dqmStore.setCurrentFolder(Form("%s/%s", dqmDirectory_.data(), "hltVertex"));
  hltVertexPlots_ = new vertexPlotEntryType();
  hltVertexPlots_->bookHistograms(dqmStore); 
  dqmStore.setCurrentFolder(Form("%s/%s", dqmDirectory_.data(), "offlineVertex"));
  offlineVertexPlots_ = new vertexPlotEntryType();
  offlineVertexPlots_->bookHistograms(dqmStore);
}

namespace
{
  std::vector<double> get_recVertex_z(const reco::VertexCollection& recVertices)
  {
    std::vector<double> recVertices_z;
    for ( auto recVertex : recVertices )
    {
      recVertices_z.push_back(recVertex.position().z());
    }
    return recVertices_z;
  }
}

void RecoVertexAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<float> genVertex_z;
  evt.getByToken(token_genVertex_z_, genVertex_z);

  const double evtWeight = 1.;

  edm::Handle<reco::VertexCollection> hltVertices;
  evt.getByToken(token_hltVertices_, hltVertices);
  std::vector<double> hltVertices_z = get_recVertex_z(*hltVertices);
  hltVertexPlots_->fillHistograms(*genVertex_z, hltVertices_z, evtWeight);

  edm::Handle<reco::VertexCollection> offlineVertices;
  evt.getByToken(token_offlineVertices_, offlineVertices);
  std::vector<double> offlineVertices_z = get_recVertex_z(*offlineVertices);
  offlineVertexPlots_->fillHistograms(*genVertex_z, offlineVertices_z, evtWeight);
}

void RecoVertexAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(RecoVertexAnalyzer);
