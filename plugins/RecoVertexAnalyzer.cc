#include "HLTrigger/TallinnHLTPFTauAnalyzer/plugins/RecoVertexAnalyzer.h"

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

  src_recVertices_ = cfg.getParameter<edm::InputTag>("srcRecVertices");
  token_recVertices_ = consumes<reco::VertexCollection>(src_recVertices_);

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

RecoVertexAnalyzer::~RecoVertexAnalyzer()
{
  delete recVertexPlots_;
}

void RecoVertexAnalyzer::beginJob()
{
  if ( !edm::Service<dqm::legacy::DQMStore>().isAvailable() ) {
    throw cms::Exception("JetToTauFakeRateAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<dqm::legacy::DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory_.data());
  recVertexPlots_ = new vertexPlotEntryType();
  recVertexPlots_->bookHistograms(dqmStore); 
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

  edm::Handle<reco::VertexCollection> recVertices;
  evt.getByToken(token_recVertices_, recVertices);
  std::vector<double> recVertices_z = get_recVertex_z(*recVertices);
  recVertexPlots_->fillHistograms(*genVertex_z, recVertices_z, evtWeight);
}

void RecoVertexAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(RecoVertexAnalyzer);
