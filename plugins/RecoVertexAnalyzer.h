#ifndef HLTTrigger_TallinnHLTPFTauAnalyzer_RecoVertexAnalyzer_h
#define HLTTrigger_TallinnHLTPFTauAnalyzer_RecoVertexAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/VertexReco/interface/Vertex.h"    // reco::Vertex
#include "DataFormats/VertexReco/interface/VertexFwd.h" // reco::VertexCollection

#include <TH1.h>
#include <TH2.h>
#include <TString.h> // TString, Form()
#include <TMath.h> // TMath::Abs(), TMath::Nint()

#include <vector>
#include <string>
#include <algorithm> // std::sort

using namespace dqm::implementation;

class RecoVertexAnalyzer : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit RecoVertexAnalyzer(const edm::ParameterSet&);
    
  // destructor
  ~RecoVertexAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag src_genVertex_z_;
  edm::EDGetTokenT<float> token_genVertex_z_;

  edm::InputTag src_hltVertices_;
  edm::EDGetTokenT<reco::VertexCollection> token_hltVertices_;

  edm::InputTag src_offlineVertices_;
  edm::EDGetTokenT<reco::VertexCollection> token_offlineVertices_;

  std::string dqmDirectory_;

  struct vertexPlotEntryType
  {
    vertexPlotEntryType()
      : me_genVertex_z_(nullptr)
      , histogram_genVertex_z_(nullptr)
      , me_recVertex0_z_(nullptr)
      , histogram_recVertex0_z_(nullptr)
      , me_delta_z_(nullptr)
      , histogram_delta_z_(nullptr)
      , me_idxRecVertex_bestMatch_(nullptr)
      , histogram_idxRecVertex_bestMatch_(nullptr)
      , me_delta_z_bestMatch_(nullptr)
      , histogram_delta_z_bestMatch_(nullptr)
      , me_numRecVertices_(nullptr)
      , histogram_numRecVertices_(nullptr)
    {}
    ~vertexPlotEntryType() 
    {}
    void bookHistograms(DQMStore& dqmStore)
    {
      me_genVertex_z_ = dqmStore.book1D("genVertex_z", "genVertex_z", 100, -25., +25.);
      histogram_genVertex_z_ = me_genVertex_z_->getTH1();
      assert(histogram_genVertex_z_);

      me_recVertex0_z_ = dqmStore.book1D("recVertex0_z", "recVertex0_z", 100, -25., +25.);
      histogram_recVertex0_z_ = me_recVertex0_z_->getTH1();
      assert(histogram_recVertex0_z_);

      me_delta_z_ = dqmStore.book1D("delta_z", "delta_z", 500, -25., +25.);
      histogram_delta_z_ = me_delta_z_->getTH1();
      assert(histogram_delta_z_);

      me_idxRecVertex_bestMatch_ = dqmStore.book1D("idxRecVertex_bestMatch", "idxRecVertex_bestMatch", 50, -0.5, +49.5);
      histogram_idxRecVertex_bestMatch_ = me_idxRecVertex_bestMatch_->getTH1();
      assert(histogram_idxRecVertex_bestMatch_);

      me_delta_z_bestMatch_ = dqmStore.book1D("delta_z_bestMatch", "delta_z_bestMatch", 500, -25., +25.);
      histogram_delta_z_bestMatch_ = me_delta_z_bestMatch_->getTH1();
      assert(histogram_delta_z_bestMatch_);

      me_numRecVertices_ = dqmStore.book1D("numRecVertices", "numRecVertices", 100, -0.5, +99.5);
      histogram_numRecVertices_ = me_numRecVertices_->getTH1();
      assert(histogram_numRecVertices_);
    }
    void fillHistograms(double genVertex_z, const std::vector<double>& recVertices_z, double evtWeight)
    {
      histogram_genVertex_z_->Fill(genVertex_z, evtWeight);

      size_t numRecVertices = recVertices_z.size();
      if ( numRecVertices > 0 ) 
      {
	histogram_recVertex0_z_->Fill(recVertices_z[0], evtWeight);
	histogram_delta_z_->Fill(recVertices_z[0] - genVertex_z, evtWeight);

	int idxRecVertex_bestMatch = -1;
	double min_delta_z = 1.e+3;
	for ( size_t idxRecVertex = 0; idxRecVertex < numRecVertices; ++idxRecVertex ) 
	{
	  double delta_z = recVertices_z[idxRecVertex] - genVertex_z;
	  if ( TMath::Abs(delta_z) < TMath::Abs(min_delta_z) ) 
	  {
	    idxRecVertex_bestMatch = idxRecVertex;
	    min_delta_z = delta_z;
	  }
	}
	assert(idxRecVertex_bestMatch >= 0 && idxRecVertex_bestMatch < (int)numRecVertices);
	histogram_idxRecVertex_bestMatch_->Fill(idxRecVertex_bestMatch, evtWeight);
	histogram_delta_z_bestMatch_->Fill(min_delta_z, evtWeight);	
      }

      histogram_numRecVertices_->Fill(numRecVertices, evtWeight);
    }
    MonitorElement* me_genVertex_z_;
    TH1* histogram_genVertex_z_;
    MonitorElement* me_recVertex0_z_;
    TH1* histogram_recVertex0_z_;
    MonitorElement* me_delta_z_;
    TH1* histogram_delta_z_;
    MonitorElement* me_idxRecVertex_bestMatch_;
    TH1* histogram_idxRecVertex_bestMatch_;
    MonitorElement* me_delta_z_bestMatch_;
    TH1* histogram_delta_z_bestMatch_;
    MonitorElement* me_numRecVertices_;
    TH1* histogram_numRecVertices_;    
  };

  vertexPlotEntryType* hltVertexPlots_;

  vertexPlotEntryType* offlineVertexPlots_;
};

#endif   
