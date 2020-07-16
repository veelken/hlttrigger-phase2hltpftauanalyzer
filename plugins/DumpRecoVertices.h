#ifndef HLTrigger_TallinnHLTPFTauAnalyzer_DumpRecoVertices_h
#define HLTrigger_TallinnHLTPFTauAnalyzer_DumpRecoVertices_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/VertexReco/interface/Vertex.h"    // reco::Vertex
#include "DataFormats/VertexReco/interface/VertexFwd.h" // reco::VertexCollection

#include <vector>
#include <string>

class DumpRecoVertices : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit DumpRecoVertices(const edm::ParameterSet&);
    
  // destructor
  ~DumpRecoVertices();
    
 private:
  void analyze(const edm::Event&, const edm::EventSetup&);

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<reco::VertexCollection> token_;
};

#endif   
