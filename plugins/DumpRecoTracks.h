#ifndef HLTrigger_TallinnHLTPFTauAnalyzer_DumpRecoTracks_h
#define HLTrigger_TallinnHLTPFTauAnalyzer_DumpRecoTracks_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"    // reco::Track            
#include "DataFormats/TrackReco/interface/TrackFwd.h" // reco::TrackCollection

#include <vector>
#include <string>

class DumpRecoTracks : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit DumpRecoTracks(const edm::ParameterSet&);
    
  // destructor
  ~DumpRecoTracks();
    
 private:
  void analyze(const edm::Event&, const edm::EventSetup&);

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<reco::TrackCollection> token_;
};

#endif   
