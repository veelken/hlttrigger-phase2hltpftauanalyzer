#ifndef HLTTrigger_TallinnHLTPFTauAnalyzer_MyPFTauSelector_h
#define HLTTrigger_TallinnHLTPFTauAnalyzer_MyPFTauSelector_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "DataFormats/TauReco/interface/PFTau.h"              // reco::PFTau
#include "DataFormats/TauReco/interface/PFTauFwd.h"           // reco::PFTauCollection
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h" // reco::PFTauDiscriminator

#include <vector>

class MyPFTauSelector : public edm::EDProducer 
{
 public:
  explicit MyPFTauSelector(const edm::ParameterSet& cfg);
  ~MyPFTauSelector();

 private:
  void produce(edm::Event& evt, const edm::EventSetup& es);

  edm::InputTag src_;
  edm::EDGetTokenT<reco::PFTauCollection> token_;
  edm::InputTag src_sumChargedIso_;
  edm::EDGetTokenT<reco::PFTauDiscriminator> token_sumChargedIso_;

  double min_pt_;
  double max_pt_;
  double min_absEta_;
  double max_absEta_;
  double min_leadTrackPt_;
  double max_leadTrackPt_;
  double min_relChargedIso_;
  double max_relChargedIso_;
  double min_absChargedIso_;
  double max_absChargedIso_;

  bool invert_;
};

#endif
