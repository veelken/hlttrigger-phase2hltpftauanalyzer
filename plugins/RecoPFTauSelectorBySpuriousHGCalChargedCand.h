#ifndef HLTrigger_TallinnHLTPFTauAnalyzer_RecoPFTauSelectorBySpuriousHGCalChargedCand_h
#define HLTrigger_TallinnHLTPFTauAnalyzer_RecoPFTauSelectorBySpuriousHGCalChargedCand_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/TauReco/interface/PFTau.h"    // reco::PFTau
#include "DataFormats/TauReco/interface/PFTauFwd.h" // reco::PFTauCollection

class RecoPFTauSelectorBySpuriousHGCalChargedCandImp
{
 public:
  typedef reco::PFTauCollection collection;

  explicit RecoPFTauSelectorBySpuriousHGCalChargedCandImp(const edm::ParameterSet& cfg, edm::ConsumesCollector&& iC);

  typename std::vector<const reco::PFTau*>::const_iterator begin() const { return selected_.begin(); }
  typename std::vector<const reco::PFTau*>::const_iterator end() const { return selected_.end(); }

  void select(const edm::Handle<reco::PFTauCollection>& pfTaus, const edm::Event& evt, const edm::EventSetup& es);

  size_t size() const { return selected_.size(); }

 private:
  std::vector<const reco::PFTau*> selected_;

  double min_trackPt_div_pfCandPt_;
  double max_trackPt_div_pfCandPt_;
};

#endif
