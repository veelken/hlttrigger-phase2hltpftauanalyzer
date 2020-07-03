#include "HLTTrigger/TallinnHLTPFTauAnalyzer/plugins/RecoPFTauSelectorBySpuriousHGCalChargedCand.h"

RecoPFTauSelectorBySpuriousHGCalChargedCandImp::RecoPFTauSelectorBySpuriousHGCalChargedCandImp(const edm::ParameterSet& cfg, edm::ConsumesCollector&& iC)
{
  min_trackPt_div_pfCandPt_ = cfg.getParameter<double>("min_trackPt_div_pfCandPt");
  max_trackPt_div_pfCandPt_ = cfg.getParameter<double>("max_trackPt_div_pfCandPt");
}

void 
RecoPFTauSelectorBySpuriousHGCalChargedCandImp::select(const edm::Handle<reco::PFTauCollection>& pfTaus, const edm::Event& evt, const edm::EventSetup& es)
{
  selected_.clear();
  //std::cout << "#pfTaus = " << pfTaus->size() << std::endl;

  for ( reco::PFTauCollection::const_iterator pfTau = pfTaus->begin(); 
        pfTau != pfTaus->end(); ++pfTau ) 
  {
    bool isSelected = false;
    const std::vector<reco::PFCandidatePtr>& pfCands = pfTau->signalPFChargedHadrCands();
    for ( std::vector<reco::PFCandidatePtr>::const_iterator pfCand = pfCands.begin();
          pfCand != pfCands.end(); ++pfCand ) 
    {
      if ( (*pfCand)->bestTrack() )
      {
        double trackPt_div_pfCandPt = (*pfCand)->bestTrack()->pt()/(*pfCand)->pt();
        if ( (min_trackPt_div_pfCandPt_ < 0. || trackPt_div_pfCandPt >= min_trackPt_div_pfCandPt_) &&
             (max_trackPt_div_pfCandPt_ < 0. || trackPt_div_pfCandPt <= max_trackPt_div_pfCandPt_) ) 
        {
          isSelected = true;
        } 
      }
    }
    if ( isSelected ) 
    {
      selected_.push_back(&(*pfTau));
    }
  }
}

#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

typedef ObjectSelector<RecoPFTauSelectorBySpuriousHGCalChargedCandImp> RecoPFTauSelectorBySpuriousHGCalChargedCand;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(RecoPFTauSelectorBySpuriousHGCalChargedCand);
