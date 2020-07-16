#include "HLTrigger/TallinnHLTPFTauAnalyzer/plugins/MyPFTauSelector.h"

#include "FWCore/Utilities/interface/InputTag.h"

MyPFTauSelector::MyPFTauSelector(const edm::ParameterSet& cfg) 
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<reco::PFTauCollection>(src_);
  src_sumChargedIso_ = cfg.getParameter<edm::InputTag>("src_sumChargedIso");
  token_sumChargedIso_ = consumes<reco::PFTauDiscriminator>(src_sumChargedIso_);

  min_pt_ = cfg.getParameter<double>("min_pt");
  max_pt_ = cfg.getParameter<double>("max_pt");
  min_absEta_ = cfg.getParameter<double>("min_absEta");
  max_absEta_ = cfg.getParameter<double>("max_absEta");
  min_leadTrackPt_ = cfg.getParameter<double>("min_leadTrackPt");
  max_leadTrackPt_ = cfg.getParameter<double>("max_leadTrackPt");
  min_relChargedIso_ = cfg.getParameter<double>("min_relChargedIso");
  max_relChargedIso_ = cfg.getParameter<double>("max_relChargedIso");
  min_absChargedIso_ = cfg.getParameter<double>("min_absChargedIso");
  max_absChargedIso_ = cfg.getParameter<double>("max_absChargedIso");

  invert_ = cfg.getParameter<bool>("invert");

  produces<reco::PFTauCollection>("");
}

MyPFTauSelector::~MyPFTauSelector()
{}

void MyPFTauSelector::produce(edm::Event& evt, const edm::EventSetup& es)
{
  std::unique_ptr<reco::PFTauCollection> pfTaus_selected(new reco::PFTauCollection());

  edm::Handle<reco::PFTauCollection> pfTaus;
  evt.getByToken(token_, pfTaus);
  
  edm::Handle<reco::PFTauDiscriminator> pfTauSumChargedIso;
  evt.getByToken(token_sumChargedIso_, pfTauSumChargedIso);

  size_t numPFTaus = pfTaus->size();
  for ( size_t idxPFTau = 0; idxPFTau < numPFTaus; ++idxPFTau ) 
  {  
    reco::PFTauRef pfTau(pfTaus, idxPFTau);
    double pfTau_absEta = std::fabs(pfTau->eta());
    double sumChargedIso = (*pfTauSumChargedIso)[pfTau];
    bool isSelected = false;
    if ( (min_pt_            < 0. || pfTau->pt()                                       >=  min_pt_                        ) &&
         (max_pt_            < 0. || pfTau->pt()                                       <=  max_pt_                        ) &&
         (min_absEta_        < 0. || pfTau_absEta                                      >=  min_absEta_                    ) &&
         (max_absEta_        < 0. || pfTau_absEta                                      <=  max_absEta_                    ) &&
         (                           pfTau->leadPFChargedHadrCand().isNonnull()                                           &&   
                                     pfTau->leadPFChargedHadrCand()->bestTrack()                                          ) && 
         (min_leadTrackPt_   < 0. || pfTau->leadPFChargedHadrCand()->bestTrack()->pt() >=  min_leadTrackPt_               ) &&
         (max_leadTrackPt_   < 0. || pfTau->leadPFChargedHadrCand()->bestTrack()->pt() <=  max_leadTrackPt_               ) &&
         (min_relChargedIso_ < 0. || sumChargedIso                                     >= (min_relChargedIso_*pfTau->pt())) &&
         (max_relChargedIso_ < 0. || sumChargedIso                                     <= (max_relChargedIso_*pfTau->pt())) &&
         (min_absChargedIso_ < 0. || sumChargedIso                                     >=  min_absChargedIso_             ) &&
	 (max_absChargedIso_ < 0. || sumChargedIso                                     <=  max_absChargedIso_             ) )
    {
      isSelected = true;
    }
    if ( (isSelected && !invert_) || (!isSelected && invert_) ) 
    {
      pfTaus_selected->push_back(*pfTau);
    }
  }

  evt.put(std::move(pfTaus_selected));
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(MyPFTauSelector);
