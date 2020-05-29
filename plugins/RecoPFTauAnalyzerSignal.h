#ifndef HLTTrigger_TallinnHLTPFTauAnalyzer_RecoPFTauAnalyzerSignal_h
#define HLTTrigger_TallinnHLTPFTauAnalyzer_RecoPFTauAnalyzerSignal_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/TauReco/interface/PFTau.h"                     // reco::PFTau
#include "DataFormats/TauReco/interface/PFTauFwd.h"                  // reco::PFTauCollection
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"        // reco::PFTauDiscriminator
#include "DataFormats/Math/interface/deltaR.h"                       // reco::deltaR
#include "DataFormats/JetReco/interface/GenJet.h"                    // reco::GenJet
#include "DataFormats/JetReco/interface/GenJetCollection.h"          // reco::GenJetCollection
#include "DataFormats/PatCandidates/interface/Tau.h"                 // pat::Tau, pat::TauCollection
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"              // JetMCTagUtils::genTauDecayMode()
#include "DataFormats/TauReco/interface/PFTau.h"                     // reco::PFTau::kOneProng0PiZero, reco::PFTau::kOneProng1PiZero,...

#include <TH1.h>     // TH1
#include <TH2.h>     // TH2
#include <TString.h> // TString, Form()
#include <TMath.h>   // TMath::Abs(), TMath::Pi()

#include <vector>    // std::vector
#include <string>    // std::string

using namespace dqm::implementation;

class RecoPFTauAnalyzerSignal : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit RecoPFTauAnalyzerSignal(const edm::ParameterSet&);
    
  // destructor
  ~RecoPFTauAnalyzerSignal();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag srcPFTaus_;
  edm::EDGetTokenT<reco::PFTauCollection> tokenPFTaus_;
  edm::InputTag srcPFTauSumChargedIso_;
  edm::EDGetTokenT<reco::PFTauDiscriminator> tokenPFTauSumChargedIso_;
  edm::InputTag srcDenominator_;
  enum { kGen, kOffline };
  int typeDenominator_;
  edm::EDGetTokenT<reco::GenJetCollection> tokenDenominator_gen_;
  edm::EDGetTokenT<pat::TauCollection> tokenDenominator_offline_;

  std::string dqmDirectory_;

  struct efficiencyPlotEntryType
  {
    efficiencyPlotEntryType(double min_pt, double max_pt, double min_absEta, double max_absEta, const std::string& decayMode, 
			    double ptThreshold, double max_relChargedIso, double max_absChargedIso)
      : me_pt_numerator_(nullptr)
      , histogram_pt_numerator_(nullptr)
      , me_pt_denominator_(nullptr)
      , histogram_pt_denominator_(nullptr)
      , me_eta_numerator_(nullptr)
      , histogram_eta_numerator_(nullptr)
      , me_eta_denominator_(nullptr)
      , histogram_eta_denominator_(nullptr)
      , me_phi_numerator_(nullptr)
      , histogram_phi_numerator_(nullptr)
      , me_phi_denominator_(nullptr)
      , histogram_phi_denominator_(nullptr)
      , me_minDeltaR_numerator_(nullptr)
      , histogram_minDeltaR_numerator_(nullptr)
      , me_minDeltaR_denominator_(nullptr)
      , histogram_minDeltaR_denominator_(nullptr)
      , me_pt_vs_absEta_numerator_(nullptr)
      , histogram_pt_vs_absEta_numerator_(nullptr)
      , me_pt_vs_absEta_denominator_(nullptr)
      , histogram_pt_vs_absEta_denominator_(nullptr)
      , min_pt_(min_pt)
      , max_pt_(max_pt)
      , min_absEta_(min_absEta)
      , max_absEta_(max_absEta)
      , decayMode_(decayMode)
      , ptThreshold_(ptThreshold)
      , max_relChargedIso_(max_relChargedIso)
      , max_absChargedIso_(max_absChargedIso)
      , dRmatch_(0.3)
    {}
    ~efficiencyPlotEntryType() 
    {}
    void bookHistograms(DQMStore& dqmStore)
    {
      TString histogramName_suffix = decayMode_.data(); 
      if      ( min_absEta_ >= 0. && max_absEta_ > 0. ) histogramName_suffix.Append(Form("_absEta%1.2fto%1.2f", min_absEta_, max_absEta_));
      else if ( min_absEta_ >= 0.                     ) histogramName_suffix.Append(Form("_absEtaGt%1.2f", min_absEta_));
      else if (                      max_absEta_ > 0. ) histogramName_suffix.Append(Form("_absEtaLt%1.2f", max_absEta_));
      if ( ptThreshold_ > 0. ) histogramName_suffix.Append(Form("_ptGt%1.0f", ptThreshold_));
      if ( max_relChargedIso_ > 0. ) histogramName_suffix.Append(Form("_relChargedIsoLt%1.2f", max_relChargedIso_));
      if ( max_absChargedIso_ > 0. ) histogramName_suffix.Append(Form("_absChargedIsoLt%1.2f", max_absChargedIso_));
      histogramName_suffix = histogramName_suffix.ReplaceAll(".", "p");

      TString histogramName_pt_numerator = Form("effL1PFTau_vs_pt_numerator_%s", histogramName_suffix.Data());
      me_pt_numerator_ = dqmStore.book1D(histogramName_pt_numerator.Data(), histogramName_pt_numerator.Data(), 16, 20., 100.);
      histogram_pt_numerator_ = me_pt_numerator_->getTH1();
      assert(histogram_pt_numerator_);
      TString histogramName_pt_denominator = Form("effL1PFTau_vs_pt_denominator_%s", histogramName_suffix.Data());
      me_pt_denominator_ = dqmStore.book1D(histogramName_pt_denominator.Data(), histogramName_pt_denominator.Data(), 16, 20., 100.);
      histogram_pt_denominator_ = me_pt_denominator_->getTH1();
      assert(histogram_pt_denominator_);

      TString histogramName_eta_numerator = Form("effL1PFTau_vs_eta_numerator_%s", histogramName_suffix.Data());
      me_eta_numerator_ = dqmStore.book1D(histogramName_eta_numerator.Data(), histogramName_eta_numerator.Data(), 30, -3., +3.);
      histogram_eta_numerator_ = me_eta_numerator_->getTH1();
      assert(histogram_eta_numerator_);
      TString histogramName_eta_denominator = Form("effL1PFTau_vs_eta_denominator_%s", histogramName_suffix.Data());
      me_eta_denominator_ = dqmStore.book1D(histogramName_eta_denominator.Data(), histogramName_eta_denominator.Data(), 30, -3., +3.);
      histogram_eta_denominator_ = me_eta_denominator_->getTH1();
      assert(histogram_eta_denominator_);

      TString histogramName_phi_numerator = Form("effL1PFTau_vs_phi_numerator_%s", histogramName_suffix.Data());
      me_phi_numerator_ = dqmStore.book1D(histogramName_phi_numerator.Data(), histogramName_phi_numerator.Data(), 18, -TMath::Pi(), +TMath::Pi());
      histogram_phi_numerator_ = me_phi_numerator_->getTH1();
      assert(histogram_phi_numerator_);
      TString histogramName_phi_denominator = Form("effL1PFTau_vs_phi_denominator_%s", histogramName_suffix.Data());
      me_phi_denominator_ = dqmStore.book1D(histogramName_phi_denominator.Data(), histogramName_phi_denominator.Data(), 18, -TMath::Pi(), +TMath::Pi());
      histogram_phi_denominator_ = me_phi_denominator_->getTH1();
      assert(histogram_phi_denominator_);

      TString histogramName_minDeltaR_numerator = Form("effL1PFTau_vs_minDeltaR_numerator_%s", histogramName_suffix.Data());
      me_minDeltaR_numerator_ = dqmStore.book1D(histogramName_minDeltaR_numerator.Data(), histogramName_minDeltaR_numerator.Data(), 50, 0., 5.); 
      histogram_minDeltaR_numerator_ = me_minDeltaR_numerator_->getTH1();
      assert(histogram_minDeltaR_numerator_);
      TString histogramName_minDeltaR_denominator = Form("effL1PFTau_vs_minDeltaR_denominator_%s", histogramName_suffix.Data());
      me_minDeltaR_denominator_ = dqmStore.book1D(histogramName_minDeltaR_denominator.Data(), histogramName_minDeltaR_denominator.Data(), 50, 0., 5.); 
      histogram_minDeltaR_denominator_ = me_minDeltaR_denominator_->getTH1();
      assert(histogram_minDeltaR_denominator_);

      const int numBins_absEta = 6;
      float binning_absEta[numBins_absEta + 1] = { 
        0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4
      };

      const int numBins_pt = 8;
      float binning_pt[numBins_pt + 1] = { 
        20., 25., 30., 35., 40., 50., 60., 80., 100.
      };

      TString histogramName_pt_vs_absEta_numerator = Form("effL1PFTau_vs_pt_vs_absEta_numerator_%s", histogramName_suffix.Data());
      me_pt_vs_absEta_numerator_ = dqmStore.book2D(histogramName_pt_vs_absEta_numerator.Data(), histogramName_pt_vs_absEta_numerator.Data(), numBins_absEta, binning_absEta, numBins_pt, binning_pt);
      histogram_pt_vs_absEta_numerator_ = dynamic_cast<TH2*>(me_pt_vs_absEta_numerator_->getTH1());
      assert(histogram_pt_vs_absEta_numerator_);
      TString histogramName_pt_vs_absEta_denominator = Form("effL1PFTau_vs_pt_vs_absEta_denominator_%s", histogramName_suffix.Data());
      me_pt_vs_absEta_denominator_ = dqmStore.book2D(histogramName_pt_vs_absEta_denominator.Data(), histogramName_pt_vs_absEta_denominator.Data(), numBins_absEta, binning_absEta, numBins_pt, binning_pt);
      histogram_pt_vs_absEta_denominator_ = dynamic_cast<TH2*>(me_pt_vs_absEta_denominator_->getTH1());
      assert(histogram_pt_vs_absEta_denominator_);
    }
    void fillHistograms(const edm::Handle<reco::PFTauCollection>& pfTaus, const reco::PFTauDiscriminator& pfTauSumChargedIso, 
                        const reco::GenJetCollection& denominatorTaus, double evtWeight)
    {
      double minDeltaR = 1.e+3;
      for ( reco::GenJetCollection::const_iterator denominatorTau1 = denominatorTaus.begin();
	    denominatorTau1 != denominatorTaus.end(); ++denominatorTau1 ) {
	for ( reco::GenJetCollection::const_iterator denominatorTau2 = denominatorTau1 + 1;
	      denominatorTau2 != denominatorTaus.end(); ++denominatorTau2 ) {
	  double dR = deltaR(denominatorTau1->eta(), denominatorTau1->phi(), denominatorTau2->eta(), denominatorTau2->phi());
	  if ( dR > 0. && dR < minDeltaR ) minDeltaR = dR;
        }
      }

      for ( auto denominatorTau : denominatorTaus )
      {
	std::string denominatorTau_decayMode = JetMCTagUtils::genTauDecayMode(denominatorTau);
	if ( decayMode_ != "all" && denominatorTau_decayMode != decayMode_ ) continue;

        bool isMatched = false;

        size_t numPFTaus = pfTaus->size();
        for ( size_t idxPFTau = 0; idxPFTau < numPFTaus; ++idxPFTau ) 
        {  
          reco::PFTauRef pfTauRef(pfTaus, idxPFTau);
          double sumChargedIso = pfTauSumChargedIso[pfTauRef];
  	  if ( (ptThreshold_       < 0. || pfTauRef->pt() >=  ptThreshold_                      ) &&
               (                           pfTauRef->leadPFChargedHadrCand().isNonnull()        &&   
                                           pfTauRef->leadPFChargedHadrCand()->bestTrack()       ) && 
               (max_relChargedIso_ < 0. || sumChargedIso  <= (max_relChargedIso_*pfTauRef->pt())) &&
	       (max_absChargedIso_ < 0. || sumChargedIso  <=  max_absChargedIso_                ) )
	  {
	    double dR = reco::deltaR(denominatorTau.eta(), denominatorTau.phi(), pfTauRef->eta(), pfTauRef->phi());
	    if ( dR < dRmatch_ ) isMatched = true;
	  }
	}

	double denominatorTau_absEta = TMath::Abs(denominatorTau.eta());
	if ( denominatorTau_absEta > min_absEta_ && denominatorTau_absEta < max_absEta_ ) 
	{
	  histogram_pt_denominator_->Fill(denominatorTau.pt(), evtWeight);
	  if ( isMatched ) histogram_pt_numerator_->Fill(denominatorTau.pt(), evtWeight);
	}
	if ( denominatorTau.pt() > min_pt_ && denominatorTau.pt() < max_pt_ ) 
	{
	  histogram_eta_denominator_->Fill(denominatorTau.eta(), evtWeight);
	  if ( isMatched ) histogram_eta_numerator_->Fill(denominatorTau.eta(), evtWeight);
	}
	if ( denominatorTau.pt() > min_pt_ && denominatorTau.pt() < max_pt_ && denominatorTau_absEta > min_absEta_ && denominatorTau_absEta < max_absEta_ )
	{
	  histogram_phi_denominator_->Fill(denominatorTau.phi(), evtWeight);
	  if ( isMatched ) histogram_phi_numerator_->Fill(denominatorTau.phi(), evtWeight);
	}
	if ( denominatorTau.pt() > min_pt_ && denominatorTau.pt() < max_pt_ && denominatorTau_absEta > min_absEta_ && denominatorTau_absEta < max_absEta_ )
	{
	  histogram_minDeltaR_denominator_->Fill(minDeltaR, evtWeight);
	  if ( isMatched ) histogram_minDeltaR_numerator_->Fill(minDeltaR, evtWeight);
	}
	histogram_pt_vs_absEta_denominator_->Fill(denominatorTau_absEta, denominatorTau.pt(), evtWeight);
	if ( isMatched ) histogram_pt_vs_absEta_numerator_->Fill(denominatorTau_absEta, denominatorTau.pt(), evtWeight);
      }
    }
    void fillHistograms(const edm::Handle<reco::PFTauCollection>& pfTaus, const reco::PFTauDiscriminator& pfTauSumChargedIso,
                        const pat::TauCollection& denominatorTaus, double evtWeight)
    {
      double minDeltaR = 1.e+3;
      for ( pat::TauCollection::const_iterator denominatorTau1 = denominatorTaus.begin();
	    denominatorTau1 != denominatorTaus.end(); ++denominatorTau1 ) {
	for ( pat::TauCollection::const_iterator denominatorTau2 = denominatorTau1 + 1;
	      denominatorTau2 != denominatorTaus.end(); ++denominatorTau2 ) {
	  double dR = deltaR(denominatorTau1->eta(), denominatorTau1->phi(), denominatorTau2->eta(), denominatorTau2->phi());
	  if ( dR > 0. && dR < minDeltaR ) minDeltaR = dR;
        }
      }

      for ( auto denominatorTau : denominatorTaus )
      {
	if ( (decayMode_ == "oneProng0Pi0"   && denominatorTau.decayMode() != reco::PFTau::kOneProng0PiZero)   ||
	     (decayMode_ == "oneProng1Pi0"   && denominatorTau.decayMode() != reco::PFTau::kOneProng1PiZero)   ||
	     (decayMode_ == "oneProng1Pi0"   && denominatorTau.decayMode() != reco::PFTau::kOneProng2PiZero)   ||
	     (decayMode_ == "threeProng0Pi0" && denominatorTau.decayMode() != reco::PFTau::kThreeProng0PiZero) ||
	     (decayMode_ == "threeProng1Pi0" && denominatorTau.decayMode() != reco::PFTau::kThreeProng1PiZero) ) continue;	       

        bool isMatched = false;

        size_t numPFTaus = pfTaus->size();
        for ( size_t idxPFTau = 0; idxPFTau < numPFTaus; ++idxPFTau ) 
        {  
          reco::PFTauRef pfTauRef(pfTaus, idxPFTau);
          double sumChargedIso = pfTauSumChargedIso[pfTauRef];
  	  if ( (ptThreshold_       < 0. || pfTauRef->pt() >=  ptThreshold_                      ) &&
               (                           pfTauRef->leadPFChargedHadrCand().isNonnull()        &&   
                                           pfTauRef->leadPFChargedHadrCand()->bestTrack()       ) && 
               (max_relChargedIso_ < 0. || sumChargedIso  <= (max_relChargedIso_*pfTauRef->pt())) &&
	       (max_absChargedIso_ < 0. || sumChargedIso  <=  max_absChargedIso_                ) )
	  {
	    double dR = reco::deltaR(denominatorTau.eta(), denominatorTau.phi(), pfTauRef->eta(), pfTauRef->phi());
	    if ( dR < dRmatch_ ) isMatched = true;
	  }
	}

	double denominatorTau_absEta = TMath::Abs(denominatorTau.eta());
	if ( denominatorTau_absEta > min_absEta_ && denominatorTau_absEta < max_absEta_ ) 
	{
	  histogram_pt_denominator_->Fill(denominatorTau.pt(), evtWeight);
	  if ( isMatched ) histogram_pt_numerator_->Fill(denominatorTau.pt(), evtWeight);
	}
	if ( denominatorTau.pt() > min_pt_ && denominatorTau.pt() < max_pt_ ) 
	{
	  histogram_eta_denominator_->Fill(denominatorTau.eta(), evtWeight);
	  if ( isMatched ) histogram_eta_numerator_->Fill(denominatorTau.eta(), evtWeight);
	}
	if ( denominatorTau.pt() > min_pt_ && denominatorTau.pt() < max_pt_ && denominatorTau_absEta > min_absEta_ && denominatorTau_absEta < max_absEta_ )
	{
	  histogram_phi_denominator_->Fill(denominatorTau.phi(), evtWeight);
	  if ( isMatched ) histogram_phi_numerator_->Fill(denominatorTau.phi(), evtWeight);
	}
	if ( denominatorTau.pt() > min_pt_ && denominatorTau.pt() < max_pt_ && denominatorTau_absEta > min_absEta_ && denominatorTau_absEta < max_absEta_ )
	{
	  histogram_minDeltaR_denominator_->Fill(minDeltaR, evtWeight);
	  if ( isMatched ) histogram_minDeltaR_numerator_->Fill(minDeltaR, evtWeight);
	}
	histogram_pt_vs_absEta_denominator_->Fill(denominatorTau_absEta, denominatorTau.pt(), evtWeight);
	if ( isMatched ) histogram_pt_vs_absEta_numerator_->Fill(denominatorTau_absEta, denominatorTau.pt(), evtWeight);
      }
    }
    MonitorElement* me_pt_numerator_;
    TH1* histogram_pt_numerator_;
    MonitorElement* me_pt_denominator_;
    TH1* histogram_pt_denominator_;
    MonitorElement* me_eta_numerator_;
    TH1* histogram_eta_numerator_;
    MonitorElement* me_eta_denominator_;
    TH1* histogram_eta_denominator_;
    MonitorElement* me_phi_numerator_;
    TH1* histogram_phi_numerator_;
    MonitorElement* me_phi_denominator_;
    TH1* histogram_phi_denominator_;
    MonitorElement* me_minDeltaR_numerator_;
    TH1* histogram_minDeltaR_numerator_;
    MonitorElement* me_minDeltaR_denominator_;
    TH1* histogram_minDeltaR_denominator_;
    MonitorElement* me_pt_vs_absEta_numerator_;
    TH2* histogram_pt_vs_absEta_numerator_;
    MonitorElement* me_pt_vs_absEta_denominator_;
    TH2* histogram_pt_vs_absEta_denominator_;
    // cuts applied to offline and generator-level taus in denominator 
    double min_pt_;
    double max_pt_;
    double min_absEta_;
    double max_absEta_;
    std::string decayMode_;
    // cuts applied to L1 trigger taus in numerator 
    double ptThreshold_;
    double max_relChargedIso_;
    double max_absChargedIso_;
    // matching between offline or generator-level taus in denominator and L1 trigger taus in numerator
    double dRmatch_; 
  };
  std::vector<efficiencyPlotEntryType*> efficiencyPlots_;
};

#endif   
