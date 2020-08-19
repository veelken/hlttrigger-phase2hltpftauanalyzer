#ifndef HLTrigger_TallinnHLTPFTauAnalyzer_BaseTauAnalyzerSignal_h
#define HLTrigger_TallinnHLTPFTauAnalyzer_BaseTauAnalyzerSignal_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/Math/interface/deltaR.h"                       // reco::deltaR
#include "DataFormats/JetReco/interface/GenJet.h"                    // reco::GenJet
#include "DataFormats/JetReco/interface/GenJetCollection.h"          // reco::GenJetCollection
#include "DataFormats/PatCandidates/interface/Tau.h"                 // pat::Tau, pat::TauCollection
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"              // JetMCTagUtils::genTauDecayMode()
#include "DataFormats/TauReco/interface/PFTau.h"                     // reco::PFTau::kOneProng0PiZero, reco::PFTau::kOneProng1PiZero,...
#include "HLTrigger/TallinnHLTPFTauAnalyzer/interface/BaseTau.h"

#include <TH1.h>     // TH1
#include <TH2.h>     // TH2
#include <TString.h> // TString, Form()
#include <TMath.h>   // TMath::Abs(), TMath::Pi()

#include <vector>    // std::vector
#include <string>    // std::string

using namespace dqm::implementation;

namespace
{
  TString getHistogramName(const TString& histogramName_old, const std::string& observableName, double observable_min_value, double observable_max_value, int precision)
  {
    TString histogramName_new = histogramName_old;
    std::string number_format = Form("1.%i", precision);
    if( observable_min_value >= 0. && observable_max_value > 0. ) 
    {
      std::string histogramName_format = std::string("_%s%") + number_format + "fto%" + number_format + "f";
      histogramName_new.Append(Form(histogramName_format.data(), observableName.data(), observable_min_value, observable_max_value));
    } 
    else if ( observable_min_value >= 0. ) 
    {
      std::string histogramName_format = std::string("_%sGt%") + number_format + "f";
      histogramName_new.Append(Form(histogramName_format.data(), observableName.data(), observable_min_value));
    }
    else if ( observable_max_value > 0. ) 
    {
      std::string histogramName_format = std::string("_%sLt%") + number_format + "f";
      histogramName_new.Append(Form(histogramName_format.data(), observableName.data(), observable_max_value));
    }
    return histogramName_new;
  }
}

class BaseTauAnalyzerSignal : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit BaseTauAnalyzerSignal(const edm::ParameterSet&);
    
  // destructor
  ~BaseTauAnalyzerSignal();
    
 protected:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  virtual std::vector<BaseTau> buildBaseTaus(const edm::Event&, const edm::EventSetup&) = 0;

  std::string moduleLabel_;

  edm::InputTag srcDenominator_;
  enum { kGen, kOffline };
  int typeDenominator_;
  edm::EDGetTokenT<reco::GenJetCollection> tokenDenominator_gen_;
  edm::EDGetTokenT<pat::TauCollection> tokenDenominator_offline_;

  double min_pt_denominator_;
  double max_pt_denominator_;
  std::vector<double> min_ptValues_numerator_;
  std::vector<double> max_ptValues_numerator_;
  size_t num_ptValues_numerator_;
  std::vector<double> min_absEtaValues_;
  std::vector<double> max_absEtaValues_;
  size_t num_absEtaValues_;
  std::vector<std::string> decayModes_;
  std::vector<double> min_leadTrackPtValues_;
  std::vector<double> max_leadTrackPtValues_;
  size_t num_leadTrackPtValues_;
  std::vector<double> min_relDiscriminatorValues_;
  std::vector<double> max_relDiscriminatorValues_;
  size_t num_relDiscriminatorValues_;
  std::vector<double> min_absDiscriminatorValues_;
  std::vector<double> max_absDiscriminatorValues_;
  size_t num_absDiscriminatorValues_;

  edm::InputTag src_evtWeight_;
  edm::EDGetTokenT<double> token_evtWeight_;

  std::string dqmDirectory_;

 private:
  struct efficiencyPlotEntryType
  {
    efficiencyPlotEntryType(double min_pt_numerator, double max_pt_numerator, double min_pt_denominator, double max_pt_denominator, 
                            double min_absEta, double max_absEta, 
                            const std::string& decayMode, 
                            double min_leadTrackPt, double max_leadTrackPt, 
                            double min_relDiscriminator, double max_relDiscriminator, double min_absDiscriminator, double max_absDiscriminator)
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
      , min_pt_numerator_(min_pt_numerator)
      , max_pt_numerator_(max_pt_numerator)
      , min_pt_denominator_(min_pt_denominator)
      , max_pt_denominator_(max_pt_denominator)
      , min_absEta_(min_absEta)
      , max_absEta_(max_absEta)
      , decayMode_(decayMode)
      , min_leadTrackPt_(min_leadTrackPt)
      , max_leadTrackPt_(max_leadTrackPt)
      , min_relDiscriminator_(min_relDiscriminator)
      , max_relDiscriminator_(max_relDiscriminator)
      , min_absDiscriminator_(min_absDiscriminator)
      , max_absDiscriminator_(max_absDiscriminator)
      , dRmatch_(0.3)
    {}
    ~efficiencyPlotEntryType() 
    {}
    void bookHistograms(DQMStore& dqmStore)
    {
      TString histogramName_suffix = decayMode_.data(); 
      histogramName_suffix = getHistogramName(histogramName_suffix, "pt_numerator", min_pt_numerator_, max_pt_numerator_, 0);
      histogramName_suffix = getHistogramName(histogramName_suffix, "absEta", min_absEta_, max_absEta_, 2);
      histogramName_suffix = getHistogramName(histogramName_suffix, "leadTrackPt", min_leadTrackPt_, max_leadTrackPt_, 0);
      histogramName_suffix = getHistogramName(histogramName_suffix, "relDiscriminator", min_relDiscriminator_, max_relDiscriminator_, 3);
      histogramName_suffix = getHistogramName(histogramName_suffix, "absDiscriminator", min_absDiscriminator_, max_absDiscriminator_, 3);
      if ( !(min_relDiscriminator_ >= 0. || max_relDiscriminator_ > 0. || min_absDiscriminator_ >= 0. || max_absDiscriminator_ > 0.) )
      {
        histogramName_suffix.Append("_noIsolation");
      }
      histogramName_suffix = histogramName_suffix.ReplaceAll(".", "p");

      TString histogramName_pt_numerator = Form("effPFTau_vs_pt_numerator_%s", histogramName_suffix.Data());
      me_pt_numerator_ = dqmStore.book1D(histogramName_pt_numerator.Data(), histogramName_pt_numerator.Data(), 16, 20., 100.);
      histogram_pt_numerator_ = me_pt_numerator_->getTH1();
      assert(histogram_pt_numerator_);
      TString histogramName_pt_denominator = Form("effPFTau_vs_pt_denominator_%s", histogramName_suffix.Data());
      me_pt_denominator_ = dqmStore.book1D(histogramName_pt_denominator.Data(), histogramName_pt_denominator.Data(), 16, 20., 100.);
      histogram_pt_denominator_ = me_pt_denominator_->getTH1();
      assert(histogram_pt_denominator_);

      TString histogramName_eta_numerator = Form("effPFTau_vs_eta_numerator_%s", histogramName_suffix.Data());
      me_eta_numerator_ = dqmStore.book1D(histogramName_eta_numerator.Data(), histogramName_eta_numerator.Data(), 30, -3., +3.);
      histogram_eta_numerator_ = me_eta_numerator_->getTH1();
      assert(histogram_eta_numerator_);
      TString histogramName_eta_denominator = Form("effPFTau_vs_eta_denominator_%s", histogramName_suffix.Data());
      me_eta_denominator_ = dqmStore.book1D(histogramName_eta_denominator.Data(), histogramName_eta_denominator.Data(), 30, -3., +3.);
      histogram_eta_denominator_ = me_eta_denominator_->getTH1();
      assert(histogram_eta_denominator_);

      TString histogramName_phi_numerator = Form("effPFTau_vs_phi_numerator_%s", histogramName_suffix.Data());
      me_phi_numerator_ = dqmStore.book1D(histogramName_phi_numerator.Data(), histogramName_phi_numerator.Data(), 18, -TMath::Pi(), +TMath::Pi());
      histogram_phi_numerator_ = me_phi_numerator_->getTH1();
      assert(histogram_phi_numerator_);
      TString histogramName_phi_denominator = Form("effPFTau_vs_phi_denominator_%s", histogramName_suffix.Data());
      me_phi_denominator_ = dqmStore.book1D(histogramName_phi_denominator.Data(), histogramName_phi_denominator.Data(), 18, -TMath::Pi(), +TMath::Pi());
      histogram_phi_denominator_ = me_phi_denominator_->getTH1();
      assert(histogram_phi_denominator_);

      TString histogramName_minDeltaR_numerator = Form("effPFTau_vs_minDeltaR_numerator_%s", histogramName_suffix.Data());
      me_minDeltaR_numerator_ = dqmStore.book1D(histogramName_minDeltaR_numerator.Data(), histogramName_minDeltaR_numerator.Data(), 50, 0., 5.); 
      histogram_minDeltaR_numerator_ = me_minDeltaR_numerator_->getTH1();
      assert(histogram_minDeltaR_numerator_);
      TString histogramName_minDeltaR_denominator = Form("effPFTau_vs_minDeltaR_denominator_%s", histogramName_suffix.Data());
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

      TString histogramName_pt_vs_absEta_numerator = Form("effPFTau_vs_pt_vs_absEta_numerator_%s", histogramName_suffix.Data());
      me_pt_vs_absEta_numerator_ = dqmStore.book2D(
        histogramName_pt_vs_absEta_numerator.Data(), histogramName_pt_vs_absEta_numerator.Data(), 
        numBins_absEta, binning_absEta, numBins_pt, binning_pt);
      histogram_pt_vs_absEta_numerator_ = dynamic_cast<TH2*>(me_pt_vs_absEta_numerator_->getTH1());
      assert(histogram_pt_vs_absEta_numerator_);
      TString histogramName_pt_vs_absEta_denominator = Form("effPFTau_vs_pt_vs_absEta_denominator_%s", histogramName_suffix.Data());
      me_pt_vs_absEta_denominator_ = dqmStore.book2D(
        histogramName_pt_vs_absEta_denominator.Data(), histogramName_pt_vs_absEta_denominator.Data(), 
        numBins_absEta, binning_absEta, numBins_pt, binning_pt);
      histogram_pt_vs_absEta_denominator_ = dynamic_cast<TH2*>(me_pt_vs_absEta_denominator_->getTH1());
      assert(histogram_pt_vs_absEta_denominator_);
    }
    void fillHistograms(const std::vector<BaseTau>& baseTaus,
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
	if ( !(decayMode_ == "all" || denominatorTau_decayMode == decayMode_) ) continue;

        bool isMatched = false;

        for ( std::vector<BaseTau>::const_iterator baseTau = baseTaus.begin();
            baseTau != baseTaus.end(); ++baseTau ) {
          double baseTau_absEta = TMath::Abs(baseTau->p4().eta());
          if ( (min_pt_numerator_     < 0. || baseTau->p4().pt()          >=  min_pt_numerator_                        ) &&
	       (max_pt_numerator_     < 0. || baseTau->p4().pt()          <=  max_pt_numerator_                        ) &&
               (min_absEta_           < 0. || baseTau_absEta              >=  min_absEta_                              ) &&
	       (max_absEta_           < 0. || baseTau_absEta              <=  max_absEta_                              ) &&
               (min_leadTrackPt_      < 0. || baseTau->leadTrackP4().pt() >=  min_leadTrackPt_                         ) &&
               (max_leadTrackPt_      < 0. || baseTau->leadTrackP4().pt() >=  max_leadTrackPt_                         ) &&
               (min_relDiscriminator_ < 0. || baseTau->discriminator()    >= (min_relDiscriminator_*baseTau->p4().pt())) &&
               (max_relDiscriminator_ < 0. || baseTau->discriminator()    <= (max_relDiscriminator_*baseTau->p4().pt())) &&
               (min_absDiscriminator_ < 0. || baseTau->discriminator()    >=  min_absDiscriminator_                    ) &&
               (max_absDiscriminator_ < 0. || baseTau->discriminator()    <=  max_absDiscriminator_                    ) )
          {
            double dR = reco::deltaR(denominatorTau.eta(), denominatorTau.phi(), baseTau->p4().eta(), baseTau->p4().phi());
	    if ( dR < dRmatch_ ) isMatched = true;
          }
	}

	double denominatorTau_absEta = TMath::Abs(denominatorTau.eta());
	if ( denominatorTau_absEta > min_absEta_ && denominatorTau_absEta < max_absEta_ ) 
	{
	  histogram_pt_denominator_->Fill(denominatorTau.pt(), evtWeight);
	  if ( isMatched ) histogram_pt_numerator_->Fill(denominatorTau.pt(), evtWeight);
	}
	if ( denominatorTau.pt() > min_pt_denominator_ && denominatorTau.pt() < max_pt_denominator_ ) 
	{
	  histogram_eta_denominator_->Fill(denominatorTau.eta(), evtWeight);
	  if ( isMatched ) histogram_eta_numerator_->Fill(denominatorTau.eta(), evtWeight);
	}
	if ( denominatorTau.pt() > min_pt_denominator_ && denominatorTau.pt() < max_pt_denominator_ && denominatorTau_absEta > min_absEta_ && denominatorTau_absEta < max_absEta_ )
	{
	  histogram_phi_denominator_->Fill(denominatorTau.phi(), evtWeight);
	  if ( isMatched ) histogram_phi_numerator_->Fill(denominatorTau.phi(), evtWeight);
	}
	if ( denominatorTau.pt() > min_pt_denominator_ && denominatorTau.pt() < max_pt_denominator_ && denominatorTau_absEta > min_absEta_ && denominatorTau_absEta < max_absEta_ )
	{
	  histogram_minDeltaR_denominator_->Fill(minDeltaR, evtWeight);
	  if ( isMatched ) histogram_minDeltaR_numerator_->Fill(minDeltaR, evtWeight);
	}
	histogram_pt_vs_absEta_denominator_->Fill(denominatorTau_absEta, denominatorTau.pt(), evtWeight);
	if ( isMatched ) histogram_pt_vs_absEta_numerator_->Fill(denominatorTau_absEta, denominatorTau.pt(), evtWeight);
      }
    }
    void fillHistograms(const std::vector<BaseTau>& baseTaus,
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
	if ( !((decayMode_ == "all"                                                                            )  ||
               (decayMode_ == "oneProng0Pi0"   && denominatorTau.decayMode() == reco::PFTau::kOneProng0PiZero  )  ||
	       (decayMode_ == "oneProng1Pi0"   && denominatorTau.decayMode() == reco::PFTau::kOneProng1PiZero  )  ||
	       (decayMode_ == "oneProng2Pi0"   && denominatorTau.decayMode() == reco::PFTau::kOneProng2PiZero  )  ||
	       (decayMode_ == "threeProng0Pi0" && denominatorTau.decayMode() == reco::PFTau::kThreeProng0PiZero)  ||
	       (decayMode_ == "threeProng1Pi0" && denominatorTau.decayMode() == reco::PFTau::kThreeProng1PiZero)) ) continue;	       

        bool isMatched = false;

        for ( std::vector<BaseTau>::const_iterator baseTau = baseTaus.begin();
            baseTau != baseTaus.end(); ++baseTau ) {
          double baseTau_absEta = TMath::Abs(baseTau->p4().eta());
          if ( (min_pt_numerator_     < 0. || baseTau->p4().pt()          >=  min_pt_numerator_                        ) &&
	       (max_pt_numerator_     < 0. || baseTau->p4().pt()          <=  max_pt_numerator_                        ) &&
               (min_absEta_           < 0. || baseTau_absEta              >=  min_absEta_                              ) &&
	       (max_absEta_           < 0. || baseTau_absEta              <=  max_absEta_                              ) &&
               (min_leadTrackPt_      < 0. || baseTau->leadTrackP4().pt() >=  min_leadTrackPt_                         ) &&
               (max_leadTrackPt_      < 0. || baseTau->leadTrackP4().pt() >=  max_leadTrackPt_                         ) &&
               (min_relDiscriminator_ < 0. || baseTau->discriminator()    >= (min_relDiscriminator_*baseTau->p4().pt())) &&
               (max_relDiscriminator_ < 0. || baseTau->discriminator()    <= (max_relDiscriminator_*baseTau->p4().pt())) &&
               (min_absDiscriminator_ < 0. || baseTau->discriminator()    >=  min_absDiscriminator_                    ) &&
               (max_absDiscriminator_ < 0. || baseTau->discriminator()    <=  max_absDiscriminator_                    ) )
          {
            double dR = reco::deltaR(denominatorTau.eta(), denominatorTau.phi(), baseTau->p4().eta(), baseTau->p4().phi());
	    if ( dR < dRmatch_ ) isMatched = true;
          }
	}

	double denominatorTau_absEta = TMath::Abs(denominatorTau.eta());
	if ( denominatorTau_absEta > min_absEta_ && denominatorTau_absEta < max_absEta_ ) 
	{
	  histogram_pt_denominator_->Fill(denominatorTau.pt(), evtWeight);
	  if ( isMatched ) histogram_pt_numerator_->Fill(denominatorTau.pt(), evtWeight);
	}
	if ( denominatorTau.pt() > min_pt_denominator_ && denominatorTau.pt() < max_pt_denominator_ ) 
	{
	  histogram_eta_denominator_->Fill(denominatorTau.eta(), evtWeight);
	  if ( isMatched ) histogram_eta_numerator_->Fill(denominatorTau.eta(), evtWeight);
	}
	if ( denominatorTau.pt() > min_pt_denominator_ && denominatorTau.pt() < max_pt_denominator_ && denominatorTau_absEta > min_absEta_ && denominatorTau_absEta < max_absEta_ )
	{
	  histogram_phi_denominator_->Fill(denominatorTau.phi(), evtWeight);
	  if ( isMatched ) histogram_phi_numerator_->Fill(denominatorTau.phi(), evtWeight);
	}
	if ( denominatorTau.pt() > min_pt_denominator_ && denominatorTau.pt() < max_pt_denominator_ && denominatorTau_absEta > min_absEta_ && denominatorTau_absEta < max_absEta_ )
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
    // (1st set of) cuts applied to HLT trigger taus in numerator 
    double min_pt_numerator_;
    double max_pt_numerator_;
    // cuts applied to offline and generator-level taus in denominator 
    double min_pt_denominator_;
    double max_pt_denominator_;
    double min_absEta_;
    double max_absEta_;
    std::string decayMode_;
    // (2nd set of) cuts applied to HLT trigger taus in numerator 
    double min_leadTrackPt_;
    double max_leadTrackPt_;
    double min_relDiscriminator_;
    double max_relDiscriminator_;
    double min_absDiscriminator_;
    double max_absDiscriminator_;
    // matching between offline or generator-level taus in denominator and L1 trigger taus in numerator
    double dRmatch_; 
  };
  std::vector<efficiencyPlotEntryType*> efficiencyPlots_;
};

#endif   
