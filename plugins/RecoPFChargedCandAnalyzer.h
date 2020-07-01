#ifndef HLTTrigger_TallinnHLTPFTauAnalyzer_RecoPFChargedCandAnalyzer_h
#define HLTTrigger_TallinnHLTPFTauAnalyzer_RecoPFChargedCandAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"    // reco::PFCandidate
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h" // reco::PFCandidateCollection

#include <TH1.h>     // TH1
#include <TString.h> // TString, Form()
#include <TMath.h>   // TMath::Abs(), TMath::Pi()

#include <vector>    // std::vector
#include <string>    // std::string

using namespace dqm::implementation;

class RecoPFChargedCandAnalyzer : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit RecoPFChargedCandAnalyzer(const edm::ParameterSet&);
    
  // destructor
  ~RecoPFChargedCandAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<reco::PFCandidateCollection> token_;

  std::string dqmDirectory_;

  struct plotEntryType
  {
    plotEntryType(double min_pt, double max_pt, double min_absEta, double max_absEta)
      : me_pfCandPt_e_(nullptr)
      , histogram_pfCandPt_e_(nullptr)
      , me_trackPt_e_(nullptr)
      , histogram_trackPt_e_(nullptr)
      , me_trackPt_div_pfCandPt_e_(nullptr)
      , histogram_trackPt_div_pfCandPt_e_(nullptr)
      , me_pfCandPt_mu_(nullptr)
      , histogram_pfCandPt_mu_(nullptr)
      , me_trackPt_mu_(nullptr)
      , histogram_trackPt_mu_(nullptr)
      , me_trackPt_div_pfCandPt_mu_(nullptr)
      , histogram_trackPt_div_pfCandPt_mu_(nullptr)
      , me_pfCandPt_h_(nullptr)
      , histogram_pfCandPt_h_(nullptr)
      , me_trackPt_h_(nullptr)
      , histogram_trackPt_h_(nullptr)
      , me_trackPt_div_pfCandPt_h_(nullptr)
      , histogram_trackPt_div_pfCandPt_h_(nullptr)
      , me_pfCandPt_other_(nullptr)
      , histogram_pfCandPt_other_(nullptr)
      , me_trackPt_other_(nullptr)
      , histogram_trackPt_other_(nullptr)
      , me_trackPt_div_pfCandPt_other_(nullptr)
      , histogram_trackPt_div_pfCandPt_other_(nullptr)
      , min_pt_(min_pt)
      , max_pt_(max_pt)
      , min_absEta_(min_absEta)
      , max_absEta_(max_absEta)
    {}
    ~plotEntryType() 
    {}
    void bookHistograms(DQMStore& dqmStore)
    {
      TString histogramName_suffix;
      if      ( min_absEta_ >= 0. && max_absEta_ > 0. ) histogramName_suffix.Append(Form("_absEta%1.2fto%1.2f", min_absEta_, max_absEta_));
      else if ( min_absEta_ >= 0.                     ) histogramName_suffix.Append(Form("_absEtaGt%1.2f", min_absEta_));
      else if (                      max_absEta_ > 0. ) histogramName_suffix.Append(Form("_absEtaLt%1.2f", max_absEta_));
      if      ( min_pt_     >= 0. && max_pt_     > 0. ) histogramName_suffix.Append(Form("_pt%1.0fto%1.0f", min_pt_, max_pt_));
      else if ( min_pt_     >= 0.                     ) histogramName_suffix.Append(Form("_ptGt%1.0f", min_pt_));
      else if (                      max_pt_     > 0. ) histogramName_suffix.Append(Form("_ptLt%1.0f", max_pt_));
      histogramName_suffix = histogramName_suffix.ReplaceAll(".", "p");

      const int numBins_pt = 19;
      float binning_pt[numBins_pt + 1] = { 
        0., 1., 2., 3., 4., 5., 6., 8., 10., 12., 15., 20., 25., 30., 35., 40., 50., 60., 80., 100.
      };

      TString histogramName_pfCandPt_e = Form("pfCandPt_e%s", histogramName_suffix.Data());
      me_pfCandPt_e_ = dqmStore.book1D(histogramName_pfCandPt_e.Data(), histogramName_pfCandPt_e.Data(), numBins_pt, binning_pt);
      histogram_pfCandPt_e_ = me_pfCandPt_e_->getTH1();
      assert(histogram_pfCandPt_e_);
      TString histogramName_trackPt_e = Form("trackPt_e%s", histogramName_suffix.Data());
      me_trackPt_e_ = dqmStore.book1D(histogramName_trackPt_e.Data(), histogramName_trackPt_e.Data(), numBins_pt, binning_pt);
      histogram_trackPt_e_ = me_trackPt_e_->getTH1();
      assert(histogram_trackPt_e_);
      TString histogramName_trackPt_div_pfCandPt_e = Form("trackPt_div_pfCandPt_e%s", histogramName_suffix.Data());
      me_trackPt_div_pfCandPt_e_ = dqmStore.book1D(histogramName_trackPt_div_pfCandPt_e.Data(), histogramName_trackPt_div_pfCandPt_e.Data(), 40, 0., 2.);
      histogram_trackPt_div_pfCandPt_e_ = me_trackPt_div_pfCandPt_e_->getTH1();
      assert(histogram_trackPt_div_pfCandPt_e_);

      TString histogramName_pfCandPt_mu = Form("pfCandPt_mu%s", histogramName_suffix.Data());
      me_pfCandPt_mu_ = dqmStore.book1D(histogramName_pfCandPt_mu.Data(), histogramName_pfCandPt_mu.Data(), numBins_pt, binning_pt);
      histogram_pfCandPt_mu_ = me_pfCandPt_mu_->getTH1();
      assert(histogram_pfCandPt_mu_);
      TString histogramName_trackPt_mu = Form("trackPt_mu%s", histogramName_suffix.Data());
      me_trackPt_mu_ = dqmStore.book1D(histogramName_trackPt_mu.Data(), histogramName_trackPt_mu.Data(), numBins_pt, binning_pt);
      histogram_trackPt_mu_ = me_trackPt_mu_->getTH1();
      assert(histogram_trackPt_mu_);
      TString histogramName_trackPt_div_pfCandPt_mu = Form("trackPt_div_pfCandPt_mu%s", histogramName_suffix.Data());
      me_trackPt_div_pfCandPt_mu_ = dqmStore.book1D(histogramName_trackPt_div_pfCandPt_mu.Data(), histogramName_trackPt_div_pfCandPt_mu.Data(), 40, 0., 2.);
      histogram_trackPt_div_pfCandPt_mu_ = me_trackPt_div_pfCandPt_mu_->getTH1();
      assert(histogram_trackPt_div_pfCandPt_mu_);

      TString histogramName_pfCandPt_h = Form("pfCandPt_h%s", histogramName_suffix.Data());
      me_pfCandPt_h_ = dqmStore.book1D(histogramName_pfCandPt_h.Data(), histogramName_pfCandPt_h.Data(), numBins_pt, binning_pt);
      histogram_pfCandPt_h_ = me_pfCandPt_h_->getTH1();
      assert(histogram_pfCandPt_h_);
      TString histogramName_trackPt_h = Form("trackPt_h%s", histogramName_suffix.Data());
      me_trackPt_h_ = dqmStore.book1D(histogramName_trackPt_h.Data(), histogramName_trackPt_h.Data(), numBins_pt, binning_pt);
      histogram_trackPt_h_ = me_trackPt_h_->getTH1();
      assert(histogram_trackPt_h_);
      TString histogramName_trackPt_div_pfCandPt_h = Form("trackPt_div_pfCandPt_h%s", histogramName_suffix.Data());
      me_trackPt_div_pfCandPt_h_ = dqmStore.book1D(histogramName_trackPt_div_pfCandPt_h.Data(), histogramName_trackPt_div_pfCandPt_h.Data(), 40, 0., 2.);
      histogram_trackPt_div_pfCandPt_h_ = me_trackPt_div_pfCandPt_h_->getTH1();
      assert(histogram_trackPt_div_pfCandPt_h_);

      TString histogramName_pfCandPt_other = Form("pfCandPt_other%s", histogramName_suffix.Data());
      me_pfCandPt_other_ = dqmStore.book1D(histogramName_pfCandPt_other.Data(), histogramName_pfCandPt_other.Data(), numBins_pt, binning_pt);
      histogram_pfCandPt_other_ = me_pfCandPt_other_->getTH1();
      assert(histogram_pfCandPt_other_);
      TString histogramName_trackPt_other = Form("trackPt_other%s", histogramName_suffix.Data());
      me_trackPt_other_ = dqmStore.book1D(histogramName_trackPt_other.Data(), histogramName_trackPt_other.Data(), numBins_pt, binning_pt);
      histogram_trackPt_other_ = me_trackPt_other_->getTH1();
      assert(histogram_trackPt_other_);
      TString histogramName_trackPt_div_pfCandPt_other = Form("trackPt_div_pfCandPt_other%s", histogramName_suffix.Data());
      me_trackPt_div_pfCandPt_other_ = dqmStore.book1D(histogramName_trackPt_div_pfCandPt_other.Data(), histogramName_trackPt_div_pfCandPt_other.Data(), 40, 0., 2.);
      histogram_trackPt_div_pfCandPt_other_ = me_trackPt_div_pfCandPt_other_->getTH1();
      assert(histogram_trackPt_div_pfCandPt_other_);
    }
    void fillHistograms(const reco::PFCandidateCollection& pfCands, double evtWeight)
    {
      for ( reco::PFCandidateCollection::const_iterator pfCand = pfCands.begin();
            pfCand != pfCands.end(); ++pfCand ) {
        double pfCand_absEta = TMath::Abs(pfCand->eta());
        if ( (min_absEta_ < 0. || pfCand_absEta       >=  min_absEta_) &&
	     (max_absEta_ < 0. || pfCand_absEta       <=  max_absEta_) &&
             (min_pt_     < 0. || pfCand->pt()        >=  min_pt_    ) &&
	     (max_pt_     < 0. || pfCand->pt()        <=  max_pt_    ) &&
             (                    pfCand->bestTrack()                ) )
        {
          double pfCand_absPdgId = TMath::Abs(pfCand->pdgId());
          TH1* histogram_pfCandPt             = nullptr;
          TH1* histogram_trackPt              = nullptr;
          TH1* histogram_trackPt_div_pfCandPt = nullptr;
          if      ( pfCand_absPdgId ==  11 )
          {
            histogram_pfCandPt             = histogram_pfCandPt_e_;
            histogram_trackPt              = histogram_trackPt_e_;
            histogram_trackPt_div_pfCandPt = histogram_trackPt_div_pfCandPt_e_;
          }
          else if ( pfCand_absPdgId ==  13 )
          {
            histogram_pfCandPt             = histogram_pfCandPt_mu_;
            histogram_trackPt              = histogram_trackPt_mu_;
            histogram_trackPt_div_pfCandPt = histogram_trackPt_div_pfCandPt_mu_;
          }
          else if ( pfCand_absPdgId == 211 )
          {
            histogram_pfCandPt             = histogram_pfCandPt_h_;
            histogram_trackPt              = histogram_trackPt_h_;
            histogram_trackPt_div_pfCandPt = histogram_trackPt_div_pfCandPt_h_;
          }
          else if ( TMath::Abs(pfCand->charge()) > 0.5 )
          {
            histogram_pfCandPt             = histogram_pfCandPt_other_;
            histogram_trackPt              = histogram_trackPt_other_;
            histogram_trackPt_div_pfCandPt = histogram_trackPt_div_pfCandPt_other_;
          }
          if ( histogram_pfCandPt && histogram_trackPt && histogram_trackPt_div_pfCandPt )
          {
            histogram_pfCandPt->Fill(pfCand->pt(), evtWeight);
            assert(pfCand->bestTrack());
            histogram_trackPt->Fill(pfCand->bestTrack()->pt(), evtWeight);
            if ( pfCand->pt() > 0.5 )
            {
              histogram_trackPt_div_pfCandPt->Fill(pfCand->bestTrack()->pt()/pfCand->pt(), evtWeight);
            }
          }
        }
      }
    }
    
    MonitorElement*me_pfCandPt_e_;
    TH1* histogram_pfCandPt_e_;
    MonitorElement* me_trackPt_e_;
    TH1* histogram_trackPt_e_;
    MonitorElement* me_trackPt_div_pfCandPt_e_;
    TH1* histogram_trackPt_div_pfCandPt_e_;
    MonitorElement* me_pfCandPt_mu_;
    TH1* histogram_pfCandPt_mu_;
    MonitorElement* me_trackPt_mu_;
    TH1* histogram_trackPt_mu_;
    MonitorElement* me_trackPt_div_pfCandPt_mu_;
    TH1* histogram_trackPt_div_pfCandPt_mu_;
    MonitorElement* me_pfCandPt_h_;
    TH1* histogram_pfCandPt_h_;
    MonitorElement* me_trackPt_h_;
    TH1* histogram_trackPt_h_;
    MonitorElement* me_trackPt_div_pfCandPt_h_;
    TH1* histogram_trackPt_div_pfCandPt_h_;
    MonitorElement* me_pfCandPt_other_;
    TH1* histogram_pfCandPt_other_;
    MonitorElement* me_trackPt_other_;
    TH1* histogram_trackPt_other_;
    MonitorElement* me_trackPt_div_pfCandPt_other_;
    TH1* histogram_trackPt_div_pfCandPt_other_;
    double min_pt_;
    double max_pt_;
    double min_absEta_;
    double max_absEta_;
  };
  std::vector<plotEntryType*> plots_;
};

#endif  
