#ifndef HLTrigger_HLTPFTauAnalyzer_GenJetAnalyzer_h
#define HLTrigger_HLTPFTauAnalyzer_GenJetAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/JetReco/interface/GenJet.h"             // reco::GenJet
#include "DataFormats/JetReco/interface/GenJetCollection.h"   // reco::GenJetCollection

#include "HLTrigger/TallinnHLTPFTauAnalyzer/interface/histogramAuxFunctions.h" // fillWithOverFlow

#include <TH1.h>     // TH1
#include <TString.h> // TString, Form()
#include <TMath.h>   // TMath::Abs(), TMath::Pi()

#include <vector>    // std::vector
#include <string>    // std::string
#include <algorithm> // std::sort()

using namespace dqm::implementation;

namespace
{
  bool isHigherPt(const reco::GenJet* genJet1, const reco::GenJet* genJet2)
  {
    return (genJet1->pt() > genJet2->pt());
  }
}

class GenJetAnalyzer : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit GenJetAnalyzer(const edm::ParameterSet&);
    
  // destructor
  ~GenJetAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<reco::GenJetCollection> token_;

  double lumiScale_;

  std::string dqmDirectory_;

  struct plotEntryType
  {
    plotEntryType(double min_pt, double max_pt, double min_absEta, double max_absEta)
      : me_leadJet_pt_(nullptr)
      , histogram_leadJet_pt_(nullptr)
      , me_leadJet_eta_(nullptr)
      , histogram_leadJet_eta_(nullptr)
      , me_subleadJet_pt_(nullptr)
      , histogram_subleadJet_pt_(nullptr)
      , me_subleadJet_eta_(nullptr)
      , histogram_subleadJet_eta_(nullptr)
      , me_numJets_(nullptr)
      , histogram_numJets_(nullptr)
      , me_jet_pt_(nullptr)
      , histogram_jet_pt_(nullptr)
      , me_jet_eta_(nullptr)
      , histogram_jet_eta_(nullptr)
      , me_HT_(nullptr)
      , histogram_HT_(nullptr)
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
      if      ( min_pt_     >= 0. && max_pt_     > 0. ) histogramName_suffix.Append(Form("_pt%1.0fto%1.0f", min_pt_, max_pt_));
      else if ( min_pt_     >= 0.                     ) histogramName_suffix.Append(Form("_ptGt%1.0f", min_pt_));
      else if (                      max_pt_     > 0. ) histogramName_suffix.Append(Form("_ptLt%1.0f", max_pt_));
      if      ( min_absEta_ >= 0. && max_absEta_ > 0. ) histogramName_suffix.Append(Form("_absEta%1.2fto%1.2f", min_absEta_, max_absEta_));
      else if ( min_absEta_ >= 0.                     ) histogramName_suffix.Append(Form("_absEtaGt%1.2f", min_absEta_));
      else if (                      max_absEta_ > 0. ) histogramName_suffix.Append(Form("_absEtaLt%1.2f", max_absEta_));      
      histogramName_suffix = histogramName_suffix.ReplaceAll(".", "p");

      TString histogramName_leadJet_pt = Form("leadJet_pt%s", histogramName_suffix.Data());
      me_leadJet_pt_ = dqmStore.book1D(histogramName_leadJet_pt.Data(), histogramName_leadJet_pt.Data(), 100, 0., 500.);
      histogram_leadJet_pt_ = me_leadJet_pt_->getTH1();
      assert(histogram_leadJet_pt_);
      TString histogramName_leadJet_eta = Form("leadJet_eta%s", histogramName_suffix.Data());
      me_leadJet_eta_ = dqmStore.book1D(histogramName_leadJet_eta.Data(), histogramName_leadJet_eta.Data(), 100, -5.0, +5.0);
      histogram_leadJet_eta_ = me_leadJet_eta_->getTH1();
      assert(histogram_leadJet_eta_);

      TString histogramName_subleadJet_pt = Form("subleadJet_pt%s", histogramName_suffix.Data());
      me_subleadJet_pt_ = dqmStore.book1D(histogramName_subleadJet_pt.Data(), histogramName_subleadJet_pt.Data(), 100, 0., 500.);
      histogram_subleadJet_pt_ = me_subleadJet_pt_->getTH1();
      assert(histogram_subleadJet_pt_);
      TString histogramName_subleadJet_eta = Form("subleadJet_eta%s", histogramName_suffix.Data());
      me_subleadJet_eta_ = dqmStore.book1D(histogramName_subleadJet_eta.Data(), histogramName_subleadJet_eta.Data(), 100, -5.0, +5.0);
      histogram_subleadJet_eta_ = me_subleadJet_eta_->getTH1();
      assert(histogram_subleadJet_eta_);

      TString histogramName_numJets = Form("numJets%s", histogramName_suffix.Data());
      me_numJets_ = dqmStore.book1D(histogramName_numJets.Data(), histogramName_numJets.Data(), 20, -0.5, 19.5);
      histogram_numJets_ = me_numJets_->getTH1();
      assert(histogram_numJets_);

      TString histogramName_jet_pt = Form("jet_pt%s", histogramName_suffix.Data());
      me_jet_pt_ = dqmStore.book1D(histogramName_jet_pt.Data(), histogramName_jet_pt.Data(), 100, 0., 500.);
      histogram_jet_pt_ = me_jet_pt_->getTH1();
      assert(histogram_jet_pt_);
      TString histogramName_jet_eta = Form("jet_eta%s", histogramName_suffix.Data());
      me_jet_eta_ = dqmStore.book1D(histogramName_jet_eta.Data(), histogramName_jet_eta.Data(), 100, -5.0, +5.0);
      histogram_jet_eta_ = me_jet_eta_->getTH1();
      assert(histogram_jet_eta_);

      TString histogramName_HT = Form("HT%s", histogramName_suffix.Data());
      me_HT_ = dqmStore.book1D(histogramName_HT.Data(), histogramName_HT.Data(), 100, 0., 1000.);
      histogram_HT_ = me_HT_->getTH1();
      assert(histogram_HT_);
    }
    void fillHistograms(const reco::GenJetCollection& genJets, double evtWeight)
    {
      std::vector<const reco::GenJet*> genjets_selected;
      for ( reco::GenJetCollection::const_iterator genJet = genJets.begin();
            genJet != genJets.end(); ++genJet ) {
        if ( (min_pt_     < 0. || genJet->pt()              > min_pt_     ) && 
             (max_pt_     < 0. || genJet->pt()              < max_pt_     ) && 
             (min_absEta_ < 0. || TMath::Abs(genJet->eta()) > min_absEta_ ) &&
             (max_absEta_ < 0. || TMath::Abs(genJet->eta()) < max_absEta_ ) ) 
	{
          genjets_selected.push_back(&(*genJet));
        }
      }
      // CV: sort genJets by decreasing pT
      std::sort(genjets_selected.begin(), genjets_selected.end(), isHigherPt);

      if ( genjets_selected.size() >= 1 )
      {
        const reco::GenJet* genJet_lead = genjets_selected[0];
        fillWithOverFlow(histogram_leadJet_pt_, genJet_lead->pt(), evtWeight);
        histogram_leadJet_eta_->Fill(genJet_lead->eta(), evtWeight);
      }
      if ( genjets_selected.size() >= 2 )
      {
        const reco::GenJet* genJet_sublead = genjets_selected[1];
        fillWithOverFlow(histogram_subleadJet_pt_, genJet_sublead->pt(), evtWeight);
        histogram_subleadJet_eta_->Fill(genJet_sublead->eta(), evtWeight);
      }

      fillWithOverFlow(histogram_numJets_, genjets_selected.size(), evtWeight);

      for ( std::vector<const reco::GenJet*>::const_iterator genJet = genjets_selected.begin();
            genJet != genjets_selected.end(); ++genJet ) {
        fillWithOverFlow(histogram_jet_pt_, (*genJet)->pt(), evtWeight);
        histogram_jet_eta_->Fill((*genJet)->eta(), evtWeight);
      }

      double HT = 0.;
      for ( reco::GenJetCollection::const_iterator genJet = genJets.begin();
            genJet != genJets.end(); ++genJet ) {
        HT += genJet->pt();
      }
      fillWithOverFlow(histogram_HT_, HT, evtWeight);
    }
    MonitorElement* me_leadJet_pt_;
    TH1* histogram_leadJet_pt_;
    MonitorElement* me_leadJet_eta_;
    TH1* histogram_leadJet_eta_;
    MonitorElement* me_subleadJet_pt_;
    TH1* histogram_subleadJet_pt_;
    MonitorElement* me_subleadJet_eta_;
    TH1* histogram_subleadJet_eta_;
    MonitorElement* me_numJets_;
    TH1* histogram_numJets_;
    MonitorElement* me_jet_pt_;
    TH1* histogram_jet_pt_;
    MonitorElement* me_jet_eta_;
    TH1* histogram_jet_eta_;
    MonitorElement* me_HT_;
    TH1* histogram_HT_;
    // cuts applied to generator-level jets
    double min_pt_;
    double max_pt_;
    double min_absEta_;
    double max_absEta_;
  };
  std::vector<plotEntryType*> plots_;
};

#endif   
