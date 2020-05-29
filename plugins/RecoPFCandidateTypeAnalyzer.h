#ifndef HLTTrigger_TallinnHLTPFTauAnalyzer_RecoPFCandidateTypeAnalyzer_h
#define HLTTrigger_TallinnHLTPFTauAnalyzer_RecoPFCandidateTypeAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "RecoTauTag/RecoTau/interface/RecoTauQualityCuts.h"                  // reco::tau::RecoTauQualityCuts
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"          // reco::PFCandidate
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"       // reco::PFCandidateCollection
#include "DataFormats/VertexReco/interface/Vertex.h"                          // reco::Vertex
#include "DataFormats/VertexReco/interface/VertexFwd.h"                       // reco::VertexCollection
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"     // PileupSummaryInfo

#include "HLTTrigger/TallinnHLTPFTauAnalyzer/interface/histogramAuxFunctions.h" // divideByBinWidth

#include <TH1.h>     // TH1
#include <TString.h> // TString, Form()

#include <vector>    // std::vector
#include <string>    // std::string

using namespace dqm::implementation;

class RecoPFCandidateTypeAnalyzer : public edm::EDAnalyzer 
{
 public:
  RecoPFCandidateTypeAnalyzer(const edm::ParameterSet& cfg);
  ~RecoPFCandidateTypeAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag src_pfCands_;
  edm::EDGetTokenT<reco::PFCandidateCollection> token_pfCands_;
  
  edm::InputTag srcVertices_;
  edm::EDGetTokenT<reco::VertexCollection> tokenVertices_;

  typedef std::vector<PileupSummaryInfo> PileupSummaryInfoCollection;
  edm::InputTag srcPileupSummaryInfo_;
  edm::EDGetTokenT<PileupSummaryInfoCollection> tokenPileupSummaryInfo_;

  reco::tau::RecoTauQualityCuts* isolationQualityCuts_dzCut_disabled_;
  reco::tau::RecoTauQualityCuts* isolationQualityCuts_dzCut_enabled_primary_;

  std::string dqmDirectory_;

  struct pfCandTypePlotEntryType
  {
    pfCandTypePlotEntryType(double min_pt, double max_absEta, const std::string& label)
      : label_(label)
      , min_pt_(min_pt)
      , max_absEta_(max_absEta)
      , me_energyFraction_vs_eta_(nullptr)
      , histogram_energyFraction_vs_eta_(nullptr)
      , me_ptFraction_vs_eta_(nullptr)
      , histogram_ptFraction_vs_eta_(nullptr)
      , me_energyFraction_vs_eta_fine_binning_(nullptr)
      , histogram_energyFraction_vs_eta_fine_binning_(nullptr)
      , me_ptFraction_vs_eta_fine_binning_(nullptr)
      , histogram_ptFraction_vs_eta_fine_binning_(nullptr)
      , me_energyFraction_vs_absEta_(nullptr)
      , histogram_energyFraction_vs_absEta_(nullptr)
      , me_ptFraction_vs_absEta_(nullptr)
      , histogram_ptFraction_vs_absEta_(nullptr)
      , me_energyFraction_vs_pt_(nullptr)
      , histogram_energyFraction_vs_pt_(nullptr)
      , me_ptFraction_vs_pt_(nullptr)
      , histogram_ptFraction_vs_pt_(nullptr)
    {}
    ~pfCandTypePlotEntryType()
    {}
    void bookHistograms(DQMStore& dqmStore)
    {
      std::string histogramName_energyFraction_vs_eta = Form("%sEnergyFraction_vs_eta", label_.data());
      me_energyFraction_vs_eta_ = dqmStore.book1D(
        histogramName_energyFraction_vs_eta.data(), histogramName_energyFraction_vs_eta.data(), 
	68, -34*0.087, +34*0.087); // 0.087 = size of one CaloTower (5 ECAL crystals) in eta
      histogram_energyFraction_vs_eta_ = me_energyFraction_vs_eta_->getTH1();
      assert(histogram_energyFraction_vs_eta_);

      std::string histogramName_ptFraction_vs_eta = Form("%sPtFraction_vs_eta", label_.data());
      me_ptFraction_vs_eta_ = dqmStore.book1D(
        histogramName_ptFraction_vs_eta.data(), histogramName_ptFraction_vs_eta.data(), 
	68, -34*0.087, +34*0.087); // 0.087 = size of one CaloTower (5 ECAL crystals) in eta
      histogram_ptFraction_vs_eta_ = me_ptFraction_vs_eta_->getTH1();
      assert(histogram_ptFraction_vs_eta_);

      std::string histogramName_energyFraction_vs_eta_fine_binning = Form("%sEnergyFraction_vs_eta_fine_binning", label_.data());
      me_energyFraction_vs_eta_fine_binning_ = dqmStore.book1D(
        histogramName_energyFraction_vs_eta_fine_binning.data(), histogramName_energyFraction_vs_eta_fine_binning.data(), 
	6000, -3.0, +3.0);
      histogram_energyFraction_vs_eta_fine_binning_ = me_energyFraction_vs_eta_fine_binning_->getTH1();
      assert(histogram_energyFraction_vs_eta_fine_binning_);

      std::string histogramName_ptFraction_vs_eta_fine_binning = Form("%sPtFraction_vs_eta_fine_binning", label_.data());
      me_ptFraction_vs_eta_fine_binning_ = dqmStore.book1D(
        histogramName_ptFraction_vs_eta_fine_binning.data(), histogramName_ptFraction_vs_eta_fine_binning.data(), 
	6000, -3.0, +3.0);
      histogram_ptFraction_vs_eta_fine_binning_ = me_ptFraction_vs_eta_fine_binning_->getTH1();
      assert(histogram_ptFraction_vs_eta_fine_binning_);

      std::string histogramName_energyFraction_vs_absEta = Form("%sEnergyFraction_vs_absEta", label_.data());
      me_energyFraction_vs_absEta_ = dqmStore.book1D(
        histogramName_energyFraction_vs_absEta.data(), histogramName_energyFraction_vs_absEta.data(), 
	34, 0., 34*0.087); // 0.087 = size of one CaloTower (5 ECAL crystals) in eta
      histogram_energyFraction_vs_absEta_ = me_energyFraction_vs_absEta_->getTH1();
      assert(histogram_energyFraction_vs_absEta_);

      std::string histogramName_ptFraction_vs_absEta = Form("%sPtFraction_vs_absEta", label_.data());
      me_ptFraction_vs_absEta_ = dqmStore.book1D(
        histogramName_ptFraction_vs_absEta.data(), histogramName_ptFraction_vs_absEta.data(), 
	34, 0., 34*0.087); // 0.087 = size of one CaloTower (5 ECAL crystals) in eta
      histogram_ptFraction_vs_absEta_ = me_ptFraction_vs_absEta_->getTH1();
      assert(histogram_ptFraction_vs_absEta_);

      const int numBins_pt = 15;
      float binning_pt[numBins_pt + 1] = { 
        0., 1., 2., 4., 6., 10., 15., 25., 40., 60., 100., 150., 250., 400., 600., 1000.
      };

      std::string histogramName_energyFraction_vs_pt = Form("%sEnergyFraction_vs_pt", label_.data());
      me_energyFraction_vs_pt_ = dqmStore.book1D(
        histogramName_energyFraction_vs_pt.data(), histogramName_energyFraction_vs_pt.data(), 
	numBins_pt, binning_pt);
      histogram_energyFraction_vs_pt_ = me_energyFraction_vs_pt_->getTH1();
      assert(histogram_energyFraction_vs_pt_);

      std::string histogramName_ptFraction_vs_pt = Form("%sPtFraction_vs_pt", label_.data());
      me_ptFraction_vs_pt_ = dqmStore.book1D(
        histogramName_ptFraction_vs_pt.data(), histogramName_ptFraction_vs_pt.data(), 
	numBins_pt, binning_pt);
      histogram_ptFraction_vs_pt_ = me_ptFraction_vs_pt_->getTH1();
      assert(histogram_ptFraction_vs_pt_);

      std::string histogramName_energyFraction_vs_numPileup = Form("%sEnergyFraction_vs_numPileup", label_.data());
      me_energyFraction_vs_numPileup_ = dqmStore.book1D(
        histogramName_energyFraction_vs_numPileup.data(), histogramName_energyFraction_vs_numPileup.data(), 
	40, 0., 400.);
      histogram_energyFraction_vs_numPileup_ = me_energyFraction_vs_numPileup_->getTH1();
      assert(histogram_energyFraction_vs_numPileup_);

      std::string histogramName_ptFraction_vs_numPileup = Form("%sPtFraction_vs_numPileup", label_.data());
      me_ptFraction_vs_numPileup_ = dqmStore.book1D(
        histogramName_ptFraction_vs_numPileup.data(), histogramName_ptFraction_vs_numPileup.data(), 
	40, 0., 400.);
      histogram_ptFraction_vs_numPileup_ = me_ptFraction_vs_numPileup_->getTH1();
      assert(histogram_ptFraction_vs_numPileup_);
    }
    void fillHistograms(const reco::PFCandidate& pfCand, int numPileup, double evtWeight)
    {
      double weight_energyFraction = pfCand.energy()*evtWeight;
      double weight_ptFraction = pfCand.pt()*evtWeight;
      if ( pfCand.pt() > min_pt_ ) 
      {
	histogram_energyFraction_vs_eta_->Fill(pfCand.eta(), weight_energyFraction);
	histogram_ptFraction_vs_eta_->Fill(pfCand.eta(), weight_ptFraction);
	histogram_energyFraction_vs_eta_fine_binning_->Fill(pfCand.eta(), weight_energyFraction);
	histogram_ptFraction_vs_eta_fine_binning_->Fill(pfCand.eta(), weight_ptFraction);
	histogram_energyFraction_vs_absEta_->Fill(TMath::Abs(pfCand.eta()), weight_energyFraction);
	histogram_ptFraction_vs_absEta_->Fill(TMath::Abs(pfCand.eta()), weight_ptFraction);
      }
      if ( TMath::Abs(pfCand.eta()) < max_absEta_ ) 
      {
	histogram_energyFraction_vs_pt_->Fill(pfCand.pt(), weight_energyFraction);
	histogram_ptFraction_vs_pt_->Fill(pfCand.pt(), weight_ptFraction);
      }
      if ( pfCand.pt() > min_pt_ && TMath::Abs(pfCand.eta()) < max_absEta_ && numPileup >= 0 ) 
      {
	histogram_energyFraction_vs_numPileup_->Fill(numPileup, weight_energyFraction);
	histogram_ptFraction_vs_numPileup_->Fill(numPileup, weight_ptFraction);
      }
    }
    void normalizeHistograms(double numEvents_processed)
    {
      if ( numEvents_processed > 0. ) 
      {
	divideByBinWidth(histogram_energyFraction_vs_eta_);
	histogram_energyFraction_vs_eta_->Scale(1./numEvents_processed);
	divideByBinWidth(histogram_ptFraction_vs_eta_);
	histogram_ptFraction_vs_eta_->Scale(1./numEvents_processed);
	divideByBinWidth(histogram_energyFraction_vs_eta_fine_binning_);
	histogram_energyFraction_vs_eta_fine_binning_->Scale(1./numEvents_processed);
	divideByBinWidth(histogram_ptFraction_vs_eta_fine_binning_);
	histogram_ptFraction_vs_eta_fine_binning_->Scale(1./numEvents_processed);
	divideByBinWidth(histogram_energyFraction_vs_absEta_);
	histogram_energyFraction_vs_absEta_->Scale(1./numEvents_processed);
	divideByBinWidth(histogram_ptFraction_vs_absEta_);
	histogram_ptFraction_vs_absEta_->Scale(1./numEvents_processed);
	divideByBinWidth(histogram_energyFraction_vs_pt_);
	histogram_energyFraction_vs_pt_->Scale(1./numEvents_processed);
	divideByBinWidth(histogram_ptFraction_vs_pt_);
	histogram_ptFraction_vs_pt_->Scale(1./numEvents_processed);
	histogram_energyFraction_vs_numPileup_->Scale(1./numEvents_processed);
	histogram_ptFraction_vs_numPileup_->Scale(1./numEvents_processed);
      }
    }
    std::string label_;
    double min_pt_;
    double max_absEta_;
    MonitorElement* me_energyFraction_vs_eta_;
    TH1* histogram_energyFraction_vs_eta_;
    MonitorElement* me_ptFraction_vs_eta_;
    TH1* histogram_ptFraction_vs_eta_;
    MonitorElement* me_energyFraction_vs_eta_fine_binning_;
    TH1* histogram_energyFraction_vs_eta_fine_binning_;
    MonitorElement* me_ptFraction_vs_eta_fine_binning_;
    TH1* histogram_ptFraction_vs_eta_fine_binning_;
    MonitorElement* me_energyFraction_vs_absEta_;
    TH1* histogram_energyFraction_vs_absEta_;
    MonitorElement* me_ptFraction_vs_absEta_;
    TH1* histogram_ptFraction_vs_absEta_;
    MonitorElement* me_energyFraction_vs_pt_;
    TH1* histogram_energyFraction_vs_pt_;
    MonitorElement* me_ptFraction_vs_pt_;
    TH1* histogram_ptFraction_vs_pt_;
    MonitorElement* me_energyFraction_vs_numPileup_;
    TH1* histogram_energyFraction_vs_numPileup_;
    MonitorElement* me_ptFraction_vs_numPileup_;
    TH1* histogram_ptFraction_vs_numPileup_;
  };
  pfCandTypePlotEntryType* pfChargedHadronPlots_;
  pfCandTypePlotEntryType* pfChargedHadronPileupPlots_;
  pfCandTypePlotEntryType* pfElectronPlots_;
  pfCandTypePlotEntryType* pfNeutralHadronPlots_;
  pfCandTypePlotEntryType* pfPhotonPlots_;
  pfCandTypePlotEntryType* pfMuonPlots_;
  
  MonitorElement* me_EventCounter_;
  TH1* histogram_EventCounter_;
};

#endif   

