#ifndef HLTTrigger_TallinnHLTPFTauAnalyzer_RecoPFTauAnalyzerBackground_h
#define HLTTrigger_TallinnHLTPFTauAnalyzer_RecoPFTauAnalyzerBackground_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/TauReco/interface/PFTau.h"                                // reco::PFTau
#include "DataFormats/TauReco/interface/PFTauFwd.h"                             // reco::PFTauCollection
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"                   // reco::PFTauDiscriminator
#include "HLTTrigger/TallinnHLTPFTauAnalyzer/interface/histogramAuxFunctions.h" // fillWithOverFlow()

#include <TH1.h>     // TH1
#include <TH2.h>     // TH2
#include <TString.h> // TString, Form()
#include <TMath.h>   // TMath::Abs(), TMath::Nint()

#include <vector>    // std::vector
#include <string>    // std::string
#include <algorithm> // std::sort

using namespace dqm::implementation;

namespace
{
  bool
  isHigherPt(const std::pair<const reco::PFTau*, double>& pfTau_wChargedIso1,
	     const std::pair<const reco::PFTau*, double>& pfTau_wChargedIso2)
  {
    return pfTau_wChargedIso1.first->pt() > pfTau_wChargedIso2.first->pt();
  }
}

class RecoPFTauAnalyzerBackground : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit RecoPFTauAnalyzerBackground(const edm::ParameterSet&);
    
  // destructor
  ~RecoPFTauAnalyzerBackground();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag srcPFTaus_;
  edm::EDGetTokenT<reco::PFTauCollection> tokenPFTaus_;
  edm::InputTag srcPFTauSumChargedIso_;
  edm::EDGetTokenT<reco::PFTauDiscriminator> tokenPFTauSumChargedIso_;

  double min_pt_;
  double max_absEta_;

  double lumiScale_;

  std::string dqmDirectory_;

  struct ratePlotEntryType
  {
    ratePlotEntryType(double min_absEta, double max_absEta, double min_leadTrackPt, double max_relChargedIso, double max_absChargedIso, double max_dz)
      : me_numPFTaus_vs_ptThreshold_(nullptr)
      , histogram_numPFTaus_vs_ptThreshold_(nullptr)
      , me_numPFTausPtGt20_(nullptr)
      , histogram_numPFTausPtGt20_(nullptr)
      , me_numPFTausPtGt25_(nullptr)
      , histogram_numPFTausPtGt25_(nullptr)
      , me_numPFTausPtGt30_(nullptr)
      , histogram_numPFTausPtGt30_(nullptr)
      , me_numPFTausPtGt35_(nullptr)
      , histogram_numPFTausPtGt35_(nullptr)
      , me_numPFTausPtGt40_(nullptr)
      , histogram_numPFTausPtGt40_(nullptr)
      , me_numPFTausPtGt45_(nullptr)
      , histogram_numPFTausPtGt45_(nullptr)
      , me_numPFTausPtGt50_(nullptr)
      , histogram_numPFTausPtGt50_(nullptr)
      , min_absEta_(min_absEta)
      , max_absEta_(max_absEta)
      , min_leadTrackPt_(min_leadTrackPt)
      , max_relChargedIso_(max_relChargedIso)
      , max_absChargedIso_(max_absChargedIso)
      , max_dz_(max_dz)
    {}
    ~ratePlotEntryType() 
    {}
    void bookHistograms(DQMStore& dqmStore)
    {
      TString histogramName_suffix;
      if      ( min_absEta_ >= 0. && max_absEta_ > 0. ) histogramName_suffix.Append(Form("_absEta%1.2fto%1.2f", min_absEta_, max_absEta_));
      else if ( min_absEta_ >= 0.                     ) histogramName_suffix.Append(Form("_absEtaGt%1.2f", min_absEta_));
      else if (                      max_absEta_ > 0. ) histogramName_suffix.Append(Form("_absEtaLt%1.2f", max_absEta_));
      if ( min_leadTrackPt_   > 0. ) histogramName_suffix.Append(Form("_leadTrackPtGt%1.0f", min_leadTrackPt_));
      if ( max_relChargedIso_ > 0. || max_absChargedIso_ > 0. )
      {
        if ( max_relChargedIso_ > 0. ) histogramName_suffix.Append(Form("_relChargedIsoLt%1.2f", max_relChargedIso_));
        if ( max_absChargedIso_ > 0. ) histogramName_suffix.Append(Form("_absChargedIsoLt%1.2f", max_absChargedIso_));
      } 
      else
      {
        histogramName_suffix.Append("_noIsolation");
      }
      histogramName_suffix = histogramName_suffix.ReplaceAll(".", "p");

      TString histogramName_numPFTaus_vs_ptThreshold = Form("numPFTaus_vs_ptThreshold%s", histogramName_suffix.Data());
      me_numPFTaus_vs_ptThreshold_ = dqmStore.book2D(histogramName_numPFTaus_vs_ptThreshold.Data(), histogramName_numPFTaus_vs_ptThreshold.Data(), 251, -0.5, 250.5, 11, -0.5, +10.5);
      histogram_numPFTaus_vs_ptThreshold_ = dynamic_cast<TH2*>(me_numPFTaus_vs_ptThreshold_->getTH1());
      assert(histogram_numPFTaus_vs_ptThreshold_);

      TString histogramName_numPFTausPtGt20 = Form("numPFTausPtGt20%s", histogramName_suffix.Data());
      me_numPFTausPtGt20_ = dqmStore.book1D(histogramName_numPFTausPtGt20.Data(), histogramName_numPFTausPtGt20.Data(), 21, -0.5, +20.5);
      histogram_numPFTausPtGt20_ = me_numPFTausPtGt20_->getTH1();
      assert(histogram_numPFTausPtGt20_);
      TString histogramName_numPFTausPtGt25 = Form("numPFTausPtGt25%s", histogramName_suffix.Data());
      me_numPFTausPtGt25_ = dqmStore.book1D(histogramName_numPFTausPtGt25.Data(), histogramName_numPFTausPtGt25.Data(), 21, -0.5, +20.5);
      histogram_numPFTausPtGt25_ = me_numPFTausPtGt25_->getTH1();
      assert(histogram_numPFTausPtGt25_);
      TString histogramName_numPFTausPtGt30 = Form("numPFTausPtGt30%s", histogramName_suffix.Data());
      me_numPFTausPtGt30_ = dqmStore.book1D(histogramName_numPFTausPtGt30.Data(), histogramName_numPFTausPtGt30.Data(), 21, -0.5, +20.5);
      histogram_numPFTausPtGt30_ = me_numPFTausPtGt30_->getTH1();
      assert(histogram_numPFTausPtGt30_);
      TString histogramName_numPFTausPtGt35 = Form("numPFTausPtGt35%s", histogramName_suffix.Data());
      me_numPFTausPtGt35_ = dqmStore.book1D(histogramName_numPFTausPtGt35.Data(), histogramName_numPFTausPtGt35.Data(), 21, -0.5, +20.5);
      histogram_numPFTausPtGt35_ = me_numPFTausPtGt35_->getTH1();
      assert(histogram_numPFTausPtGt35_);
      TString histogramName_numPFTausPtGt40 = Form("numPFTausPtGt40%s", histogramName_suffix.Data());
      me_numPFTausPtGt40_ = dqmStore.book1D(histogramName_numPFTausPtGt40.Data(), histogramName_numPFTausPtGt40.Data(), 21, -0.5, +20.5);
      histogram_numPFTausPtGt40_ = me_numPFTausPtGt40_->getTH1();
      assert(histogram_numPFTausPtGt40_);
      TString histogramName_numPFTausPtGt45 = Form("numPFTausPtGt45%s", histogramName_suffix.Data());
      me_numPFTausPtGt45_ = dqmStore.book1D(histogramName_numPFTausPtGt45.Data(), histogramName_numPFTausPtGt45.Data(), 21, -0.5, +20.5);
      histogram_numPFTausPtGt45_ = me_numPFTausPtGt45_->getTH1();
      assert(histogram_numPFTausPtGt45_);
      TString histogramName_numPFTausPtGt50 = Form("numPFTausPtGt50%s", histogramName_suffix.Data());
      me_numPFTausPtGt50_ = dqmStore.book1D(histogramName_numPFTausPtGt50.Data(), histogramName_numPFTausPtGt50.Data(), 21, -0.5, +20.5);
      histogram_numPFTausPtGt50_ = me_numPFTausPtGt50_->getTH1();
      assert(histogram_numPFTausPtGt50_);
    }
    void fillHistograms(const std::vector<std::pair<const reco::PFTau*, double>>& pfTaus_wChargedIso, double evtWeight)
    {
      std::vector<const reco::PFTau*> pfTaus_passingAbsEta;
      for ( std::vector<std::pair<const reco::PFTau*, double>>::const_iterator pfTau_wChargedIso = pfTaus_wChargedIso.begin();
            pfTau_wChargedIso != pfTaus_wChargedIso.end(); ++pfTau_wChargedIso ) {
        const reco::PFTau* pfTau = pfTau_wChargedIso->first;
        double pfTau_absEta = TMath::Abs(pfTau->eta());
        if ( (min_absEta_        < 0. || pfTau_absEta                                      >=  min_absEta_                    ) &&
	     (max_absEta_        < 0. || pfTau_absEta                                      <=  max_absEta_                    ) &&
             (                           pfTau->leadPFChargedHadrCand().isNonnull()                                           &&   
                                         pfTau->leadPFChargedHadrCand()->bestTrack()                                          ) && 
             (min_leadTrackPt_   < 0. || pfTau->leadPFChargedHadrCand()->bestTrack()->pt() >=  min_leadTrackPt_               ) && 
	     (max_relChargedIso_ < 0. || pfTau_wChargedIso->second                         <= (max_relChargedIso_*pfTau->pt())) &&
	     (max_absChargedIso_ < 0. || pfTau_wChargedIso->second                         <=  max_absChargedIso_             ) )
	{
	  pfTaus_passingAbsEta.push_back(pfTau_wChargedIso->first);
	}
      }

      int numPFTausPtGt20 = 0;
      int numPFTausPtGt25 = 0;
      int numPFTausPtGt30 = 0;
      int numPFTausPtGt35 = 0;
      int numPFTausPtGt40 = 0;
      int numPFTausPtGt45 = 0;
      int numPFTausPtGt50 = 0;
      for ( std::vector<const reco::PFTau*>::const_iterator pfTau = pfTaus_passingAbsEta.begin();
	    pfTau != pfTaus_passingAbsEta.end(); ++pfTau ) {
	if ( (*pfTau)->pt() > 20. ) ++numPFTausPtGt20;
	if ( (*pfTau)->pt() > 25. ) ++numPFTausPtGt25;
	if ( (*pfTau)->pt() > 30. ) ++numPFTausPtGt30;
	if ( (*pfTau)->pt() > 35. ) ++numPFTausPtGt35;
	if ( (*pfTau)->pt() > 40. ) ++numPFTausPtGt40;
	if ( (*pfTau)->pt() > 45. ) ++numPFTausPtGt45;
	if ( (*pfTau)->pt() > 50. ) ++numPFTausPtGt50;
      }
      fillWithOverFlow(histogram_numPFTausPtGt20_, numPFTausPtGt20, evtWeight);
      fillWithOverFlow(histogram_numPFTausPtGt25_, numPFTausPtGt25, evtWeight);
      fillWithOverFlow(histogram_numPFTausPtGt30_, numPFTausPtGt30, evtWeight);
      fillWithOverFlow(histogram_numPFTausPtGt35_, numPFTausPtGt35, evtWeight);
      fillWithOverFlow(histogram_numPFTausPtGt40_, numPFTausPtGt40, evtWeight);
      fillWithOverFlow(histogram_numPFTausPtGt45_, numPFTausPtGt45, evtWeight);
      fillWithOverFlow(histogram_numPFTausPtGt50_, numPFTausPtGt50, evtWeight);

      TAxis* xAxis = histogram_numPFTaus_vs_ptThreshold_->GetXaxis();
      TAxis* yAxis = histogram_numPFTaus_vs_ptThreshold_->GetYaxis();
      int numBinsX = xAxis->GetNbins();
      bool max_numPFTaus_passingPt_isZero = false;
      for ( int idxBin = 1; idxBin <= numBinsX; ++idxBin )
      {
	double ptThreshold = xAxis->GetBinCenter(idxBin);

	int max_numPFTaus_passingPt = 0;
        if ( !max_numPFTaus_passingPt_isZero ) 
        {
	  for ( std::vector<const reco::PFTau*>::const_iterator pfTau_zVtxRef = pfTaus_passingAbsEta.begin(); 
	        pfTau_zVtxRef != pfTaus_passingAbsEta.end(); ++pfTau_zVtxRef ) {
	    int numPFTaus_passingPt = 0;
	    for ( std::vector<const reco::PFTau*>::const_iterator pfTau_toMatch = pfTaus_passingAbsEta.begin(); 
		  pfTau_toMatch != pfTaus_passingAbsEta.end(); ++pfTau_toMatch ) {
	      if ( (*pfTau_toMatch)->pt() > ptThreshold )
	      {
                if ( (*pfTau_zVtxRef)->leadPFChargedHadrCand().isNonnull() && (*pfTau_zVtxRef)->leadPFChargedHadrCand()->bestTrack() &&
                     (*pfTau_toMatch)->leadPFChargedHadrCand().isNonnull() && (*pfTau_toMatch)->leadPFChargedHadrCand()->bestTrack() )
	        {  
                  //double dz = TMath::Abs((*pfTau_zVtxRef)->leadPFChargedHadrCand()->vertex().z() - (*pfTau_toMatch)->leadPFChargedHadrCand()->vertex().z());
                  double dz = TMath::Abs((*pfTau_zVtxRef)->leadPFChargedHadrCand()->bestTrack()->vertex().z() - (*pfTau_toMatch)->leadPFChargedHadrCand()->bestTrack()->vertex().z());
		  if ( dz < max_dz_ ) 
	          {
		    ++numPFTaus_passingPt;
		  }
	        }
	      }
	    }
	    if ( numPFTaus_passingPt > max_numPFTaus_passingPt ) 
	    {
              max_numPFTaus_passingPt = numPFTaus_passingPt;
            }
	  }
	}

	if ( max_numPFTaus_passingPt > yAxis->GetXmax() ) max_numPFTaus_passingPt = TMath::Nint(yAxis->GetXmax() - 0.5);
	histogram_numPFTaus_vs_ptThreshold_->Fill(ptThreshold, max_numPFTaus_passingPt, evtWeight);

	if ( max_numPFTaus_passingPt == 0 ) 
	{
	  max_numPFTaus_passingPt_isZero = true;
	}
      }
    }
    void scaleHistograms(double sf)
    {
      histogram_numPFTaus_vs_ptThreshold_->Scale(sf);
      histogram_numPFTausPtGt20_->Scale(sf);
      histogram_numPFTausPtGt25_->Scale(sf);
      histogram_numPFTausPtGt30_->Scale(sf);
      histogram_numPFTausPtGt35_->Scale(sf);
      histogram_numPFTausPtGt40_->Scale(sf);   
      histogram_numPFTausPtGt45_->Scale(sf);
      histogram_numPFTausPtGt50_->Scale(sf);  
    }
    MonitorElement* me_numPFTaus_vs_ptThreshold_;
    TH2* histogram_numPFTaus_vs_ptThreshold_;
    MonitorElement* me_numPFTausPtGt20_;
    TH1* histogram_numPFTausPtGt20_;
    MonitorElement* me_numPFTausPtGt25_;
    TH1* histogram_numPFTausPtGt25_;
    MonitorElement* me_numPFTausPtGt30_;
    TH1* histogram_numPFTausPtGt30_;
    MonitorElement* me_numPFTausPtGt35_;
    TH1* histogram_numPFTausPtGt35_;
    MonitorElement* me_numPFTausPtGt40_;
    TH1* histogram_numPFTausPtGt40_;   
    MonitorElement* me_numPFTausPtGt45_;
    TH1* histogram_numPFTausPtGt45_;
    MonitorElement* me_numPFTausPtGt50_;
    TH1* histogram_numPFTausPtGt50_;  
    double min_absEta_;    
    double max_absEta_;    
    double min_leadTrackPt_;
    double max_relChargedIso_;
    double max_absChargedIso_;
    double max_dz_;
  };
  std::vector<ratePlotEntryType*> ratePlots_;

  unsigned long numEvents_processed_;
};

#endif   
