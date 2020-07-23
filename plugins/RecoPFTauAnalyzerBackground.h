#ifndef HLTrigger_TallinnHLTPFTauAnalyzer_RecoPFTauAnalyzerBackground_h
#define HLTrigger_TallinnHLTPFTauAnalyzer_RecoPFTauAnalyzerBackground_h

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
#include "HLTrigger/TallinnHLTPFTauAnalyzer/interface/histogramAuxFunctions.h" // fillWithOverFlow()

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
    }
    MonitorElement* me_numPFTaus_vs_ptThreshold_;
    TH2* histogram_numPFTaus_vs_ptThreshold_;
    double min_absEta_;    
    double max_absEta_;    
    double min_leadTrackPt_;
    double max_relChargedIso_;
    double max_absChargedIso_;
    double max_dz_;
  };
  std::vector<ratePlotEntryType*> ratePlots_;
};

#endif   
