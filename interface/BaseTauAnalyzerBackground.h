#ifndef HLTrigger_TallinnHLTPFTauAnalyzer_TauAnalyzerBackgroundBase_h
#define HLTrigger_TallinnHLTPFTauAnalyzer_TauAnalyzerBackgroundBase_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "HLTrigger/TallinnHLTPFTauAnalyzer/interface/BaseTau.h"
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

class BaseTauAnalyzerBackground : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit BaseTauAnalyzerBackground(const edm::ParameterSet&);
    
  // destructor
  virtual ~BaseTauAnalyzerBackground();
    
 protected:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  virtual std::vector<BaseTau> buildBaseTaus(const edm::Event&, const edm::EventSetup&) = 0;

  std::string moduleLabel_;

  double min_pt_;
  double max_pt_;
  std::vector<double> min_absEtaValues_;
  std::vector<double> max_absEtaValues_;
  size_t num_absEtaValues_;
  std::vector<double> min_leadTrackPtValues_;
  std::vector<double> max_leadTrackPtValues_;
  size_t num_leadTrackPtValues_;
  std::vector<double> min_relDiscriminatorValues_;
  std::vector<double> max_relDiscriminatorValues_;
  size_t num_relDiscriminatorValues_;
  std::vector<double> min_absDiscriminatorValues_;
  std::vector<double> max_absDiscriminatorValues_;
  size_t num_absDiscriminatorValues_;
  std::vector<double> min_dzValues_;
  std::vector<double> max_dzValues_;
  size_t num_dzValues_;

  edm::InputTag src_evtWeight_;
  edm::EDGetTokenT<double> token_evtWeight_;

  std::string dqmDirectory_;

 private:
  struct ratePlotEntryType
  {
    ratePlotEntryType(double min_absEta, double max_absEta, 
                      double min_leadTrackPt, double max_leadTrackPt, 
                      double min_relDiscriminator, double max_relDiscriminator, double min_absDiscriminator, double max_absDiscriminator, 
                      double min_dz, double max_dz)
      : me_numPFTaus_vs_ptThreshold_(nullptr)
      , histogram_numPFTaus_vs_ptThreshold_(nullptr)
      , min_absEta_(min_absEta)
      , max_absEta_(max_absEta)
      , min_leadTrackPt_(min_leadTrackPt)
      , max_leadTrackPt_(min_leadTrackPt)
      , min_relDiscriminator_(min_relDiscriminator)
      , max_relDiscriminator_(max_relDiscriminator)
      , min_absDiscriminator_(min_absDiscriminator)
      , max_absDiscriminator_(max_absDiscriminator)
      , min_dz_(min_dz)
      , max_dz_(max_dz)
    {}
    ~ratePlotEntryType() 
    {}
    void bookHistograms(DQMStore& dqmStore)
    {
      TString histogramName_suffix;
      histogramName_suffix = getHistogramName(histogramName_suffix, "absEta", min_absEta_, max_absEta_, 2);
      histogramName_suffix = getHistogramName(histogramName_suffix, "leadTrackPt", min_leadTrackPt_, max_leadTrackPt_, 0);
      histogramName_suffix = getHistogramName(histogramName_suffix, "relDiscriminator", min_relDiscriminator_, max_relDiscriminator_, 3);
      histogramName_suffix = getHistogramName(histogramName_suffix, "absDiscriminator", min_absDiscriminator_, max_absDiscriminator_, 3);
      if ( !(min_relDiscriminator_ >= 0. || max_relDiscriminator_ > 0. || min_absDiscriminator_ >= 0. || max_absDiscriminator_ > 0.) )
      {
        histogramName_suffix.Append("_noIsolation");
      }
      histogramName_suffix = getHistogramName(histogramName_suffix, "dz", min_dz_, max_dz_, 1);
      histogramName_suffix = histogramName_suffix.ReplaceAll(".", "p");

      TString histogramName_numPFTaus_vs_ptThreshold = Form("numPFTaus_vs_ptThreshold%s", histogramName_suffix.Data());
      me_numPFTaus_vs_ptThreshold_ = dqmStore.book2D(
        histogramName_numPFTaus_vs_ptThreshold.Data(), histogramName_numPFTaus_vs_ptThreshold.Data(), 
        251, -0.5, 250.5, 11, -0.5, +10.5);
      histogram_numPFTaus_vs_ptThreshold_ = dynamic_cast<TH2*>(me_numPFTaus_vs_ptThreshold_->getTH1());
      assert(histogram_numPFTaus_vs_ptThreshold_);
    }
    void fillHistograms(const std::vector<BaseTau>& baseTaus, double evtWeight)
    {
      std::vector<const BaseTau*> baseTaus_passingAbsEta;
      for ( std::vector<BaseTau>::const_iterator baseTau = baseTaus.begin();
            baseTau != baseTaus.end(); ++baseTau ) {
        double baseTau_absEta = TMath::Abs(baseTau->p4().eta());
        if ( (min_absEta_           < 0. || baseTau_absEta              >=  min_absEta_                              ) &&
	     (max_absEta_           < 0. || baseTau_absEta              <=  max_absEta_                              ) &&
             (min_leadTrackPt_      < 0. || baseTau->leadTrackP4().pt() >=  min_leadTrackPt_                         ) &&
             (max_leadTrackPt_      < 0. || baseTau->leadTrackP4().pt() >=  max_leadTrackPt_                         ) &&
             (min_relDiscriminator_ < 0. || baseTau->discriminator()    >= (min_relDiscriminator_*baseTau->p4().pt())) &&
             (max_relDiscriminator_ < 0. || baseTau->discriminator()    <= (max_relDiscriminator_*baseTau->p4().pt())) &&
             (min_absDiscriminator_ < 0. || baseTau->discriminator()    >=  min_absDiscriminator_                    ) &&
             (max_absDiscriminator_ < 0. || baseTau->discriminator()    <=  max_absDiscriminator_                    ) )
        {
	  baseTaus_passingAbsEta.push_back(&(*baseTau));
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
	  for ( std::vector<const BaseTau*>::const_iterator baseTau_zVtxRef = baseTaus_passingAbsEta.begin(); 
	        baseTau_zVtxRef != baseTaus_passingAbsEta.end(); ++baseTau_zVtxRef ) {
	    int numPFTaus_passingPt = 0;
	    for ( std::vector<const BaseTau*>::const_iterator baseTau_toMatch = baseTaus_passingAbsEta.begin(); 
		  baseTau_toMatch != baseTaus_passingAbsEta.end(); ++baseTau_toMatch ) {
	      if ( (*baseTau_toMatch)->p4().pt() > ptThreshold )
	      {
                double dz = TMath::Abs((*baseTau_zVtxRef)->zVtx() - (*baseTau_toMatch)->zVtx());
		if ( dz < max_dz_ ) 
	        {
		  ++numPFTaus_passingPt;
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
    double max_leadTrackPt_;
    double min_relDiscriminator_;
    double max_relDiscriminator_;
    double min_absDiscriminator_;
    double max_absDiscriminator_;
    double min_dz_;
    double max_dz_;
  };
  std::vector<ratePlotEntryType*> ratePlots_;
};

#endif   
