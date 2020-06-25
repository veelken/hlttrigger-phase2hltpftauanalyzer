#ifndef HLTTrigger_TallinnHLTPFTauAnalyzer_RecoPFTauPairAnalyzer_h
#define HLTTrigger_TallinnHLTPFTauAnalyzer_RecoPFTauPairAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"                       // reco::CandidateView
#include "DataFormats/Candidate/interface/Candidate.h"                          // reco::Candidate
#include "DataFormats/Phase2HLTPFTaus/interface/PFTauPair.h"                    // reco::PFTauPair
#include "DataFormats/Phase2HLTPFTaus/interface/PFTauPairFwd.h"                 // reco::PFTauPairCollection
#include "HLTTrigger/TallinnHLTPFTauAnalyzer/interface/histogramAuxFunctions.h" // fillWithOverFlow()

#include <TH1.h>     // TH1
#include <TH2.h>     // TH2
#include <TString.h> // TString, Form()
#include <TMath.h>   // TMath::Abs(), TMath::Nint()

#include <vector>    // std::vector
#include <string>    // std::string
#include <algorithm> // std::sort

using namespace dqm::legacy;

class RecoPFTauPairAnalyzer : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit RecoPFTauPairAnalyzer(const edm::ParameterSet&);
    
  // destructor
  ~RecoPFTauPairAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag srcPFTauPairs_;
  edm::EDGetTokenT<reco::PFTauPairCollection> tokenPFTauPairs_;
  edm::InputTag srcRefTaus_;
  edm::EDGetTokenT<reco::CandidateView> tokenRefTaus_;

  double min_refTau_pt_;
  double max_refTau_pt_;
  double min_refTau_absEta_;
  double max_refTau_absEta_;
  double dRmatch_;

  double lumiScale_;

  std::string dqmDirectory_;

  struct efficiency_or_ratePlotEntryType
  {
    efficiency_or_ratePlotEntryType(double min_absEta, double max_absEta, double max_relChargedIso, double max_absChargedIso, double max_dz)
      : me_efficiency_or_rate_(nullptr)
      , histogram_efficiency_or_rate_(nullptr)
      , me_denominator_(nullptr)
      , histogram_denominator_(nullptr)
      , min_absEta_(min_absEta)
      , max_absEta_(max_absEta)
      , max_relChargedIso_(max_relChargedIso)
      , max_absChargedIso_(max_absChargedIso)
      , max_dz_(max_dz)
    {}
    ~efficiency_or_ratePlotEntryType() 
    {}
    void bookHistograms(DQMStore& dqmStore)
    {
      TString histogramName_suffix;
      if      ( min_absEta_ >= 0. && max_absEta_ > 0. ) histogramName_suffix.Append(Form("_absEta%1.2fto%1.2f", min_absEta_, max_absEta_));
      else if ( min_absEta_ >= 0.                     ) histogramName_suffix.Append(Form("_absEtaGt%1.2f", min_absEta_));
      else if (                      max_absEta_ > 0. ) histogramName_suffix.Append(Form("_absEtaLt%1.2f", max_absEta_));
      if ( max_relChargedIso_ > 0. ) histogramName_suffix.Append(Form("_relChargedIsoLt%1.2f", max_relChargedIso_));
      if ( max_absChargedIso_ > 0. ) histogramName_suffix.Append(Form("_absChargedIsoLt%1.2f", max_absChargedIso_));
      histogramName_suffix = histogramName_suffix.ReplaceAll(".", "p");

      TString histogramName_efficiency_or_rate = Form("efficiency_or_rate%s", histogramName_suffix.Data());
      me_efficiency_or_rate_ = dqmStore.book2D(histogramName_efficiency_or_rate.Data(), histogramName_efficiency_or_rate.Data(), 101, -0.5, 100.5, 101, -0.5, 100.5);
      histogram_efficiency_or_rate_ = dynamic_cast<TH2*>(me_efficiency_or_rate_->getTH1());
      assert(histogram_efficiency_or_rate_);

      TString histogramName_denominator = Form("denominator%s", histogramName_suffix.Data());
      me_denominator_ = dqmStore.book1D(histogramName_denominator.Data(), histogramName_denominator.Data(), 1, -0.5, +0.5);
      histogram_denominator_ = me_denominator_->getTH1();
      assert(histogram_denominator_);
    }
    void fillHistograms(const std::vector<const reco::PFTauPair*>& pfTauPairs, double evtWeight)
    {
      std::vector<const reco::PFTauPair*> pfTauPairs_passingAbsEta;
      for ( std::vector<const reco::PFTauPair*>::const_iterator pfTauPair = pfTauPairs.begin(); 
	    pfTauPair != pfTauPairs.end(); ++pfTauPair ) 
      {
	const reco::PFTauRef& leadPFTau = (*pfTauPair)->leadPFTau();
        double leadPFTau_sumChargedIso = (*pfTauPair)->leadPFTau_sumChargedIso();
	const reco::PFTauRef& subleadPFTau = (*pfTauPair)->subleadPFTau();
        double subleadPFTau_sumChargedIso = (*pfTauPair)->subleadPFTau_sumChargedIso();
	if ( ( min_absEta_ < 0.                                                         || 
	      (TMath::Abs(leadPFTau->eta())    >=  min_absEta_                            && 
	       TMath::Abs(subleadPFTau->eta()) >=  min_absEta_                          ) ) &&
	     ( max_absEta_ < 0.                                                         || 
	      (TMath::Abs(leadPFTau->eta())    <=  max_absEta_                            && 
	       TMath::Abs(subleadPFTau->eta()) <=  max_absEta_                          ) ) &&
	     (max_relChargedIso_ < 0. || 
	      (leadPFTau_sumChargedIso         <= (max_relChargedIso_*leadPFTau->pt()   ) &&
	       subleadPFTau_sumChargedIso      <= (max_relChargedIso_*subleadPFTau->pt()))) &&
	     (max_absChargedIso_ < 0. || 
	      (leadPFTau_sumChargedIso         <=  max_absChargedIso_                     &&
	       subleadPFTau_sumChargedIso      <=  max_absChargedIso_                    )) )
	{
	  if ( leadPFTau->leadPFChargedHadrCand().isNonnull()    && leadPFTau->leadPFChargedHadrCand()->bestTrack()    &&
               subleadPFTau->leadPFChargedHadrCand().isNonnull() && subleadPFTau->leadPFChargedHadrCand()->bestTrack() )
	  {
            double dz = (*pfTauPair)->dz();
            if ( max_dz_ < 0. || dz < max_dz_ ) 
	    {
	      pfTauPairs_passingAbsEta.push_back(*pfTauPair);
	    }
          }
	}
      }

      TAxis* xAxis = histogram_efficiency_or_rate_->GetXaxis();
      int numBinsX = xAxis->GetNbins();
      TAxis* yAxis = histogram_efficiency_or_rate_->GetYaxis();
      int numBinsY = yAxis->GetNbins();
      bool max_numPFTauPairs_passingPt_isZero = false;
      for ( int idxBinX = 1; idxBinX <= numBinsX; ++idxBinX )
      {
	double leadPFTau_ptThreshold = xAxis->GetBinCenter(idxBinX);
	int max_numPFTauPairs_passingPt = 0;
	if ( !max_numPFTauPairs_passingPt_isZero ) 
        {
	  for ( int idxBinY = 1; idxBinY <= numBinsY; ++idxBinY )
          {
	    double subleadPFTau_ptThreshold = yAxis->GetBinCenter(idxBinY);
	    int numPFTauPairs_passingPt = 0;
	    for ( std::vector<const reco::PFTauPair*>::const_iterator pfTauPair = pfTauPairs_passingAbsEta.begin(); 
		  pfTauPair != pfTauPairs_passingAbsEta.end(); ++pfTauPair ) {
	      if ( (*pfTauPair)->leadPFTau()->pt()    > leadPFTau_ptThreshold    &&
		   (*pfTauPair)->subleadPFTau()->pt() > subleadPFTau_ptThreshold )
	      {
	        ++numPFTauPairs_passingPt;
	      }
	    }
	    if ( numPFTauPairs_passingPt >= 1 )
	    {
	      histogram_efficiency_or_rate_->Fill(leadPFTau_ptThreshold, subleadPFTau_ptThreshold, evtWeight);
	    }
	    if ( numPFTauPairs_passingPt > max_numPFTauPairs_passingPt ) 
	    {
	      max_numPFTauPairs_passingPt = numPFTauPairs_passingPt;
	    }
  	  }
	  if ( max_numPFTauPairs_passingPt == 0 ) 
	  {
	    max_numPFTauPairs_passingPt_isZero = true;
	  }
	}
      }

      histogram_denominator_->Fill(0., evtWeight);
    }
    void normalizeHistograms()
    {
      if ( histogram_denominator_->Integral() > 0. ) 
      {
	histogram_efficiency_or_rate_->Scale(1./histogram_denominator_->Integral());
      }
    }
    MonitorElement* me_efficiency_or_rate_;
    TH2* histogram_efficiency_or_rate_;
    MonitorElement* me_denominator_;
    TH1* histogram_denominator_;
    double min_absEta_;    
    double max_absEta_;   
    double max_relChargedIso_;
    double max_absChargedIso_;
    double max_dz_;
  };
  std::vector<efficiency_or_ratePlotEntryType*> efficiency_or_ratePlots_;
};

#endif   
