#ifndef HLTTrigger_TallinnHLTPFTauAnalyzer_L1VertexAnalyzer_h
#define HLTTrigger_TallinnHLTPFTauAnalyzer_L1VertexAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/VertexReco/interface/Vertex.h"                                 // reco::Vertex
#include "DataFormats/VertexReco/interface/VertexFwd.h"                              // reco::VertexCollection
#include "DataFormats/TrackReco/interface/Track.h"                                   // reco::Track
#include "DataFormats/TrackReco/interface/TrackFwd.h"                                // reco::TrackCollection
#include "DataFormats/TrackReco/interface/TrackBase.h"                               // reco::TrackBase::Point
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"                 // reco::PFCandidate
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"              // reco::PFCandidateCollection  
#include "DataFormats/JetReco/interface/GenJet.h"                                    // reco::GenJet
#include "DataFormats/JetReco/interface/GenJetCollection.h"                          // reco::GenJetCollection
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"                        // reco::GenParticle
#include "DataFormats/Math/interface/deltaR.h"                                       // deltaR

#include "HLTTrigger/TallinnHLTPFTauAnalyzer/interface/GenChargedHadronToTrackMatch.h" // GenChargedHadronToRecoTrackMatch

#include <TH1.h>
#include <TString.h> // TString, Form()
#include <TMath.h> // TMath::Abs()

#include <vector>
#include <string>

using namespace dqm::implementation;

namespace 
{
  // auxiliary struct for generator-level charged hadrons produced in tau decays
  class GenChargedHadron_and_genTau_decayMode
  {
   public:
    GenChargedHadron_and_genTau_decayMode(const reco::GenParticle* genChargedHadron, const std::string& genTau_decayMode)
      : genChargedHadron_(genChargedHadron)
      , genTau_decayMode_(genTau_decayMode)
    {
      assert(genChargedHadron_);
    }
    ~GenChargedHadron_and_genTau_decayMode()
    {}

    const reco::GenParticle* genChargedHadron()        const { return genChargedHadron_;                    }          
    double                   genChargedHadron_pt()     const { return genChargedHadron_->pt();              }
    double                   genChargedHadron_eta()    const { return genChargedHadron_->eta();             }
    double                   genChargedHadron_absEta() const { return TMath::Abs(genChargedHadron_->eta()); }
    double                   genChargedHadron_phi()    const { return genChargedHadron_->phi();             }
    const std::string&       genTau_decayMode()        const { return genTau_decayMode_;                    }

   private:
    const reco::GenParticle* genChargedHadron_;
    std::string genTau_decayMode_;
  };

  // auxiliary struct for matches between generator-level charged hadrons produced in tau decays and tracks reconstructed offline or on HLT level
  template <class T1, class T2>
  class GenChargedHadronToTrackMatch_and_genTau_decayMode
  {
   public:
    GenChargedHadronToTrackMatch_and_genTau_decayMode(const T2& genChargedHadronToTrackMatch, const std::string& genTau_decayMode, double dR)
      : genChargedHadronToTrackMatch_(genChargedHadronToTrackMatch)
      , genTau_decayMode_(genTau_decayMode)
      , dR_(dR)
    {}
    ~GenChargedHadronToTrackMatch_and_genTau_decayMode()
    {}

    const reco::Candidate* genChargedHadron()             const { return genChargedHadronToTrackMatch_.genChargedHadron();        }
    bool                   hasGenChargedHadron()          const { return genChargedHadronToTrackMatch_.hasGenChargedHadron();     }
    double                 genChargedHadron_pt()          const { return genChargedHadronToTrackMatch_.genChargedHadron_pt();     }
    double                 genChargedHadron_eta()         const { return genChargedHadronToTrackMatch_.genChargedHadron_eta();    }
    double                 genChargedHadron_absEta()      const { return genChargedHadronToTrackMatch_.genChargedHadron_absEta(); }
    double                 genChargedHadron_phi()         const { return genChargedHadronToTrackMatch_.genChargedHadron_phi();    }
    const T1*              recTrack()                     const { return genChargedHadronToTrackMatch_.recTrack();                }
    bool                   hasRecTrack()                  const { return genChargedHadronToTrackMatch_.hasRecTrack();             }
    double                 recTrack_pt()                  const { return genChargedHadronToTrackMatch_.recTrack_pt();             }
    double                 recTrack_eta()                 const { return genChargedHadronToTrackMatch_.recTrack_eta();            }
    double                 recTrack_absEta()              const { return genChargedHadronToTrackMatch_.recTrack_absEta();         }
    double                 recTrack_phi()                 const { return genChargedHadronToTrackMatch_.recTrack_phi();            }
    const std::string&     genTau_decayMode()             const { return genTau_decayMode_;                                       }
    double                 dR()                           const { return dR_;                                                     }
    const T2&              genChargedHadronToTrackMatch() const { return genChargedHadronToTrackMatch_;                           }
    
  private:
    T2 genChargedHadronToTrackMatch_;
    std::string genTau_decayMode_;
    double dR_; // used to sort (rank) matches
  };
  
  typedef GenChargedHadronToTrackMatch_and_genTau_decayMode<reco::Track, GenChargedHadronToRecoTrackMatch> GenChargedHadronToRecoTrackMatch_and_genTau_decayMode;

  template <class T>
  void printGenChargedHadronToTrackMatches(const std::vector<const T*>& genChargedHadronToTrackMatches)
  {
    size_t idx = 0;
    for ( auto match : genChargedHadronToTrackMatches )
    {
      std::cout << " genChargedHadron (#" << idx << "):" 
		<< " pT = " << match->genChargedHadron_pt() << ","
		<< " eta = " << match->genChargedHadron_eta() << ","
		<< " phi = " << match->genChargedHadron_phi() << ","
		<< " which is matched to recTrack:";
      if ( match->hasRecTrack() )
      {
	std::cout << " pT = " << match->recTrack_pt() << ","
		  << " eta = " << match->recTrack_eta() << ","
		  << " phi = " << match->recTrack_phi() 
		  << " (dR = " << match->dR() << ")";
      }
      else
      {
	std::cout << " N/A";
      }
      std::cout << std::endl;
      ++idx;
    }
  }

  template <class T>
  void printGenChargedHadronToTrackMatches(const std::vector<T>& genChargedHadronToTrackMatches)
  {
    std::vector<const T*> genChargedHadronToTrackMatches_ptr;
    for ( typename std::vector<T>::const_iterator match = genChargedHadronToTrackMatches.begin();
	  match != genChargedHadronToTrackMatches.end(); ++match )
    {
      genChargedHadronToTrackMatches_ptr.push_back(&(*match));
    }
    printGenChargedHadronToTrackMatches(genChargedHadronToTrackMatches_ptr);
  }
}

class RecoTrackAnalyzer : public edm::EDAnalyzer 
{
 public:
  RecoTrackAnalyzer(const edm::ParameterSet& cfg);
  ~RecoTrackAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag src_genTaus_;
  edm::EDGetTokenT<reco::GenJetCollection> token_genTaus_;

  enum { kGenVtx, kRecVtx };
  int vtxMode_;

  edm::InputTag src_genVertex_position_;
  edm::EDGetTokenT<reco::TrackBase::Point> token_genVertex_position_;

  edm::InputTag src_offlineVertices_;
  edm::EDGetTokenT<reco::VertexCollection> token_offlineVertices_;

  edm::InputTag src_offlineTracks_;
  edm::EDGetTokenT<reco::TrackCollection> token_offlineTracks_;

  edm::InputTag src_offlinePFCands_;
  edm::EDGetTokenT<reco::PFCandidateCollection> token_offlinePFCands_;

  edm::InputTag src_hltVertices_;
  edm::EDGetTokenT<reco::VertexCollection> token_hltVertices_;

  edm::InputTag src_hltTracks_;
  edm::EDGetTokenT<reco::TrackCollection> token_hltTracks_;

  edm::InputTag src_hltPFCands_;
  edm::EDGetTokenT<reco::PFCandidateCollection> token_hltPFCands_;

  std::string dqmDirectory_;

  struct efficiencyPlotEntryType
  {
    efficiencyPlotEntryType(const std::string& recTrack_type, 
			    double genChargedHadron_min_pt, double genChargedHadron_max_pt, double genChargedHadron_min_absEta, double genChargedHadron_max_absEta, 
			    const std::string& genTau_decayMode)
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
      , me_resolution_(nullptr)
      , histogram_resolution_(nullptr)
      , recTrack_type_(recTrack_type)
      , genChargedHadron_min_pt_(genChargedHadron_min_pt)
      , genChargedHadron_max_pt_(genChargedHadron_max_pt)
      , genChargedHadron_min_absEta_(genChargedHadron_min_absEta)
      , genChargedHadron_max_absEta_(genChargedHadron_max_absEta)
      , genTau_decayMode_(genTau_decayMode)
    {
      if ( recTrack_type_.length() > 0 ) 
      {	
	recTrack_type_capitalized_ = recTrack_type_;
	recTrack_type_capitalized_[0] = toupper(recTrack_type_capitalized_[0]);	
      }
      if ( genTau_decayMode_.length() > 0 ) 
      {	
	genTau_decayMode_capitalized_ = genTau_decayMode_;
	genTau_decayMode_capitalized_[0] = toupper(genTau_decayMode_capitalized_[0]);	
      }
    }
    ~efficiencyPlotEntryType() 
    {}
    void bookHistograms(DQMStore& dqmStore)
    {
      TString histogramName_suffix;
      if ( genChargedHadron_min_absEta_ >= 0. && genChargedHadron_max_absEta_ > 0. ) 
      { 
	histogramName_suffix.Append(Form("_absEta%1.2fto%1.2f", genChargedHadron_min_absEta_, genChargedHadron_max_absEta_));
      }
      else if ( genChargedHadron_min_absEta_ >= 0. ) 
      {
	histogramName_suffix.Append(Form("_absEtaGt%1.2f", genChargedHadron_min_absEta_));
      }
      else if ( genChargedHadron_max_absEta_ > 0. ) 
      {
	histogramName_suffix.Append(Form("_absEtaLt%1.2f", genChargedHadron_max_absEta_));
      }
      if ( genChargedHadron_min_pt_ >= 0. && genChargedHadron_max_pt_ > 0. ) 
      {
	histogramName_suffix.Append(Form("_pt%1.2fto%1.2f", genChargedHadron_min_pt_, genChargedHadron_max_pt_));
      }
      else if ( genChargedHadron_min_pt_ >= 0. ) 
      {
	histogramName_suffix.Append(Form("_ptGt%1.2f", genChargedHadron_min_pt_));
      }
      else if ( genChargedHadron_max_pt_ > 0. ) 
      {
	histogramName_suffix.Append(Form("_ptLt%1.2f", genChargedHadron_max_pt_));
      }
      if ( genTau_decayMode_ != "all" ) 
      {
	histogramName_suffix.Append(Form("_gen%sTau", genTau_decayMode_capitalized_.data()));      
      }
      histogramName_suffix = histogramName_suffix.ReplaceAll(".", "p");

      const int numBins_ptL = 19;
      float binning_ptL[numBins_ptL + 1] = { 
        0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12.5, 15., 17.5, 20., 25., 30., 40., 60., 100.
      };

      TString histogramName_pt_numerator = Form("eff%s_vs_pt_numerator%s", recTrack_type_capitalized_.data(), histogramName_suffix.Data());
      me_pt_numerator_ = dqmStore.book1D(histogramName_pt_numerator.Data(), histogramName_pt_numerator.Data(), numBins_ptL, binning_ptL);
      histogram_pt_numerator_ = me_pt_numerator_->getTH1();
      assert(histogram_pt_numerator_);
      TString histogramName_pt_denominator = Form("eff%s_vs_pt_denominator%s", recTrack_type_capitalized_.data(), histogramName_suffix.Data());
      me_pt_denominator_ = dqmStore.book1D(histogramName_pt_denominator.Data(), histogramName_pt_denominator.Data(), numBins_ptL, binning_ptL);
      histogram_pt_denominator_ = me_pt_denominator_->getTH1();
      assert(histogram_pt_denominator_);

      TString histogramName_eta_numerator = Form("eff%s_vs_eta_numerator%s", recTrack_type_capitalized_.data(), histogramName_suffix.Data());
      me_eta_numerator_ = dqmStore.book1D(histogramName_eta_numerator.Data(), histogramName_eta_numerator.Data(), 30, -3., +3.);
      histogram_eta_numerator_ = me_eta_numerator_->getTH1();
      assert(histogram_eta_numerator_);
      TString histogramName_eta_denominator = Form("eff%s_vs_eta_denominator%s", recTrack_type_capitalized_.data(), histogramName_suffix.Data());
      me_eta_denominator_ = dqmStore.book1D(histogramName_eta_denominator.Data(), histogramName_eta_denominator.Data(), 30, -3., +3.);
      histogram_eta_denominator_ = me_eta_denominator_->getTH1();
      assert(histogram_eta_denominator_);

      TString histogramName_phi_numerator = Form("eff%s_vs_phi_numerator%s", recTrack_type_capitalized_.data(), histogramName_suffix.Data());
      me_phi_numerator_ = dqmStore.book1D(histogramName_phi_numerator.Data(), histogramName_phi_numerator.Data(), 18, -TMath::Pi(), +TMath::Pi());
      histogram_phi_numerator_ = me_phi_numerator_->getTH1();
      assert(histogram_phi_numerator_);
      TString histogramName_phi_denominator = Form("eff%s_vs_phi_denominator%s", recTrack_type_capitalized_.data(), histogramName_suffix.Data());
      me_phi_denominator_ = dqmStore.book1D(histogramName_phi_denominator.Data(), histogramName_phi_denominator.Data(), 18, -TMath::Pi(), +TMath::Pi());
      histogram_phi_denominator_ = me_phi_denominator_->getTH1();
      assert(histogram_phi_denominator_);

      TString histogramName_minDeltaR_numerator = Form("eff%s_vs_minDeltaR_numerator%s", recTrack_type_capitalized_.data(), histogramName_suffix.Data());
      me_minDeltaR_numerator_ = dqmStore.book1D(histogramName_minDeltaR_numerator.Data(), histogramName_minDeltaR_numerator.Data(), 200, 0., 0.05);
      histogram_minDeltaR_numerator_ = me_minDeltaR_numerator_->getTH1();
      assert(histogram_minDeltaR_numerator_);
      TString histogramName_minDeltaR_denominator = Form("eff%s_vs_minDeltaR_denominator%s", recTrack_type_capitalized_.data(), histogramName_suffix.Data());
      me_minDeltaR_denominator_ = dqmStore.book1D(histogramName_minDeltaR_denominator.Data(), histogramName_minDeltaR_denominator.Data(), 200, 0., 0.05);
      histogram_minDeltaR_denominator_ = me_minDeltaR_denominator_->getTH1();
      assert(histogram_minDeltaR_denominator_);

      const int numBins_absEta = 6;
      float binning_absEta[numBins_absEta + 1] = { 
        0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4
      };

      const int numBins_ptS = 11;
      float binning_ptS[numBins_ptS + 1] = { 
        0., 2., 4., 6., 8., 10., 15., 20., 30., 40., 60., 100.
      };

      TString histogramName_pt_vs_absEta_numerator = Form("eff%s_vs_pt_vs_absEta_numerator%s", recTrack_type_capitalized_.data(), histogramName_suffix.Data());
      me_pt_vs_absEta_numerator_ = dqmStore.book2D(histogramName_pt_vs_absEta_numerator.Data(), histogramName_pt_vs_absEta_numerator.Data(), numBins_absEta, binning_absEta, numBins_ptS, binning_ptS);
      histogram_pt_vs_absEta_numerator_ = dynamic_cast<TH2*>(me_pt_vs_absEta_numerator_->getTH1());
      assert(histogram_pt_vs_absEta_numerator_);
      TString histogramName_pt_vs_absEta_denominator = Form("eff%s_vs_pt_vs_absEta_denominator%s", recTrack_type_capitalized_.data(), histogramName_suffix.Data());
      me_pt_vs_absEta_denominator_ = dqmStore.book2D(histogramName_pt_vs_absEta_denominator.Data(), histogramName_pt_vs_absEta_denominator.Data(), numBins_absEta, binning_absEta, numBins_ptS, binning_ptS);
      histogram_pt_vs_absEta_denominator_ = dynamic_cast<TH2*>(me_pt_vs_absEta_denominator_->getTH1());
      assert(histogram_pt_vs_absEta_denominator_);

      TString histogramName_resolution = Form("resolution%s%s", recTrack_type_capitalized_.data(), histogramName_suffix.Data());
      me_resolution_ = dqmStore.book1D(histogramName_resolution.Data(), histogramName_resolution.Data(), 100, -0.05, +0.05);  
      histogram_resolution_ = me_resolution_->getTH1();
      assert(histogram_resolution_);
    }
    template <class T>
    void fillHistograms(const std::vector<const T*>& cleanedGenChargedHadronToTrackMatches, double evtWeight)
    {
      //std::cout << "<efficiencyPlotEntryType::fillHistograms (evtWeight = " << evtWeight << ")>:" << std::endl;

      double minDeltaR = 1.e+3;
      for ( typename std::vector<const T*>::const_iterator match1 = cleanedGenChargedHadronToTrackMatches.begin();
	    match1 != cleanedGenChargedHadronToTrackMatches.end(); ++match1 ) {
	for ( typename std::vector<const T*>::const_iterator match2 = match1 + 1;
	      match2 != cleanedGenChargedHadronToTrackMatches.end(); ++match2 ) {
	  double dR = deltaR((*match1)->genChargedHadron_eta(), (*match1)->genChargedHadron_phi(), (*match2)->genChargedHadron_eta(), (*match2)->genChargedHadron_phi());
	  if ( dR > 0. && dR < minDeltaR ) minDeltaR = dR;
        }
      }

      //size_t idx = 0;
      for ( auto match : cleanedGenChargedHadronToTrackMatches )
      {
	//std::cout << " genChargedHadron (#" << idx << "):" 
	//	    << " pT = " << match->genChargedHadron_pt() << ","
	//	    << " eta = " << match->genChargedHadron_eta() << ","
	//	    << " phi = " << match->genChargedHadron_phi() << std::endl;
	//++idx;
	  
	if ( genTau_decayMode_ != "all" && match->genTau_decayMode() != genTau_decayMode_ ) continue;

	if ( match->genChargedHadron_pt()     > genChargedHadron_min_pt_     && match->genChargedHadron_pt()     < genChargedHadron_max_pt_     && 
	     match->genChargedHadron_absEta() > genChargedHadron_min_absEta_ && match->genChargedHadron_absEta() < genChargedHadron_max_absEta_ )
	{
	  if ( match->hasRecTrack() )
	  {
	    //std::cout << " filling numerator histograms" << std::endl; 
	    histogram_pt_numerator_->Fill(match->genChargedHadron_pt(), evtWeight);
	    histogram_eta_numerator_->Fill(match->genChargedHadron_eta(), evtWeight);
	    histogram_phi_numerator_->Fill(match->genChargedHadron_phi(), evtWeight);
	    histogram_minDeltaR_numerator_->Fill(minDeltaR, evtWeight);
	    histogram_pt_vs_absEta_numerator_->Fill(match->genChargedHadron_absEta(), match->genChargedHadron_pt(), evtWeight);
	  }

	  //std::cout << " filling denominator histograms" << std::endl; 
	  histogram_pt_denominator_->Fill(match->genChargedHadron_pt(), evtWeight);
	  histogram_eta_denominator_->Fill(match->genChargedHadron_eta(), evtWeight);
	  histogram_phi_denominator_->Fill(match->genChargedHadron_phi(), evtWeight);
	  histogram_minDeltaR_denominator_->Fill(minDeltaR, evtWeight);
	  histogram_pt_vs_absEta_denominator_->Fill(match->genChargedHadron_absEta(), match->genChargedHadron_pt(), evtWeight);

	  if ( match->hasRecTrack() )
  	  {
	    histogram_resolution_->Fill((1./match->recTrack_pt()) - (1./match->genChargedHadron_pt()), evtWeight);
	  }
	}
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
    MonitorElement* me_resolution_;
    TH1* histogram_resolution_;
    std::string recTrack_type_;
    std::string recTrack_type_capitalized_;
    double genChargedHadron_min_pt_;
    double genChargedHadron_max_pt_;
    double genChargedHadron_min_absEta_;
    double genChargedHadron_max_absEta_;
    std::string genTau_decayMode_;
    std::string genTau_decayMode_capitalized_;
  };
  std::vector<efficiencyPlotEntryType*> efficiencyPlots_offlineTracks_woQualityCuts_;
  std::vector<efficiencyPlotEntryType*> efficiencyPlots_offlineTracks_wQualityCuts_;
  std::vector<efficiencyPlotEntryType*> efficiencyPlots_offlinePFCandTracks_woQualityCuts_;
  std::vector<efficiencyPlotEntryType*> efficiencyPlots_offlinePFCandTracks_wQualityCuts_;
  std::vector<efficiencyPlotEntryType*> efficiencyPlots_hltTracks_woQualityCuts_;
  std::vector<efficiencyPlotEntryType*> efficiencyPlots_hltTracks_wQualityCuts_;
  std::vector<efficiencyPlotEntryType*> efficiencyPlots_hltPFCandTracks_woQualityCuts_;
  std::vector<efficiencyPlotEntryType*> efficiencyPlots_hltPFCandTracks_wQualityCuts_;

  bool debug_;
};

#endif   

