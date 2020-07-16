#ifndef HLTrigger_HLTPFTauAnalyzer_RecoPFCandidateConeEnResponseAnalyzer_h
#define HLTrigger_HLTPFTauAnalyzer_RecoPFCandidateConeEnResponseAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"            // reco::PFCandidate
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"         // reco::PFCandidateCollection
#include "DataFormats/TrackReco/interface/Track.h"                              // reco::Track
#include "DataFormats/TrackReco/interface/TrackFwd.h"                           // reco::TrackCollection
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"                   // reco::GenParticle
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"                // reco::GenParticleCollection
#include "DataFormats/Common/interface/View.h"                                  // edm::View
#include "DataFormats/Candidate/interface/Candidate.h"                          // reco::Candidate
#include "DataFormats/JetReco/interface/GenJet.h"                               // reco::GenJet
#include "DataFormats/JetReco/interface/GenJetCollection.h"                     // reco::GenJetCollection
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"                         // JetMCTagUtils::genTauDecayMode()
#include "DataFormats/Math/interface/deltaR.h"                                  // reco::deltaR

#include "HLTrigger/TallinnHLTPFTauAnalyzer/interface/histogramAuxFunctions.h" // fillWithOverFlow2d

#include <TH1.h>     // TH1
#include <TH2.h>     // TH2
#include <TString.h> // TString, Form()
#include <TMath.h>   // TMath::Abs(), TMath::Pi()

#include <vector>    // std::vector
#include <string>    // std::string

using namespace dqm::implementation;

class RecoPFCandidateConeEnResponseAnalyzer : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit RecoPFCandidateConeEnResponseAnalyzer(const edm::ParameterSet&);
    
  // destructor
  ~RecoPFCandidateConeEnResponseAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag srcPFCandidates_;
  edm::EDGetTokenT<reco::PFCandidateCollection> tokenPFCandidates_;
  edm::InputTag srcTracks_;
  edm::EDGetTokenT<reco::TrackCollection> tokenTracks_;
  edm::InputTag srcGenParticles_;
  edm::EDGetTokenT<edm::View<reco::Candidate>> tokenGenParticles_;
  edm::InputTag srcGenHadTaus_;
  edm::EDGetTokenT<reco::GenJetCollection> tokenGenHadTaus_;

  std::vector<double> dRcones_;

  std::string dqmDirectory_;

  struct responsePlotEntryType
  {
    responsePlotEntryType(const std::string& decayMode, double dRcone)
      : me_sumChargedPFCands_(nullptr)
      , histogram_sumChargedPFCands_(nullptr)
      , me_sumAllPFCands_(nullptr)
      , histogram_sumAllPFCands_(nullptr)
      , me_leadChargedPFCand_(nullptr)
      , histogram_leadChargedPFCand_(nullptr)
      , me_sumTracks_(nullptr)
      , histogram_sumTracks_(nullptr)
      , me_leadTrack_(nullptr)
      , histogram_leadTrack_(nullptr)
      , decayMode_(decayMode)
      , dRcone_(dRcone)
    {}
    ~responsePlotEntryType() 
    {}
    void bookHistograms(DQMStore& dqmStore)
    {
      TString histogramName_suffix = Form("%s_dRconeEq%1.2f", decayMode_.data(), dRcone_);
      histogramName_suffix = histogramName_suffix.ReplaceAll(".", "p");

      TString histogramName_sumChargedPFCands = Form("sumChargedPFCands_%s", histogramName_suffix.Data());
      me_sumChargedPFCands_ = dqmStore.book2D(histogramName_sumChargedPFCands.Data(), histogramName_sumChargedPFCands.Data(), 24, -2.4, +2.4, 40, 0., 4.);
      histogram_sumChargedPFCands_ = dynamic_cast<TH2*>(me_sumChargedPFCands_->getTH1());
      assert(histogram_sumChargedPFCands_);
      TString histogramName_sumAllPFCands = Form("sumAllPFCands_%s", histogramName_suffix.Data());
      me_sumAllPFCands_ = dqmStore.book2D(histogramName_sumAllPFCands.Data(), histogramName_sumAllPFCands.Data(), 24, -2.4, +2.4, 40, 0., 4.);
      histogram_sumAllPFCands_ = dynamic_cast<TH2*>(me_sumAllPFCands_->getTH1());
      assert(histogram_sumAllPFCands_);
      TString histogramName_leadChargedPFCand = Form("leadChargedPFCand_%s", histogramName_suffix.Data());
      me_leadChargedPFCand_ = dqmStore.book2D(histogramName_leadChargedPFCand.Data(), histogramName_leadChargedPFCand.Data(), 24, -2.4, +2.4, 40, 0., 4.);
      histogram_leadChargedPFCand_ = dynamic_cast<TH2*>(me_leadChargedPFCand_->getTH1());
      assert(histogram_leadChargedPFCand_);

      TString histogramName_sumTracks = Form("sumTracks_%s", histogramName_suffix.Data());
      me_sumTracks_ = dqmStore.book2D(histogramName_sumTracks.Data(), histogramName_sumTracks.Data(), 24, -2.4, +2.4, 40, 0., 4.);
      histogram_sumTracks_ = dynamic_cast<TH2*>(me_sumTracks_->getTH1());
      assert(histogram_sumTracks_);
      TString histogramName_leadTrack = Form("leadTrack_%s", histogramName_suffix.Data());
      me_leadTrack_ = dqmStore.book2D(histogramName_leadTrack.Data(), histogramName_leadTrack.Data(), 24, -2.4, +2.4, 40, 0., 4.);
      histogram_leadTrack_ = dynamic_cast<TH2*>(me_leadTrack_->getTH1());
      assert(histogram_leadTrack_);
    }
    void fillHistograms(const reco::PFCandidateCollection& pfCands, const reco::TrackCollection& tracks, const edm::View<reco::Candidate>& genParticles,
                        const reco::GenJetCollection& genHadTaus, double evtWeight)
    {
      for ( reco::GenJetCollection::const_iterator genHadTau = genHadTaus.begin();
            genHadTau != genHadTaus.end(); ++genHadTau ) {
        std::string genHadTau_decayMode = JetMCTagUtils::genTauDecayMode(*genHadTau);
	if ( decayMode_ != "all" && genHadTau_decayMode != decayMode_ ) continue;
        if ( !(genHadTau->pt() > 20. && TMath::Abs(genHadTau->eta()) < 2.4) ) continue;

        double sumChargedPFCandPt  = 0.;
        double sumAllPFCandPt      = 0.;
        double leadChargedPFCandPt = 0.;
        for ( reco::PFCandidateCollection::const_iterator pfCand = pfCands.begin();
              pfCand != pfCands.end(); ++pfCand ) {
          double dR = reco::deltaR(genHadTau->eta(), genHadTau->phi(), pfCand->eta(), pfCand->phi());
          if ( dR <= dRcone_ ) 
          {
            if ( TMath::Abs(pfCand->charge()) > 0.5 ) 
            {
              sumChargedPFCandPt += pfCand->pt();
              if ( pfCand->pt() > leadChargedPFCandPt ) leadChargedPFCandPt = pfCand->pt();
            }
            sumAllPFCandPt += pfCand->pt();
          }
        }

        double sumTrackPt  = 0.;
        double leadTrackPt = 0.;
        for ( reco::TrackCollection::const_iterator track = tracks.begin();
              track != tracks.end(); ++track ) {
          double dR = reco::deltaR(genHadTau->eta(), genHadTau->phi(), track->eta(), track->phi());
          if ( dR <= dRcone_ ) 
          {
            sumTrackPt += track->pt();
            if ( track->pt() > leadTrackPt ) leadTrackPt = track->pt();
          }
        }

        double sumChargedGenParticlePt  = 0.;
        double sumAllGenParticlePt      = 0.;
        double leadChargedGenParticlePt = 0.;
        for ( edm::View<reco::Candidate>::const_iterator genParticle = genParticles.begin();
              genParticle != genParticles.end(); ++genParticle ) {
          double dR = reco::deltaR(genHadTau->eta(), genHadTau->phi(), genParticle->eta(), genParticle->phi());
          if ( dR <= dRcone_ ) 
          {
            if ( TMath::Abs(genParticle->charge()) > 0.5 ) 
            {
              sumChargedGenParticlePt += genParticle->pt();
              if ( genParticle->pt() > leadChargedGenParticlePt ) leadChargedGenParticlePt = genParticle->pt();
            }
            sumAllGenParticlePt += genParticle->pt();
          }
        }

        if ( sumChargedGenParticlePt > 5. )
        {
          fillWithOverFlow2d(histogram_sumChargedPFCands_, genHadTau->eta(), sumChargedPFCandPt/sumChargedGenParticlePt, evtWeight);
          fillWithOverFlow2d(histogram_sumTracks_, genHadTau->eta(), sumTrackPt/sumChargedGenParticlePt, evtWeight);
        }
        if ( sumAllGenParticlePt > 5. ) 
        {
          fillWithOverFlow2d(histogram_sumAllPFCands_, genHadTau->eta(), sumAllPFCandPt/sumAllGenParticlePt, evtWeight);
        }
        if ( leadChargedGenParticlePt > 5. )
        {
          fillWithOverFlow2d(histogram_leadChargedPFCand_, genHadTau->eta(), leadChargedPFCandPt/leadChargedGenParticlePt, evtWeight);
          fillWithOverFlow2d(histogram_leadTrack_, genHadTau->eta(), leadTrackPt/leadChargedGenParticlePt, evtWeight);
        }
      } 
    }
    MonitorElement* me_sumChargedPFCands_;
    TH2* histogram_sumChargedPFCands_;
    MonitorElement* me_sumAllPFCands_;
    TH2* histogram_sumAllPFCands_;
    MonitorElement* me_leadChargedPFCand_;
    TH2* histogram_leadChargedPFCand_;
    MonitorElement* me_sumTracks_;
    TH2* histogram_sumTracks_;
    MonitorElement* me_leadTrack_;
    TH2* histogram_leadTrack_;
    std::string decayMode_;
    double dRcone_; 
  };
  std::vector<responsePlotEntryType*> responsePlots_;
};

#endif   
