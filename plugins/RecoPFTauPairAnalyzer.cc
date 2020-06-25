#include "HLTTrigger/TallinnHLTPFTauAnalyzer/plugins/RecoPFTauPairAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"      // edm::Handle
#include "DataFormats/Math/interface/deltaR.h"        // reco::deltaR

#include "TMath.h"   // TMath::Abs()

#include <iostream>
#include <iomanip>
#include <algorithm> // std::sort

RecoPFTauPairAnalyzer::RecoPFTauPairAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
  , min_refTau_pt_(-1.)
  , max_refTau_pt_(-1.)
  , min_refTau_absEta_(-1.)
  , max_refTau_absEta_(-1.)
  , dRmatch_(0.3)
{
  srcPFTauPairs_ = cfg.getParameter<edm::InputTag>("srcPFTauPairs");
  tokenPFTauPairs_ = consumes<reco::PFTauPairCollection>(srcPFTauPairs_);
  srcRefTaus_ = cfg.getParameter<edm::InputTag>("srcRefTaus");
  if ( srcRefTaus_.label() != "" )
  {
    tokenRefTaus_ = consumes<reco::CandidateView>(srcRefTaus_);
    min_refTau_pt_ = cfg.getParameter<double>("min_refTau_pt");
    max_refTau_pt_ = cfg.getParameter<double>("max_refTau_pt");
    min_refTau_absEta_ = cfg.getParameter<double>("min_refTau_absEta");
    max_refTau_absEta_ = cfg.getParameter<double>("max_refTau_absEta");
  }

  lumiScale_ = cfg.getParameter<double>("lumiScale"); // for background: rate in Hz corresponding to one MC event

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

RecoPFTauPairAnalyzer::~RecoPFTauPairAnalyzer()
{
  for ( auto efficiency_or_ratePlot : efficiency_or_ratePlots_ ) 
  {
    delete efficiency_or_ratePlot;
  }
}

void RecoPFTauPairAnalyzer::beginJob()
{
  if ( !edm::Service<dqm::legacy::DQMStore>().isAvailable() ) {
    throw cms::Exception("RecoPFTauPairAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<dqm::legacy::DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory_.data());
  
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(-1.,  1.4,   0.40, -1., 0.4)); // vLoose
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(-1.,  1.4,   0.20, -1., 0.4)); // Loose
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(-1.,  1.4,   0.10, -1., 0.4)); // Medium
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(-1.,  1.4,   0.05, -1., 0.4)); // Tight
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(-1.,  1.4,   0.02, -1., 0.4)); // vTight
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(-1.,  1.4,   0.01, -1., 0.4)); // vvTight

  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType( 1.4, 2.172, 0.40, -1., 0.4)); // vLoose
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType( 1.4, 2.172, 0.20, -1., 0.4)); // Loose
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType( 1.4, 2.172, 0.10, -1., 0.4)); // Medium
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType( 1.4, 2.172, 0.05, -1., 0.4)); // Tight
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType( 1.4, 2.172, 0.02, -1., 0.4)); // vTight
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType( 1.4, 2.172, 0.01, -1., 0.4)); // vvTight

  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType( 1.4, 2.4,   0.40, -1., 0.4)); // vLoose
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType( 1.4, 2.4,   0.20, -1., 0.4)); // Loose
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType( 1.4, 2.4,   0.10, -1., 0.4)); // Medium
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType( 1.4, 2.4,   0.05, -1., 0.4)); // Tight
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType( 1.4, 2.4,   0.02, -1., 0.4)); // vTight
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType( 1.4, 2.4,   0.01, -1., 0.4)); // vvTight

  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(-1.,  2.172, 0.40, -1., 0.4)); // vLoose
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(-1.,  2.172, 0.20, -1., 0.4)); // Loose
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(-1.,  2.172, 0.10, -1., 0.4)); // Medium
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(-1.,  2.172, 0.05, -1., 0.4)); // Tight
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(-1.,  2.172, 0.02, -1., 0.4)); // vTight
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(-1.,  2.172, 0.01, -1., 0.4)); // vvTight

  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(-1.,  2.4,   0.40, -1., 0.4)); // vLoose
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(-1.,  2.4,   0.20, -1., 0.4)); // Loose
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(-1.,  2.4,   0.10, -1., 0.4)); // Medium
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(-1.,  2.4,   0.05, -1., 0.4)); // Tight
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(-1.,  2.4,   0.02, -1., 0.4)); // vTight
  efficiency_or_ratePlots_.push_back(new efficiency_or_ratePlotEntryType(-1.,  2.4,   0.01, -1., 0.4)); // vvTight

  for ( auto efficiency_or_ratePlot : efficiency_or_ratePlots_ ) 
  {
    efficiency_or_ratePlot->bookHistograms(dqmStore);
  }
}

namespace
{
  bool
  isGenMatched(const reco::PFTau& pfTau, const std::vector<const reco::Candidate*>& refTaus, double dRmatch)
  {
    for ( std::vector<const reco::Candidate*>::const_iterator refTau = refTaus.begin();
	  refTau != refTaus.end(); ++refTau ) {
      double dR = reco::deltaR(pfTau.eta(), pfTau.phi(), (*refTau)->eta(), (*refTau)->phi());
      if ( dR < dRmatch ) return true;
    }
    return false;
  }
}

void RecoPFTauPairAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<reco::PFTauPairCollection> pfTauPairs;
  evt.getByToken(tokenPFTauPairs_, pfTauPairs);

  std::vector<const reco::PFTauPair*> pfTauPairs_matched;
  if ( srcRefTaus_.label() != "" ) 
  {
    edm::Handle<reco::CandidateView> refTaus;
    evt.getByToken(tokenRefTaus_, refTaus);

    std::vector<const reco::Candidate*> refTaus_passingAbsEtaAndPt;
    for ( reco::CandidateView::const_iterator refTau = refTaus->begin();
	  refTau != refTaus->end(); ++refTau ) {
      double refTau_absEta = TMath::Abs(refTau->eta());
      if ( refTau->pt() > min_refTau_pt_ && refTau->pt() < max_refTau_pt_ && refTau_absEta > min_refTau_absEta_ && refTau_absEta < max_refTau_absEta_ )
      {
	refTaus_passingAbsEtaAndPt.push_back(&(*refTau));
      }
    }
    if ( !(refTaus_passingAbsEtaAndPt.size() >= 2) ) return;

    for ( reco::PFTauPairCollection::const_iterator pfTauPair = pfTauPairs->begin();
          pfTauPair != pfTauPairs->end(); ++pfTauPair ) {
      if ( isGenMatched(*pfTauPair->leadPFTau(),    refTaus_passingAbsEtaAndPt, dRmatch_) &&
           isGenMatched(*pfTauPair->subleadPFTau(), refTaus_passingAbsEtaAndPt, dRmatch_) ) {
        pfTauPairs_matched.push_back(&(*pfTauPair));
      }
    }
  } 
  else 
  {
    for ( reco::PFTauPairCollection::const_iterator pfTauPair = pfTauPairs->begin();
          pfTauPair != pfTauPairs->end(); ++pfTauPair ) {
      pfTauPairs_matched.push_back(&(*pfTauPair));
    }
  }

  const double evtWeight = lumiScale_;

  for ( auto efficiency_or_ratePlot : efficiency_or_ratePlots_ ) 
  {
    efficiency_or_ratePlot->fillHistograms(pfTauPairs_matched, evtWeight);
  }
}

void RecoPFTauPairAnalyzer::endJob()
{
  for ( auto efficiency_or_ratePlot : efficiency_or_ratePlots_ ) 
  {
    efficiency_or_ratePlot->normalizeHistograms();
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(RecoPFTauPairAnalyzer);
