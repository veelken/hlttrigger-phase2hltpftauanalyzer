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
  srcPFTaus_ = cfg.getParameter<edm::InputTag>("srcPFTaus");
  tokenPFTaus_ = consumes<reco::PFTauCollection>(srcPFTaus_);
  srcPFTauSumChargedIso_ = cfg.getParameter<edm::InputTag>("srcPFTauSumChargedIso");
  tokenPFTauSumChargedIso_ = consumes<reco::PFTauDiscriminator>(srcPFTauSumChargedIso_);
  srcRefTaus_ = cfg.getParameter<edm::InputTag>("srcRefTaus");
  if ( srcRefTaus_.label() != "" )
  {
    tokenRefTaus_ = consumes<reco::CandidateView>(srcRefTaus_);
    min_refTau_pt_ = cfg.getParameter<double>("min_refTau_pt");
    max_refTau_pt_ = cfg.getParameter<double>("max_refTau_pt");
    min_refTau_absEta_ = cfg.getParameter<double>("min_refTau_absEta");
    max_refTau_absEta_ = cfg.getParameter<double>("max_refTau_absEta");
  }

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
  isHigherPt(const std::pair<const reco::PFTau*, double>& pfTau_wChargedIso1,
	     const std::pair<const reco::PFTau*, double>& pfTau_wChargedIso2)
  {
    return pfTau_wChargedIso1.first->pt() > pfTau_wChargedIso2.first->pt();
  }

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
  edm::Handle<reco::PFTauCollection> pfTaus;
  evt.getByToken(tokenPFTaus_, pfTaus);
  edm::Handle<reco::PFTauDiscriminator> pfTauSumChargedIso;
  evt.getByToken(tokenPFTauSumChargedIso_, pfTauSumChargedIso);

  std::vector<std::pair<const reco::PFTau*, double>> pfTaus_wChargedIso_sorted;
  size_t numPFTaus = pfTaus->size();
  for ( size_t idxPFTau = 0; idxPFTau < numPFTaus; ++idxPFTau ) 
  { 
    reco::PFTauRef pfTauRef(pfTaus, idxPFTau);
    double sumChargedIso = (*pfTauSumChargedIso)[pfTauRef];
    pfTaus_wChargedIso_sorted.push_back(std::pair<const reco::PFTau*, double>(pfTauRef.get(), sumChargedIso));
  }
  std::sort(pfTaus_wChargedIso_sorted.begin(), pfTaus_wChargedIso_sorted.end(), isHigherPt);

  reco::PFTauPairCollection pfTauPairs;
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

    for ( std::vector<std::pair<const reco::PFTau*, double>>::const_iterator leadPFTau_wChargedIso = pfTaus_wChargedIso_sorted.begin();
           leadPFTau_wChargedIso != pfTaus_wChargedIso_sorted.end(); ++leadPFTau_wChargedIso ) {
      if ( !isGenMatched(*leadPFTau_wChargedIso->first, refTaus_passingAbsEtaAndPt, dRmatch_) ) continue;
      for ( std::vector<std::pair<const reco::PFTau*, double>>::const_iterator subleadPFTau_wChargedIso = leadPFTau_wChargedIso + 1;
            subleadPFTau_wChargedIso != pfTaus_wChargedIso_sorted.end(); ++subleadPFTau_wChargedIso ) {
        if ( !isGenMatched(*subleadPFTau_wChargedIso->first, refTaus_passingAbsEtaAndPt, dRmatch_) ) continue;
	pfTauPairs.push_back(reco::PFTauPair(
          leadPFTau_wChargedIso->first, leadPFTau_wChargedIso->second, 
          subleadPFTau_wChargedIso->first, subleadPFTau_wChargedIso->second));
      }
    }
  } 
  else
  {
    for ( std::vector<std::pair<const reco::PFTau*, double>>::const_iterator leadPFTau_wChargedIso = pfTaus_wChargedIso_sorted.begin();
           leadPFTau_wChargedIso != pfTaus_wChargedIso_sorted.end(); ++leadPFTau_wChargedIso ) {
      for ( std::vector<std::pair<const reco::PFTau*, double>>::const_iterator subleadPFTau_wChargedIso = leadPFTau_wChargedIso + 1;
            subleadPFTau_wChargedIso != pfTaus_wChargedIso_sorted.end(); ++subleadPFTau_wChargedIso ) {
	pfTauPairs.push_back(reco::PFTauPair(
          leadPFTau_wChargedIso->first, leadPFTau_wChargedIso->second, 
          subleadPFTau_wChargedIso->first, subleadPFTau_wChargedIso->second));
      }
    }
  }

  const double evtWeight = 1.;

  for ( auto efficiency_or_ratePlot : efficiency_or_ratePlots_ ) 
  {
    efficiency_or_ratePlot->fillHistograms(pfTauPairs, evtWeight);
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
