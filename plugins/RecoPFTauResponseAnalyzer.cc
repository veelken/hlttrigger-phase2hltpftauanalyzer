#include "HLTrigger/TallinnHLTPFTauAnalyzer/plugins/RecoPFTauResponseAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

RecoPFTauResponseAnalyzer::RecoPFTauResponseAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  srcPFTaus_ = cfg.getParameter<edm::InputTag>("srcPFTaus");
  tokenPFTaus_ = consumes<reco::PFTauCollection>(srcPFTaus_);
  srcPFTauSumChargedIso_ = cfg.getParameter<edm::InputTag>("srcPFTauSumChargedIso");
  tokenPFTauSumChargedIso_ = consumes<reco::PFTauDiscriminator>(srcPFTauSumChargedIso_);
  srcRefTaus_ = cfg.getParameter<edm::InputTag>("srcRefTaus");
  std::string typeRefTaus_string = cfg.getParameter<std::string>("typeRefTaus");
  if ( typeRefTaus_string == "gen"     ) 
  {
    typeRefTaus_ = kGen;
    tokenRefTaus_gen_ = consumes<reco::GenJetCollection>(srcRefTaus_);
  }
  else if ( typeRefTaus_string == "offline" ) 
  {
    typeRefTaus_ = kOffline;
    tokenRefTaus_offline_ = consumes<pat::TauCollection>(srcRefTaus_);
  }
  else
  {
    throw cms::Exception("TallinnL1PFTauResponseAnalyzer") 
      << " Invalid Configuration parameter 'typeRefTaus' = " << typeRefTaus_string << " !!\n";;
  }

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

RecoPFTauResponseAnalyzer::~RecoPFTauResponseAnalyzer()
{
  for ( auto responsePlot : responsePlots_ ) 
  {
    delete responsePlot;
  }
}

void RecoPFTauResponseAnalyzer::beginJob()
{
  if ( !edm::Service<dqm::legacy::DQMStore>().isAvailable() ) {
    throw cms::Exception("RecoPFTauResponseAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<dqm::legacy::DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory_.data());

  std::vector<std::string> decayModes = { "oneProng0Pi0", "oneProng1Pi0", "oneProng2Pi0", "threeProng0Pi0", "threeProng1Pi0", "all" };
  std::vector<double> min_ptValues = { 20., 25., 30., 35., 40., 45., 50. };
  std::vector<double> max_ptValues = { -1., -1., -1., -1., -1., -1., -1. }; 
  assert(min_ptValues.size() == max_ptValues.size());
  size_t numPtRanges = min_ptValues.size();
  std::vector<double> min_absEtaValues = { -1.,   1.4,   1.4, -1.,    -1.  };
  std::vector<double> max_absEtaValues = {  1.4,  2.172, 2.4,  2.172,  2.4 };
  assert(min_absEtaValues.size() == max_absEtaValues.size());
  size_t numAbsEtaRanges = min_absEtaValues.size();
  std::vector<double> min_leadTrackPtValues = { 1., 2., 5. };
  for ( size_t idxPtRange = 0; idxPtRange < numPtRanges; ++idxPtRange )
  {
    double min_pt = min_ptValues[idxPtRange];
    double max_pt = max_ptValues[idxPtRange];
    for ( size_t idxAbsEtaRange = 0; idxAbsEtaRange < numAbsEtaRanges; ++idxAbsEtaRange )
    {
      double min_absEta = min_absEtaValues[idxAbsEtaRange];
      double max_absEta = max_absEtaValues[idxAbsEtaRange];
      for ( auto decayMode : decayModes )
      {
        for ( auto min_leadTrackPt : min_leadTrackPtValues )
        {
	  TString dqmDirectory = dqmDirectory_.data();
          if ( min_pt >= 0. && max_pt > 0. ) 
	  { 
	    dqmDirectory.Append(Form("/pt%1.0fto%1.0f", min_pt, max_pt));
          }
          else if ( min_pt >= 0. ) 
          {
            dqmDirectory.Append(Form("/ptGt%1.0f", min_pt));
          }
          else if ( max_pt > 0. ) 
          {
            dqmDirectory.Append(Form("/ptLt%1.0f", max_pt));
          }
	  if ( min_absEta >= 0. && max_absEta > 0. ) 
  	  { 
	   dqmDirectory.Append(Form("/absEta%1.2fto%1.2f", min_absEta, max_absEta));
          }
          else if ( min_absEta >= 0. ) 
          {
            dqmDirectory.Append(Form("/absEtaGt%1.2f", min_absEta));
          }
          else if ( max_absEta > 0. ) 
          {
            dqmDirectory.Append(Form("/absEtaLt%1.2f", max_absEta));
          }
          //std::string decayMode_capitalized = decayMode;
          //decayMode_capitalized[0] = toupper(decayMode_capitalized[0]);	
          //dqmDirectory.Append(Form("/gen%sTau", decayMode_capitalized.data()));
          dqmDirectory.Append(Form("/%s", decayMode.data()));
          dqmDirectory = dqmDirectory.ReplaceAll(".", "p");

          dqmStore.setCurrentFolder(dqmDirectory.Data());
          responsePlotEntryType* responsePlots_vLoose  = new responsePlotEntryType(min_pt, max_pt, min_absEta, max_absEta, decayMode, min_leadTrackPt, 0.40, -1.); // vLoose
  	  responsePlots_vLoose->bookHistograms(dqmStore);
          responsePlots_.push_back(responsePlots_vLoose);
          responsePlotEntryType* responsePlots_Loose   = new responsePlotEntryType(min_pt, max_pt, min_absEta, max_absEta, decayMode, min_leadTrackPt, 0.20, -1.); // Loose
	  responsePlots_Loose->bookHistograms(dqmStore);
          responsePlots_.push_back(responsePlots_Loose);
          responsePlotEntryType* responsePlots_Medium  = new responsePlotEntryType(min_pt, max_pt, min_absEta, max_absEta, decayMode, min_leadTrackPt, 0.10, -1.); // Medium
  	  responsePlots_Medium->bookHistograms(dqmStore);
          responsePlots_.push_back(responsePlots_Medium);
          responsePlotEntryType* responsePlots_Tight   = new responsePlotEntryType(min_pt, max_pt, min_absEta, max_absEta, decayMode, min_leadTrackPt, 0.05, -1.); // Tight
	  responsePlots_Tight->bookHistograms(dqmStore);
          responsePlots_.push_back(responsePlots_Tight);
          responsePlotEntryType* responsePlots_vTight  = new responsePlotEntryType(min_pt, max_pt, min_absEta, max_absEta, decayMode, min_leadTrackPt, 0.02, -1.); // vTight
	  responsePlots_vTight->bookHistograms(dqmStore);
          responsePlots_.push_back(responsePlots_vTight);
          responsePlotEntryType* responsePlots_vvTight = new responsePlotEntryType(min_pt, max_pt, min_absEta, max_absEta, decayMode, min_leadTrackPt, 0.01, -1.); // vvTight
	  responsePlots_vvTight->bookHistograms(dqmStore);
          responsePlots_.push_back(responsePlots_vvTight);
        }
      }
    }
  }
}

void RecoPFTauResponseAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<reco::PFTauCollection> pfTaus;
  evt.getByToken(tokenPFTaus_, pfTaus);
  
  edm::Handle<reco::PFTauDiscriminator> pfTauSumChargedIso;
  evt.getByToken(tokenPFTauSumChargedIso_, pfTauSumChargedIso);

  const double evtWeight = 1.;

  if ( typeRefTaus_ == kGen )
  {
    edm::Handle<reco::GenJetCollection> refTaus_gen;
    evt.getByToken(tokenRefTaus_gen_, refTaus_gen);

    for ( auto responsePlot : responsePlots_ ) 
    {    
      responsePlot->fillHistograms(pfTaus, *pfTauSumChargedIso, *refTaus_gen, evtWeight);
    }
  }

  if ( typeRefTaus_ == kOffline )
  {
    edm::Handle<pat::TauCollection> refTaus_offline;
    evt.getByToken(tokenRefTaus_offline_, refTaus_offline);

    for ( auto responsePlot : responsePlots_ ) 
    {    
      responsePlot->fillHistograms(pfTaus, *pfTauSumChargedIso, *refTaus_offline, evtWeight);
    }
  }
}

void RecoPFTauResponseAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(RecoPFTauResponseAnalyzer);
