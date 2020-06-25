#include "HLTTrigger/TallinnHLTPFTauAnalyzer/plugins/RecoPFTauAnalyzerSignal.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

RecoPFTauAnalyzerSignal::RecoPFTauAnalyzerSignal(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  srcPFTaus_ = cfg.getParameter<edm::InputTag>("srcPFTaus");
  tokenPFTaus_ = consumes<reco::PFTauCollection>(srcPFTaus_);
  srcPFTauSumChargedIso_ = cfg.getParameter<edm::InputTag>("srcPFTauSumChargedIso");
  tokenPFTauSumChargedIso_ = consumes<reco::PFTauDiscriminator>(srcPFTauSumChargedIso_);
  srcDenominator_ = cfg.getParameter<edm::InputTag>("srcDenominator");
  std::string typeDenominator_string = cfg.getParameter<std::string>("typeDenominator");
  if ( typeDenominator_string == "gen" ) 
  {
    typeDenominator_ = kGen;
    tokenDenominator_gen_ = consumes<reco::GenJetCollection>(srcDenominator_);
  }
  else if ( typeDenominator_string == "offline" ) 
  {
    typeDenominator_ = kOffline;
    tokenDenominator_offline_ = consumes<pat::TauCollection>(srcDenominator_);
  }
  else
  {
    throw cms::Exception("TallinnL1PFTauAnalyzerSignal") 
      << " Invalid Configuration parameter 'typeDenominator' = " << typeDenominator_string << " !!\n";;
  }

  lumiScale_ = cfg.getParameter<double>("lumiScale");

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

RecoPFTauAnalyzerSignal::~RecoPFTauAnalyzerSignal()
{
  for ( auto efficiencyPlot : efficiencyPlots_ ) 
  {
    delete efficiencyPlot;
  }
}

void RecoPFTauAnalyzerSignal::beginJob()
{
  if ( !edm::Service<dqm::legacy::DQMStore>().isAvailable() ) {
    throw cms::Exception("RecoPFTauAnalyzerSignal") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<dqm::legacy::DQMStore>());

  std::vector<std::string> decayModes = { "oneProng0Pi0", "oneProng1Pi0", "oneProng2Pi0", "threeProng0Pi0", "threeProng1Pi0", "all" };
  std::vector<double> min_absEtaValues = { -1.,   1.4,   1.4, -1.,    -1.  };
  std::vector<double> max_absEtaValues = {  1.4,  2.172, 2.4,  2.172,  2.4 };
  assert(min_absEtaValues.size() == max_absEtaValues.size());
  size_t numAbsEtaRanges = min_absEtaValues.size();
  std::vector<double> ptThresholds = { 20., 25., 30., 35., 40., 45., 50. };
  std::vector<double> min_leadTrackPtValues = { 1., 2., 5. };
  for ( size_t idxAbsEtaRange = 0; idxAbsEtaRange < numAbsEtaRanges; ++idxAbsEtaRange )
  {
    double min_absEta = min_absEtaValues[idxAbsEtaRange];
    double max_absEta = max_absEtaValues[idxAbsEtaRange];
    for ( auto decayMode : decayModes )
    {
      for ( auto ptThreshold : ptThresholds )
      {
        for ( auto min_leadTrackPt : min_leadTrackPtValues )
        {
  	  TString dqmDirectory = dqmDirectory_.data();
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
          efficiencyPlotEntryType* efficiencyPlots_noIsolation  = new efficiencyPlotEntryType(45., 1.e+3, min_absEta, max_absEta, decayMode, 
            ptThreshold, min_leadTrackPt, -1.,   -1.); // no isolation cut applied
	  efficiencyPlots_noIsolation->bookHistograms(dqmStore);
          efficiencyPlotEntryType* efficiencyPlots_vLoose  = new efficiencyPlotEntryType(45., 1.e+3, min_absEta, max_absEta, decayMode, 
            ptThreshold, min_leadTrackPt,  0.40, -1.); // vLoose
	  efficiencyPlots_vLoose->bookHistograms(dqmStore);
	  efficiencyPlots_.push_back(efficiencyPlots_vLoose);
	  efficiencyPlotEntryType* efficiencyPlots_Loose   = new efficiencyPlotEntryType(45., 1.e+3, min_absEta, max_absEta, decayMode, 
            ptThreshold, min_leadTrackPt,  0.20, -1.); // Loose
	  efficiencyPlots_Loose->bookHistograms(dqmStore);
	  efficiencyPlots_.push_back(efficiencyPlots_Loose);
	  efficiencyPlotEntryType* efficiencyPlots_Medium  = new efficiencyPlotEntryType(45., 1.e+3, min_absEta, max_absEta, decayMode, 
            ptThreshold, min_leadTrackPt,  0.10, -1.); // Medium
	  efficiencyPlots_Medium->bookHistograms(dqmStore);
	  efficiencyPlots_.push_back(efficiencyPlots_Medium);
	  efficiencyPlotEntryType* efficiencyPlots_Tight   = new efficiencyPlotEntryType(45., 1.e+3, min_absEta, max_absEta, decayMode, 
            ptThreshold, min_leadTrackPt,  0.05, -1.); // Tight
	  efficiencyPlots_Tight->bookHistograms(dqmStore);
	  efficiencyPlots_.push_back(efficiencyPlots_Tight);
	  efficiencyPlotEntryType* efficiencyPlots_vTight  = new efficiencyPlotEntryType(45., 1.e+3, min_absEta, max_absEta, decayMode, 
            ptThreshold, min_leadTrackPt,  0.02, -1.); // vTight
	  efficiencyPlots_vTight->bookHistograms(dqmStore);
	  efficiencyPlots_.push_back(efficiencyPlots_vTight);
	  efficiencyPlotEntryType* efficiencyPlots_vvTight = new efficiencyPlotEntryType(45., 1.e+3, min_absEta, max_absEta, decayMode,
            ptThreshold, min_leadTrackPt,  0.01, -1.); // vvTight
	  efficiencyPlots_vvTight->bookHistograms(dqmStore);
	  efficiencyPlots_.push_back(efficiencyPlots_vvTight);
        }
      }
    }
  }
}

void RecoPFTauAnalyzerSignal::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<reco::PFTauCollection> pfTaus;
  evt.getByToken(tokenPFTaus_, pfTaus);
  
  edm::Handle<reco::PFTauDiscriminator> pfTauSumChargedIso;
  evt.getByToken(tokenPFTauSumChargedIso_, pfTauSumChargedIso);
  
  const double evtWeight = lumiScale_;

  if ( typeDenominator_ == kGen )
  {
    edm::Handle<reco::GenJetCollection> denominatorTaus_gen;
    evt.getByToken(tokenDenominator_gen_, denominatorTaus_gen);

    for ( auto efficiencyPlot : efficiencyPlots_ ) 
    {    
      efficiencyPlot->fillHistograms(pfTaus, *pfTauSumChargedIso, *denominatorTaus_gen, evtWeight);
    }
  }

  if ( typeDenominator_ == kOffline )
  {
    edm::Handle<pat::TauCollection> denominatorTaus_offline;
    evt.getByToken(tokenDenominator_offline_, denominatorTaus_offline);

    for ( auto efficiencyPlot : efficiencyPlots_ ) 
    {    
      efficiencyPlot->fillHistograms(pfTaus, *pfTauSumChargedIso, *denominatorTaus_offline, evtWeight);
    }
  }
}

void RecoPFTauAnalyzerSignal::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(RecoPFTauAnalyzerSignal);
