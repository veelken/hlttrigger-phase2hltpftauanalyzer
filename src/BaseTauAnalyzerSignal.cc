#include "HLTrigger/TallinnHLTPFTauAnalyzer/interface/BaseTauAnalyzerSignal.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Common/interface/Handle.h"

#include "HLTrigger/TallinnHLTPFTauAnalyzer/interface/generalAuxFunctions.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

typedef std::vector<double> vdouble;
typedef std::vector<std::string> vstring;

namespace
{
  void
  checkArray(const std::string& arrayName, const vdouble& array_min_values, const vdouble& array_max_values)
  {
    if ( array_min_values.size() != array_max_values.size() )
      throw cms::Exception("BaseTauAnalyzerSignal") 
        << " Mismatch in number of 'min_" << arrayName << "' and 'max_" << arrayName << "' Configuration parameters !!\n";
  }
}

BaseTauAnalyzerSignal::BaseTauAnalyzerSignal(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
{
  //std::cout << "<BaseTauAnalyzerSignal::BaseTauAnalyzerSignal (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;

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

  min_pt_denominator_ = cfg.getParameter<double>("min_pt_denominator");
  max_pt_denominator_ = cfg.getParameter<double>("max_pt_denominator");
  //std::cout << "min_pt_denominator = " << min_pt_denominator_ << std::endl;
  //std::cout << "max_pt_denominator = " << max_pt_denominator_ << std::endl;

  min_ptValues_numerator_ = cfg.getParameter<vdouble>("min_pt_numerator");
  max_ptValues_numerator_ = cfg.getParameter<vdouble>("max_pt_numerator");
  checkArray("pt_numerator", min_ptValues_numerator_, max_ptValues_numerator_);
  num_ptValues_numerator_ = min_ptValues_numerator_.size();
  //std::cout << "num_ptValues_numerator = " << num_ptValues_numerator_ << std::endl;
  //std::cout << " min_ptValues_numerator = " << format_vdouble(min_ptValues_numerator_) << std::endl;
  //std::cout << " max_ptValues_numerator = " << format_vdouble(max_ptValues_numerator_) << std::endl;

  min_absEtaValues_ = cfg.getParameter<vdouble>("min_absEta");
  max_absEtaValues_ = cfg.getParameter<vdouble>("max_absEta");
  checkArray("absEta", min_absEtaValues_, max_absEtaValues_);
  num_absEtaValues_ = min_absEtaValues_.size();
  //std::cout << "num_absEtaValues = " << num_absEtaValues_ << std::endl;
  //std::cout << " min_absEtaValues = " << format_vdouble(min_absEtaValues_) << std::endl;
  //std::cout << " max_absEtaValues = " << format_vdouble(max_absEtaValues_) << std::endl;

  decayModes_ = cfg.getParameter<vstring>("decayModes");
  //std::cout << "decayModes = " << format_vstring(decayModes_) << std::endl;

  min_leadTrackPtValues_ = cfg.getParameter<vdouble>("min_leadTrackPt");
  max_leadTrackPtValues_ = cfg.getParameter<vdouble>("max_leadTrackPt");
  checkArray("leadTrackPt", min_leadTrackPtValues_, max_leadTrackPtValues_);
  num_leadTrackPtValues_ = min_leadTrackPtValues_.size();
  //std::cout << "num_leadTrackPtValues = " << num_leadTrackPtValues_ << std::endl;
  //std::cout << " min_leadTrackPtValues = " << format_vdouble(min_leadTrackPtValues_) << std::endl;
  //std::cout << " max_leadTrackPtValues = " << format_vdouble(max_leadTrackPtValues_) << std::endl;

  if ( cfg.exists("min_relDiscriminator") && cfg.exists("max_relDiscriminator") ) // used to apply cut on charged-isolation pT sum relative to tau pT
  {
    min_relDiscriminatorValues_ = cfg.getParameter<vdouble>("min_relDiscriminator");
    max_relDiscriminatorValues_ = cfg.getParameter<vdouble>("max_relDiscriminator");
    checkArray("relDiscriminator", min_relDiscriminatorValues_, max_relDiscriminatorValues_);
    num_relDiscriminatorValues_ = min_relDiscriminatorValues_.size();
    //std::cout << "num_relDiscriminatorValues = " << num_relDiscriminatorValues_ << std::endl;
    //std::cout << " min_relDiscriminatorValues = " << format_vdouble(min_relDiscriminatorValues_) << std::endl;
    //std::cout << " max_relDiscriminatorValues = " << format_vdouble(max_relDiscriminatorValues_) << std::endl;
  }
  if ( cfg.exists("min_relDiscriminator") && cfg.exists("max_relDiscriminator") ) // used to apply cut on DeepTau tau ID discriminant, independent of tau pT
  {
    min_absDiscriminatorValues_ = cfg.getParameter<vdouble>("min_absDiscriminator");
    max_absDiscriminatorValues_ = cfg.getParameter<vdouble>("max_absDiscriminator");
    checkArray("absDiscriminator", min_absDiscriminatorValues_, max_absDiscriminatorValues_);
    num_absDiscriminatorValues_ = min_absDiscriminatorValues_.size();
    //std::cout << "num_absDiscriminatorValues = " << num_absDiscriminatorValues_ << std::endl;
    //std::cout << " min_absDiscriminatorValues = " << format_vdouble(min_absDiscriminatorValues_) << std::endl;
    //std::cout << " max_absDiscriminatorValues = " << format_vdouble(max_absDiscriminatorValues_) << std::endl;
  }
  if ( (max_relDiscriminatorValues_.size() + max_absDiscriminatorValues_.size()) == 0 )
    throw cms::Exception("BaseTauAnalyzerSignal") 
      << " No '' or '' Configuration parameters defined !!\n";

  src_evtWeight_ = cfg.getParameter<edm::InputTag>("src_evtWeight");
  token_evtWeight_ = consumes<double>(src_evtWeight_);
  //std::cout << "src_evtWeight = " << src_evtWeight_.label() << std::endl;

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

BaseTauAnalyzerSignal::~BaseTauAnalyzerSignal()
{
  for ( auto efficiencyPlot : efficiencyPlots_ ) 
  {
    delete efficiencyPlot;
  }
}

void BaseTauAnalyzerSignal::beginJob()
{
  //std::cout << "<BaseTauAnalyzerSignal::beginJob>:" << std::endl;

  if ( !edm::Service<dqm::legacy::DQMStore>().isAvailable() ) {
    throw cms::Exception("BaseTauAnalyzerSignal") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<dqm::legacy::DQMStore>());

  for ( size_t idx_pt_numerator = 0; idx_pt_numerator < num_ptValues_numerator_; ++idx_pt_numerator )
  {
    double min_pt_numerator = min_ptValues_numerator_[idx_pt_numerator];
    double max_pt_numerator = max_ptValues_numerator_[idx_pt_numerator];
    for ( size_t idx_absEta = 0; idx_absEta < num_absEtaValues_; ++idx_absEta )
    {
      double min_absEta = min_absEtaValues_[idx_absEta];
      double max_absEta = max_absEtaValues_[idx_absEta];
      for ( auto decayMode : decayModes_ )
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
        //std::cout << "dqmDirectory = " << dqmDirectory << std::endl;      
	dqmStore.setCurrentFolder(dqmDirectory.Data());

        for ( size_t idx_leadTrackPt = 0; idx_leadTrackPt < num_leadTrackPtValues_; ++idx_leadTrackPt )
        {
          double min_leadTrackPt = min_leadTrackPtValues_[idx_leadTrackPt];
          double max_leadTrackPt = max_leadTrackPtValues_[idx_leadTrackPt];
          for ( size_t idx_relDiscriminator = 0; idx_relDiscriminator < num_relDiscriminatorValues_; ++idx_relDiscriminator )
          {
            double min_relDiscriminator = min_relDiscriminatorValues_[idx_relDiscriminator];
            double max_relDiscriminator = max_relDiscriminatorValues_[idx_relDiscriminator];
            efficiencyPlotEntryType* efficiencyPlot = new efficiencyPlotEntryType(
              min_pt_numerator, max_pt_numerator, min_pt_denominator_, max_pt_denominator_, 
              min_absEta, max_absEta,
              decayMode,
              min_leadTrackPt, max_leadTrackPt, 
              min_relDiscriminator, max_relDiscriminator, -1., -1.);
            efficiencyPlot->bookHistograms(dqmStore);
            efficiencyPlots_.push_back(efficiencyPlot);
          } // idx_relDiscriminator
          for ( size_t idx_absDiscriminator = 0; idx_absDiscriminator < num_absDiscriminatorValues_; ++idx_absDiscriminator )
          {
            double min_absDiscriminator = min_absDiscriminatorValues_[idx_absDiscriminator];
            double max_absDiscriminator = max_absDiscriminatorValues_[idx_absDiscriminator];
            efficiencyPlotEntryType* efficiencyPlot = new efficiencyPlotEntryType(
              min_pt_numerator, max_pt_numerator, min_pt_denominator_, max_pt_denominator_, 
              min_absEta, max_absEta,
              decayMode,
              min_leadTrackPt, max_leadTrackPt, 
              -1., -1., min_absDiscriminator, max_absDiscriminator);
            efficiencyPlot->bookHistograms(dqmStore);
            efficiencyPlots_.push_back(efficiencyPlot);
          } // idx_absDiscriminator
        } // idx_leadTrackPt
      } // decayMode
    } // idx_absEta
  } // idx_pt_numerator
}

namespace
{
  bool
  isHigherPt(const BaseTau& tau1, const BaseTau& tau2)
  {
    return tau1.p4().pt() > tau2.p4().pt();
  }
}

void BaseTauAnalyzerSignal::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  std::vector<BaseTau> baseTaus = buildBaseTaus(evt, es);
  
  // CV: sort taus by decreasing pT
  std::sort(baseTaus.begin(), baseTaus.end(), isHigherPt);
  
  edm::Handle<double> evtWeight;
  evt.getByToken(token_evtWeight_, evtWeight);

  if ( typeDenominator_ == kGen )
  {
    edm::Handle<reco::GenJetCollection> denominatorTaus_gen;
    evt.getByToken(tokenDenominator_gen_, denominatorTaus_gen);

    for ( auto efficiencyPlot : efficiencyPlots_ ) 
    {    
      efficiencyPlot->fillHistograms(baseTaus, *denominatorTaus_gen, *evtWeight);
    }
  }

  if ( typeDenominator_ == kOffline )
  {
    edm::Handle<pat::TauCollection> denominatorTaus_offline;
    evt.getByToken(tokenDenominator_offline_, denominatorTaus_offline);

    for ( auto efficiencyPlot : efficiencyPlots_ ) 
    {    
      efficiencyPlot->fillHistograms(baseTaus, *denominatorTaus_offline, *evtWeight);
    }
  }
}

void BaseTauAnalyzerSignal::endJob()
{}
