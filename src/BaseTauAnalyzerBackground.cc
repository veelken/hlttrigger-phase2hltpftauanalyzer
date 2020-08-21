#include "HLTrigger/TallinnHLTPFTauAnalyzer/interface/BaseTauAnalyzerBackground.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

typedef std::vector<double> vdouble;

namespace
{
  void
  checkArray(const std::string& arrayName, const vdouble& array_min_values, const vdouble& array_max_values)
  {
    if ( array_min_values.size() != array_max_values.size() )
      throw cms::Exception("BaseTauAnalyzerBackground") 
        << " Mismatch in number of 'min_" << arrayName << "' and 'max_" << arrayName << "' Configuration parameters !!\n";
  }
}

BaseTauAnalyzerBackground::BaseTauAnalyzerBackground(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
  , num_absEtaValues_(0)
  , num_leadTrackPtValues_(0)
  , num_relDiscriminatorValues_(0)
  , num_absDiscriminatorValues_(0)
  , num_dzValues_(0)
{
  min_pt_ = cfg.getParameter<double>("min_pt");
  max_pt_ = cfg.getParameter<double>("max_pt");

  min_absEtaValues_ = cfg.getParameter<vdouble>("min_absEta");
  max_absEtaValues_ = cfg.getParameter<vdouble>("max_absEta");
  checkArray("absEta", min_absEtaValues_, max_absEtaValues_);
  num_absEtaValues_ = min_absEtaValues_.size();

  min_leadTrackPtValues_ = cfg.getParameter<vdouble>("min_leadTrackPt");
  max_leadTrackPtValues_ = cfg.getParameter<vdouble>("max_leadTrackPt");
  checkArray("leadTrackPt", min_leadTrackPtValues_, max_leadTrackPtValues_);
  num_leadTrackPtValues_ = min_leadTrackPtValues_.size();

  if ( cfg.exists("min_relDiscriminator") && cfg.exists("max_relDiscriminator") ) // used to apply cut on charged-isolation pT sum relative to tau pT
  {
    min_relDiscriminatorValues_ = cfg.getParameter<vdouble>("min_relDiscriminator");
    max_relDiscriminatorValues_ = cfg.getParameter<vdouble>("max_relDiscriminator");
    checkArray("relDiscriminator", min_relDiscriminatorValues_, max_relDiscriminatorValues_);
    num_relDiscriminatorValues_ = min_relDiscriminatorValues_.size();
  }
  if ( cfg.exists("min_relDiscriminator") && cfg.exists("max_relDiscriminator") ) // used to apply cut on DeepTau tau ID discriminant, independent of tau pT
  {
    min_absDiscriminatorValues_ = cfg.getParameter<vdouble>("min_absDiscriminator");
    max_absDiscriminatorValues_ = cfg.getParameter<vdouble>("max_absDiscriminator");
    checkArray("absDiscriminator", min_absDiscriminatorValues_, max_absDiscriminatorValues_);
    num_absDiscriminatorValues_ = min_absDiscriminatorValues_.size();
  }
  if ( (max_relDiscriminatorValues_.size() + max_absDiscriminatorValues_.size()) == 0 )
    throw cms::Exception("BaseTauAnalyzerBackground") 
      << " No '' or '' Configuration parameters defined !!\n";

  min_dzValues_ = cfg.getParameter<vdouble>("min_dzValues");
  max_dzValues_ = cfg.getParameter<vdouble>("max_dzValues"); 
  checkArray("dz", min_dzValues_, max_dzValues_);
  num_dzValues_ = min_dzValues_.size();

  src_evtWeight_ = cfg.getParameter<edm::InputTag>("src_evtWeight");
  token_evtWeight_ = consumes<double>(src_evtWeight_);

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

BaseTauAnalyzerBackground::~BaseTauAnalyzerBackground()
{
  for ( auto ratePlot : ratePlots_ ) 
  {
    delete ratePlot;
  }
}

void BaseTauAnalyzerBackground::beginJob()
{
  if ( !edm::Service<dqm::legacy::DQMStore>().isAvailable() ) {
    throw cms::Exception("BaseTauAnalyzerBackground") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<dqm::legacy::DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory_.data());

  for ( size_t idx_absEta = 0; idx_absEta < num_absEtaValues_; ++idx_absEta )
  {
    double min_absEta = min_absEtaValues_[idx_absEta];
    double max_absEta = max_absEtaValues_[idx_absEta];
    for ( size_t idx_leadTrackPt = 0; idx_leadTrackPt < num_leadTrackPtValues_; ++idx_leadTrackPt )
    {
      double min_leadTrackPt = min_leadTrackPtValues_[idx_leadTrackPt];
      double max_leadTrackPt = max_leadTrackPtValues_[idx_leadTrackPt];
      for ( size_t idx_dz = 0; idx_dz < num_dzValues_; ++idx_dz )
      {
        double min_dz = min_dzValues_[idx_dz];
        double max_dz = max_dzValues_[idx_dz];
        for ( size_t idx_relDiscriminator = 0; idx_relDiscriminator < num_relDiscriminatorValues_; ++idx_relDiscriminator )
        {
          double min_relDiscriminator = min_relDiscriminatorValues_[idx_relDiscriminator];
          double max_relDiscriminator = max_relDiscriminatorValues_[idx_relDiscriminator];
          ratePlots_.push_back(new ratePlotEntryType(
            min_absEta, max_absEta,
            min_leadTrackPt, max_leadTrackPt, 
            min_relDiscriminator, max_relDiscriminator, -1., -1.,
            min_dz, max_dz));
        } // idx_relDiscriminator
        for ( size_t idx_absDiscriminator = 0; idx_absDiscriminator < num_absDiscriminatorValues_; ++idx_absDiscriminator )
        {
          double min_absDiscriminator = min_absDiscriminatorValues_[idx_absDiscriminator];
          double max_absDiscriminator = max_absDiscriminatorValues_[idx_absDiscriminator];
          ratePlots_.push_back(new ratePlotEntryType(
            min_absEta, max_absEta,
            min_leadTrackPt, max_leadTrackPt, 
            -1., -1., min_absDiscriminator, max_absDiscriminator, 
            min_dz, max_dz));
        } // idx_absDiscriminator
      } // idx_dz
    } // idx_leadTrackPt
  } // idx_absEta

  for ( auto ratePlot : ratePlots_ ) 
  {
    ratePlot->bookHistograms(dqmStore);
  }
}

namespace
{
  bool
  isHigherPt(const BaseTau& tau1, const BaseTau& tau2)
  {
    return tau1.p4().pt() > tau2.p4().pt();
  }
}

void BaseTauAnalyzerBackground::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  std::vector<BaseTau> baseTaus = buildBaseTaus(evt, es);
  
  // CV: sort taus by decreasing pT
  std::sort(baseTaus.begin(), baseTaus.end(), isHigherPt);

  edm::Handle<double> evtWeight;
  evt.getByToken(token_evtWeight_, evtWeight);

  for ( auto ratePlot : ratePlots_ ) 
  {
    ratePlot->fillHistograms(baseTaus, *evtWeight);
  }
}

void BaseTauAnalyzerBackground::endJob()
{}
