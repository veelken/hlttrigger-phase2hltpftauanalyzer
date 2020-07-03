#include "HLTTrigger/TallinnHLTPFTauAnalyzer/plugins/DumpRecoPFTaus.h"

#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

DumpRecoPFTaus::DumpRecoPFTaus(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
  , debug_(cfg.getUntrackedParameter<bool>("debug", false))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  token_ = consumes<reco::PFTauCollection>(src_);
  src_sumChargedIso_ = cfg.getParameter<edm::InputTag>("src_sumChargedIso");
  token_sumChargedIso_ = consumes<reco::PFTauDiscriminator>(src_sumChargedIso_);
  src_discriminators_ = cfg.getParameter<std::vector<edm::InputTag>>("src_discriminators");
  for ( std::vector<edm::InputTag>::const_iterator src_discriminator = src_discriminators_.begin();
        src_discriminator != src_discriminators_.end(); ++src_discriminator ) {
    token_discriminators_.push_back(consumes<reco::PFTauDiscriminator>(*src_discriminator));
  }
}

DumpRecoPFTaus::~DumpRecoPFTaus()
{}

namespace
{
  bool isSelected(const reco::PFTau& tau, double sumChargedIso, const std::map<std::string, double>& discriminator_values, bool debug)
  {
    const double min_PFTau_pt                   = 20.;
    const double max_PFTau_absEta               =  2.4;
    const double min_leadPFChargedHadron_pt     =  5.;
    const double max_leadPFChargedHadron_absEta =  2.4;
    const double max_leadPFChargedHadron_dz     =  0.2;
    const double max_chargedIso                 = 1.e+3;
    const double max_chargedRelIso              = 0.05;
    bool retVal = true;
    if ( !(tau.pt() > min_PFTau_pt) )
    {
      if ( debug )
      {
	std::cout << " FAILS min_PFTau_pt cut." << std::endl;
      }
      retVal = false;
    }
    if ( !(std::fabs(tau.eta()) < max_PFTau_absEta) )
    {
      if ( debug )
      {
	std::cout << " FAILS max_PFTau_absEta cut." << std::endl;
      }
      retVal = false;
    }
    if ( !(tau.leadPFChargedHadrCand().isNonnull() && 
           tau.leadPFChargedHadrCand()->pt() > min_leadPFChargedHadron_pt) )
    {
      if ( debug )
      {
	std::cout << " FAILS min_leadPFChargedHadron_pt cut." << std::endl;
      }
      retVal = false;
    }
    if ( !(tau.leadPFChargedHadrCand().isNonnull() &&
	   std::fabs(tau.leadPFChargedHadrCand()->eta()) < max_leadPFChargedHadron_absEta) )
    {
      if ( debug )
      {
	std::cout << " FAILS max_leadPFChargedHadron_absEta cut." << std::endl;
      }
      retVal = false;
    }
    if ( !(tau.leadPFChargedHadrCand().isNonnull() && tau.leadPFChargedHadrCand()->bestTrack() &&
           std::fabs(tau.leadPFChargedHadrCand()->bestTrack()->dz(tau.leadPFChargedHadrCand()->vertex())) < max_leadPFChargedHadron_dz) )
    { 
      if ( debug )
      {
	std::cout << " FAILS max_leadPFChargedHadron_dz cut." << std::endl;
      }
      retVal = false;
    }
    if ( !(sumChargedIso < max_chargedIso) )
    {
      if ( debug )
      {
	std::cout << " FAILS max_chargedIso cut." << std::endl;
      }
      retVal = false;
    }
    if ( !(sumChargedIso < max_chargedRelIso*tau.pt()) )
    {
      if ( debug )
      {
	std::cout << " FAILS max_chargedRelIso cut." << std::endl;
      }
      retVal = false;
    }
    for ( std::map<std::string, double>::const_iterator discriminator_value = discriminator_values.begin();
          discriminator_value != discriminator_values.end(); ++discriminator_value ) {
      if ( discriminator_value->second < 0.5 )
      {
	std::cout << " FAILS '" << discriminator_value->first << "' discriminator." << std::endl;
      }
      retVal = false;
    }
    if ( debug && retVal )
    {
      std::cout << " PASSES selection." << std::endl;
    }
    return retVal;
  }
}

void DumpRecoPFTaus::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  std::cout << "<DumpRecoPFTaus::analyze (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
  std::cout << " src = " << src_ << std::endl;

  edm::Handle<reco::PFTauCollection> taus;
  evt.getByToken(token_, taus);
  
  edm::Handle<reco::PFTauDiscriminator> sumChargedIso;
  evt.getByToken(token_sumChargedIso_, sumChargedIso);

  size_t numTaus = taus->size();
  for ( size_t idxTau = 0; idxTau < numTaus; ++idxTau ) 
  {
    reco::PFTauRef tauRef(taus, idxTau);
    double sumChargedIso_value = (*sumChargedIso)[tauRef];
    std::cout << "PFTau #" << idxTau << ": " << " pT = " << tauRef->pt() << ", eta = " << tauRef->eta() << ", phi = " << tauRef->phi() << ","
              << " decayMode = " << tauRef->decayMode() << ", mass = " << tauRef->mass() << ", chargedIso = " << sumChargedIso_value << std::endl;
    std::cout << "lead. ChargedPFCand:";
    if ( tauRef->leadPFChargedHadrCand().isNonnull() )
    {
      std::cout << " pT = " << tauRef->leadPFChargedHadrCand()->pt() << ", eta = " << tauRef->leadPFChargedHadrCand()->eta() << ", phi = " << tauRef->leadPFChargedHadrCand()->phi();
    }
    else 
    {
      std::cout << " N/A";
    }
    std::cout << std::endl;
    std::cout << "lead. Track:";
    if ( tauRef->leadPFChargedHadrCand().isNonnull() && tauRef->leadPFChargedHadrCand()->bestTrack() ) 
    {
      const reco::Track* leadingTrack = tauRef->leadPFChargedHadrCand()->bestTrack();
      std::cout << " pT = " << leadingTrack->pt() << ", eta = " << leadingTrack->eta() << ", phi = " << leadingTrack->phi();
    }
    else 
    {
      std::cout << " N/A";
    }
    std::cout << std::endl;
    std::map<std::string, double> discriminator_values;
    size_t numDiscriminators = token_discriminators_.size();
    for ( size_t idxDiscriminator = 0; idxDiscriminator < numDiscriminators; ++idxDiscriminator ) 
    {
      edm::Handle<reco::PFTauDiscriminator> discriminators;
      evt.getByToken(token_discriminators_[idxDiscriminator], discriminators);
      discriminator_values[src_discriminators_[idxDiscriminator].label()] = (*discriminators)[tauRef];
    }
    if ( debug_ ) 
    {
      if ( !isSelected(*tauRef, sumChargedIso_value, discriminator_values, debug_) )
      {
	std::cout << "--> CHECK !!" << std::endl;
      }
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DumpRecoPFTaus);





