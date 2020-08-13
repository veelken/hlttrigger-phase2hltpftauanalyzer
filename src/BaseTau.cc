#include "HLTrigger/TallinnHLTPFTauAnalyzer/interface/BaseTau.h"

BaseTau::BaseTau(const reco::Candidate::LorentzVector& p4, const reco::Candidate::LorentzVector& leadTrackP4, double discriminator, double zVtx)
 : p4_(p4)
 , leadTrackP4_(leadTrackP4)
 , discriminator_(discriminator)
 , zVtx_(zVtx)
{}

BaseTau::~BaseTau()
{}

const reco::Candidate::LorentzVector& 
BaseTau::p4() const
{
  return p4_;
}

const reco::Candidate::LorentzVector& 
BaseTau::leadTrackP4() const
{
  return leadTrackP4_;
}

double 
BaseTau::discriminator() const
{
  return discriminator_;
}

double 
BaseTau::zVtx() const
{
  return zVtx_;
}
