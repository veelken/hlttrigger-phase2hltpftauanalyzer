#ifndef HLTrigger_TallinnHLTPFTauAnalyzer_BaseTau_h
#define HLTrigger_TallinnHLTPFTauAnalyzer_BaseTau_h

#include "DataFormats/Candidate/interface/Candidate.h"

class BaseTau
{
 public:
  BaseTau(const reco::Candidate::LorentzVector& p4, const reco::Candidate::LorentzVector& leadTrackP4, double discriminator, double zVtx);
  ~BaseTau();

  const reco::Candidate::LorentzVector& p4() const;
  const reco::Candidate::LorentzVector& leadTrackP4() const;
  double discriminator() const;
  double zVtx() const;

 private:
  reco::Candidate::LorentzVector p4_;
  reco::Candidate::LorentzVector leadTrackP4_;
  double discriminator_;
  double zVtx_;
};

#endif
