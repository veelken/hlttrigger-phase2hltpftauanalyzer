#ifndef HLTTrigger_TallinnHLTPFTauAnalyzer_GenChargedHadronToTrackMatch_h
#define HLTTrigger_TallinnHLTPFTauAnalyzer_GenChargedHadronToTrackMatch_h

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" // reco::GenParticle
#include "DataFormats/TrackReco/interface/Track.h"            // reco::Track

class GenChargedHadronToTrackMatchBase
{
 public:
  GenChargedHadronToTrackMatchBase(const reco::Candidate* genChargedHadron);
  virtual ~GenChargedHadronToTrackMatchBase();

  const reco::Candidate* genChargedHadron() const;
  bool hasGenChargedHadron() const;
  double genChargedHadron_pt() const;
  double genChargedHadron_eta() const;
  double genChargedHadron_absEta() const;
  double genChargedHadron_phi() const;

  bool hasRecTrack() const;
  double recTrack_pt() const;
  double recTrack_eta() const;
  double recTrack_absEta() const;
  double recTrack_phi() const;

 protected:
  const reco::Candidate* genChargedHadron_;
  bool hasGenChargedHadron_;
  double genChargedHadron_pt_;
  double genChargedHadron_eta_;
  double genChargedHadron_absEta_;
  double genChargedHadron_phi_;
  
  bool hasRecTrack_;
  double recTrack_pt_;
  double recTrack_eta_;
  double recTrack_absEta_;
  double recTrack_phi_;
};

class GenChargedHadronToRecoTrackMatch : public GenChargedHadronToTrackMatchBase
{
 public:
  GenChargedHadronToRecoTrackMatch(const reco::Candidate* genChargedHadron, const reco::Track* recTrack);
  ~GenChargedHadronToRecoTrackMatch();
  
  const reco::Track* recTrack() const;

  bool isOverlap(const GenChargedHadronToTrackMatchBase* other);

 private:
  const reco::Track* recTrack_;
};

bool isOverlap(const GenChargedHadronToRecoTrackMatch& match1, const GenChargedHadronToRecoTrackMatch& match2);

#endif
