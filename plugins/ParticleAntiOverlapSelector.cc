#include "HLTTrigger/TallinnHLTPFTauAnalyzer/plugins/ParticleAntiOverlapSelector.h"

#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

#include "DataFormats/PatCandidates/interface/Tau.h"                 // pat::Tau
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"     // pat::PackedCandidate
#include "DataFormats/TauReco/interface/PFTau.h"                     // reco::PFTau
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h" // reco::PFCandidate
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"        // reco::GenParticle
#include "DataFormats/JetReco/interface/GenJet.h"                    // reco::GenJet

typedef ObjectSelector<ParticleAntiOverlapSelector<pat::Tau>> PATTauAntiOverlapSelector;
typedef ObjectSelector<ParticleAntiOverlapSelector<pat::PackedCandidate>> PackedCandidateAntiOverlapSelector;
typedef ObjectSelector<ParticleAntiOverlapSelector<reco::PFTau>> PFTauAntiOverlapSelector;
typedef ObjectSelector<ParticleAntiOverlapSelector<reco::PFCandidate>> PFCandidateAntiOverlapSelector;
typedef ObjectSelector<ParticleAntiOverlapSelector<reco::GenParticle>> GenParticleAntiOverlapSelector;
typedef ObjectSelector<ParticleAntiOverlapSelector<reco::GenJet>> GenJetAntiOverlapSelector;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(PATTauAntiOverlapSelector);
DEFINE_FWK_MODULE(PackedCandidateAntiOverlapSelector);
DEFINE_FWK_MODULE(PFTauAntiOverlapSelector);
DEFINE_FWK_MODULE(PFCandidateAntiOverlapSelector);
DEFINE_FWK_MODULE(GenParticleAntiOverlapSelector);
DEFINE_FWK_MODULE(GenJetAntiOverlapSelector);
