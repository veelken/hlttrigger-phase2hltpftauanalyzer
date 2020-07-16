#include "HLTrigger/TallinnHLTPFTauAnalyzer/plugins/ParticleAntiOverlapSelector.h"

#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

#include "DataFormats/PatCandidates/interface/Tau.h"                 // pat::Tau
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"     // pat::PackedCandidate
#include "DataFormats/TauReco/interface/PFTau.h"                     // reco::PFTau
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h" // reco::PFCandidate
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"        // reco::GenParticle
#include "DataFormats/JetReco/interface/GenJet.h"                    // reco::GenJet
#include "DataFormats/PatCandidates/interface/Jet.h"                 // pat::Jet
#include "DataFormats/JetReco/interface/PFJet.h"                     // reco::PFJet

typedef ObjectSelector<ParticleAntiOverlapSelector<pat::Tau>> PATTauAntiOverlapSelector;
typedef ObjectSelector<ParticleAntiOverlapSelector<pat::PackedCandidate>> PackedCandidateAntiOverlapSelector;
typedef ObjectSelector<ParticleAntiOverlapSelector<reco::PFTau>> PFTauAntiOverlapSelector;
typedef ObjectSelector<ParticleAntiOverlapSelector<reco::PFCandidate>> PFCandidateAntiOverlapSelector;
typedef ObjectSelector<ParticleAntiOverlapSelector<reco::GenParticle>> GenParticleAntiOverlapSelector;
typedef ObjectSelector<ParticleAntiOverlapSelector<reco::GenJet>> GenJetAntiOverlapSelector;
typedef ObjectSelector<ParticleAntiOverlapSelector<pat::Jet>> PATJetAntiOverlapSelector;
typedef ObjectSelector<ParticleAntiOverlapSelector<reco::PFJet>> PFJetAntiOverlapSelector;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(PATTauAntiOverlapSelector);
DEFINE_FWK_MODULE(PackedCandidateAntiOverlapSelector);
DEFINE_FWK_MODULE(PFTauAntiOverlapSelector);
DEFINE_FWK_MODULE(PFCandidateAntiOverlapSelector);
DEFINE_FWK_MODULE(GenParticleAntiOverlapSelector);
DEFINE_FWK_MODULE(GenJetAntiOverlapSelector);
DEFINE_FWK_MODULE(PATJetAntiOverlapSelector);
DEFINE_FWK_MODULE(PFJetAntiOverlapSelector);
