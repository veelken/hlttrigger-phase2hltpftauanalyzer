#ifndef HLTTrigger_TallinnHLTPFTauAnalyzer_DumpGenParticles_h
#define HLTTrigger_TallinnHLTPFTauAnalyzer_DumpGenParticles_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"    // reco::GenParticle
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h" // reco::GenParticleCollection

#include <vector>
#include <string>

class DumpGenParticles : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit DumpGenParticles(const edm::ParameterSet&);
    
  // destructor
  ~DumpGenParticles();
    
 private:
  void analyze(const edm::Event&, const edm::EventSetup&);

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<reco::GenParticleCollection> token_;
};

#endif   
