#ifndef HLTTrigger_TallinnHLTPFTauAnalyzer_GenVertexProducer_h
#define HLTTrigger_TallinnHLTPFTauAnalyzer_GenVertexProducer_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"    // reco::GenParticle
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h" // reco::GenParticleCollection

#include <vector>

class GenVertexProducer : public edm::EDProducer 
{
 public:
  explicit GenVertexProducer(const edm::ParameterSet& cfg);
  ~GenVertexProducer();

 private:
  void produce(edm::Event& evt, const edm::EventSetup& es);

  edm::InputTag src_;
  edm::EDGetTokenT<reco::GenParticleCollection> token_;

  typedef std::vector<int> vint;
  vint pdgIds_;
};

#endif
