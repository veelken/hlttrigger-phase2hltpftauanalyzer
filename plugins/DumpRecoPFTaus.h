#ifndef HLTTrigger_TallinnHLTPFTauAnalyzer_DumpTallinRecoPFTaus_h
#define HLTTrigger_TallinnHLTPFTauAnalyzer_DumpTallinRecoPFTaus_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TauReco/interface/PFTau.h"              // reco::PFTau
#include "DataFormats/TauReco/interface/PFTauFwd.h"           // reco::PFTauCollection
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h" // reco::PFTauDiscriminator

#include <vector>
#include <string>

class DumpRecoPFTaus : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit DumpRecoPFTaus(const edm::ParameterSet&);
    
  // destructor
  ~DumpRecoPFTaus();
    
 private:
  void analyze(const edm::Event&, const edm::EventSetup&);

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<reco::PFTauCollection> token_;
  edm::InputTag src_sumChargedIso_;
  edm::EDGetTokenT<reco::PFTauDiscriminator> token_sumChargedIso_;
  std::vector<edm::InputTag> src_discriminators_;
  std::vector<edm::EDGetTokenT<reco::PFTauDiscriminator>> token_discriminators_;

  bool debug_;
};

#endif   
