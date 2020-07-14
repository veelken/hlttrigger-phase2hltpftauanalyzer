#ifndef HLTTrigger_TallinnHLTPFTauAnalyzer_GenPtHatAnalyzer_h
#define HLTTrigger_TallinnHLTPFTauAnalyzer_GenPtHatAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include <TH1.h>     // TH1
#include <TString.h> // TString, Form()

#include <vector>    // std::vector
#include <string>    // std::string

using namespace dqm::implementation;

class GenPtHatAnalyzer : public edm::EDAnalyzer 
{
 public:
  GenPtHatAnalyzer(const edm::ParameterSet& cfg);
  ~GenPtHatAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<GenEventInfoProduct> token_;

  double lumiScale_;

  std::string dqmDirectory_;

  MonitorElement* me_genPtHat_;
  TH1* histogram_genPtHat_;
};

#endif   

