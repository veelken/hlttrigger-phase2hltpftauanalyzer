#ifndef HLTTrigger_TallinnHLTPFTauAnalyzer_RhoAnalyzer_h
#define HLTTrigger_TallinnHLTPFTauAnalyzer_RhoAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include <TH1.h>     // TH1
#include <TString.h> // TString, Form()

#include <vector>    // std::vector
#include <string>    // std::string

using namespace dqm::implementation;

class RhoAnalyzer : public edm::EDAnalyzer 
{
 public:
  RhoAnalyzer(const edm::ParameterSet& cfg);
  ~RhoAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag src_rho_;
  edm::EDGetTokenT<double> token_rho_;

  edm::InputTag src_rhoNeutral_;
  edm::EDGetTokenT<double> token_rhoNeutral_;

  std::string dqmDirectory_;

  MonitorElement* me_rho_;
  TH1* histogram_rho_;
  
  MonitorElement* me_rhoNeutral_;
  TH1* histogram_rhoNeutral_;
};

#endif   

