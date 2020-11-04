#ifndef HLTrigger_TallinnHLTPFTauAnalyzer_EvtWeightAnalyzer_h
#define HLTrigger_TallinnHLTPFTauAnalyzer_EvtWeightAnalyzer_h

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

class EvtWeightAnalyzer : public edm::EDAnalyzer 
{
 public:
  EvtWeightAnalyzer(const edm::ParameterSet& cfg);
  ~EvtWeightAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag src_;
  edm::EDGetTokenT<double> token_;

  int numBins_evtWeight_;
  double xMin_evtWeight_;
  double xMax_evtWeight_;

  int numBins_log10evtWeight_;
  double xMin_log10evtWeight_;
  double xMax_log10evtWeight_;

  std::string dqmDirectory_;

  MonitorElement* me_evtWeight_;
  TH1* histogram_evtWeight_;

  MonitorElement* me_log10evtWeight_;
  TH1* histogram_log10evtWeight_;
};

#endif   

