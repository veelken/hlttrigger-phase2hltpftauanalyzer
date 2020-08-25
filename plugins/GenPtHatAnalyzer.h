#ifndef HLTrigger_TallinnHLTPFTauAnalyzer_GenPtHatAnalyzer_h
#define HLTrigger_TallinnHLTPFTauAnalyzer_GenPtHatAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include <TH1.h>     // TH1
#include <TH2.h>     // TH2
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

  edm::InputTag src_genEventInfo_;
  edm::EDGetTokenT<GenEventInfoProduct> token_genEventInfo_;

  edm::InputTag src_genJets_;
  edm::EDGetTokenT<reco::GenJetCollection> token_genJets_;

  edm::InputTag src_pileupSummaryInfo_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> token_pileupSummaryInfo_;

  double lumiScale_;

  std::string dqmDirectory_;

  MonitorElement* me_genPtHat_hardscatter_;
  TH1* histogram_genPtHat_hardscatter_;

  MonitorElement* me_leadGenJetPt_vs_genPtHat_;
  TH2* histogram_leadGenJetPt_vs_genPtHat_;

  MonitorElement* me_subleadGenJetPt_vs_genPtHat_;
  TH2* histogram_subleadGenJetPt_vs_genPtHat_;

  MonitorElement* me_genPtHat_pileup_;
  TH1* histogram_genPtHat_pileup_;
};

#endif   

