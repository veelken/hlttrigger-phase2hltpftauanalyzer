#ifndef HLTTrigger_HLTPFTauAnalyzer_RecoMEtResolutionAnalyzer_h
#define HLTTrigger_HLTPFTauAnalyzer_RecoMEtResolutionAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h" 
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/METReco/interface/PFMET.h"     // reco::PFMET
#include "DataFormats/METReco/interface/PFMETFwd.h"  // reco::PFMETCollection
#include "DataFormats/METReco/interface/GenMET.h"    // reco::GenMET
#include "DataFormats/METReco/interface/GenMETFwd.h" // reco::GenMETCollection

#include <TH1.h>     // TH1
#include <TString.h> // TString, Form()
#include <TMath.h>   // TMath::Abs(), TMath::Pi()

#include <vector>    // std::vector
#include <string>    // std::string

using namespace dqm::implementation;

class RecoMEtResolutionAnalyzer : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit RecoMEtResolutionAnalyzer(const edm::ParameterSet&);
    
  // destructor
  ~RecoMEtResolutionAnalyzer();
    
 private:
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag srcRecMEt_;
  edm::EDGetTokenT<reco::PFMETCollection> tokenRecMEt_;
  edm::InputTag srcGenMEt_;
  edm::EDGetTokenT<reco::GenMETCollection> tokenGenMEt_;

  std::string dqmDirectory_;

  MonitorElement* me_recMEt_Pt_;
  TH1* histogram_recMEt_Pt_;
  MonitorElement* me_genMEt_Pt_;
  TH1* histogram_genMEt_Pt_;
  MonitorElement* me_deltaMEt_Px_;
  TH1* histogram_deltaMEt_Px_;
  MonitorElement* me_deltaMEt_Py_;
  TH1* histogram_deltaMEt_Py_;
};

#endif   
