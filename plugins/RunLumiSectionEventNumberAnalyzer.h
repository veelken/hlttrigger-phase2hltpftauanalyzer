#ifndef HLTTrigger_TallinnHLTPFTauAnalyzer_RunLumiSectionEventNumberAnalyzer_h
#define HLTTrigger_TallinnHLTPFTauAnalyzer_RunLumiSectionEventNumberAnalyzer_h

/** \class PrintRunLumiSectionEventNumber
 *
 * Write run + luminosity section + event numbers
 * of processed events to output stream
 *
 * \author Christian Veelken, Tallinn
 *
 */

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <iostream>
#include <string>

class RunLumiSectionEventNumberAnalyzer : public edm::EDAnalyzer
{
 public:
  // constructor 
  explicit RunLumiSectionEventNumberAnalyzer(const edm::ParameterSet&);
    
  // destructor
  virtual ~RunLumiSectionEventNumberAnalyzer();
    
 private:
  void analyze(const edm::Event&, const edm::EventSetup&);

  std::ostream* outputStream_;

  int cfgError_;

  bool isOutputFile_;

  std::string separator_;
};

#endif   
