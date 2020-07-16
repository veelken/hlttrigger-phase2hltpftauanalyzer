#ifndef HLTrigger_TallinnHLTPFTauAnalyzer_RunLumiSectionEventNumberFilter_h
#define HLTrigger_TallinnHLTPFTauAnalyzer_RunLumiSectionEventNumberFilter_h

/** \class RunLumiSectionEventNumberFilter
 *
 * Select events based on run + luminosity section + event number pairs
 * written (a three columns separated by white-space character) into an ASCII file
 * 
 * \author Christian Veelken, Tallinn
 *
 */

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>
#include <map>

class RunLumiSectionEventNumberFilter : public edm::EDFilter
{
 public:
  // constructor 
  explicit RunLumiSectionEventNumberFilter(const edm::ParameterSet&);
    
  // destructor
  virtual ~RunLumiSectionEventNumberFilter();
    
 private:
  void beginJob();

  bool filter(edm::Event&, const edm::EventSetup&);

//--- read ASCII file containing run and event numbers
  void readRunLumiSectionEventNumberFile();
  
  std::string runLumiSectionEventNumberFileName_;

  std::string separator_;

  typedef std::set<edm::EventNumber_t> eventNumberSet;
  typedef std::map<edm::LuminosityBlockNumber_t, eventNumberSet> lumiSectionEventNumberMap;
  std::map<edm::RunNumber_t, lumiSectionEventNumberMap> runLumiSectionEventNumbers_;
  
  typedef std::map<edm::EventNumber_t, int> matchedEventNumbersMap;
  typedef std::map<edm::LuminosityBlockNumber_t, matchedEventNumbersMap> matchedLumiSectionEventNumberMap;
  std::map<edm::RunNumber_t, matchedLumiSectionEventNumberMap> matchedRunLumiSectionEventNumbers_;

  int cfgError_;

  long numEventsProcessed_;
  long numEventsToBeSelected_;
  long numEventsSelected_;
};

#endif   
