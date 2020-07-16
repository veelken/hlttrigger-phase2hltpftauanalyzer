#include "HLTrigger/TallinnHLTPFTauAnalyzer/plugins/RunLumiSectionEventNumberAnalyzer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>
#include <iomanip>
#include <fstream>

namespace
{
  std::ostream* getOutputOptions(const edm::ParameterSet& cfg, bool& isOutputFile, int& error)
  {
    std::ostream* outputStream = 0;    
    if ( cfg.exists("output") ) {
      std::string output = cfg.getParameter<std::string>("output");    
      if ( output == "std::cout" ) {
	outputStream = &std::cout;
	isOutputFile = false;
      } else if ( output == "std::cerr" ) {
	outputStream = &std::cerr;
	isOutputFile = false;
      } else if ( output != "" ) {
	outputStream = new std::ofstream(output.data(), std::ios::out);
	isOutputFile = true;
      } else {
	edm::LogError ("getOutputOptions") 
	  << " Invalid Configuration Parameter output = " << output << " --> skipping !!";
	error = 1;
      }
    } else {
      outputStream = &std::cout;
      isOutputFile = false;
    }    
    return outputStream;
  }
}

RunLumiSectionEventNumberAnalyzer::RunLumiSectionEventNumberAnalyzer(const edm::ParameterSet& cfg)
  : cfgError_(0)
{
  //std::cout << "<RunLumiSectionEventNumberAnalyzer::RunLumiSectionEventNumberAnalyzer>:" << std::endl;

  outputStream_ = getOutputOptions(cfg, isOutputFile_, cfgError_);

  separator_ = cfg.exists("separator") ? 
    cfg.getParameter<std::string>("separator") : ":";
}

RunLumiSectionEventNumberAnalyzer::~RunLumiSectionEventNumberAnalyzer()
{
//--- close output file
  if ( isOutputFile_ ) delete outputStream_;
}

void RunLumiSectionEventNumberAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
//--- check that configuration parameters contain no errors
  if ( cfgError_ ) {
    edm::LogError ("analyze") << " Error in Configuration ParameterSet --> skipping !!";
    return;
  }

//--- retrieve run and event numbers from the event
  edm::RunNumber_t runNumber = evt.id().run();
  edm::LuminosityBlockNumber_t lumiSectionNumber = evt.luminosityBlock();
  edm::EventNumber_t eventNumber = evt.id().event();

  *outputStream_ << runNumber << separator_ << lumiSectionNumber << separator_ << eventNumber << std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(RunLumiSectionEventNumberAnalyzer);

