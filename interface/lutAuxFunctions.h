#ifndef HLTrigger_TallinnHLTPFTauAnalyzer_lutAuxFunctions_h
#define HLTrigger_TallinnHLTPFTauAnalyzer_lutAuxFunctions_h

#include <FWCore/ParameterSet/interface/ParameterSet.h> // edm::ParameterSet

#include "HLTrigger/TallinnHLTPFTauAnalyzer/interface/LocalFileInPath.h" // LocalFileInPath

// forward declarations
class TFile;
class TH1;
class TH2;

// define auxiliary functions
TFile *
openFile(const LocalFileInPath & fileName);

TH1 *
loadTH1(TFile * inputFile,
        const std::string & histogramName);

TH2 *
loadTH2(TFile * inputFile,
        const std::string & histogramName);

#endif // HLTrigger_TallinnHLTPFTauAnalyzer_lutAuxFunctions_h
