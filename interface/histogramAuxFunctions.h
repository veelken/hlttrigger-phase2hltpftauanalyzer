#ifndef HLTrigger_TallinnHLTPFTauAnalyzer_histogramAuxFunctions_h
#define HLTrigger_TallinnHLTPFTauAnalyzer_histogramAuxFunctions_h

#include <TH1.h> // TH1
#include <TH2.h> // TH2

void
fill(TH1* histogram,
     double x,
     double evtWeight);

void
fillWithOverFlow(TH1* histogram,
                 double x,
                 double evtWeight);

void
fill2d(TH2* histogram,
       double x,
       double y,
       double evtWeight);

void
fillWithOverFlow2d(TH2* histogram,
                   double x,
                   double y,
                   double evtWeight);

void 
divideByBinWidth(TH1* histogram);

#endif // HLTrigger_TallinnHLTPFTauAnalyzer_histogramAuxFunctions_h
