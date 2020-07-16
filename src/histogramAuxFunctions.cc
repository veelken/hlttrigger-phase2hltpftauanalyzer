#include "HLTrigger/TallinnHLTPFTauAnalyzer/interface/histogramAuxFunctions.h"

void
fill(TH1* histogram,
     double x,
     double evtWeight)
{
  const TAxis* xAxis = histogram->GetXaxis();
  int idxBin = xAxis->FindBin(x);
  int numBins = xAxis->GetNbins();

  if( idxBin >= 1 && idxBin <= numBins )
  {
    histogram->Fill(x, evtWeight);
  }
}

void
fillWithOverFlow(TH1 * histogram,
                 double x,
                 double evtWeight)
{
  const TAxis* xAxis = histogram->GetXaxis();
  int idxBin = xAxis->FindBin(x);
  int numBins = xAxis->GetNbins();
  if ( idxBin < 1       ) idxBin = 1;
  if ( idxBin > numBins ) idxBin = numBins;
  double binCenter = xAxis->GetBinCenter(idxBin);

  histogram->Fill(binCenter, evtWeight);
}

void
fill2d(TH2* histogram,
       double x,
       double y,
       double evtWeight)
{
  const TAxis* xAxis = histogram->GetXaxis();
  int idxBinX = xAxis->FindBin(x);
  int numBinsX = xAxis->GetNbins();

  const TAxis* yAxis = histogram->GetYaxis();
  int idxBinY = yAxis->FindBin(y);
  int numBinsY = yAxis->GetNbins();

  if ( idxBinX >= 1 && idxBinX <= numBinsX && idxBinY >= 1 && idxBinY <= numBinsY )
  {
    histogram->Fill(x, y, evtWeight);
  }  
}

void
fillWithOverFlow2d(TH2 * histogram,
                   double x,
                   double y,
                   double evtWeight)
{
  const TAxis* xAxis = histogram->GetXaxis();
  int idxBinX = xAxis->FindBin(x);
  int numBinsX = xAxis->GetNbins();
  if ( idxBinX < 1        ) idxBinX = 1;
  if ( idxBinX > numBinsX ) idxBinX = numBinsX;
  double binCenterX = xAxis->GetBinCenter(idxBinX);

  const TAxis* yAxis = histogram->GetYaxis();
  int idxBinY = yAxis->FindBin(y);
  int numBinsY = yAxis->GetNbins();
  if ( idxBinY < 1        ) idxBinY = 1;
  if ( idxBinY > numBinsY ) idxBinY = numBinsY;
  double binCenterY = yAxis->GetBinCenter(idxBinY);

  histogram->Fill(binCenterX, binCenterY, evtWeight);
}

void divideByBinWidth(TH1* histogram)
{
  if ( !histogram ) return;
  TAxis* xAxis = histogram->GetXaxis();
  int numBins = xAxis->GetNbins();
  for ( int iBin = 1; iBin <= numBins; ++iBin ) {
    double binContent = histogram->GetBinContent(iBin);
    double binError = histogram->GetBinError(iBin);
    double binWidth = xAxis->GetBinWidth(iBin);
    histogram->SetBinContent(iBin, binContent/binWidth);
    histogram->SetBinError(iBin, binError/binWidth);
  }
}
