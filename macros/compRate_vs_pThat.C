
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <math.h>
#include <limits>

TH1* loadHistogram(TFile* inputFile, const std::string& histogramName)
{
  TH1* histogram = dynamic_cast<TH1*>(inputFile->Get(histogramName.data()));
  if ( !histogram ) {
    std::cerr << "Failed to load histogram = " << histogramName << " from file = " << inputFile->GetName() << " !!" << std::endl;
    assert(0);
  }
  if ( !histogram->GetSumw2N() ) histogram->Sumw2();
  return histogram;
}

double compIntegral(TH1* histogram, double minX, double maxX)
{
  const TAxis* xAxis = histogram->GetXaxis();
  int numBins = xAxis->GetNbins();
  double integral = 0.;
  for ( int idxBin = 1; idxBin <= numBins; ++idxBin ) 
  {
    double binCenter = xAxis->GetBinCenter(idxBin);
    double binContent = histogram->GetBinContent(idxBin);
    if ( binCenter > minX && binCenter < maxX ) 
    {
      integral += binContent;
    }
  }
  return integral;
}

void compRate_vs_pThat()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  std::string inputFilePath = "/hdfs/local/veelken/Phase2HLT/rate/2020Jul29/";
  std::string inputFileName = "hadd_QCD_all.root";
  std::string inputFileName_full = inputFilePath;
  if ( inputFileName_full.find_last_of("/") != (inputFileName_full.size() - 1) ) inputFileName_full.append("/");
  inputFileName_full.append(inputFileName);
  TFile* inputFile = new TFile(inputFileName_full.data());
  if ( !inputFile ) {
    std::cerr << "Failed to open input file = '" << inputFileName_full << "' !!" << std::endl;
    assert(0);
  }

  std::string histogramName_pThat = "GenPtHatAnalyzer/genPtHat";
  TH1* histogram_pThat = loadHistogram(inputFile, histogramName_pThat);

  std::vector<double> binning_pThat = { 30., 50., 80., 120., 170., 300., 470. };
  int numBins_pThat = binning_pThat.size() - 1;
  for ( int idxBin = 0; idxBin < numBins_pThat; ++ idxBin ) 
  {
    double min_pThat = binning_pThat[idxBin];
    double max_pThat = binning_pThat[idxBin + 1];
    double rate = compIntegral(histogram_pThat, min_pThat, max_pThat);
    std::cout << "rate (" << min_pThat << " < pT_hat < " << max_pThat << " GeV): " << rate << " Hz" << std::endl;
  }

  delete inputFile;
}
