
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <math.h>
#include <limits>

TFile* openFile(const std::string& inputFilePath, const std::string& inputFileName)
{
  std::string inputFileName_full = inputFilePath;
  if ( inputFileName_full.find_last_of("/") != (inputFileName_full.size() - 1) ) inputFileName_full.append("/");
  inputFileName_full.append(inputFileName);
  TFile* inputFile = new TFile(inputFileName_full.data());
  if ( !inputFile ) {
    std::cerr << "Failed to open input file = '" << inputFileName_full << "' !!" << std::endl;
    assert(0);
  }
  return inputFile;
}

double compIntegral(const TH1* histogram, bool includeUnderflowBin = false, bool includeOverflowBin = false)
{
  int firstBin = 1;
  if ( includeUnderflowBin ) firstBin -= 1;
  int lastBin = histogram->GetNbinsX();
  if ( includeOverflowBin ) lastBin += 1;
  double integral = 0.;
  for ( int idxBin = firstBin; idxBin <= lastBin; ++idxBin ) {
    integral += histogram->GetBinContent(idxBin);
  }
  return integral;
}

TH1* loadHistogram(TFile* inputFile, const std::string& histogramName)
{
  TH1* histogram = dynamic_cast<TH1*>(inputFile->Get(histogramName.data()));
  if ( !histogram ) {
    std::cerr << "Failed to load histogram = " << histogramName << " from file = " << inputFile->GetName() << " !!" << std::endl;
    assert(0);
  }
  if ( !histogram->GetSumw2N() ) histogram->Sumw2();
  double integral = compIntegral(histogram);
  if ( integral > 0. ) {
    histogram->Scale(1./integral);
  }
  return histogram;
}

void showHistograms(double canvasSizeX, double canvasSizeY,
                    TH1* histogram1, const std::string& legendEntry1,
                    TH1* histogram2, const std::string& legendEntry2,
                    TH1* histogram3, const std::string& legendEntry3,
                    TH1* histogram4, const std::string& legendEntry4,
                    TH1* histogram5, const std::string& legendEntry5,
                    TH1* histogram6, const std::string& legendEntry6,
                    int colors[], int lineStyles[], 
                    double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, 
                    std::vector<std::string>& labelTextLines, double labelTextSize,
                    double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
                    double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
                    bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
                    const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  
  canvas->SetTopMargin(0.05);
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.14);
  canvas->SetRightMargin(0.05);

  canvas->SetLogy(useLogScale);
  
  canvas->SetGridx(1);
  canvas->SetGridy(1);

  if ( !histogram1 ) {
    std::cerr << "<showHistograms>: histogram1 = NULL --> skipping !!" << std::endl;
    return;
  }

  histogram1->SetTitle("");
  histogram1->SetStats(false);
  histogram1->SetMinimum(yMin);
  histogram1->SetMaximum(yMax);

  TAxis* xAxis = histogram1->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleSize(0.045);
  xAxis->SetTitleOffset(xAxisOffset);  
  if ( xMax > xMin ) {
    std::cout << "limiting x-axis range to " << xMin << ".." << xMax << std::endl;
    xAxis->SetRangeUser(xMin, xMax);
  }

  TAxis* yAxis = histogram1->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleSize(0.045);
  yAxis->SetTitleOffset(yAxisOffset);

  histogram1->SetLineColor(colors[0]);
  histogram1->SetLineWidth(2);
  histogram1->SetLineStyle(lineStyles[0]);
  histogram1->Draw("hist");

  if ( histogram2 ) {
    histogram2->SetLineColor(colors[1]);
    histogram2->SetLineWidth(2);
    histogram2->SetLineStyle(lineStyles[1]);
    histogram2->Draw("histsame");
  }

  if ( histogram3 ) {
    histogram3->SetLineColor(colors[2]);
    histogram3->SetLineWidth(2);
    histogram3->SetLineStyle(lineStyles[2]);
    histogram3->Draw("histsame");
  }

  if ( histogram4 ) {
    histogram4->SetLineColor(colors[3]);
    histogram4->SetLineWidth(2);
    histogram4->SetLineStyle(lineStyles[3]);
    histogram4->Draw("histsame");
  }

  if ( histogram5 ) {
    histogram5->SetLineColor(colors[4]);
    histogram5->SetLineWidth(2);
    histogram5->SetLineStyle(lineStyles[4]);
    histogram5->Draw("histsame");
  }

  if ( histogram6 ) {
    histogram6->SetLineColor(colors[5]);
    histogram6->SetLineWidth(2);
    histogram6->SetLineStyle(lineStyles[5]);
    histogram6->Draw("histsame");
  }

  TLegend* legend = 0;
  if ( legendEntry1 != "" ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, "", "brNDC"); 
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetTextSize(legendTextSize);
    legend->AddEntry(histogram1, legendEntry1.data(), "l");
    if ( histogram2 ) legend->AddEntry(histogram2, legendEntry2.data(), "l");
    if ( histogram3 ) legend->AddEntry(histogram3, legendEntry3.data(), "l");
    if ( histogram4 ) legend->AddEntry(histogram4, legendEntry4.data(), "l");
    if ( histogram5 ) legend->AddEntry(histogram5, legendEntry5.data(), "l");
    if ( histogram6 ) legend->AddEntry(histogram6, legendEntry6.data(), "l");
    legend->Draw();
  }

  TPaveText* label = 0;
  if ( labelTextLines.size() > 0 ) {
    label = new TPaveText(labelPosX, labelPosY, labelPosX + labelSizeX, labelPosY + labelSizeY, "brNDC");
    for ( std::vector<std::string>::const_iterator labelTextLine = labelTextLines.begin();
          labelTextLine != labelTextLines.end(); ++labelTextLine ) {
      label->AddText(labelTextLine->data());
    }
    label->SetFillColor(10);
    label->SetBorderSize(0);
    label->SetTextColor(1);
    label->SetTextAlign(12);
    label->SetTextSize(labelTextSize);
    label->Draw();
  }

  histogram1->Draw("axissame");

  canvas->Update();
  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  
  delete label;
  delete legend;
  delete canvas;  
}

double compIntegral_within_Range(const TH1* histogram, double xMin, double xMax)
{
  const TAxis* xAxis = histogram->GetXaxis();
  double integral = 0.;
  int numBinsX = xAxis->GetNbins();
  for ( int idxBinX = 1; idxBinX <= numBinsX; ++idxBinX ) { 
    double binCenter = xAxis->GetBinCenter(idxBinX);
    double binContent = histogram->GetBinContent(idxBinX);
    if ( binCenter >= xMin && binCenter < xMax ) {
      integral += binContent;
    }
  }
  return integral;
}

void makeVertexPlots()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  std::string inputFilePath = Form("%s/src/HLTTrigger/TallinnHLTPFTauAnalyzer/test/", gSystem->Getenv("CMSSW_BASE"));
  std::string inputFileName_signal = "analyzePFTaus_signal_2020Jun26.root";
  TFile* inputFile_signal = openFile(inputFilePath, inputFileName_signal);
  std::string inputFileName_background = "analyzePFTaus_background_2020Jun26.root";
  //TFile* inputFile_background = openFile(inputFilePath, inputFileName_background);
  TFile* inputFile_background = nullptr;

  std::vector<std::string> vertexAlgos;
  vertexAlgos.push_back("offlinePrimaryVertices");
  vertexAlgos.push_back("hltPhase2PixelVertices");
  vertexAlgos.push_back("hltPhase2TrimmedPixelVertices");

  std::vector<std::string> observables;
  observables.push_back("delta_z");
  observables.push_back("delta_z_bestMatch");
  observables.push_back("idxRecVertex_bestMatch");

  std::map<std::string, int> rebin; // key = observable
  rebin["delta_z"]                               =   1;
  rebin["delta_z_bestMatch"]                     =   1;
  rebin["idxRecVertex_bestMatch"]                =   1;

  std::map<std::string, double> xMin; // key = observable
  xMin["delta_z"]                                = -10.;
  xMin["delta_z_bestMatch"]                      = -10.;
  xMin["idxRecVertex_bestMatch"]                 =  -0.5;

  std::map<std::string, double> xMax; // key = observable
  xMax["delta_z"]                                = +10.;
  xMax["delta_z_bestMatch"]                      = +10.;
  xMax["idxRecVertex_bestMatch"]                 =   9.5;
  
  std::map<std::string, std::string> xAxisTitles; // key = observable
  xAxisTitles["delta_z"]                         = "z_{vtx}^{rec[0]} - z_{vtx}^{gen} [cm]";
  xAxisTitles["delta_z_bestMatch"]               = "min_{i}(z_{vtx}^{rec[i]} - z_{vtx}^{gen}) [GeV]";
  xAxisTitles["idxRecVertex_bestMatch"]          = "i";

  std::map<std::string, std::string> dqmDirectories; // key = vertexAlgo
  dqmDirectories["offlinePrimaryVertices"]        = "RecoVertexAnalyzer/offlinePrimaryVertices";
  dqmDirectories["hltPhase2PixelVertices"]        = "RecoVertexAnalyzer/hltPhase2PixelVertices";
  dqmDirectories["hltPhase2TrimmedPixelVertices"] = "RecoVertexAnalyzer/hltPhase2TrimmedPixelVertices";  
  
  int colors[6]       = {  1,  2,  8,  4,  6,  7 };
  int lineStyles[6]   = {  1,  1,  1,  1,  1,  1 };
  int markerStyles[6] = { 22, 32, 20, 24, 21, 25 };

  enum { kSignal, kBackground };
  for ( int idxSignal_or_Background = kSignal; idxSignal_or_Background <= kSignal; ++idxSignal_or_Background ) { // CV: signal only
    TFile* inputFile = nullptr;
    if      ( idxSignal_or_Background == kSignal     ) inputFile = inputFile_signal;
    else if ( idxSignal_or_Background == kBackground ) inputFile = inputFile_background;
    else assert(0);
    for ( std::vector<std::string>::const_iterator observable = observables.begin();
	  observable != observables.end(); ++observable ) {
      std::map<std::string, TH1*> histograms;          // key = vertexAlgo
      std::map<std::string, TH1*> histograms_rebinned; // key = vertexAlgo
      for ( std::vector<std::string>::const_iterator vertexAlgo = vertexAlgos.begin();
	    vertexAlgo != observables.end(); ++vertexAlgo ) {
        std::string histogramName = Form("%s/%s", dqmDirectories[*vertexAlgo].data(), observable->data()); 
        TH1* histogram = loadHistogram(inputFile, histogramName);
        histograms[*vertexAlgo] = histogram;
        TH1* histogram_rebinned = ( rebin[*observable] > 1 ) ? histogram->Rebin(rebin[*observable]) : histogram;
        histograms_rebinned[*vertexAlgo] = histogram_rebinned;
      
        if ( (*observable) == "delta_z" ) {
          std::cout << "Vertex collection = '" << (*vertexAlgo) << "':" << std::endl;
	  std::cout << " Efficiency for dz < 0.2 cm = " << compIntegral_within_Range(histogram, -0.2, +0.2) << std::endl;
          std::cout << " Efficiency for dz < 0.4 cm = " << compIntegral_within_Range(histogram, -0.4, +0.4) << std::endl;
        }
      }

      std::vector<std::string> labelTextLines;
      std::string outputFileName;
      if      ( idxSignal_or_Background == kSignal     ) outputFileName = Form("makeVertexPlots_signal_%s.png",     observable->data());
      else if ( idxSignal_or_Background == kBackground ) outputFileName = Form("makeVertexPlots_background_%s.png", observable->data());
      else assert(0);
      showHistograms(1150, 850,
		     histograms_rebinned["offlinePrimaryVertices"],        "offline",
		     histograms_rebinned["hltPhase2PixelVertices"],        "online",
		     histograms_rebinned["hltPhase2TrimmedPixelVertices"], "online(trimmed)",
		     0, "",
		     0, "",
		     0, "",
		     colors, lineStyles, 
		     0.050, 0.62, 0.74, 0.31, 0.18, 
		     labelTextLines, 0.050,
		     0.70, 0.62, 0.23, 0.06, 
		     xMin[*observable], xMax[*observable], xAxisTitles[*observable], 1.2, 
		     true, 1.e-4, 1.99, "Events", 1.4, 
		     outputFileName);
    }
  }

  delete inputFile_signal;
  //delete inputFile_background;
}

