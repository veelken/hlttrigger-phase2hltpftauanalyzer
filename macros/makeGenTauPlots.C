
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
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

TH1* loadHistogram(TFile* inputFile, const std::string& histogramName, bool doNormalize)
{
  TH1* histogram = dynamic_cast<TH1*>(inputFile->Get(histogramName.data()));
  if ( !histogram ) {
    std::cerr << "Failed to load histogram = " << histogramName << " from file = " << inputFile->GetName() << " !!" << std::endl;
    assert(0);
  }
  if ( !histogram->GetSumw2N() ) histogram->Sumw2();
  if ( doNormalize ) {
    double integral = compIntegral(histogram);
    if ( integral > 0. ) {
      histogram->Scale(1./integral);
    }
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
                    int colors[], int markerStyles[], int lineStyles[], 
		    double topMargin, double leftMargin, double bottomMargin, double rightMargin,  
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
  
  canvas->SetTopMargin(topMargin);
  canvas->SetLeftMargin(leftMargin);
  canvas->SetBottomMargin(bottomMargin);
  canvas->SetRightMargin(rightMargin);

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

  //const std::string drawOption = "hist";
  //const std::string drawOption_legend = "l";
  const std::string drawOption = "e1p";
  const std::string drawOption_legend = "p";

  histogram1->SetMarkerColor(colors[0]);
  histogram1->SetMarkerSize(2);
  histogram1->SetMarkerStyle(markerStyles[0]);
  histogram1->SetLineColor(colors[0]);
  histogram1->SetLineWidth(2);
  histogram1->SetLineStyle(lineStyles[0]);
  histogram1->Draw(drawOption.data());

  if ( histogram2 ) {
    histogram2->SetMarkerColor(colors[1]);
    histogram2->SetMarkerSize(2);
    histogram2->SetMarkerStyle(markerStyles[1]);
    histogram2->SetLineColor(colors[1]);
    histogram2->SetLineWidth(2);
    histogram2->SetLineStyle(lineStyles[1]);
    histogram2->Draw(Form("%ssame", drawOption.data()));
  }

  if ( histogram3 ) {
    histogram3->SetMarkerColor(colors[2]);
    histogram3->SetMarkerSize(2);
    histogram3->SetMarkerStyle(markerStyles[2]);
    histogram3->SetLineColor(colors[2]);
    histogram3->SetLineWidth(2);
    histogram3->SetLineStyle(lineStyles[2]);
    histogram3->Draw(Form("%ssame", drawOption.data()));
  }

  if ( histogram4 ) {
    histogram4->SetMarkerColor(colors[3]);
    histogram4->SetMarkerSize(2);
    histogram4->SetMarkerStyle(markerStyles[3]);
    histogram4->SetLineColor(colors[3]);
    histogram4->SetLineWidth(2);
    histogram4->SetLineStyle(lineStyles[3]);
    histogram4->Draw(Form("%ssame", drawOption.data()));
  }

  if ( histogram5 ) {
    histogram5->SetMarkerColor(colors[4]);
    histogram5->SetMarkerSize(2);
    histogram5->SetMarkerStyle(markerStyles[4]);
    histogram5->SetLineColor(colors[4]);
    histogram5->SetLineWidth(2);
    histogram5->SetLineStyle(lineStyles[4]);
    histogram5->Draw(Form("%ssame", drawOption.data()));
  }

  if ( histogram6 ) {
    histogram6->SetMarkerColor(colors[5]);
    histogram6->SetMarkerSize(2);
    histogram6->SetMarkerStyle(markerStyles[5]);
    histogram6->SetLineColor(colors[5]);
    histogram6->SetLineWidth(2);
    histogram6->SetLineStyle(lineStyles[5]);
    histogram6->Draw(Form("%ssame", drawOption.data()));
  }

  TLegend* legend = 0;
  if ( legendEntry1 != "" ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, "", "brNDC"); 
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetTextSize(legendTextSize);
    legend->AddEntry(histogram1, legendEntry1.data(), drawOption_legend.data());
    if ( histogram2 ) legend->AddEntry(histogram2, legendEntry2.data(), drawOption_legend.data());
    if ( histogram3 ) legend->AddEntry(histogram3, legendEntry3.data(), drawOption_legend.data());
    if ( histogram4 ) legend->AddEntry(histogram4, legendEntry4.data(), drawOption_legend.data());
    if ( histogram5 ) legend->AddEntry(histogram5, legendEntry5.data(), drawOption_legend.data());
    if ( histogram6 ) legend->AddEntry(histogram6, legendEntry6.data(), drawOption_legend.data());
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

void makeGenTauPlots()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  std::string inputFilePath = Form("%s/src/L1Trigger/TallinnL1PFTauAnalyzer/test/", gSystem->Getenv("CMSSW_BASE"));
  std::string inputFileName = "GenHadTauAnalyzer_signal_qqH_2019Oct16.root";
  TFile* inputFile = openFile(inputFilePath, inputFileName);

  std::vector<std::string> observables;
  observables.push_back("pt");
  observables.push_back("eta");

  std::vector<std::string> tauTypes;
  tauTypes.push_back("genTau");
  tauTypes.push_back("genTau_offlineMatched");
  tauTypes.push_back("genHadTau");
  tauTypes.push_back("genHadTau_offlineMatched");
  //tauTypes.push_back("offlineHadTau");

  std::map<std::string, std::string> dqmDirectories; // key = tauType
  dqmDirectories["genTau"]                   = "DQMData/GenTauAnalyzer";
  dqmDirectories["genTau_offlineMatched"]    = "DQMData/GenTauAnalyzer_offlineMatched";
  dqmDirectories["genHadTau"]                = "DQMData/GenHadTauAnalyzer";
  dqmDirectories["genHadTau_offlineMatched"] = "DQMData/GenHadTauAnalyzer_offlineMatched";
  dqmDirectories["offlineHadTau"]            = "DQMData/PATTauAnalyzer";

  std::map<std::string, int> canvasSizesX; // key = tauType
  canvasSizesX["genTau"]                    = 1150;
  canvasSizesX["genTau_offlineMatched"]     = 1150;
  canvasSizesX["genHadTau"]                 = 1150;
  canvasSizesX["genHadTau_offlineMatched"]  = 1150;
  canvasSizesX["offlineHadTau"]             = 1150;

  std::map<std::string, int> canvasSizesY; // key = tauType
  canvasSizesY["genTau"]                    = 1050;
  canvasSizesY["genTau_offlineMatched"]     = 1050;
  canvasSizesY["genHadTau"]                 =  850;
  canvasSizesY["genHadTau_offlineMatched"]  =  850;
  canvasSizesY["offlineHadTau"]             =  850;

  std::map<std::string, double> leftMargins; // key = tauType
  leftMargins["genTau"]                     =    0.17;
  leftMargins["genTau_offlineMatched"]      =    0.17;
  leftMargins["genHadTau"]                  =    0.14;
  leftMargins["genHadTau_offlineMatched"]   =    0.14;
  leftMargins["offlineHadTau"]              =    0.14;

  std::map<std::string, double> bottomMargins; // key = tauType
  bottomMargins["genTau"]                   =    0.11;
  bottomMargins["genTau_offlineMatched"]    =    0.11;
  bottomMargins["genHadTau"]                =    0.14;
  bottomMargins["genHadTau_offlineMatched"] =    0.14;
  bottomMargins["offlineHadTau"]            =    0.14;

  std::map<std::string, int> rebin; // key = observable
  rebin["pt"]           =   1;
  rebin["eta"]          =   1;
  rebin["deltaR"]       =   1;
  
  std::map<std::string, double> xMin; // key = observable
  xMin["pt"]            =   0.;
  xMin["eta"]           =  -3.;
  xMin["deltaR"]        =   0.;

  std::map<std::string, double> xMax; // key = observable
  xMax["pt"]            = 200.;
  xMax["eta"]           =  +3.;
  xMax["deltaR"]        =   5.;
  
  std::map<std::string, std::string> xAxisTitles; // key = observable
  xAxisTitles["pt"]     = "p_{T} [GeV]";
  xAxisTitles["eta"]    = "#eta";
  xAxisTitles["deltaR"] = "#Delta R";

  std::map<std::string, double> xAxisOffsets; // key = tauType
  xAxisOffsets["genTau"]                   =    1.1;
  xAxisOffsets["genTau_offlineMatched"]    =    1.1;
  xAxisOffsets["genHadTau"]                =    1.2;
  xAxisOffsets["genHadTau_offlineMatched"] =    1.2;
  xAxisOffsets["offlineHadTau"]            =    1.2;

  std::map<std::string, double> yMin; // key = observable
  yMin["pt"]            =   0.;
  yMin["eta"]           =   0.;
  yMin["deltaR"]        =   0.;

  std::map<std::string, double> yMax; // key = observable
  yMax["pt"]            =   0.10;
  yMax["eta"]           =   0.08;
  yMax["deltaR"]        =   0.1;

  std::map<std::string, std::string> yAxisTitles; // key = observable
  yAxisTitles["pt"]     = "#frac{dN}{dp_{T}} [1/GeV]";
  yAxisTitles["eta"]    = "#frac{dN}{d#eta}";
  yAxisTitles["deltaR"] = "#frac{dN}{d#Delta R}";
  
  std::map<std::string, double> yAxisOffsets; // key = tauType
  yAxisOffsets["genTau"]                   =    1.7;
  yAxisOffsets["genTau_offlineMatched"]    =    1.7;
  yAxisOffsets["genHadTau"]                =    1.4;
  yAxisOffsets["genHadTau_offlineMatched"] =    1.4;
  yAxisOffsets["offlineHadTau"]            =    1.4;

  std::map<std::string, std::map<std::string, std::string>> legendEntries; // key = tauType, "leading" / "subleading"
  legendEntries["genTau"]["leading"]                      = "lead. #tau";
  legendEntries["genTau"]["subleading"]                   = "sublead. #tau";
  legendEntries["genTau_offlineMatched"]["leading"]       = "lead. #tau";
  legendEntries["genTau_offlineMatched"]["subleading"]    = "sublead. #tau";
  legendEntries["genHadTau"]["leading"]                   = "lead. #tau_{h}";
  legendEntries["genHadTau"]["subleading"]                = "sublead. #tau_{h}";
  legendEntries["genHadTau_offlineMatched"]["leading"]    = "lead. #tau_{h}";
  legendEntries["genHadTau_offlineMatched"]["subleading"] = "sublead. #tau_{h}";
  legendEntries["offlineHadTau"]["leading"]               = "lead. #tau_{h}";
  legendEntries["offlineHadTau"]["subleading"]            = "sublead. #tau_{h}";
  
  std::map<std::string, double> legendPosX; // key = tauType
  legendPosX["genTau"]                   =    0.68;
  legendPosX["genTau_offlineMatched"]    =    0.68;
  legendPosX["genHadTau"]                =    0.71;
  legendPosX["genHadTau_offlineMatched"] =    0.71;
  legendPosX["offlineHadTau"]            =    0.71;

  std::map<std::string, double> legendPosY; // key = tauType
  legendPosY["genTau"]                   =    0.74;
  legendPosY["genTau_offlineMatched"]    =    0.74;
  legendPosY["genHadTau"]                =    0.74;
  legendPosY["genHadTau_offlineMatched"] =    0.74;
  legendPosY["offlineHadTau"]            =    0.74;

  int colors[6] = { 1, 8, 2, 4, 6, 7 };
  int lineStyles[6] = { 1, 1, 1, 1, 1, 1 };
  int markerStyles[6] = { 26, 23, 20, 24, 21, 25 };

  std::vector<std::string> labelTextLines;

  for ( std::vector<std::string>::const_iterator tauType = tauTypes.begin();
	tauType != tauTypes.end(); ++tauType ) {
    for ( std::vector<std::string>::const_iterator observable = observables.begin();
	  observable != observables.end(); ++observable ) {
      std::string histogramName_leading = Form("%s/leadingTau_%s_all_absEtaLt2p40_pt20p00to1000p00", dqmDirectories[*tauType].data(), observable->data()); 
      TH1* histogram_leading = loadHistogram(inputFile, histogramName_leading, true);
      TH1* histogram_leading_rebinned = (TH1*)histogram_leading->Clone(Form("%s_rebinned", histogram_leading->GetName()));
      if ( rebin[*observable] > 1 ) histogram_leading_rebinned->Rebin(rebin[*observable]);

      std::string histogramName_subleading = Form("%s/subleadingTau_%s_all_absEtaLt2p40_pt20p00to1000p00", dqmDirectories[*tauType].data(), observable->data()); 
      TH1* histogram_subleading = loadHistogram(inputFile, histogramName_subleading, true);
      TH1* histogram_subleading_rebinned = (TH1*)histogram_subleading->Clone(Form("%s_rebinned", histogram_subleading->GetName()));
      if ( rebin[*observable] > 1 ) histogram_leading_rebinned->Rebin(rebin[*observable]);

      std::string outputFileName = Form("makeGenTauPlots_%s_%s.png", observable->data(), tauType->data());
      showHistograms(canvasSizesX[*tauType], canvasSizesY[*tauType], 
		     histogram_leading_rebinned,    legendEntries[*tauType]["leading"],
		     histogram_subleading_rebinned, legendEntries[*tauType]["subleading"],
		     0, "",
		     0, "",
		     0, "",
		     0, "",
		     colors, markerStyles, lineStyles, 
		     0.05, leftMargins[*tauType], bottomMargins[*tauType], 0.05,
		     0.050, legendPosX[*tauType], legendPosY[*tauType], 0.22, 0.18, 
		     labelTextLines, 0.050,
		     0.70, 0.62, 0.23, 0.06, 
		     xMin[*observable], xMax[*observable], xAxisTitles[*observable], xAxisOffsets[*tauType], 
		     false, yMin[*observable], yMax[*observable], yAxisTitles[*observable], yAxisOffsets[*tauType], 
		     outputFileName);

      delete histogram_leading_rebinned;
      delete histogram_subleading_rebinned;
    }

    std::string histogramName_deltaR = Form("%s/%s", dqmDirectories[*tauType].data(), "deltaR"); 
    TH1* histogram_deltaR = loadHistogram(inputFile, histogramName_deltaR, true);
    TH1* histogram_deltaR_rebinned = (TH1*)histogram_deltaR->Clone(Form("%s_rebinned", histogram_deltaR->GetName()));
    if ( rebin["deltaR"] > 1 ) histogram_deltaR_rebinned->Rebin(rebin["deltaR"]);

    std::string outputFileName_deltaR = Form("makeGenTauPlots_%s_%s.png", "deltaR", tauType->data());
    showHistograms(1150, 850,
		   histogram_deltaR_rebinned, "",
		   0, "",
		   0, "",
		   0, "",
		   0, "",
		   0, "",
		   colors, markerStyles, lineStyles, 
		   0.05, 0.14, 0.14, 0.05,
		   0.050, 0.71, 0.74, 0.22, 0.18, 
		   labelTextLines, 0.050,
		   0.70, 0.62, 0.23, 0.06, 
		   xMin["deltaR"], xMax["deltaR"], xAxisTitles["deltaR"], 1.2, 
		   false, yMin["deltaR"], yMax["deltaR"], yAxisTitles["deltaR"], 1.4, 
		   outputFileName_deltaR);
    
    delete histogram_deltaR_rebinned;
  }

  delete inputFile;
}

