
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
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

void showHistograms(double canvasSizeX, double canvasSizeY,
                    TH1* histogram_minbias, const std::string& legendEntry_minbias,
                    TH1* histogram_qcd1, const std::string& legendEntry_qcd1,
                    TH1* histogram_qcd2, const std::string& legendEntry_qcd2,
                    TH1* histogram_qcd3, const std::string& legendEntry_qcd3,
                    TH1* histogram_qcd4, const std::string& legendEntry_qcd4,
                    TH1* histogram_qcd5, const std::string& legendEntry_qcd5,
                    TH1* histogram_qcd6, const std::string& legendEntry_qcd6,
                    int colors[], int fillStyles[],
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
  canvas->SetBottomMargin(0.12);
  canvas->SetRightMargin(0.05);

  canvas->SetLogy(useLogScale);
  
  canvas->SetGridx(1);
  canvas->SetGridy(1);

  if ( !histogram_minbias ) {
    std::cerr << "<showHistograms>: histogram_minbias = NULL --> skipping !!" << std::endl;
    return;
  }

  THStack* histogramStack_qcd = new THStack("histogramStack", "histogramStack");

  histogramStack_qcd->SetTitle("");
  histogramStack_qcd->SetMinimum(yMin);
  histogramStack_qcd->SetMaximum(yMax);

  if ( histogram_qcd6 ) {
    histogram_qcd6->SetStats(false);
    histogram_qcd6->SetFillColor(colors[5]);
    histogram_qcd6->SetFillStyle(fillStyles[5]);
    histogram_qcd6->SetLineColor(colors[5]);
    histogram_qcd6->SetLineWidth(2);
    histogram_qcd6->SetLineStyle(1);
    histogramStack_qcd->Add(histogram_qcd6);
  }

  if ( histogram_qcd5 ) {
    histogram_qcd5->SetStats(false);
    histogram_qcd5->SetFillColor(colors[4]);
    histogram_qcd5->SetFillStyle(fillStyles[4]);
    histogram_qcd5->SetLineColor(colors[4]);
    histogram_qcd5->SetLineWidth(2);
    histogram_qcd5->SetLineStyle(1);
    histogramStack_qcd->Add(histogram_qcd5);
  }

  if ( histogram_qcd4 ) {
    histogram_qcd4->SetStats(false);
    histogram_qcd4->SetFillColor(colors[3]);
    histogram_qcd4->SetFillStyle(fillStyles[3]);
    histogram_qcd4->SetLineColor(colors[3]);
    histogram_qcd4->SetLineWidth(2);
    histogram_qcd4->SetLineStyle(1);
    histogramStack_qcd->Add(histogram_qcd4);
  }

  if ( histogram_qcd3 ) {
    histogram_qcd3->SetStats(false);
    histogram_qcd3->SetFillColor(colors[2]);
    histogram_qcd3->SetFillStyle(fillStyles[2]);
    histogram_qcd3->SetLineColor(colors[2]);
    histogram_qcd3->SetLineWidth(2);
    histogram_qcd3->SetLineStyle(1);
    histogramStack_qcd->Add(histogram_qcd3);
  }

  if ( histogram_qcd2 ) {
    histogram_qcd2->SetStats(false);
    histogram_qcd2->SetFillColor(colors[1]);
    histogram_qcd2->SetFillStyle(fillStyles[1]);
    histogram_qcd2->SetLineColor(colors[1]);
    histogram_qcd2->SetLineWidth(2);
    histogram_qcd2->SetLineStyle(1);
    histogramStack_qcd->Add(histogram_qcd2);
  }

  if ( histogram_qcd1 ) {
    histogram_qcd1->SetStats(false);
    histogram_qcd1->SetFillColor(colors[0]);
    histogram_qcd1->SetFillStyle(fillStyles[0]);
    histogram_qcd1->SetLineColor(colors[0]);
    histogram_qcd1->SetLineWidth(2);
    histogram_qcd1->SetLineStyle(1);
    histogramStack_qcd->Add(histogram_qcd1);
  }

  histogramStack_qcd->Draw("hist");

  std::cout << "histogramStack_qcd->GetHistogram() = " << histogramStack_qcd->GetHistogram() << std::endl;

  TAxis* xAxis = histogramStack_qcd->GetHistogram()->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleSize(0.045);
  xAxis->SetTitleOffset(xAxisOffset);  
  if ( xMax > xMin ) {
    std::cout << "limiting x-axis range to " << xMin << ".." << xMax << std::endl;
    xAxis->SetRangeUser(xMin, xMax);
  }

  TAxis* yAxis = histogramStack_qcd->GetHistogram()->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleSize(0.045);
  yAxis->SetTitleOffset(yAxisOffset);

  histogram_minbias->SetTitle("");
  histogram_minbias->SetStats(false);
  histogram_minbias->SetLineColor(1);
  histogram_minbias->SetLineWidth(3);
  histogram_minbias->Draw("histsame");

  TLegend* legend = 0;
  if ( legendEntry_minbias != "" ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, "", "brNDC"); 
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetTextSize(legendTextSize);
    legend->AddEntry(histogram_minbias, legendEntry_minbias.data(), "l");
    if ( histogram_qcd1 ) legend->AddEntry(histogram_qcd1, legendEntry_qcd1.data(), "f");
    if ( histogram_qcd2 ) legend->AddEntry(histogram_qcd2, legendEntry_qcd2.data(), "f");
    if ( histogram_qcd3 ) legend->AddEntry(histogram_qcd3, legendEntry_qcd3.data(), "f");
    if ( histogram_qcd4 ) legend->AddEntry(histogram_qcd4, legendEntry_qcd4.data(), "f");
    if ( histogram_qcd5 ) legend->AddEntry(histogram_qcd5, legendEntry_qcd5.data(), "f");
    if ( histogram_qcd6 ) legend->AddEntry(histogram_qcd6, legendEntry_qcd6.data(), "f");
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

  histogram_minbias->Draw("axissame");

  canvas->Update();
  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  
  delete histogramStack_qcd;
  delete label;
  delete legend;
  delete canvas;  
}

void makeGenJetPlots()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  std::string inputFilePath = "/hdfs/local/veelken/Phase2HLT/rate/2020Aug04v3/";

  std::vector<std::string> processes;
  processes.push_back("minbias");
  processes.push_back("QCD");

  std::map<std::string, std::string> inputFileNames; // key = process
  inputFileNames["minbias"] = "hadd_minbias_all.root";
  inputFileNames["QCD"]     = "hadd_QCD_all.root";

  std::map<std::string, TFile*> inputFiles; // key = process
  for ( std::vector<std::string>::const_iterator process = processes.begin();
	process != processes.end(); ++process ) {
    std::string inputFileName_full = inputFilePath;
    if ( inputFileName_full.find_last_of("/") != (inputFileName_full.size() - 1) ) inputFileName_full.append("/");
    inputFileName_full.append(inputFileNames[*process]);
    TFile* inputFile = new TFile(inputFileName_full.data());
    if ( !inputFile ) {
      std::cerr << "Failed to open input file = '" << inputFileName_full << "' !!" << std::endl;
      assert(0);
    }
    inputFiles[*process] = inputFile;
  }

  std::vector<std::string> ptRanges;
  ptRanges.push_back("ptGt20");
  ptRanges.push_back("ptGt25");
  ptRanges.push_back("ptGt30");

  std::vector<std::string> absEtaRanges;
  absEtaRanges.push_back("absEtaLt2p40");

  std::vector<std::string> observables;
  observables.push_back("leadJet_pt");
  observables.push_back("subleadJet_pt");
  observables.push_back("numJets");
  observables.push_back("HT");
  
  std::map<std::string, std::string> xAxisTitles; // key = observable
  xAxisTitles["leadJet_pt"]    = "p_{T} [GeV]";
  xAxisTitles["subleadJet_pt"] = "p_{T} [GeV]";
  xAxisTitles["numJets"]       = "N_{jet}";

  std::vector<std::string> ptHatRanges = { "qcd1", "qcd2", "qcd3", "qcd4", "qcd5", "qcd6" };

  std::map<std::string, std::string> sampleNames_qcd;   // key = ptHatRange
  sampleNames_qcd["qcd1"] = "qcd_pt30to50";
  sampleNames_qcd["qcd2"] = "qcd_pt50to80";
  sampleNames_qcd["qcd3"] = "qcd_pt80to120";
  sampleNames_qcd["qcd4"] = "qcd_pt120to170";
  sampleNames_qcd["qcd5"] = "qcd_pt170to300";
  sampleNames_qcd["qcd6"] = "qcd_ptGt300";

  std::map<std::string, std::string> legendEntries_qcd; // key = ptHatRange
  legendEntries_qcd["qcd1"] = " 30 < #hat{p}_{T} <  50";
  legendEntries_qcd["qcd2"] = " 50 < #hat{p}_{T} <  80";
  legendEntries_qcd["qcd3"] = " 80 < #hat{p}_{T} < 120";
  legendEntries_qcd["qcd4"] = "120 < #hat{p}_{T} < 170";
  legendEntries_qcd["qcd5"] = "170 < #hat{p}_{T} < 300";
  legendEntries_qcd["qcd6"] = "300 < #hat{p}_{T} < 470";

  std::string dqmDirectory = "GenJetAnalyzer";

  std::vector<std::string> labelTextLines;

  int colors[6]     = {    7,    6,    4,   28,    8,    2 };
  int fillStyles[6] = { 1001, 1001, 1001, 1001, 1001, 1001 };

  for ( std::vector<std::string>::const_iterator ptRange = ptRanges.begin();
	ptRange != ptRanges.end(); ++ptRange ) {
    for ( std::vector<std::string>::const_iterator absEtaRange = absEtaRanges.begin();
	  absEtaRange != absEtaRanges.end(); ++absEtaRange ) {
      for ( std::vector<std::string>::const_iterator observable = observables.begin();
	    observable != observables.end(); ++observable ) {
        std::string histogramName_minbias = Form(("%s/" + dqmDirectory + "/%s/%s/%s__%s_%s").data(), 
          "minbias", ptRange->data(), absEtaRange->data(), 
          observable->data(), ptRange->data(), absEtaRange->data());
        TH1* histogram_minbias = loadHistogram(inputFiles["minbias"], histogramName_minbias);
        std::cout << "minbias: " << histogram_minbias->Integral() << std::endl;
        std::map<std::string, TH1*> histograms_qcd; // key = ptHatRange
        for ( std::vector<std::string>::const_iterator ptHatRange = ptHatRanges.begin();
              ptHatRange != ptHatRanges.end(); ++ptHatRange ) {
          std::string histogramName_qcd = Form(("%s/" + dqmDirectory + "/%s/%s/%s__%s_%s").data(), 
            sampleNames_qcd[*ptHatRange].data(), ptRange->data(), absEtaRange->data(), 
            observable->data(), ptRange->data(), absEtaRange->data());
          TH1* histogram_qcd = loadHistogram(inputFiles["QCD"], histogramName_qcd);
          histogram_qcd->Scale(1./200);
          std::cout << "QCD(" << legendEntries_qcd[*ptHatRange] << "): " << histogram_qcd->Integral() << std::endl;
          histograms_qcd[*ptHatRange] = histogram_qcd;
        }
        std::string outputFileName = Form("makeGenJetPlots_%s_for_%s_and_%s.png", 
          observable->data(), ptRange->data(), absEtaRange->data());
        showHistograms(1150, 1150,
                       histogram_minbias,      "Minbias",
                       histograms_qcd["qcd6"], legendEntries_qcd["qcd6"],
                       histograms_qcd["qcd5"], legendEntries_qcd["qcd5"],
                       histograms_qcd["qcd4"], legendEntries_qcd["qcd4"],
                       histograms_qcd["qcd3"], legendEntries_qcd["qcd3"],
                       histograms_qcd["qcd2"], legendEntries_qcd["qcd2"],
                       histograms_qcd["qcd1"], legendEntries_qcd["qcd1"],
		       colors, fillStyles, 
		       0.040, 0.63, 0.63, 0.25, 0.31,
		       labelTextLines, 0.050,
		       0.63, 0.66, 0.26, 0.07, 
		       -1., -1., xAxisTitles[*observable], 1.2, 
		       true, 1.e0, 1.e+6, "Rate", 1.4, 
		       outputFileName);
      }
    }
  }

  for ( std::map<std::string, TFile*>::const_iterator inputFile = inputFiles.begin();
        inputFile != inputFiles.end(); ++inputFile ) {
    delete inputFile->second;
  }
}

