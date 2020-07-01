
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
  std::cout << "<loadHistogram>: Loading histogram = " << histogramName << " from file = " << inputFile->GetName() << "." << std::endl;
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
  TAxis* xAxis = histogram->GetXaxis();
  int numBins = histogram->GetNbinsX();
  for ( int idxBin = 1; idxBin <= numBins; ++idxBin ) {
    double binContent = histogram->GetBinContent(idxBin);
    double binError = histogram->GetBinError(idxBin);
    double binEdgeLo =  xAxis->GetBinLowEdge(idxBin);
    double binEdgeHi =  xAxis->GetBinUpEdge(idxBin);
    double binWidth = binEdgeHi - binEdgeLo;
    std::cout << "bin #" << idxBin (x = " << binEdgeLo << ".." << binEdgeHi << "): " << binContent << " +/- " << binError << std::endl;
    assert(binWidth > 0.);
    histogram->SetBinContent(idxBin, binContent/binWidth);
    histogram->SetBinError(idxBin, binError/binWidth);
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

void makePFChargedCandPlots()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  std::string inputFilePath = Form("%s/src/HLTTrigger/TallinnHLTPFTauAnalyzer/test/", gSystem->Getenv("CMSSW_BASE"));
  std::string inputFileName = "analyzePFChargedCands_backgrounds_2020Jul01.root";
  TFile* inputFile = openFile(inputFilePath, inputFileName);

  std::vector<std::string> ptRanges;
  ptRanges.push_back("ptLt1");
  ptRanges.push_back("pt1to2");
  ptRanges.push_back("pt2to5");
  ptRanges.push_back("pt5to10");
  ptRanges.push_back("pt10to20");
  ptRanges.push_back("ptGt20");

  std::vector<std::string> absEtaRanges;
  absEtaRanges.push_back("absEtaLt1p40");
  //absEtaRanges.push_back("absEta1p40to2p17");
  absEtaRanges.push_back("absEta1p40to2p40");
  //absEtaRanges.push_back("absEtaLt2p17");
  absEtaRanges.push_back("absEtaLt2p40");

  std::vector<std::string> particleTypes;
  particleTypes.push_back("e");
  particleTypes.push_back("mu");
  particleTypes.push_back("h");
  particleTypes.push_back("other");

  std::map<std::string, std::string> legendEntries; // key = particleType
  legendEntries["e"]     = "e";
  legendEntries["mu"]    = "#mu";
  legendEntries["h"]     = "h^{#pm}";
  legendEntries[othere"] = "Other";

  std::vector<std::string> observables;
  observables.push_back("pfCandPt");
  observables.push_back("trackPt");
  observables.push_back("trackPt_div_pfCandPt");

  std::map<std::string, int> rebin;   // key = observable
  rebin["pfCandPt"]                   = 1;
  rebin["trackPt"]                    = 1;
  rebin["trackPt_div_pfCandPt"]       = 1;

  std::map<std::string, double> xMin; // key = observable
  xMin["pfCandPt"]                    =   0.;
  xMin["trackPt"]                     =   0.;
  xMin["trackPt_div_pfCandPt"]        =   0.;

  std::map<std::string, double> xMax; // key = observable
  xMax["pfCandPt"]                    = 100.;
  xMax["trackPt"]                     = 100.;
  xMax["trackPt_div_pfCandPt"]        =   2.;
  
  std::map<std::string, std::string> xAxisTitles; // key = observable
  xAxisTitles["pfCandPt"]             = "p_{T}^{pf} [GeV]";
  xAxisTitles["trackPt"]              = "p_{T}^{trk} [GeV]";
  xAxisTitles["trackPt_div_pfCandPt"] = "R = p_{T}^{trk} / p_{T}^{pf}";

  std::map<std::string, std::string> yAxisTitles; // key = observable
  yAxisTitles["pfCandPt"]             = "dN/dp_{T}^{pf} [GeV^{-1}]";
  yAxisTitles["trackPt"]              = "dN/dp_{T}^{trk} [GeV^{-1}]";
  yAxisTitles["trackPt_div_pfCandPt"] = "dN/dR";

  std::string dqmDirectory = "recoPFChargedCandAnalyzer";
  
  int colors[6]       = {  1,  2,  8,  4,  6,  7 };
  int lineStyles[6]   = {  1,  1,  1,  1,  1,  1 };
  int markerStyles[6] = { 22, 32, 20, 24, 21, 25 };

  for ( std::vector<std::string>::const_iterator observable = observables.begin();
	observable != observables.end(); ++observable ) {      
    for ( std::vector<std::string>::const_iterator ptRange = ptRanges.begin();
          ptRange != ptRanges.end(); ++ptRange ) {
      for ( std::vector<std::string>::const_iterator absEtaRange = absEtaRanges.begin();
            absEtaRange != absEtaRanges.end(); ++absEtaRange ) {
        for ( std::vector<std::string>::const_iterator particleType = particleTypes.begin();
	      particleType != particleTypes.end(); ++particleType ) {	    
	  std::map<std::string, TH1*> histograms; // key = particleType
	  std::string histogramName = Form("%s/%s_%s_%s_%s", dqmDirectory.data()
            observable->data(), particleType->data(), ptRange->data(), absEtaRange->data());
          TH1* histogram = loadHistogram(inputFile, histogramName);
          histograms[*particleType] = histogram;
	}
        for ( std::vector<std::string>::const_iterator particleType = particleTypes.begin();
	      particleType != particleTypes.end(); ++particleType ) {	  
          const TH1* histogram = histograms[*particleType];
	  std::vector<std::string> labelTextLines1;
	  labelTextLines1.push_back(Form("Mean = %1.3f", histogram->GetMean()));
	  labelTextLines1.push_back(Form("RMS = %1.3f", histogram->GetRMS()));
	  std::string outputFileName1 = Form("makePFChargedCandPlots_%s_%s_%s_%s.png", 
            observable->data(), particleType->data(), ptRange->data(), absEtaRange->data());
	  showHistograms(1150, 850,
		  	 histogram, "",
			 nullptr, "",
		         nullptr, "",	     
			 nullptr, "",
			 nullptr, "",
			 nullptr, "",
                         colors, lineStyles, 
			 0.050, 0.70, 0.74, 0.23, 0.18, 
			 labelTextLines1, 0.050,
			 0.71, 0.76, 0.23, 0.15, 
			 xMin[*observable], xMax[*observable], xAxisTitles[*observable], 1.2, 
			 true, 1.e-4, 1.99, yAxisTitles[*observable], 1.4, 
			 outputFileName1);
	  }
        }
        const TH1* histogram_e     = histograms["e"];
        const TH1* histogram_mu    = histograms["mu"];
        const TH1* histogram_h     = histograms["h"];
        const TH1* histogram_other = histograms["other"];
        std::vector<std::string> labelTextLines2;
	std::string outputFileName2 = Form("makePFChargedCandPlots_%s_all_%s_%s.png", 
          observable->data(), ptRange->data(), absEtaRange->data());
	showHistograms(1150, 850,
                       histogram_e,     legendEntry["e"],
                       histogram_mu,    legendEntry["mu"],
                       histogram_h,     legendEntry["h"],
                       histogram_other, legendEntry["other"],
                       nullptr, "",
		       nullptr, "",
                       colors, lineStyles,
                       0.050, 0.70, 0.71, 0.23, 0.21, 
                       labelTextLines2, 0.050,
                       0.71, 0.76, 0.23, 0.15, 
                       xMin[*observable], xMax[*observable], xAxisTitles[*observable], 1.2, 
                       true, 1.e-4, 1.99, yAxisTitles[*observable], 1.4, 
                       outputFileName2);
      }
    }
  }

  delete inputFile;
}

