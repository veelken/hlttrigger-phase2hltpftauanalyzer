
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

void makeResponsePlots()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  std::string inputFilePath = Form("%s/src/HLTTrigger/TallinnHLTPFTauAnalyzer/test/", gSystem->Getenv("CMSSW_BASE"));
  std::string inputFileName = "analyzePFTauResponse_signal_2020Jun18.root";
  TFile* inputFile = openFile(inputFilePath, inputFileName);

  std::vector<std::string> pfAlgos;
  pfAlgos.push_back("PFTau");
  pfAlgos.push_back("HpsPFTau");

  std::vector<std::string> vertexOptions;
  vertexOptions.push_back("WithOfflineVertices");
  vertexOptions.push_back("WithOnlineVertices");
  vertexOptions.push_back("WithOnlineVerticesTrimmed");
  
  std::vector<std::string> refTauTypes;
  refTauTypes.push_back("GenHadTaus");
  refTauTypes.push_back("OfflineTaus");

  std::vector<std::string> ptRanges;
  ptRanges.push_back("ptGt20");
  //ptRanges.push_back("ptGt25");
  ptRanges.push_back("ptGt30");
  //ptRanges.push_back("ptGt35");
  ptRanges.push_back("ptGt40");
  //ptRanges.push_back("ptGt45");
  ptRanges.push_back("ptGt50");

  std::vector<std::string> absEtaRanges;
  absEtaRanges.push_back("absEtaLt1p40");
  absEtaRanges.push_back("absEta1p40to2p17");
  absEtaRanges.push_back("absEta1p40to2p40");
  absEtaRanges.push_back("absEtaLt2p17");
  absEtaRanges.push_back("absEtaLt2p40");

  std::vector<std::string> decayModes;
  //decayModes.push_back("oneProng0Pi0");
  //decayModes.push_back("oneProng1Pi0");
  //decayModes.push_back("oneProng2Pi0");
  //decayModes.push_back("threeProng0Pi0");
  //decayModes.push_back("threeProng1Pi0");
  decayModes.push_back("all");

  std::vector<std::string> isolationWPs;
  //isolationWPs.push_back("relChargedIsoLt0p40");
  //isolationWPs.push_back("relChargedIsoLt0p20");
  //isolationWPs.push_back("relChargedIsoLt0p10");
  isolationWPs.push_back("relChargedIsoLt0p05");
  //isolationWPs.push_back("relChargedIsoLt0p02");
  //isolationWPs.push_back("relChargedIsoLt0p01");

  std::vector<std::string> observables;
  observables.push_back("response");

  std::map<std::string, int> rebin;   // key = observable
  rebin["response"]                      = 1;

  std::map<std::string, double> xMin; // key = observable
  xMin["response"]                       = 0.;

  std::map<std::string, double> xMax; // key = observable
  xMax["response"]                       = 2.;
  
  std::map<std::string, std::map<std::string, std::string>> xAxisTitles; // key = refTauType, observable
  xAxisTitles["GenHadTaus"]["response"]  = "p_{T}^{online} / p_{T}^{gen}";
  xAxisTitles["OfflineTaus"]["response"] = "p_{T}^{online} / p_{T}^{offline}";

  std::string dqmDirectory = "%sResponseAnalyzer%s";
  
  int colors[6] = { 1, 2, 8, 4, 6, 7 };
  int lineStyles[6] = { 1, 1, 1, 1, 1, 1 };
  int markerStyles[6] = { 22, 32, 20, 24, 21, 25 };

  // TallinnL1PFTaus 
  for ( std::vector<std::string>::const_iterator pfAlgo = pfAlgos.begin();
	pfAlgo != pfAlgos.end(); ++pfAlgo ) {
    for ( std::vector<std::string>::const_iterator vertexOption = vertexOptions.begin();
	  vertexOption != vertexOptions.end(); ++vertexOption ) {
      for ( std::vector<std::string>::const_iterator observable = observables.begin();
	    observable != observables.end(); ++observable ) {      
        for ( std::vector<std::string>::const_iterator ptRange = ptRanges.begin();
	      ptRange != ptRanges.end(); ++ptRange ) {
          for ( std::vector<std::string>::const_iterator absEtaRange = absEtaRanges.begin();
	        absEtaRange != absEtaRanges.end(); ++absEtaRange ) {
  	    for ( std::vector<std::string>::const_iterator isolationWP = isolationWPs.begin();
	          isolationWP != isolationWPs.end(); ++isolationWP ) {	    
	      std::map<std::string, std::map<std::string, TH1*>> histograms; // keys = refTauType, decayMode
	      for ( std::vector<std::string>::const_iterator refTauType = refTauTypes.begin();
	  	    refTauType != refTauTypes.end(); ++refTauType ) {
	        for ( std::vector<std::string>::const_iterator decayMode = decayModes.begin();
		      decayMode != decayModes.end(); ++decayMode ) {
	          std::string histogramName = Form((dqmDirectory + "_wrt%s/%s/%s/%s/%s_%s_%s_%s_%s").data(), 
                    pfAlgo->data(), vertexOption->data(), refTauType->data(), ptRange->data(), absEtaRange->data(), decayMode->data(), 
                    observable->data(), decayMode->data(), ptRange->data(), absEtaRange->data(), isolationWP->data());
		  TH1* histogram = loadHistogram(inputFile, histogramName);
		  histograms[*refTauType][*decayMode] = histogram;
	        }
	      }

	      for ( std::vector<std::string>::const_iterator refTauType = refTauTypes.begin();
		    refTauType != refTauTypes.end(); ++refTauType ) {
	        TH1* histogram = histograms[*refTauType]["all"]; 
	        std::vector<std::string> labelTextLines;
	        labelTextLines.push_back(Form("Mean = %1.3f", histogram->GetMean()));
	        labelTextLines.push_back(Form("RMS = %1.3f", histogram->GetRMS()));
	        std::string outputFileName = Form("makeResponsePlots_%s%s_wrt%s_%s_%s_%s_%s.png", 
                  pfAlgo->data(), vertexOption->data(), refTauType->data(), observable->data(), ptRange->data(), absEtaRange->data(), isolationWP->data());
	        showHistograms(1150, 850,
		  	       histogram, "",
			       0, "",
			       0, "",	     
			       0, "",
			       0, "",
			       0, "",
			       colors, lineStyles, 
			       0.050, 0.70, 0.74, 0.23, 0.18, 
			       labelTextLines, 0.050,
			       0.71, 0.76, 0.23, 0.15, 
			       xMin[*observable], xMax[*observable], xAxisTitles[*refTauType][*observable], 1.2, 
			       true, 1.e-4, 1.99, "Events", 1.4, 
			       outputFileName);
	      }
/*
	      for ( std::vector<std::string>::const_iterator refTauType = refTauTypes.begin();
		    refTauType != refTauTypes.end(); ++refTauType ) {
	        TH1* histogram_oneProng0Pi0   = histograms[*refTauType]["oneProng0Pi0"];
	        TH1* histogram_oneProng1Pi0   = histograms[*refTauType]["oneProng1Pi0"];
	        TH1* histogram_oneProng2Pi0   = histograms[*refTauType]["oneProng2Pi0"];
	        TH1* histogram_threeProng0Pi0 = histograms[*refTauType]["threeProng0Pi0"];
	        TH1* histogram_threeProng1Pi0 = histograms[*refTauType]["threeProng1Pi0"];
	        std::vector<std::string> labelTextLines;
	        std::string outputFileName = Form("makeResponsePlots_%s%s_wrt%s_vs_decayModes_%s_%s_%s_%s.png", 
                  pfAlgo->data(), vertexOption->data(), refTauType->data(), observable->data(), ptRange->data(), absEtaRange->data(), isolationWP->data());
	        showHistograms(1150, 850,
		  	       histogram_oneProng0Pi0,   "h^{#pm}",
			       histogram_oneProng1Pi0,   "h^{#pm}#pi^{0}",
			       histogram_oneProng2Pi0,   "h^{#pm}#pi^{0}#pi^{0}",
			       histogram_threeProng0Pi0, "h^{#pm}h^{#mp}h^{#pm}",
			       histogram_threeProng1Pi0, "h^{#pm}h^{#mp}h^{#pm}#pi^{0}",
			       0, "",
			       colors, lineStyles, 
			       0.050, 0.73, 0.67, 0.20, 0.26, 
			       labelTextLines, 0.050,
			       0.70, 0.62, 0.23, 0.06, 
			       xMin[*observable], xMax[*observable], xAxisTitles[*refTauType][*observable], 1.2, 
			       true, 1.e-4, 1.99, "Events", 1.4, 
			       outputFileName);
              }
 */
	    }
	  }
	}
      }
    }
  }

  delete inputFile;
}

