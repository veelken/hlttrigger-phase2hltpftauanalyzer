
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
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

TH2* loadHistogram2d(TFile* inputFile, const std::string& histogramName)
{
  TH2* histogram = dynamic_cast<TH2*>(inputFile->Get(histogramName.data()));
  if ( !histogram ) {
    std::cerr << "Failed to load histogram = " << histogramName << " from file = " << inputFile->GetName() << " !!" << std::endl;
    assert(0);
  }
  if ( !histogram->GetSumw2N() ) histogram->Sumw2();
  return histogram;
}

double compRate(TH2* histogram2d, int idxBinX, int minNumTaus)
{
  const TAxis* yAxis = histogram2d->GetYaxis();
  int numBinsY = yAxis->GetNbins();
  double integral = 0.;
  double integral_failed = 0.;
  for ( int idxBinY = 1; idxBinY <= numBinsY; ++idxBinY ) { 
    double binContent = histogram2d->GetBinContent(idxBinX, idxBinY);
    integral += binContent;
    double numTaus = yAxis->GetBinCenter(idxBinY);
    if ( numTaus < minNumTaus ) integral_failed += binContent;
  }
  double integral_passed = integral - integral_failed;
  double max_rate = 2.8e+7; // bunch-crossing frequency of colliding bunches = 28 MHz
  //rate = max_rate*(integral_passed/integral);
  double rate = integral_passed;
  if ( rate > max_rate ) rate = max_rate;
  return rate;
}

TH1* makeRateHistogram(TH2* histogram2d, int minNumTaus)
{
  std::string histogramName_rate = Form("%s_rate_minNumTausEq%i", histogram2d->GetName(), minNumTaus);
  TH1* histogram_rate = histogram2d->ProjectionX(histogramName_rate.data());
  histogram_rate->Reset();
  const TAxis* xAxis = histogram2d->GetXaxis();
  int numBinsX = xAxis->GetNbins();
  for ( int idxBinX = 1; idxBinX <= numBinsX; ++idxBinX ) { 
    double rate = compRate(histogram2d, idxBinX, minNumTaus);
    histogram_rate->SetBinContent(idxBinX, rate);
    histogram_rate->SetBinError(idxBinX, 0.); // CV: error computation not implemented (needed) yet
  }
  return histogram_rate;
}

TH1* sumHistograms(const std::map<std::string, TH1*>& histograms, const std::vector<std::string>& processes)
{
  TH1* histogram_sum = nullptr;
  for ( std::vector<std::string>::const_iterator process = processes.begin();
        process != processes.end(); ++process ) {
    if ( histograms.find(*process) == histograms.end() ) {
      std::cerr << "Error in <sumHistograms>: Failed to find histogram for process = '" << (*process) << "' !!" << std::endl;
      assert(0);
    }
    const TH1* histogram = histograms.find(*process)->second;
    if ( histogram_sum ) {
      histogram_sum->Add(histogram);
    } else {
      std::string histogramName_sum = Form("%s_sum", histogram->GetName());
      histogram_sum = (TH1*)histogram->Clone(histogramName_sum.data());
      if ( !histogram_sum->GetSumw2N() ) histogram_sum->Sumw2();
    }
  }
  if ( !histogram_sum ) {
    std::cerr << "Error in <sumHistograms>: Failed to sum histograms (#processes = " << processes.size() << ") !!" << std::endl;
    assert(0);
  }
  return histogram_sum;
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
  canvas->SetBottomMargin(0.12);
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

TH1* getHistogram(std::map<std::string, TH1*>& histograms, const std::string& isolationWP)
{
  TH1* histogram = nullptr;
  if ( isolationWP != "" && histograms.find(isolationWP) != histograms.end() ) 
  {
    histogram = histograms.find(isolationWP)->second;
  }
  return histogram;
}

std::string getLegendEntry(std::map<std::string, std::string>& legendEntries, const std::string& isolationWP)
{
  std::string legendEntry = "";
  if ( isolationWP != "" && legendEntries.find(isolationWP) != legendEntries.end() ) 
  {
    legendEntry = legendEntries.find(isolationWP)->second;
  }
  return legendEntry;
}

void makeRatePlots()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  std::string inputFilePath = "/hdfs/local/veelken/Phase2HLT/rate/2020Jul29/";

  std::vector<std::string> processes;
  processes.push_back("minbias");
  processes.push_back("QCD");
  processes.push_back("DY");
  processes.push_back("W");

  std::map<std::string, std::string> inputFileNames; // key = process
  inputFileNames["minbias"] = "hadd_minbias_all.root";
  inputFileNames["QCD"]     = "hadd_QCD_all.root";
  inputFileNames["DY"]      = "hadd_DY_all.root";
  inputFileNames["W"]       = "hadd_W_all.root";

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

  std::vector<std::string> pfAlgos;
  //pfAlgos.push_back("PFTau");
  pfAlgos.push_back("HpsPFTau");

  std::vector<std::string> vertexOptions;
  vertexOptions.push_back("8HitsMaxDeltaZWithOfflineVertices");
  vertexOptions.push_back("8HitsMaxDeltaZToLeadTrackWithOfflineVertices");
  vertexOptions.push_back("8HitsMaxDeltaZWithOnlineVertices");
  vertexOptions.push_back("8HitsMaxDeltaZToLeadTrackWithOnlineVertices");
  //vertexOptions.push_back("8HitsMaxDeltaZWithOnlineVerticesTrimmed");
  //vertexOptions.push_back("8HitsMaxDeltaZToLeadTrackWithOnlineVerticesTrimmed");
/*
  vertexOptions.push_back("5HitsMaxDeltaZWithOfflineVertices");
  vertexOptions.push_back("5HitsMaxDeltaZToLeadTrackWithOfflineVertices");
  vertexOptions.push_back("5HitsMaxDeltaZWithOnlineVertices");
  vertexOptions.push_back("5HitsMaxDeltaZToLeadTrackWithOnlineVertices");
  //vertexOptions.push_back("5HitsMaxDeltaZWithOnlineVerticesTrimmed");
  //vertexOptions.push_back("5HitsMaxDeltaZToLeadTrackWithOnlineVerticesTrimmed");
  vertexOptions.push_back("3HitsMaxDeltaZWithOfflineVertices");
  vertexOptions.push_back("3HitsMaxDeltaZToLeadTrackWithOfflineVertices");
  vertexOptions.push_back("3HitsMaxDeltaZWithOnlineVertices");
  vertexOptions.push_back("3HitsMaxDeltaZToLeadTrackWithOnlineVertices");
  //vertexOptions.push_back("3HitsMaxDeltaZWithOnlineVerticesTrimmed");
  //vertexOptions.push_back("3HitsMaxDeltaZToLeadTrackWithOnlineVerticesTrimmed");
 */
  std::map<std::string, std::string> srcVertices; // key = vertexOption
  srcVertices["8HitsMaxDeltaZWithOfflineVertices"]                  = "offlinePrimaryVertices";
  srcVertices["8HitsMaxDeltaZToLeadTrackWithOfflineVertices"]       = "offlinePrimaryVertices";
  srcVertices["8HitsMaxDeltaZWithOnlineVertices"]                   = "hltPhase2PixelVertices";
  srcVertices["8HitsMaxDeltaZToLeadTrackWithOnlineVertices"]        = "hltPhase2PixelVertices";
  srcVertices["8HitsMaxDeltaZWithOnlineVerticesTrimmed"]            = "hltPhase2TrimmedPixelVertices";
  srcVertices["8HitsMaxDeltaZToLeadTrackWithOnlineVerticesTrimmed"] = "hltPhase2TrimmedPixelVertices";
  srcVertices["5HitsMaxDeltaZWithOfflineVertices"]                  = "offlinePrimaryVertices";
  srcVertices["5HitsMaxDeltaZToLeadTrackWithOfflineVertices"]       = "offlinePrimaryVertices";
  srcVertices["5HitsMaxDeltaZWithOnlineVertices"]                   = "hltPhase2PixelVertices";
  srcVertices["5HitsMaxDeltaZToLeadTrackWithOnlineVertices"]        = "hltPhase2PixelVertices";
  srcVertices["5HitsMaxDeltaZWithOnlineVerticesTrimmed"]            = "hltPhase2TrimmedPixelVertices";
  srcVertices["5HitsMaxDeltaZToLeadTrackWithOnlineVerticesTrimmed"] = "hltPhase2TrimmedPixelVertices";
  srcVertices["3HitsMaxDeltaZWithOfflineVertices"]                  = "offlinePrimaryVertices";
  srcVertices["3HitsMaxDeltaZToLeadTrackWithOfflineVertices"]       = "offlinePrimaryVertices";
  srcVertices["3HitsMaxDeltaZWithOnlineVertices"]                   = "hltPhase2PixelVertices";
  srcVertices["3HitsMaxDeltaZToLeadTrackWithOnlineVertices"]        = "hltPhase2PixelVertices";
  srcVertices["3HitsMaxDeltaZWithOnlineVerticesTrimmed"]            = "hltPhase2TrimmedPixelVertices";
  srcVertices["3HitsMaxDeltaZToLeadTrackWithOnlineVerticesTrimmed"] = "hltPhase2TrimmedPixelVertices";

  std::vector<std::string> tauIdOptions;
  tauIdOptions.push_back("sumChargedIso");
  tauIdOptions.push_back("deepTau");

  std::vector<std::string> l1MatchingOptions;
  l1MatchingOptions.push_back("");            // CV: no matching of HLT taus to L1 taus
  l1MatchingOptions.push_back("MatchedToL1"); 

  std::vector<std::string> absEtaRanges;
  absEtaRanges.push_back("absEtaLt1p40");
  //absEtaRanges.push_back("absEta1p40to2p17");
  absEtaRanges.push_back("absEta1p40to2p40");
  //absEtaRanges.push_back("absEtaLt2p17");
  absEtaRanges.push_back("absEtaLt2p40");

  std::vector<std::string> min_leadTrackPtValues;
  min_leadTrackPtValues.push_back("leadTrackPtGt1");
  min_leadTrackPtValues.push_back("leadTrackPtGt2");
  min_leadTrackPtValues.push_back("leadTrackPtGt5");

  std::map<std::string, std::vector<std::string>> isolationWPs; // key = tauIdOption
  isolationWPs["sumChargedIso"].push_back("noIsolation");
  isolationWPs["sumChargedIso"].push_back("relDiscriminatorLt0p400");
  isolationWPs["sumChargedIso"].push_back("relDiscriminatorLt0p200");
  isolationWPs["sumChargedIso"].push_back("relDiscriminatorLt0p100");
  isolationWPs["sumChargedIso"].push_back("relDiscriminatorLt0p050");
  isolationWPs["deepTau"].push_back("absDiscriminatorGt0p260");
  isolationWPs["deepTau"].push_back("absDiscriminatorGt0p425");
  isolationWPs["deepTau"].push_back("absDiscriminatorGt0p598");
  isolationWPs["deepTau"].push_back("absDiscriminatorGt0p785");
  isolationWPs["deepTau"].push_back("absDiscriminatorGt0p883");
  isolationWPs["deepTau"].push_back("absDiscriminatorGt0p931");
  
  std::map<std::string, std::map<std::string, std::string>> legendEntries_vs_isolationWPs; // key = tauIdOption, isolationWP
  legendEntries_vs_isolationWPs["sumChargedIso"]["noIsolation"]             = "No Isolation";
  legendEntries_vs_isolationWPs["sumChargedIso"]["relDiscriminatorLt0p400"] = "I_{ch} < 0.40*p_{T}";
  legendEntries_vs_isolationWPs["sumChargedIso"]["relDiscriminatorLt0p200"] = "I_{ch} < 0.20*p_{T}";
  legendEntries_vs_isolationWPs["sumChargedIso"]["relDiscriminatorLt0p100"] = "I_{ch} < 0.10*p_{T}";
  legendEntries_vs_isolationWPs["sumChargedIso"]["relDiscriminatorLt0p050"] = "I_{ch} < 0.05*p_{T}";
  legendEntries_vs_isolationWPs["deepTau"]["absDiscriminatorGt0p260"]       = "D > 0.260";
  legendEntries_vs_isolationWPs["deepTau"]["absDiscriminatorGt0p425"]       = "D > 0.425";
  legendEntries_vs_isolationWPs["deepTau"]["absDiscriminatorGt0p598"]       = "D > 0.598";
  legendEntries_vs_isolationWPs["deepTau"]["absDiscriminatorGt0p785"]       = "D > 0.785";
  legendEntries_vs_isolationWPs["deepTau"]["absDiscriminatorGt0p883"]       = "D > 0.883";
  legendEntries_vs_isolationWPs["deepTau"]["absDiscriminatorGt0p931"]       = "D > 0.931";

  std::map<std::string, std::string> legendEntries_vs_leadTrackPt; // key = min_leadTrackPt
  legendEntries_vs_leadTrackPt["leadTrackPtGt1"] = "lead. Track p_{T} > 1 GeV";
  legendEntries_vs_leadTrackPt["leadTrackPtGt2"] = "lead. Track p_{T} > 2 GeV";
  legendEntries_vs_leadTrackPt["leadTrackPtGt5"] = "lead. Track p_{T} > 5 GeV";

  std::map<std::string, std::string> legendEntries_vs_l1MatchingOption; // key = l1MatchingOption
  legendEntries_vs_l1MatchingOption[""]            = "wo. L1 Matching";
  legendEntries_vs_l1MatchingOption["MatchedToL1"] = "w. L1 Matching";

  std::map<std::string, std::string> legendEntries_vs_processes; // key = process
  legendEntries_vs_processes["minbias"] = "Minimum Bias";
  legendEntries_vs_processes["QCD"]     = "QCD";
  legendEntries_vs_processes["DY"]      = "Drell-Yan";
  legendEntries_vs_processes["W"]       = "W+jets";

  std::string dqmDirectory = "%sAnalyzerBackground%s%s_%s";

  std::vector<std::string> labelTextLines;

  int colors[6]     = { 1, 2, 8, 4, 6, 7 };
  int lineStyles[6] = { 1, 1, 1, 1, 1, 1 };

  typedef std::map<std::string, TH1*>              string_to_TH1Map1;
  typedef std::map<std::string, string_to_TH1Map1> string_to_TH1Map2;
  typedef std::map<std::string, string_to_TH1Map2> string_to_TH1Map3;
  typedef std::map<std::string, string_to_TH1Map3> string_to_TH1Map4;
  typedef std::map<std::string, string_to_TH1Map4> string_to_TH1Map5;
  typedef std::map<std::string, string_to_TH1Map5> string_to_TH1Map6;
  typedef std::map<std::string, string_to_TH1Map6> string_to_TH1Map7;
  typedef std::map<std::string, string_to_TH1Map7> string_to_TH1Map8;
  string_to_TH1Map8 histograms_rateSingleTau_vs_processes;        // key = pfAlgo, vertexOption, tauIdOption, l1MatchingOption, absEtaRange, min_leadTrackPt, isolationWP, process
  string_to_TH1Map8 histograms_rateDoubleTau_vs_processes;        // key = pfAlgo, vertexOption, tauIdOption, l1MatchingOption, absEtaRange, min_leadTrackPt, isolationWP, process
  string_to_TH1Map7 histograms_rateSingleTau_vs_isolationWPs;     // key = pfAlgo, vertexOption, tauIdOption, l1MatchingOption, absEtaRange, min_leadTrackPt, isolationWP
  string_to_TH1Map7 histograms_rateDoubleTau_vs_isolationWPs;     // key = pfAlgo, vertexOption, tauIdOption, l1MatchingOption, absEtaRange, min_leadTrackPt, isolationWP
  string_to_TH1Map7 histograms_rateSingleTau_vs_leadTrackPt;      // key = pfAlgo, vertexOption, tauIdOption, l1MatchingOption, absEtaRange, isolationWP, min_leadTrackPt
  string_to_TH1Map7 histograms_rateDoubleTau_vs_leadTrackPt;      // key = pfAlgo, vertexOption, tauIdOption, l1MatchingOption, absEtaRange, isolationWP, min_leadTrackPt
  string_to_TH1Map7 histograms_rateSingleTau_vs_l1MatchingOption; // key = pfAlgo, vertexOption, tauIdOption, absEtaRange, min_leadTrackPt, isolationWP, l1MatchingOption
  string_to_TH1Map7 histograms_rateDoubleTau_vs_l1MatchingOption; // key = pfAlgo, vertexOption, tauIdOption, absEtaRange, min_leadTrackPt, isolationWP, l1MatchingOption
  
  for ( std::vector<std::string>::const_iterator pfAlgo = pfAlgos.begin();
	pfAlgo != pfAlgos.end(); ++pfAlgo ) {
    for ( std::vector<std::string>::const_iterator vertexOption = vertexOptions.begin();
	  vertexOption != vertexOptions.end(); ++vertexOption ) {
      for ( std::vector<std::string>::const_iterator tauIdOption = tauIdOptions.begin();
            tauIdOption != tauIdOptions.end(); ++tauIdOption ) {
        for ( std::vector<std::string>::const_iterator absEtaRange = absEtaRanges.begin();
	      absEtaRange != absEtaRanges.end(); ++absEtaRange ) {
          for ( std::vector<std::string>::const_iterator l1MatchingOption = l1MatchingOptions.begin();
	        l1MatchingOption != l1MatchingOptions.end(); ++l1MatchingOption ) {        
            for ( std::vector<std::string>::const_iterator min_leadTrackPt = min_leadTrackPtValues.begin();
                  min_leadTrackPt != min_leadTrackPtValues.end(); ++min_leadTrackPt ) {              
              for ( std::vector<std::string>::const_iterator isolationWP = isolationWPs[*tauIdOption].begin();
	            isolationWP != isolationWPs[*tauIdOption].end(); ++isolationWP ) {
                for ( std::vector<std::string>::const_iterator process = processes.begin();
	              process != processes.end(); ++process ) {
                  std::string histogram2dName = Form(("%s/%s/" + dqmDirectory + "/numPFTaus_vs_ptThreshold_%s_%s_%s").data(), 
                    process->data(), srcVertices[*vertexOption].data(), 
                    pfAlgo->data(), vertexOption->data(), l1MatchingOption->data(), tauIdOption->data(),
                    absEtaRange->data(), min_leadTrackPt->data(), isolationWP->data());
                  TH2* histogram2d = loadHistogram2d(inputFiles[*process], histogram2dName);

                  TH1* histogram_rateSingleTau = makeRateHistogram(histogram2d, 1);
                  histograms_rateSingleTau_vs_processes[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption]
                    [*absEtaRange][*min_leadTrackPt][*isolationWP][*process] = histogram_rateSingleTau;
                  TH1* histogram_rateDoubleTau = makeRateHistogram(histogram2d, 2);
                  histograms_rateDoubleTau_vs_processes[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption]
                    [*absEtaRange][*min_leadTrackPt][*isolationWP][*process] = histogram_rateDoubleTau;
                } // process
                histograms_rateSingleTau_vs_isolationWPs[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption]
                  [*absEtaRange][*min_leadTrackPt][*isolationWP] = sumHistograms(
                    histograms_rateSingleTau_vs_processes[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption][*absEtaRange][*min_leadTrackPt][*isolationWP], processes);
                histograms_rateDoubleTau_vs_isolationWPs[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption]
                  [*absEtaRange][*min_leadTrackPt][*isolationWP] = sumHistograms(
                    histograms_rateDoubleTau_vs_processes[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption][*absEtaRange][*min_leadTrackPt][*isolationWP], processes);
                histograms_rateSingleTau_vs_leadTrackPt[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption]
                  [*absEtaRange][*isolationWP][*min_leadTrackPt] = sumHistograms(
                    histograms_rateSingleTau_vs_processes[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption][*absEtaRange][*min_leadTrackPt][*isolationWP], processes);              
                histograms_rateDoubleTau_vs_leadTrackPt[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption]
                  [*absEtaRange][*isolationWP][*min_leadTrackPt] = sumHistograms(
                    histograms_rateDoubleTau_vs_processes[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption][*absEtaRange][*min_leadTrackPt][*isolationWP], processes);
                histograms_rateSingleTau_vs_l1MatchingOption[*pfAlgo][*vertexOption][*tauIdOption]
                  [*absEtaRange][*min_leadTrackPt][*isolationWP][*l1MatchingOption] = sumHistograms(
                    histograms_rateSingleTau_vs_processes[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption][*absEtaRange][*min_leadTrackPt][*isolationWP], processes); 
                histograms_rateDoubleTau_vs_l1MatchingOption[*pfAlgo][*vertexOption][*tauIdOption]
                  [*absEtaRange][*min_leadTrackPt][*isolationWP][*l1MatchingOption] = sumHistograms(
                    histograms_rateDoubleTau_vs_processes[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption][*absEtaRange][*min_leadTrackPt][*isolationWP], processes); 
              } // isolationWP
            } // min_leadTrackPt
          } // l1MatchingOption

          for ( std::vector<std::string>::const_iterator l1MatchingOption = l1MatchingOptions.begin();
	        l1MatchingOption != l1MatchingOptions.end(); ++l1MatchingOption ) {
            for ( std::vector<std::string>::const_iterator min_leadTrackPt = min_leadTrackPtValues.begin();
                  min_leadTrackPt != min_leadTrackPtValues.end(); ++min_leadTrackPt ) {

              std::string isolationWP1 = ( isolationWPs[*tauIdOption].size() >= 1 ) ? isolationWPs[*tauIdOption][0] : "";
              std::string isolationWP2 = ( isolationWPs[*tauIdOption].size() >= 2 ) ? isolationWPs[*tauIdOption][1] : "";
              std::string isolationWP3 = ( isolationWPs[*tauIdOption].size() >= 3 ) ? isolationWPs[*tauIdOption][2] : "";
              std::string isolationWP4 = ( isolationWPs[*tauIdOption].size() >= 4 ) ? isolationWPs[*tauIdOption][3] : "";
              std::string isolationWP5 = ( isolationWPs[*tauIdOption].size() >= 5 ) ? isolationWPs[*tauIdOption][4] : "";
              std::string isolationWP6 = ( isolationWPs[*tauIdOption].size() >= 6 ) ? isolationWPs[*tauIdOption][5] : "";

              string_to_TH1Map1 histograms1 = histograms_rateSingleTau_vs_isolationWPs[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption][*absEtaRange][*min_leadTrackPt];
              std::string outputFileName1 = Form("makeRatePlots_SingleTau_%s%s%s_%s_%s_%s_vs_isolationWP.png", 
                pfAlgo->data(), vertexOption->data(), l1MatchingOption->data(),
                tauIdOption->data(), absEtaRange->data(), min_leadTrackPt->data());
              showHistograms(1150, 1150,
                             getHistogram(histograms1, isolationWP1), getLegendEntry(legendEntries_vs_isolationWPs[*tauIdOption], isolationWP1),
                             getHistogram(histograms1, isolationWP2), getLegendEntry(legendEntries_vs_isolationWPs[*tauIdOption], isolationWP2),
                             getHistogram(histograms1, isolationWP3), getLegendEntry(legendEntries_vs_isolationWPs[*tauIdOption], isolationWP3),
                             getHistogram(histograms1, isolationWP4), getLegendEntry(legendEntries_vs_isolationWPs[*tauIdOption], isolationWP4),
                             getHistogram(histograms1, isolationWP5), getLegendEntry(legendEntries_vs_isolationWPs[*tauIdOption], isolationWP5),
                             getHistogram(histograms1, isolationWP6), getLegendEntry(legendEntries_vs_isolationWPs[*tauIdOption], isolationWP6),
		             colors, lineStyles, 
  		             0.040, 0.66, 0.66, 0.23, 0.28,
		             labelTextLines, 0.050,
		             0.63, 0.66, 0.26, 0.07, 
		             -1., -1., "HLT #tau p_{T} Threshold [GeV]", 1.2, 
		             true, 1.e0, 1.e+6, "Single #tau Trigger Rate [Hz]", 1.4, 
		             outputFileName1);

              string_to_TH1Map1 histograms2 = histograms_rateDoubleTau_vs_isolationWPs[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption][*absEtaRange][*min_leadTrackPt];
              std::string outputFileName2 = Form("makeRatePlots_DoubleTau_%s%s%s_%s_%s_%s_vs_isolationWP.png", 
                pfAlgo->data(), vertexOption->data(), l1MatchingOption->data(),
                tauIdOption->data(), absEtaRange->data(), min_leadTrackPt->data());
              showHistograms(1150, 1150,
                             getHistogram(histograms2, isolationWP1), getLegendEntry(legendEntries_vs_isolationWPs[*tauIdOption], isolationWP1),
                             getHistogram(histograms2, isolationWP2), getLegendEntry(legendEntries_vs_isolationWPs[*tauIdOption], isolationWP2),
                             getHistogram(histograms2, isolationWP3), getLegendEntry(legendEntries_vs_isolationWPs[*tauIdOption], isolationWP3),
                             getHistogram(histograms2, isolationWP4), getLegendEntry(legendEntries_vs_isolationWPs[*tauIdOption], isolationWP4),
                             getHistogram(histograms2, isolationWP5), getLegendEntry(legendEntries_vs_isolationWPs[*tauIdOption], isolationWP5),
                             getHistogram(histograms2, isolationWP6), getLegendEntry(legendEntries_vs_isolationWPs[*tauIdOption], isolationWP6),
		             colors, lineStyles, 
		             0.040, 0.66, 0.66, 0.23, 0.28,
		             labelTextLines, 0.050,
		             0.63, 0.66, 0.26, 0.07, 
		             -1., -1., "HLT #tau p_{T} Threshold [GeV]", 1.2, 
		             true, 1.e0, 1.e+8, "Double #tau Trigger Rate [Hz]", 1.4, 
		             outputFileName2);
            } // min_leadTrackPt
          } // l1MatchingOption

          for ( std::vector<std::string>::const_iterator l1MatchingOption = l1MatchingOptions.begin();
	        l1MatchingOption != l1MatchingOptions.end(); ++l1MatchingOption ) {
            for ( std::vector<std::string>::const_iterator isolationWP = isolationWPs[*tauIdOption].begin();
	          isolationWP != isolationWPs[*tauIdOption].end(); ++isolationWP ) {
              string_to_TH1Map1 histograms3 = histograms_rateSingleTau_vs_leadTrackPt[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption][*absEtaRange][*isolationWP];
              std::string outputFileName3 = Form("makeRatePlots_SingleTau_%s%s%s_%s_%s_%s_vs_leadTrackPt.png", 
                pfAlgo->data(), vertexOption->data(), l1MatchingOption->data(),
                tauIdOption->data(), absEtaRange->data(), isolationWP->data());
              showHistograms(1150, 1150,
                             histograms3["leadTrackPtGt1"], legendEntries_vs_leadTrackPt["leadTrackPtGt1"],
	      	             histograms3["leadTrackPtGt2"], legendEntries_vs_leadTrackPt["leadTrackPtGt2"],
		             histograms3["leadTrackPtGt5"], legendEntries_vs_leadTrackPt["leadTrackPtGt5"],
		             nullptr, "",
		             nullptr, "",
		             nullptr, "",
		             colors, lineStyles, 
		             0.040, 0.47, 0.79, 0.42, 0.13,
		             labelTextLines, 0.050,
		             0.63, 0.66, 0.26, 0.07, 
		             -1., -1., "HLT #tau p_{T} Threshold [GeV]", 1.2, 
		             true, 1.e0, 1.e+8, "Single #tau Trigger Rate [Hz]", 1.4, 
		             outputFileName3);

              string_to_TH1Map1 histograms4 = histograms_rateDoubleTau_vs_leadTrackPt[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption][*absEtaRange][*isolationWP];
              std::string outputFileName4 = Form("makeRatePlots_DoubleTau_%s%s%s_%s_%s_%s_vs_leadTrackPt.png", 
                pfAlgo->data(), vertexOption->data(), l1MatchingOption->data(),
                tauIdOption->data(), absEtaRange->data(), isolationWP->data());
              showHistograms(1150, 1150,
                             histograms4["leadTrackPtGt1"], legendEntries_vs_leadTrackPt["leadTrackPtGt1"],
	         	     histograms4["leadTrackPtGt2"], legendEntries_vs_leadTrackPt["leadTrackPtGt2"],
	 	             histograms4["leadTrackPtGt5"], legendEntries_vs_leadTrackPt["leadTrackPtGt5"],
		             nullptr, "",
		             nullptr, "",
		             nullptr, "",
		             colors, lineStyles, 
		             0.040, 0.47, 0.79, 0.42, 0.13,
		             labelTextLines, 0.050,
		             0.63, 0.66, 0.26, 0.07, 
		             -1., -1., "HLT #tau p_{T} Threshold [GeV]", 1.2, 
		             true, 1.e0, 1.e+8, "Double #tau Trigger Rate [Hz]", 1.4, 
		             outputFileName4);
            } // isolationWP
          } // l1MatchingOption

          for ( std::vector<std::string>::const_iterator min_leadTrackPt = min_leadTrackPtValues.begin();
                min_leadTrackPt != min_leadTrackPtValues.end(); ++min_leadTrackPt ) {
            for ( std::vector<std::string>::const_iterator isolationWP = isolationWPs[*tauIdOption].begin();
	          isolationWP != isolationWPs[*tauIdOption].end(); ++isolationWP ) {
              string_to_TH1Map1 histograms5 = histograms_rateSingleTau_vs_l1MatchingOption[*pfAlgo][*vertexOption][*tauIdOption][*absEtaRange][*min_leadTrackPt][*isolationWP];
              std::string outputFileName5 = Form("makeRatePlots_SingleTau_%s%s_%s_%s_%s_%s_vs_l1MatchingOption.png", 
                pfAlgo->data(), vertexOption->data(), 
                tauIdOption->data(), absEtaRange->data(), min_leadTrackPt->data(), isolationWP->data());
              showHistograms(1150, 1150,
                             histograms5[""],            legendEntries_vs_l1MatchingOption[""],
	  	             histograms5["MatchedToL1"], legendEntries_vs_l1MatchingOption["MatchedToL1"],
		             nullptr, "",
		             nullptr, "",
		             nullptr, "",
		             nullptr, "",
		             colors, lineStyles, 
		             0.040, 0.47, 0.79, 0.42, 0.13,
		             labelTextLines, 0.050,
		             0.63, 0.66, 0.26, 0.07, 
		             -1., -1., "HLT #tau p_{T} Threshold [GeV]", 1.2, 
		             true, 1.e0, 1.e+8, "Single #tau Trigger Rate [Hz]", 1.4, 
		             outputFileName5);

              string_to_TH1Map1 histograms6 = histograms_rateDoubleTau_vs_l1MatchingOption[*pfAlgo][*vertexOption][*tauIdOption][*absEtaRange][*min_leadTrackPt][*isolationWP];
              std::string outputFileName6 = Form("makeRatePlots_DoubleTau_%s%s_%s_%s_%s_%s_vs_l1MatchingOption.png", 
                pfAlgo->data(), vertexOption->data(),
                tauIdOption->data(), absEtaRange->data(), min_leadTrackPt->data(), isolationWP->data());
              showHistograms(1150, 1150,
	  	             histograms6[""],            legendEntries_vs_l1MatchingOption[""],
		             histograms6["MatchedToL1"], legendEntries_vs_l1MatchingOption["MatchedToL1"],
		             nullptr, "",
		             nullptr, "",
		             nullptr, "",
		             nullptr, "",
		             colors, lineStyles, 
		             0.040, 0.47, 0.79, 0.42, 0.13,
		             labelTextLines, 0.050,
		             0.63, 0.66, 0.26, 0.07, 
		             -1., -1., "HLT #tau p_{T} Threshold [GeV]", 1.2, 
		             true, 1.e0, 1.e+8, "Double #tau Trigger Rate [Hz]", 1.4, 
		             outputFileName6);
            } // isolationWP
          } // min_leadTrackPt

          for ( std::vector<std::string>::const_iterator l1MatchingOption = l1MatchingOptions.begin();
	        l1MatchingOption != l1MatchingOptions.end(); ++l1MatchingOption ) {
            for ( std::vector<std::string>::const_iterator min_leadTrackPt = min_leadTrackPtValues.begin();
                  min_leadTrackPt != min_leadTrackPtValues.end(); ++min_leadTrackPt ) {
              for ( std::vector<std::string>::const_iterator isolationWP = isolationWPs[*tauIdOption].begin();
	            isolationWP != isolationWPs[*tauIdOption].end(); ++isolationWP ) {
                string_to_TH1Map1 histograms7 = histograms_rateSingleTau_vs_processes[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption][*absEtaRange][*min_leadTrackPt][*isolationWP];
                std::string outputFileName7 = Form("makeRatePlots_SingleTau_%s%s%s_%s_%s_%s_%s_vs_processes.png", 
                  pfAlgo->data(), vertexOption->data(), l1MatchingOption->data(), 
                  tauIdOption->data(), absEtaRange->data(), min_leadTrackPt->data(), isolationWP->data());
                showHistograms(1150, 1150,
                               histograms7["minbias"], legendEntries_vs_processes["minbias"],
                               histograms7["QCD"],     legendEntries_vs_processes["QCD"],
                               histograms7["DY"],      legendEntries_vs_processes["DY"],
                               histograms7["W"],       legendEntries_vs_processes["W"],
                               nullptr, "",
		               nullptr, "",
		               colors, lineStyles, 
		               0.040, 0.65, 0.66, 0.24, 0.28,
		               labelTextLines, 0.050,
		               0.63, 0.66, 0.26, 0.07, 
		               -1., -1., "HLT #tau p_{T} Threshold [GeV]", 1.2, 
		               true, 1.e-1, 1.e+6, "Single #tau Trigger Rate [Hz]", 1.4, 
		               outputFileName7);

                string_to_TH1Map1 histograms8 = histograms_rateDoubleTau_vs_processes[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption][*absEtaRange][*min_leadTrackPt][*isolationWP];
                std::string outputFileName8 = Form("makeRatePlots_DoubleTau_%s%s%s_%s_%s_%s_%s_vs_processes.png", 
                  pfAlgo->data(), vertexOption->data(), l1MatchingOption->data(), 
                  tauIdOption->data(), absEtaRange->data(), min_leadTrackPt->data(), isolationWP->data());
                showHistograms(1150, 1150,
                               histograms8["minbias"], legendEntries_vs_processes["minbias"],
                               histograms8["QCD"],     legendEntries_vs_processes["QCD"],
                               histograms8["DY"],      legendEntries_vs_processes["DY"],
                               histograms8["W"],       legendEntries_vs_processes["W"],
                               nullptr, "",
		               nullptr, "",
		               colors, lineStyles, 
		               0.040, 0.65, 0.66, 0.24, 0.28,
		               labelTextLines, 0.050,
		               0.63, 0.66, 0.26, 0.07, 
		               -1., -1., "HLT #tau p_{T} Threshold [GeV]", 1.2, 
		               true, 1.e-1, 1.e+6, "Double #tau Trigger Rate [Hz]", 1.4, 
		               outputFileName8);
              } // isolationWP
            } // min_leadTrackPt
          } // l1MatchingOption
        } // absEtaRange
      } // tauIdOption
    } // vertexOption
  } // pfAlgo

  for ( std::map<std::string, TFile*>::const_iterator inputFile = inputFiles.begin();
        inputFile != inputFiles.end(); ++inputFile ) {
    delete inputFile->second;
  }
}

