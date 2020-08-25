
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TF1.h>
#include <TColor.h>
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

double square(double x)
{
  return x*x;
}

std::pair<double, double> 
compRate(TH2* histogram2d, int idxBinX, int minNumTaus)
{
  const TAxis* yAxis = histogram2d->GetYaxis();
  int numBinsY = yAxis->GetNbins();
  double integral_passed = 0.;
  double integralErr2_passed = 0.;
  for ( int idxBinY = 1; idxBinY <= numBinsY; ++idxBinY ) { 
    double binContent = histogram2d->GetBinContent(idxBinX, idxBinY);
    double binError = histogram2d->GetBinError(idxBinX, idxBinY);
    double numTaus = yAxis->GetBinCenter(idxBinY);
    if ( numTaus >= minNumTaus ) 
    {
      integral_passed += binContent;
      integralErr2_passed += square(binError);
    }
  }
  double rate = integral_passed;
  double rateErr = TMath::Sqrt(integralErr2_passed);
  const double max_rate = 2.8e+7; // bunch-crossing frequency of colliding bunches = 28 MHz
  if ( rate > max_rate ) 
  { 
    std::cerr << "Warning in <compRate>: computed rate = " << rate << " exceeds bunch-crossing frequency of 28 MHz !!" << std::endl;
    rate = max_rate;
    rateErr = 0.;
  }
  return std::pair<double, double>(rate, rateErr);
}

TH1* makeRateHistogram(TH2* histogram2d, int minNumTaus)
{
  std::string histogramName_rate = Form("%s_rate_minNumTausEq%i", histogram2d->GetName(), minNumTaus);
  TH1* histogram_rate = histogram2d->ProjectionX(histogramName_rate.data());
  histogram_rate->Reset();
  const TAxis* xAxis = histogram2d->GetXaxis();
  int numBinsX = xAxis->GetNbins();
  for ( int idxBinX = 1; idxBinX <= numBinsX; ++idxBinX ) { 
    std::pair<double, double> rate_and_rateErr = compRate(histogram2d, idxBinX, minNumTaus);
    histogram_rate->SetBinContent(idxBinX, rate_and_rateErr.first);
    histogram_rate->SetBinError(idxBinX, rate_and_rateErr.second);
  }
  return histogram_rate;
}

TH1* sumHistograms(const std::vector<TH1*>& histograms)
{
  TH1* histogram_sum = nullptr;
  for ( std::vector<TH1*>::const_iterator histogram = histograms.begin();
        histogram != histograms.end(); ++histogram ) {
    if ( histogram_sum ) {
      histogram_sum->Add(*histogram);
    } else {
      std::string histogramName_sum = Form("%s_sum", (*histogram)->GetName());
      histogram_sum = (TH1*)(*histogram)->Clone(histogramName_sum.data());
      if ( !histogram_sum->GetSumw2N() ) histogram_sum->Sumw2();
    }
  }
  return histogram_sum;
}

void setStyle_uncertainty(TH1* histogram)
{
  const int color_int = 12;
  const double alpha = 0.40;
  TColor* color = gROOT->GetColor(color_int);
  static int newColor_int = -1;
  static TColor* newColor = 0;
  if ( !newColor ) {
    newColor_int = gROOT->GetListOfColors()->GetSize() + 1;
    newColor = new TColor(newColor_int, color->GetRed(), color->GetGreen(), color->GetBlue(), "", alpha);
  }
  histogram->SetLineColor(newColor_int);
  histogram->SetLineWidth(0);
  histogram->SetFillColor(newColor_int);
  histogram->SetFillStyle(1001);
}

void showHistograms(double canvasSizeX, double canvasSizeY,
                    TH1* histogram_minbias_unstitched, const std::string& legendEntry_minbias_unstitched,
                    TH1* histogram_minbias_stitched, const std::string& legendEntry_minbias_stitched,
                    TH1* histogram_qcd1_stitched, const std::string& legendEntry_qcd1_stitched,
                    TH1* histogram_qcd2_stitched, const std::string& legendEntry_qcd2_stitched,
                    TH1* histogram_qcd3_stitched, const std::string& legendEntry_qcd3_stitched,
                    TH1* histogram_qcd4_stitched, const std::string& legendEntry_qcd4_stitched,
                    TH1* histogram_qcd5_stitched, const std::string& legendEntry_qcd5_stitched,
                    TH1* histogram_qcd6_stitched, const std::string& legendEntry_qcd6_stitched,
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
  canvas->Draw();

  TPad* topPad = new TPad("topPad", "topPad", 0.00, 0.35, 1.00, 1.00);
  topPad->SetFillColor(10);
  topPad->SetTopMargin(0.015);
  topPad->SetLeftMargin(0.17);
  topPad->SetBottomMargin(0.03);
  topPad->SetRightMargin(0.03);
  topPad->SetLogy(useLogScale);
  topPad->SetGridx(1);
  topPad->SetGridy(1);

  TPad* bottomPad = new TPad("bottomPad", "bottomPad", 0.00, 0.01, 1.00, 0.35);
  bottomPad->SetFillColor(10);
  bottomPad->SetTopMargin(0.02);
  bottomPad->SetLeftMargin(0.17);
  bottomPad->SetBottomMargin(0.26);
  bottomPad->SetRightMargin(0.03);
  bottomPad->SetLogy(false);

  canvas->cd();
  topPad->Draw();
  topPad->cd();

  if ( !(histogram_minbias_unstitched && histogram_minbias_stitched) ) {
    std::cerr << "<showHistograms>: histogram_minbias = NULL --> skipping !!" << std::endl;
    return;
  }

  THStack* histogramStack = new THStack("histogramStack", "histogramStack");

  histogramStack->SetTitle("");
  histogramStack->SetMinimum(yMin);
  histogramStack->SetMaximum(yMax);

  if ( histogram_qcd6_stitched ) {
    histogram_qcd6_stitched->SetStats(false);
    histogram_qcd6_stitched->SetFillColor(colors[5]);
    histogram_qcd6_stitched->SetFillStyle(fillStyles[5]);
    histogram_qcd6_stitched->SetLineColor(colors[5]);
    histogram_qcd6_stitched->SetLineWidth(2);
    histogram_qcd6_stitched->SetLineStyle(1);
    histogramStack->Add(histogram_qcd6_stitched);
  }

  if ( histogram_qcd5_stitched ) {
    histogram_qcd5_stitched->SetStats(false);
    histogram_qcd5_stitched->SetFillColor(colors[4]);
    histogram_qcd5_stitched->SetFillStyle(fillStyles[4]);
    histogram_qcd5_stitched->SetLineColor(colors[4]);
    histogram_qcd5_stitched->SetLineWidth(2);
    histogram_qcd5_stitched->SetLineStyle(1);
    histogramStack->Add(histogram_qcd5_stitched);
  }

  if ( histogram_qcd4_stitched ) {
    histogram_qcd4_stitched->SetStats(false);
    histogram_qcd4_stitched->SetFillColor(colors[3]);
    histogram_qcd4_stitched->SetFillStyle(fillStyles[3]);
    histogram_qcd4_stitched->SetLineColor(colors[3]);
    histogram_qcd4_stitched->SetLineWidth(2);
    histogram_qcd4_stitched->SetLineStyle(1);
    histogramStack->Add(histogram_qcd4_stitched);
  }

  if ( histogram_qcd3_stitched ) {
    histogram_qcd3_stitched->SetStats(false);
    histogram_qcd3_stitched->SetFillColor(colors[2]);
    histogram_qcd3_stitched->SetFillStyle(fillStyles[2]);
    histogram_qcd3_stitched->SetLineColor(colors[2]);
    histogram_qcd3_stitched->SetLineWidth(2);
    histogram_qcd3_stitched->SetLineStyle(1);
    histogramStack->Add(histogram_qcd3_stitched);
  }

  if ( histogram_qcd2_stitched ) {
    histogram_qcd2_stitched->SetStats(false);
    histogram_qcd2_stitched->SetFillColor(colors[1]);
    histogram_qcd2_stitched->SetFillStyle(fillStyles[1]);
    histogram_qcd2_stitched->SetLineColor(colors[1]);
    histogram_qcd2_stitched->SetLineWidth(2);
    histogram_qcd2_stitched->SetLineStyle(1);
    histogramStack->Add(histogram_qcd2_stitched);
  }

  if ( histogram_qcd1_stitched ) {
    histogram_qcd1_stitched->SetStats(false);
    histogram_qcd1_stitched->SetFillColor(colors[0]);
    histogram_qcd1_stitched->SetFillStyle(fillStyles[0]);
    histogram_qcd1_stitched->SetLineColor(colors[0]);
    histogram_qcd1_stitched->SetLineWidth(2);
    histogram_qcd1_stitched->SetLineStyle(1);
    histogramStack->Add(histogram_qcd1_stitched);
  }

  if ( histogram_minbias_stitched ) {
    histogram_minbias_stitched->SetStats(false);
    histogram_minbias_stitched->SetFillColor(10);
    histogram_minbias_stitched->SetFillStyle(1001);
    histogram_minbias_stitched->SetLineColor(1);
    histogram_minbias_stitched->SetLineWidth(2);
    histogram_minbias_stitched->SetLineStyle(1);
    histogramStack->Add(histogram_minbias_stitched);
  }

  histogramStack->Draw("hist");

  std::vector<TH1*> histograms_stitched;
  histograms_stitched.push_back(histogram_minbias_stitched);
  if ( histogram_qcd1_stitched ) histograms_stitched.push_back(histogram_qcd1_stitched);
  if ( histogram_qcd2_stitched ) histograms_stitched.push_back(histogram_qcd2_stitched);
  if ( histogram_qcd3_stitched ) histograms_stitched.push_back(histogram_qcd3_stitched);
  if ( histogram_qcd4_stitched ) histograms_stitched.push_back(histogram_qcd4_stitched);
  if ( histogram_qcd5_stitched ) histograms_stitched.push_back(histogram_qcd5_stitched);
  if ( histogram_qcd6_stitched ) histograms_stitched.push_back(histogram_qcd6_stitched);
  TH1* histogramErr_stitched = sumHistograms(histograms_stitched);
  setStyle_uncertainty(histogramErr_stitched);
  histogramErr_stitched->Draw("e2same");
  
  //std::cout << "histogramStack->GetHistogram() = " << histogramStack->GetHistogram() << std::endl;

  TAxis* xAxis_top = histogramStack->GetHistogram()->GetXaxis();
  xAxis_top->SetTitle(xAxisTitle.data());
  xAxis_top->SetTitleOffset(1.20);
  xAxis_top->SetLabelColor(10);
  xAxis_top->SetTitleColor(10);
  xAxis_top->SetNdivisions(505);
  if ( xMax > xMin ) {
    std::cout << "limiting x-axis range to " << xMin << ".." << xMax << std::endl;
    xAxis_top->SetRangeUser(xMin, xMax);
  }

  TAxis* yAxis_top = histogramStack->GetHistogram()->GetYaxis();
  yAxis_top->SetTitle(yAxisTitle.data());
  yAxis_top->SetTitleOffset(yAxisOffset);
  yAxis_top->SetTitleSize(0.060);
  yAxis_top->SetLabelSize(0.045);
  yAxis_top->SetTickLength(0.02);  

  histogram_minbias_unstitched->SetTitle("");
  histogram_minbias_unstitched->SetStats(false);
  histogram_minbias_unstitched->SetMarkerColor(1);
  histogram_minbias_unstitched->SetMarkerStyle(20);
  histogram_minbias_unstitched->SetMarkerSize(1);
  histogram_minbias_unstitched->SetLineColor(1);
  histogram_minbias_unstitched->SetLineWidth(2);
  histogram_minbias_unstitched->Draw("e1psame");

  TLegend* legend = 0;
  if ( legendEntry_minbias_stitched != "" ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, "", "brNDC"); 
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetTextSize(legendTextSize);
    if ( legendEntry_minbias_unstitched != "" ) legend->AddEntry(histogram_minbias_unstitched, legendEntry_minbias_unstitched.data(), "p");
    legend->AddEntry(histogram_minbias_stitched, legendEntry_minbias_stitched.data(), "f");
    if ( histogram_qcd1_stitched ) legend->AddEntry(histogram_qcd1_stitched, legendEntry_qcd1_stitched.data(), "f");
    if ( histogram_qcd2_stitched ) legend->AddEntry(histogram_qcd2_stitched, legendEntry_qcd2_stitched.data(), "f");
    if ( histogram_qcd3_stitched ) legend->AddEntry(histogram_qcd3_stitched, legendEntry_qcd3_stitched.data(), "f");
    if ( histogram_qcd4_stitched ) legend->AddEntry(histogram_qcd4_stitched, legendEntry_qcd4_stitched.data(), "f");
    if ( histogram_qcd5_stitched ) legend->AddEntry(histogram_qcd5_stitched, legendEntry_qcd5_stitched.data(), "f");
    if ( histogram_qcd6_stitched ) legend->AddEntry(histogram_qcd6_stitched, legendEntry_qcd6_stitched.data(), "f");
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

  TF1* line1 = new TF1("line1", "150", xAxis_top->GetXmin(), xAxis_top->GetXmax());
  line1->SetLineStyle(3);
  line1->SetLineWidth(2);
  line1->SetLineColor(kRed);
  line1->Draw("same");

  histogram_minbias_unstitched->Draw("e1psame");
  histogram_minbias_unstitched->Draw("axissame");

  canvas->cd();
  bottomPad->Draw();
  bottomPad->cd();

  TH1* histogramRatio = (TH1*)histogram_minbias_unstitched->Clone("histogramRatio");
  if ( !histogramRatio->GetSumw2N() ) histogramRatio->Sumw2();
  histogramRatio->SetTitle("");
  histogramRatio->SetStats(false);
  histogramRatio->SetMinimum(0.01);
  histogramRatio->SetMaximum(1.99);
  histogramRatio->SetMarkerColor(histogram_minbias_unstitched->GetMarkerColor());
  histogramRatio->SetMarkerStyle(histogram_minbias_unstitched->GetMarkerStyle());
  histogramRatio->SetMarkerSize(histogram_minbias_unstitched->GetMarkerSize());
  histogramRatio->SetLineColor(histogram_minbias_unstitched->GetLineColor());
  
  TH1* histogramRatioUncertainty = (TH1*)histogram_minbias_unstitched->Clone("histogramRatioUncertainty");
  if ( !histogramRatioUncertainty->GetSumw2N() ) histogramRatioUncertainty->Sumw2();
  histogramRatioUncertainty->SetMarkerColor(10);
  histogramRatioUncertainty->SetMarkerSize(0);
  histogramRatioUncertainty->SetFillColor(histogramErr_stitched->GetFillColor());
  histogramRatioUncertainty->SetFillStyle(histogramErr_stitched->GetFillStyle());

  int numBins_bottom = histogramRatio->GetNbinsX();
  for ( int idxBin = 1; idxBin <= numBins_bottom; ++idxBin ) {
    double binContent_unstitched = histogram_minbias_unstitched->GetBinContent(idxBin);
    double binError_unstitched = histogram_minbias_unstitched->GetBinError(idxBin);
    double binContent_stitched = histogramErr_stitched->GetBinContent(idxBin);
    double binError_stitched = histogramErr_stitched->GetBinError(idxBin);

    if ( binContent_stitched > 0. ) {
      histogramRatio->SetBinContent(idxBin, binContent_unstitched/binContent_stitched);
      histogramRatio->SetBinError(idxBin, binError_unstitched/binContent_stitched);

      histogramRatioUncertainty->SetBinContent(idxBin, 1.);
      histogramRatioUncertainty->SetBinError(idxBin, binError_stitched/binContent_stitched);
    }
  }

  TAxis* xAxis_bottom = histogramRatio->GetXaxis();
  xAxis_bottom->SetTitle(xAxis_top->GetTitle());
  xAxis_bottom->SetLabelColor(1);
  xAxis_bottom->SetTitleColor(1);
  xAxis_bottom->SetTitleOffset(1.30);
  xAxis_bottom->SetTitleSize(0.100);
  xAxis_bottom->SetTitleFont(xAxis_top->GetTitleFont());
  xAxis_bottom->SetLabelOffset(0.02);
  xAxis_bottom->SetLabelSize(0.085);
  xAxis_bottom->SetTickLength(0.065);
  xAxis_bottom->SetNdivisions(505);
  if ( xMax > xMin ) {
    xAxis_bottom->SetRangeUser(xMin, xMax);
  }

  TAxis* yAxis_bottom = histogramRatio->GetYaxis();
  yAxis_bottom->SetTitle("Ratio");
  yAxis_bottom->SetLabelColor(1);
  yAxis_bottom->SetTitleColor(1);
  yAxis_bottom->SetTitleOffset(0.85);
  yAxis_bottom->SetTitleFont(yAxis_top->GetTitleFont());
  yAxis_bottom->SetNdivisions(505);
  yAxis_bottom->CenterTitle();
  yAxis_bottom->SetTitleSize(0.095);
  yAxis_bottom->SetLabelSize(0.075);
  yAxis_bottom->SetTickLength(0.040);  

  histogramRatio->Draw("ep");

  TF1* line2 = new TF1("line2", "1", xAxis_bottom->GetXmin(), xAxis_bottom->GetXmax());
  line2->SetLineStyle(3);
  line2->SetLineWidth(2);
  line2->SetLineColor(kBlack);
  line2->Draw("same");

  histogramRatioUncertainty->Draw("e2same");  
  line2->Draw("same");
  histogramRatio->Draw("epsame");

  canvas->Update();

  canvas->Update();
  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  //canvas->Print(std::string(outputFileName_plot).append(".png").data());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  
  delete histogramStack;
  delete label;
  delete legend;
  delete line1;
  delete topPad;
  delete histogramRatio;
  delete histogramRatioUncertainty;
  delete line2;
  delete bottomPad;
  delete canvas;  
}

void makeRatePlots_stacked()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  std::string inputFilePath = "/hdfs/local/veelken/Phase2HLT/rate/2020Aug22v2/";

  std::vector<std::string> processes;
  processes.push_back("minbias");
  processes.push_back("qcd_pt30to50");
  processes.push_back("qcd_pt50to80");
  processes.push_back("qcd_pt80to120");
  processes.push_back("qcd_pt120to170");
  processes.push_back("qcd_pt170to300");

  std::map<std::string, std::string> inputFileNames; // key = process
  inputFileNames["minbias"]        = "hadd_stage1_minbias_hps_offlinePrimaryVertices_dz_wrt_primaryVertex_8Hits.root";
  inputFileNames["qcd_pt30to50"]   = "hadd_stage1_qcd_pt30to50_hps_offlinePrimaryVertices_dz_wrt_primaryVertex_8Hits.root";
  inputFileNames["qcd_pt50to80"]   = "hadd_stage1_qcd_pt50to80_hps_offlinePrimaryVertices_dz_wrt_primaryVertex_8Hits.root";
  inputFileNames["qcd_pt80to120"]  = "hadd_stage1_qcd_pt80to120_hps_offlinePrimaryVertices_dz_wrt_primaryVertex_8Hits.root";
  inputFileNames["qcd_pt120to170"] = "hadd_stage1_qcd_pt120to170_hps_offlinePrimaryVertices_dz_wrt_primaryVertex_8Hits.root";
  inputFileNames["qcd_pt170to300"] = "hadd_stage1_qcd_pt170to300_hps_offlinePrimaryVertices_dz_wrt_primaryVertex_8Hits.root";

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
  pfAlgos.push_back("HpsPFTau");

  std::vector<std::string> vertexOptions;
  vertexOptions.push_back("8HitsMaxDeltaZWithOfflineVertices");

  std::map<std::string, std::string> srcVertices; // key = vertexOption
  srcVertices["8HitsMaxDeltaZWithOfflineVertices"]                  = "offlinePrimaryVertices";
  srcVertices["8HitsMaxDeltaZToLeadTrackWithOfflineVertices"]       = "offlinePrimaryVertices";
  srcVertices["8HitsMaxDeltaZWithOnlineVertices"]                   = "hltPhase2PixelVertices";
  srcVertices["8HitsMaxDeltaZToLeadTrackWithOnlineVertices"]        = "hltPhase2PixelVertices";
  srcVertices["8HitsMaxDeltaZWithOnlineVerticesTrimmed"]            = "hltPhase2TrimmedPixelVertices";
  srcVertices["8HitsMaxDeltaZToLeadTrackWithOnlineVerticesTrimmed"] = "hltPhase2TrimmedPixelVertices";
  
  std::vector<std::string> tauIdOptions;
  tauIdOptions.push_back("recoSumChargedIso");
  tauIdOptions.push_back("patSumChargedIso");
  tauIdOptions.push_back("patDeepTau");

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
  min_leadTrackPtValues.push_back("leadTrackPtGt5");

  std::map<std::string, std::vector<std::string>> isolationWPs; // key = tauIdOption
  isolationWPs["recoSumChargedIso"].push_back("noIsolation");
  isolationWPs["recoSumChargedIso"].push_back("relDiscriminatorLt0p400");
  isolationWPs["recoSumChargedIso"].push_back("relDiscriminatorLt0p200");
  isolationWPs["recoSumChargedIso"].push_back("relDiscriminatorLt0p100");
  isolationWPs["recoSumChargedIso"].push_back("relDiscriminatorLt0p050");
  isolationWPs["patSumChargedIso"] = isolationWPs["recoSumChargedIso"];
  isolationWPs["patDeepTau"].push_back("absDiscriminatorGt0p260");
  isolationWPs["patDeepTau"].push_back("absDiscriminatorGt0p425");
  isolationWPs["patDeepTau"].push_back("absDiscriminatorGt0p598");
  isolationWPs["patDeepTau"].push_back("absDiscriminatorGt0p785");
  isolationWPs["patDeepTau"].push_back("absDiscriminatorGt0p883");
  isolationWPs["patDeepTau"].push_back("absDiscriminatorGt0p931");

  std::map<std::string, std::string> legendEntries; // key = process
  legendEntries["minbias"]        = "Minbias";
  legendEntries["qcd_pt30to50"]   = " 30 < #hat{p}_{T} <  50";
  legendEntries["qcd_pt50to80"]   = " 50 < #hat{p}_{T} <  80";
  legendEntries["qcd_pt80to120"]  = " 80 < #hat{p}_{T} < 120";
  legendEntries["qcd_pt120to170"] = "120 < #hat{p}_{T} < 170";
  legendEntries["qcd_pt170to300"] = "170 < #hat{p}_{T} < 300";

  std::map<std::string, std::string> dqmDirectories; // key = process
  dqmDirectories["minbias"]        = "minbias/%s/%s/%sAnalyzerBackground%s%s_%s";
  dqmDirectories["qcd_pt30to50"]   = "QCD/%s/%s/%sAnalyzerBackground%s%s_%s";
  dqmDirectories["qcd_pt50to80"]   = "QCD/%s/%s/%sAnalyzerBackground%s%s_%s";
  dqmDirectories["qcd_pt80to120"]  = "QCD/%s/%s/%sAnalyzerBackground%s%s_%s";
  dqmDirectories["qcd_pt120to170"] = "QCD/%s/%s/%sAnalyzerBackground%s%s_%s";
  dqmDirectories["qcd_pt170to300"] = "QCD/%s/%s/%sAnalyzerBackground%s%s_%s";

  std::vector<std::string> labelTextLines;

  int colors[6]     = {    7,    6,    4,   28,    8,    2 };
  int fillStyles[6] = { 1001, 1001, 1001, 1001, 1001, 1001 };

  typedef std::map<std::string, TH1*>              string_to_TH1Map1;
  typedef std::map<std::string, string_to_TH1Map1> string_to_TH1Map2;
  typedef std::map<std::string, string_to_TH1Map2> string_to_TH1Map3;
  typedef std::map<std::string, string_to_TH1Map3> string_to_TH1Map4;
  typedef std::map<std::string, string_to_TH1Map4> string_to_TH1Map5;
  typedef std::map<std::string, string_to_TH1Map5> string_to_TH1Map6;
  typedef std::map<std::string, string_to_TH1Map6> string_to_TH1Map7;
  typedef std::map<std::string, string_to_TH1Map7> string_to_TH1Map8;
  string_to_TH1Map7 histograms_rateSingleTau_minbias_unstitched;    // key = pfAlgo, vertexOption, tauIdOption, l1MatchingOption, absEtaRange, min_leadTrackPt, isolationWP
  string_to_TH1Map8 histograms_rateSingleTau_vs_processes_stitched; // key = pfAlgo, vertexOption, tauIdOption, l1MatchingOption, absEtaRange, min_leadTrackPt, isolationWP, process
  string_to_TH1Map7 histograms_rateDoubleTau_minbias_unstitched;    // key = pfAlgo, vertexOption, tauIdOption, l1MatchingOption, absEtaRange, min_leadTrackPt, isolationWP
  string_to_TH1Map8 histograms_rateDoubleTau_vs_processes_stitched; // key = pfAlgo, vertexOption, tauIdOption, l1MatchingOption, absEtaRange, min_leadTrackPt, isolationWP, process

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

                std::string histogram2dName_minbias_unstitched = Form((dqmDirectories["minbias"] + "/numPFTaus_vs_ptThreshold_%s_%s_%s_dzLt0p2").data(), 
                  srcVertices[*vertexOption].data(), "lumiScale",
                  pfAlgo->data(), vertexOption->data(), l1MatchingOption->data(), tauIdOption->data(),
                  absEtaRange->data(), min_leadTrackPt->data(), isolationWP->data());
                TH2* histogram2d_minbias_unstitched = loadHistogram2d(inputFiles["minbias"], histogram2dName_minbias_unstitched);

                TH1* histogram_rateSingleTau_minbias_unstitched = makeRateHistogram(histogram2d_minbias_unstitched, 1);
                histograms_rateSingleTau_minbias_unstitched[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption]
                  [*absEtaRange][*min_leadTrackPt][*isolationWP] = histogram_rateSingleTau_minbias_unstitched;
                TH1* histogram_rateDoubleTau_minbias_unstitched = makeRateHistogram(histogram2d_minbias_unstitched, 2);
                histograms_rateDoubleTau_minbias_unstitched[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption]
                  [*absEtaRange][*min_leadTrackPt][*isolationWP] = histogram_rateDoubleTau_minbias_unstitched;

                for ( std::vector<std::string>::const_iterator process = processes.begin();
	              process != processes.end(); ++process ) {
                  std::string histogram2dName_process_stitched = Form((dqmDirectories[*process] + "/numPFTaus_vs_ptThreshold_%s_%s_%s_dzLt0p2").data(), 
                    srcVertices[*vertexOption].data(), "stitchingWeight",
                    pfAlgo->data(), vertexOption->data(), l1MatchingOption->data(), tauIdOption->data(),
                    absEtaRange->data(), min_leadTrackPt->data(), isolationWP->data());
                  TH2* histogram2d_process_stitched = loadHistogram2d(inputFiles[*process], histogram2dName_process_stitched);

                  TH1* histogram_rateSingleTau_process_stitched = makeRateHistogram(histogram2d_process_stitched, 1);
                  histograms_rateSingleTau_vs_processes_stitched[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption]
                    [*absEtaRange][*min_leadTrackPt][*isolationWP][*process] = histogram_rateSingleTau_process_stitched;
                  TH1* histogram_rateDoubleTau_process_stitched = makeRateHistogram(histogram2d_process_stitched, 2);
                  histograms_rateDoubleTau_vs_processes_stitched[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption]
                    [*absEtaRange][*min_leadTrackPt][*isolationWP][*process] = histogram_rateDoubleTau_process_stitched;
                } // process

                TH1* histogram1_minbias_unstitched = histograms_rateSingleTau_minbias_unstitched[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption]
                  [*absEtaRange][*min_leadTrackPt][*isolationWP];
                string_to_TH1Map1 histograms1_processes_stitched = histograms_rateSingleTau_vs_processes_stitched[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption]
                  [*absEtaRange][*min_leadTrackPt][*isolationWP];
                std::string outputFileName1 = Form("makeRatePlots_stacked_SingleTau_%s%s%s_%s_%s_%s_%s.png", 
                  pfAlgo->data(), vertexOption->data(), l1MatchingOption->data(), 
                    tauIdOption->data(), absEtaRange->data(), min_leadTrackPt->data(), isolationWP->data());
                showHistograms(1150, 1600,
                  histogram1_minbias_unstitched,                    "",
                  histograms1_processes_stitched["minbias"],        legendEntries["minbias"],
                  histograms1_processes_stitched["qcd_pt30to50"],   legendEntries["qcd_pt30to50"],
                  histograms1_processes_stitched["qcd_pt50to80"],   legendEntries["qcd_pt50to80"],
                  histograms1_processes_stitched["qcd_pt80to120"],  legendEntries["qcd_pt80to120"],
                  histograms1_processes_stitched["qcd_pt120to170"], legendEntries["qcd_pt120to170"],
                  histograms1_processes_stitched["qcd_pt170to300"], legendEntries["qcd_pt170to300"],
                  nullptr, "",
		  colors, fillStyles, 
		  0.045, 0.62, 0.57, 0.28, 0.39,
		  labelTextLines, 0.050,
		  0.63, 0.66, 0.26, 0.07, 
		  0., 250., "HLT #tau p_{T} Threshold [GeV]", 1.25, 
		  true, 1.e+1, 1.e+7, "Single #tau Trigger Rate [Hz]", 1.40, 
		  outputFileName1);
                  
                TH1* histogram2_minbias_unstitched = histograms_rateDoubleTau_minbias_unstitched[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption]
                  [*absEtaRange][*min_leadTrackPt][*isolationWP];
                string_to_TH1Map1 histograms2_processes_stitched = histograms_rateDoubleTau_vs_processes_stitched[*pfAlgo][*vertexOption][*tauIdOption][*l1MatchingOption]
                  [*absEtaRange][*min_leadTrackPt][*isolationWP];
                std::string outputFileName2 = Form("makeRatePlots_stacked_DoubleTau_%s%s%s_%s_%s_%s_%s.png", 
                  pfAlgo->data(), vertexOption->data(), l1MatchingOption->data(), 
                  tauIdOption->data(), absEtaRange->data(), min_leadTrackPt->data(), isolationWP->data());
                showHistograms(1150, 1600,
                  histogram2_minbias_unstitched,                    "",
                  histograms2_processes_stitched["minbias"],        legendEntries["minbias"],
                  histograms2_processes_stitched["qcd_pt30to50"],   legendEntries["qcd_pt30to50"],
                  histograms2_processes_stitched["qcd_pt50to80"],   legendEntries["qcd_pt50to80"],
                  histograms2_processes_stitched["qcd_pt80to120"],  legendEntries["qcd_pt80to120"],
                  histograms2_processes_stitched["qcd_pt120to170"], legendEntries["qcd_pt120to170"],
                  histograms2_processes_stitched["qcd_pt170to300"], legendEntries["qcd_pt170to300"],
                  nullptr, "",
		  colors, fillStyles, 
		  0.045, 0.62, 0.57, 0.28, 0.39,
		  labelTextLines, 0.050,
		  0.63, 0.66, 0.26, 0.07, 
		  0., 100., "HLT #tau p_{T} Threshold [GeV]", 1.25, 
		  true, 1.e+1, 1.e+6, "Double #tau Trigger Rate [Hz]", 1.40, 
		  outputFileName2);
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

