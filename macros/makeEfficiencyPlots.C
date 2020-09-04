
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
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

void dumpGraph(TGraph* graph)
{
  int numPoints = graph->GetN();
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
    double x, y;
    graph->GetPoint(idxPoint, x, y);
    std::cout << "point #" << idxPoint << ": x = " << x << ", y = " << y << std::endl;
  }
}

TGraph* makeEfficiencyGraph(TH1* histogram_numerator, TH1* histogram_denominator)
{
  //std::cout << "<makeEfficiencyGraph>:" << std::endl;
  assert(histogram_numerator->GetNbinsX() == histogram_denominator->GetNbinsX());
  TAxis* xAxis = histogram_numerator->GetXaxis();
  int numBins = xAxis->GetNbins();
  for ( int idxBin = 1; idxBin <= numBins; ++idxBin ) { 
    double binContent_numerator = histogram_numerator->GetBinContent(idxBin);
    double binContent_denominator = histogram_denominator->GetBinContent(idxBin);
    if( binContent_numerator > binContent_denominator ) {
      std::cerr << "Error in <makeEfficiencyGraph>: numerator = " << binContent_numerator << " exceeds denominator = " << binContent_denominator 
		<< " @ x = " << histogram_denominator->GetBinCenter(idxBin) << " !!" << std::endl;
      assert(0);
    }
    //std::cout << "bin #" << idxBin << " (x=" << xAxis->GetBinLowEdge(idxBin) << ".." << xAxis->GetBinUpEdge(idxBin) << "):" 
    //          << " numerator = " << binContent_numerator << " +/- " << histogram_numerator->GetBinError(idxBin) << ","
    //          << " denominator = " << binContent_denominator << " +/- " << histogram_denominator->GetBinError(idxBin) 
    //          << " (ratio = " << binContent_numerator/binContent_denominator << ")" << std::endl;
  }
  histogram_numerator->GetXaxis()->SetRange(1, numBins);
  histogram_denominator->GetXaxis()->SetRange(1, numBins);
  TString graphName_efficiency = TString(histogram_numerator->GetName()).ReplaceAll("_numerator", "");
  TGraphAsymmErrors* graph_efficiency = new TGraphAsymmErrors(numBins);
  graph_efficiency->Divide(histogram_numerator, histogram_denominator, "w");
  //dumpGraph(graph_efficiency);
  graph_efficiency->SetName(graphName_efficiency.Data());
  return graph_efficiency;
}

/**
 * @brief Integral of Crystal Ball function for fitting trigger efficiency turn-on curves (code from Pascal Paganini)
 * @param m: pT of reconstructed hadronic tau candidates;
 *           the other parameters refer to the shape of the Crystal Ball function, cf. Section 6.2.3 of AN-2016/027
 * @return efficiency for passing trigger, per hadronic tau leg
 */

double
integralCrystalBall(double m,
                    double m0,
                    double sigma,
                    double alpha,
                    double n,
                    double norm)
{
  //std::cout << "<integralCrystalBall>:" << std::endl;
  //std::cout << " m     = " << m << std::endl;
  //std::cout << " m0    = " << m0 << std::endl;
  //std::cout << " sigma = " << sigma << std::endl;
  //std::cout << " alpha = " << alpha << std::endl;
  //std::cout << " n     = " << n << std::endl;
  //std::cout << " norm  = " << norm << std::endl;

  const double sqrtPiOver2 = 1.2533141373;
  const double sqrt2 = 1.4142135624;

  const double sig = std::fabs(static_cast<double>(sigma));

  double t = (m - m0) / sig;
  if ( alpha < 0 )
  {
    t = -t;
  }

  const double absAlpha = std::fabs(alpha / sig);
  const double a = std::pow(n / absAlpha, n) * std::exp(-0.5 * absAlpha * absAlpha);
  const double b = absAlpha - n / absAlpha;

  if ( a >= std::numeric_limits<double>::max() )
  {
    return -1.;
  }

  double ApproxErf;
  double arg = absAlpha / sqrt2;
  if      ( arg >  5. ) ApproxErf =  1;
  else if ( arg < -5. ) ApproxErf = -1;
  else                  ApproxErf = std::erf(arg);

  const double leftArea  = (1 + ApproxErf) * sqrtPiOver2;
  const double rightArea = (a / std::pow(absAlpha - b, n - 1)) / (n - 1);
  const double area = leftArea + rightArea;

  double retVal = 0.;
  if ( t <= absAlpha )
  {
    arg = t / sqrt2;
    if     (arg >  5.) ApproxErf =  1;
    else if(arg < -5.) ApproxErf = -1;
    else               ApproxErf = std::erf(arg);
    retVal = norm * (1 + ApproxErf) * sqrtPiOver2 / area;
  }
  else
  {
    retVal = norm * (leftArea +  a * (1 / std::pow(t-b, n - 1) - 1 / std::pow(absAlpha - b, n - 1)) / (1 - n)) / area;
  }

  //std:: cout << "--> returning " << retVal << std::endl;
  return retVal;
}

Double_t integralCrystalBall_fcn(Double_t* x, Double_t* par)
{
  return integralCrystalBall(x[0], par[0], par[1], par[2], par[3], par[4]);
}

TF1* makeEfficiencyFit(TGraph* graph, const std::string& fitFunctionName, double xMin, double xMax)
{
  TF1* fitFunction = 0;

  int numPoints = graph->GetN();
  if ( numPoints >= 2 )
  {
    TGraph* graph_inverse = new TGraph(numPoints);
    for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) { 
      double x, y;
      graph->GetPoint(idxPoint, x, y);
      graph_inverse->SetPoint(idxPoint, y, x);
    }
    double m0 = graph_inverse->Eval(0.5);
    if ( m0 < 20. ) m0 = 20.;

    fitFunction = new TF1(fitFunctionName.data(), integralCrystalBall_fcn, xMin, xMax, 5);
    fitFunction->SetParameter(0, m0);
    fitFunction->SetParameter(1, 5.);
    fitFunction->SetParameter(2, 0.1);
    fitFunction->SetParameter(3, 2.);
    fitFunction->SetParameter(4, 1.);
    graph->Fit(fitFunction, "QN");

    fitFunction->SetLineColor(graph->GetLineColor());
    fitFunction->SetLineWidth(graph->GetLineWidth());
    fitFunction->SetLineStyle(graph->GetLineStyle());
  
    delete graph_inverse;
  }
  
  return fitFunction;
}

void showGraphs(double canvasSizeX, double canvasSizeY,
		TGraph* graph1, const std::string& legendEntry1,
		TGraph* graph2, const std::string& legendEntry2,
		TGraph* graph3, const std::string& legendEntry3,
		TGraph* graph4, const std::string& legendEntry4,
		TGraph* graph5, const std::string& legendEntry5,
		TGraph* graph6, const std::string& legendEntry6,
		bool addFitFunctions,
		int colors[], int markerStyles[], int lineStyles[], 
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

  if ( !graph1 ) {
    std::cerr << "<showGraphs>: graph1 = NULL --> skipping !!" << std::endl;
    return;
  }

  TH1* dummyHistogram = new TH1F("dummyHistogram", "dummyHistogram", 10, xMin, xMax);
  dummyHistogram->SetTitle("");
  dummyHistogram->SetStats(false);
  dummyHistogram->SetMinimum(yMin);
  dummyHistogram->SetMaximum(yMax);

  TAxis* xAxis = dummyHistogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleSize(0.045);
  xAxis->SetTitleOffset(xAxisOffset);  
  if ( xMax > xMin ) {
    std::cout << "limiting x-axis range to " << xMin << ".." << xMax << std::endl;
    xAxis->SetRangeUser(xMin, xMax);
  }

  TAxis* yAxis = dummyHistogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleSize(0.045);
  yAxis->SetTitleOffset(yAxisOffset);

  dummyHistogram->Draw();
  //dummyHistogram->Draw("axis");
  
  graph1->SetMarkerColor(colors[0]);
  graph1->SetMarkerSize(2);
  graph1->SetMarkerStyle(markerStyles[0]);
  graph1->SetLineColor(colors[0]);
  graph1->SetLineWidth(2);
  graph1->SetLineStyle(lineStyles[0]);
  graph1->Draw("P");
  TF1* fitFunction1 = 0;
  if ( addFitFunctions )
  {
    fitFunction1 = makeEfficiencyFit(graph1, "fitFunction1", xMin, xMax);
    if ( fitFunction1 )
    {
      fitFunction1->SetLineColor(colors[0]);
      fitFunction1->Draw("Lsame");
      graph1->Draw("P");
    }
  }

  TF1* fitFunction2 = 0;
  if ( graph2 ) {
    graph2->SetMarkerColor(colors[1]);
    graph2->SetMarkerSize(2);
    graph2->SetMarkerStyle(markerStyles[1]);
    graph2->SetLineColor(colors[1]);
    graph2->SetLineWidth(2);
    graph2->SetLineStyle(lineStyles[1]);
    graph2->Draw("P");
    if ( addFitFunctions )
    {
      fitFunction2 = makeEfficiencyFit(graph2, "fitFunction2", xMin, xMax);
      if ( fitFunction2 )
      {
	fitFunction2->SetLineColor(colors[1]);
        fitFunction2->Draw("Lsame");
        graph2->Draw("P");
      }
    }
  }

  TF1* fitFunction3 = 0;
  if ( graph3 ) {
    graph3->SetMarkerColor(colors[2]);
    graph3->SetMarkerSize(2);
    graph3->SetMarkerStyle(markerStyles[2]);
    graph3->SetLineColor(colors[2]);
    graph3->SetLineWidth(2);
    graph3->SetLineStyle(lineStyles[2]);
    graph3->Draw("P");
    if ( addFitFunctions )
    {
      fitFunction3 = makeEfficiencyFit(graph3, "fitFunction3", xMin, xMax);
      if ( fitFunction3 )
      {
	fitFunction3->SetLineColor(colors[2]);
        fitFunction3->Draw("Lsame");
        graph3->Draw("P");
      }
    }
  }

  TF1* fitFunction4 = 0;
  if ( graph4 ) {
    graph4->SetMarkerColor(colors[3]);
    graph4->SetMarkerSize(2);
    graph4->SetMarkerStyle(markerStyles[3]);
    graph4->SetLineColor(colors[3]);
    graph4->SetLineWidth(2);
    graph4->SetLineStyle(lineStyles[3]);
    graph4->Draw("P");
    if ( addFitFunctions )
    {
      fitFunction4 = makeEfficiencyFit(graph4, "fitFunction4", xMin, xMax);
      if ( fitFunction4 )
      {
	fitFunction4->SetLineColor(colors[3]);
	fitFunction4->Draw("Lsame");
	graph4->Draw("P");
      }
    }
  }

  TF1* fitFunction5 = 0;
  if ( graph5 ) {
    graph5->SetMarkerColor(colors[4]);
    graph5->SetMarkerSize(2);
    graph5->SetMarkerStyle(markerStyles[4]);
    graph5->SetLineColor(colors[4]);
    graph5->SetLineWidth(2);
    graph5->SetLineStyle(lineStyles[4]);
    graph5->Draw("P");
    if ( addFitFunctions )
    {
      fitFunction5 = makeEfficiencyFit(graph5, "fitFunction5", xMin, xMax);
      if ( fitFunction5 )
      {
	fitFunction5->SetLineColor(colors[4]);
	fitFunction5->Draw("Lsame");
	graph5->Draw("P");
      }
    }
  }

  TF1* fitFunction6 = 0;
  if ( graph6 ) {
    graph6->SetMarkerColor(colors[5]);
    graph6->SetMarkerSize(2);
    graph6->SetMarkerStyle(markerStyles[5]);
    graph6->SetLineColor(colors[5]);
    graph6->SetLineWidth(2);
    graph6->SetLineStyle(lineStyles[5]);
    graph6->Draw("P");
    if ( addFitFunctions )
    {
      fitFunction6 = makeEfficiencyFit(graph6, "fitFunction6", xMin, xMax);
      if ( fitFunction6 )
      {
	fitFunction6->SetLineColor(colors[5]);
	fitFunction6->Draw("Lsame");
	graph6->Draw("P");
      }
    }
  }

  TLegend* legend = 0;
  if ( legendEntry1 != "" ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, "", "brNDC"); 
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetTextSize(legendTextSize);
    legend->AddEntry(graph1, legendEntry1.data(), "p");
    if ( graph2 ) legend->AddEntry(graph2, legendEntry2.data(), "p");
    if ( graph3 ) legend->AddEntry(graph3, legendEntry3.data(), "p");
    if ( graph4 ) legend->AddEntry(graph4, legendEntry4.data(), "p");
    if ( graph5 ) legend->AddEntry(graph5, legendEntry5.data(), "p");
    if ( graph6 ) legend->AddEntry(graph6, legendEntry6.data(), "p");
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

  dummyHistogram->Draw("axissame");

  canvas->RedrawAxis();

  canvas->Update();
  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  //canvas->Print(std::string(outputFileName_plot).append(".png").data());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").data());

  delete dummyHistogram;
  delete fitFunction1;
  delete fitFunction2;
  delete fitFunction3;
  delete fitFunction4;
  delete fitFunction5;
  delete fitFunction6;
  delete label;
  delete legend;
  delete canvas;  
}

TGraph* getGraph(std::map<std::string, TGraph*>& graphs, const std::string& isolationWP)
{
  TGraph* graph = nullptr;
  if ( isolationWP != "" && graphs.find(isolationWP) != graphs.end() ) 
  {
    graph = graphs.find(isolationWP)->second;
  }
  return graph;
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

std::vector<std::string> getLabelTextLines(const std::string& ptThreshold)
{
  std::vector<std::string> labelTextLines;
  if ( ptThreshold == "ptGt20" ) labelTextLines.push_back("p_{T} > 20 GeV");
  if ( ptThreshold == "ptGt30" ) labelTextLines.push_back("p_{T} > 30 GeV");
  if ( ptThreshold == "ptGt40" ) labelTextLines.push_back("p_{T} > 40 GeV");
  else assert(0);
  return labelTextLines;
}

void makeEfficiencyPlots()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  std::string inputFilePath = "/hdfs/local/veelken/Phase2HLT/efficiency/2020Sep03/";
  std::string inputFileName = "hadd_qqH_htt_all.root";
  std::string inputFileName_full = inputFilePath;
  if ( inputFileName_full.find_last_of("/") != (inputFileName_full.size() - 1) ) inputFileName_full.append("/");
  inputFileName_full.append(inputFileName);
  TFile* inputFile = new TFile(inputFileName_full.data());
  if ( !inputFile ) {
    std::cerr << "Failed to open input file = '" << inputFileName_full << "' !!" << std::endl;
    assert(0);
  }

  std::string l1_pfAlgo = "WithStripsPF"; // 'WithStripsPF' or 'WithoutStripsPF'
  std::string l1_ptThreshold = "ptGt20";
  std::string l1_min_leadTrackPt = "leadTrackPtGt5";
  std::string l1_isolationWP = "relChargedIsoLt0p05";

  std::vector<std::string> hlt_pfAlgos;
  //hlt_pfAlgos.push_back("PFTau");
  hlt_pfAlgos.push_back("HpsPFTau");

  std::vector<std::string> hlt_vertexOptions;
  hlt_vertexOptions.push_back("8HitsMaxDeltaZWithOfflineVertices");
  //hlt_vertexOptions.push_back("8HitsMaxDeltaZToLeadTrackWithOfflineVertices");
  //hlt_vertexOptions.push_back("8HitsMaxDeltaZWithOnlineVertices");
  //hlt_vertexOptions.push_back("8HitsMaxDeltaZToLeadTrackWithOnlineVertices");
  //hlt_vertexOptions.push_back("8HitsMaxDeltaZWithOnlineVerticesTrimmed");
  //hlt_vertexOptions.push_back("8HitsMaxDeltaZToLeadTrackWithOnlineVerticesTrimmed");
/*
  hlt_vertexOptions.push_back("5HitsMaxDeltaZWithOfflineVertices");
  hlt_vertexOptions.push_back("5HitsMaxDeltaZToLeadTrackWithOfflineVertices");
  hlt_vertexOptions.push_back("5HitsMaxDeltaZWithOnlineVertices");
  hlt_vertexOptions.push_back("5HitsMaxDeltaZToLeadTrackWithOnlineVertices");
  //hlt_vertexOptions.push_back("5HitsMaxDeltaZWithOnlineVerticesTrimmed");
  //hlt_vertexOptions.push_back("5HitsMaxDeltaZToLeadTrackWithOnlineVerticesTrimmed");
  hlt_vertexOptions.push_back("3HitsMaxDeltaZWithOfflineVertices");
  hlt_vertexOptions.push_back("3HitsMaxDeltaZToLeadTrackWithOfflineVertices");
  hlt_vertexOptions.push_back("3HitsMaxDeltaZWithOnlineVertices");
  hlt_vertexOptions.push_back("3HitsMaxDeltaZToLeadTrackWithOnlineVertices");
  //hlt_vertexOptions.push_back("3HitsMaxDeltaZWithOnlineVerticesTrimmed");
  //hlt_vertexOptions.push_back("3HitsMaxDeltaZToLeadTrackWithOnlineVerticesTrimmed");
 */

  std::map<std::string, std::string> hlt_srcVertices; // key = vertexOption
  hlt_srcVertices["8HitsMaxDeltaZWithOfflineVertices"]                  = "offlinePrimaryVertices";
  hlt_srcVertices["8HitsMaxDeltaZToLeadTrackWithOfflineVertices"]       = "offlinePrimaryVertices";
  hlt_srcVertices["8HitsMaxDeltaZWithOnlineVertices"]                   = "hltPhase2PixelVertices";
  hlt_srcVertices["8HitsMaxDeltaZToLeadTrackWithOnlineVertices"]        = "hltPhase2PixelVertices";
  hlt_srcVertices["8HitsMaxDeltaZWithOnlineVerticesTrimmed"]            = "hltPhase2TrimmedPixelVertices";
  hlt_srcVertices["8HitsMaxDeltaZToLeadTrackWithOnlineVerticesTrimmed"] = "hltPhase2TrimmedPixelVertices";
  hlt_srcVertices["5HitsMaxDeltaZWithOfflineVertices"]                  = "offlinePrimaryVertices";
  hlt_srcVertices["5HitsMaxDeltaZToLeadTrackWithOfflineVertices"]       = "offlinePrimaryVertices";
  hlt_srcVertices["5HitsMaxDeltaZWithOnlineVertices"]                   = "hltPhase2PixelVertices";
  hlt_srcVertices["5HitsMaxDeltaZToLeadTrackWithOnlineVertices"]        = "hltPhase2PixelVertices";
  hlt_srcVertices["5HitsMaxDeltaZWithOnlineVerticesTrimmed"]            = "hltPhase2TrimmedPixelVertices";
  hlt_srcVertices["5HitsMaxDeltaZToLeadTrackWithOnlineVerticesTrimmed"] = "hltPhase2TrimmedPixelVertices";
  hlt_srcVertices["3HitsMaxDeltaZWithOfflineVertices"]                  = "offlinePrimaryVertices";
  hlt_srcVertices["3HitsMaxDeltaZToLeadTrackWithOfflineVertices"]       = "offlinePrimaryVertices";
  hlt_srcVertices["3HitsMaxDeltaZWithOnlineVertices"]                   = "hltPhase2PixelVertices";
  hlt_srcVertices["3HitsMaxDeltaZToLeadTrackWithOnlineVertices"]        = "hltPhase2PixelVertices";
  hlt_srcVertices["3HitsMaxDeltaZWithOnlineVerticesTrimmed"]            = "hltPhase2TrimmedPixelVertices";
  hlt_srcVertices["3HitsMaxDeltaZToLeadTrackWithOnlineVerticesTrimmed"] = "hltPhase2TrimmedPixelVertices";

  std::vector<std::string> hlt_tauIdOptions;
  hlt_tauIdOptions.push_back("recoSumChargedIso");
  hlt_tauIdOptions.push_back("patSumChargedIso");
  hlt_tauIdOptions.push_back("patDeepTau");

  std::vector<std::string> l1MatchingOptions;
  l1MatchingOptions.push_back("");            // CV: no matching of HLT taus to L1 taus
  //l1MatchingOptions.push_back("MatchedToL1"); 

  std::vector<std::string> observables;
  observables.push_back("pt");
  observables.push_back("eta");
  //observables.push_back("phi");
  //observables.push_back("minDeltaR");

  std::vector<std::string> absEtaRanges;
  //absEtaRanges.push_back("absEtaLt1p40");
  //absEtaRanges.push_back("absEta1p40to2p17");
  //absEtaRanges.push_back("absEta1p40to2p40");
  //absEtaRanges.push_back("absEtaLt2p17");
  absEtaRanges.push_back("absEtaLt2p40");

  std::vector<std::string> decayModes;
  //decayModes.push_back("oneProng0Pi0");
  //decayModes.push_back("oneProng1Pi0");
  //decayModes.push_back("oneProng2Pi0");
  //decayModes.push_back("threeProng0Pi0");
  //decayModes.push_back("threeProng1Pi0");
  decayModes.push_back("all");

  std::vector<std::string> hlt_ptThresholds;
  //hlt_ptThresholds.push_back("pt_numeratorGt20");
  //hlt_ptThresholds.push_back("pt_numeratorGt25");
  hlt_ptThresholds.push_back("pt_numeratorGt30");
  //hlt_ptThresholds.push_back("pt_numeratorGt35");
  //hlt_ptThresholds.push_back("pt_numeratorGt40");
  //hlt_ptThresholds.push_back("pt_numeratorGt45");
  //hlt_ptThresholds.push_back("pt_numeratorGt50");

  std::vector<std::string> hlt_min_leadTrackPtValues;
  //hlt_min_leadTrackPtValues.push_back("leadTrackPtGt1");
  //hlt_min_leadTrackPtValues.push_back("leadTrackPtGt2");
  hlt_min_leadTrackPtValues.push_back("leadTrackPtGt5");

  std::map<std::string, std::vector<std::string>> hlt_isolationWPs; // key = hlt_tauIdOption
  hlt_isolationWPs["recoSumChargedIso"].push_back("noIsolation");
  hlt_isolationWPs["recoSumChargedIso"].push_back("relDiscriminatorLt0p400");
  hlt_isolationWPs["recoSumChargedIso"].push_back("relDiscriminatorLt0p200");
  hlt_isolationWPs["recoSumChargedIso"].push_back("relDiscriminatorLt0p100");
  hlt_isolationWPs["recoSumChargedIso"].push_back("relDiscriminatorLt0p050");
  hlt_isolationWPs["patSumChargedIso"] = hlt_isolationWPs["recoSumChargedIso"];
  hlt_isolationWPs["patDeepTau"].push_back("absDiscriminatorGt0p260");
  hlt_isolationWPs["patDeepTau"].push_back("absDiscriminatorGt0p425");
  hlt_isolationWPs["patDeepTau"].push_back("absDiscriminatorGt0p598");
  hlt_isolationWPs["patDeepTau"].push_back("absDiscriminatorGt0p785");
  hlt_isolationWPs["patDeepTau"].push_back("absDiscriminatorGt0p883");
  hlt_isolationWPs["patDeepTau"].push_back("absDiscriminatorGt0p931");

  std::map<std::string, double> xMin; // key = observable
  xMin["pt"]        =   0.;
  xMin["eta"]       =  -3.0;
  xMin["phi"]       = -TMath::Pi();
  xMin["minDeltaR"] =   0.;
  
  std::map<std::string, double> xMax; // key = observable
  xMax["pt"]        = 100.;
  xMax["eta"]       =  +3.0;
  xMax["phi"]       = +TMath::Pi();
  xMax["minDeltaR"] =   5.;

  std::map<std::string, std::string> xAxisTitles; // key = observable
  //xAxisTitles["pt"]        = "True #tau_{h} p_{T} [GeV]";
  //xAxisTitles["eta"]       = "True #tau_{h} #eta";
  //xAxisTitles["phi"]       = "True #tau_{h} #phi";
  xAxisTitles["pt"]        = "Offline #tau_{h} p_{T} [GeV]";
  xAxisTitles["eta"]       = "Offline #tau_{h} #eta";
  xAxisTitles["phi"]       = "Offline #tau_{h} #phi";
  xAxisTitles["minDeltaR"] = "#Delta R";
  
  std::map<std::string, std::map<std::string, std::string>> legendEntries_vs_isolationWPs; // key = hlt_tauIdOption, hlt_isolationWP
  legendEntries_vs_isolationWPs["recoSumChargedIso"]["noIsolation"]             = "No Isolation";
  legendEntries_vs_isolationWPs["recoSumChargedIso"]["relDiscriminatorLt0p400"] = "I_{ch} < 0.40*p_{T}";
  legendEntries_vs_isolationWPs["recoSumChargedIso"]["relDiscriminatorLt0p200"] = "I_{ch} < 0.20*p_{T}";
  legendEntries_vs_isolationWPs["recoSumChargedIso"]["relDiscriminatorLt0p100"] = "I_{ch} < 0.10*p_{T}";
  legendEntries_vs_isolationWPs["recoSumChargedIso"]["relDiscriminatorLt0p050"] = "I_{ch} < 0.05*p_{T}";
  legendEntries_vs_isolationWPs["patSumChargedIso"] = legendEntries_vs_isolationWPs["recoSumChargedIso"];
  legendEntries_vs_isolationWPs["patDeepTau"]["absDiscriminatorGt0p260"]        = "D > 0.260";
  legendEntries_vs_isolationWPs["patDeepTau"]["absDiscriminatorGt0p425"]        = "D > 0.425";
  legendEntries_vs_isolationWPs["patDeepTau"]["absDiscriminatorGt0p598"]        = "D > 0.598";
  legendEntries_vs_isolationWPs["patDeepTau"]["absDiscriminatorGt0p785"]        = "D > 0.785";
  legendEntries_vs_isolationWPs["patDeepTau"]["absDiscriminatorGt0p883"]        = "D > 0.883";
  legendEntries_vs_isolationWPs["patDeepTau"]["absDiscriminatorGt0p931"]        = "D > 0.931";

  std::map<std::string, std::string> legendEntries_vs_leadTrackPt; // key = min_leadTrackPt
  legendEntries_vs_leadTrackPt["leadTrackPtGt1"] = "lead. Track p_{T} > 1 GeV";
  legendEntries_vs_leadTrackPt["leadTrackPtGt2"] = "lead. Track p_{T} > 2 GeV";
  legendEntries_vs_leadTrackPt["leadTrackPtGt5"] = "lead. Track p_{T} > 5 GeV";

  std::map<std::string, std::string> legendEntries_vs_decayModes; // key = decayMode
  legendEntries_vs_decayModes["oneProng0Pi0"]   = "h^{#pm}";
  legendEntries_vs_decayModes["oneProng1Pi0"]   = "h^{#pm}#pi^{0}";
  legendEntries_vs_decayModes["oneProng2Pi0"]   = "h^{#pm}#pi^{0}#pi^{0}";
  legendEntries_vs_decayModes["threeProng0Pi0"] = "h^{#pm}h^{#mp}h^{#pm}";
  legendEntries_vs_decayModes["threeProng1Pi0"] = "h^{#pm}h^{#mp}h^{#pm}#pi^{0}";
  legendEntries_vs_decayModes["all"]            = "all";

  std::string hlt_dqmDirectory = "%sAnalyzerSignal%s%s_%s";
  std::string l1_dqmDirectory = "L1HPSPFTauAnalyzerSignal"; 

  int colors[6]       = {  1,  2,  8,  4,  6,  7 };
  int lineStyles[6]   = {  1,  1,  1,  1,  1,  1 };
  int markerStyles[6] = { 22, 32, 20, 24, 21, 25 };

  typedef std::map<std::string, TGraph*>              string_to_TGraphMap1;
  typedef std::map<std::string, string_to_TGraphMap1> string_to_TGraphMap2;
  typedef std::map<std::string, string_to_TGraphMap2> string_to_TGraphMap3;
  typedef std::map<std::string, string_to_TGraphMap3> string_to_TGraphMap4;
  typedef std::map<std::string, string_to_TGraphMap4> string_to_TGraphMap5;
  typedef std::map<std::string, string_to_TGraphMap5> string_to_TGraphMap6;
  typedef std::map<std::string, string_to_TGraphMap6> string_to_TGraphMap7;
  typedef std::map<std::string, string_to_TGraphMap7> string_to_TGraphMap8;
  typedef std::map<std::string, string_to_TGraphMap8> string_to_TGraphMap9;
  typedef std::map<std::string, string_to_TGraphMap9> string_to_TGraphMap10;
 
  // CV: compute L1 efficiency
  string_to_TGraphMap2 graphs_l1_efficiency; // key = observable, absEtaRange
  for ( std::vector<std::string>::const_iterator observable = observables.begin();
	observable != observables.end(); ++observable ) {      
    for ( std::vector<std::string>::const_iterator absEtaRange = absEtaRanges.begin();
	  absEtaRange != absEtaRanges.end(); ++absEtaRange ) {
      std::string histogramName_numerator = Form("%s%s_wrtGenHadTaus/%s/%s/effL1PFTau_vs_%s_numerator_all_%s_%s_%s", 
        l1_dqmDirectory.data(), l1_pfAlgo.data(), absEtaRange->data(), "all", 
        observable->data(), absEtaRange->data(), l1_ptThreshold.data(), l1_isolationWP.data());
      //std::string histogramName_numerator = Form("%s%s_wrtOfflineTaus/%s/%s/effL1PFTau_vs_%s_numerator_all_%s_%s_%s", 
      //  l1_dqmDirectory.data(), l1_pfAlgo.data(), absEtaRange->data(), "all", 
      //  observable->data(), absEtaRange->data(), l1_ptThreshold.data(), l1_isolationWP.data());
      TH1* histogram_numerator = loadHistogram(inputFile, histogramName_numerator);
      std::string histogramName_denominator = Form("%s%s_wrtGenHadTaus/%s/%s/effL1PFTau_vs_%s_denominator_all_%s_%s_%s", 
        l1_dqmDirectory.data(), l1_pfAlgo.data(), absEtaRange->data(), "all", 
        observable->data(), absEtaRange->data(), l1_ptThreshold.data(), l1_isolationWP.data());
      //std::string histogramName_denominator = Form("%s%s_wrtOfflineTaus/%s/%s/effL1PFTau_vs_%s_denominator_all_%s_%s_%s", 
      //  l1_dqmDirectory.data(), l1_pfAlgo.data(), absEtaRange->data(), "all", 
      //  observable->data(), absEtaRange->data(), l1_ptThreshold.data(), l1_isolationWP.data());
      TH1* histogram_denominator = loadHistogram(inputFile, histogramName_denominator);
      TGraph* graph_efficiency = makeEfficiencyGraph(histogram_numerator, histogram_denominator);
      graphs_l1_efficiency[*observable][*absEtaRange] = graph_efficiency;
    } // absEtaRange
  } // observable

  // CV: compute HLT efficiency
  string_to_TGraphMap9  graphs_hlt_efficiency_vs_isolationWPs; // key = hlt_pfAlgo, hlt_vertexOption, hlt_tauIdOption, l1MatchingOption, 
                                                               //       observable, absEtaRange, hlt_ptThreshold, hlt_min_leadTrackPt, hlt_isolationWP
  string_to_TGraphMap9  graphs_hlt_efficiency_vs_leadTrackPt;  // key = hlt_pfAlgo, hlt_vertexOption, hlt_tauIdOption, l1MatchingOption, 
                                                               //       observable, absEtaRange, hlt_ptThreshold, hlt_isolationWP, hlt_min_leadTrackPt
  string_to_TGraphMap10 graphs_hlt_efficiency_vs_decayModes;   // key = hlt_pfAlgo, hlt_vertexOption, hlt_tauIdOption, l1MatchingOption, 
                                                               //       observable, absEtaRange, hlt_ptThreshold, hlt_min_leadTrackPt, hlt_isolationWP, decayMode
  for ( std::vector<std::string>::const_iterator hlt_pfAlgo = hlt_pfAlgos.begin();
	hlt_pfAlgo != hlt_pfAlgos.end(); ++hlt_pfAlgo ) {
    for ( std::vector<std::string>::const_iterator hlt_vertexOption = hlt_vertexOptions.begin();
	  hlt_vertexOption != hlt_vertexOptions.end(); ++hlt_vertexOption ) {
      for ( std::vector<std::string>::const_iterator hlt_tauIdOption = hlt_tauIdOptions.begin();
            hlt_tauIdOption != hlt_tauIdOptions.end(); ++hlt_tauIdOption ) {
        for ( std::vector<std::string>::const_iterator observable = observables.begin();
	      observable != observables.end(); ++observable ) {      
          for ( std::vector<std::string>::const_iterator absEtaRange = absEtaRanges.begin();
	        absEtaRange != absEtaRanges.end(); ++absEtaRange ) {
            for ( std::vector<std::string>::const_iterator hlt_ptThreshold = hlt_ptThresholds.begin();
	          hlt_ptThreshold != hlt_ptThresholds.end(); ++hlt_ptThreshold ) {      
              for ( std::vector<std::string>::const_iterator l1MatchingOption = l1MatchingOptions.begin();
	            l1MatchingOption != l1MatchingOptions.end(); ++l1MatchingOption ) { 
                for ( std::vector<std::string>::const_iterator hlt_min_leadTrackPt = hlt_min_leadTrackPtValues.begin();
                      hlt_min_leadTrackPt != hlt_min_leadTrackPtValues.end(); ++hlt_min_leadTrackPt ) {
                  for ( std::vector<std::string>::const_iterator hlt_isolationWP = hlt_isolationWPs[*hlt_tauIdOption].begin();
	                hlt_isolationWP != hlt_isolationWPs[*hlt_tauIdOption].end(); ++hlt_isolationWP ) {
                    std::string histogramName_numerator = Form(("%s/" + hlt_dqmDirectory + "_wrtGenHadTaus/%s/all/effPFTau_vs_%s_numerator_all_%s_%s_%s_%s").data(), 
                      hlt_srcVertices[*hlt_vertexOption].data(), hlt_pfAlgo->data(), hlt_vertexOption->data(), l1MatchingOption->data(), hlt_tauIdOption->data(), absEtaRange->data(),
                      observable->data(), hlt_ptThreshold->data(), absEtaRange->data(), hlt_min_leadTrackPt->data(), hlt_isolationWP->data());
	            //std::string histogramName_numerator = Form(("%s/" + hlt_dqmDirectory + "_wrtOfflineTaus/%s/all/effPFTau_vs_%s_numerator_all_%s_%s_%s_%s").data(), 
                    //  hlt_srcVertices[*hlt_vertexOption].data(), hlt_pfAlgo->data(), hlt_vertexOption->data(), l1MatchingOption->data(), hlt_tauIdOption->data(), absEtaRange->data(),
                    //  observable->data(), hlt_ptThreshold->data(), absEtaRange->data(), hlt_min_leadTrackPt->data(), hlt_isolationWP->data());
                    TH1* histogram_numerator = loadHistogram(inputFile, histogramName_numerator);
	            std::string histogramName_denominator = Form(("%s/" + hlt_dqmDirectory + "_wrtGenHadTaus/%s/all/effPFTau_vs_%s_denominator_all_%s_%s_%s_%s").data(), 
                      hlt_srcVertices[*hlt_vertexOption].data(), hlt_pfAlgo->data(), hlt_vertexOption->data(), l1MatchingOption->data(), hlt_tauIdOption->data(), absEtaRange->data(),
                      observable->data(), hlt_ptThreshold->data(), absEtaRange->data(), hlt_min_leadTrackPt->data(), hlt_isolationWP->data());
	            //std::string histogramName_denominator = Form(("%s/" + hlt_dqmDirectory + "_wrtOfflineTaus/%s/all/effPFTau_vs_%s_denominator_all_%s_%s_%s_%s").data(), 
                    //  hlt_srcVertices[*hlt_vertexOption].data(), hlt_pfAlgo->data(), hlt_vertexOption->data(), l1MatchingOption->data(), hlt_tauIdOption->data(), absEtaRange->data(),
                    //  observable->data(), hlt_ptThreshold->data(), absEtaRange->data(), hlt_min_leadTrackPt->data(), hlt_isolationWP->data());
                    TH1* histogram_denominator = loadHistogram(inputFile, histogramName_denominator);
	            TGraph* graph_efficiency = makeEfficiencyGraph(histogram_numerator, histogram_denominator);
	            graphs_hlt_efficiency_vs_isolationWPs[*hlt_pfAlgo][*hlt_vertexOption][*hlt_tauIdOption][*l1MatchingOption]
                      [*observable][*absEtaRange][*hlt_ptThreshold][*hlt_min_leadTrackPt][*hlt_isolationWP] = graph_efficiency;
                    graphs_hlt_efficiency_vs_leadTrackPt[*hlt_pfAlgo][*hlt_vertexOption][*hlt_tauIdOption][*l1MatchingOption]
                      [*observable][*absEtaRange][*hlt_ptThreshold][*hlt_isolationWP][*hlt_min_leadTrackPt] = graph_efficiency;
                    for ( std::vector<std::string>::const_iterator decayMode = decayModes.begin();
	  	          decayMode != decayModes.end(); ++decayMode ) {
                      std::string histogramName_numerator = Form(("%s/" + hlt_dqmDirectory + "_wrtGenHadTaus/%s/%s/effPFTau_vs_%s_numerator_%s_%s_%s_%s_%s").data(), 
                        hlt_srcVertices[*hlt_vertexOption].data(), 
                        hlt_pfAlgo->data(), hlt_vertexOption->data(), l1MatchingOption->data(), hlt_tauIdOption->data(), 
                        absEtaRange->data(), decayMode->data(),
                        observable->data(), decayMode->data(), hlt_ptThreshold->data(), absEtaRange->data(), hlt_min_leadTrackPt->data(), hlt_isolationWP->data());
	              //std::string histogramName_numerator = Form(("%s/" + hlt_dqmDirectory + "_wrtOfflineTaus/%s/%s/effPFTau_vs_%s_numerator_%s_%s_%s_%s_%s").data(), 
                      //  hlt_srcVertices[*hlt_vertexOption].data(), 
                      //  hlt_pfAlgo->data(), hlt_vertexOption->data(), l1MatchingOption->data(), hlt_tauIdOption->data(), 
                      //  absEtaRange->data(), decayMode->data(),
                      //  observable->data(), decayMode->data(), hlt_ptThreshold->data(), absEtaRange->data(), hlt_min_leadTrackPt->data(), hlt_isolationWP->data());
                      TH1* histogram_numerator = loadHistogram(inputFile, histogramName_numerator);
	              std::string histogramName_denominator = Form(("%s/" + hlt_dqmDirectory + "_wrtGenHadTaus/%s/%s/effPFTau_vs_%s_denominator_%s_%s_%s_%s_%s").data(), 
                        hlt_srcVertices[*hlt_vertexOption].data(), 
                        hlt_pfAlgo->data(), hlt_vertexOption->data(), l1MatchingOption->data(), hlt_tauIdOption->data(), 
                        absEtaRange->data(), decayMode->data(),
                        observable->data(), decayMode->data(), hlt_ptThreshold->data(), absEtaRange->data(), hlt_min_leadTrackPt->data(), hlt_isolationWP->data());
	              //std::string histogramName_denominator = Form(("%s/" + hlt_dqmDirectory + "_wrtOfflineTaus/%s/%s/effPFTau_vs_%s_denominator_%s_%s_%s_%s_%s").data(), 
                      //  hlt_srcVertices[*hlt_vertexOption].data(), 
                      //  hlt_pfAlgo->data(), hlt_vertexOption->data(), l1MatchingOption->data(), hlt_tauIdOption->data(), 
                      //  absEtaRange->data(), decayMode->data(),
                      //  observable->data(), decayMode->data(), hlt_ptThreshold->data(), absEtaRange->data(), hlt_min_leadTrackPt->data(), hlt_isolationWP->data());
                      TH1* histogram_denominator = loadHistogram(inputFile, histogramName_denominator);
	              TGraph* graph_efficiency = makeEfficiencyGraph(histogram_numerator, histogram_denominator);
	              graphs_hlt_efficiency_vs_decayModes[*hlt_pfAlgo][*hlt_vertexOption][*hlt_tauIdOption][*l1MatchingOption]
                        [*observable][*absEtaRange][*hlt_ptThreshold][*hlt_min_leadTrackPt][*hlt_isolationWP][*decayMode] = graph_efficiency;
                    } // decayMode
                  } // hlt_isolationWP
                } // hlt_min_leadTrackPt
              } // l1MatchingOption
            } // hlt_ptThreshold
          } // absEtaRange
        } // observable
      } // hlt_tauIdOption
    } // hlt_vertexOption
  } // hlt_pfAlgo

  // CV: make plots
  for ( std::vector<std::string>::const_iterator hlt_pfAlgo = hlt_pfAlgos.begin();
	hlt_pfAlgo != hlt_pfAlgos.end(); ++hlt_pfAlgo ) {
    for ( std::vector<std::string>::const_iterator hlt_vertexOption = hlt_vertexOptions.begin();
	  hlt_vertexOption != hlt_vertexOptions.end(); ++hlt_vertexOption ) {
      for ( std::vector<std::string>::const_iterator hlt_tauIdOption = hlt_tauIdOptions.begin();
            hlt_tauIdOption != hlt_tauIdOptions.end(); ++hlt_tauIdOption ) {
        for ( std::vector<std::string>::const_iterator observable = observables.begin();
	      observable != observables.end(); ++observable ) {      
          for ( std::vector<std::string>::const_iterator absEtaRange = absEtaRanges.begin();
	        absEtaRange != absEtaRanges.end(); ++absEtaRange ) {
            for ( std::vector<std::string>::const_iterator hlt_ptThreshold = hlt_ptThresholds.begin();
  	          hlt_ptThreshold != hlt_ptThresholds.end(); ++hlt_ptThreshold ) {      
              for ( std::vector<std::string>::const_iterator l1MatchingOption = l1MatchingOptions.begin();
	            l1MatchingOption != l1MatchingOptions.end(); ++l1MatchingOption ) {
                for ( std::vector<std::string>::const_iterator hlt_min_leadTrackPt = hlt_min_leadTrackPtValues.begin();
                      hlt_min_leadTrackPt != hlt_min_leadTrackPtValues.end(); ++hlt_min_leadTrackPt ) {

                  std::string hlt_isolationWP1 = ( hlt_isolationWPs[*hlt_tauIdOption].size() >= 1 ) ? hlt_isolationWPs[*hlt_tauIdOption][0] : "";
                  std::string hlt_isolationWP2 = ( hlt_isolationWPs[*hlt_tauIdOption].size() >= 2 ) ? hlt_isolationWPs[*hlt_tauIdOption][1] : "";
                  std::string hlt_isolationWP3 = ( hlt_isolationWPs[*hlt_tauIdOption].size() >= 3 ) ? hlt_isolationWPs[*hlt_tauIdOption][2] : "";
                  std::string hlt_isolationWP4 = ( hlt_isolationWPs[*hlt_tauIdOption].size() >= 4 ) ? hlt_isolationWPs[*hlt_tauIdOption][3] : "";
                  std::string hlt_isolationWP5 = ( hlt_isolationWPs[*hlt_tauIdOption].size() >= 5 ) ? hlt_isolationWPs[*hlt_tauIdOption][4] : "";
                  std::string hlt_isolationWP6 = ( hlt_isolationWPs[*hlt_tauIdOption].size() >= 6 ) ? hlt_isolationWPs[*hlt_tauIdOption][5] : "";

  	          string_to_TGraphMap1 graphs1 = graphs_hlt_efficiency_vs_isolationWPs[*hlt_pfAlgo][*hlt_vertexOption][*hlt_tauIdOption][*l1MatchingOption]
                    [*observable][*absEtaRange][*hlt_ptThreshold][*hlt_min_leadTrackPt];                 
	          bool addFitFunctions1 = false;
                  double legendPosX = 0.71;
	          if ( (*observable) == "pt" ) 
	          {
	            addFitFunctions1 = true;
	          }
                  if ( (*observable) == "eta" ) 
	          {
                    legendPosX = 0.43;
	          }
	          std::vector<std::string> labelTextLines1 = getLabelTextLines(*hlt_ptThreshold);
                  std::string outputFileName1 = Form("makeEfficiencyPlots_%s%s%s_%s_%s_%s_%s_%s_vs_isolationWP.png", 
                    hlt_pfAlgo->data(), hlt_vertexOption->data(), l1MatchingOption->data(), 
                    hlt_tauIdOption->data(), observable->data(), absEtaRange->data(), hlt_ptThreshold->data(), hlt_min_leadTrackPt->data());
	          showGraphs(1150, 850,
                             getGraph(graphs1, hlt_isolationWP1), getLegendEntry(legendEntries_vs_isolationWPs[*hlt_tauIdOption], hlt_isolationWP1),
                             getGraph(graphs1, hlt_isolationWP2), getLegendEntry(legendEntries_vs_isolationWPs[*hlt_tauIdOption], hlt_isolationWP2),
                             getGraph(graphs1, hlt_isolationWP3), getLegendEntry(legendEntries_vs_isolationWPs[*hlt_tauIdOption], hlt_isolationWP3),
                             getGraph(graphs1, hlt_isolationWP4), getLegendEntry(legendEntries_vs_isolationWPs[*hlt_tauIdOption], hlt_isolationWP4),
                             getGraph(graphs1, hlt_isolationWP5), getLegendEntry(legendEntries_vs_isolationWPs[*hlt_tauIdOption], hlt_isolationWP5),
                             getGraph(graphs1, hlt_isolationWP6), getLegendEntry(legendEntries_vs_isolationWPs[*hlt_tauIdOption], hlt_isolationWP6),
  		             addFitFunctions1,
		             colors, markerStyles, lineStyles, 
		             0.040, legendPosX, 0.17, 0.23, 0.28, 
		             labelTextLines1, 0.050,
		             0.17, 0.85, 0.26, 0.05, 
		             xMin[*observable], xMax[*observable], xAxisTitles[*observable], 1.2, 
		             false, 0., 1.09, "Efficiency", 1.4, 
		             outputFileName1);
                } // hlt_min_leadTrackPt
              } // l1MatchingOption
/*
              for ( std::vector<std::string>::const_iterator l1MatchingOption = l1MatchingOptions.begin();
	            l1MatchingOption != l1MatchingOptions.end(); ++l1MatchingOption ) {
                for ( std::vector<std::string>::const_iterator hlt_isolationWP = hlt_isolationWPs[*hlt_tauIdOption].begin();
	              hlt_isolationWP != hlt_isolationWPs[*hlt_tauIdOption].end(); ++hlt_isolationWP ) {
                  string_to_TGraphMap1 graphs2 = graphs_hlt_efficiency_vs_leadTrackPt[*hlt_pfAlgo][*hlt_vertexOption][*hlt_tauIdOption][*l1MatchingOption]
                    [*observable][*absEtaRange][*hlt_ptThreshold][*hlt_isolationWP];
  	          bool addFitFunctions2 = false;
                  double legendPosX = 0.61;
	          if ( (*observable) == "pt" ) 
	          {
	            addFitFunctions2 = true;
	          }
                  if ( (*observable) == "eta" ) 
	          {
                    legendPosX = 0.38;
	          }
	          std::vector<std::string> labelTextLines2 = getLabelTextLines(*hlt_ptThreshold);
                  std::string outputFileName2 = Form("makeEfficiencyPlots_%s%s%s_%s_%s_%s_%s_%s_vs_leadTrackPt.png", 
                    hlt_pfAlgo->data(), hlt_vertexOption->data(), l1MatchingOption->data(), 
                    hlt_tauIdOption->data(), observable->data(), absEtaRange->data(), hlt_ptThreshold->data(), hlt_isolationWP->data());
	          showGraphs(1150, 850,
                             graphs2["leadTrackPtGt1"], legendEntries_vs_leadTrackPt["leadTrackPtGt1"],
  		             graphs2["leadTrackPtGt2"], legendEntries_vs_leadTrackPt["leadTrackPtGt2"],
		             graphs2["leadTrackPtGt5"], legendEntries_vs_leadTrackPt["leadTrackPtGt5"],
                             nullptr, "",
                             nullptr, "",
                             nullptr, "",
		             addFitFunctions2,
		             colors, markerStyles, lineStyles, 
		             0.040, legendPosX, 0.17, 0.33, 0.14, 
		             labelTextLines2, 0.050,
		             0.17, 0.85, 0.26, 0.05, 
		             xMin[*observable], xMax[*observable], xAxisTitles[*observable], 1.2, 
		             false, 0., 1.09, "Efficiency", 1.4, 
		             outputFileName2);
                } // hlt_isolationWP
              } // l1MatchingOption

              for ( std::vector<std::string>::const_iterator l1MatchingOption = l1MatchingOptions.begin();
	            l1MatchingOption != l1MatchingOptions.end(); ++l1MatchingOption ) {
                for ( std::vector<std::string>::const_iterator hlt_min_leadTrackPt = hlt_min_leadTrackPtValues.begin();
                      hlt_min_leadTrackPt != hlt_min_leadTrackPtValues.end(); ++hlt_min_leadTrackPt ) {
    	          for ( std::vector<std::string>::const_iterator hlt_isolationWP = hlt_isolationWPs[*hlt_tauIdOption].begin();
	                hlt_isolationWP != hlt_isolationWPs[*hlt_tauIdOption].end(); ++hlt_isolationWP ) {  
	            string_to_TGraphMap1 graphs3 = graphs_hlt_efficiency_vs_decayModes[*hlt_pfAlgo][*hlt_vertexOption][*hlt_tauIdOption][*l1MatchingOption]
                      [*observable][*absEtaRange][*hlt_ptThreshold][*hlt_min_leadTrackPt][*hlt_isolationWP];
	            std::vector<std::string> labelTextLines3 = getLabelTextLines(*hlt_ptThreshold);
	            std::string outputFileName3 = Form("makeEfficiencyPlots_%s%s%s_%s_%s_%s_%s_%s_%s_vs_decayMode.png", 
                      hlt_pfAlgo->data(), hlt_vertexOption->data(), l1MatchingOption->data(), 
                      hlt_tauIdOption->data(), observable->data(), absEtaRange->data(), hlt_ptThreshold->data(), hlt_min_leadTrackPt->data(), hlt_isolationWP->data());
	            showGraphs(1150, 850,
  		               graphs3["oneProng0Pi0"],   legendEntries_vs_decayModes["oneProng0Pi0"],
		               graphs3["oneProng1Pi0"],   legendEntries_vs_decayModes["oneProng1Pi0"],
		               graphs3["oneProng2Pi0"],   legendEntries_vs_decayModes["oneProng2Pi0"],
		               graphs3["threeProng0Pi0"], legendEntries_vs_decayModes["threeProng0Pi0"],
		               graphs3["threeProng1Pi0"], legendEntries_vs_decayModes["threeProng1Pi0"],
		               0, "",
		               false,
		               colors, markerStyles, lineStyles, 
		               0.045, 0.69, 0.17, 0.23, 0.28, 
		               labelTextLines3, 0.045,
		               0.17, 0.85, 0.26, 0.05, 
		               xMin[*observable], xMax[*observable], xAxisTitles[*observable], 1.2, 
		               false, 0., 1.09, "Efficiency", 1.4, 
		               outputFileName3);
                  } // hlt_isolationWP
                } // hlt_min_leadTrackPt
              } // l1MatchingOption

              // CV: make final plot that shows trigger efficiency in steps
              for ( std::vector<std::string>::const_iterator l1MatchingOption = l1MatchingOptions.begin();
	            l1MatchingOption != l1MatchingOptions.end(); ++l1MatchingOption ) {
                for ( std::vector<std::string>::const_iterator hlt_min_leadTrackPt = hlt_min_leadTrackPtValues.begin();
                      hlt_min_leadTrackPt != hlt_min_leadTrackPtValues.end(); ++hlt_min_leadTrackPt ) {
    	          for ( std::vector<std::string>::const_iterator hlt_isolationWP = hlt_isolationWPs[*hlt_tauIdOption].begin();
	                hlt_isolationWP != hlt_isolationWPs[*hlt_tauIdOption].end(); ++hlt_isolationWP ) {  
                    std::vector<TGraph*> graphs;
                    std::vector<std::string> legendEntries;
                    if ( (*l1MatchingOption) == "MatchedToL1" ) 
                    {
                      graphs.push_back(graphs_l1_efficiency[*observable][*absEtaRange]);
                      legendEntries.push_back("L1");
                    }
                    // CV: track finding
                    graphs.push_back(graphs_hlt_efficiency_vs_isolationWPs[*hlt_pfAlgo][*hlt_vertexOption][*hlt_tauIdOption][*l1MatchingOption]
                      [*observable][*absEtaRange][*hlt_ptThreshold]["leadTrackPtGt1"]["noIsolation"]);
                    legendEntries.push_back("Track finding");
                    // CV: lead. track pT > 5 GeV
                    graphs.push_back(graphs_hlt_efficiency_vs_isolationWPs[*hlt_pfAlgo][*hlt_vertexOption][*hlt_tauIdOption][*l1MatchingOption]
                      [*observable][*absEtaRange][*hlt_ptThreshold]["leadTrackPtGt5"]["noIsolation"]);
                    legendEntries.push_back(legendEntries_vs_leadTrackPt["leadTrackPtGt5"]);
                    if ( (*hlt_tauIdOption) == "recoSumChargedIso" || (*hlt_tauIdOption) == "patSumChargedIso" ) 
                    {
                      // CV: charged isolation < 0.05 * tau pT
                      graphs.push_back(graphs_hlt_efficiency_vs_isolationWPs[*hlt_pfAlgo][*hlt_vertexOption][*hlt_tauIdOption][*l1MatchingOption]
                        [*observable][*absEtaRange][*hlt_ptThreshold]["leadTrackPtGt5"]["relChargedIsoLt0p05"]);
                      legendEntries.push_back(legendEntries_vs_isolationWPs[*hlt_tauIdOption]["relDiscriminatorLt0p050"]);
                    }
                    else if ( (*hlt_tauIdOption) == "patDeepTau" )
                    {
                      // CV: VLoose DeepTau WP
                      graphs.push_back(graphs_hlt_efficiency_vs_isolationWPs[*hlt_pfAlgo][*hlt_vertexOption][*hlt_tauIdOption][*l1MatchingOption]
                        [*observable][*absEtaRange][*hlt_ptThreshold]["leadTrackPtGt5"]["absDiscriminatorGt0p598"]);
                      legendEntries.push_back(legendEntries_vs_isolationWPs[*hlt_tauIdOption]["absDiscriminatorGt0p598"]);
                    }
                    else assert(0);
                    TGraph* graph4 = ( graphs.size() >= 4 ) ? graphs[3] : nullptr;
                    const std::string& legendEntry4 = ( legendEntries.size() >= 4 ) ? legendEntries[3] : "";
  	            bool addFitFunctions4 = false;
                    double legendPosX = 0.57;
	            if ( (*observable) == "pt" ) 
	            {
	              addFitFunctions4 = true;
	            }
                    if ( (*observable) == "eta" ) 
	            {
                      legendPosX = 0.36;
	            }
                    double legendPosY  = 0.17;
                    double legendSizeY = 0.14;
                    if ( graph4 ) 
                    {
                      legendPosY  = 0.155;
                      legendSizeY = 0.18;
                    }
	            std::vector<std::string> labelTextLines4 = getLabelTextLines(*hlt_ptThreshold);
                    std::string outputFileName4 = Form("makeEfficiencyPlots_%s%s%s_%s_%s_%s_%s_final.png", 
                      hlt_pfAlgo->data(), hlt_vertexOption->data(), l1MatchingOption->data(), 
                      hlt_tauIdOption->data(), observable->data(), absEtaRange->data(), hlt_ptThreshold->data());
	            showGraphs(1150, 850,
                               graphs[0], legendEntries[0],
                               graphs[1], legendEntries[1],
                               graphs[2], legendEntries[2],
                               graph4,    legendEntry4,
                               nullptr, "",
                               nullptr, "",
		               addFitFunctions4,
		               colors, markerStyles, lineStyles, 
		               0.040, legendPosX, legendPosY, 0.33, legendSizeY,  
		               labelTextLines4, 0.050,
		               0.17, 0.85, 0.26, 0.05, 
		               xMin[*observable], xMax[*observable], xAxisTitles[*observable], 1.2, 
		               false, 0., 1.09, "Efficiency", 1.4, 
  		               outputFileName4);
	          } // hlt_isolationWP
  	        } // hlt_min_leadTrackPt
              } // l1MatchingOption
 */
            } // hlt_ptThreshold
          } // absEtaRange
        } // observable
      } // hlt_tauIdOption
    } // hlt_vertexOption
  } // hlt_pfAlgo

  delete inputFile;
}

