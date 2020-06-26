
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
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

TGraph* makeEfficiencyGraph(TH1* histogram_numerator, TH1* histogram_denominator)
{
  std::cout << "<makeEfficiencyGraph>:" << std::endl;
  std::cout << " numerator = " << histogram_numerator->GetName() << ", #bins = " << histogram_numerator->GetNbinsX() << std::endl;
  std::cout << " denominator = " << histogram_denominator->GetName() << ", #bins = " << histogram_denominator->GetNbinsX() << std::endl;
  assert(histogram_numerator->GetNbinsX() == histogram_denominator->GetNbinsX());
  int numBinsX = histogram_numerator->GetNbinsX();
  for ( int idxBinX = 1; idxBinX <= numBinsX; ++idxBinX ) { 
    double binContent_numerator = histogram_numerator->GetBinContent(idxBinX);
    double binContent_denominator = histogram_denominator->GetBinContent(idxBinX);
    std::cout << " bin #" << idxBinX << ": numerator = " << binContent_numerator << ", denominator = " << binContent_denominator << std::endl;
    if( binContent_numerator > binContent_denominator ) {
      std::cerr << "Error in <makeEfficiencyGraph>: numerator = " << binContent_numerator << " exceeds denominator = " << binContent_denominator 
		<< " @ x = " << histogram_denominator->GetBinCenter(idxBinX) << " !!" << std::endl;
      assert(0);
    }
  }
  TString graphName_efficiency = TString(histogram_numerator->GetName()).ReplaceAll("_numerator", "");
  TGraphAsymmErrors* graph_efficiency = new TGraphAsymmErrors(histogram_numerator, histogram_denominator, "w");
  graph_efficiency->SetName(graphName_efficiency.Data());
  return graph_efficiency;
}

void showGraphs(double canvasSizeX, double canvasSizeY,
		TGraph* graph1, const std::string& legendEntry1,
		TGraph* graph2, const std::string& legendEntry2,
		TGraph* graph3, const std::string& legendEntry3,
		TGraph* graph4, const std::string& legendEntry4,
		TGraph* graph5, const std::string& legendEntry5,
		TGraph* graph6, const std::string& legendEntry6,
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

  dummyHistogram->Draw("axis");

  graph1->SetMarkerColor(colors[0]);
  graph1->SetMarkerSize(2);
  graph1->SetMarkerStyle(markerStyles[0]);
  graph1->SetLineColor(colors[0]);
  graph1->SetLineWidth(2);
  graph1->SetLineStyle(lineStyles[0]);
  graph1->Draw("P");

  if ( graph2 ) {
    graph2->SetMarkerColor(colors[1]);
    graph2->SetMarkerSize(2);
    graph2->SetMarkerStyle(markerStyles[1]);
    graph2->SetLineColor(colors[1]);
    graph2->SetLineWidth(2);
    graph2->SetLineStyle(lineStyles[1]);
    graph2->Draw("P");
  }

  if ( graph3 ) {
    graph3->SetMarkerColor(colors[2]);
    graph3->SetMarkerSize(2);
    graph3->SetMarkerStyle(markerStyles[2]);
    graph3->SetLineColor(colors[2]);
    graph3->SetLineWidth(2);
    graph3->SetLineStyle(lineStyles[2]);
    graph3->Draw("P");
  }

  if ( graph4 ) {
    graph4->SetMarkerColor(colors[3]);
    graph4->SetMarkerSize(2);
    graph4->SetMarkerStyle(markerStyles[3]);
    graph4->SetLineColor(colors[3]);
    graph4->SetLineWidth(2);
    graph4->SetLineStyle(lineStyles[3]);
    graph4->Draw("P");
  }

  if ( graph5 ) {
    graph5->SetMarkerColor(colors[4]);
    graph5->SetMarkerSize(2);
    graph5->SetMarkerStyle(markerStyles[4]);
    graph5->SetLineColor(colors[4]);
    graph5->SetLineWidth(2);
    graph5->SetLineStyle(lineStyles[4]);
    graph5->Draw("P");
  }

  if ( graph6 ) {
    graph6->SetMarkerColor(colors[5]);
    graph6->SetMarkerSize(2);
    graph6->SetMarkerStyle(markerStyles[5]);
    graph6->SetLineColor(colors[5]);
    graph6->SetLineWidth(2);
    graph6->SetLineStyle(lineStyles[5]);
    graph6->Draw("P");
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
  delete dummyHistogram;
}

void makeTrackingEfficiencyPlots()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  std::string inputFilePath = Form("%s/src/HLTTrigger/TallinnHLTPFTauAnalyzer/test/", gSystem->Getenv("CMSSW_BASE"));
  std::string inputFileName = "analyzePFTaus_signal_2020Jun24.root";
  std::string inputFileName_full = inputFilePath;
  if ( inputFileName_full.find_last_of("/") != (inputFileName_full.size() - 1) ) inputFileName_full.append("/");
  inputFileName_full.append(inputFileName);
  TFile* inputFile = new TFile(inputFileName_full.data());
  if ( !inputFile ) {
    std::cerr << "Failed to open input file = '" << inputFileName_full << "' !!" << std::endl;
    assert(0);
  }

  std::vector<std::string> recTrack_types;
  recTrack_types.push_back("hltTrack");
  recTrack_types.push_back("hltPFCandTrack");
  recTrack_types.push_back("offlineTrack");
  recTrack_types.push_back("offlinePFCandTrack");

  std::vector<std::string> recTrack_options;
  recTrack_options.push_back("woQualityCuts");
  recTrack_options.push_back("wQualityCuts");

  std::vector<std::string> vtxModes;
  vtxModes.push_back("GenVertex");
  vtxModes.push_back("RecVertex");
  
  std::vector<std::string> observables;
  observables.push_back("pt");
  observables.push_back("eta");
  observables.push_back("phi");
  observables.push_back("minDeltaR");

  std::vector<std::string> absEtaRanges;
  absEtaRanges.push_back("absEtaLt1p40");
  absEtaRanges.push_back("absEta1p40to2p17");
  //absEtaRanges.push_back("absEta1p40to2p40");
  absEtaRanges.push_back("absEtaLt2p17");
  //absEtaRanges.push_back("absEtaLt2p40");

  std::vector<std::string> decayModes;
  //decayModes.push_back("oneProng0Pi0");
  //decayModes.push_back("oneProng1Pi0");
  //decayModes.push_back("oneProng2Pi0");
  //decayModes.push_back("threeProng0Pi0");
  //decayModes.push_back("threeProng1Pi0");
  decayModes.push_back("all");

  std::map<std::string, double> xMin; // key = observable
  xMin["pt"]                                            =   0.;
  xMin["eta"]                                           =  -3.0;
  xMin["phi"]                                           =  -TMath::Pi();
  xMin["minDeltaR"]                                     =   0.;

  std::map<std::string, double> xMax; // key = observable
  xMax["pt"]                                            = 100.;
  xMax["eta"]                                           =  +3.0;
  xMax["phi"]                                           =  +TMath::Pi();
  xMax["minDeltaR"]                                     =  0.03;

  std::map<std::string, std::string> xAxisTitles; // key = observable
  xAxisTitles["pt"]                                     = "p_{T} [GeV]";
  xAxisTitles["eta"]                                    = "#eta";
  xAxisTitles["phi"]                                    = "#phi";
  xAxisTitles["minDeltaR"]                              = "min(#Delta R)";

  std::map<std::string, std::string> legendEntries_vs_recTrack_types; // key = recTrack_type
  legendEntries_vs_recTrack_types["hltTrack"]           = "Online Track";
  legendEntries_vs_recTrack_types["hltPFCandTrack"]     = "Online PFlow";
  legendEntries_vs_recTrack_types["offlineTrack"]       = "Offline Track";
  legendEntries_vs_recTrack_types["offlinePFCandTrack"] = "Offline PFlow";

  std::map<std::string, std::string> legendEntries_vs_recTrack_options; // key = recTrack_option
  legendEntries_vs_recTrack_options["woQualityCuts"]    = "Without Quality Cuts";
  legendEntries_vs_recTrack_options["wQualityCuts"]     = "With Quality Cuts";

  std::map<std::string, std::string> legendEntries_vs_decayModes; // key = decayMode
  legendEntries_vs_decayModes["oneProng0Pi0"]           = "h^{#pm}";
  legendEntries_vs_decayModes["oneProng1Pi0"]           = "h^{#pm}#pi^{0}";
  legendEntries_vs_decayModes["oneProng2Pi0"]           = "h^{#pm}#pi^{0}#pi^{0}";
  legendEntries_vs_decayModes["threeProng0Pi0"]         = "h^{#pm}h^{#mp}h^{#pm}";
  legendEntries_vs_decayModes["threeProng1Pi0"]         = "h^{#pm}h^{#mp}h^{#pm}#pi^{0}";
  legendEntries_vs_decayModes["all"]                    = "all";

  std::string dqmDirectory = "recoTrackAnalyzer";

  int colors[6] = { 1, 2, 8, 4, 6, 7 };
  int lineStyles[6] = { 1, 1, 1, 1, 1, 1 };
  int markerStyles[6] = { 22, 32, 20, 24, 21, 25 };

  typedef std::map<std::string, TGraph*>               string_to_graph_map_1;
  typedef std::map<std::string, string_to_graph_map_1> string_to_graph_map_2;
  typedef std::map<std::string, string_to_graph_map_2> string_to_graph_map_3;
  typedef std::map<std::string, string_to_graph_map_3> string_to_graph_map_4;
  typedef std::map<std::string, string_to_graph_map_4> string_to_graph_map_5;
  typedef std::map<std::string, string_to_graph_map_5> string_to_graph_map_6;
  string_to_graph_map_6 graphs_efficiency; // key = recTrack_type, recTrack_option, vtxMode, observable, absEtaRange, decayMode

  for ( std::vector<std::string>::const_iterator vtxMode = vtxModes.begin();
	vtxMode != vtxModes.end(); ++vtxMode ) {
    for ( std::vector<std::string>::const_iterator recTrack_type = recTrack_types.begin();
	  recTrack_type != recTrack_types.end(); ++recTrack_type ) {
      for ( std::vector<std::string>::const_iterator recTrack_option = recTrack_options.begin();
	    recTrack_option != recTrack_options.end(); ++recTrack_option ) {
        for ( std::vector<std::string>::const_iterator observable = observables.begin();
	      observable != observables.end(); ++observable ) {
  	  for ( std::vector<std::string>::const_iterator absEtaRange = absEtaRanges.begin();
	        absEtaRange != absEtaRanges.end(); ++absEtaRange ) {
	    for ( std::vector<std::string>::const_iterator decayMode = decayModes.begin();
		  decayMode != decayModes.end(); ++decayMode ) {
	      std::string recTrack_type_capitalized = *recTrack_type;
	      recTrack_type_capitalized[0] = toupper(recTrack_type_capitalized[0]);
	      std::string decayMode_capitalized = *decayMode;
	      decayMode_capitalized[0] = toupper(decayMode_capitalized[0]);
	      std::string dqmDirectory_full = Form("%sWrt%s/%s/gen%sTau/%s_%s", 
                dqmDirectory.data(), vtxMode->data(), absEtaRange->data(), decayMode_capitalized.data(), recTrack_type->data(), recTrack_option->data());	    
	      std::string histogramName_numerator = Form("%s/eff%s_%s_vs_%s_numerator_%s_pt1p00to1000p00", 
	        dqmDirectory_full.data(), recTrack_type_capitalized.data(), recTrack_option->data(), observable->data(), absEtaRange->data());
	      if ( (*decayMode) != "all" ) histogramName_numerator.append(Form("_gen%sTau", decayMode_capitalized.data()));
	      TH1* histogram_numerator = loadHistogram(inputFile, histogramName_numerator);
	      std::string histogramName_denominator = Form("%s/eff%s_%s_vs_%s_denominator_%s_pt1p00to1000p00", 
	        dqmDirectory_full.data(), recTrack_type_capitalized.data(), recTrack_option->data(), observable->data(), absEtaRange->data());
	      if ( (*decayMode) != "all" ) histogramName_denominator.append(Form("_gen%sTau", decayMode_capitalized.data()));
	      TH1* histogram_denominator = loadHistogram(inputFile, histogramName_denominator);
	      TGraph* graph_efficiency = makeEfficiencyGraph(histogram_numerator, histogram_denominator);
	      graphs_efficiency[*recTrack_type][*recTrack_option][*vtxMode][*observable][*absEtaRange][*decayMode] = graph_efficiency;
	    }
	  }
	}
      }
    }
  }

  // 1st set of plots: comparing efficiencies of different algorithms (HLT Track, HLT PFlow, offline Track, offline PFlow)
  for ( std::vector<std::string>::const_iterator vtxMode = vtxModes.begin();
	vtxMode != vtxModes.end(); ++vtxMode ) {
    for ( std::vector<std::string>::const_iterator recTrack_option = recTrack_options.begin();
	  recTrack_option != recTrack_options.end(); ++recTrack_option ) {
      for ( std::vector<std::string>::const_iterator absEtaRange = absEtaRanges.begin();
	    absEtaRange != absEtaRanges.end(); ++absEtaRange ) {
        for ( std::vector<std::string>::const_iterator observable = observables.begin();
	      observable != observables.end(); ++observable ) {
	  TGraph* graph_hltTrack = graphs_efficiency["hltTrack"][*recTrack_option][*vtxMode][*observable][*absEtaRange]["all"];
	  assert(graph_hltTrack);
	  TGraph* graph_hltPFlow = graphs_efficiency["hltPFCandTrack"][*recTrack_option][*vtxMode][*observable][*absEtaRange]["all"];
	  assert(graph_hltPFlow);
	  TGraph* graph_offlineTrack = graphs_efficiency["offlineTrack"][*recTrack_option][*vtxMode][*observable][*absEtaRange]["all"];
	  assert(graph_offlineTrack);
	  TGraph* graph_offlinePFlow = graphs_efficiency["offlinePFCandTrack"][*recTrack_option][*vtxMode][*observable][*absEtaRange]["all"];
	  assert(graph_offlinePFlow);
	
  	  double legendPosX = 0.68;
          double legendPosY = 0.17;
          if ( (*observable) == "minDeltaR" ) {
    	    legendPosX = 0.18;
	    legendPosY = 0.17;
          }
  	  std::vector<std::string> labelTextLines;
	  std::string outputFileName = Form("makeTrackingEfficiencyPlots_vs_algo_and_%s_%s_wrt%s_%s.png", 
            observable->data(), recTrack_option->data(), vtxMode->data(), absEtaRange->data());
	  showGraphs(1150, 850,
		     graph_hltTrack,     legendEntries_vs_recTrack_types["hltTrack"],
		     graph_hltPFlow,     legendEntries_vs_recTrack_types["hltPFCandTrack"],
		     graph_offlineTrack, legendEntries_vs_recTrack_types["offlineTrack"],
		     graph_offlinePFlow, legendEntries_vs_recTrack_types["offlinePFCandTrack"],
		     0, "",
		     0, "",
		     colors, markerStyles, lineStyles, 
		     0.045, legendPosX, legendPosY, 0.25, 0.26, 
		     labelTextLines, 0.050,
		     0.63, 0.65, 0.26, 0.07, 
		     xMin[*observable], xMax[*observable], xAxisTitles[*observable], 1.2, 
		     false, 0., 1.09, "Efficiency", 1.4, 
		     outputFileName);
	}
      }
    }
  }
	
  // 2nd set of plots: comparing efficiencies obtained before and after quality cuts are applied
  for ( std::vector<std::string>::const_iterator vtxMode = vtxModes.begin();
	vtxMode != vtxModes.end(); ++vtxMode ) {
    for ( std::vector<std::string>::const_iterator recTrack_type = recTrack_types.begin();
	  recTrack_type != recTrack_types.end(); ++recTrack_type ) {
      for ( std::vector<std::string>::const_iterator absEtaRange = absEtaRanges.begin();
	    absEtaRange != absEtaRanges.end(); ++absEtaRange ) {
        for ( std::vector<std::string>::const_iterator observable = observables.begin();
	      observable != observables.end(); ++observable ) {
	  TGraph* graph_woQualityCuts = graphs_efficiency[*recTrack_type]["woQualityCuts"][*vtxMode][*observable][*absEtaRange]["all"];
	  assert(graph_woQualityCuts);
	  TGraph* graph_wQualityCuts = graphs_efficiency[*recTrack_type]["wQualityCuts"][*vtxMode][*observable][*absEtaRange]["all"];
	  assert(graph_wQualityCuts);
	
  	  double legendPosX = 0.68;
          double legendPosY = 0.17;
          if ( (*observable) == "minDeltaR" ) {
    	    legendPosX = 0.18;
	    legendPosY = 0.17;
          }
	  std::vector<std::string> labelTextLines;
	  std::string outputFileName = Form("makeTrackingEfficiencyPlots_vs_recTrackOption_and_%s_%s_wrt%s_%s.png", 
            observable->data(), recTrack_type->data(), vtxMode->data(), absEtaRange->data());
	  showGraphs(1150, 850,
		     graph_woQualityCuts, legendEntries_vs_recTrack_options["woQualityCuts"],
		     graph_wQualityCuts,  legendEntries_vs_recTrack_options["wQualityCuts"],
		     0, "",
		     0, "",
		     0, "",
		     0, "",
		     colors, markerStyles, lineStyles, 
		     0.045, legendPosX, legendPosY, 0.32, 0.13, 
		     labelTextLines, 0.050,
		     0.63, 0.65, 0.26, 0.07, 
		     xMin[*observable], xMax[*observable], xAxisTitles[*observable], 1.2, 
		     false, 0., 1.09, "Efficiency", 1.4, 
		     outputFileName);
	}
      }
    }
  }
/*
  // 3rd set of plots: comparing efficiencies for different generator-level tau decay modes
  for ( std::vector<std::string>::const_iterator vtxMode = vtxModes.begin();
	vtxMode != vtxModes.end(); ++vtxMode ) {
    for ( std::vector<std::string>::const_iterator recTrack_type = recTrack_types.begin();
	  recTrack_type != recTrack_types.end(); ++recTrack_type ) {
      for ( std::vector<std::string>::const_iterator recTrack_option = recTrack_options.begin();
	    recTrack_option != recTrack_options.end(); ++recTrack_option ) {
        for ( std::vector<std::string>::const_iterator absEtaRange = absEtaRanges.begin();
	      absEtaRange != absEtaRanges.end(); ++absEtaRange ) {
  	  for ( std::vector<std::string>::const_iterator observable = observables.begin();
		observable != observables.end(); ++observable ) {
	    TGraph* graph_oneProng0Pi0 = graphs_efficiency[*recTrack_type][*recTrack_option][*vtxMode][*observable][*absEtaRange]["oneProng0Pi0"];
	    assert(graph_oneProng0Pi0);
	    TGraph* graph_oneProng1Pi0 = graphs_efficiency[*recTrack_type][*recTrack_option][*vtxMode][*observable][*absEtaRange]["oneProng1Pi0"];
	    assert(graph_oneProng1Pi0);
	    TGraph* graph_oneProng2Pi0 = graphs_efficiency[*recTrack_type][*recTrack_option][*vtxMode][*observable][*absEtaRange]["oneProng2Pi0"];
	    assert(graph_oneProng2Pi0);
	    TGraph* graph_threeProng0Pi0 = graphs_efficiency[*recTrack_type][*recTrack_option][*vtxMode][*observable][*absEtaRange]["threeProng0Pi0"];
	    assert(graph_threeProng0Pi0);
	    TGraph* graph_threeProng1Pi0 = graphs_efficiency[*recTrack_type][*recTrack_option][*vtxMode][*observable][*absEtaRange]["threeProng1Pi0"];
	    assert(graph_threeProng1Pi0);
	    
	    std::vector<std::string> labelTextLines;
	    std::string outputFileName = Form("makeTrackingEfficiencyPlots_vs_decayMode_and_%s_%s_%s_wrt%s_%s.png", 
              observable->data(), recTrack_type->data(), recTrack_option->data(), vtxMode->data(), absEtaRange->data());
	    showGraphs(1150, 850,
		       graph_oneProng0Pi0,   legendEntries_vs_decayModes["oneProng0Pi0"],
		       graph_oneProng1Pi0,   legendEntries_vs_decayModes["oneProng1Pi0"],
		       graph_oneProng2Pi0,   legendEntries_vs_decayModes["oneProng2Pi0"],
		       graph_threeProng0Pi0, legendEntries_vs_decayModes["threeProng0Pi0"],
		       graph_threeProng1Pi0, legendEntries_vs_decayModes["threeProng1Pi0"],
		       0, "",
		       colors, markerStyles, lineStyles, 
		       0.045, 0.68, 0.17, 0.25, 0.26, 
		       labelTextLines, 0.050,
		       0.63, 0.65, 0.26, 0.07, 
		       xMin[*observable], xMax[*observable], xAxisTitles[*observable], 1.2, 
		       false, 0., 1.09, "Efficiency", 1.4, 
		       outputFileName);
	  }
	}
      }
    }
  }
 */

  delete inputFile;
}

