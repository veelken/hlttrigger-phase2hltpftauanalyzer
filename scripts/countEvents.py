#!/usr/bin/env python

import getpass
import os

from HLTrigger.TallinnHLTPFTauAnalyzer.tools.jobTools import getInputFileNames
from HLTrigger.TallinnHLTPFTauAnalyzer.tools.safe_root import ROOT

inputFilePaths = [
  ##'/hdfs/cms/store/user/rdewanje/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5/HLTConfig_VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5_wOfflineVtx_wL1_2FM/',
  ##'/hdfs/cms/store/user/rdewanje/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5/HLTConfig_VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5_wOnlineVtx_wL1_2FM/',
  ##'/hdfs/cms/store/user/rdewanje/MinBias_TuneCP5_14TeV-pythia8/HLTConfig_MinBias_TuneCP5_14TeV-pythia8_wOfflineVtx_wL1/',
  ##'/hdfs/cms/store/user/rdewanje/MinBias_TuneCP5_14TeV-pythia8/HLTConfig_MinBias_TuneCP5_14TeV-pythia8_wOnlineVtx_wL1/',
  '/hdfs/cms/store/user/rdewanje/QCD_Pt_30to50_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_30to50_TuneCP5_14TeV_pythia8_wOfflineVtx_wL1/', # CV: sample still growing !!
  ##'/hdfs/cms/store/user/rdewanje/QCD_Pt_30to50_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_30to50_TuneCP5_14TeV_pythia8_wOnlineVtx_wL1/',
  ##'/hdfs/cms/store/user/rdewanje/QCD_Pt_50to80_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_50to80_TuneCP5_14TeV_pythia8_wOfflineVtx_wL1/',
  ##'/hdfs/cms/store/user/rdewanje/QCD_Pt_50to80_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_50to80_TuneCP5_14TeV_pythia8_wOnlineVtx_wL1/', 
  ##'/hdfs/cms/store/user/rdewanje/QCD_Pt_80to120_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_80to120_TuneCP5_14TeV_pythia8_wOfflineVtx_wL1/',
  ##'/hdfs/cms/store/user/rdewanje/QCD_Pt_80to120_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_80to120_TuneCP5_14TeV_pythia8_wOnlineVtx_wL1/',
  ##'/hdfs/cms/store/user/rdewanje/QCD_Pt_120to170_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_120to170_TuneCP5_14TeV_pythia8_wOfflineVtx_wL1/',
  ##'/hdfs/cms/store/user/rdewanje/QCD_Pt_120to170_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_120to170_TuneCP5_14TeV_pythia8_wOnlineVtx_wL1/',
  ##'/hdfs/cms/store/user/rdewanje/QCD_Pt_170to300_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_170to300_TuneCP5_14TeV_pythia8_wOfflineVtx_wL1/',
  ##'/hdfs/cms/store/user/rdewanje/QCD_Pt_170to300_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_170to300_TuneCP5_14TeV_pythia8_wOnlineVtx_wL1/', 
  ##'/hdfs/cms/store/user/rdewanje/QCD_Pt_300to470_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_300to470_TuneCP5_14TeV_pythia8_wOfflineVtx_wL1/',
  ##'/hdfs/cms/store/user/rdewanje/QCD_Pt_300to470_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_300to470_TuneCP5_14TeV_pythia8_wOnlineVtx_wL1/', 
  ##'/hdfs/cms/store/user/rdewanje/DYJetsToLL_M-10to50_TuneCP5_14TeV-madgraphMLM-pythia8/HLTConfig_DYJetsToLL_M-10to50_TuneCP5_14TeV-madgraphMLM-pythia8_wOfflineVtx_wL1/',
  ##'/hdfs/cms/store/user/rdewanje/DYJetsToLL_M-10to50_TuneCP5_14TeV-madgraphMLM-pythia8/HLTConfig_DYJetsToLL_M-10to50_TuneCP5_14TeV-madgraphMLM-pythia8_wOnlineVtx_wL1/',
  ##'/hdfs/cms/store/user/rdewanje/DYToLL_M-50_TuneCP5_14TeV-pythia8/HLTConfig_DYToLL_M-50_TuneCP5_14TeV-pythia8_wOfflineVtx_wL1/',
  ##'/hdfs/cms/store/user/rdewanje/DYToLL_M-50_TuneCP5_14TeV-pythia8/HLTConfig_DYToLL_M-50_TuneCP5_14TeV-pythia8_wOnlineVtx_wL1/',
  ##'/hdfs/cms/store/user/rdewanje/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/HLTConfig_WJetsToLNu_TuneCP5_14TeV_pythia8_wOfflineVtx_wL1/',
  ##'/hdfs/cms/store/user/rdewanje/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/HLTConfig_WJetsToLNu_TuneCP5_14TeV_pythia8_wOnlineVtx_wL1/'
]

events_treeName = 'Events'

for inputFilePath in inputFilePaths:
  print("Searching for input files in path = '%s'" % inputFilePath)
  if os.path.exists(inputFilePath):
    inputFileNames = getInputFileNames(inputFilePath)
    print("Found %i input files." % len(inputFileNames))
    numEvents = 0
    for idx, inputFileName in enumerate(inputFileNames):
      if (idx + 1) % 100 == 0:
        print(" Processing %ith input file..." % (idx + 1))
      inputFile = ROOT.TFile.Open(inputFileName, 'read')
      events_tree = inputFile.Get(events_treeName)
      numEvents += events_tree.GetEntries()
    print(" numEvents = %i" % numEvents)
    print("")
  else:
    print("Path = '%s' does not exist --> skipping !!" % inputFilePath)
