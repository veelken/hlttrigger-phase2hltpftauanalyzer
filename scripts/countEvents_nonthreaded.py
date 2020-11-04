#!/usr/bin/env python

import getpass
import os

from HLTrigger.TallinnHLTPFTauAnalyzer.tools.jobTools import getInputFileNames
from HLTrigger.TallinnHLTPFTauAnalyzer.tools.safe_root import ROOT

inputFilePaths = [
  #'/hdfs/cms/store/user/sbhowmik/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5/VBFHToTauTau_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1113_20201009/',
  #'/hdfs/cms/store/user/sbhowmik/MinBias_TuneCP5_14TeV-pythia8/MinBias_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1113_20201009/',
  #'/hdfs/cms/store/user/sbhowmik/QCD_Pt_30to50_TuneCP5_14TeV_pythia8/QCDPt30to50_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1113_20201009/',
  #'/hdfs/cms/store/user/sbhowmik/QCD_Pt_50to80_TuneCP5_14TeV_pythia8/QCDPt50to80_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1113_20201009/',
  #'/hdfs/cms/store/user/sbhowmik/QCD_Pt_80to120_TuneCP5_14TeV_pythia8/QCDPt80to120_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1113_20201009/',
  #'/hdfs/cms/store/user/sbhowmik/QCD_Pt_120to170_TuneCP5_14TeV_pythia8/QCDPt120to170_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1113_20201009/',
  #'/hdfs/cms/store/user/sbhowmik/QCD_Pt_170to300_TuneCP5_14TeV_pythia8/QCDPt170to300_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1113_20201009/',
  #'/hdfs/cms/store/user/sbhowmik/QCD_Pt_300to470_TuneCP5_14TeV_pythia8/QCDPt300to470_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1113_20201009/',
  #'/hdfs/cms/store/user/sbhowmik/QCD_Pt_470to600_TuneCP5_14TeV_pythia8/QCDPt470to600_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1113_20201009/',   
  #'/hdfs/cms/store/user/sbhowmik/QCD_Pt_600oInf_TuneCP5_14TeV_pythia8/QCDPt600toInf_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1113_20201009/',
  '/hdfs/cms/store/user/sbhowmik/GluGluHToTauTau_M125_14TeV_powheg_pythia8_TuneCP5/GluGluHToTauTau_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1113_20201009/',
  '/hdfs/cms/store/user/sbhowmik/DYToLL_M-50_TuneCP5_14TeV-pythia8/DYJetsToLL_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1113_20201009/',
  '/hdfs/cms/store/user/sbhowmik/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/WJetsToLNu_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1113_20201009/',
  '/hdfs/cms/store/user/sbhowmik/TTTo2L2Nu_TuneCP5_14TeV-powheg-pythia8/TTTo2L2Nu_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1113_20201009/',
  '/hdfs/cms/store/user/sbhowmik/TTToSemiLepton_TuneCP5_14TeV-powheg-pythia8/TTToSemiLepton_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1113_20201009/',
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
            if not inputFile or inputFile.IsZombie():
                print("Failed to open input file = '%s' --> skipping !!" % inputFileName)
                continue
            events_tree = inputFile.Get(events_treeName)
            numEvents += events_tree.GetEntries()
        print(" numEvents = %i" % numEvents)
        print("")
    else:
        print("Path = '%s' does not exist --> skipping !!" % inputFilePath)
