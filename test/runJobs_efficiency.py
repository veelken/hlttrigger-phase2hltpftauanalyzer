#!/usr/bin/env python

import os

from HLTTrigger.TallinnHLTPFTauAnalyzer.tools import getInputFileNames

signal_samples = {
  'qqH_htt' : {
    'inputFilePath' : {
      'offlinePrimaryVertices' : '/hdfs/cms/store/user/rdewanje/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5/HLTConfig_w_offlineVtxCollection_HGCalFix_VBFHTT_Phase2HLTTDRWinter20_PU200_CMSSW_11_1_0_pre8/200627_142633/',
      'hltPhase2PixelVertices' : '/hdfs/cms/store/user/rdewanje/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5/HLTConfig_w_onlineVtxCollection_HGCalFix_VBFHTT_Phase2HLTTDRWinter20_PU200_CMSSW_11_1_0_pre8/200627_141415/'  
    },
    'numJobs' : 32,
    'process' : "qqH_htt"
  }
}

version = "2020Jun31"

def runCommand(command):
  #print("executing command = '%s'" % command)
  os.system(command)

def buildConfigFile(cfgFile_original, cfgFile_modified, inputFileNames, process, srcVertices, outputFileName):
  rmCommand   = 'rm -f %s' % cfgFile_modified
  runCommand(rmCommand)
 
  sedCommand  = 'sed'
  sedCommand += ' "s/##inputFilePath/inputFilePath/; s/\$inputFilePath/None/;'
  sedCommand += '  s/##inputFileNames/inputFileNames/; s/\$inputFileNames/%s/;' % [ inputFileName.replace("/", "\/") for inputFileName in inputFileNames ]
  sedCommand += '  s/##processName/processName/; s/\$processName/%s/;' % process
  sedCommand += '  s/##srcVertices/srcVertices/; s/\$srcVertices/%s/;' % srcVertices
  sedCommand += '  s/##outputFileName/outputFileName/; s/\$outputFileName/%s/"' % outputFileName
  sedCommand += ' %s > %s' % (cfgFile_original, cfgFile_modified)
  runCommand(sedCommand)
 
makeCommands = []
outputFileNames = []
for sampleName, sample in signal_samples.items():
  for srcVertices in [ "offlinePrimaryVertices", "hltPhase2PixelVertices" ]:
    print("processing sample = '%s': srcVertices = '%s'" % (sampleName, srcVertices)) 
    inputFilePath = sample['inputFilePath'][srcVertices]
    print(" inputFilePath = '%s'" % inputFilePath)
    inputFileNames = getInputFileNames.getInputFileNames(inputFilePath)
    numInputFiles = len(inputFileNames)
    print("Found %i input files." % numInputFiles)
    numJobs = sample['numJobs']
    for jobId in range(numJobs):
      idxFirstFile = jobId*numInputFiles/numJobs
      idxLastFile = (jobId + 1)*numInputFiles/numJobs - 1
      inputFileNames_job = inputFileNames[idxFirstFile:idxLastFile + 1]
      #print("job #%i: inputFiles = %s" % (jobId, inputFileNames_job))
      cfgFileName_modified = "analyzePFTaus_signal_%s_%s_%i_cfg.py" % (sampleName, srcVertices, jobId)
      outputFileName = "analyzePFTaus_signal_%s_%s_%s_%i.root" % (sampleName, srcVertices, version, jobId)
      outputFileNames.append(outputFileName)
      buildConfigFile("analyzePFTaus_signal_cfg.py", cfgFileName_modified, inputFileNames_job, sample['process'], srcVertices, outputFileName)
      logFileName = cfgFileName_modified.replace("_cfg.py", ".log")
      makeCommands.append("%s:" % outputFileName)
      rmCommand = 'rm -f %s' % logFileName
      makeCommands.append("\t%s" % rmCommand)
      cmsRunCommand = 'cmsRun %s >& %s' % (cfgFileName_modified, logFileName)
      makeCommands.append("\t%s" % cmsRunCommand)

outputFileName_hadd = "analyzePFTaus_signal_all_%s.root" % version
makeCommands.append("%s: %s" % (outputFileName_hadd, " ".join(outputFileNames)))
rmCommand = 'rm -f %s' % outputFileName_hadd
makeCommands.append("\t%s" % rmCommand)
if len(outputFileNames) >= 2:
  haddCommand = 'hadd %s %s' % (outputFileName_hadd, " ".join(outputFileNames))
  makeCommands.append("\t%s" % haddCommand)
elif len(outputFileNames) == 1:
  cpCommand = 'cp %s %s' % (outputFileNames[0], outputFileName_hadd)
  makeCommands.append("\t%s" % cpCommand)
else:
  raise ValueError("No samples defined, so there is nothing to do !!")

print("Building Makefile: version = '%s'" % version)
makeFileName = "Makefile_runJobs_efficiency"
makeFile = open(makeFileName, "w") 
makeFile.write("all: %s\n" % outputFileName_hadd)
for makeCommand in makeCommands:
    makeFile.write("%s\n" % makeCommand)
makeFile.close()

print("Finished building config files. Now execute 'make -j 16 -f %s'." % makeFileName)
