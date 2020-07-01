#!/usr/bin/env python

import os

from HLTTrigger.TallinnHLTPFTauAnalyzer.tools import getInputFileNames

# CV: define instantaneous luminosity for HL-LHC running period taken from the twiki
#       https://twiki.cern.ch/twiki/bin/viewauth/CMS/HighLevelTriggerPhase2#Rate_calculations
instLuminosity = 7.5e+34 # cm^-2 s^-1

# CV: define conversion factor from pb to cm^2, taken from wikipedia
#       https://en.wikipedia.org/wiki/Barn_(unit)
conversionFactor = 1.e-36

background_samples = {
  # CV: minbias cross-section taken to be 75mb, 
  #     resulting in 200 pileup interaction per bunch-crossing for instantaneous luminosity of 7.5e+34 cm^-2 s^-1
  #    (assuming that 70% of all bunches are colliding, i.e. collision frequence of 28 MHz)
  'minbias' : {
    'inputFilePath' : {
      'offlinePrimaryVertices' : '/hdfs/cms/store/user/rdewanje/MinBias_TuneCP5_14TeV-pythia8/HLTConfig_w_offlineVtxCollection_HGCalFix_MINBIAS_Phase2HLTTDRWinter20_PU200_CMSSW_11_1_0_pre8/200627_142511/',
      'hltPhase2PixelVertices' : '/hdfs/cms/store/user/rdewanje/MinBias_TuneCP5_14TeV-pythia8/HLTConfig_w_onlineVtxCollection_HGCalFix_MINBIAS_Phase2HLTTDRWinter20_PU200_CMSSW_11_1_0_pre8/200627_135114/'  
    },
    'numJobs' : 32, 
    'crossSection' : 75.0e+9/200, # pb
    'process' : "minbias"
  },
  # CV: cross-sections for QCD, Drell-Yan, and W+jets production taken from the twiki
  #       https://twiki.cern.ch/twiki/bin/viewauth/CMS/HighLevelTriggerPhase2#Rate_calculations
##  'qcd_pt15to20' : {
##    'inputFilePath' : "",
##    'crossSection' : 923300000.0,
##    'process' : "QCD"
##  },
##  'qcd_pt20to30' : {
##    'inputFilePath' : "",
##    'crossSection' : 436000000.0,
##    'process' : "QCD"
##  },
##  'qcd_pt30to50' : {
##    'inputFilePath' : "",
##    'crossSection' : 118400000.0,
##    'process' : "QCD"
##  },
##  'qcd_pt50to80' : {
##    'inputFilePath' : "",
##    'crossSection' : 17650000.0,
##    'process' : "QCD"
##  },
##  'qcd_pt80to120' : {
##    'inputFilePath' : "",
##    'crossSection' : 2671000.0,
##    'process' : "QCD"
##  },
##  'qcd_pt120to170' : {
##    'inputFilePath' : "",
##    'crossSection' : 469700.0 ,
##    'process' : "QCD"
##  },
##  'qcd_pt170to300' : {
##    'inputFilePath' : "",
##    'crossSection' : 121700.0,
##    'process' : "QCD"
##  },  
##  'qcd_ptGt300' : {
##    'inputFilePath' : "",
##    'crossSection' : 9171.0,
##    'process' : "QCD"
##  },
##  'dy_mass10to50' : {
##    'inputFilePath' : "",
##    'crossSection' : 16880.0,
##    'process' : "DY"
##  },
##  'dy_massGt50' : {
##    'inputFilePath' : "",
##    'crossSection' : 5795.0,
##    'process' : "DY"
##  },
##  'w' : {
##    'inputFilePath' : "",
##    'crossSection' : 56990.0,
##    'process' : "W"
##  },
}

version = "2020Jun31"

def runCommand(command):
  #print("executing command = '%s'" % command)
  os.system(command)

def buildConfigFile(cfgFile_original, cfgFile_modified, inputFileNames, process, lumiScale, srcVertices, outputFileName):
  rmCommand   = 'rm -f %s' % cfgFile_modified
  runCommand(rmCommand)
 
  sedCommand  = 'sed'
  sedCommand += ' "s/##inputFilePath/inputFilePath/; s/\$inputFilePath/None/;'
  sedCommand += '  s/##inputFileNames/inputFileNames/; s/\$inputFileNames/%s/;' % [ inputFileName.replace("/", "\/") for inputFileName in inputFileNames ]
  sedCommand += '  s/##processName/processName/; s/\$processName/%s/;' % process
  sedCommand += '  s/##lumiScale/lumiScale/; s/\$lumiScale/%s/;' % lumiScale
  sedCommand += '  s/##srcVertices/srcVertices/; s/\$srcVertices/%s/;' % srcVertices
  sedCommand += '  s/##outputFileName/outputFileName/; s/\$outputFileName/%s/"' % outputFileName
  sedCommand += ' %s > %s' % (cfgFile_original, cfgFile_modified)
  runCommand(sedCommand)
 
makeCommands = []
outputFileNames = []
for sampleName, sample in background_samples.items(): 
  for srcVertices in [ "offlinePrimaryVertices", "hltPhase2PixelVertices" ]:
    print("processing sample = '%s': srcVertices = '%s'" % (sampleName, srcVertices)) 
    inputFilePath = sample['inputFilePath'][srcVertices]
    print(" inputFilePath = '%s'" % inputFilePath)
    inputFileNames = getInputFileNames.getInputFileNames(inputFilePath)
    numInputFiles = len(inputFileNames)
    print("Found %i input files." % numInputFiles)
    numJobs = sample['numJobs']
    lumiScale = sample['crossSection']*conversionFactor*instLuminosity/numJobs # rate in Hz corresponding to one MC event
    # CV: lumiScale needs to be divided by numJobs,
    #     as rate normalization implemented in RecoPFTauAnalyzerBackground::endJob method assumes that whole MC sample is analyzed using single cmsRun job
    print(" lumiScale = %1.2f" % lumiScale)
    for jobId in range(numJobs):
      idxFirstFile = jobId*numInputFiles/numJobs
      idxLastFile = (jobId + 1)*numInputFiles/numJobs - 1
      inputFileNames_job = inputFileNames[idxFirstFile:idxLastFile + 1]
      #print("job #%i: inputFiles = %s" % (jobId, inputFileNames_job))
      cfgFileName_modified = "analyzePFTaus_background_%s_%s_%i_cfg.py" % (sampleName, srcVertices, jobId)
      outputFileName = "analyzePFTaus_background_%s_%s_%s_%i.root" % (sampleName, srcVertices, version, jobId)
      outputFileNames.append(outputFileName)
      buildConfigFile("analyzePFTaus_background_cfg.py", cfgFileName_modified, inputFileNames_job, sample['process'], lumiScale, srcVertices, outputFileName)
      logFileName = cfgFileName_modified.replace("_cfg.py", ".log")
      makeCommands.append("%s:" % outputFileName)
      rmCommand = 'rm -f %s' % logFileName
      makeCommands.append("\t%s" % rmCommand)
      cmsRunCommand = 'cmsRun %s >& %s' % (cfgFileName_modified, logFileName)
      makeCommands.append("\t%s" % cmsRunCommand)

outputFileName_hadd = "analyzePFTaus_background_all_%s.root" % version
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
makeFileName = "Makefile_runJobs_rate"
makeFile = open(makeFileName, "w") 
makeFile.write("all: %s\n" % outputFileName_hadd)
for makeCommand in makeCommands:
    makeFile.write("%s\n" % makeCommand)
makeFile.close()

print("Finished building config files. Now execute 'make -j 16 -f %s'." % makeFileName)
