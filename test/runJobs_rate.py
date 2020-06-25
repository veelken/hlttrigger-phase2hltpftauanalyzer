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
    'inputFilePath' : "/hdfs/cms/store/user/rdewanje/MinBias_TuneCP5_14TeV-pythia8/HLTConfig_Christian_MINBIAS_Phase2HLTTDRWinter20_PU200_CMSSW_11_1_0_pre6_numCores8_maxMem16kMB_T2_EE_Estonia_blacklist/200612_212910/",
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

version = "2020Jun24"

def runCommand(command):
  print("executing command = '%s'" % command)
  os.system(command)

def buildConfigFile(cfgFile_original, cfgFile_modified, inputFilePath, process, lumiScale, outputFileName):
  rmCommand   = 'rm -f %s' % cfgFile_modified
  runCommand(rmCommand)
 
  sedCommand  = 'sed'
  sedCommand += ' "s/##inputFilePath/inputFilePath/; s/\$inputFilePath/%s/;' % inputFilePath.replace("/", "\/")
  sedCommand += '  s/##processName/processName/; s/\$processName/%s/;' % process
  sedCommand += '  s/##lumiScale/lumiScale/; s/\$lumiScale/%s/;' % lumiScale
  sedCommand += '  s/##outputFileName/outputFileName/; s/\$outputFileName/%s/"' % outputFileName
  sedCommand += ' %s > %s' % (cfgFile_original, cfgFile_modified)
  runCommand(sedCommand)
 
shellCommands = []
outputFileNames = []
for sampleName, sample in background_samples.items():
  print("processing sample = '%s'" % sampleName)
  inputFilePath = sample['inputFilePath']
  print(" inputFilePath = '%s'" % inputFilePath)
  lumiScale = sample['crossSection']*conversionFactor*instLuminosity # rate in Hz corresponding to one MC event
  print(" lumiScale = %1.2f" % lumiScale)
  cfgFileName_modified = "analyzePFTaus_background_%s_cfg.py" % sampleName
  outputFileName = "analyzePFTaus_background_%s_%s.root" % (sampleName, version)
  outputFileNames.append(outputFileName)
  buildConfigFile("analyzePFTaus_background_cfg.py", cfgFileName_modified, inputFilePath, sample['process'], lumiScale, outputFileName)
  logFileName = cfgFileName_modified.replace("_cfg.py", ".log")
  rmCommand = 'rm -f %s' % logFileName
  shellCommands.append(rmCommand)
  cmsRunCommand = 'cmsRun %s >& %s' % (cfgFileName_modified, logFileName)
  shellCommands.append(cmsRunCommand)

outputFileName_hadd = "analyzePFTaus_background_all_%s.root" % version
rmCommand = 'rm -f %s' % outputFileName_hadd
shellCommands.append(rmCommand)
if len(outputFileNames) >= 2:
  haddCommand = 'hadd %s %s' % (outputFileName_hadd, " ".join(outputFileNames))
  shellCommands.append(haddCommand)
elif len(outputFileNames) == 1:
  cpCommand = 'cp %s %s' % (outputFileNames[0], outputFileName_hadd)
  shellCommands.append(cpCommand)
else:
  raise ValueError("No samples defined, so there is nothing to do !!")

shellScriptFileName = "runJobs_rate.sh"
shellScriptFile = open(shellScriptFileName, "w") 
for shellCommand in shellCommands:
    shellScriptFile.write("%s\n" % shellCommand)
shellScriptFile.close()

print("Finished building config files. Now execute 'source %s'." % shellScriptFileName)
