#!/usr/bin/env python

import getpass
import os

from HLTTrigger.TallinnHLTPFTauAnalyzer.tools.jobTools import getInputFileNames, build_sbatchManagerFile, build_Makefile

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
      'offlinePrimaryVertices' : '/hdfs/cms/store/user/rdewanje/MinBias_TuneCP5_14TeV-pythia8/HLTConfig_w_offlineVtxCollection_HGCalFix_MINBIAS_Phase2HLTTDRWinter20_PU200_CMSSW_11_1_0_pre8_isoFix/200630_131646/',
      'hltPhase2PixelVertices' : '/hdfs/cms/store/user/rdewanje/MinBias_TuneCP5_14TeV-pythia8/HLTConfig_w_onlineVtxCollection_HGCalFix_MINBIAS_Phase2HLTTDRWinter20_PU200_CMSSW_11_1_0_pre8_isoFix/200630_131531/'  
    },
    'numJobs' : 10, 
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

run_algorithms = [ "hps" ]
run_srcVertices = [ "offlinePrimaryVertices", "hltPhase2PixelVertices" ]
run_isolation_maxDeltaZOptions = [ "primaryVertex", "leadTrack" ]
run_isolation_minTrackHits = [ 3, 5, 8 ]

version = "2020Jul01"

configDir  = os.path.join("/home",       getpass.getuser(), "Phase2HLT/rate", version)
outputDir  = os.path.join("/hdfs/local", getpass.getuser(), "Phase2HLT/rate", version)
workingDir = os.getcwd()
cmsswDir   = os.getenv('CMSSW_BASE')

def run_command(command):
  #print("executing command = '%s'" % command)
  os.system(command)

run_command('mkdir -p %s' % configDir)
run_command('mkdir -p %s' % outputDir)

def build_cfgFile(cfgFileName_original, cfgFileName_modified, 
                  inputFileNames, process, lumiScale, 
                  srcVertices, algorithm, isolation_maxDeltaZOption, isolation_minTrackHits, 
                  outputFileName):
  print("Building configFile = '%s'" % cfgFileName_modified)

  rmCommand   = 'rm -f %s' % cfgFileName_modified
  run_command(rmCommand)
 
  sedCommand  = 'sed'
  sedCommand += ' "s/##inputFilePath/inputFilePath/; s/\$inputFilePath/None/;'
  sedCommand += '  s/##inputFileNames/inputFileNames/; s/\$inputFileNames/%s/;' % [ inputFileName.replace("/", "\/") for inputFileName in inputFileNames ]
  sedCommand += '  s/##processName/processName/; s/\$processName/%s/;' % process
  sedCommand += '  s/##lumiScale/lumiScale/; s/\$lumiScale/%s/;' % lumiScale
  sedCommand += '  s/##srcVertices/srcVertices/; s/\$srcVertices/%s/;' % srcVertices
  sedCommand += '  s/##algorithms/algorithms/; s/\$algorithm/%s/;' % algorithm
  sedCommand += '  s/##isolation_maxDeltaZOptions/isolation_maxDeltaZOptions/; s/\$isolation_maxDeltaZOption/%s/;' % isolation_maxDeltaZOption
  sedCommand += '  s/##isolation_minTrackHits/isolation_minTrackHits/; s/\$isolation_minTrackHits/%s/;' % isolation_minTrackHits
  sedCommand += '  s/##outputFileName/outputFileName/; s/\$outputFileName/%s/"' % outputFileName
  sedCommand += ' %s > %s' % (cfgFileName_original, cfgFileName_modified)
  run_command(sedCommand)

jobOptions = {} # key = process + algorithm + isolation_maxDeltaZOption + isolation_minTrackHits (all separated by underscore)
for sampleName, sample in background_samples.items(): 
  process = sample['process']
  for srcVertices in run_srcVertices:
    print("processing sample = '%s': srcVertices = '%s'" % (sampleName, srcVertices)) 
    inputFilePath = sample['inputFilePath'][srcVertices]
    print(" inputFilePath = '%s'" % inputFilePath)
    inputFileNames = getInputFileNames(inputFilePath)
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
      for algorithm in run_algorithms:
        for isolation_maxDeltaZOption in run_isolation_maxDeltaZOptions:
          for isolation_minTrackHits in run_isolation_minTrackHits:
            job_key = '%s_%s_%s_%s_%iHits' % (process, algorithm, srcVertices, isolation_maxDeltaZOption, isolation_minTrackHits)
            if not job_key in jobOptions.keys():
              jobOptions[job_key] = []        
            cfgFileName_modified = os.path.join(configDir, "analyzePFTaus_background_%s_%s_%s_dz_wrt_%s_%iHits_%i_cfg.py" % \
              (sampleName, algorithm, srcVertices, isolation_maxDeltaZOption, isolation_minTrackHits, jobId))
            outputFileName = "analyzePFTaus_background_%s_%s_%s_dz_wrt_%s_%iHits_%i.root" % \
              (sampleName, algorithm, srcVertices, isolation_maxDeltaZOption, isolation_minTrackHits, jobId)
            build_cfgFile(
              "analyzePFTaus_background_cfg.py", cfgFileName_modified, 
              inputFileNames_job, sample['process'], lumiScale, 
              srcVertices, algorithm, isolation_maxDeltaZOption, isolation_minTrackHits, 
              outputFileName)
            logFileName = cfgFileName_modified.replace("_cfg.py", ".log")
            jobOptions[job_key].append({
              'inputFileNames' : inputFileNames_job,
              'cfgFileName'    : cfgFileName_modified,
              'outputFilePath' : outputDir,
              'outputFileName' : outputFileName,
              'logFileName'    : logFileName,
            })

sbatchManagerFileName = os.path.join(configDir, "sbatch_analyzePFTaus_background.py")
jobOptions_sbatchManager = []
for job_key, jobs in jobOptions.items():
  jobOptions_sbatchManager.extend(jobs)
build_sbatchManagerFile(sbatchManagerFileName, jobOptions_sbatchManager, workingDir, cmsswDir, version)

jobOptions_Makefile_sbatch = []
jobOptions_Makefile_sbatch.append({
  'target'          : "phony",
  'dependencies'    : [],
  'commands'        : [ 'python %s' % sbatchManagerFileName ],
  'outputFileNames' : [ os.path.join(job['outputFilePath'], job['outputFileName']) for job in jobOptions_sbatchManager ],
})
makeFileName_sbatch = os.path.join(configDir, "Makefile_sbatch")
build_Makefile(makeFileName_sbatch, jobOptions_Makefile_sbatch)

jobOptions_Makefile_hadd = []
for job_key, jobs in jobOptions.items():
  inputFileNames = [ os.path.join(job['outputFilePath'], job['outputFileName']) for job in jobs ]
  outputFileName = "hadd_%s.root" % job_key
  commands = []
  commands.append('rm -f %s' % outputFileName)
  commands.append('hadd %s %s' % (outputFileName, " ".join(inputFileNames)))
  commands.append('cp -f %s %s' % (outputFileName, os.path.join(outputDir, outputFileName)))
  commands.append('sleep 5s')
  commands.append('rm -f %s' % outputFileName)
  jobOptions_Makefile_hadd.append({
    'target'          : os.path.join(outputDir, outputFileName),
    'dependencies'    : inputFileNames,
    'commands'        : commands,
    'outputFileNames' : [ os.path.join(outputDir, outputFileName) ],
  })
for sampleName, sample in background_samples.items(): 
  process = sample['process']
  inputFileNames = []
  for job in jobOptions_Makefile_hadd:
    for outputFileName_job in job['outputFileNames']:
      if outputFileName_job.find(sampleName) != -1:
        inputFileNames.append(outputFileName_job)
  outputFileName = "hadd_%s_all.root" % process
  commands = []
  commands.append('rm -f %s' % outputFileName)
  commands.append('hadd %s %s' % (outputFileName, " ".join(inputFileNames)))
  commands.append('cp -f %s %s' % (outputFileName, os.path.join(outputDir, outputFileName)))
  commands.append('sleep 5s')
  commands.append('rm -f %s' % outputFileName)
  jobOptions_Makefile_hadd.append({
    'target'          : os.path.join(outputDir, outputFileName),
    'dependencies'    : inputFileNames,
    'commands'        : commands,
    'outputFileNames' : [ os.path.join(outputDir, outputFileName) ],
  })
makeFileName_hadd = os.path.join(configDir, "Makefile_hadd")
build_Makefile(makeFileName_hadd, jobOptions_Makefile_hadd)

message  = "Finished building config files."
message += " Now execute 'make -f %s' to submit the jobs to the batch system." % makeFileName_sbatch
message += " Once all batch jobs have finished processing, execute 'make -j 4 -f %s' to merge the output files." % makeFileName_hadd
print(message)
