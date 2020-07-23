#!/usr/bin/env python

import getpass
import os

from HLTrigger.TallinnHLTPFTauAnalyzer.tools.jobTools import getInputFileNames, build_sbatchManagerFile, build_Makefile

signal_samples = {
  'qqH_htt' : {
    'samples' : {
      'offlinePrimaryVertices' : { 
        'inputFilePath' : '/hdfs/cms/store/user/rdewanje/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5/HLTConfig_VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5_wOfflineVtx_wL1_2FM/',
        'numEvents' : 140108
      },
      'hltPhase2PixelVertices' : {
        'inputFilePath' : '/hdfs/cms/store/user/rdewanje/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5/HLTConfig_VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5_wOnlineVtx_wL1_2FM/',
        'numEvents' : 144510
      }
    },
    'numJobs' : 10,
    'process' : "qqH_htt"
  }
}

run_hlt_algorithms = [ "hps" ]
run_hlt_srcVertices = [ "offlinePrimaryVertices", "hltPhase2PixelVertices" ]
run_hlt_isolation_maxDeltaZOptions = [ "primaryVertex", "leadTrack" ]
##run_hlt_isolation_minTrackHits = [ 3, 5, 8 ]
run_hlt_isolation_minTrackHits = [ 8 ]
l1_useStrips = True
##cfgFileName_original = "analyzePFTaus_signal_cfg.py"
cfgFileName_original = "analyzePFTaus_and_L1HPSPFTaus_signal_cfg.py"

version = "2020Jul23"

configDir  = os.path.join("/home",       getpass.getuser(), "Phase2HLT/efficiency", version)
outputDir  = os.path.join("/hdfs/local", getpass.getuser(), "Phase2HLT/efficiency", version)
workingDir = os.getcwd()
cmsswDir   = os.getenv('CMSSW_BASE')

def run_command(command):
  #print("executing command = '%s'" % command)
  os.system(command)

run_command('mkdir -p %s' % configDir)
run_command('mkdir -p %s' % outputDir)

def build_cfgFile(cfgFile_original, cfgFile_modified, 
                  inputFileNames, process, 
                  hlt_srcVertices, hlt_algorithm, hlt_isolation_maxDeltaZOption, hlt_isolation_minTrackHits, l1_useStrips,
                  outputFileName):
  print("Building configFile = '%s'" % cfgFileName_modified)

  rmCommand   = 'rm -f %s' % cfgFileName_modified
  run_command(rmCommand)
 
  sedCommand  = 'sed'
  sedCommand += ' "s/##inputFilePath/inputFilePath/; s/\$inputFilePath/None/;'
  sedCommand += '  s/##inputFileNames/inputFileNames/; s/\$inputFileNames/%s/;' % [ inputFileName.replace("/", "\/") for inputFileName in inputFileNames ]
  sedCommand += '  s/##processName/processName/; s/\$processName/%s/;' % process
  sedCommand += '  s/##hlt_srcVertices/hlt_srcVertices/; s/\$hlt_srcVertices/%s/;' % hlt_srcVertices
  sedCommand += '  s/##hlt_algorithms/hlt_algorithms/; s/\$hlt_algorithm/%s/;' % hlt_algorithm
  sedCommand += '  s/##hlt_isolation_maxDeltaZOptions/hlt_isolation_maxDeltaZOptions/; s/\$hlt_isolation_maxDeltaZOption/%s/;' % hlt_isolation_maxDeltaZOption
  sedCommand += '  s/##hlt_isolation_minTrackHits/hlt_isolation_minTrackHits/; s/\$hlt_isolation_minTrackHits/%s/;' % hlt_isolation_minTrackHits
  sedCommand += '  s/##l1_useStrips/l1_useStrips/; s/\$l1_useStrips/%s/;' % l1_useStrips
  sedCommand += '  s/##outputFileName/outputFileName/; s/\$outputFileName/%s/"' % outputFileName
  sedCommand += ' %s > %s' % (cfgFile_original, cfgFile_modified)
  run_command(sedCommand)
 
jobOptions = {} # key = process + algorithm + isolation_maxDeltaZOption + isolation_minTrackHits (all separated by underscore)
for sampleName, sample in signal_samples.items():
  process = sample['process']
  for hlt_srcVertices in run_hlt_srcVertices:
    print("processing sample = '%s': hlt_srcVertices = '%s'" % (sampleName, hlt_srcVertices)) 
    inputFilePath = sample['samples'][hlt_srcVertices]['inputFilePath']
    print(" inputFilePath = '%s'" % inputFilePath)
    inputFileNames = None
    numInputFiles = None
    if os.path.exists(inputFilePath):
      inputFileNames = getInputFileNames(inputFilePath)
      numInputFiles = len(inputFileNames)
      print("Found %i input files." % numInputFiles)
    else:
      print("Path = '%s' does not exist --> skipping !!" % inputFilePath)
      continue
    numJobs = sample['numJobs']
    for jobId in range(numJobs):
      idxFirstFile = jobId*numInputFiles/numJobs
      idxLastFile = (jobId + 1)*numInputFiles/numJobs - 1
      inputFileNames_job = inputFileNames[idxFirstFile:idxLastFile + 1]
      #print("job #%i: inputFiles = %s" % (jobId, inputFileNames_job))
      for hlt_algorithm in run_hlt_algorithms:
        for hlt_isolation_maxDeltaZOption in run_hlt_isolation_maxDeltaZOptions:
          for hlt_isolation_minTrackHits in run_hlt_isolation_minTrackHits:
            job_key = '%s_%s_%s_%s_%iHits' % (process, hlt_algorithm, hlt_srcVertices, hlt_isolation_maxDeltaZOption, hlt_isolation_minTrackHits)
            if not job_key in jobOptions.keys():
              jobOptions[job_key] = []        
            cfgFileName_modified = os.path.join(configDir, "analyzePFTaus_signal_%s_%s_%s_dz_wrt_%s_%iHits_%i_cfg.py" % \
              (sampleName, hlt_algorithm, hlt_srcVertices, hlt_isolation_maxDeltaZOption, hlt_isolation_minTrackHits, jobId))
            outputFileName = "analyzePFTaus_signal_%s_%s_%s_dz_wrt_%s_%iHits_%i.root" % \
              (sampleName, hlt_algorithm, hlt_srcVertices, hlt_isolation_maxDeltaZOption, hlt_isolation_minTrackHits, jobId)
            build_cfgFile(
              cfgFileName_original, cfgFileName_modified, 
              inputFileNames_job, sample['process'], 
              hlt_srcVertices, hlt_algorithm, hlt_isolation_maxDeltaZOption, hlt_isolation_minTrackHits, l1_useStrips,
              outputFileName)
            logFileName = cfgFileName_modified.replace("_cfg.py", ".log")
            jobOptions[job_key].append({
              'inputFileNames' : inputFileNames_job,
              'cfgFileName'    : cfgFileName_modified,
              'outputFilePath' : outputDir,
              'outputFileName' : outputFileName,
              'logFileName'    : logFileName,
            })

sbatchManagerFileName = os.path.join(configDir, "sbatch_analyzePFTaus_signal.py")
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
  commands.append('mv %s %s' % (outputFileName, os.path.join(outputDir, outputFileName)))
  jobOptions_Makefile_hadd.append({
    'target'          : os.path.join(outputDir, outputFileName),
    'dependencies'    : inputFileNames,
    'commands'        : commands,
    'outputFileNames' : [ os.path.join(outputDir, outputFileName) ],
  })
for sampleName, sample in signal_samples.items(): 
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
  commands.append('sleep 30s')
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
