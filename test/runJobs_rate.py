#!/usr/bin/env python

import getpass
import os

from HLTrigger.TallinnHLTPFTauAnalyzer.tools.jobTools import getInputFileNames, build_sbatchManagerFile, build_Makefile

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
##  'minbias' : {
##    'inputFilePath' : {
##      'offlinePrimaryVertices' : '',
##      'hltPhase2PixelVertices' : ''  
##    },
##    'numJobs' : 10, 
##    'crossSection' : 75.0e+9/200, # pb
##    'process' : "minbias"
##  },
  # CV: cross-sections for QCD, Drell-Yan, and W+jets production taken from the twiki
  #       https://twiki.cern.ch/twiki/bin/viewauth/CMS/HighLevelTriggerPhase2#Rate_calculations
##  'qcd_pt15to20' : {
##    'inputFilePath' : {
##      'offlinePrimaryVertices' : '',
##      'hltPhase2PixelVertices' : ''  
##    },
##    'numJobs' : 10, 
##    'crossSection' : 923300000.0,
##    'process' : "QCD"
##  },
##  'qcd_pt20to30' : {
##    'inputFilePath' : {
##      'offlinePrimaryVertices' : '',
##      'hltPhase2PixelVertices' : ''  
##    },
##    'numJobs' : 10, 
##    'crossSection' : 436000000.0,
##    'process' : "QCD"
##  },
  'qcd_pt30to50' : {
    'inputFilePath' : {
      'offlinePrimaryVertices' : '/hdfs/cms/store/user/rdewanje/QCD_Pt_30to50_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_30to50_TuneCP5_14TeV_pythia8_wOfflineVtx/200710_232606/',
      'hltPhase2PixelVertices' : '/hdfs/cms/store/user/rdewanje/QCD_Pt_30to50_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_30to50_TuneCP5_14TeV_pythia8_wOnlineVtx/200710_233400/'  
    },
    'numJobs' : 10, 
    'crossSection' : 118400000.0,
    'process' : "QCD"
  },
  'qcd_pt50to80' : {
    'inputFilePath' : {
      'offlinePrimaryVertices' : '/hdfs/cms/store/user/rdewanje/QCD_Pt_50to80_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_50to80_TuneCP5_14TeV_pythia8_wOfflineVtx/200710_232636/',
      'hltPhase2PixelVertices' : '/hdfs/cms/store/user/rdewanje/QCD_Pt_50to80_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_50to80_TuneCP5_14TeV_pythia8_wOnlineVtx/200710_230414/'  
    },
    'numJobs' : 10, 
    'crossSection' : 17650000.0,
    'process' : "QCD"
  },
  'qcd_pt80to120' : {
    'inputFilePath' : {
      'offlinePrimaryVertices' : '/hdfs/cms/store/user/rdewanje/QCD_Pt_80to120_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_80to120_TuneCP5_14TeV_pythia8_wOfflineVtx/200710_232521/',
      'hltPhase2PixelVertices' : '/hdfs/cms/store/user/rdewanje/QCD_Pt_80to120_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_80to120_TuneCP5_14TeV_pythia8_wOnlineVtx/200710_233333/'  
    },
    'numJobs' : 10, 
    'crossSection' : 2671000.0,
    'process' : "QCD"
  },
  'qcd_pt120to170' : {
    'inputFilePath' : {
      'offlinePrimaryVertices' : '/hdfs/cms/store/user/rdewanje/QCD_Pt_120to170_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_120to170_TuneCP5_14TeV_pythia8_wOfflineVtx/200710_232405/',
      'hltPhase2PixelVertices' : '/hdfs/cms/store/user/rdewanje/QCD_Pt_120to170_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_120to170_TuneCP5_14TeV_pythia8_wOnlineVtx/200710_233305/'  
    },
    'numJobs' : 10, 
    'crossSection' : 469700.0 ,
    'process' : "QCD"
  },
  'qcd_pt170to300' : {
    'inputFilePath' : {
      'offlinePrimaryVertices' : '/hdfs/cms/store/user/rdewanje/QCD_Pt_170to300_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_170to300_TuneCP5_14TeV_pythia8_wOfflineVtx/200710_232451/',
      'hltPhase2PixelVertices' : '/hdfs/cms/store/user/rdewanje/QCD_Pt_170to300_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_170to300_TuneCP5_14TeV_pythia8_wOnlineVtx/200710_233114/'  
    },
    'numJobs' : 10, 
    'crossSection' : 121700.0,
    'process' : "QCD"
  },  
  'qcd_ptGt300' : {
    'inputFilePath' : {
      'offlinePrimaryVertices' : '/hdfs/cms/store/user/rdewanje/QCD_Pt_300to470_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_300to470_TuneCP5_14TeV_pythia8_wOfflineVtx/200710_232252/',
      'hltPhase2PixelVertices' : '/hdfs/cms/store/user/rdewanje/QCD_Pt_300to470_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_300to470_TuneCP5_14TeV_pythia8_wOnlineVtx/200710_233141/'  
    },
    'numJobs' : 10, 
    'crossSection' : 9171.0,
    'process' : "QCD"
  },
  'dy_mass10to50' : {
    'inputFilePath' : {
      'offlinePrimaryVertices' : '/hdfs/cms/store/user/rdewanje/DYJetsToLL_M-10to50_TuneCP5_14TeV-madgraphMLM-pythia8/HLTConfig_DYJetsToLL_M-10to50_TuneCP5_14TeV-madgraphMLM-pythia8_wOfflineVtx/200710_232002/',
      'hltPhase2PixelVertices' : '/hdfs/cms/store/user/rdewanje/DYJetsToLL_M-10to50_TuneCP5_14TeV-madgraphMLM-pythia8/HLTConfig_DYJetsToLL_M-10to50_TuneCP5_14TeV-madgraphMLM-pythia8_wOnlineVtx/200710_232851/'  
    },
    'numJobs' : 10, 
    'crossSection' : 16880.0,
    'process' : "DY"
  },
  'dy_massGt50' : {
    'inputFilePath' : {
      'offlinePrimaryVertices' : '/hdfs/cms/store/user/rdewanje/DYToLL_M-50_TuneCP5_14TeV-pythia8/HLTConfig_DYToLL_M-50_TuneCP5_14TeV-pythia8_wOfflineVtx/200710_231855/',
      'hltPhase2PixelVertices' : '/hdfs/cms/store/user/rdewanje/DYToLL_M-50_TuneCP5_14TeV-pythia8/HLTConfig_DYToLL_M-50_TuneCP5_14TeV-pythia8_wOnlineVtx/200710_232823/'  
    },
    'numJobs' : 10, 
    'crossSection' : 5795.0,
    'process' : "DY"
  },
##  'w' : {
##    'inputFilePath' : {
##      'offlinePrimaryVertices' : '',
##      'hltPhase2PixelVertices' : ''  
##    },
##    'numJobs' : 10, 
##    'crossSection' : 56990.0,
##    'process' : "W"
##  },
}

run_algorithms = [ "hps" ]
run_srcVertices = [ "offlinePrimaryVertices", "hltPhase2PixelVertices" ]
run_isolation_maxDeltaZOptions = [ "primaryVertex", "leadTrack" ]
##run_isolation_minTrackHits = [ 3, 5, 8 ]
run_isolation_minTrackHits = [ 8 ]
l1_useStrips = True
##cfgFileName_original = "analyzePFTaus_signal_cfg.py"
cfgFileName_original = "analyzePFTaus_and_L1HPSPFTaus_signal_cfg.py"

version = "2020Jul17"

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
                  hlt_srcVertices, hlt_algorithm, hlt_isolation_maxDeltaZOption, hlt_isolation_minTrackHits, l1_useStrips, 
                  outputFileName):
  print("Building configFile = '%s'" % cfgFileName_modified)

  rmCommand   = 'rm -f %s' % cfgFileName_modified
  run_command(rmCommand)
 
  sedCommand  = 'sed'
  sedCommand += ' "s/##inputFilePath/inputFilePath/; s/\$inputFilePath/None/;'
  sedCommand += '  s/##inputFileNames/inputFileNames/; s/\$inputFileNames/%s/;' % [ inputFileName.replace("/", "\/") for inputFileName in inputFileNames ]
  sedCommand += '  s/##processName/processName/; s/\$processName/%s/;' % process
  sedCommand += '  s/##lumiScale/lumiScale/; s/\$lumiScale/%s/;' % lumiScale
  sedCommand += '  s/##hlt_srcVertices/hlt_srcVertices/; s/\$hlt_srcVertices/%s/;' % hlt_srcVertices
  sedCommand += '  s/##hlt_algorithms/hlt_algorithms/; s/\$hlt_algorithm/%s/;' % hlt_algorithm
  sedCommand += '  s/##hlt_isolation_maxDeltaZOptions/hlt_isolation_maxDeltaZOptions/; s/\$hlt_isolation_maxDeltaZOption/%s/;' % hlt_isolation_maxDeltaZOption
  sedCommand += '  s/##hlt_isolation_minTrackHits/hlt_isolation_minTrackHits/; s/\$hlt_isolation_minTrackHits/%s/;' % hlt_isolation_minTrackHits
  sedCommand += '  s/##l1_useStrips/l1_useStrips/; s/\$l1_useStrips/%s/;' % l1_useStrips
  sedCommand += '  s/##outputFileName/outputFileName/; s/\$outputFileName/%s/"' % outputFileName
  sedCommand += ' %s > %s' % (cfgFileName_original, cfgFileName_modified)
  run_command(sedCommand)

jobOptions = {} # key = sampleName + hlt_algorithm + hlt_isolation_maxDeltaZOption + hlt_isolation_minTrackHits (all separated by underscore)
for sampleName, sample in background_samples.items(): 
  process = sample['process']
  for hlt_srcVertices in run_hlt_srcVertices:
    print("processing sample = '%s': hlt_srcVertices = '%s'" % (sampleName, hlt_srcVertices)) 
    inputFilePath = sample['inputFilePath'][hlt_srcVertices]
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
      for hlt_algorithm in run_hlt_algorithms:
        for hlt_isolation_maxDeltaZOption in run_hlt_isolation_maxDeltaZOptions:
          for hlt_isolation_minTrackHits in run_hlt_isolation_minTrackHits:
            job_key = '%s_%s_%s_dz_wrt_%s_%iHits' % (sampleName, hlt_algorithm, hlt_srcVertices, hlt_isolation_maxDeltaZOption, hlt_isolation_minTrackHits)
            if not job_key in jobOptions.keys():
              jobOptions[job_key] = []        
            cfgFileName_modified = os.path.join(configDir, "analyzePFTaus_background_%s_%i_cfg.py" % \
              (job_key, jobId))
            outputFileName = "analyzePFTaus_background_%s_%i.root" % \
              (job_key, jobId)
            build_cfgFile(
              cfgFileName_original, cfgFileName_modified, 
              inputFileNames_job, sample['process'], lumiScale, 
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

jobOptions_Makefile_hadd_stage1 = {} # key = outputFileName
for sampleName, sample in background_samples.items(): 
  for job_key, jobs in jobOptions.items():
    if job_key.find(sampleName) != -1:
      inputFileNames = [ os.path.join(job['outputFilePath'], job['outputFileName']) for job in jobs ]
      outputFileName = "hadd_stage1_%s.root" % job_key
      if not outputFileName in jobOptions_Makefile_hadd_stage1.keys():
        commands = []
        commands.append('rm -f %s' % outputFileName)
        commands.append('hadd %s %s' % (outputFileName, " ".join(inputFileNames)))
        commands.append('cp -f %s %s' % (outputFileName, os.path.join(outputDir, outputFileName)))
        commands.append('sleep 5s')
        commands.append('rm -f %s' % outputFileName)
        jobOptions_Makefile_hadd_stage1[outputFileName]= {
          'target'          : os.path.join(outputDir, outputFileName),
          'dependencies'    : inputFileNames,
          'commands'        : commands,
          'outputFileNames' : [ os.path.join(outputDir, outputFileName) ],
        }
jobOptions_Makefile_hadd_stage2 = {} # key = outputFileName
processes = []
for sampleName, sample in background_samples.items():
  if not sample['process'] in processes:
    processes.append(sample['process'])
print("processes = ", processes)
for process in processes:  
  for hlt_srcVertices in run_hlt_srcVertices:
    for hlt_algorithm in run_hlt_algorithms:
      for hlt_isolation_maxDeltaZOption in run_hlt_isolation_maxDeltaZOptions:
        for hlt_isolation_minTrackHits in run_hlt_isolation_minTrackHits:
          inputFileNames = []
          for sampleName, sample in background_samples.items():
            if sample['process'] == process:
              for job in jobOptions_Makefile_hadd_stage1.values():
                for outputFileName_job in job['outputFileNames']:
                  job_key = '%s_%s_%s_dz_wrt_%s_%iHits' % (sampleName, hlt_algorithm, hlt_srcVertices, hlt_isolation_maxDeltaZOption, hlt_isolation_minTrackHits)
                  if outputFileName_job.find(job_key) != -1:
                    inputFileNames.append(outputFileName_job)
          job_key_hadd_stage2 = '%s_%s_%s_dz_wrt_%s_%iHits' % (process, hlt_algorithm, hlt_srcVertices, hlt_isolation_maxDeltaZOption, hlt_isolation_minTrackHits)
          outputFileName = "hadd_stage2_%s.root" % job_key_hadd_stage2
          if not outputFileName in jobOptions_Makefile_hadd_stage2.keys():
            commands = []
            commands.append('rm -f %s' % outputFileName)
            commands.append('hadd %s %s' % (outputFileName, " ".join(inputFileNames)))
            commands.append('cp -f %s %s' % (outputFileName, os.path.join(outputDir, outputFileName)))
            commands.append('sleep 5s')
            commands.append('rm -f %s' % outputFileName)
            jobOptions_Makefile_hadd_stage2[outputFileName]= {
              'target'          : os.path.join(outputDir, outputFileName),
              'dependencies'    : inputFileNames,
              'commands'        : commands,
              'outputFileNames' : [ os.path.join(outputDir, outputFileName) ],
            }
jobOptions_Makefile_hadd_stage3 = {} # key = outputFileName
for sampleName, sample in background_samples.items(): 
  process = sample['process']
  inputFileNames = []
  for job in jobOptions_Makefile_hadd_stage2.values():
    for outputFileName_job in job['outputFileNames']:
      if outputFileName_job.find(process) != -1:
        inputFileNames.append(outputFileName_job)
  outputFileName = "hadd_%s_all.root" % process
  if not outputFileName in jobOptions_Makefile_hadd_stage3.keys():
    commands = []
    commands.append('rm -f %s' % outputFileName)
    commands.append('hadd %s %s' % (outputFileName, " ".join(inputFileNames)))
    commands.append('cp -f %s %s' % (outputFileName, os.path.join(outputDir, outputFileName)))
    commands.append('sleep 30s')
    commands.append('rm -f %s' % outputFileName)
    jobOptions_Makefile_hadd_stage3[outputFileName] = {
      'target'          : os.path.join(outputDir, outputFileName),
      'dependencies'    : inputFileNames,
      'commands'        : commands,
      'outputFileNames' : [ os.path.join(outputDir, outputFileName) ],
    }
makeFileName_hadd = os.path.join(configDir, "Makefile_hadd")
jobOptions_Makefile_hadd = []
jobOptions_Makefile_hadd.extend([ job for job in jobOptions_Makefile_hadd_stage1.values() ])
jobOptions_Makefile_hadd.extend([ job for job in jobOptions_Makefile_hadd_stage2.values() ])
jobOptions_Makefile_hadd.extend([ job for job in jobOptions_Makefile_hadd_stage3.values() ])
build_Makefile(makeFileName_hadd, jobOptions_Makefile_hadd)

message  = "Finished building config files."
message += " Now execute 'make -f %s' to submit the jobs to the batch system." % makeFileName_sbatch
message += " Once all batch jobs have finished processing, execute 'make -j 4 -f %s' to merge the output files." % makeFileName_hadd
print(message)
