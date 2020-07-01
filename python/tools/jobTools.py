from HLTTrigger.TallinnHLTPFTauAnalyzer.tools.common import logging

import os
import re
import subprocess

inputFile_regex = r"step3_RAW2DIGI_RECO_[a-zA-Z0-9-_]+.root"
inputFile_matcher = re.compile(inputFile_regex)

def getInputFileNames(inputFilePath):
    inputFileNames = []
    files = os.listdir(inputFilePath)
    for file in files:
        if os.path.isdir(os.path.join(inputFilePath, file)):
            inputFileNames.extend(getInputFileNames(os.path.join(inputFilePath, file)))
        else:
            # check if name of inputFile matches regular expression
            if inputFile_matcher.match(file):
                inputFileNames.append("file:%s" % os.path.join(inputFilePath, file))
    return inputFileNames

#----------------------------------------------------------------------------------------------------
# CV: copied from tthAnalysis/HiggsToTauTau/python/jobTools.py
def run_cmd(command, do_not_log = False, stdout_file = None, stderr_file = None,
            return_stderr = False):
  """Runs given commands and logs stdout and stderr to files
  """
  p = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
  stdout, stderr = p.communicate()
  # Remove trailing newline
  stdout = stdout.rstrip('\n')
  stderr = stderr.rstrip('\n')

  if stdout_file:
    stdout_file.write(command + "\n")
    stdout_file.write('%s\n' % stdout)
  if stderr_file:
    stderr_file.write('%s\n' % stderr)

  if not do_not_log:
    logging.debug("Executed command: '%s'" % command)
    logging.debug("stdout: '%s'" % stdout)
    logging.debug("stderr: '%s'" % stderr)

  if return_stderr:
    return stdout, stderr
  return stdout

def get_log_version(list_of_log_files):
  """
  :param list_of_log_files: tuple of strings, List of log files that need to have the same version number
  :return: tuple of strings, List of log files with the same version number
  Instead of passing log files one-by-one, the more consistent way of dealing with these things is
  to loop over a set of log files that represents a single joint iteration of the jobs.
  """
  if all(map(lambda path: not os.path.exists(path), list_of_log_files)):
    # if none of the files exist, then retain the path names
    return list_of_log_files
  # loop over version numbers until none of the paths exist
  version_idx = 1
  while True:
    list_of_log_files_versioned = tuple(map(lambda path: "%s.%i" % (path, version_idx), list_of_log_files))
    if all(map(lambda path: not os.path.exists(path), list_of_log_files_versioned)):
      return list_of_log_files_versioned
    else:
      # some log files already exist -> increase the version number
      version_idx += 1
#----------------------------------------------------------------------------------------------------

def build_sbatchManagerFile(sbatchManagerFileName, jobOptions, workingDir, cmsswDir, version):
  print("Building sbatchManagerFile = '%s'" % sbatchManagerFileName)

  sbatchManagerName = "Phase2HLT_rate_%s" % version

  lines_sbatchManager = []
  lines_sbatchManager.append("from HLTTrigger.TallinnHLTPFTauAnalyzer.tools.sbatchManager import sbatchManager")
  lines_sbatchManager.append("m = sbatchManager('%s', verbose = False, dry_run = False, use_home = False, min_file_size = 20000, max_num_submittedJobs = 5000)" % sbatchManagerName)
  lines_sbatchManager.append("m.setWorkingDir('%s')" % workingDir)
  lines_sbatchManager.append("m.setcmssw_base_dir('%s')" % cmsswDir)
  lines_sbatchManager.append("m.log_completion = False")
  for job in jobOptions:
    lines_sbatchManager.append("m.submitJob(")
    lines_sbatchManager.append("  inputFiles             = %s," % job['inputFileNames'])
    lines_sbatchManager.append("  executable             = 'cmsRun',")
    lines_sbatchManager.append("  command_line_parameter = '%s'," % job['cfgFileName'])
    lines_sbatchManager.append("  outputFilePath         = '%s'," % job['outputFilePath'])
    lines_sbatchManager.append("  outputFiles            = [ '%s' ]," % job['outputFileName'])
    lines_sbatchManager.append("  scriptFile             = '%s'," % job['logFileName'].replace(".log", ".sh"))
    lines_sbatchManager.append("  logFile                = '%s'," % job['logFileName'])
    lines_sbatchManager.append("  skipIfOutputFileExists = False,")
    lines_sbatchManager.append("  job_template_file      = '%s'," % os.path.join(cmsswDir, "src/HLTTrigger/TallinnHLTPFTauAnalyzer/python/templates/sbatch-node.sh.template"))
    lines_sbatchManager.append(")")
  lines_sbatchManager.append("m.waitForJobs()")
  sbatchManagerFile = open(sbatchManagerFileName, "w") 
  for line in lines_sbatchManager:
    sbatchManagerFile.write("%s\n" % line)
  sbatchManagerFile.close()

def build_Makefile(makeFileName, jobOptions):
  print("Building Makefile = '%s'" % makeFileName)

  lines_Makefile = []
  lines_Makefile.append(".DEFAULT_GOAL := all")
  lines_Makefile.append("SHELL := /bin/bash")
  lines_Makefile.append("all: %s" % " ".join([ job['target'] for job in jobOptions ]))
  for job in jobOptions:
    lines_Makefile.append("%s: %s" % (job['target'], " ".join(job['dependencies'])))
    for command in job['commands']:
      lines_Makefile.append("\t%s" % command)
  lines_Makefile.append("clean:")
  outputFileNames = []
  for job in jobOptions:
    outputFileNames.extend(job['outputFileNames'])
  for outputFileName in outputFileNames:
    lines_Makefile.append("\trm -f %s" % outputFileName)

  makeFile = open(makeFileName, "w") 
  for line in lines_Makefile:
    makeFile.write("%s\n" % line)
  makeFile.close()
 
