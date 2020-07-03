#----------------------------------------------------------------------------------------------------
# CV: The code in this file has been copied from 
#       https://github.com/HEP-KBFI/tth-htt/tree/master/python
#----------------------------------------------------------------------------------------------------

import logging
import os
import subprocess
import sys

#--------------------------------------------------------------------------------
# CV: copied from tthAnalysis/HiggsToTauTau/python/analysisTools.py
def createFile(fileName, lines, nofNewLines = 2):
    """Auxiliary function to write new config file,
       containg the lines given as argument.
    """
    content = "\n".join(lines)
    content += nofNewLines * "\n"
    with open(fileName, "w") as f:
      f.write(content)
#--------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# CV: copied from tthAnalysis/HiggsToTauTau/python/common.py
logging.basicConfig(
  stream = sys.stdout,
  level  = logging.INFO,
  format = '%(asctime)s - %(levelname)s: %(message)s',
)
#----------------------------------------------------------------------------------------------------

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
