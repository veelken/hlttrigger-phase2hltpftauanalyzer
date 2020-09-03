#!/usr/bin/env python

import getpass
import os
import time
import threading
import Queue
import json

from HLTrigger.TallinnHLTPFTauAnalyzer.tools.jobTools import getInputFileNames
from HLTrigger.TallinnHLTPFTauAnalyzer.tools.safe_root import ROOT

inputFilePaths = {
  'VBFHToTauTau_wOfflineVtx'        : '/hdfs/cms/store/user/rdewanje/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5/HLTConfig_VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5_wOfflineVtx_wDeepTau4/',
  'MinBias_wOfflineVtx'             : '/hdfs/cms/store/user/rdewanje/MinBias_TuneCP5_14TeV-pythia8/HLTConfig_MinBias_TuneCP5_14TeV-pythia8_wOfflineVtx_wDeepTau4/',
  'QCD_Pt_30to50_wOfflineVtx'       : '/hdfs/cms/store/user/rdewanje/QCD_Pt_30to50_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_30to50_TuneCP5_14TeV_pythia8_wOfflineVtx_wDeepTau4/',
  'QCD_Pt_50to80_wOfflineVtx'       : '/hdfs/cms/store/user/rdewanje/QCD_Pt_50to80_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_50to80_TuneCP5_14TeV_pythia8_wOfflineVtx_wDeepTau4/',
  'QCD_Pt_80to120_wOfflineVtx'      : '/hdfs/cms/store/user/rdewanje/QCD_Pt_80to120_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_80to120_TuneCP5_14TeV_pythia8_wOfflineVtx_wDeepTau4/',
  'QCD_Pt_120to170_wOfflineVtx'     : '/hdfs/cms/store/user/rdewanje/QCD_Pt_120to170_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_120to170_TuneCP5_14TeV_pythia8_wOfflineVtx_wDeepTau4/',
  'QCD_Pt_170to300_wOfflineVtx'     : '/hdfs/cms/store/user/rdewanje/QCD_Pt_170to300_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_170to300_TuneCP5_14TeV_pythia8_wOfflineVtx_wDeepTau4/',
  'QCD_Pt_300to470_wOfflineVtx'     : '/hdfs/cms/store/user/rdewanje/QCD_Pt_300to470_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_300to470_TuneCP5_14TeV_pythia8_wOfflineVtx_wDeepTau4/',
  'DYJetsToLL_M-10to50_wOfflineVtx' : '/hdfs/cms/store/user/rdewanje/DYJetsToLL_M-10to50_TuneCP5_14TeV-madgraphMLM-pythia8/HLTConfig_DYJetsToLL_M-10to50_TuneCP5_14TeV-madgraphMLM-pythia8_wOfflineVtx_wDeepTau4/',
  'DYToLL_M-50_wOfflineVtx'         : '/hdfs/cms/store/user/rdewanje/DYToLL_M-50_TuneCP5_14TeV-pythia8/HLTConfig_DYToLL_M-50_TuneCP5_14TeV-pythia8_wOfflineVtx_wDeepTau4/',
  'WJetsToLNu_wOfflineVtx'          : '/hdfs/cms/store/user/rdewanje/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/HLTConfig_WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8_wOfflineVtx_wDeepTau4'
}

events_treeName = 'Events'

#--------------------------------------------------------------------------------
# CV: code for multi-threading based on example code taken from
#       https://www.tutorialspoint.com/python/python_multithreading.htm
#     and
#       https://docs.python.org/2/library/queue.html
#--------------------------------------------------------------------------------

numThreads = 16

queueLock = threading.Lock()
workQueue = Queue.Queue()

class WorkerThread(threading.Thread):
    def __init__(self, threadID, threadName, workQueue):
        threading.Thread.__init__(self)
        self.threadID   = threadID
        self.threadName = threadName
        self.workQueue  = workQueue
        self.numEvents  = {} # key = sampleName
        self.start_time = {} # key = sampleName
        self.end_time   = {} # key = sampleName
    def run(self): 
        queueLock.acquire()
        ( sampleName, inputFilePath ) = workQueue.get()
        queueLock.release()
        print("%s: %s" % (sampleName, time.ctime(time.time())))
        print("%s: Searching for input files in path = '%s'" % (sampleName, inputFilePath))
        numEvents = 0
        start_time = time.time() 
        end_time = None
        if os.path.exists(inputFilePath):
            inputFileNames = getInputFileNames(inputFilePath)
            print("%s: Found %i input files." % (sampleName, len(inputFileNames)))
            for idx, inputFileName in enumerate(inputFileNames):
                if (idx + 1) % 100 == 0:
                    print("%s: Processing %ith input file..." % (sampleName, (idx + 1)))
                inputFile = ROOT.TFile.Open(inputFileName, 'read')
                if not inputFile or inputFile.IsZombie():
                    print("%s: Failed to open input file = '%s' --> skipping !!" % (sampleName, inputFileName))
                    continue
                events_tree = inputFile.Get(events_treeName)
                numEvents += events_tree.GetEntries()      
            end_time = time.time()
            print("%s: %s" % (sampleName, time.ctime(time.time())))
            print("%s: Finished processing %i input files." % (sampleName, len(inputFileNames)))
            print(" (numEvents = %i, run-time = %s seconds)" % (numEvents, end_time - start_time))
        else:
            print("%s: Path = '%s' does not exist --> skipping !!" % (sampleName, inputFilePath))
            end_time = time.time()
        self.numEvents[sampleName]  = numEvents
        self.start_time[sampleName] = start_time
        self.end_time[sampleName]   = end_time
        workQueue.task_done()
    def getNumEvents(self):
        return self.numEvents

# Create the threads
threads = []
for threadID in range(1, numThreads + 1):
    thread = WorkerThread(threadID, "Thread-%i" % threadID, workQueue)
    thread.daemon = True
    threads.append(thread)

# Fill the queue
queueLock.acquire()
for sampleName, inputFilePath in inputFilePaths.items():
    workQueue.put(( sampleName, inputFilePath ))
queueLock.release()

for thread in threads:
    thread.start()

# Wait for all threads to complete
workQueue.join()

# Collect results
numEventsDict = {}
for thread in threads:
    numEventsDict.update(thread.getNumEvents()) 

# Print results
for sampleName, inputFilePath in inputFilePaths.items():
    numEvents = "N/A"
    if sampleName in numEventsDict.keys():
        numEvents = numEventsDict[sampleName]
    print("Path = '%s': numEvents = %s" % (inputFilePath, numEvents))

jsonFile = open('countEvents.json', 'w')
json.dump(numEventsDict, jsonFile)
jsonFile.close()
