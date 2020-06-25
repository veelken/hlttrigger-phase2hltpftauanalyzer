import FWCore.ParameterSet.Config as cms

import os
import re

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
