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
