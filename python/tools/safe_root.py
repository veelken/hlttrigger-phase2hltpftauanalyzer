#----------------------------------------------------------------------------------------------------
# CV: The code in this file has been copied from 
#       https://github.com/HEP-KBFI/tth-htt/blob/master/python/safe_root.py
#----------------------------------------------------------------------------------------------------

import sys
argv = sys.argv
sys.argv = []
import ROOT
ROOT.gSystem.ResetSignals()
sys.argv = argv
