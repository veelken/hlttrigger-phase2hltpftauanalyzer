#--------------------------------------------------------------------------------
# CV: copied from tthAnalysis/HiggsToTauTau/python/safe_root.py
import sys
argv = sys.argv
sys.argv = []
import ROOT
ROOT.gSystem.ResetSignals()
sys.argv = argv
#--------------------------------------------------------------------------------
