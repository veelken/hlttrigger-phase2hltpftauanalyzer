
def get_suffix(hlt_srcVertices, hlt_isolation_maxDeltaZOption, hlt_isolation_minTrackHits):
    suffix = "%iHits" % hlt_isolation_minTrackHits
    if hlt_isolation_maxDeltaZOption == "primaryVertex":
      suffix += "MaxDeltaZ"
    elif hlt_isolation_maxDeltaZOption == "leadTrack":
      suffix += "MaxDeltaZToLeadTrack"
    else:
      raise ValueError("Invalid parameter hlt_isolation_maxDeltaZOption = '%s' !!" % hlt_isolation_maxDeltaZOption)
    if hlt_srcVertices == "offlinePrimaryVertices":
      suffix += "WithOfflineVertices"
    elif hlt_srcVertices == "hltPhase2PixelVertices":
      suffix += "WithOnlineVertices"
    elif hlt_srcVertices == "hltPhase2TrimmedPixelVertices":
      suffix += "WithOnlineVerticesTrimmed"
    else:
      raise ValueError("Invalid parameter hlt_srcVertices = '%s' !!" % hlt_srcVertices)  
    return suffix
