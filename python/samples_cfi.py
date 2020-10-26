
signal_samples = {
  'qqH_htt' : {
    'samples' : {
      'offlinePrimaryVertices' : { 
        'inputFilePath' : '/hdfs/cms/store/user/sbhowmik/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5/VBFHToTauTau_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1113_20201009/',
        'numEvents' : 296308
      },
      'hltPhase2PixelVertices' : {
        'inputFilePath' : '/hdfs/cms/store/user/rdewanje/VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5/HLTConfig_VBFHToTauTau_M125_14TeV_powheg_pythia8_correctedGridpack_tuneCP5_wOnlineVtx_wL1_2FM/',
        'numEvents' : 300000
      }
    },
    'numJobs' : 10,
    'process' : "qqH_htt"
  }
}

background_samples = {
  # CV: minbias cross-section taken to be 75mb, 
  #     resulting in 200 pileup interaction per bunch-crossing for instantaneous luminosity of 7.5e+34 cm^-2 s^-1
  #    (assuming that 70% of all bunches are colliding, i.e. collision frequence of 28 MHz)
  'minbias' : {
    'samples' : {
      'offlinePrimaryVertices' : { 
        'inputFilePath' : '/hdfs/cms/store/user/sbhowmik/MinBias_TuneCP5_14TeV-pythia8/MinBias_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1113_20201009/',
        'numEvents' : 981592
      },
      'hltPhase2PixelVertices' : {
        'inputFilePath' : '/hdfs/cms/store/user/rdewanje/MinBias_TuneCP5_14TeV-pythia8/HLTConfig_MinBias_TuneCP5_14TeV-pythia8_wOnlineVtx_wL1/',
        'numEvents' : 1000000
      }
    },
    'numJobs' : 10, 
    'crossSection' : 80.0e+9/200, # pb
    'process' : "minbias"
  },
  # CV: cross-sections for QCD, Drell-Yan, and W+jets production taken from the twiki
  #       https://twiki.cern.ch/twiki/bin/viewauth/CMS/HighLevelTriggerPhase2#Rate_calculations
##  'qcd_pt15to20' : {
##    'samples' : {
##      'offlinePrimaryVertices' : { 
##        'inputFilePath' : 
##        'numEvents' : 
##      },
##      'hltPhase2PixelVertices' : {
##        'inputFilePath' : 
##        'numEvents' : 
##      }
##    },
##    'numJobs' : 10, 
##    'crossSection' : 923300000.0/200,
##    'process' : "QCD"
##  },
##  'qcd_pt20to30' : {
##    'samples' : {
##      'offlinePrimaryVertices' : { 
##        'inputFilePath' : 
##        'numEvents' : 
##      },
##      'hltPhase2PixelVertices' : {
##        'inputFilePath' : 
##        'numEvents' : 
##      }
##    },
##    'numJobs' : 10, 
##    'crossSection' : 436000000.0/200,
##    'process' : "QCD"
##  },
  'qcd_pt30to50' : {
    'samples' : {
      'offlinePrimaryVertices' : { 
        'inputFilePath' : '/hdfs/cms/store/user/sbhowmik/QCD_Pt_30to50_TuneCP5_14TeV_pythia8/QCDPt30to50_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1113_20201009/',
        'numEvents' : 483498
      },
      'hltPhase2PixelVertices' : {
        'inputFilePath' : '/hdfs/cms/store/user/rdewanje/QCD_Pt_30to50_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_30to50_TuneCP5_14TeV_pythia8_wOnlineVtx_wL1/',
        'numEvents' : 487066
      }
    },
    'numJobs' : 10, 
    'crossSection' : 118400000.0/200,
    ##'kFactor_pT_hat' : 1.15,
    'kFactor_pT_hat' : 1.0,
    'process' : "QCD"
  },
  'qcd_pt50to80' : {
    'samples' : {
      'offlinePrimaryVertices' : { 
        'inputFilePath' : '/hdfs/cms/store/user/sbhowmik/QCD_Pt_50to80_TuneCP5_14TeV_pythia8/QCDPt50to80_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1113_20201009/',
        'numEvents' : 300000
      },
      'hltPhase2PixelVertices' : {
        'inputFilePath' : '/hdfs/cms/store/user/rdewanje/QCD_Pt_50to80_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_50to80_TuneCP5_14TeV_pythia8_wOnlineVtx_wL1/',
        'numEvents' : 235400
      }
    },
    'numJobs' : 10, 
    'crossSection' : 17650000.0/200,
    ##'kFactor_pT_hat' : 1.15,
    'kFactor_pT_hat' : 1.0,
    'process' : "QCD"
  },
  'qcd_pt80to120' : {
    'samples' : {
      'offlinePrimaryVertices' : { 
        'inputFilePath' : '/hdfs/cms/store/user/sbhowmik/QCD_Pt_80to120_TuneCP5_14TeV_pythia8/QCDPt80to120_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1113_20201009/',
        'numEvents' : 100000
      },
      'hltPhase2PixelVertices' : {
        'inputFilePath' : '/hdfs/cms/store/user/rdewanje/QCD_Pt_80to120_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_80to120_TuneCP5_14TeV_pythia8_wOnlineVtx_wL1/',
        'numEvents' : 100000
      }
    },
    'numJobs' : 10, 
    'crossSection' : 2671000.0/200,
    ##'kFactor_pT_hat' : 1.15,
    'kFactor_pT_hat' : 1.0,
    'process' : "QCD"
  },
  'qcd_pt120to170' : {
    'samples' : {
      'offlinePrimaryVertices' : { 
        'inputFilePath' : '/hdfs/cms/store/user/sbhowmik/QCD_Pt_120to170_TuneCP5_14TeV_pythia8/QCDPt120to170_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1113_20201009/',
        'numEvents' : 49601
      },
      'hltPhase2PixelVertices' : {
        'inputFilePath' : '/hdfs/cms/store/user/rdewanje/QCD_Pt_120to170_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_120to170_TuneCP5_14TeV_pythia8_wOnlineVtx_wL1/',
        'numEvents' : 50000
      }
    },
    'numJobs' : 10, 
    'crossSection' : 469700.0/200,
    ##'kFactor_pT_hat' : 1.10,
    'kFactor_pT_hat' : 1.0,
    'process' : "QCD"
  },
  'qcd_pt170to300' : {
    'samples' : {
      'offlinePrimaryVertices' : { 
        'inputFilePath' : '/hdfs/cms/store/user/sbhowmik/QCD_Pt_170to300_TuneCP5_14TeV_pythia8/QCDPt170to300_Phase2HLTTDRSummer20ReRECOMiniAOD_CMSSW_1113_20201009/',
        'numEvents' : 50000
      },
      'hltPhase2PixelVertices' : {
        'inputFilePath' : '/hdfs/cms/store/user/rdewanje/QCD_Pt_170to300_TuneCP5_14TeV_pythia8/HLTConfig_QCD_Pt_170to300_TuneCP5_14TeV_pythia8_wOnlineVtx_wL1/',
        'numEvents' : 50000
      }
    },
    'numJobs' : 10, 
    'crossSection' : 121700.0/200,
    ##'kFactor_pT_hat' : 1.10,
    'kFactor_pT_hat' : 1.0,
    'process' : "QCD"
  },  
  'dy_mass10to50' : {
    'samples' : {
      'offlinePrimaryVertices' : { 
        'inputFilePath' : '',
        'numEvents' : XXX,
      },
      'hltPhase2PixelVertices' : {
        'inputFilePath' : '/hdfs/cms/store/user/rdewanje/DYJetsToLL_M-10to50_TuneCP5_14TeV-madgraphMLM-pythia8/HLTConfig_DYJetsToLL_M-10to50_TuneCP5_14TeV-madgraphMLM-pythia8_wOnlineVtx_wL1/',
        'numEvents' : 96923
      }
    },
    'numJobs' : 10, 
    'crossSection' : 16880.0,
    'process' : "DY"
  },
  'dy_massGt50' : {
    'samples' : {
      'offlinePrimaryVertices' : { 
        'inputFilePath' : '',
        'numEvents' : XXX, 
      },
      'hltPhase2PixelVertices' : {
        'inputFilePath' : '/hdfs/cms/store/user/rdewanje/DYToLL_M-50_TuneCP5_14TeV-pythia8/HLTConfig_DYToLL_M-50_TuneCP5_14TeV-pythia8_wOnlineVtx_wL1/',
        'numEvents' : 10000
      }
    },
    'numJobs' : 10, 
    'crossSection' : 5795.0,
    'process' : "DY"
  },
  'w' : {
    'samples' : {
      'offlinePrimaryVertices' : { 
        'inputFilePath' : '',
        'numEvents' : XXX,
      },
      'hltPhase2PixelVertices' : {
        'inputFilePath' : '/hdfs/cms/store/user/rdewanje/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/HLTConfig_WJetsToLNu_TuneCP5_14TeV_pythia8_wOnlineVtx_wL1/',
        'numEvents' : 63619
      }
    },
    'numJobs' : 10, 
    'crossSection' : 56990.0,
    'process' : "W"
  }
}
