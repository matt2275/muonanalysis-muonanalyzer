'''Author: g. karathanasis. georgios.karathanasis@cern.ch
cfg to run tag and probe ntuple for muon POG. It runs both on AOD and miniAOD
Modified by Andre Frankenthal (a.franken@cern.ch) -- September 2020
usage: cmsRun run_muonAnalyzer_cfg.py option1=value1 option2=value2

Example: cmsRun test/run_analysis_muonAnalyzer_cfg.py isMC=True isFullAOD=True mcType=EE isHIUPC=True

'''

from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms

options = VarParsing('python')

# options.register('resonance', 'Z',
    # VarParsing.multiplicity.singleton,
    # VarParsing.varType.string,
    # "Set resonance ('Z'/'JPsi')"
# )

options.register('mcType', 'TauTau',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Set resonance ('MuMu'/'TauTau'/'EE'/'Other')"
)


options.register('isFullAOD', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Set to False for MiniAOD datatier"
)

options.register('isMC', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Set to True for MC"
)

options.register('globalTag', '',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Set global tag"
)

options.register('reportEvery', 1000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report frequency"
)


options.register('numThreads', 1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Number of CMSSW threads" 
)

# this parameter is added for Jet Branches (ID varies for different era)
options.register('era', 'Run2018',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "era"
)

options.register('includeJets', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Set to False to exclude jets information in output ntuples"
)

options.register('fromCRAB', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Is config run from CRAB"
)

# options.register('isStandAlone', True,
    # VarParsing.multiplicity.singleton,
    # VarParsing.varType.bool,
    # "run StandAlone Muon Analyzers"
# )

options.register('isHIUPC', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "run HIUPC Muon Analyzers"
)

# options.register('FileNum', 0,
    # VarParsing.multiplicity.singleton,
    # VarParsing.varType.int,
    # "Number of Mu File" 
# )
options.parseArguments()

# defaults


if options._beenSet['globalTag'] and options.globalTag != '':
    globaltag = options.globalTag
else:
    # globaltag = '102X_dataRun2_v11' if not options.isMC else '102X_upgrade2018_realistic_v15'
    globaltag = '103X_dataRun2_Prompt_LowPtPhotonReg_v1' if not options.isMC else '102X_upgrade2018_realistic_v15'
    # globaltag = '103X_upgrade2018_realistic_HI_v1' if not options.isMC else '102X_upgrade2018_realistic_v15'

# Run local test if no input files provided
if len(options.inputFiles) == 0:
    if options.isMC:
       if options.mcType == "TauTau":
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/ggTauTau_TuneCP5_5p02TeV_amcatnlo_pythia8/AODSIM/NoPUlowPtPhotonReg_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v2/280000/0578D671-11C5-8843-88C5-67E784ACB9E0.root')
          options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/ggTauTau_TuneCP5_5p02TeV_SuperChic_pythia8/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v2/10000/05D69AFE-D2D0-6647-803A-E2F4FA2D87F1.root')
          # options.inputFiles.append('/store/group/comm_luminosity/gkrintir/GammaGammatoTauTau_5p02TeV_gammaUPCChFF-pLHE-v3/GammaGammatoTauTau_5p02TeV_gammaUPCChFF_AODSIM-v3/220906_115337/0000/HIN-HINPbPbAutumn18DR-00184_1.root')
          
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/ggTauTau_TuneCP5_5p02TeV_SuperChic_pythia8/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v2/10000/05D69AFE-D2D0-6647-803A-E2F4FA2D87F1.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/ggTauTau_TuneCP5_5p02TeV_SuperChic_pythia8/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v2/10000/32AD5A11-4F91-E24A-BE84-0C6B07455E05.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/ggTauTau_TuneCP5_5p02TeV_SuperChic_pythia8/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v2/10000/42076A67-5F1A-CA43-BEB7-038D1EC71E69.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/ggTauTau_TuneCP5_5p02TeV_SuperChic_pythia8/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v2/10000/A40236B3-F858-BF45-8FA1-51878A8955FB.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/ggTauTau_TuneCP5_5p02TeV_SuperChic_pythia8/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v2/10000/B53B304D-0905-3145-8E18-B5A9A374A2FB.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/ggTauTau_TuneCP5_5p02TeV_SuperChic_pythia8/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v2/10000/E17D3970-920C-EE45-8718-27D43F97882A.root')
       if options.mcType == "MuMu":
           # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/02762AC2-E426-CA46-A9BE-3B67E3ABC372.root')
           # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/0A46C879-0011-BD44-BED3-E05053DAFE28.root')
           
           
           # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/0A46C879-0011-BD44-BED3-E05053DAFE28.root')
           options.inputFiles.append('/store/user/ajofrehe/gtauSim/mumuFSR/output/gammaUPCmumuFSR/gammaUPCmumuFSR-lowPtPhoton/230223_095036/0000/ggmumuFSR-EDFF-gammaUPC2018_1.root')
           # options.inputFiles.append('/store/user/ajofrehe/gtauSim/mumuFSR/output/gammaUPCmumuFSR/gammaUPCmumuFSR-lowPtPhoton/230223_095036/0000/ggmumuFSR-EDFF-gammaUPC2018_10.root')
           # options.inputFiles.append('/store/user/ajofrehe/gtauSim/mumuFSR/output/gammaUPCmumuFSR/gammaUPCmumuFSR-lowPtPhoton/230223_095036/0000/ggmumuFSR-EDFF-gammaUPC2018_100.root')
           # options.inputFiles.append('/store/user/ajofrehe/gtauSim/mumuFSR/output/gammaUPCmumuFSR/gammaUPCmumuFSR-lowPtPhoton/230223_095036/0000/ggmumuFSR-EDFF-gammaUPC2018_101.root')
           # options.inputFiles.append('/store/user/ajofrehe/gtauSim/mumuFSR/output/gammaUPCmumuFSR/gammaUPCmumuFSR-lowPtPhoton/230223_095036/0000/ggmumuFSR-EDFF-gammaUPC2018_102.root')
           # options.inputFiles.append('/store/user/ajofrehe/gtauSim/mumuFSR/output/gammaUPCmumuFSR/gammaUPCmumuFSR-lowPtPhoton/230223_095036/0000/ggmumuFSR-EDFF-gammaUPC2018_103.root')      
           
           # options.inputFiles.append('file:/eos/user/a/ajofrehe/4Pranati/AOD/DR20{:02d}.root'.format(options.FileNum))           
 

           # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/02762AC2-E426-CA46-A9BE-3B67E3ABC372.root')
           # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/02D64F67-B2B4-494E-B9B1-5D7B0DD22244.root')
           # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/04B4C8C5-DE0D-6146-AC1C-F27DAD4EE6C0.root')
           # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/04CD57C1-5BBA-A94A-848D-8B036B7786BE.root')
           # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/0659B805-B026-6444-B784-989A8663E582.root') 
        
# Files for MuMu SL public gen Efficiency

          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/02762AC2-E426-CA46-A9BE-3B67E3ABC372.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/02D64F67-B2B4-494E-B9B1-5D7B0DD22244.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/04B4C8C5-DE0D-6146-AC1C-F27DAD4EE6C0.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/04CD57C1-5BBA-A94A-848D-8B036B7786BE.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/0659B805-B026-6444-B784-989A8663E582.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/06A5A5A1-FCAD-FE4B-A154-CAAD5CDDFF2D.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/06CB403B-96BC-864E-A353-3FDF7253749D.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/0771E778-2FD1-5E48-8102-21CD5E012D3E.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/078C31D4-71A6-D445-BD23-B4E5CF9A5EA5.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/07C74492-647B-404F-B83E-AC109F8A18D5.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/08279F99-F25B-ED4D-823A-A469C234124D.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/094D9AD3-BDDF-B642-9B0C-57731D218724.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/097D219B-3FF7-1B44-8C51-B1CF1732699E.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/09C3768B-D62B-8F40-B46F-6BAB48827981.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/0A46C879-0011-BD44-BED3-E05053DAFE28.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/0AEB3D71-E017-D840-BCA0-6619AB2CA778.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/0BE2782D-B1DB-4D49-BD0E-594B2220634F.root')


# files for MuMu SL private ( No Pt Cut ) gen efficiency


          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_1.root')
          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_10.root')
          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_100.root')
          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_101.root')
          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_102.root')
          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_103.root')
          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_104.root')
          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_105.root')
          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_106.root')
          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_107.root')
          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_108.root')
          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_109.root')
          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_11.root')
          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_110.root')
          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_111.root')
          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_112.root')
          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_113.root')
          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_114.root')
          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_115.root')
          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_116.root')
          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_117.root')
          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_118.root')
          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_119.root')
          # options.inputFiles.append('/store/user/shuaiy/RiceHIN/STARlight/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/GammaGamma2MuMu_GENSIM_woPtCut_PbPb5TeV_v1/STARlight_GammaGamma2MuMu_Reco_woPtCut_PbPb5TeV_v1/191111_221615/0000/step3_RECO_12.root')
        
       if options.mcType == "EE":
          options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_STARlight/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v5/270000/00F70C94-633E-EB45-807C-2FA338E5401F.root') 
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_SuperChic/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v4/110000/0059A52B-5A4C-AF40-8CFF-9C7442A4343E.root')

# files for EE SL gen Efficiency Study

          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_STARlight/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v5/270000/00F70C94-633E-EB45-807C-2FA338E5401F.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_STARlight/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v5/270000/01142839-02B6-0B45-A57B-2FAAAC715A67.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_STARlight/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v5/270000/011CCA66-CA75-8C45-AF06-9AB715FC1453.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_STARlight/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v5/270000/01DD01C9-9A2A-D24E-B3D3-7F7E10976201.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_STARlight/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v5/270000/02CC5C5B-2A60-D34D-AEA9-A95610F06868.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_STARlight/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v5/270000/02E37F85-4C9E-B045-AE4D-85B7483DB31D.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_STARlight/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v5/270000/03309A71-3815-F848-9F3A-9A4DB908EDB6.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_STARlight/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v5/270000/0446E44C-48D1-C746-AD18-3B5C42A133B4.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_STARlight/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v5/270000/06AC07E7-23CD-384D-81E6-80CBAD44695C.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_STARlight/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v5/270000/07684F47-B7E0-CD43-9E8E-39F5948396B2.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_STARlight/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v5/270000/084020B0-7363-EF44-9D39-09CE038E7320.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_STARlight/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v5/270000/0890569E-4C91-234D-9524-2D45144320F9.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_STARlight/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v5/270000/090A46FE-4015-DB48-B61B-33BED7BE0C7E.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_STARlight/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v5/270000/0BA0E6EF-1949-294F-8DFA-183259867D5B.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_STARlight/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v5/270000/0C080B0F-B9D0-494F-91B1-9505C1548BB0.root')
          
          
# files for EE SC gen Efficiency study

          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_SuperChic/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v4/110000/0059A52B-5A4C-AF40-8CFF-9C7442A4343E.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_SuperChic/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v4/110000/00A94507-8E92-7A40-981D-95C2EC5F7882.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_SuperChic/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v4/110000/0160E78D-C745-A441-870E-20AA5E996410.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_SuperChic/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v4/110000/01DA896A-1209-994E-AF9A-6774C5DF2922.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_SuperChic/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v4/110000/0284C703-6E73-384F-8EB4-E5C93BF3CACA.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_SuperChic/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v4/110000/032F56BF-BFB9-FF4A-9C83-CD0AA79ABA71.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_SuperChic/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v4/110000/045419DC-6C7A-244C-80B1-E312BE178393.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_SuperChic/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v4/110000/0513C325-AFAB-074E-B9D9-B74A401CA666.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_SuperChic/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v4/110000/07AB4411-753C-934C-8E4C-2255F9068812.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_SuperChic/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v4/110000/0806613B-650F-EB4A-98AE-D40B9655A77C.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_SuperChic/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v4/110000/0B682A08-0749-2946-8F4B-F9B78C48204F.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_SuperChic/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v4/110000/0D276E25-CCAF-6149-9967-ECE5089465C1.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_SuperChic/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v4/110000/0E20DA10-DBC2-FF41-B851-A2E0B1472DA2.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_SuperChic/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v4/110000/0ECC7A22-BCBD-F94D-BA1A-891A26F0A911.root')
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_SuperChic/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v4/110000/1038B3B2-4EEB-D04A-AF18-1561E9E9E26D.root')
          
       if options.mcType == "JPsi":
          # options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/JPsi_pThat-2_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/AODSIM/mva98_103X_upgrade2018_realistic_HI_v11-v1/120000/00004340-2236-EB40-8E63-BF5F66A147DC.root') 
          options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/JPsi_pThat-2_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/AODSIM/mva98_103X_upgrade2018_realistic_HI_v11-v1/120000/002F1011-2C55-BA46-A9BE-D684848356FB.root') 
          options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/JPsi_pThat-2_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/AODSIM/mva98_103X_upgrade2018_realistic_HI_v11-v1/120000/0039314F-C1B2-344B-A4B3-D597D225DF24.root') 
          
       if options.mcType == "Other":
          options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/ggBBbar_4f_TuneCP5_5p02TeV_MG5_aMCatNLO_pythia8/AODSIM/NoPUlowPtPhotonReg_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2_ext1-v2/230000/03E94DC0-FE38-4E49-9445-CA999E866E6C.root')              
    else:
        # options.inputFiles.append('/store/hidata/HIRun2018A/HIForward/AOD/04Apr2019-v1/00000/36696797-7247-4F47-843D-A2D3630801C5.root')   
        options.inputFiles.append('/store/hidata/HIRun2018A/HIForward/AOD/ForLByL-v2/240000/AF5F01C7-6411-A442-96B0-D9672D8830E8.root')
        
        # options.inputFiles.append('/store/hidata/HIRun2018A/HIForward/AOD/ForLByL-v2/20000/03158A59-174B-A948-8CBA-5684193F9B80.root')
        # options.inputFiles.append('/store/hidata/HIRun2018A/HIForward/AOD/ForLByL-v2/20000/03C086F2-733A-A44B-BC2F-BE755070E42B.root')
        # options.inputFiles.append('/store/hidata/HIRun2018A/HIForward/AOD/ForLByL-v2/20000/040644DA-861F-2B41-A894-212103D23CB8.root')
        # options.inputFiles.append('/store/hidata/HIRun2018A/HIForward/AOD/ForLByL-v2/20000/06C83806-A438-5849-87DF-643363E6BEE1.root')
        # options.inputFiles.append('/store/hidata/HIRun2018A/HIForward/AOD/ForLByL-v2/20000/0712DDCB-1349-3B49-96D6-12F6713F9F96.root')


options.outputFile=""

if options.outputFile == ".root":
    if options.isMC:
        options.outputFile="_" + options.mcType + "_mc" +options.outputFile
        # options.outputFile="_" + options.mcType + "_mc_20{:02d}".format(options.FileNum) +options.outputFile
    else:
        options.outputFile="_data"+ options.outputFile
        
    if options.isFullAOD:
        options.outputFile="output_full"+ options.outputFile
        # options.outputFile="/eos/user/m/mnickel/CrabOut_MultiTag/DiMuon_Arash/output_full"+ options.outputFile
    else:
        options.outputFile+="output_mini"+ options.outputFile

print ('Output Root file:')
print (options.outputFile)


process = cms.Process("MuonAnalysis")



# process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
    # ignoreTotal = cms.untracked.int32(1),
    # showMallocInfo = cms.untracked.bool(True),
    # monitorPssAndPrivate = cms.untracked.bool(True)
# )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,globaltag, '')

# only used for J/Psi centrality 

# if(options.isMC):
    # print('\n\033[31m~*~ USING CENTRALITY TABLE FOR Hydjet Drum5F ~*~\033[0m\n')
    # process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
    # process.GlobalTag.toGet.extend([
        # cms.PSet(record = cms.string("HeavyIonRcd"),
            # # tag = cms.string("CentralityTable_HFtowers200_HydjetDrum5Ev8_v1030pre5x02_mc"),
            # tag = cms.string("CentralityTable_HFtowers200_HydjetDrum5F_v1032x01_mc"),
            # connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
            # label = cms.untracked.string("HFtowers")
            # ),
        # ])
# else:
    # print('\n\033[31m~*~ USING CENTRALITY TABLE FOR PbPb 2018 DATA ~*~\033[0m\n')
    # process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
    # process.GlobalTag.toGet.extend([
        # cms.PSet(record = cms.string("HeavyIonRcd"),
            # tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run2v1031x02_offline"),
            # connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
            # label = cms.untracked.string("HFtowers")
            # ),
        # ])

# process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
# process.centralityBin.Centrality = cms.InputTag("hiCentrality")
# process.centralityBin.centralityVariable = cms.string("HFtowers")

# process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))
# process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000));

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(options.inputFiles),
        secondaryFileNames=cms.untracked.vstring(),
        inputCommands=cms.untracked.vstring(
            'keep *',
            'drop *_ctppsPixelClusters_*_*'
        )
)

if options.includeJets:
    # for b-tagging
    process.load("RecoBTag.ImpactParameter.impactParameter_cff")
    process.load("RecoBTag.SecondaryVertex.secondaryVertex_cff")
    process.load("RecoBTag.SoftLepton.softLepton_cff")
    process.load("RecoBTag.Combined.combinedMVA_cff")
    process.load("RecoBTag.CTagging.cTagging_cff")
    process.load("RecoBTag.Combined.deepFlavour_cff")
    process.load("JetMETCorrections.Configuration.JetCorrectors_cff")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    numberOfThreads = cms.untracked.uint32(options.numThreads)
)

from MuonAnalysis.MuonAnalyzer.tools.ntuple_tools import *
process = muonAnalysis_customizeHIUPCStandAloneFullAOD_Analysis(process)  
if not options.isMC:
   process.muon.jetCorrector = cms.InputTag("ak4PFCHSL1FastL2L3ResidualCorrector")
   

process.muon.isMC = options.isMC
process.muon.MCType = options.mcType
process.muon.includeJets = options.includeJets
process.muon.era = options.era

if options.isHIUPC:
   process.load('RecoHI.ZDCRecHit.QWZDC2018Producer_cfi')
   process.load('RecoHI.ZDCRecHit.QWZDC2018RecHit_cfi')

###############################################################################
#Recover peripheral primary vertices
#https://twiki.cern.ch/twiki/bin/view/CMS/HITracking2018PbPb#Peripheral%20Vertex%20Recovery
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesRecovery_cfi")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")


# Trigger matching
muonSrc = "muons" if options.isFullAOD else "slimmedMuons"
from MuonAnalysis.MuonAssociators.muonL1Match_cfi import muonL1Match as _muonL1Match
process.muonL1Info = _muonL1Match.clone(
    src = cms.InputTag(muonSrc),
    useMB2InOverlap = cms.bool(True),
    useStage2L1 = cms.bool(True),
    preselection = cms.string(""),
    matched = cms.InputTag("gmtStage2Digis:Muon:")
)
process.muonL1InfoByQ = process.muonL1Info.clone(
    sortBy = cms.string("quality"),
    sortByQuality  = cms.bool(True),
    sortByDeltaPhi = cms.bool(False),
    sortByDeltaEta = cms.bool(False),
    sortByPt       = cms.bool(False)
)

from MuonAnalysis.MuonAnalyzer.hltInfo_HI_Analysis_cff import getHLTInfo, selectTriggers
hltInfo = getHLTInfo(options.era)  # add HIUPC section and remove resonance
excludeDSA = (not options.isFullAOD)
process.muon.triggerPaths = cms.vstring(selectTriggers(hltInfo['triggerPaths'], True, False, excludeDSA))
process.muon.tagFilters = cms.vstring(selectTriggers(hltInfo['tagFilters'], not options.isFullAOD, True, excludeDSA))
process.muon.probeFilters = cms.vstring(selectTriggers(hltInfo['probeFilters'], not options.isFullAOD, True, excludeDSA))

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleListDrawer",
  maxEventsToPrint = cms.untracked.int32(100),
  printVertex = cms.untracked.bool(True),
  printOnlyHardInteraction = cms.untracked.bool(False), # Print only status=3 particles. This will not work for Pythia8, which does not have any such particles.
  src = cms.InputTag("genParticles")
)

# Standard selectors
from MuonAnalysis.MuonAnalyzer.selectorInfo_cff import getSelectorNamesAndBits
selectorNames, selectorBits = getSelectorNamesAndBits(options.era, options.isFullAOD)
process.muon.probeSelectorNames = cms.vstring(selectorNames)
process.muon.probeSelectorBits = cms.vuint32(selectorBits)

# if options.includeJets:
    # if not options.isMC:
        # process.analysis_step = cms.Path(
            # process.offlinePrimaryVerticesRecovery +
            # process.muonL1Info +
            # process.muonL1InfoByQ +
            # process.ak4PFCHSL1FastL2L3ResidualCorrectorChain +
            # process.muSequence
        # )
    # else:
        # process.analysis_step = cms.Path(
            # process.muonL1Info +
            # process.muonL1InfoByQ +
            # process.ak4PFCHSL1FastL2L3CorrectorChain +
            # process.muSequence
        # )
# else:
if not options.isMC and options.isHIUPC:
    if options.includeJets:
        process.analysis_step = cms.Path(
            process.offlinePrimaryVerticesRecovery +
            process.muonL1Info +
            process.muonL1InfoByQ +
            process.zdcdigi +
            process.QWzdcreco +
            process.ak4PFCHSL1FastL2L3ResidualCorrectorChain +
            process.muSequence
        )
    else:
        process.analysis_step = cms.Path(
            process.offlinePrimaryVerticesRecovery +
            process.muonL1Info +
            process.muonL1InfoByQ +
            process.zdcdigi +
            process.QWzdcreco +
            # process.hiCentrality +
            # process.centralityBin + 
            process.muSequence
        )
else:
    if options.includeJets:
        process.analysis_step = cms.Path(
            process.offlinePrimaryVerticesRecovery +
            process.muonL1Info +
            process.muonL1InfoByQ +
            process.ak4PFCHSL1FastL2L3CorrectorChain +
            process.muSequence
        )
    else:
        process.analysis_step = cms.Path(
            # process.printTree+
            process.offlinePrimaryVerticesRecovery +
            process.muonL1Info +
            process.muonL1InfoByQ +
            # process.hiCentrality +
            # process.centralityBin + 
            process.muSequence
        )


process.TFileService = cms.Service("TFileService",
        fileName = cms.string(options.outputFile)
)
process.endjob_step = cms.EndPath(process.endOfProcess)

# process.fevt = cms.OutputModule("PoolOutputModule",
#     outputCommands = cms.untracked.vstring(),
#     fileName = cms.untracked.string("edm_output.root")
# )

process.schedule = cms.Schedule(process.analysis_step, process.endjob_step)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
