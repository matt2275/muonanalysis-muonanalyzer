'''Author: g. karathanasis. georgios.karathanasis@cern.ch
cfg to run tag and probe ntuple for muon POG. It runs both on AOD and miniAOD
Modified by Andre Frankenthal (a.franken@cern.ch) -- September 2020
usage: cmsRun run_muonAnalyzer_cfg.py option1=value1 option2=value2
'''

from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms

options = VarParsing('python')

options.register('resonance', 'Z',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Set resonance ('Z'/'JPsi')"
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

options.parseArguments()

# defaults

if options.fromCRAB == True:
    options.maxEvents = -1
else:
    options.maxEvents = -1

if options._beenSet['globalTag'] and options.globalTag != '':
    globaltag = options.globalTag
else:
    # globaltag = '103X_dataRun2_v6' if not options.isMC else '103X_upgrade2018_realistic_HI_v11'
    # globaltag = '102X_upgrade2018_realistic_v15' if not options.isMC else '102X_upgrade2018_realistic_v15'    
    globaltag = '106X_dataRun2_v20' if not options.isMC else '102X_upgrade2018_realistic_v15'
    #globaltag = '102X_dataRun2_v11' if not options.isMC else '102X_upgrade2018_realistic_v15'
    #globaltag = '102X_dataRun2_v11' if not options.isMC else '103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2'


# Run local test if no input files provided
if len(options.inputFiles) == 0:
    if options.resonance == 'Z':
        if options.isFullAOD:
            if options.isMC:
                #options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/021657EA-7BFE-B548-AF55-7FF2F732A065.root'),
                #options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/102X_upgrade2018_realistic_v15-v1/80002/FF9AF238-78B6-CF48-BC7C-05025D85A45C.root')
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/021657EA-7BFE-B548-AF55-7FF2F732A065.root')
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/02BEC7FE-546C-0C4A-8D9B-F5BF313236E5.root')
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/06249AB6-2DEE-B640-8AFF-F3B1494864F1.root')
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/06787DBF-220D-954D-96C0-041445C3EEBF.root')
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/0A0EDFFD-419A-784F-9090-C286669EFF95.root')
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/0A2E5082-C970-0A4A-A549-50E3A410BABB.root')
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/0D23C693-BAAD-474D-922A-CF1EDE3B311C.root')
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/0F31BCBA-DAE0-6C4C-9986-C5629CC8A966.root')
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/10A284FE-0894-AC46-BA69-0FE7FC4B20E7.root')
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/10A78652-4801-4248-A188-AF61FDADFF1D.root')
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/11EC1AB9-FF3D-FC49-9383-0A1A1B48EDD5.root')
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/124FB222-00E9-B04E-BB3A-E5A441323BC5.root')
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/1283ECB7-636A-A842-A79E-D0B934FA9024.root')
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/137AC9C7-A3FF-5F40-8E34-25165DC82CFE.root')
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/15D1C459-96DC-0947-B09E-DEAB97C39C18.root')
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/167DC01C-0F17-024C-BB40-AF1A246F714E.root')
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/173651D4-398F-1246-8EFA-DB3557B1437A.root')
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/18B6179B-0DD7-FE4D-A4A9-BDAE5F77C63C.root')
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/19416B30-A13B-7246-9247-8C1A7298C71C.root')
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/1992C353-9D12-2944-B8FA-24B0497C43C1.root')
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/1C2C8436-5196-CE44-A08A-79FD066F0576.root')
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/1D2159DA-CA9C-5747-BEAC-7B1F5157E0C5.root')
                options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/DYToMuMu_M-20_13TeV_pythia8/GEN-SIM-RECO/102X_upgrade2018_realistic_v15-v2/100000/1D73D9D6-3163-3244-BABE-6AA478A30B6F.root')
            else:
                #options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/60004/FB123080-071C-F64D-BAFD-F2F292F7FC64.root')
                #options.inputFiles.append('/store/data/Run2018D/SingleMuon/AOD/12Nov2019_UL2018-v6/230000/1217FEF3-168A-7C46-96EA-F5D8B6C6BC4A.root')
                #options.inputFiles.append('/store/data/Run2018D/SingleMuon/AOD/12Nov2019_UL2018-v8/120000/002F4F96-E55E-8141-B8B8-B6345AB90AA2.root')
                #options.inputFiles.append('/store/data/Run2017E/SingleMuon/AOD/09Aug2019_UL2017-v1/260003/7417C8DE-DCB6-C947-9B0D-99D8B1EA4A49.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/60004/FB123080-071C-F64D-BAFD-F2F292F7FC64.root')
                options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/001498AA-6457-0B42-BEFA-9C72A6A91D2E.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/001D1DB9-405A-D14D-B915-55F071A76BF3.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/0087BDD1-CE79-084C-8A72-CE64EE3B2B3C.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/0089CA7E-648B-CF4C-A59A-4906AD3F26D8.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/008C5121-EC97-BD4C-8461-A0FA43E28DD7.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/008D63AC-B709-704C-A058-1AE819034C56.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/00B09A1D-4AC0-3943-AB47-A81AC75D158F.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/00E49621-DF30-F34E-B84C-332A2AE0B781.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/015821DF-68E2-BB4E-9B8E-D4E129C265AE.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/015A3917-A98F-4D47-9FA0-AA14F8AB6BC8.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/016D3250-5DA8-2145-A862-862ED8F24CFA.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/0195DE92-8133-0249-8A1D-676A8931EA1F.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/0199FEB3-0DFB-A04C-BC05-4C1A7E3DBC2A.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/01BE8D26-F46A-9F40-969E-8D68602D0A5F.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/01F0D35D-7039-BE4F-BF6F-4DD6558F2CD0.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/0237CC50-72C1-D14E-BC9B-18E584E67AA3.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/0250BEDD-6910-6144-8AF5-DB8EA9D79057.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/028A254D-4BCF-A248-BC5D-D491CF441B3C.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/0299B1F1-BF3A-AC4A-93D7-3C731584CE39.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/02A68556-4DD1-3A4B-A184-4CD14A941DEF.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/02EA41F6-F852-9543-8CF1-D1E85D044C34.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/0307E886-1A51-AA41-B550-09E7B6692DC9.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/03110E4E-F57F-5B46-B74E-7F834725B5F0.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/0315D0D3-41D3-1F4D-9F97-1B20A2A9A2B8.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/03BC2DAA-E33F-A849-847D-30C32BEA6231.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/03CF3DFA-EC9F-1F42-ADE1-48420D8945E9.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/03D355DD-113B-784C-B6AD-00F4F6294935.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/03EA2407-B40B-514A-8E5F-F35E27EA9F90.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/04130A67-F541-2D4B-8B78-C56BF0F5516B.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/043D375C-D4FD-FE4A-B573-3EE6E797EF51.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/0501A525-D845-B54F-AB99-103F1276ABD5.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/0554562D-2CBF-CC4D-92AD-57A027E16EC9.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/059351F5-2999-9544-AEBE-DC049F6B2016.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/05CFF9AF-DD02-4549-8130-B516F9A8835E.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/05EB3967-0FFE-7842-9C88-3BAA5C5B8F68.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/06598F25-85A5-044D-A0BD-E3076728AA5A.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/065B2DD7-EDEB-0248-9038-BB15FDA396B8.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/069993E0-3848-C14C-B80B-56443C227C21.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/0733AC52-7ADB-9046-8FFF-4BAAB2CD90DB.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/0769C291-BC1A-E143-AA17-8FE1B3EA032A.root')
                # options.inputFiles.append('/store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00000/07AC2FCD-D220-CF40-89AE-333E8BBED0C1.root')
        else:
            if options.isMC:
                options.inputFiles.append('/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_Zpt-150toInf_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/100000/FAACE9E0-1D0E-204E-9960-078F095EA34C.root')
            else:
                options.inputFiles.append('/store/data/Run2018D/SingleMuon/MINIAOD/PromptReco-v2/000/325/159/00000/BE9BB28C-8AEC-1B4B-A7BD-AD1C9A0D67A8.root')
    elif options.resonance == 'JPsi':
        if options.isFullAOD:
            if options.isMC:
                options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/02762AC2-E426-CA46-A9BE-3B67E3ABC372.root')
                #options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/ggBBbar_4f_TuneCP5_5p02TeV_MG5_aMCatNLO_pythia8/AODSIM/NoPUlowPtPhotonReg_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2_ext1-v2/230000/03FF4C20-EBBF-2944-BD7B-6B9D0CBABCDC.root')
                #options.inputFiles.append('/store/mc/RunIIAutumn18DRPremix/JpsiToMuMu_JpsiPt8_TuneCP5_13TeV-pythia8/AODSIM/102X_upgrade2018_realistic_v15-v1/270001/FFF2FC1D-18CB-7244-9663-4E36963494B7.root')
            else:
                #options.inputFiles.append('/store/data/Run2018A/Charmonium/AOD/17Sep2018-v1/100001/07679496-4DEF-1B44-BA04-768765A80599.root')
                #options.inputFiles.append('/store/hidata/HIRun2015/HIForward/AOD/02May2016-v1/00000/0038049F-0D25-E611-B57F-F01FAFD691F4.root')
                options.inputFiles.append('/store/hidata/HIRun2015/HIForward/AOD/02May2016-v1/00000/046290BB-E417-E611-B342-F01FAFE15CBD.root')
                #options.inputFiles.append('/store/hidata/HIRun2018A/HIForward/AOD/ForLByL-v2/20000/072336D7-34FB-BE4A-816B-49E276923F5E.root')

if options.outputFile=="":
    options.outputFile="output"
    if options.isMC:
        options.outputFile+="_mc"
    else:
        options.outputFile+="_data" 
    if options.isFullAOD:
        options.outputFile+="_full"
    else:
        options.outputFile+="_mini"
    options.outputFile+=".root"


process = cms.Process("MuonAnalysis")

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

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

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
if options.isFullAOD:
    if options.resonance == 'Z':
        process = muonAnalysis_customizeStandAloneFullAOD_Z(process)  # Aug 3 2022 Changed to Standalone Z
    else:
        process = muonAnalysis_customizeFullAOD_JPsi(process)
    if not options.isMC:
        process.muon.jetCorrector = cms.InputTag(
            "ak4PFCHSL1FastL2L3ResidualCorrector")
else:
    if options.resonance == 'Z':
        process = muonAnalysis_customizeMiniAOD_Z(process)
    else:
        process = muonAnalysis_customizeMiniAOD(process)

process.muon.isMC = options.isMC
process.muon.includeJets = options.includeJets
process.muon.era = options.era

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

#### 
#### Adding Old Framework things 
####


## SELECT WHAT DATASET YOU'RE RUNNING ON
TRIGGER="SingleMu"
#TRIGGER="DoubleMu"

## ==== Fast Filters ====
process.goodVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 25 && position.Rho <= 2"),
    filter = cms.bool(True),
)
process.noScraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

process.load("HLTrigger.HLTfilters.triggerResultsFilter_cfi")


if TRIGGER == "SingleMu":
#    process.triggerResultsFilter.triggerConditions = cms.vstring( 'HLT_Mu45_eta2p1_v*', 'HLT_Mu50_v*',
#                                                                  'HLT_IsoMu27_v*',   'HLT_IsoMu24_v*',   'HLT_IsoMu22_v*',   'HLT_IsoMu20_v*',
#                                                                  'HLT_IsoTkMu27_v*', 'HLT_IsoTkMu24_v*', 'HLT_IsoTkMu22_v*', 'HLT_IsoTkMu20_v*'  )
    process.triggerResultsFilter.triggerConditions = cms.vstring('HLT_IsoMu24_v*')
elif TRIGGER == "DoubleMu":
    process.triggerResultsFilter.triggerConditions = cms.vstring( 'HLT_Mu8_v*', 'HLT_Mu17_v*',
                                                                  'HLT_Mu8_TrkIsoVVL_v*', 'HLT_Mu17_TrkIsoVVL_v*',
                                                                  'HLT_Mu17_TkMu8_v*', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*' )
else:
    raise RuntimeError, "TRIGGER must be 'SingleMu' or 'DoubleMu'"

process.triggerResultsFilter.l1tResults = "gtStage2Digis"
process.triggerResultsFilter.throw = False
process.triggerResultsFilter.hltResults = cms.InputTag("TriggerResults","","HLT")

#decomment when you have it
#process.triggerResultsFilterFake = process.triggerResultsFilter.clone(
#    triggerConditions = cms.vstring( 'HLT_Mu40_v*', 'HLT_Mu5_v*', 'HLT_Mu12_v*', 'HLT_Mu24_v*')
#)

process.fastFilter     = cms.Sequence(process.goodVertexFilter + process.noScraping + process.triggerResultsFilter)
#process.fastFilterFake = cms.Sequence(process.goodVertexFilter + process.noScraping + process.triggerResultsFilterFake)

##    __  __                       
##   |  \/  |_   _  ___  _ __  ___ 
##   | |\/| | | | |/ _ \| '_ \/ __|
##   | |  | | |_| | (_) | | | \__ \
##   |_|  |_|\__,_|\___/|_| |_|___/
##                                 
## ==== Merge CaloMuons and Tracks into the collection of reco::Muons  ====
from RecoMuon.MuonIdentification.calomuons_cfi import calomuons;
process.mergedMuons = cms.EDProducer("CaloMuonMerger",
    mergeTracks = cms.bool(True),
    mergeCaloMuons = cms.bool(False), # AOD
    muons     = cms.InputTag("muons"), 
    caloMuons = cms.InputTag("calomuons"),
    tracks    = cms.InputTag("generalTracks"),
    minCaloCompatibility = calomuons.minCaloCompatibility,
    ## Apply some minimal pt cut
    muonsCut     = cms.string("pt > 3 && track.isNonnull"),
    caloMuonsCut = cms.string("pt > 3"),
    tracksCut    = cms.string("pt > 3"),
)

## ==== Trigger matching
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
## with some customization
process.muonMatchHLTL2.maxDeltaR = 0.3 # Zoltan tuning - it was 0.5
process.muonMatchHLTL3.maxDeltaR = 0.1
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *
changeRecoMuonInput(process, "mergedMuons")
useL1Stage2Candidates(process)
appendL1MatchingAlgo(process)
#addHLTL1Passthrough(process)
from MuonAnalysis.TagAndProbe.common_variables_cff import *
process.load("MuonAnalysis.TagAndProbe.common_modules_cff")

process.tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("pt > 27 && "+MuonIDFlags.Tight2012.value()+
                     " && !triggerObjectMatchesByCollection('hltIterL3MuonCandidates').empty()"),
)
process.pseudoTag = cms.EDFilter("MuonSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("pt > 27 && isGlobalMuon && numberOfMatchedStations >= 2 && pfIsolationR04().sumChargedHadronPt/pt < 0.2")
)
if TRIGGER == "DoubleMu":
    process.tagMuons.cut = ("pt > 6 && (isGlobalMuon || isTrackerMuon) && isPFMuon "+
                            " && !triggerObjectMatchesByCollection('hltL3MuonCandidates').empty()"+
                            " && pfIsolationR04().sumChargedHadronPt/pt < 0.2")

process.oneTag  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tagMuons"), minNumber = cms.uint32(1))

process.probeMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("track.isNonnull"),  # no real cut now
)

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    #cut = cms.string('60 < mass < 140 && abs(daughter(0).vz - daughter(1).vz) < 4'),
    cut = cms.string('70 < mass < 115 && abs(daughter(0).charge - daughter(1).charge) < 10'),
    decay = cms.string('tagMuons@+ probeMuons@-')
)
process.onePair = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tpPairs"), minNumber = cms.uint32(1))

from MuonAnalysis.TagAndProbe.muon.tag_probe_muon_extraIso_cff import ExtraIsolationVariables

from MuonAnalysis.TagAndProbe.puppiIso_cfi import load_fullPFpuppiIsolation
process.fullPuppIsolationSequence = load_fullPFpuppiIsolation(process)
from MuonAnalysis.TagAndProbe.puppiIso_cff import PuppiIsolationVariables

process.tpTree = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("None"),
    # probe variables: all useful ones
    variables = cms.PSet(
        AllVariables,
        ExtraIsolationVariables,
        PuppiIsolationVariables,
        isoTrk03Abs = cms.InputTag("probeMuonsIsoValueMaps","probeMuonsIsoFromDepsTk"),
        isoTrk03Rel = cms.InputTag("probeMuonsIsoValueMaps","probeMuonsRelIsoFromDepsTk"),
        dxyBS = cms.InputTag("muonDxyPVdzmin","dxyBS"),
        dxyPVdzmin = cms.InputTag("muonDxyPVdzmin","dxyPVdzmin"),
        dzPV = cms.InputTag("muonDxyPVdzmin","dzPV"),
        JetPtRatio= cms.InputTag("AddLeptonJetRelatedVariables","JetPtRatio"),
        JetPtRel= cms.InputTag("AddLeptonJetRelatedVariables","JetPtRel"),
        JetNDauCharged= cms.InputTag("AddLeptonJetRelatedVariables","JetNDauCharged"),
        JetBTagCSV= cms.InputTag("AddLeptonJetRelatedVariables","JetBTagCSV"),
        miniIsoCharged = cms.InputTag("muonMiniIsoCharged","miniIso"),
        activity_miniIsoCharged = cms.InputTag("muonMiniIsoCharged","activity"),
        miniIsoPUCharged = cms.InputTag("muonMiniIsoPUCharged","miniIso"),
        activity_miniIsoPUCharged = cms.InputTag("muonMiniIsoPUCharged","activity"),
        miniIsoNeutrals = cms.InputTag("muonMiniIsoNeutrals","miniIso"),
        activity_miniIsoNeutrals = cms.InputTag("muonMiniIsoNeutrals","activity"),
        miniIsoPhotons = cms.InputTag("muonMiniIsoPhotons","miniIso"),
        activity_miniIsoPhotons = cms.InputTag("muonMiniIsoPhotons","activity"),
        nSplitTk  = cms.InputTag("splitTrackTagger"),
        mt  = cms.InputTag("probeMetMt","mt"),
        CutBasedIdGlobalHighPt_new = cms.InputTag("muonHighPt","highPtIDNew"),

    ),
    flags = cms.PSet(
       TrackQualityFlags,
       MuonIDFlags,
       HighPtTriggerFlags,
       HighPtTriggerFlagsDebug,
    ),
    tagVariables = cms.PSet(
     #   TriggerVariables, 
     #   MVAIsoVariablesPlainTag, 
     #   pt = cms.string("pt"),
     #   eta = cms.string("eta"),
     #   phi = cms.string("phi"),
     #   combRelIso = cms.string("(isolationR03.emEt + isolationR03.hadEt + isolationR03.sumPt)/pt"),
     #   chargedHadIso04 = cms.string("pfIsolationR04().sumChargedHadronPt"),
     #   neutralHadIso04 = cms.string("pfIsolationR04().sumNeutralHadronEt"),
     #   photonIso04 = cms.string("pfIsolationR04().sumPhotonEt"),
     #   combRelIsoPF04dBeta = IsolationVariables.combRelIsoPF04dBeta,
     #   combRelIsoPF03dBeta = IsolationVariables.combRelIsoPF03dBeta,
     #   dzPV = cms.InputTag("muonDxyPVdzminTags","dzPV"),
        AllVariables,
        # ExtraIsolationVariables,
        nVertices   = cms.InputTag("nverticesModule"),
        # isoTrk03Abs = cms.InputTag("probeMuonsIsoValueMaps","probeMuonsIsoFromDepsTk"),
        # isoTrk03Rel = cms.InputTag("probeMuonsIsoValueMaps","probeMuonsRelIsoFromDepsTk"),
        dxyBS = cms.InputTag("muonDxyPVdzminTags","dxyBS"),
        dxyPVdzmin = cms.InputTag("muonDxyPVdzminTags","dxyPVdzmin"),
        dzPV = cms.InputTag("muonDxyPVdzminTags","dzPV"),
        # nSplitTk  = cms.InputTag("splitTrackTagger"),
        l1rate = cms.InputTag("l1rate"),
        bx     = cms.InputTag("l1rate","bx"),
        #mu17ps = cms.InputTag("l1hltprescale","HLTMu17TotalPrescale"), 
        #mu8ps  = cms.InputTag("l1hltprescale","HLTMu8TotalPrescale"), 
        instLumi = cms.InputTag("addEventInfo", "instLumi"),
        met = cms.InputTag("tagMetMt","met"),
        mt  = cms.InputTag("tagMetMt","mt"),
        CutBasedIdGlobalHighPt_new = cms.InputTag("muonHighPtTags","highPtIDNew"),
    ),
    tagFlags = cms.PSet(
        HighPtTriggerFlags,
        HighPtTriggerFlagsDebug,
        ),
    pairVariables = cms.PSet(
        nJets30 = cms.InputTag("njets30Module"),
        dz      = cms.string("daughter(0).vz - daughter(1).vz"),
        pt      = cms.string("pt"), 
        rapidity = cms.string("rapidity"),
        deltaR   = cms.string("deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)"), 
        probeMultiplicity = cms.InputTag("probeMultiplicity"),
        probeMultiplicity_TMGM = cms.InputTag("probeMultiplicityTMGM"),
        probeMultiplicity_Pt10_M60140 = cms.InputTag("probeMultiplicityPt10M60140"),
        ## New TuneP variables
        newTuneP_probe_pt            = cms.InputTag("newTunePVals", "pt"),
        newTuneP_probe_sigmaPtOverPt = cms.InputTag("newTunePVals", "ptRelError"),
        newTuneP_probe_trackType     = cms.InputTag("newTunePVals", "trackType"),
        newTuneP_mass                = cms.InputTag("newTunePVals", "mass"),
    ),
    pairFlags = cms.PSet(
        BestZ = cms.InputTag("bestPairByZMass"),
    ),
    isMC           = cms.bool(False),
    addRunLumiInfo = cms.bool(True),
)
if TRIGGER == "DoubleMu":
    for K,F in MuonIDFlags.parameters_().iteritems():
        setattr(process.tpTree.tagFlags, K, F)


process.load("MuonAnalysis.TagAndProbe.muon.tag_probe_muon_extraIso_cfi")
process.load("PhysicsTools.PatAlgos.recoLayer0.pfParticleSelectionForIso_cff")

process.miniIsoSeq = cms.Sequence(
    process.pfParticleSelectionForIsoSequence +
    process.muonMiniIsoCharged + 
    process.muonMiniIsoPUCharged + 
    process.muonMiniIsoNeutrals + 
    process.muonMiniIsoPhotons 
)

# process.load("JetMETCorrections.Configuration.JetCorrectionProducersAllAlgos_cff")
# process.ak4PFCHSJetsL1L2L3 = process.ak4PFCHSJetsL1.clone( correctors = ['ak4PFCHSL1FastL2L3'] )

process.extraProbeVariablesSeq = cms.Sequence(
    process.probeMuonsIsoSequence +
    process.computeCorrectedIso + 
    process.splitTrackTagger +
    process.muonDxyPVdzmin + 
    process.muonHighPt + 
    process.probeMetMt + process.tagMetMt +
    process.miniIsoSeq +
    # process.ak4PFCHSJetsL1L2L3 +
    process.ak4PFCHSL1FastL2L3CorrectorChain * process.AddLeptonJetRelatedVariables +
    process.fullPuppIsolationSequence 
)

process.tnpSimpleSequence = cms.Sequence(
    process.tagMuons +
    process.oneTag     +
    process.probeMuons +
    process.tpPairs    +
    process.onePair    +
    process.nverticesModule +
    process.njets30Module +
    process.extraProbeVariablesSeq +
    process.probeMultiplicities + 
    process.addEventInfo +
    process.l1rate +
    #process.l1hltprescale + 
    process.bestPairByZMass + 
    process.newTunePVals +
    process.muonDxyPVdzminTags +
    process.muonHighPtTags + 
    process.tpTree
)

process.tagAndProbe = cms.Path( 
    process.fastFilter +
    process.mergedMuons                 *
    process.patMuonsWithTriggerSequence +
    process.tnpSimpleSequence
)

##    _____               _    _             
##   |_   _| __ __ _  ___| | _(_)_ __   __ _ 
##     | || '__/ _` |/ __| |/ / | '_ \ / _` |
##     | || | | (_| | (__|   <| | | | | (_| |
##     |_||_|  \__,_|\___|_|\_\_|_| |_|\__, |
##                                     |___/ 

## Then make another collection for standalone muons, using standalone track to define the 4-momentum
process.muonsSta = cms.EDProducer("RedefineMuonP4FromTrack",
    src   = cms.InputTag("muons"),
    track = cms.string("outer"),
)
## Match to trigger, to measure the efficiency of HLT tracking
from PhysicsTools.PatAlgos.tools.helpers import *
process.patMuonsWithTriggerSequenceSta = cloneProcessingSnippet(process, process.patMuonsWithTriggerSequence, "Sta")
process.patMuonsWithTriggerSequenceSta.replace(process.patTriggerFullSta, process.patTriggerFull)
process.patTriggerSta.src = 'patTriggerFull'
process.muonMatchHLTL2Sta.maxDeltaR = 0.5
process.muonMatchHLTL3Sta.maxDeltaR = 0.5
massSearchReplaceAnyInputTag(process.patMuonsWithTriggerSequenceSta, "mergedMuons", "muonsSta")
#massSearchReplaceAnyInputTag(process.patMuonsWithTriggerSequenceSta, "muons", "muonsSta")
## Define probes and T&P pairs
process.probeMuonsSta = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTriggerSta"),
    cut = cms.string("pt > 3 && outerTrack.isNonnull"), # no real cut now
)
process.pseudoProbeSta = cms.EDFilter("MuonSelector",
    src = cms.InputTag("muonsSta"),
    cut = cms.string("pt > 3 && outerTrack.isNonnull"),
)


process.tpPairsSta = process.tpPairs.clone(decay = "tagMuons@+ probeMuonsSta@-", cut = '70 < mass < 115')

process.onePairSta = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tpPairsSta"), minNumber = cms.uint32(1))

process.pseudoPairsSta = process.tpPairsSta.clone(decay = "pseudoTag@+ pseudoProbeSta@-")
process.onePseudoPairSta = process.onePairSta.clone(src = 'pseudoPairsSta')
process.fastPseudoTnPSta = cms.Sequence(process.pseudoTag + process.muonsSta + process.pseudoProbeSta + process.pseudoPairsSta + process.onePseudoPairSta)

process.staToTkMatch.maxDeltaR     = 0.3
process.staToTkMatch.maxDeltaPtRel = 2.
process.staToTkMatchNoZ.maxDeltaR     = 0.3
process.staToTkMatchNoZ.maxDeltaPtRel = 2.

process.load("MuonAnalysis.TagAndProbe.tracking_reco_info_cff")

process.tpTreeSta = process.tpTree.clone(
    tagProbePairs = "tpPairsSta",
    arbitration   = "OneProbe",
    variables = cms.PSet(
        KinematicVariables, 
        StaOnlyVariables,
        ## track matching variables
        tk_deltaR     = cms.InputTag("staToTkMatch","deltaR"),
        tk_deltaEta   = cms.InputTag("staToTkMatch","deltaEta"),
        tk_deltaR_NoZ   = cms.InputTag("staToTkMatchNoZ","deltaR"),
        tk_deltaEta_NoZ = cms.InputTag("staToTkMatchNoZ","deltaEta"),
    ),
    flags = cms.PSet(
        outerValidHits = cms.string("outerTrack.numberOfValidHits > 0"),
        TM  = cms.string("isTrackerMuon"),
        Glb = cms.string("isGlobalMuon"),
        Tk  = cms.string("track.isNonnull"),
        StaTkSameCharge = cms.string("outerTrack.isNonnull && innerTrack.isNonnull && (outerTrack.charge == innerTrack.charge)"),
    ),
    tagVariables = cms.PSet(
        pt = cms.string("pt"),
        eta = cms.string("eta"),
        phi = cms.string("phi"),
        nVertices = cms.InputTag("nverticesModule"),
        combRelIso = cms.string("(isolationR03.emEt + isolationR03.hadEt + isolationR03.sumPt)/pt"),
        chargedHadIso04 = cms.string("pfIsolationR04().sumChargedHadronPt"),
        neutralHadIso04 = cms.string("pfIsolationR04().sumNeutralHadronEt"),
        photonIso04 = cms.string("pfIsolationR04().sumPhotonEt"),
        combRelIsoPF04dBeta = IsolationVariables.combRelIsoPF04dBeta,
        l1rate = cms.InputTag("l1rate"),
        bx     = cms.InputTag("l1rate","bx"),
        instLumi = cms.InputTag("addEventInfo", "instLumi"),
    ),
    pairVariables = cms.PSet(
        nJets30 = cms.InputTag("njets30ModuleSta"),
        dz      = cms.string("daughter(0).vz - daughter(1).vz"),
        pt      = cms.string("pt"), 
        rapidity = cms.string("rapidity"),
        deltaR   = cms.string("deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)"), 
    ),
    pairFlags = cms.PSet(),
)
process.njets30ModuleSta = process.njets30Module.clone(pairs = "tpPairsSta")

process.tnpSimpleSequenceSta = cms.Sequence(
    process.tagMuons +
    process.oneTag     +
    process.probeMuonsSta   +
    process.tpPairsSta      +
    process.onePairSta      +
    process.nverticesModule +
    process.staToTkMatchSequenceZ +
    process.njets30ModuleSta +
    process.addEventInfo +
    process.l1rate +
    process.tpTreeSta
)

## Add extra RECO-level info
if False:
    process.tnpSimpleSequenceSta.replace(process.tpTreeSta, process.tkClusterInfo+process.tpTreeSta)
    process.tpTreeSta.tagVariables.nClustersStrip = cms.InputTag("tkClusterInfo","siStripClusterCount")
    process.tpTreeSta.tagVariables.nClustersPixel = cms.InputTag("tkClusterInfo","siPixelClusterCount")
    process.tnpSimpleSequenceSta.replace(process.tpTreeSta, process.tkLogErrors+process.tpTreeSta)
    process.tpTreeSta.tagVariables.nLogErrFirst = cms.InputTag("tkLogErrors","firstStep")
    process.tpTreeSta.tagVariables.nLogErrPix   = cms.InputTag("tkLogErrors","pixelSteps")
    process.tpTreeSta.tagVariables.nLogErrAny   = cms.InputTag("tkLogErrors","anyStep")

if False: 
    process.tracksNoMuonSeeded = cms.EDFilter("TrackSelector",
              src = cms.InputTag("generalTracks"),
              cut = cms.string(" || ".join("isAlgoInMask('%s')" % a for a in [
              'initialStep', 'lowPtTripletStep', 'pixelPairStep', 'detachedTripletStep',
              'mixedTripletStep', 'pixelLessStep', 'tobTecStep', 'jetCoreRegionalStep',
              'lowPtQuadStep', 'highPtTripletStep', 'detachedQuadStep' ] ) )
    )
    process.pCutTracks0 = process.pCutTracks.clone(src = 'tracksNoMuonSeeded')
    process.tkTracks0 = process.tkTracks.clone(src = 'pCutTracks0')
    process.tkTracksNoZ0 = process.tkTracksNoZ.clone(src = 'tkTracks0')
    process.preTkMatchSequenceZ.replace(
            process.tkTracksNoZ, process.tkTracksNoZ +
            process.tracksNoMuonSeeded + process.pCutTracks0 + process.tkTracks0 + process.tkTracksNoZ0)
    process.staToTkMatch0 = process.staToTkMatch.clone(matched = 'tkTracks0')
    process.staToTkMatchNoZ0 = process.staToTkMatchNoZ.clone(matched = 'tkTracksNoZ0')
    process.staToTkMatchSequenceZ.replace( process.staToTkMatch, process.staToTkMatch + process.staToTkMatch0 )
    process.staToTkMatchSequenceZ.replace( process.staToTkMatchNoZ, process.staToTkMatchNoZ + process.staToTkMatchNoZ0 )
    process.tpTreeSta.variables.tk0_deltaR     = cms.InputTag("staToTkMatch0","deltaR")
    process.tpTreeSta.variables.tk0_deltaEta   = cms.InputTag("staToTkMatch0","deltaEta")
    process.tpTreeSta.variables.tk0_deltaR_NoZ   = cms.InputTag("staToTkMatchNoZ0","deltaR")
    process.tpTreeSta.variables.tk0_deltaEta_NoZ = cms.InputTag("staToTkMatchNoZ0","deltaEta")

process.tagAndProbeSta = cms.Path( 
    process.fastFilter +
    process.fastPseudoTnPSta +
    process.mergedMuons * process.patMuonsWithTriggerSequence +
    process.muonsSta                       +
    process.patMuonsWithTriggerSequenceSta +
    process.tnpSimpleSequenceSta
)

##### 
##### End of Tracking
##### 

process.OldFilter     = cms.Sequence(
    process.fastFilter +
    process.fastPseudoTnPSta +
    process.mergedMuons * process.patMuonsWithTriggerSequence +
#    process.patMuonsWithTriggerSequence  +
    process.muonsSta                       +
    process.patMuonsWithTriggerSequenceSta +
    process.tagMuons +
    process.oneTag     +
    process.probeMuonsSta   +
    process.tpPairsSta      +
    process.onePairSta   
    )

####
#### Edd of Old Framework Things
####

###
###


###
###
from MuonAnalysis.MuonAnalyzer.hltInfo_SATest_cff import getHLTInfo, selectTriggers # Aug 3 2022 changed to SATest
hltInfo = getHLTInfo(options.resonance, options.era)
excludeDSA = (not options.isFullAOD)
process.muon.triggerPaths = cms.vstring(selectTriggers(hltInfo['triggerPaths'], True, False, excludeDSA))
process.muon.tagFilters = cms.vstring(selectTriggers(hltInfo['tagFilters'], not options.isFullAOD, True, excludeDSA))
process.muon.probeFilters = cms.vstring(selectTriggers(hltInfo['probeFilters'], not options.isFullAOD, True, excludeDSA))

# Standard selectors
from MuonAnalysis.MuonAnalyzer.selectorInfo_cff import getSelectorNamesAndBits
selectorNames, selectorBits = getSelectorNamesAndBits(options.era, options.isFullAOD)
process.muon.probeSelectorNames = cms.vstring(selectorNames)
process.muon.probeSelectorBits = cms.vuint32(selectorBits)

if options.includeJets:
    if not options.isMC:
        process.analysis_step = cms.Path(
            process.muonL1Info +
            process.muonL1InfoByQ +
            process.ak4PFCHSL1FastL2L3ResidualCorrectorChain +
            process.muSequence
        )
    else:
        process.analysis_step = cms.Path(
            process.muonL1Info +
            process.muonL1InfoByQ +
            process.ak4PFCHSL1FastL2L3CorrectorChain +
            process.muSequence
        )
else:
    process.analysis_step = cms.Path(
        process.OldFilter +
        process.muonL1Info +
        process.muonL1InfoByQ +
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

process.schedule = cms.Schedule(process.tagAndProbeSta,process.analysis_step, process.endjob_step)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
