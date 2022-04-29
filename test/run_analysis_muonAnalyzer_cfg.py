'''Author: g. karathanasis. georgios.karathanasis@cern.ch
cfg to run tag and probe ntuple for muon POG. It runs both on AOD and miniAOD
Modified by Andre Frankenthal (a.franken@cern.ch) -- September 2020
usage: cmsRun run_muonAnalyzer_cfg.py option1=value1 option2=value2
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

options.register('isHIUPC', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "run HIUPC Muon Analyzers"
)
options.parseArguments()

# defaults

if options._beenSet['globalTag'] and options.globalTag != '':
    globaltag = options.globalTag
else:
    # globaltag = '102X_dataRun2_v11' if not options.isMC else '102X_upgrade2018_realistic_v15'
    globaltag = '103X_dataRun2_Prompt_LowPtPhotonReg_v1' if not options.isMC else '102X_upgrade2018_realistic_v15'

# Run local test if no input files provided
if len(options.inputFiles) == 0:
    if options.isMC:
       if options.mcType == "TauTau":
          options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/ggTauTau_TuneCP5_5p02TeV_amcatnlo_pythia8/AODSIM/NoPUlowPtPhotonReg_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v2/280000/0578D671-11C5-8843-88C5-67E784ACB9E0.root')
       if options.mcType == "MuMu":
           options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/02762AC2-E426-CA46-A9BE-3B67E3ABC372.root')
           options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/02762AC2-E426-CA46-A9BE-3B67E3ABC372.root')
           options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/02D64F67-B2B4-494E-B9B1-5D7B0DD22244.root')
           options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/04B4C8C5-DE0D-6146-AC1C-F27DAD4EE6C0.root')
           options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/04CD57C1-5BBA-A94A-848D-8B036B7786BE.root')
           options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/GammaGammatoMuMu_5p02TeV_STARlight/AODSIM/NoPU_103X_upgrade2018_realistic_HI_v11-v2/280000/0659B805-B026-6444-B784-989A8663E582.root')
       if options.mcType == "EE":
          options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/QEDGammaGamma_5p02TeV_STARlight/AODSIM/NoPUlowPtPhotonReg_LbyL_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2-v5/270000/00F70C94-633E-EB45-807C-2FA338E5401F.root')
       if options.mcType == "Other":
          options.inputFiles.append('/store/himc/HINPbPbAutumn18DR/ggBBbar_4f_TuneCP5_5p02TeV_MG5_aMCatNLO_pythia8/AODSIM/NoPUlowPtPhotonReg_103X_upgrade2018_realistic_HI_LowPtPhotonReg_v2_ext1-v2/230000/03E94DC0-FE38-4E49-9445-CA999E866E6C.root')              
    else:
        # options.inputFiles.append('/store/hidata/HIRun2018A/HIForward/AOD/04Apr2019-v1/00000/36696797-7247-4F47-843D-A2D3630801C5.root')   
        options.inputFiles.append('/store/hidata/HIRun2018A/HIForward/AOD/ForLByL-v2/240000/AF5F01C7-6411-A442-96B0-D9672D8830E8.root')   


options.outputFile=""

if options.outputFile == ".root":
    if options.isMC:
        options.outputFile="_" + options.mcType + "_mc" + options.outputFile
    else:
        options.outputFile="_data"+ options.outputFile
        
    if options.isFullAOD:
        options.outputFile="output_full"+ options.outputFile
    else:
        options.outputFile+="output_mini"+ options.outputFile

print ('Output Root file:')
print (options.outputFile)


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
  maxEventsToPrint = cms.untracked.int32(1000),
  printVertex = cms.untracked.bool(False),
  printOnlyHardInteraction = cms.untracked.bool(False), # Print only status=3 particles. This will not work for Pythia8, which does not have any such particles.
  src = cms.InputTag("genParticles")
)

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
  if not options.isMC and options.isHIUPC:
    process.analysis_step = cms.Path(
        process.muonL1Info +
        process.muonL1InfoByQ +
        process.zdcdigi +
        process.QWzdcreco +
        process.muSequence
    )
  else:
    process.analysis_step = cms.Path(

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

process.schedule = cms.Schedule(process.analysis_step, process.endjob_step)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
