'''author gkarathanasis
option for AOD run'''

import FWCore.ParameterSet.Config as cms

muon = cms.EDAnalyzer('HIUPC_Analysis_GenEfficiency',
        isMC=cms.bool(False),
        includeJets=cms.bool(False),
        era = cms.string('dummy'), # updated in run_muonAnalyzer_cfg.py
        genEventInfo = cms.InputTag('generator'),
        pileupInfo=cms.InputTag('addPileupInfo'),
        Rho=cms.InputTag('fixedGridRhoFastjetAll'),
        beamSpot=cms.InputTag('offlineBeamSpot'),
        vertices=cms.InputTag("offlinePrimaryVerticesRecovery"),
        muons=cms.InputTag("muons"),
        electrons=cms.InputTag("gedGsfElectrons"),
        photons=cms.InputTag("gedPhotons"),
        tracks=cms.InputTag("generalTracks"),
        CaloTowers = cms.InputTag("towerMaker"),
        # Centrality= cms.InputTag("hiCentrality"),
        CentralitySrc    = cms.InputTag("hiCentrality"),
        CentralityBinSrc = cms.InputTag("centralityBin","HFtowers"),
        RecHits = cms.InputTag("QWzdcreco"),        
        SAmuons=cms.InputTag("standAloneMuons"),
        dSAmuons=cms.InputTag("displacedStandAloneMuons"),
        dGlmuons=cms.InputTag("displacedGlobalMuons"),
        staCosmic=cms.InputTag("cosmicMuons"),
        triggerResults=cms.InputTag("TriggerResults::HLT"),
        triggerObjects=cms.InputTag('hltTriggerSummaryAOD'),
        l1Matches = cms.InputTag("muonL1Info"),
        l1MatchesQuality = cms.InputTag("muonL1Info", "quality"),
        l1MatchesDeltaR = cms.InputTag("muonL1Info", "deltaR"),
        l1MatchesByQ = cms.InputTag("muonL1InfoByQ"),
        l1MatchesByQQuality = cms.InputTag("muonL1InfoByQ", "quality"),
        l1MatchesByQDeltaR = cms.InputTag("muonL1InfoByQ", "deltaR"),
        muonSimInfo = cms.InputTag("muonSimClassifier"),
        triggerPaths=cms.vstring(), # updated in run_muonAnalyzer_cfg.py
        tagFilters=cms.vstring(), # updated in run_muonAnalyzer_cfg.py
        probeFilters=cms.vstring(), # updated in run_muonAnalyzer_cfg.py
        probeSelectorNames=cms.vstring(), # updated in run_muonAnalyzer_cfg.py
        probeSelectorBits=cms.vuint32(), # updated in run_muonAnalyzer_cfg.py
        gen = cms.InputTag("genParticles"),
        rhoJetsNC = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
        PFCands = cms.InputTag("particleFlow"),
        jets = cms.InputTag("ak4PFJetsCHS"),
        jetCorrector = cms.InputTag("ak4PFCHSL1FastL2L3Corrector"),
        genJets = cms.InputTag("ak4GenJets"),
        deepCSVProbb = cms.InputTag("pfDeepCSVJetTags:probb"),
        deepCSVProbbb = cms.InputTag("pfDeepCSVJetTags:probbb"),
        deepFlavProbb = cms.InputTag("pfDeepFlavourJetTags:probb"),
        deepFlavProbbb = cms.InputTag("pfDeepFlavourJetTags:probbb"),
        keepMuons = cms.bool(True),
        keepElectrons = cms.bool(True),
        keepTracks = cms.bool(True),
        keepPFcands = cms.bool(True),
        keepPhotons = cms.bool(True),
        keepCaloTowers = cms.bool(True),
        keepZDC = cms.bool(True),
        debug = cms.int32(0),
        MCType = cms.string(""), # filled in test_run file                                                      
)

fullAODSequence=cms.Sequence(muon)
