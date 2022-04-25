'''author gkarathanasis
option for AOD run'''

import FWCore.ParameterSet.Config as cms

muon = cms.EDAnalyzer('StandAloneMuonFullAODAnalyzer',
        isMC=cms.bool(False),
        includeJets=cms.bool(False),
        era = cms.string('dummy'), # updated in run_muonAnalyzer_cfg.py
        genEventInfo = cms.InputTag('generator'),
        pileupInfo=cms.InputTag('addPileupInfo'),
        Rho=cms.InputTag('fixedGridRhoFastjetAll'),
        beamSpot=cms.InputTag('offlineBeamSpot'),
        vertices=cms.InputTag("offlinePrimaryVertices"),
        muons=cms.InputTag("muons"),
        tracks=cms.InputTag("generalTracks"),
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
        trgDRwindow= cms.double(0.1), # dr winwow hlt mu/offline
        tagQuality = cms.uint32(0),
        tagSelection = cms.string("pt()>5"),
        probeHPurity = cms.bool(False),
        probeSelection = cms.string("pt()>2"),
        probeSelectionSA = cms.string("pt()>2"),
        muonOnly = cms.bool(False), # allow only reco or pat Muon for probes
        probeMuonSelection = cms.string("pt()>0"), #string for probe (reco or pat Muon)
        pairMassMin = cms.double(2.8),
        pairMassMax = cms.double(3.4), # 9999.0 for high mass C&C
        pairDz = cms.double(10.1),
        RequireVtxCreation = cms.bool(True),
        minSVtxProb = cms.double(0),
        maxDzProbeTrkMuon = cms.double(0.01), # max Dz(mu1,mu2)
        maxRelPtProbeTrkMuon = cms.double(1.0),# max [pt(mu)-pt(trk)]/pt(trk) for probe/offline
        maxDRProbeTrkMuon =  cms.double(0.03), # max DR for probe/offline
        minPtTrkSA=cms.double(10.0),  # min Pt for trk matching for SA muons
        maxDzProbeTrkSAMuon = cms.double(1.0), # max Dz(trk,SAmu)
        maxRelPtProbeTrkSAMuon = cms.double(1.0),# max [pt(SAmu)-pt(trk)]/pt(trk) for probe/offline
        maxDRProbeTrkSAMuon =  cms.double(0.3), # max DR for trk / SA muon
        maxDRProbeTrkDSA =  cms.double(0.5), # max DR for general track and dSA
        momPdgId= cms.uint32(443),
        genRecoDrMatch = cms.double(0.03),
        saveStandAloneTree = cms.bool(True),
        saveTnPTree = cms.bool(True),
        debug = cms.int32(0),
        propM1 = cms.PSet(
            useStation2 = cms.bool(False),
            useTrack = cms.string("tracker"),
            useState = cms.string("atVertex"),  # in AOD
            useSimpleGeometry = cms.bool(True), # use just one cylinder and two planes, not all the fancy chambers  
        ),                                                        
)

fullAODSequence=cms.Sequence(muon)
