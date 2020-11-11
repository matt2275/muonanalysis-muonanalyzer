'''author: g. karathanasis
parameters for miniAOD'''

import FWCore.ParameterSet.Config as cms


Path=["HLT_Mu8_v","HLT_Mu17_v","HLT_Mu19_v","HLT_Mu20_v","HLT_IsoMu20_v","HLT_IsoMu24_v","HLT_Mu50_v"]  #paths for tag muon

TagTriggerMatching = [
  # "HLT_Mu8_v",
  # "HLT_Mu17_v",
  # "HLT_Mu19_v",
  # "HLT_Mu20_v",
  # "HLT_IsoMu20_v",
  "HLT_IsoMu24_v",
  "HLT_Mu50_v",
  "HLT_OldMu100_v",
  "HLT_TkMu100_v",
  "hltL1fL1sMu22L1Filtered0",
  "hltL2fL1sSingleMu22L1f0L2Filtered10Q",
  "hltL3fL1sSingleMu22L1f0L2f10QL3Filtered24Q",
  "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07"
]

ProbeTriggerMatching = [
  # "HLT_Mu8_v",
  # "HLT_Mu17_v",
  # "HLT_Mu19_v",
  # "HLT_Mu20_v",
  # "HLT_IsoMu20_v",
  "HLT_IsoMu24_v",
  "HLT_Mu50_v",
  "HLT_OldMu100_v",
  "HLT_TkMu100_v",
  "hltL1fL1sMu22L1Filtered0",
  "hltL2fL1sSingleMu22L1f0L2Filtered10Q",
  "hltL3fL1sSingleMu22L1f0L2f10QL3Filtered24Q",
  "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07"
]

muon = cms.EDAnalyzer('MuonMiniAODAnalyzer',
           beamSpot=cms.InputTag('offlineBeamSpot'),
           vertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
           muons=cms.InputTag("slimmedMuons"),
           triggerResults=cms.InputTag("TriggerResults::HLT"),
           l1Matches = cms.InputTag("muonL1Info"),
           l1MatchesQuality = cms.InputTag("muonL1Info", "quality"),
           l1MatchesDeltaR = cms.InputTag("muonL1Info", "deltaR"),
           l1MatchesByQ = cms.InputTag("muonL1InfoByQ"),
           l1MatchesByQQuality = cms.InputTag("muonL1InfoByQ", "quality"),
           l1MatchesByQDeltaR = cms.InputTag("muonL1InfoByQ", "deltaR"),
           PFCands=cms.InputTag("packedPFCandidates"),
           lostTracks=cms.InputTag("lostTracks"),
           gen = cms.InputTag("prunedGenParticles"),
           HLTPaths=cms.vstring(Path),
           TagPathsOrFilters=cms.vstring(TagTriggerMatching),
           ProbePathsOrFilters=cms.vstring(ProbeTriggerMatching),
           tagQuality = cms.uint32(0), # quality of tag muon following muonSelector convention
           tagSelection = cms.string("pt()>0"), # string to pass cuts on tag
           ProbeHPurity = cms.bool(True), # skips non High purity probes
           probeSelection = cms.string("pt()>0"), #string for probe
           pairMassMin = cms.double(60.0), # min mass of mu pair
           pairMassMax = cms.double(140.0), # max mss of mu pair
           pairDz = cms.double(1e99), #max Dz of mu1,mu2
           RequireVtxCreation = cms.bool(False), # if true kills pairs w/o vtx
           minSVtxProb = cms.double(-0.01), # min prob of mu pair
           maxDzProbeTrkMuon = cms.double(0.01), # max Dz(mu1,mu2)
           maxRelPtProbeTrkMuon = cms.double(1.0),# max [pt(mu)-pt(trk)]/pt(trk) for probe/offline
           maxDRProbeTrkMuon =  cms.double(0.03), # max DR for probe/offline
           momPdgId = cms.uint32(23),
           genRecoDrMatch= cms.double(0.03)
           
)

miniAODSequence=cms.Sequence(muon)
