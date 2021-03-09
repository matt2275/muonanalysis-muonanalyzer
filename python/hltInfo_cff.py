"""
Dictionary for trigger information

hltInfoAll[resonance][era][type]:
- resonance: Z or JPsi
- era: Run2016, Run2017, Run2018
- type:
  - triggerPaths: HLT paths accepted in a givin event
  - tagFilters: HLT filters matched to the tag muon (MiniAOD: path names are also supported),
                the tag should pass at least one of the filters to be stored
  - probeFilters: HLT filters matched to the probe muon (MiniAOD: path names are also supported)
* Recommended naming convention:
  - path name: HLT_X_v, e.g. HLT_IsoMu24_v, HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v
  - filter name: full filter name, e.g. hltL3fL1sSingleMu22L1f0L2f10QL3Filtered24Q
"""

hltInfoAll = {
  'Z': {
    'Run2016': {
      'triggerPaths': [
        "HLT_IsoMu24_v",
        "HLT_IsoTkMu24_v",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
        "HLT_TkMu17_v",
        "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
        "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
        "HLT_Mu50_v",
        "HLT_TkMu50_v",
        "HLT_Mu17_TrkIsoVVL_v",
        "HLT_Mu17_v",
      ],
      'tagFilters': [
        'HLT_IsoMu24_v',
        'hltL2fL1sMu22L1f0L2Filtered10Q',
        'hltL3fL1sMu22L1f0L2f10QL3Filtered24Q',
        'hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09',
        'HLT_IsoTkMu24_v',
        'hltL3fL1sMu22f0TkFiltered24Q',
        'hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09',
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v',
        'hltL2pfL1sDoubleMu114ORDoubleMu125L1f0L2PreFiltered0',
        'hltL2fL1sDoubleMu114L1f0L2Filtered10OneMu',
        'hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8',
        'hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17',
        'hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4',
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v',
        'hltL2pfL1sDoubleMu114ORDoubleMu125L1f0L2PreFiltered0',
        'hltL2fL1sDoubleMu114L1f0L2Filtered10OneMu',
        'hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8',
        'hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17',
        'hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4',
        'hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2',
        'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v',
        'hltL2fL1sDoubleMu114L1f0OneMuL2Filtered10',
        'hltL3fL1sDoubleMu114L1f0L2f10L3Filtered17',
        'hltDiMuonGlbFiltered17TrkFiltered8',
        'hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4',
        'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v',
        'hltL2fL1sDoubleMu114L1f0OneMuL2Filtered10',
        'hltL3fL1sDoubleMu114L1f0L2f10L3Filtered17',
        'hltDiMuonGlbFiltered17TrkFiltered8',
        'hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4',
        'hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2',
        'HLT_TkMu17_v',
        'hltL3fL1sMu10lqTkFiltered17Q',
        'HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v',
        'hltL3fL1sDoubleMu114TkFiltered17Q',
        'hltDiTkMuonTkFiltered17TkFiltered8',
        'hltDiMuonTrk17Trk8RelTrkIsoFiltered0p4',
        'HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v',
        'hltL3fL1sDoubleMu114TkFiltered17Q',
        'hltDiTkMuonTkFiltered17TkFiltered8',
        'hltDiMuonTrk17Trk8RelTrkIsoFiltered0p4',
        'hltDiMuonTrk17Trk8RelTrkIsoFiltered0p4DzFiltered0p2',
        'HLT_Mu50_v',
        'hltL2fL1sMu22Or25L1f0L2Filtered10Q',
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q',
        'HLT_TkMu50_v',
        'hltL3fL1sMu25f0TkFiltered50Q',
        'HLT_Mu17_TrkIsoVVL_v',
        'hltL2fL1sMu10lqL1f0L2Filtered10',
        'hltL3fL1sMu10lqL1f0L2f10L3Filtered17',
        'hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4',
        'HLT_Mu17_v',
        'hltL2fL1sMu10lqL1f0L2Filtered10',
        'hltL3fL1sMu10lqL1f0L2f10L3Filtered17',
      ],
      'probeFilters': [
        'HLT_IsoMu24_v',
        'hltL2fL1sMu22L1f0L2Filtered10Q',
        'hltL3fL1sMu22L1f0L2f10QL3Filtered24Q',
        'hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09',
        'HLT_IsoTkMu24_v',
        'hltL3fL1sMu22f0TkFiltered24Q',
        'hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09',
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v',
        'hltL2pfL1sDoubleMu114ORDoubleMu125L1f0L2PreFiltered0',
        'hltL2fL1sDoubleMu114L1f0L2Filtered10OneMu',
        'hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8',
        'hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17',
        'hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4',
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v',
        'hltL2pfL1sDoubleMu114ORDoubleMu125L1f0L2PreFiltered0',
        'hltL2fL1sDoubleMu114L1f0L2Filtered10OneMu',
        'hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8',
        'hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17',
        'hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4',
        'hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2',
        'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v',
        'hltL2fL1sDoubleMu114L1f0OneMuL2Filtered10',
        'hltL3fL1sDoubleMu114L1f0L2f10L3Filtered17',
        'hltDiMuonGlbFiltered17TrkFiltered8',
        'hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4',
        'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v',
        'hltL2fL1sDoubleMu114L1f0OneMuL2Filtered10',
        'hltL3fL1sDoubleMu114L1f0L2f10L3Filtered17',
        'hltDiMuonGlbFiltered17TrkFiltered8',
        'hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4',
        'hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2',
        'HLT_TkMu17_v',
        'hltL3fL1sMu10lqTkFiltered17Q',
        'HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v',
        'hltL3fL1sDoubleMu114TkFiltered17Q',
        'hltDiTkMuonTkFiltered17TkFiltered8',
        'hltDiMuonTrk17Trk8RelTrkIsoFiltered0p4',
        'HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v',
        'hltL3fL1sDoubleMu114TkFiltered17Q',
        'hltDiTkMuonTkFiltered17TkFiltered8',
        'hltDiMuonTrk17Trk8RelTrkIsoFiltered0p4',
        'hltDiMuonTrk17Trk8RelTrkIsoFiltered0p4DzFiltered0p2',
        'HLT_Mu50_v',
        'hltL2fL1sMu22Or25L1f0L2Filtered10Q',
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q',
        'HLT_TkMu50_v',
        'hltL3fL1sMu25f0TkFiltered50Q',
        'HLT_Mu17_TrkIsoVVL_v',
        'hltL2fL1sMu10lqL1f0L2Filtered10',
        'hltL3fL1sMu10lqL1f0L2f10L3Filtered17',
        'hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4',
        'HLT_Mu17_v',
        'hltL2fL1sMu10lqL1f0L2Filtered10',
        'hltL3fL1sMu10lqL1f0L2f10L3Filtered17',
      ]
    },
    'Run2017': {
      'triggerPaths': [
        "HLT_IsoMu27_v",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",
        "HLT_Mu50_v",
        "HLT_OldMu100_v",
        "HLT_TkMu100_v",
        "HLT_Mu17_TrkIsoVVL_v",
        "HLT_Mu17_v",
      ],
      'tagFilters': [
        'HLT_IsoMu27_v',
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q',
        'hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07',
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v',
        'hltL3fL1DoubleMu155fPreFiltered8',
        'hltL3fL1DoubleMu155fFiltered17',
        'hltDiMuon178RelTrkIsoFiltered0p4',
        'hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2',
        'hltDiMuon178Mass8Filtered',
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v',
        'hltL3fL1DoubleMu155fPreFiltered8',
        'hltL3fL1DoubleMu155fFiltered17',
        'hltDiMuon178RelTrkIsoFiltered0p4',
        'hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2',
        'hltDiMuon178Mass3p8Filtered',
        'HLT_Mu50_v',
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q',
        'HLT_OldMu100_v',
        'hltL2fOldL1sMu22or25L1f0L2Filtered10Q',
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered100Q',
        'HLT_TkMu100_v',
        'hltL3fL1sMu25f0TkFiltered100Q',
        'HLT_Mu17_TrkIsoVVL_v',
        'hltL3fL1sMu10lqL1f0L2f10L3Filtered17',
        'hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4',
        'HLT_Mu17_v',
        'hltL3fL1sMu10lqL1f0L2f10L3Filtered17',
      ],
      'probeFilters': [
        'HLT_IsoMu27_v',
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q',
        'hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07',
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v',
        'hltL3fL1DoubleMu155fPreFiltered8',
        'hltL3fL1DoubleMu155fFiltered17',
        'hltDiMuon178RelTrkIsoFiltered0p4',
        'hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2',
        'hltDiMuon178Mass8Filtered',
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v',
        'hltL3fL1DoubleMu155fPreFiltered8',
        'hltL3fL1DoubleMu155fFiltered17',
        'hltDiMuon178RelTrkIsoFiltered0p4',
        'hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2',
        'hltDiMuon178Mass3p8Filtered',
        'HLT_Mu50_v',
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q',
        'HLT_OldMu100_v',
        'hltL2fOldL1sMu22or25L1f0L2Filtered10Q',
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered100Q',
        'HLT_TkMu100_v',
        'hltL3fL1sMu25f0TkFiltered100Q',
        'HLT_Mu17_TrkIsoVVL_v',
        'hltL3fL1sMu10lqL1f0L2f10L3Filtered17',
        'hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4',
        'HLT_Mu17_v',
        'hltL3fL1sMu10lqL1f0L2f10L3Filtered17',
      ]
    },
    'Run2018': {
      'triggerPaths': [
        "HLT_IsoMu24_v",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",
        "HLT_Mu50_v",
        "HLT_OldMu100_v",
        "HLT_TkMu100_v",
        "HLT_Mu17_TrkIsoVVL_v",
        "HLT_Mu17_v",
      ],
      'tagFilters': [
        'HLT_IsoMu24_v',
        'hltL3fL1sSingleMu22L1f0L2f10QL3Filtered24Q',
        'hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07',
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v',
        'hltL3fL1DoubleMu155fPreFiltered8',
        'hltL3fL1DoubleMu155fFiltered17',
        'hltDiMuon178RelTrkIsoFiltered0p4',
        'hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2',
        'hltDiMuon178Mass3p8Filtered',
        'HLT_Mu50_v',
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q',
        'HLT_OldMu100_v',
        'hltL2fOldL1sMu22or25L1f0L2Filtered10Q',
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered100Q',
        'HLT_TkMu100_v',
        'hltL3fL1sMu25f0TkFiltered100Q',
        'HLT_Mu17_TrkIsoVVL_v',
        'hltL3fL1sMu15DQlqL1f0L2f10L3Filtered17',
        'hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4',
        'HLT_Mu17_v',
        'hltL3fL1sMu15DQlqL1f0L2f10L3Filtered17',
      ],
      'probeFilters': [
        'HLT_IsoMu24_v',
        'hltL3fL1sSingleMu22L1f0L2f10QL3Filtered24Q',
        'hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07',
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v',
        'hltL3fL1DoubleMu155fPreFiltered8',
        'hltL3fL1DoubleMu155fFiltered17',
        'hltDiMuon178RelTrkIsoFiltered0p4',
        'hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2',
        'hltDiMuon178Mass3p8Filtered',
        'HLT_Mu50_v',
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q',
        'HLT_OldMu100_v',
        'hltL2fOldL1sMu22or25L1f0L2Filtered10Q',
        'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered100Q',
        'HLT_TkMu100_v',
        'hltL3fL1sMu25f0TkFiltered100Q',
        'HLT_Mu17_TrkIsoVVL_v',
        'hltL3fL1sMu15DQlqL1f0L2f10L3Filtered17',
        'hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4',
        'HLT_Mu17_v',
        'hltL3fL1sMu15DQlqL1f0L2f10L3Filtered17',
      ]
    }
  },
  'JPsi': {
    'Run2016': {
      'triggerPaths': [
      ],
      'tagFilters': [
      ],
      'probeFilters': [
      ]
    },
    'Run2017': {
      'triggerPaths': [
      ],
      'tagFilters': [
      ],
      'probeFilters': [
      ]
    },
    'Run2018': {
      'triggerPaths': [
        "HLT_Mu7p5_Track7_Jpsi",
        "HLT_Mu7p5_Track3p5_Jpsi",
        "HLT_Mu7p5_Track2_Jpsi"
      ],
      'tagFilters': [
        'hltL3fLMu7p5TrackL3Filtered7p5',
      ],
      'probeFilters': [
        'hltL3fLMu7p5TrackL3Filtered7p5',
      ]
    }
  }
}

def check_size(info, max_size = 100, keys = []):
  for k, v in info.items():
    assert isinstance(v, dict) or isinstance(v, list), v
    if isinstance(v, dict):
      keys.append(k)
      check_size(v, max_size, keys)
    else:
      if len(v) > max_size:
        loc = '|'.join(keys)
        print "WARNING either put less than {} paths/filters in hltInfoAll|{}|{}, or increase the quota from NtupleContent.h/.cc".format(max_size+1, loc, k)
        exit()
    if len(keys) > 0 and k == info.keys()[-1]:
      keys.pop()

def getShortEraForHLT(era):
  if 'Run2016' in era:
    return 'Run2016'
  elif 'Run2017' in era:
    return 'Run2017'
  else:  # any other era will use the trigger list for Run2018
    return 'Run2018'

def selectTriggers(trgList, keepPaths = True, keepFilters = True):
  assert (keepFilters or keepPaths)
  if keepFilters and keepPaths:
    return trgList

  out = []
  for trg in trgList:
    if keepFilters and not trg.startswith('hlt'):
      continue
    if keepPaths and not trg.startswith('HLT_'):
      continue
    out.append(trg)
  return out

def getHLTInfo(resonance, era):
  check_size(hltInfoAll)
  return hltInfoAll[resonance][getShortEraForHLT(era)]

