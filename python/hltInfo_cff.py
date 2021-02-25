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
"""

hltInfoAll = {
  'Z': {
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
        "HLT_Mu8_v",
        "HLT_Mu17_v",
        "HLT_Mu19_v",
        "HLT_Mu20_v",
        "HLT_IsoMu20_v",
        "HLT_IsoMu24_v",
        "HLT_Mu50_v",
      ],
      'tagFilters': [
        "HLT_Mu8_v",
        "HLT_Mu17_v",
        "HLT_Mu19_v",
        "HLT_Mu20_v",
        "HLT_IsoMu20_v",
        "HLT_IsoMu24_v",
        "HLT_Mu50_v",
        "HLT_OldMu100_v",
        "HLT_TkMu100_v",
        "hltL1fL1sMu22L1Filtered0",
        "hltL2fL1sSingleMu22L1f0L2Filtered10Q",
        "hltL3fL1sSingleMu22L1f0L2f10QL3Filtered24Q",
        "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07",
      ],
      'probeFilters': [
        "HLT_Mu8_v",
        "HLT_Mu17_v",
        "HLT_Mu19_v",
        "HLT_Mu20_v",
        "HLT_IsoMu20_v",
        "HLT_IsoMu24_v",
        "HLT_Mu50_v",
        "HLT_OldMu100_v",
        "HLT_TkMu100_v",
        "hltL1fL1sMu22L1Filtered0",
        "hltL2fL1sSingleMu22L1f0L2Filtered10Q",
        "hltL3fL1sSingleMu22L1f0L2f10QL3Filtered24Q",
        "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07",
      ]
    },
    'Run2018': {
      'triggerPaths': [
        "HLT_Mu8_v",
        "HLT_Mu17_v",
        "HLT_Mu19_v",
        "HLT_Mu20_v",
        "HLT_IsoMu20_v",
        "HLT_IsoMu24_v",
        "HLT_Mu50_v",
      ],
      'tagFilters': [
        "HLT_Mu8_v",
        "HLT_Mu17_v",
        "HLT_Mu19_v",
        "HLT_Mu20_v",
        "HLT_IsoMu20_v",
        "HLT_IsoMu24_v",
        "HLT_Mu50_v",
        "HLT_OldMu100_v",
        "HLT_TkMu100_v",
        "hltL1fL1sMu22L1Filtered0",
        "hltL2fL1sSingleMu22L1f0L2Filtered10Q",
        "hltL3fL1sSingleMu22L1f0L2f10QL3Filtered24Q",
        "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07",
      ],
      'probeFilters': [
        "HLT_Mu8_v",
        "HLT_Mu17_v",
        "HLT_Mu19_v",
        "HLT_Mu20_v",
        "HLT_IsoMu20_v",
        "HLT_IsoMu24_v",
        "HLT_Mu50_v",
        "HLT_OldMu100_v",
        "HLT_TkMu100_v",
        "hltL1fL1sMu22L1Filtered0",
        "hltL2fL1sSingleMu22L1f0L2Filtered10Q",
        "hltL3fL1sSingleMu22L1f0L2f10QL3Filtered24Q",
        "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07",
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

def getHLTInfo(resonance, era):
  check_size(hltInfoAll)
  return hltInfoAll[resonance][getShortEraForHLT(era)]

