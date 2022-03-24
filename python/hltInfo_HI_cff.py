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

        "HLT_HIUPCSingleMuNotHF2Pixel_SingleTrack_v",
      ],
      'tagFilters': [
 
        'hltL1MuOpenNotHF2L1Filtered0',

      ],
      'probeFilters': [
        'hltL1MuOpenNotHF2L1Filtered0',
      ]
    },
    'Run2017': {
      'triggerPaths': [

        "HLT_HIUPCSingleMuNotHF2Pixel_SingleTrack_v",
      ],
      'tagFilters': [
 
        'hltL1MuOpenNotHF2L1Filtered0',

      ],
      'probeFilters': [
        'hltL1MuOpenNotHF2L1Filtered0',
      ]
    },
    'Run2018': {
      'triggerPaths': [
        "HLT_HIUPC_SingleMuOpen_NotMBHF2AND_v",
        "HLT_HIUPC_SingleMu0_NotMBHF2AND_v1",
      ],
      'tagFilters': [
        'hltL1sSingleMuOpenNotMBHF2AND',
        'hltL1sSingleMu0NotMBHF2AND',
      ],
      'probeFilters': [
        'hltL1sSingleMuOpenNotMBHF2AND',
        'hltL1sSingleMu0NotMBHF2AND',

      ]
    }
  },
  'JPsi': {
    'Run2016': {
      'triggerPaths': [

        "HLT_HIUPCSingleMuNotHF2Pixel_SingleTrack_v",
      ],
      'tagFilters': [
 
        'hltL1MuOpenNotHF2L1Filtered0',

      ],
      'probeFilters': [
        'hltL1MuOpenNotHF2L1Filtered0',
      ]
    },
    'Run2017': {
      'triggerPaths': [

        "HLT_HIUPCSingleMuNotHF2Pixel_SingleTrack_v",
      ],
      'tagFilters': [
 
        'hltL1MuOpenNotHF2L1Filtered0',

      ],
      'probeFilters': [
        'hltL1MuOpenNotHF2L1Filtered0',
      ]
    },
    'Run2018': {
      'triggerPaths': [
        "HLT_HIUPC_SingleMuOpen_NotMBHF2AND_v",
        "HLT_HIUPC_SingleMu0_NotMBHF2AND_v1",
      ],
      'tagFilters': [
        'hltL1sSingleMuOpenNotMBHF2AND',
        'hltL1sSingleMu0NotMBHF2AND',
      ],
      'probeFilters': [
        'hltL1sSingleMuOpenNotMBHF2AND',
        'hltL1sSingleMu0NotMBHF2AND',

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

def selectTriggers(trgList, keepPaths = True, keepFilters = True, excludeDSA = False):
  assert (keepFilters or keepPaths)
  if keepFilters and keepPaths and not excludeDSA:
    return trgList

  out = []
  for trg in trgList:
    if trg in out:
        continue
    if excludeDSA and ('NoVtx' in trg or 'NoVertex' in trg):
        continue
    if keepFilters and trg.startswith('hlt'):
      out.append(trg)
    if keepPaths and trg.startswith('HLT_'):
      out.append(trg)
  return out

def getHLTInfo(resonance, era):
  check_size(hltInfoAll)
  return hltInfoAll[resonance][getShortEraForHLT(era)]

