# MuonAnalysis-MuonAnalyzer

Package to run tag/probe ntuples for Muon POG on both AOD and miniAOD format.

Requires CMSSW_10_6_X or higher.

## Setup
```
cmsrel CMSSW_10_6_18 
cd CMSSW_10_6_18/src
cmsenv
git cms-init
git clone https://gitlab.cern.ch/cms-muonPOG/muonanalysis-muonanalyzer.git MuonAnalysis/MuonAnalyzer
scram b -j 8
```

## Usage
```
cmsRun MuonAnalysis/MuonAnalyzer/test/run_muonAnalyzer_cfg.py
```
