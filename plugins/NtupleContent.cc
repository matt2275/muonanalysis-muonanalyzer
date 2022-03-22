#include "NtupleContent.h"

NtupleContent::NtupleContent() {}

NtupleContent::~NtupleContent() {}

void NtupleContent::SetTree(TTree *mytree) { t1 = mytree; }

void NtupleContent::CreateBranches(const std::vector<std::string> &HLTs,
                                   const std::vector<std::string> &selectorNames) {
  // General
  t1->Branch("run", &run);
  t1->Branch("event", &event);
  t1->Branch("lumi", &ls);
  t1->Branch("fromFullAOD", &fromFullAOD);
  t1->Branch("genWeight", &genWeight);
  t1->Branch("BSpot_x", &BSpot_x);
  t1->Branch("BSpot_y", &BSpot_y);
  t1->Branch("BSpot_z", &BSpot_z);
  t1->Branch("pv_x", &pv_x);
  t1->Branch("pv_y", &pv_y);
  t1->Branch("pv_z", &pv_z);
  t1->Branch("nVertices", &nvertices);
  t1->Branch("nTrueInteractions", &trueNumInteractions);
  t1->Branch("nPUInteractions", &puNumInteractions);
  t1->Branch("rho", &Rho);
  t1->Branch("nmuons", &nmuons);
  t1->Branch("ntag", &ntag);
  t1->Branch("npairs", &npairs);
  t1->Branch("pair_probeMultiplicity", &pair_probeMultiplicity);
  t1->Branch("genmu1_pt", &genmu1_pt);
  t1->Branch("genmu1_eta", &genmu1_eta);
  t1->Branch("genmu1_phi", &genmu1_phi);
  t1->Branch("genmu1_charge", &genmu1_charge);
  t1->Branch("genmu2_pt", &genmu2_pt);
  t1->Branch("genmu2_eta", &genmu2_eta);
  t1->Branch("genmu2_phi", &genmu2_phi);
  t1->Branch("genmu2_charge", &genmu2_charge);
  t1->Branch("genMass", &genMass);
  t1->Branch("genmuFSfromHP1_pt", &genmuFSfromHP1_pt);
  t1->Branch("genmuFSfromHP1_eta", &genmuFSfromHP1_eta);
  t1->Branch("genmuFSfromHP1_phi", &genmuFSfromHP1_phi);
  t1->Branch("genmuFSfromHP1_charge", &genmuFSfromHP1_charge);
  t1->Branch("genmuFSfromHP2_pt", &genmuFSfromHP2_pt);
  t1->Branch("genmuFSfromHP2_eta", &genmuFSfromHP2_eta);
  t1->Branch("genmuFSfromHP2_phi", &genmuFSfromHP2_phi);
  t1->Branch("genmuFSfromHP2_charge", &genmuFSfromHP2_charge);
  t1->Branch("genMassFSfromHP", &genMassFSfromHP);
  for (unsigned int ihlt = 0; ihlt < HLTs.size(); ihlt++)
    t1->Branch(TString(HLTs[ihlt]), &trigger[ihlt]);
  t1->Branch("genJets_pt", &genJets_pt);
  t1->Branch("genJets_eta", &genJets_eta);
  t1->Branch("genJets_phi", &genJets_phi);
  t1->Branch("genJets_mass", &genJets_mass);
  t1->Branch("jets_pt", &jets_pt);
  t1->Branch("jets_eta", &jets_eta);
  t1->Branch("jets_phi", &jets_phi);
  t1->Branch("jets_mass", &jets_mass);
  t1->Branch("jets_isTight", &jets_isTight);
  t1->Branch("jets_isTightLepVeto", &jets_isTightLepVeto);
  t1->Branch("jets_bTag_deepCSV", &jets_bTag_deepCSV);
  t1->Branch("jets_bTag_deepFlav", &jets_bTag_deepFlav);
  t1->Branch("nTightJets", &nTightJets);
  t1->Branch("nTightLepVetoJets", &nTightLepVetoJets);
  // Tag specific
  t1->Branch("tag_pt", &tag_pt);
  t1->Branch("tag_eta", &tag_eta);
  t1->Branch("tag_phi", &tag_phi);
  t1->Branch("tag_charge", &tag_charge);
  t1->Branch("tag_pterr", &tag_pterr);
  t1->Branch("tag_dxy", &tag_dxy);
  t1->Branch("tag_dz", &tag_dz);
  t1->Branch("tag_isPF", &tag_isPF);
  t1->Branch("tag_isSA", &tag_isSA);
  t1->Branch("tag_isdSA", &tag_isdSA);
  t1->Branch("tag_isTracker", &tag_isTracker);
  t1->Branch("tag_isGlobal", &tag_isGlobal);
  t1->Branch("tag_isLoose", &tag_isLoose);
  t1->Branch("tag_isMedium", &tag_isMedium);
  t1->Branch("tag_isTight", &tag_isTight);
  t1->Branch("tag_isSoft", &tag_isSoft);
  t1->Branch("tag_isHighPt", &tag_isHighPt);
  t1->Branch("tag_relIso04", &tag_relIso04);
  t1->Branch("tag_miniIso", &tag_miniIso);
  t1->Branch("tag_miniIsoCharged", &tag_miniIsoCharged);
  t1->Branch("tag_miniIsoPhotons", &tag_miniIsoPhotons);
  t1->Branch("tag_miniIsoNeutrals", &tag_miniIsoPhotons);
  t1->Branch("tag_isMatchedGen", &tag_isMatchedGen);
  t1->Branch("tag_minDR", &tag_minDR);
  t1->Branch("tag_ptRel_minDR", &tag_ptRel_minDR);
  t1->Branch("tag_iso03_sumPt", &tag_iso03_sumPt);
  t1->Branch("tag_pfIso03_charged", &tag_pfIso03_charged);
  t1->Branch("tag_pfIso03_neutral", &tag_pfIso03_neutral);
  t1->Branch("tag_pfIso03_photon", &tag_pfIso03_photon);
  t1->Branch("tag_pfIso03_sumPU", &tag_pfIso03_sumPU);
  t1->Branch("tag_pfIso04_charged", &tag_pfIso04_charged);
  t1->Branch("tag_pfIso04_neutral", &tag_pfIso04_neutral);
  t1->Branch("tag_pfIso04_photon", &tag_pfIso04_photon);
  t1->Branch("tag_pfIso04_sumPU", &tag_pfIso04_sumPU);
  t1->Branch("tag_tuneP_pt", &tag_tuneP_pt);
  t1->Branch("tag_tuneP_pterr", &tag_tuneP_pterr);
  t1->Branch("tag_nsegments", &tag_nsegments);
  // Probe specific
  t1->Branch("iprobe", &iprobe);
  t1->Branch("probe_pt", &probe_pt);
  t1->Branch("probe_eta", &probe_eta);
  t1->Branch("probe_phi", &probe_phi);
  t1->Branch("probe_charge", &probe_charge);
  t1->Branch("probe_pterr", &probe_pterr);
  t1->Branch("probe_dxy", &probe_dxy);
  t1->Branch("probe_dz", &probe_dz);
  t1->Branch("probe_isPF", &probe_isPF);
  t1->Branch("probe_isSA", &probe_isSA);
  t1->Branch("probe_isTracker", &probe_isTracker);
  t1->Branch("probe_isGlobal", &probe_isGlobal);
  t1->Branch("probe_isLoose", &probe_isLoose);
  t1->Branch("probe_isMedium", &probe_isMedium);
  t1->Branch("probe_isTight", &probe_isTight);
  t1->Branch("probe_isSoft", &probe_isSoft);
  t1->Branch("probe_isHighPt", &probe_isHighPt);
  t1->Branch("probe_isArbitratedTracker", &probe_isArbitratedTracker);
  t1->Branch("probe_isMuMatched", &probe_isMuMatched);
  t1->Branch("probe_isdSA", &probe_isdSA);
  t1->Branch("probe_isdGlobal", &probe_isdGlobal);
  t1->Branch("probe_isCosmic", &probe_isCosmic);
  t1->Branch("probe_ncosmic", &probe_ncosmic);
  t1->Branch("probe_cosmic_minDR", &probe_cosmic_minDR);
  //  t1->Branch("probe_isGood", &probe_isGood);
  t1->Branch("probe_isHighPurity", &probe_isHighPurity);
  t1->Branch("probe_validFraction", &probe_validFraction);
  t1->Branch("probe_trkChi2", &probe_trkChi2);
  t1->Branch("probe_positionChi2", &probe_positionChi2);
  t1->Branch("probe_trkKink", &probe_trkKink);
  // t1->Branch("probe_segmentCompatibility", &probe_segmentCompatibility);
  t1->Branch("probe_trackerLayers", &probe_trackerLayers);
  t1->Branch("probe_pixelLayers", &probe_pixelLayers);
  t1->Branch("probe_muonStations", &probe_muonStations);
  t1->Branch("probe_muonHits", &probe_muonHits);
  t1->Branch("probe_DTHits", &probe_DTHits);
  t1->Branch("probe_CSCHits", &probe_CSCHits);
  t1->Branch("probe_relIso04", &probe_relIso04);
  t1->Branch("probe_miniIso", &probe_miniIso);
  t1->Branch("probe_miniIsoCharged", &probe_miniIsoCharged);
  t1->Branch("probe_miniIsoPhotons", &probe_miniIsoPhotons);
  t1->Branch("probe_miniIsoNeutrals", &probe_miniIsoPhotons);
  t1->Branch("probe_isMatchedGen", &probe_isMatchedGen);
  t1->Branch("probe_minDR", &probe_minDR);
  t1->Branch("probe_ptRel_minDR", &probe_ptRel_minDR);
  t1->Branch("probe_iso03_sumPt", &probe_iso03_sumPt);
  t1->Branch("probe_pfIso03_charged", &probe_pfIso03_charged);
  t1->Branch("probe_pfIso03_neutral", &probe_pfIso03_neutral);
  t1->Branch("probe_pfIso03_photon", &probe_pfIso03_photon);
  t1->Branch("probe_pfIso03_sumPU", &probe_pfIso03_sumPU);
  t1->Branch("probe_pfIso04_charged", &probe_pfIso04_charged);
  t1->Branch("probe_pfIso04_neutral", &probe_pfIso04_neutral);
  t1->Branch("probe_pfIso04_photon", &probe_pfIso04_photon);
  t1->Branch("probe_pfIso04_sumPU", &probe_pfIso04_sumPU);
  t1->Branch("probe_pixelHits", &probe_pixelHits);
  t1->Branch("probe_matchedStations", &probe_matchedStations);
  t1->Branch("probe_expectedMatchedStations", &probe_expectedMatchedStations);
  t1->Branch("probe_RPCLayers", &probe_RPCLayers);
  t1->Branch("probe_stationMask", &probe_stationMask);
  t1->Branch("probe_nShowers", &probe_nShowers);
  t1->Branch("probe_tuneP_pt", &probe_tuneP_pt);
  t1->Branch("probe_tuneP_pterr", &probe_tuneP_pterr);
  t1->Branch("probe_tuneP_muonHits", &probe_tuneP_muonHits);
  t1->Branch("probe_nsegments", &probe_nsegments);

  t1->Branch("l1pt", &l1pt);
  t1->Branch("l1q", &l1q);
  t1->Branch("l1dr", &l1dr);
  t1->Branch("l1ptByQ", &l1ptByQ);
  t1->Branch("l1qByQ", &l1qByQ);
  t1->Branch("l1drByQ", &l1drByQ);

  t1->Branch("tag_l1pt", &tag_l1pt);
  t1->Branch("tag_l1q", &tag_l1q);
  t1->Branch("tag_l1dr", &tag_l1dr);
  t1->Branch("tag_l1ptByQ", &tag_l1ptByQ);
  t1->Branch("tag_l1qByQ", &tag_l1qByQ);
  t1->Branch("tag_l1drByQ", &tag_l1drByQ);

  t1->Branch("probe_dsa_segmentMatches", &probe_dsa_segmentMatches);
  t1->Branch("probe_dsa_nsegments", &probe_dsa_nsegments);
  t1->Branch("probe_dsa_muonStations", &probe_dsa_muonStations);
  t1->Branch("probe_dsa_muonHits", &probe_dsa_muonHits);
  t1->Branch("probe_dsa_DTHits", &probe_dsa_DTHits);
  t1->Branch("probe_dsa_CSCHits", &probe_dsa_CSCHits);
  t1->Branch("probe_dsa_pterr", &probe_dsa_pterr);
  t1->Branch("probe_dsa_dxy", &probe_dsa_dxy);
  t1->Branch("probe_dsa_dz", &probe_dsa_dz);
  t1->Branch("probe_dsa_trkChi2", &probe_dsa_trkChi2);
  t1->Branch("probe_dsa_pt", &probe_dsa_pt);
  t1->Branch("probe_dsa_eta", &probe_dsa_eta);
  t1->Branch("probe_dsa_phi", &probe_dsa_phi);
  t1->Branch("probe_dsa_outerEta", &probe_dsa_outerEta);
  t1->Branch("probe_dsa_outerPhi", &probe_dsa_outerPhi);
  t1->Branch("probe_dsa_minDR", &probe_dsa_minDR);
  t1->Branch("probe_dsa_charge", &probe_dsa_charge);

  t1->Branch("tag_dsa_segmentMatches", &tag_dsa_segmentMatches);
  t1->Branch("tag_dsa_nsegments", &tag_dsa_nsegments);
  t1->Branch("tag_dsa_muonStations", &tag_dsa_muonStations);
  t1->Branch("tag_dsa_muonHits", &tag_dsa_muonHits);
  t1->Branch("tag_dsa_DTHits", &tag_dsa_DTHits);
  t1->Branch("tag_dsa_CSCHits", &tag_dsa_CSCHits);
  t1->Branch("tag_dsa_pterr", &tag_dsa_pterr);
  t1->Branch("tag_dsa_dxy", &tag_dsa_dxy);
  t1->Branch("tag_dsa_dz", &tag_dsa_dz);
  t1->Branch("tag_dsa_trkChi2", &tag_dsa_trkChi2);
  t1->Branch("tag_dsa_pt", &tag_dsa_pt);
  t1->Branch("tag_dsa_eta", &tag_dsa_eta);
  t1->Branch("tag_dsa_phi", &tag_dsa_phi);
  t1->Branch("tag_dsa_outerEta", &tag_dsa_outerEta);
  t1->Branch("tag_dsa_outerPhi", &tag_dsa_outerPhi);
  t1->Branch("tag_dsa_minDR", &tag_dsa_minDR);
  t1->Branch("tag_dsa_charge", &tag_dsa_charge);

  // selectors for probe
  for (unsigned int isel = 0; isel < selectorNames.size(); ++isel) {
    t1->Branch(TString("probe_" + selectorNames[isel]), &probe_selectors[isel]);
  }

  t1->Branch("probe_dgl_segmentMatches", &probe_dgl_segmentMatches);
  t1->Branch("probe_dgl_nsegments", &probe_dgl_nsegments);
  t1->Branch("probe_dgl_muonStations", &probe_dgl_muonStations);
  t1->Branch("probe_dgl_muonHits", &probe_dgl_muonHits);
  t1->Branch("probe_dgl_totalHits", &probe_dgl_totalHits);
  t1->Branch("probe_dgl_outerTrackerHits", &probe_dgl_outerTrackerHits);
  t1->Branch("probe_dgl_trackerHits", &probe_dgl_trackerHits);
  t1->Branch("probe_dgl_DTHits", &probe_dgl_DTHits);
  t1->Branch("probe_dgl_CSCHits", &probe_dgl_CSCHits);
  t1->Branch("probe_dgl_pterr", &probe_dgl_pterr);
  t1->Branch("probe_dgl_dxy", &probe_dgl_dxy);
  t1->Branch("probe_dgl_dz", &probe_dgl_dz);
  t1->Branch("probe_dgl_trkChi2", &probe_dgl_trkChi2);
  t1->Branch("probe_dgl_pt", &probe_dgl_pt);
  t1->Branch("probe_dgl_eta", &probe_dgl_eta);
  t1->Branch("probe_dgl_phi", &probe_dgl_phi);
  t1->Branch("probe_dgl_charge", &probe_dgl_charge);
  t1->Branch("probe_dgl_minDR", &probe_dgl_minDR);

  // Pair specific
  t1->Branch("pair_pt", &pair_pt);
  t1->Branch("pair_eta", &pair_eta);
  t1->Branch("pair_phi", &pair_phi);
  t1->Branch("pair_mass", &pair_mass);
  t1->Branch("pair_fit_mass", &pair_fit_mass);
  t1->Branch("pair_svprob", &pair_svprob);
  t1->Branch("pair_normalchi2", &pair_normalchi2);
  t1->Branch("pair_dz", &pair_dz);
  t1->Branch("pair_dR", &pair_dR);
  t1->Branch("pair_drM1", &pair_drM1);
  t1->Branch("pair_rank_vtx_prob", &pair_rank_vtx_prob);
  t1->Branch("pair_rank_dz_PV_SV", &pair_rank_dz_PV_SV);
  t1->Branch("pair_rank_dPhi_muons", &pair_rank_dPhi_muons);
  t1->Branch("pair_rank_dM_Z_Mmumu", &pair_rank_dM_Z_Mmumu);

  t1->Branch("pair_tuneP_pt", &pair_tuneP_pt);
  t1->Branch("pair_tuneP_eta", &pair_tuneP_eta);
  t1->Branch("pair_tuneP_phi", &pair_tuneP_phi);
  t1->Branch("pair_tuneP_mass", &pair_tuneP_mass);
  t1->Branch("pair_tuneP_fit_mass", &pair_tuneP_fit_mass);
  t1->Branch("pair_tuneP_svprob", &pair_tuneP_svprob);
  t1->Branch("pair_tuneP_normalchi2", &pair_tuneP_normalchi2);
  t1->Branch("pair_tuneP_dz", &pair_tuneP_dz);
  t1->Branch("pair_tuneP_dR", &pair_tuneP_dR);

  t1->Branch("tag_simType", &tag_simType);
  t1->Branch("tag_simExtType", &tag_simExtType);
  t1->Branch("tag_simFlavour", &tag_simFlavour);
  t1->Branch("tag_simHeaviestMotherFlavour", &tag_simHeaviestMotherFlavour);
  t1->Branch("tag_simPdgId", &tag_simPdgId);
  t1->Branch("tag_simMotherPdgId", &tag_simMotherPdgId);
  t1->Branch("tag_simBX", &tag_simBX);
  t1->Branch("tag_simProdRho", &tag_simProdRho);
  t1->Branch("tag_simProdZ", &tag_simProdZ);
  t1->Branch("tag_simPt", &tag_simPt);
  t1->Branch("tag_simEta", &tag_simEta);
  t1->Branch("tag_simPhi", &tag_simPhi);

  t1->Branch("probe_simType", &probe_simType);
  t1->Branch("probe_simExtType", &probe_simExtType);
  t1->Branch("probe_simFlavour", &probe_simFlavour);
  t1->Branch("probe_simHeaviestMotherFlavour", &probe_simHeaviestMotherFlavour);
  t1->Branch("probe_simPdgId", &probe_simPdgId);
  t1->Branch("probe_simMotherPdgId", &probe_simMotherPdgId);
  t1->Branch("probe_simBX", &probe_simBX);
  t1->Branch("probe_simProdRho", &probe_simProdRho);
  t1->Branch("probe_simProdZ", &probe_simProdZ);
  t1->Branch("probe_simPt", &probe_simPt);
  t1->Branch("probe_simEta", &probe_simEta);
  t1->Branch("probe_simPhi", &probe_simPhi);
}

void NtupleContent::CreateExtraTrgBranches(const std::vector<std::string> &HLTs, bool isTag = false) {
  for (unsigned int ihlt = 0; ihlt < HLTs.size(); ihlt++) {
    if (isTag) {
      t1->Branch(TString("tag_" + HLTs[ihlt]), &tag_trg[ihlt]);
      t1->Branch(TString("tag_" + HLTs[ihlt] + "_pt"), &tag_trg_pt[ihlt]);
      t1->Branch(TString("tag_" + HLTs[ihlt] + "_eta"), &tag_trg_eta[ihlt]);
      t1->Branch(TString("tag_" + HLTs[ihlt] + "_phi"), &tag_trg_phi[ihlt]);
      t1->Branch(TString("tag_" + HLTs[ihlt] + "_dr"), &tag_trg_dr[ihlt]);
    } else {
      t1->Branch(TString("probe_" + HLTs[ihlt]), &probe_trg[ihlt]);
      t1->Branch(TString("probe_" + HLTs[ihlt] + "_pt"), &probe_trg_pt[ihlt]);
      t1->Branch(TString("probe_" + HLTs[ihlt] + "_eta"), &probe_trg_eta[ihlt]);
      t1->Branch(TString("probe_" + HLTs[ihlt] + "_phi"), &probe_trg_phi[ihlt]);
      t1->Branch(TString("probe_" + HLTs[ihlt] + "_dr"), &probe_trg_dr[ihlt]);
    }
  }
}

void NtupleContent::ClearBranches() {
  run = -1;
  event = -1;
  ls = -1;
  genWeight = -99;
  BSpot_x = -99;
  BSpot_y = -99;
  BSpot_z = -99;
  pv_x = -99;
  pv_y = -99;
  pv_z = -99;
  nvertices = 0;
  trueNumInteractions = -1.0;
  puNumInteractions = -1;
  Rho = -1;
  nmuons = 0;
  ntag = 0;
  npairs = 0;
  pair_probeMultiplicity = 0;

  for (unsigned int itrg = 0; itrg < NTRIGGERMAX; itrg++) {
    trigger[itrg] = false;
    tag_trg[itrg] = false;
    tag_trg_pt[itrg] = -99;
    tag_trg_eta[itrg] = -99;
    tag_trg_phi[itrg] = -99;
    tag_trg_dr[itrg] = 99;
    probe_trg[itrg] = false;
    probe_trg_pt[itrg] = -99;
    probe_trg_eta[itrg] = -99;
    probe_trg_phi[itrg] = -99;
    probe_trg_dr[itrg] = 99;
  }

  for (unsigned int isel = 0; isel < 100; isel++) {
    probe_selectors[isel] = false;
  }

  // Gens
  genmu1_pt = 0;
  genmu1_eta = -99;
  genmu1_phi = -99;
  genmu1_charge = 0;
  genmu2_pt = 0;
  genmu2_eta = -99;
  genmu2_phi = -99;
  genmu2_charge = 0;
  genMass = -99;

  genmuFSfromHP1_pt = -99;
  genmuFSfromHP1_eta = -99;
  genmuFSfromHP1_phi = -99;
  genmuFSfromHP1_charge = -99;
  genmuFSfromHP2_pt = -99;
  genmuFSfromHP2_eta = -99;
  genmuFSfromHP2_phi = -99;
  genmuFSfromHP2_charge = -99;
  genMassFSfromHP = -99;

  trg_filter.clear();
  trg_pt.clear();
  trg_eta.clear();
  trg_phi.clear();
  prb_filter.clear();
  prb_pt.clear();
  prb_eta.clear();
  prb_phi.clear();

  genJets_pt.clear();
  genJets_eta.clear();
  genJets_phi.clear();
  genJets_mass.clear();

  jets_pt.clear();
  jets_eta.clear();
  jets_phi.clear();
  jets_mass.clear();
  jets_isTight.clear();
  jets_isTightLepVeto.clear();
  jets_bTag_deepCSV.clear();
  jets_bTag_deepFlav.clear();
  nTightJets = 0;
  nTightLepVetoJets = 0;

  tag_pt = 0;
  tag_eta = -99;
  tag_phi = -99;
  tag_charge = -99;
  tag_pterr = 0;
  tag_dxy = -99;
  tag_dz = -99;
  tag_isPF = false;
  tag_isSA = false;
  tag_isdSA = false;
  tag_isTracker = false;
  tag_isGlobal = false;
  tag_isLoose = false;
  tag_isMedium = false;
  tag_isTight = false;
  tag_isSoft = false;
  tag_isHighPt = false;
  tag_relIso04 = -99;
  tag_miniIso = -1.;
  tag_miniIsoCharged = 0.;
  tag_miniIsoPhotons = 0.;
  tag_miniIsoNeutrals = 0.;
  tag_isMatchedGen = false;
  tag_minDR = 0.;
  tag_ptRel_minDR = 0.;
  tag_iso03_sumPt = -99;
  tag_pfIso03_charged = -99;
  tag_pfIso03_neutral = -99;
  tag_pfIso03_photon = -99;
  tag_pfIso03_sumPU = -99;
  tag_pfIso04_charged = -99;
  tag_pfIso04_neutral = -99;
  tag_pfIso04_photon = -99;
  tag_pfIso04_sumPU = -99;
  tag_tuneP_pt = -99;
  tag_tuneP_pterr = -99;
  tag_nsegments = -99;

  iprobe = 0;
  probe_pt = 0;
  probe_eta = -99;
  probe_phi = -99;
  probe_charge = -99;
  probe_isLoose = false;
  probe_isMedium = false;
  probe_isTight = false;
  probe_isSoft = false;
  probe_isHighPt = false;
  probe_isArbitratedTracker = false;
  probe_isMuMatched = false;
  probe_isPF = false;
  probe_isSA = false;
  probe_isTracker = false;
  probe_isGlobal = false;
  probe_isdSA = false;
  probe_isdGlobal = false;
  probe_isCosmic = false;
  probe_ncosmic = -99;
  probe_cosmic_minDR = +99;
  probe_isGood = false;
  probe_isHighPurity = false;
  probe_validFraction = -99;
  probe_trkChi2 = -99;
  probe_positionChi2 = -99;
  probe_trkKink = -99;
  probe_segmentCompatibility = -99;
  probe_trackerLayers = -99;
  probe_pixelLayers = -99;
  probe_muonStations = -99;
  probe_muonHits = -99;
  probe_DTHits = -99;
  probe_CSCHits = -99;
  probe_pterr = 0;
  probe_dxy = -99;
  probe_dz = -99;
  probe_relIso04 = -99;
  probe_miniIso = -1.;
  probe_miniIsoCharged = 0.;
  probe_miniIsoPhotons = 0.;
  probe_miniIsoNeutrals = 0.;
  probe_isMatchedGen = false;
  probe_minDR = 0.;
  probe_ptRel_minDR = 0.;
  probe_iso03_sumPt = -99;
  probe_pfIso03_charged = -99;
  probe_pfIso03_neutral = -99;
  probe_pfIso03_photon = -99;
  probe_pfIso03_sumPU = -99;
  probe_pfIso04_charged = -99;
  probe_pfIso04_neutral = -99;
  probe_pfIso04_photon = -99;
  probe_pfIso04_sumPU = -99;
  probe_pixelHits = -99;
  probe_matchedStations = -99;
  probe_expectedMatchedStations = -99;
  probe_RPCLayers = -99;
  probe_stationMask = 0;
  probe_nShowers = -99;
  probe_tuneP_pt = -99;
  probe_tuneP_pterr = -99;
  probe_tuneP_muonHits = -99;
  probe_nsegments = -99;

  l1pt = -99;
  l1q = -99;
  l1dr = 99;
  l1ptByQ = -99;
  l1qByQ = -99;
  l1drByQ = 99;

  tag_l1pt = -99;
  tag_l1q = -99;
  tag_l1dr = 99;
  tag_l1ptByQ = -99;
  tag_l1qByQ = -99;
  tag_l1drByQ = 99;

  probe_dsa_segmentMatches = -99;
  probe_dsa_nsegments = -99;
  probe_dsa_muonStations = -99;
  probe_dsa_muonHits = -99;
  probe_dsa_DTHits = -99;
  probe_dsa_CSCHits = -99;
  probe_dsa_pterr = 0;
  probe_dsa_dxy = -99;
  probe_dsa_dz = -99;
  probe_dsa_trkChi2 = -99;
  probe_dsa_pt = 0;
  probe_dsa_eta = -99;
  probe_dsa_phi = -99;
  probe_dsa_outerEta = -99;
  probe_dsa_outerPhi = -99;
  probe_dsa_minDR = +99;
  probe_dsa_charge = -99;

  tag_dsa_segmentMatches = -99;
  tag_dsa_nsegments = -99;
  tag_dsa_muonStations = -99;
  tag_dsa_muonHits = -99;
  tag_dsa_DTHits = -99;
  tag_dsa_CSCHits = -99;
  tag_dsa_pterr = 0;
  tag_dsa_dxy = -99;
  tag_dsa_dz = -99;
  tag_dsa_trkChi2 = -99;
  tag_dsa_pt = 0;
  tag_dsa_eta = -99;
  tag_dsa_phi = -99;
  tag_dsa_outerEta = -99;
  tag_dsa_outerPhi = -99;
  tag_dsa_minDR = +99;
  tag_dsa_charge = -99;

  probe_dgl_segmentMatches = -99;
  probe_dgl_nsegments = -99;
  probe_dgl_muonStations = -99;
  probe_dgl_muonHits = -99;
  probe_dgl_outerTrackerHits = -99;
  probe_dgl_trackerHits = -99;
  probe_dgl_totalHits = -99;
  probe_dgl_DTHits = -99;
  probe_dgl_CSCHits = -99;
  probe_dgl_pterr = 0;
  probe_dgl_dxy = -99;
  probe_dgl_dz = -99;
  probe_dgl_trkChi2 = -99;
  probe_dgl_pt = 0;
  probe_dgl_eta = -99;
  probe_dgl_phi = -99;
  probe_dgl_charge = -99;
  probe_dgl_minDR = +99;

  pair_pt = 0;
  pair_mass = 0;
  pair_eta = -99;
  pair_phi = -99;
  pair_fit_mass = 0;
  pair_svprob = 0;
  pair_normalchi2 = 0;
  pair_dz = -99;
  pair_dR = -99;
  pair_drM1 = -99;
  pair_rank_vtx_prob = -1;
  pair_rank_dz_PV_SV = -1;
  pair_rank_dPhi_muons = -1;
  pair_rank_dM_Z_Mmumu = -1;

  pair_tuneP_pt = -99;
  pair_tuneP_mass = -99;
  pair_tuneP_eta = -99;
  pair_tuneP_phi = -99;
  pair_tuneP_fit_mass = -99;
  pair_tuneP_svprob = -99;
  pair_tuneP_normalchi2 = -99;
  pair_tuneP_dz = -99;
  pair_tuneP_dR = -99;

  tag_simType = -99;
  tag_simExtType = -99;
  tag_simFlavour = -99;
  tag_simHeaviestMotherFlavour = -99;
  tag_simPdgId = -99;
  tag_simMotherPdgId = -99;
  tag_simBX = -99;
  tag_simProdRho = -99;
  tag_simProdZ = -99;
  tag_simPt = -99;
  tag_simEta = -99;
  tag_simPhi = -99;

  probe_simType = -99;
  probe_simExtType = -99;
  probe_simFlavour = -99;
  probe_simHeaviestMotherFlavour = -99;
  probe_simPdgId = -99;
  probe_simMotherPdgId = -99;
  probe_simBX = -99;
  probe_simProdRho = -99;
  probe_simProdZ = -99;
  probe_simPt = -99;
  probe_simEta = -99;
  probe_simPhi = -99;
}
