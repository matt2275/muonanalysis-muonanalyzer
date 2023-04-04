#include "Analysis_NtupleContent.h"

Analysis_NtupleContent::Analysis_NtupleContent() {}

Analysis_NtupleContent::~Analysis_NtupleContent() {}

  void Analysis_NtupleContent::SetTreeVariables(bool keepMuons_ ,bool keepElectrons_ ,bool keepTracks_ ,bool keepPFcands_ 
  ,bool keepPhotons_ ,bool keepCaloTowers_ ,bool keepZDC_ ){
    keepMuons = keepMuons_ ;
    keepElectrons = keepElectrons_ ;
    keepTracks = keepTracks_ ;
    keepPFcands = keepPFcands_;
    keepPhotons = keepPhotons_;
    keepCaloTowers = keepCaloTowers_;
    keepZDC = keepZDC_;
  }

void Analysis_NtupleContent::SetTree(TTree *mytree) { t1 = mytree; }

void Analysis_NtupleContent::CreateBranches(const std::vector<std::string> &HLTs,
                                   const std::vector<std::string> &selectorNames) {
  // General
  t1->Branch("run", &run);
  t1->Branch("event", &event);
  t1->Branch("ls", &ls);
  t1->Branch("fromFullAOD", &fromFullAOD);
  t1->Branch("genWeight", &genWeight);
  t1->Branch("BSpot_x", &BSpot_x);
  t1->Branch("BSpot_y", &BSpot_y);
  t1->Branch("BSpot_z", &BSpot_z);
  t1->Branch("pv_x", &pv_x);
  t1->Branch("pv_y", &pv_y);
  t1->Branch("pv_z", &pv_z);
  t1->Branch("nVertices", &nvertices);
  t1->Branch("hasValidVertex", &hasValidVertex);
  t1->Branch("hasFakeVertex", &hasFakeVertex);
  t1->Branch("nTrueInteractions", &trueNumInteractions);
  t1->Branch("nPUInteractions", &puNumInteractions);
  t1->Branch("rho", &Rho);
  t1->Branch("ntag", &ntag);
  t1->Branch("ntag_electron", &ntag_electron);
  t1->Branch("npairs", &npairs);
  t1->Branch("npairs_3prong", &npairs_3prong);
  t1->Branch("tag_isMuon", &tag_isMuon);
  t1->Branch("tag_isElectron", &tag_isElectron);
  t1->Branch("probe_isMuon", &probe_isMuon );
  t1->Branch("probe_isElectron", &probe_isElectron );
  t1->Branch("probe_isPion", &probe_isPion );
  t1->Branch("probe_isOther", &probe_isOther );
  t1->Branch("probe_typeSum", &probe_typeSum );
  
  t1->Branch("gentrk1_pt", &gentrk1_pt);
  t1->Branch("gentrk1_eta", &gentrk1_eta);
  t1->Branch("gentrk1_phi", &gentrk1_phi);
  t1->Branch("gentrk1_charge", &gentrk1_charge);
  t1->Branch("gentrk1_M", &gentrk1_M);
  t1->Branch("gentrk1_vtx_x", &gentrk1_vtx_x);
  t1->Branch("gentrk1_vtx_y", &gentrk1_vtx_y);
  t1->Branch("gentrk1_vtx_z", &gentrk1_vtx_z);
  t1->Branch("gentrk1_pdgId", &gentrk1_pdgId);  
  t1->Branch("gentrk2_pt", &gentrk2_pt);
  t1->Branch("gentrk2_eta", &gentrk2_eta);
  t1->Branch("gentrk2_phi", &gentrk2_phi);
  t1->Branch("gentrk2_charge", &gentrk2_charge);
  t1->Branch("gentrk2_M", &gentrk2_M);
  t1->Branch("gentrk2_vtx_x", &gentrk2_vtx_x);
  t1->Branch("gentrk2_vtx_y", &gentrk2_vtx_y);
  t1->Branch("gentrk2_vtx_z", &gentrk2_vtx_z);
  t1->Branch("gentrk2_pdgId", &gentrk2_pdgId); 
  t1->Branch("genFinalMass", &genFinalMass);
  t1->Branch("genFinalPt", &genFinalPt);
  t1->Branch("genFinalEta", &genFinalEta);
  t1->Branch("genFinalPhi", &genFinalPhi);
  t1->Branch("genFinalE", &genFinalE);
  t1->Branch("genFinalPx", &genFinalPx);
  t1->Branch("genFinalPy", &genFinalPy);
  t1->Branch("genFinalPz", &genFinalPz);  
  

  t1->Branch("gentau1_pt", &gentau1_pt);
  t1->Branch("gentau1_eta", &gentau1_eta);
  t1->Branch("gentau1_phi", &gentau1_phi);
  t1->Branch("gentau1_charge", &gentau1_charge);
  t1->Branch("gentau1_vtx_x", &gentau1_vtx_x);
  t1->Branch("gentau1_vtx_y", &gentau1_vtx_y);
  t1->Branch("gentau1_vtx_z", &gentau1_vtx_z);
  t1->Branch("gentau1_pdgId", &gentau1_pdgId);  
  t1->Branch("gentau1_decay_el", &gentau1_decay_el);
  t1->Branch("gentau1_decay_mu", &gentau1_decay_mu); 
  t1->Branch("gentau1_decay_1prong", &gentau1_decay_1prong); 
  t1->Branch("gentau1_decay_3prong", &gentau1_decay_3prong); 
  t1->Branch("gentau1_decay_other", &gentau1_decay_other); 
  t1->Branch("gentau1_decay", &gentau1_decay); 
  t1->Branch("gentau1_N_pi0", &gentau1_N_pi0);   
  t1->Branch("gentau2_pt", &gentau2_pt);
  t1->Branch("gentau2_eta", &gentau2_eta);
  t1->Branch("gentau2_phi", &gentau2_phi);
  t1->Branch("gentau2_charge", &gentau2_charge);
  t1->Branch("gentau2_vtx_x", &gentau2_vtx_x);
  t1->Branch("gentau2_vtx_y", &gentau2_vtx_y);
  t1->Branch("gentau2_vtx_z", &gentau2_vtx_z);
  t1->Branch("gentau2_pdgId", &gentau2_pdgId);
  t1->Branch("gentau2_decay_el", &gentau2_decay_el);
  t1->Branch("gentau2_decay_mu", &gentau2_decay_mu); 
  t1->Branch("gentau2_decay_1prong", &gentau2_decay_1prong); 
  t1->Branch("gentau2_decay_3prong", &gentau2_decay_3prong); 
  t1->Branch("gentau2_decay_other", &gentau2_decay_other); 
  t1->Branch("gentau2_decay", &gentau2_decay);  
  t1->Branch("gentau2_N_pi0", &gentau2_N_pi0);   
  t1->Branch("genDiTauMass", &genDiTauMass);
  t1->Branch("genDiTauPt", &genDiTauPt);
  t1->Branch("genDiTauEta", &genDiTauEta);
  t1->Branch("genDiTauPhi", &genDiTauPhi);
  t1->Branch("genDiTauE", &genDiTauE);
  t1->Branch("genDiTauPx", &genDiTauPx);
  t1->Branch("genDiTauPy", &genDiTauPy);
  t1->Branch("genDiTauPz", &genDiTauPz);
  
  t1->Branch("nGen", &nGen );
  t1->Branch("gen_pdgId", &gen_pdgId );
  t1->Branch("gen_pT", &gen_pT );
  t1->Branch("gen_eta", &gen_eta );
  t1->Branch("gen_phi", &gen_phi ); 
  t1->Branch("gen_charge", &gen_charge );
  t1->Branch("gen_mass", &gen_mass );  
  t1->Branch("gen_vtx_x", &gen_vtx_x );
  t1->Branch("gen_vtx_y", &gen_vtx_y );
  t1->Branch("gen_vtx_z", &gen_vtx_z );
  t1->Branch("gen_status", &gen_status );
  t1->Branch("gen_E", &gen_E );
  t1->Branch("gen_Et", &gen_Et );  
  
  
  t1->Branch("nGen_tau1", &nGen_tau1 );
  t1->Branch("gen_pdgId_tau1", &gen_pdgId_tau1 );
  t1->Branch("gen_pT_tau1", &gen_pT_tau1 );
  t1->Branch("gen_eta_tau1", &gen_eta_tau1 );
  t1->Branch("gen_phi_tau1", &gen_phi_tau1 ); 
  t1->Branch("gen_charge_tau1", &gen_charge_tau1 );
  t1->Branch("gen_mass_tau1", &gen_mass_tau1 );  
  t1->Branch("gen_vtx_x_tau1", &gen_vtx_x_tau1 );
  t1->Branch("gen_vtx_y_tau1", &gen_vtx_y_tau1 );
  t1->Branch("gen_vtx_z_tau1", &gen_vtx_z_tau1 );
  t1->Branch("gen_status_tau1", &gen_status_tau1 );
  t1->Branch("gen_E_tau1", &gen_E_tau1 );
  t1->Branch("gen_Et_tau1", &gen_Et_tau1 );

  t1->Branch("nGen_tau2", &nGen_tau2 );
  t1->Branch("gen_pdgId_tau2", &gen_pdgId_tau2 );
  t1->Branch("gen_pT_tau2", &gen_pT_tau2 );
  t1->Branch("gen_eta_tau2", &gen_eta_tau2 );
  t1->Branch("gen_phi_tau2", &gen_phi_tau2 ); 
  t1->Branch("gen_charge_tau2", &gen_charge_tau2 );
  t1->Branch("gen_mass_tau2", &gen_mass_tau2 );  
  t1->Branch("gen_vtx_x_tau2", &gen_vtx_x_tau2 );
  t1->Branch("gen_vtx_y_tau2", &gen_vtx_y_tau2 );
  t1->Branch("gen_vtx_z_tau2", &gen_vtx_z_tau2 );
  t1->Branch("gen_status_tau2", &gen_status_tau2 );
  t1->Branch("gen_E_tau2", &gen_E_tau2 );
  t1->Branch("gen_Et_tau2", &gen_Et_tau2 );

  t1->Branch("gentrk1_match_dr", &gentrk1_match_dr );
  t1->Branch("gentrk1_match_dphi", &gentrk1_match_dphi );
  t1->Branch("gentrk1_match_diff_eta", &gentrk1_match_diff_eta );
  t1->Branch("gentrk1_match_diff_pt", &gentrk1_match_diff_pt ); 
  t1->Branch("gentrk1_diff_vtx_x", &gentrk1_diff_vtx_x );
  t1->Branch("gentrk1_diff_vtx_y", &gentrk1_diff_vtx_y );
  t1->Branch("gentrk1_diff_vtx_z", & gentrk1_diff_vtx_z);
  t1->Branch("gentrk1_isTag", &gentrk1_isTag );
  
  t1->Branch("gentrk2_match_dr", &gentrk2_match_dr );
  t1->Branch("gentrk2_match_dphi", &gentrk2_match_dphi );
  t1->Branch("gentrk2_match_diff_eta", &gentrk2_match_diff_eta );
  t1->Branch("gentrk2_match_diff_pt", &gentrk2_match_diff_pt ); 
  t1->Branch("gentrk2_diff_vtx_x", &gentrk2_diff_vtx_x );
  t1->Branch("gentrk2_diff_vtx_y", &gentrk2_diff_vtx_y );
  t1->Branch("gentrk2_diff_vtx_z", &gentrk2_diff_vtx_z );  
  t1->Branch("gentrk2_isTag", &gentrk2_isTag );
  
  
  for (unsigned int ihlt = 0; ihlt < HLTs.size(); ihlt++)
    t1->Branch(TString(HLTs[ihlt]), &trigger[ihlt]);

  // Tag specific
  t1->Branch("tag_pt", &tag_pt);
  t1->Branch("tag_eta", &tag_eta);
  t1->Branch("tag_phi", &tag_phi);
  t1->Branch("tag_charge", &tag_charge);
  t1->Branch("tag_pterr", &tag_pterr);
  t1->Branch("tag_dxy", &tag_dxy);
  t1->Branch("tag_dz", &tag_dz);
  t1->Branch("tag_vtx_x", &tag_vtx_x);
  t1->Branch("tag_vtx_y", &tag_vtx_y);
  t1->Branch("tag_vtx_z", &tag_vtx_z);
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
  t1->Branch("tag_pfIso04_charged", &tag_pfIso04_charged);
  t1->Branch("tag_pfIso04_neutral", &tag_pfIso04_neutral);
  t1->Branch("tag_pfIso04_photon", &tag_pfIso04_photon);
  t1->Branch("tag_pfIso04_sumPU", &tag_pfIso04_sumPU);
  t1->Branch("tag_tuneP_pt", &tag_tuneP_pt);
  t1->Branch("tag_tuneP_pterr", &tag_tuneP_pterr);
  t1->Branch("tag_nsegments", &tag_nsegments);
  t1->Branch("tag_hasTrackMatch", &tag_hasTrackMatch);
  t1->Branch("tag_TrackMatchDR", &tag_TrackMatchDR);
  
  // Probe specific
  t1->Branch("iprobe", &iprobe);
  t1->Branch("probe_pt", &probe_pt);
  t1->Branch("probe_eta", &probe_eta);
  t1->Branch("probe_phi", &probe_phi);
  t1->Branch("probe_charge", &probe_charge);
  t1->Branch("probe_vtx_x", &probe_vtx_x);
  t1->Branch("probe_vtx_y", &probe_vtx_y);
  t1->Branch("probe_vtx_z", &probe_vtx_z);
  // t1->Branch("probe_pterr", &probe_pterr);
  // t1->Branch("probe_dxy", &probe_dxy);
  // t1->Branch("probe_dz", &probe_dz);
  // t1->Branch("probe_isPF", &probe_isPF);
  // t1->Branch("probe_isSA", &probe_isSA);
  // t1->Branch("probe_isTracker", &probe_isTracker);
  // t1->Branch("probe_isGlobal", &probe_isGlobal);
  // t1->Branch("probe_isLoose", &probe_isLoose);
  // t1->Branch("probe_isMedium", &probe_isMedium);
  // t1->Branch("probe_isTight", &probe_isTight);
  // t1->Branch("probe_isSoft", &probe_isSoft);
  // t1->Branch("probe_isHighPt", &probe_isHighPt);
  // t1->Branch("probe_isArbitratedTracker", &probe_isArbitratedTracker);
  // t1->Branch("probe_isMuMatched", &probe_isMuMatched);
  // t1->Branch("probe_isdSA", &probe_isdSA);
  // t1->Branch("probe_isdGlobal", &probe_isdGlobal);
  t1->Branch("probe_isCosmic", &probe_isCosmic);
  t1->Branch("probe_ncosmic", &probe_ncosmic);
  t1->Branch("probe_cosmic_minDR", &probe_cosmic_minDR);
  // //  t1->Branch("probe_isGood", &probe_isGood);
  t1->Branch("probe_isHighPurity", &probe_isHighPurity);
  // t1->Branch("probe_validFraction", &probe_validFraction);
  // t1->Branch("probe_trkChi2", &probe_trkChi2);
  // t1->Branch("probe_positionChi2", &probe_positionChi2);
  // t1->Branch("probe_trkKink", &probe_trkKink);
  // // t1->Branch("probe_segmentCompatibility", &probe_segmentCompatibility);
  // t1->Branch("probe_trackerLayers", &probe_trackerLayers);
  // t1->Branch("probe_pixelLayers", &probe_pixelLayers);
  // t1->Branch("probe_muonStations", &probe_muonStations);
  // t1->Branch("probe_muonHits", &probe_muonHits);
  // t1->Branch("probe_DTHits", &probe_DTHits);
  // t1->Branch("probe_CSCHits", &probe_CSCHits);
  t1->Branch("probe_relIso04", &probe_relIso04);
  t1->Branch("probe_miniIso", &probe_miniIso);
  t1->Branch("probe_miniIsoCharged", &probe_miniIsoCharged);
  t1->Branch("probe_miniIsoPhotons", &probe_miniIsoPhotons);
  t1->Branch("probe_miniIsoNeutrals", &probe_miniIsoPhotons);
  t1->Branch("probe_isMatchedGen", &probe_isMatchedGen);
  // t1->Branch("probe_minDR", &probe_minDR);
  // t1->Branch("probe_ptRel_minDR", &probe_ptRel_minDR);
  // t1->Branch("probe_iso03_sumPt", &probe_iso03_sumPt);
  // t1->Branch("probe_pfIso04_charged", &probe_pfIso04_charged);
  // t1->Branch("probe_pfIso04_neutral", &probe_pfIso04_neutral);
  // t1->Branch("probe_pfIso04_photon", &probe_pfIso04_photon);
  // t1->Branch("probe_pfIso04_sumPU", &probe_pfIso04_sumPU);
  // t1->Branch("probe_pixelHits", &probe_pixelHits);
  // t1->Branch("probe_matchedStations", &probe_matchedStations);
  // t1->Branch("probe_expectedMatchedStations", &probe_expectedMatchedStations);
  // t1->Branch("probe_RPCLayers", &probe_RPCLayers);
  // t1->Branch("probe_stationMask", &probe_stationMask);
  // t1->Branch("probe_nShowers", &probe_nShowers);
  // t1->Branch("probe_tuneP_pt", &probe_tuneP_pt);
  // t1->Branch("probe_tuneP_pterr", &probe_tuneP_pterr);
  // t1->Branch("probe_tuneP_muonHits", &probe_tuneP_muonHits);
  // t1->Branch("probe_nsegments", &probe_nsegments);
  t1->Branch("probe_hasMuonMatch", &probe_hasMuonMatch );
  t1->Branch("probe_MuonMatchDR", &probe_MuonMatchDR );
  t1->Branch("probe_hasElectronMatch", &probe_hasElectronMatch );
  t1->Branch("probe_ElectronMatchDR", &probe_ElectronMatchDR );
  t1->Branch("probe_hasPFCMatch", &probe_hasPFCMatch );
  t1->Branch("probe_PFCMatchDR", &probe_PFCMatchDR );
  t1->Branch("probe_pfcID", &probe_pfcID );

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


  // selectors for probe
  for (unsigned int isel = 0; isel < selectorNames.size(); ++isel) {
    t1->Branch(TString("probe_" + selectorNames[isel]), &probe_selectors[isel]);
  }


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
  t1->Branch("pair_dphi", &pair_dphi);
  t1->Branch("pair_rank_vtx_prob", &pair_rank_vtx_prob);
  t1->Branch("pair_rank_dPhi_muons", &pair_rank_dPhi_muons);
  t1->Branch("pair_rank_Mass_Mmumu", &pair_rank_Mass_Mmumu);

  
  t1->Branch("refit_vtx_x", &refit_vtx_x);
  t1->Branch("refit_vtx_y", &refit_vtx_y);
  t1->Branch("refit_vtx_z", &refit_vtx_z);
  t1->Branch("tag_refit_transverse_IP", &tag_transverse_IP);
  t1->Branch("probe_refit_transverse_IP", &probe_transverse_IP);
  t1->Branch("tag_transverse_IP", &tag_transverse_IP);
  t1->Branch("probe_transverse_IP", &probe_transverse_IP);


  t1->Branch("pair_tuneP_pt", &pair_tuneP_pt);
  t1->Branch("pair_tuneP_eta", &pair_tuneP_eta);
  t1->Branch("pair_tuneP_phi", &pair_tuneP_phi);
  t1->Branch("pair_tuneP_mass", &pair_tuneP_mass);
  t1->Branch("pair_tuneP_fit_mass", &pair_tuneP_fit_mass);
  t1->Branch("pair_tuneP_svprob", &pair_tuneP_svprob);
  t1->Branch("pair_tuneP_normalchi2", &pair_tuneP_normalchi2);
  t1->Branch("pair_tuneP_dz", &pair_tuneP_dz);
  t1->Branch("pair_tuneP_dR", &pair_tuneP_dR);
  t1->Branch("pair_tuneP_dR", &pair_tuneP_dphi);
  
  t1->Branch("refit_tuneP_vtx_x", &refit_tuneP_vtx_x);
  t1->Branch("refit_tuneP_vtx_y", &refit_tuneP_vtx_y);
  t1->Branch("refit_tuneP_vtx_z", &refit_tuneP_vtx_z);
  t1->Branch("tag_refit_tuneP_transverse_IP", &tag_tuneP_transverse_IP);
  t1->Branch("probe_refit_tuneP_transverse_IP", &probe_tuneP_transverse_IP);
  t1->Branch("tag_tuneP_transverse_IP", &tag_tuneP_transverse_IP);
  t1->Branch("probe_tuneP_transverse_IP", &probe_tuneP_transverse_IP);
  
  //commented out to save space, added useful details for HIUPC collisions
    if(keepMuons){
    t1->Branch("nMu", &nMu);
    t1->Branch("mu_isProbe", &mu_isProbe);
    t1->Branch("mu_isTag", &mu_isTag);    
    t1->Branch("muPt", &muPt);
    t1->Branch("muEta", &muEta);
    t1->Branch("muPhi", &muPhi);
    t1->Branch("muCharge", &muCharge);
    t1->Branch("muType", &muType);
    t1->Branch("muIsGood", &muIsGood);
    
    t1->Branch("muIsGlobal", &muIsGlobal);
    t1->Branch("muIsTracker", &muIsTracker);
    t1->Branch("muIsPF", &muIsPF);
    t1->Branch("muIsSTA", &muIsSTA);

    t1->Branch("muD0", &muD0);
    t1->Branch("muDz", &muDz);
    t1->Branch("muIP3D", &muIP3D);
    t1->Branch("muD0Err", &muD0Err);
    t1->Branch("muDzErr", &muDzErr);
    t1->Branch("muIP3DErr", &muIP3DErr);
    t1->Branch("muChi2NDF", &muChi2NDF);
    t1->Branch("muInnerD0", &muInnerD0);
    t1->Branch("muInnerDz", &muInnerDz);
    
    t1->Branch("muInnerD0Err", &muInnerD0Err);
    t1->Branch("muInnerDzErr", &muInnerDzErr);
    t1->Branch("muInnerPt", &muInnerPt);
    t1->Branch("muInnerPtErr", &muInnerPtErr);
    t1->Branch("muInnerEta", &muInnerEta);

    t1->Branch("muTrkLayers", &muTrkLayers);
    t1->Branch("muPixelLayers", &muPixelLayers);
    t1->Branch("muPixelHits", &muPixelHits);
    t1->Branch("muMuonHits", &muMuonHits);
    t1->Branch("muTrkQuality", &muTrkQuality);
    t1->Branch("muStations", &muStations);
    t1->Branch("muIsoTrk", &muIsoTrk);
    t1->Branch("muPFChIso", &muPFChIso);
    t1->Branch("muPFPhoIso", &muPFPhoIso);
    t1->Branch("muPFNeuIso", &muPFNeuIso);
    t1->Branch("muPFPUIso", &muPFPUIso);
    t1->Branch("muIDSoft", &muIDSoft);
    t1->Branch("muIDLoose", &muIDLoose);
    t1->Branch("muIDMedium", &muIDMedium);
    t1->Branch("muIDMediumPrompt", &muIDMediumPrompt);
    t1->Branch("muIDTight", &muIDTight);
    t1->Branch("muIDGlobalHighPt", &muIDGlobalHighPt);
    t1->Branch("muIDTrkHighPt", &muIDTrkHighPt);
    t1->Branch("muIDInTime", &muIDInTime);

    }
    
    if(keepTracks){
    t1->Branch("nTrk", &nTrk);
    t1->Branch("trkisTag", &trkisTag);
    t1->Branch("trkisProbe", &trkisProbe);
    t1->Branch("trkPt", &trkPt);
    t1->Branch("trkP", &trkP);
    t1->Branch("trkEta", &trkEta);
    t1->Branch("trkPhi", &trkPhi);
    t1->Branch("trkcharge", &trkcharge);
    t1->Branch("trkvx", &trkvx);
    t1->Branch("trkvy", &trkvy);
    t1->Branch("trkvz", &trkvz);
    t1->Branch("trknormchi2", &trknormchi2);
    t1->Branch("trkchi2", &trkchi2);
    t1->Branch("trkd0", &trkd0);
    t1->Branch("trkdxy", &trkdxy);
    t1->Branch("trkdz", &trkdz);
    // t1->Branch("trkdxyBS", &trkdxyBS);
    // t1->Branch("trkdzBS", &trkdzBS);
    t1->Branch("trkdxyError", &trkdxyError);
    t1->Branch("trkdzError", &trkdzError);
    t1->Branch("trkValidHits", &trkValidHits);
    t1->Branch("trkMissHits", &trkMissHits);
    t1->Branch("trkPurity", &trkPurity);
    }


   if(keepElectrons){
   t1->Branch("nEle", &nEle);  
   t1->Branch("ele_isProbe", &ele_isProbe);
   t1->Branch("ele_isTag", &ele_isTag);
   t1->Branch("eleCharge", &eleCharge);
   t1->Branch("eleChargeConsistent", &eleChargeConsistent);
   t1->Branch("eleSCPixCharge", &eleSCPixCharge);
   // t1->Branch("eleCtfCharge", &eleCtfCharge);
   t1->Branch("eleEn", &eleEn);
   t1->Branch("eleD0", &eleD0);
   t1->Branch("eleDz", &eleDz);
   t1->Branch("eleIP3D", &eleIP3D);
   t1->Branch("eleD0Err", &eleD0Err);
   t1->Branch("eleDzErr", &eleDzErr);
   t1->Branch("eleIP3DErr", &eleIP3DErr);
   t1->Branch("eleTrkPt", &eleTrkPt);
   t1->Branch("eleTrkEta", &eleTrkEta);
   t1->Branch("eleTrkPhi", &eleTrkPhi);
   t1->Branch("eleTrkCharge", &eleTrkCharge);
   t1->Branch("eleTrkPtErr", &eleTrkPtErr);
   t1->Branch("eleTrkChi2", &eleTrkChi2);
   t1->Branch("eleTrkNdof", &eleTrkNdof);
   t1->Branch("eleTrkNormalizedChi2", &eleTrkNormalizedChi2);
   t1->Branch("eleTrkValidHits", &eleTrkValidHits);
   t1->Branch("eleTrkLayers", &eleTrkLayers);

   t1->Branch("elePt", &elePt);
   t1->Branch("eleEta", &eleEta);
   t1->Branch("elePhi", &elePhi);
   t1->Branch("eleSCEn", &eleSCEn);
   t1->Branch("eleESEn", &eleESEn);
   t1->Branch("eleSCEta", &eleSCEta);
   t1->Branch("eleSCPhi", &eleSCPhi);
   t1->Branch("eleSCRawEn", &eleSCRawEn);
   t1->Branch("eleSCEtaWidth", &eleSCEtaWidth);
   t1->Branch("eleSCPhiWidth", &eleSCPhiWidth);
   t1->Branch("eleHoverE", &eleHoverE);
   t1->Branch("eleHoverEBc", &eleHoverEBc);
   t1->Branch("eleEoverP", &eleEoverP);
   t1->Branch("eleEoverPInv", &eleEoverPInv);
   t1->Branch("eleEcalE", &eleEcalE);
   t1->Branch("elePAtVtx", &elePAtVtx);
   t1->Branch("elePAtSC", &elePAtSC);
   t1->Branch("elePAtCluster", &elePAtCluster);
   t1->Branch("elePAtSeed", &elePAtSeed);
   t1->Branch("eleBrem", &eleBrem);
   t1->Branch("eledEtaAtVtx", &eledEtaAtVtx);
   t1->Branch("eledPhiAtVtx", &eledPhiAtVtx);
   t1->Branch("eledEtaSeedAtVtx", &eledEtaSeedAtVtx);
   t1->Branch("eleSigmaIEtaIEta", &eleSigmaIEtaIEta);
   t1->Branch("eleSigmaIEtaIEta_2012", &eleSigmaIEtaIEta_2012);
   t1->Branch("eleSigmaIPhiIPhi", &eleSigmaIPhiIPhi);
   t1->Branch("eleMissHits", &eleMissHits);
   }
   
   if(keepPFcands){
  t1->Branch("nPFCands", &nPFCands);
  t1->Branch("pfcand_isProbe", &pfcand_isProbe);  
  t1->Branch("pfcand_pdgId", &pfcand_pdgId);
  t1->Branch("pfcand_charge", &pfcand_charge);  
  t1->Branch("pfcand_pt", &pfcand_pt);  
  t1->Branch("pfcand_eta", &pfcand_eta);  
  t1->Branch("pfcand_phi", &pfcand_phi);    
  t1->Branch("pfcand_vtx_x", &pfcand_vtx_x); 
  t1->Branch("pfcand_vtx_y", &pfcand_vtx_y);    
  t1->Branch("pfcand_vtx_z", &pfcand_vtx_z);      
   }
  if(keepPhotons){
  t1->Branch("nPho", &nPho);
  t1->Branch("phoE", &phoE);
  t1->Branch("phoEt", &phoEt);
  t1->Branch("phoEta", &phoEta);
  t1->Branch("phoPhi", &phoPhi);

  t1->Branch("phoEcorrStdEcal", &phoEcorrStdEcal);
  t1->Branch("phoEcorrPhoEcal", &phoEcorrPhoEcal);
  t1->Branch("phoEcorrRegr1", &phoEcorrRegr1);
  t1->Branch("phoEcorrRegr2", &phoEcorrRegr2);
  t1->Branch("phoEcorrErrStdEcal", &phoEcorrErrStdEcal);
  t1->Branch("phoEcorrErrPhoEcal", &phoEcorrErrPhoEcal);
  t1->Branch("phoEcorrErrRegr1", &phoEcorrErrRegr1);
  t1->Branch("phoEcorrErrRegr2", &phoEcorrErrRegr2);

  t1->Branch("phoSCE", &phoSCE);
  t1->Branch("phoSCEt", &phoSCEt);
  t1->Branch("phoSCRawE", &phoSCRawE);
  t1->Branch("phoSCEta", &phoSCEta);
  t1->Branch("phoSCPhi", &phoSCPhi);
  t1->Branch("phoSCEtaWidth", &phoSCEtaWidth);
  t1->Branch("phoSCPhiWidth", &phoSCPhiWidth);
  t1->Branch("phoSCBrem", &phoSCBrem);
  t1->Branch("phoSCnHits", &phoSCnHits);
  t1->Branch("phoSCflags", &phoSCflags);
  t1->Branch("phoSCinClean", &phoSCinClean);
  t1->Branch("phoSCinUnClean", &phoSCinUnClean);
  t1->Branch("phoSCnBC", &phoSCnBC);
  t1->Branch("phoESEn", &phoESEn);
    
  t1->Branch("phoIsPFPhoton", &phoIsPFPhoton);
  t1->Branch("phoIsStandardPhoton", &phoIsStandardPhoton);
  t1->Branch("phoHasPixelSeed", &phoHasPixelSeed);
  t1->Branch("phoHasConversionTracks", &phoHasConversionTracks);
    // tree_->Branch("phoEleVeto", &phoEleVeto);        // TODO: not available in reco::
  t1->Branch("phoHadTowerOverEm", &phoHadTowerOverEm);
  t1->Branch("phoHoverE", &phoHoverE);
  t1->Branch("phoHoverEValid", &phoHoverEValid);
  t1->Branch("phoSigmaIEtaIEta", &phoSigmaIEtaIEta);
  t1->Branch("phoR9", &phoR9);
  }
  
  if(keepCaloTowers){
  t1->Branch("nTower", &nTower);
  t1->Branch("maxHFp", &maxHFp);
  t1->Branch("maxHFm", &maxHFm);
  t1->Branch("CaloTower_hadE", &CaloTower_hadE);
  t1->Branch("CaloTower_emE", &CaloTower_emE);
  t1->Branch("CaloTower_e", &CaloTower_e);
  t1->Branch("CaloTower_et", &CaloTower_et);
  t1->Branch("CaloTower_eta", &CaloTower_eta);
  t1->Branch("CaloTower_phi", &CaloTower_phi);
  t1->Branch("CaloTower_HF_AboveThreshold", &CaloTower_HF_AboveThreshold );
  t1->Branch("CaloTower_Had_AboveThreshold", &CaloTower_Had_AboveThreshold );
  t1->Branch("CaloTower_EM_AboveThreshold", &CaloTower_EM_AboveThreshold );
  t1->Branch("CaloTower_HF_inNoiseRegion", &CaloTower_HF_inNoiseRegion );
  t1->Branch("CaloTower_Had_inNoiseRegion", &CaloTower_Had_inNoiseRegion );
  t1->Branch("CaloTower_EM_inNoiseRegion", &CaloTower_EM_inNoiseRegion );
  }
  if(keepZDC){
  t1->Branch("ZDC_n",&ZDC_n);
  t1->Branch("ZDC_e",ZDC_e,"ZDC_e[ZDC_n]/F");
  t1->Branch("ZDC_zside",ZDC_zside,"ZDC_zside[ZDC_n]/I");
  t1->Branch("ZDC_section",ZDC_section,"ZDC_section[ZDC_n]/I");
  t1->Branch("ZDC_channel",ZDC_channel,"ZDC_channel[ZDC_n]/I");
  t1->Branch("ZDC_PM_Total_Energy",&ZDC_PM_Total_Energy );  
  t1->Branch("ZDC_P_Total_Energy",&ZDC_P_Total_Energy );
  t1->Branch("ZDC_P_ECal_Energy",&ZDC_P_ECal_Energy );
  t1->Branch("ZDC_P_HCal_Energy",&ZDC_P_HCal_Energy );          
  t1->Branch("ZDC_M_Total_Energy",&ZDC_M_Total_Energy );
  t1->Branch("ZDC_M_ECal_Energy",&ZDC_M_ECal_Energy );
  t1->Branch("ZDC_M_HCal_Energy",&ZDC_M_HCal_Energy); 
  
  }
}


void Analysis_NtupleContent::SetTree_GenVtxStudy(TTree *mytree) { t2 = mytree; }

void Analysis_NtupleContent::SetTree_Test(TTree *mytree) {
   t3 = mytree; 
   t3->Branch("run", &run);
   }
void Analysis_NtupleContent::CreateBranches_GenVtxStudy() {
 
 
  t2->Branch("CutThrough_Num", &CutThrough_Num);
  
  t2->Branch("run", &run);
  t2->Branch("event", &event);
  t2->Branch("ls", &ls);
  t2->Branch("fromFullAOD", &fromFullAOD);
  t2->Branch("genWeight", &genWeight);
  t2->Branch("npairs", &npairs);
  t2->Branch("nvertices", &nvertices);
  t2->Branch("iprobe", &iprobe);
  t2->Branch("BSpot_x", &BSpot_x);
  t2->Branch("BSpot_y", &BSpot_y);
  t2->Branch("BSpot_z", &BSpot_z);
  t2->Branch("pv_x", &pv_x);
  t2->Branch("pv_y", &pv_y);
  t2->Branch("pv_z", &pv_z);
  t2->Branch("hasValidVertex", &hasValidVertex);
  t2->Branch("hasFakeVertex", &hasFakeVertex);
 
  t2->Branch("gentrk1_pt", &gentrk1_pt);
  t2->Branch("gentrk1_eta", &gentrk1_eta);
  t2->Branch("gentrk1_phi", &gentrk1_phi);
  t2->Branch("gentrk1_charge", &gentrk1_charge);
  t2->Branch("gentrk1_vtx_x", &gentrk1_vtx_x);
  t2->Branch("gentrk1_vtx_y", &gentrk1_vtx_y);
  t2->Branch("gentrk1_vtx_z", &gentrk1_vtx_z);
  t2->Branch("gentrk1_px", &gentrk1_px);
  t2->Branch("gentrk1_py", &gentrk1_py);
  t2->Branch("gentrk1_pz", &gentrk1_pz);
  t2->Branch("gentrk1_pdgId", &gentrk1_pdgId);  
  t2->Branch("gentrk2_pt", &gentrk2_pt);
  t2->Branch("gentrk2_eta", &gentrk2_eta);
  t2->Branch("gentrk2_phi", &gentrk2_phi);
  t2->Branch("gentrk2_charge", &gentrk2_charge);
  t2->Branch("gentrk2_vtx_x", &gentrk2_vtx_x);
  t2->Branch("gentrk2_vtx_y", &gentrk2_vtx_y);
  t2->Branch("gentrk2_vtx_z", &gentrk2_vtx_z);
  t2->Branch("gentrk2_px", &gentrk2_px);
  t2->Branch("gentrk2_py", &gentrk2_py);
  t2->Branch("gentrk2_pz", &gentrk2_pz);
  t2->Branch("gentrk2_pdgId", &gentrk2_pdgId); 
  t2->Branch("genFinalMass", &genFinalMass);
  t2->Branch("genFinalPt", &genFinalPt);
  t2->Branch("genFinalEta", &genFinalEta);
  t2->Branch("genFinalPhi", &genFinalPhi);
  t2->Branch("genFinalE", &genFinalE);
  t2->Branch("genFinalPx", &genFinalPx);
  t2->Branch("genFinalPy", &genFinalPy);
  t2->Branch("genFinalPz", &genFinalPz);  
  
  
  t2->Branch("indep_pt", &indep_pt);
  t2->Branch("indep_eta", &indep_eta);
  t2->Branch("indep_phi", &indep_phi);
  t2->Branch("indep_mass", &indep_mass);
  t2->Branch("gentau1_gamma_pz", &gentau1_gamma_pz);
  t2->Branch("gentau2_gamma_pz", &gentau2_gamma_pz);
  t2->Branch("gentau1_pt", &gentau1_pt);
  t2->Branch("gentau1_eta", &gentau1_eta);
  t2->Branch("gentau1_phi", &gentau1_phi);
  t2->Branch("gentau1_charge", &gentau1_charge);
  t2->Branch("gentau1_vtx_x", &gentau1_vtx_x);
  t2->Branch("gentau1_vtx_y", &gentau1_vtx_y);
  t2->Branch("gentau1_vtx_z", &gentau1_vtx_z);
  t2->Branch("gentau1_px", &gentau1_px);
  t2->Branch("gentau1_py", &gentau1_py);
  t2->Branch("gentau1_pz", &gentau1_pz);
  t2->Branch("gentau1_pdgId", &gentau1_pdgId);  
  t2->Branch("gentau1_decay_el", &gentau1_decay_el);
  t2->Branch("gentau1_decay_mu", &gentau1_decay_mu); 
  t2->Branch("gentau1_decay_1prong", &gentau1_decay_1prong); 
  t2->Branch("gentau1_decay_3prong", &gentau1_decay_3prong); 
  t2->Branch("gentau1_decay_other", &gentau1_decay_other); 
  t2->Branch("gentau1_decay", &gentau1_decay);
  t2->Branch("gentau1_N_pi0", &gentau1_N_pi0);  
  t2->Branch("gentau2_pt", &gentau2_pt);
  t2->Branch("gentau2_eta", &gentau2_eta);
  t2->Branch("gentau2_phi", &gentau2_phi);
  t2->Branch("gentau2_charge", &gentau2_charge);
  t2->Branch("gentau2_vtx_x", &gentau2_vtx_x);
  t2->Branch("gentau2_vtx_y", &gentau2_vtx_y);
  t2->Branch("gentau2_vtx_z", &gentau2_vtx_z);
  t2->Branch("gentau2_px", &gentau2_px);
  t2->Branch("gentau2_py", &gentau2_py);
  t2->Branch("gentau2_pz", &gentau2_pz);
  t2->Branch("gentau2_pdgId", &gentau2_pdgId);
  t2->Branch("gentau2_decay_el", &gentau2_decay_el);
  t2->Branch("gentau2_decay_mu", &gentau2_decay_mu); 
  t2->Branch("gentau2_decay_1prong", &gentau2_decay_1prong); 
  t2->Branch("gentau2_decay_3prong", &gentau2_decay_3prong); 
  t2->Branch("gentau2_decay_other", &gentau2_decay_other); 
  t2->Branch("gentau2_decay", &gentau2_decay); 
  t2->Branch("gentau2_N_pi0", &gentau2_N_pi0);   
  t2->Branch("genDiTauMass", &genDiTauMass);
  t2->Branch("genDiTauPt", &genDiTauPt);
  t2->Branch("genDiTauEta", &genDiTauEta);
  t2->Branch("genDiTauPhi", &genDiTauPhi);
  t2->Branch("genDiTauE", &genDiTauE);
  t2->Branch("genDiTauPx", &genDiTauPx);
  t2->Branch("genDiTauPy", &genDiTauPy);
  t2->Branch("genDiTauPz", &genDiTauPz);
  
  
  t2->Branch("nGen_tau1", &nGen_tau1 );
  t2->Branch("gen_pdgId_tau1", &gen_pdgId_tau1 );
  t2->Branch("gen_pT_tau1", &gen_pT_tau1 );
  t2->Branch("gen_eta_tau1", &gen_eta_tau1 );
  t2->Branch("gen_phi_tau1", &gen_phi_tau1 ); 
  t2->Branch("gen_charge_tau1", &gen_charge_tau1 );
  t2->Branch("gen_mass_tau1", &gen_mass_tau1 );  
  t2->Branch("gen_vtx_x_tau1", &gen_vtx_x_tau1 );
  t2->Branch("gen_vtx_y_tau1", &gen_vtx_y_tau1 );
  t2->Branch("gen_vtx_z_tau1", &gen_vtx_z_tau1 );
  t2->Branch("gen_status_tau1", &gen_status_tau1 );
  t2->Branch("gen_E_tau1", &gen_E_tau1 );
  t2->Branch("gen_Et_tau1", &gen_Et_tau1 );

  t2->Branch("nGen_tau2", &nGen_tau2 );
  t2->Branch("gen_pdgId_tau2", &gen_pdgId_tau2 );
  t2->Branch("gen_pT_tau2", &gen_pT_tau2 );
  t2->Branch("gen_eta_tau2", &gen_eta_tau2 );
  t2->Branch("gen_phi_tau2", &gen_phi_tau2 ); 
  t2->Branch("gen_charge_tau2", &gen_charge_tau2 );
  t2->Branch("gen_mass_tau2", &gen_mass_tau2 );  
  t2->Branch("gen_vtx_x_tau2", &gen_vtx_x_tau2 );
  t2->Branch("gen_vtx_y_tau2", &gen_vtx_y_tau2 );
  t2->Branch("gen_vtx_z_tau2", &gen_vtx_z_tau2 );
  t2->Branch("gen_status_tau2", &gen_status_tau2 );
  t2->Branch("gen_E_tau2", &gen_E_tau2 );
  t2->Branch("gen_Et_tau2", &gen_Et_tau2 );

  // t2->Branch("genDiTau_pair_svprob", &genDiTau_pair_svprob);
  // t2->Branch("genDiTau_pair_fit_mass", &genDiTau_pair_fit_mass);
  // t2->Branch("genDiTau_pair_normalchi2", &genDiTau_pair_normalchi2);
  t2->Branch("genDiTau_vtx_x", &genDiTau_vtx_x);
  t2->Branch("genDiTau_vtx_y", &genDiTau_vtx_y);
  t2->Branch("genDiTau_vtx_z", &genDiTau_vtx_z);


  // t2->Branch("genFinal_pair_svprob", &genFinal_pair_svprob);
  // t2->Branch("genFinal_pair_fit_mass", &genFinal_pair_fit_mass);
  // t2->Branch("genFinal_pair_normalchi2", &genFinal_pair_normalchi2);  
  t2->Branch("genFinal_vtx_x", &genFinal_vtx_x);
  t2->Branch("genFinal_vtx_y", &genFinal_vtx_y);
  t2->Branch("genFinal_vtx_z", &genFinal_vtx_z);


  t2->Branch("genDiTau_vtx_gentau1_IP", &genDiTau_vtx_gentau1_IP);
  t2->Branch("genDiTau_vtx_gentau2_IP", &genDiTau_vtx_gentau2_IP);
  t2->Branch("genDiTau_vtx_gentrk1_IP", &genDiTau_vtx_gentrk1_IP);
  t2->Branch("genDiTau_vtx_gentrk2_IP", &genDiTau_vtx_gentrk2_IP);
  
  t2->Branch("genFinal_vtx_gentau1_IP", &genFinal_vtx_gentau1_IP);
  t2->Branch("genFinal_vtx_gentau2_IP", &genFinal_vtx_gentau2_IP);
  t2->Branch("genFinal_vtx_gentrk1_IP", &genFinal_vtx_gentrk1_IP);
  t2->Branch("genFinal_vtx_gentrk2_IP", &genFinal_vtx_gentrk2_IP);

  // t2->Branch("genDiTau_refit_vtx_gentau1_IP", &genDiTau_refit_vtx_gentau1_IP);
  // t2->Branch("genDiTau_refit_vtx_gentau2_IP", &genDiTau_refit_vtx_gentau2_IP);
  // t2->Branch("genDiTau_refit_vtx_gentrk1_IP", &genDiTau_refit_vtx_gentrk1_IP);
  // t2->Branch("genDiTau_refit_vtx_gentrk2_IP", &genDiTau_refit_vtx_gentrk2_IP);
  
  // t2->Branch("genFinal_refit_vtx_gentau1_IP", &genFinal_refit_vtx_gentau1_IP);
  // t2->Branch("genFinal_refit_vtx_gentau2_IP", &genFinal_refit_vtx_gentau2_IP);
  // t2->Branch("genFinal_refit_vtx_gentrk1_IP", &genFinal_refit_vtx_gentrk1_IP);
  // t2->Branch("genFinal_refit_vtx_gentrk2_IP", &genFinal_refit_vtx_gentrk2_IP);
  

}
void Analysis_NtupleContent::CreateExtraTrgBranches(const std::vector<std::string> &HLTs, bool isTag = false) {
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


void Analysis_NtupleContent::CreateExtraTrgBranches_3prong(const std::vector<std::string> &HLTs, bool isTag = false) {
  for (unsigned int ihlt = 0; ihlt < HLTs.size(); ihlt++) {
    if (isTag) {
      t4->Branch(TString("tag_" + HLTs[ihlt]), &tag_trg[ihlt]);
      t4->Branch(TString("tag_" + HLTs[ihlt] + "_pt"), &tag_trg_pt[ihlt]);
      t4->Branch(TString("tag_" + HLTs[ihlt] + "_eta"), &tag_trg_eta[ihlt]);
      t4->Branch(TString("tag_" + HLTs[ihlt] + "_phi"), &tag_trg_phi[ihlt]);
      t4->Branch(TString("tag_" + HLTs[ihlt] + "_dr"), &tag_trg_dr[ihlt]);
    }
    // else {
      // t1->Branch(TString("probe_" + HLTs[ihlt]), &probe_trg[ihlt]);
      // t1->Branch(TString("probe_" + HLTs[ihlt] + "_pt"), &probe_trg_pt[ihlt]);
      // t1->Branch(TString("probe_" + HLTs[ihlt] + "_eta"), &probe_trg_eta[ihlt]);
      // t1->Branch(TString("probe_" + HLTs[ihlt] + "_phi"), &probe_trg_phi[ihlt]);
      // t1->Branch(TString("probe_" + HLTs[ihlt] + "_dr"), &probe_trg_dr[ihlt]);
    // }
  }
}


  void Analysis_NtupleContent::SetTree_3ProngStudy(TTree *mytree) { t4 = mytree; }
  void Analysis_NtupleContent::CreateBranches_3ProngStudy(const std::vector<std::string> &HLTs,
                                   const std::vector<std::string> &selectorNames){
     
  t4->Branch("run", &run);
  t4->Branch("event", &event);
  t4->Branch("ls", &ls);
  t4->Branch("fromFullAOD", &fromFullAOD);
  t4->Branch("genWeight", &genWeight);
  t4->Branch("BSpot_x", &BSpot_x);
  t4->Branch("BSpot_y", &BSpot_y);
  t4->Branch("BSpot_z", &BSpot_z);
  t4->Branch("pv_x", &pv_x);
  t4->Branch("pv_y", &pv_y);
  t4->Branch("pv_z", &pv_z);
  t4->Branch("nVertices", &nvertices);
  t4->Branch("hasValidVertex", &hasValidVertex);
  t4->Branch("hasFakeVertex", &hasFakeVertex);
  t4->Branch("nTrueInteractions", &trueNumInteractions);
  t4->Branch("nPUInteractions", &puNumInteractions);
  t4->Branch("rho", &Rho);
  t4->Branch("ntag", &ntag);
  t4->Branch("ntag_electron", &ntag_electron);
  t4->Branch("npairs", &npairs);
  t4->Branch("npairs_3prong", &npairs_3prong);
  t4->Branch("tag_isMuon", &tag_isMuon);
  t4->Branch("tag_isElectron", &tag_isElectron);
  t4->Branch("iprobe_3prong", &iprobe_3prong);
  
  t4->Branch("gentrk1_pt", &gentrk1_pt);
  t4->Branch("gentrk1_eta", &gentrk1_eta);
  t4->Branch("gentrk1_phi", &gentrk1_phi);
  t4->Branch("gentrk1_charge", &gentrk1_charge);
  t4->Branch("gentrk1_M", &gentrk1_M);
  t4->Branch("gentrk1_vtx_x", &gentrk1_vtx_x);
  t4->Branch("gentrk1_vtx_y", &gentrk1_vtx_y);
  t4->Branch("gentrk1_vtx_z", &gentrk1_vtx_z);
  t4->Branch("gentrk1_pdgId", &gentrk1_pdgId);  
  t4->Branch("gentrk2_pt", &gentrk2_pt);
  t4->Branch("gentrk2_eta", &gentrk2_eta);
  t4->Branch("gentrk2_phi", &gentrk2_phi);
  t4->Branch("gentrk2_charge", &gentrk2_charge);
  t4->Branch("gentrk2_M", &gentrk2_M);
  t4->Branch("gentrk2_vtx_x", &gentrk2_vtx_x);
  t4->Branch("gentrk2_vtx_y", &gentrk2_vtx_y);
  t4->Branch("gentrk2_vtx_z", &gentrk2_vtx_z);
  t4->Branch("gentrk2_pdgId", &gentrk2_pdgId); 
  t4->Branch("genFinalMass", &genFinalMass);
  t4->Branch("genFinalPt", &genFinalPt);
  t4->Branch("genFinalEta", &genFinalEta);
  t4->Branch("genFinalPhi", &genFinalPhi);
  t4->Branch("genFinalE", &genFinalE);
  t4->Branch("genFinalPx", &genFinalPx);
  t4->Branch("genFinalPy", &genFinalPy);
  t4->Branch("genFinalPz", &genFinalPz);  
  

  t4->Branch("gentau1_pt", &gentau1_pt);
  t4->Branch("gentau1_eta", &gentau1_eta);
  t4->Branch("gentau1_phi", &gentau1_phi);
  t4->Branch("gentau1_charge", &gentau1_charge);
  t4->Branch("gentau1_vtx_x", &gentau1_vtx_x);
  t4->Branch("gentau1_vtx_y", &gentau1_vtx_y);
  t4->Branch("gentau1_vtx_z", &gentau1_vtx_z);
  t4->Branch("gentau1_pdgId", &gentau1_pdgId);  
  t4->Branch("gentau1_decay_el", &gentau1_decay_el);
  t4->Branch("gentau1_decay_mu", &gentau1_decay_mu); 
  t4->Branch("gentau1_decay_1prong", &gentau1_decay_1prong); 
  t4->Branch("gentau1_decay_3prong", &gentau1_decay_3prong); 
  t4->Branch("gentau1_decay_other", &gentau1_decay_other); 
  t4->Branch("gentau1_decay", &gentau1_decay); 
  t4->Branch("gentau1_N_pi0", &gentau1_N_pi0);   
  t4->Branch("gentau2_pt", &gentau2_pt);
  t4->Branch("gentau2_eta", &gentau2_eta);
  t4->Branch("gentau2_phi", &gentau2_phi);
  t4->Branch("gentau2_charge", &gentau2_charge);
  t4->Branch("gentau2_vtx_x", &gentau2_vtx_x);
  t4->Branch("gentau2_vtx_y", &gentau2_vtx_y);
  t4->Branch("gentau2_vtx_z", &gentau2_vtx_z);
  t4->Branch("gentau2_pdgId", &gentau2_pdgId);
  t4->Branch("gentau2_decay_el", &gentau2_decay_el);
  t4->Branch("gentau2_decay_mu", &gentau2_decay_mu); 
  t4->Branch("gentau2_decay_1prong", &gentau2_decay_1prong); 
  t4->Branch("gentau2_decay_3prong", &gentau2_decay_3prong); 
  t4->Branch("gentau2_decay_other", &gentau2_decay_other); 
  t4->Branch("gentau2_decay", &gentau2_decay);  
  t4->Branch("gentau2_N_pi0", &gentau2_N_pi0);   
  t4->Branch("genDiTauMass", &genDiTauMass);
  t4->Branch("genDiTauPt", &genDiTauPt);
  t4->Branch("genDiTauEta", &genDiTauEta);
  t4->Branch("genDiTauPhi", &genDiTauPhi);
  t4->Branch("genDiTauE", &genDiTauE);
  t4->Branch("genDiTauPx", &genDiTauPx);
  t4->Branch("genDiTauPy", &genDiTauPy);
  t4->Branch("genDiTauPz", &genDiTauPz);
  
  t4->Branch("nGen", &nGen );
  t4->Branch("gen_pdgId", &gen_pdgId );
  t4->Branch("gen_pT", &gen_pT );
  t4->Branch("gen_eta", &gen_eta );
  t4->Branch("gen_phi", &gen_phi ); 
  t4->Branch("gen_charge", &gen_charge );
  t4->Branch("gen_mass", &gen_mass );  
  t4->Branch("gen_vtx_x", &gen_vtx_x );
  t4->Branch("gen_vtx_y", &gen_vtx_y );
  t4->Branch("gen_vtx_z", &gen_vtx_z );
  t4->Branch("gen_status", &gen_status );
  t4->Branch("gen_E", &gen_E );
  t4->Branch("gen_Et", &gen_Et );  
  
  
  t4->Branch("nGen_tau1", &nGen_tau1 );
  t4->Branch("gen_pdgId_tau1", &gen_pdgId_tau1 );
  t4->Branch("gen_pT_tau1", &gen_pT_tau1 );
  t4->Branch("gen_eta_tau1", &gen_eta_tau1 );
  t4->Branch("gen_phi_tau1", &gen_phi_tau1 ); 
  t4->Branch("gen_charge_tau1", &gen_charge_tau1 );
  t4->Branch("gen_mass_tau1", &gen_mass_tau1 );  
  t4->Branch("gen_vtx_x_tau1", &gen_vtx_x_tau1 );
  t4->Branch("gen_vtx_y_tau1", &gen_vtx_y_tau1 );
  t4->Branch("gen_vtx_z_tau1", &gen_vtx_z_tau1 );
  t4->Branch("gen_status_tau1", &gen_status_tau1 );
  t4->Branch("gen_E_tau1", &gen_E_tau1 );
  t4->Branch("gen_Et_tau1", &gen_Et_tau1 );

  t4->Branch("nGen_tau2", &nGen_tau2 );
  t4->Branch("gen_pdgId_tau2", &gen_pdgId_tau2 );
  t4->Branch("gen_pT_tau2", &gen_pT_tau2 );
  t4->Branch("gen_eta_tau2", &gen_eta_tau2 );
  t4->Branch("gen_phi_tau2", &gen_phi_tau2 ); 
  t4->Branch("gen_charge_tau2", &gen_charge_tau2 );
  t4->Branch("gen_mass_tau2", &gen_mass_tau2 );  
  t4->Branch("gen_vtx_x_tau2", &gen_vtx_x_tau2 );
  t4->Branch("gen_vtx_y_tau2", &gen_vtx_y_tau2 );
  t4->Branch("gen_vtx_z_tau2", &gen_vtx_z_tau2 );
  t4->Branch("gen_status_tau2", &gen_status_tau2 );
  t4->Branch("gen_E_tau2", &gen_E_tau2 );
  t4->Branch("gen_Et_tau2", &gen_Et_tau2 );

  t4->Branch("gentrk1_match_dr", &gentrk1_match_dr );
  t4->Branch("gentrk1_match_dphi", &gentrk1_match_dphi );
  t4->Branch("gentrk1_match_diff_eta", &gentrk1_match_diff_eta );
  t4->Branch("gentrk1_match_diff_pt", &gentrk1_match_diff_pt ); 
  t4->Branch("gentrk1_diff_vtx_x", &gentrk1_diff_vtx_x );
  t4->Branch("gentrk1_diff_vtx_y", &gentrk1_diff_vtx_y );
  t4->Branch("gentrk1_diff_vtx_z", & gentrk1_diff_vtx_z);
  t4->Branch("gentrk1_isTag", &gentrk1_isTag );
  
  t4->Branch("gentrk2_match_dr", &gentrk2_match_dr );
  t4->Branch("gentrk2_match_dphi", &gentrk2_match_dphi );
  t4->Branch("gentrk2_match_diff_eta", &gentrk2_match_diff_eta );
  t4->Branch("gentrk2_match_diff_pt", &gentrk2_match_diff_pt ); 
  t4->Branch("gentrk2_diff_vtx_x", &gentrk2_diff_vtx_x );
  t4->Branch("gentrk2_diff_vtx_y", &gentrk2_diff_vtx_y );
  t4->Branch("gentrk2_diff_vtx_z", &gentrk2_diff_vtx_z );  
  t4->Branch("gentrk2_isTag", &gentrk2_isTag );
  
  
  for (unsigned int ihlt = 0; ihlt < HLTs.size(); ihlt++)
    t4->Branch(TString(HLTs[ihlt]), &trigger[ihlt]);

  // Tag specific
  t4->Branch("tag_pt", &tag_pt);
  t4->Branch("tag_eta", &tag_eta);
  t4->Branch("tag_phi", &tag_phi);
  t4->Branch("tag_charge", &tag_charge);
  t4->Branch("tag_pterr", &tag_pterr);
  t4->Branch("tag_dxy", &tag_dxy);
  t4->Branch("tag_dz", &tag_dz);
  t4->Branch("tag_vtx_x", &tag_vtx_x);
  t4->Branch("tag_vtx_y", &tag_vtx_y);
  t4->Branch("tag_vtx_z", &tag_vtx_z);
  t4->Branch("tag_isPF", &tag_isPF);
  t4->Branch("tag_isSA", &tag_isSA);
  t4->Branch("tag_isdSA", &tag_isdSA);
  t4->Branch("tag_isTracker", &tag_isTracker);
  t4->Branch("tag_isGlobal", &tag_isGlobal);
  t4->Branch("tag_isLoose", &tag_isLoose);
  t4->Branch("tag_isMedium", &tag_isMedium);
  t4->Branch("tag_isTight", &tag_isTight);
  t4->Branch("tag_isSoft", &tag_isSoft);
  t4->Branch("tag_isHighPt", &tag_isHighPt);
  t4->Branch("tag_relIso04", &tag_relIso04);
  t4->Branch("tag_miniIso", &tag_miniIso);
  t4->Branch("tag_miniIsoCharged", &tag_miniIsoCharged);
  t4->Branch("tag_miniIsoPhotons", &tag_miniIsoPhotons);
  t4->Branch("tag_miniIsoNeutrals", &tag_miniIsoPhotons);
  t4->Branch("tag_isMatchedGen", &tag_isMatchedGen);
  t4->Branch("tag_minDR", &tag_minDR);
  t4->Branch("tag_ptRel_minDR", &tag_ptRel_minDR);
  t4->Branch("tag_iso03_sumPt", &tag_iso03_sumPt);
  t4->Branch("tag_pfIso04_charged", &tag_pfIso04_charged);
  t4->Branch("tag_pfIso04_neutral", &tag_pfIso04_neutral);
  t4->Branch("tag_pfIso04_photon", &tag_pfIso04_photon);
  t4->Branch("tag_pfIso04_sumPU", &tag_pfIso04_sumPU);
  t4->Branch("tag_tuneP_pt", &tag_tuneP_pt);
  t4->Branch("tag_tuneP_pterr", &tag_tuneP_pterr);
  t4->Branch("tag_nsegments", &tag_nsegments);
  t4->Branch("tag_hasTrackMatch", &tag_hasTrackMatch);
  t4->Branch("tag_TrackMatchDR", &tag_TrackMatchDR);
 

// New 3 prong variables 


  
  t4->Branch("tau_3prong_trk1_pt", &tau_3prong_trk1_pt );
  t4->Branch("tau_3prong_trk1_eta", &tau_3prong_trk1_eta );
  t4->Branch("tau_3prong_trk1_phi", &tau_3prong_trk1_phi );
  t4->Branch("tau_3prong_trk1_charge", &tau_3prong_trk1_charge );
  t4->Branch("tau_3prong_trk1_vtx_x", &tau_3prong_trk1_vtx_x );
  t4->Branch("tau_3prong_trk1_vtx_y", &tau_3prong_trk1_vtx_y );
  t4->Branch("tau_3prong_trk1_vtx_z", &tau_3prong_trk1_vtx_z );
        
  t4->Branch("tau_3prong_trk2_pt", &tau_3prong_trk2_pt );
  t4->Branch("tau_3prong_trk2_eta", &tau_3prong_trk2_eta );
  t4->Branch("tau_3prong_trk2_phi", &tau_3prong_trk2_phi );
   t4->Branch("tau_3prong_trk2_charge", &tau_3prong_trk2_charge );
  t4->Branch("tau_3prong_trk2_vtx_x", &tau_3prong_trk2_vtx_x);
  t4->Branch("tau_3prong_trk2_vtx_y", &tau_3prong_trk2_vtx_y );
  t4->Branch("tau_3prong_trk2_vtx_z", &tau_3prong_trk2_vtx_z);
     
  t4->Branch("tau_3prong_trk3_pt", &tau_3prong_trk3_pt );
  t4->Branch("tau_3prong_trk3_eta", &tau_3prong_trk3_eta );
  t4->Branch("tau_3prong_trk3_phi", &tau_3prong_trk3_phi );
   t4->Branch("tau_3prong_trk3_charge", &tau_3prong_trk3_charge );
  t4->Branch("tau_3prong_trk3_vtx_x", &tau_3prong_trk3_vtx_x );
  t4->Branch("tau_3prong_trk3_vtx_y", &tau_3prong_trk3_vtx_y );
  t4->Branch("tau_3prong_trk3_vtx_z", &tau_3prong_trk3_vtx_z );
        
  t4->Branch("tau_3prong_total_pt", &tau_3prong_total_pt );
  t4->Branch("tau_3prong_total_eta", &tau_3prong_total_eta );
  t4->Branch("tau_3prong_total_phi", &tau_3prong_total_phi );
  t4->Branch("tau_3prong_total_M", &tau_3prong_total_M );
  t4->Branch("tau_3prong_total_pt_sum", &tau_3prong_total_pt_sum );
   t4->Branch("tau_3prong_total_charge", &tau_3prong_total_charge );
   t4->Branch("tau_3prong_total_Dz", &tau_3prong_total_Dz );
  t4->Branch("tau_3prong_total_vtx_x", &tau_3prong_total_vtx_x );
  t4->Branch("tau_3prong_total_vtx_y", &tau_3prong_total_vtx_y );
  t4->Branch("tau_3prong_total_vtx_z", &tau_3prong_total_vtx_z );
  t4->Branch("tau_3prong_total_vtx_prob", &tau_3prong_total_vtx_prob );
  t4->Branch("tau_3prong_total_vtx_chi2", &tau_3prong_total_vtx_chi2 );
  t4->Branch("tau_3prong_DR_1_2", &tau_3prong_DR_1_2 );
  t4->Branch("tau_3prong_DR_1_3", &tau_3prong_DR_1_3 );
  t4->Branch("tau_3prong_DR_2_3", &tau_3prong_DR_2_3 );
  t4->Branch("tau_3prong_DR_tag_1", &tau_3prong_DR_tag_1 );
  t4->Branch("tau_3prong_DR_tag_2", &tau_3prong_DR_tag_2 );
  t4->Branch("tau_3prong_DR_tag_3", &tau_3prong_DR_tag_3 );
  t4->Branch("tau_3prong_total_area", &tau_3prong_total_area );
  
  
 
  t4->Branch("pair_rank_vtx_prob_3prong", &pair_rank_vtx_prob_3prong );
  t4->Branch("pair_rank_dPhi_muons_3prong", &pair_rank_dPhi_muons_3prong );
  t4->Branch("pair_rank_Mass_Mmumu_3prong", &pair_rank_Mass_Mmumu_3prong );


  t4->Branch("pair_pt", &pair_pt);
  t4->Branch("pair_eta", &pair_eta);
  t4->Branch("pair_phi", &pair_phi);
  t4->Branch("pair_mass", &pair_mass);
  t4->Branch("pair_fit_mass", &pair_fit_mass);
  t4->Branch("pair_svprob", &pair_svprob);
  t4->Branch("pair_normalchi2", &pair_normalchi2);
  t4->Branch("pair_dz", &pair_dz);
  t4->Branch("pair_dR", &pair_dR);
  t4->Branch("pair_dphi", &pair_dphi);
  t4->Branch("pair_rank_vtx_prob_3prong", &pair_rank_vtx_prob_3prong);
  t4->Branch("pair_rank_dPhi_muons_3prong", &pair_rank_dPhi_muons_3prong);
  t4->Branch("pair_rank_Mass_Mmumu_3prong", &pair_rank_Mass_Mmumu_3prong);

 //commented out to save space, added useful details for HIUPC collisions
    if(keepMuons){
    t4->Branch("nMu", &nMu);
    t4->Branch("mu_isProbe", &mu_isProbe);
    t4->Branch("mu_isTag", &mu_isTag);    
    t4->Branch("muPt", &muPt);
    t4->Branch("muEta", &muEta);
    t4->Branch("muPhi", &muPhi);
    t4->Branch("muCharge", &muCharge);
    t4->Branch("muType", &muType);
    t4->Branch("muIsGood", &muIsGood);
    
    t4->Branch("muIsGlobal", &muIsGlobal);
    t4->Branch("muIsTracker", &muIsTracker);
    t4->Branch("muIsPF", &muIsPF);
    t4->Branch("muIsSTA", &muIsSTA);

    t4->Branch("muD0", &muD0);
    t4->Branch("muDz", &muDz);
    t4->Branch("muIP3D", &muIP3D);
    t4->Branch("muD0Err", &muD0Err);
    t4->Branch("muDzErr", &muDzErr);
    t4->Branch("muIP3DErr", &muIP3DErr);
    t4->Branch("muChi2NDF", &muChi2NDF);
    t4->Branch("muInnerD0", &muInnerD0);
    t4->Branch("muInnerDz", &muInnerDz);
    
    t4->Branch("muInnerD0Err", &muInnerD0Err);
    t4->Branch("muInnerDzErr", &muInnerDzErr);
    t4->Branch("muInnerPt", &muInnerPt);
    t4->Branch("muInnerPtErr", &muInnerPtErr);
    t4->Branch("muInnerEta", &muInnerEta);

    t4->Branch("muTrkLayers", &muTrkLayers);
    t4->Branch("muPixelLayers", &muPixelLayers);
    t4->Branch("muPixelHits", &muPixelHits);
    t4->Branch("muMuonHits", &muMuonHits);
    t4->Branch("muTrkQuality", &muTrkQuality);
    t4->Branch("muStations", &muStations);
    t4->Branch("muIsoTrk", &muIsoTrk);
    t4->Branch("muPFChIso", &muPFChIso);
    t4->Branch("muPFPhoIso", &muPFPhoIso);
    t4->Branch("muPFNeuIso", &muPFNeuIso);
    t4->Branch("muPFPUIso", &muPFPUIso);
    t4->Branch("muIDSoft", &muIDSoft);
    t4->Branch("muIDLoose", &muIDLoose);
    t4->Branch("muIDMedium", &muIDMedium);
    t4->Branch("muIDMediumPrompt", &muIDMediumPrompt);
    t4->Branch("muIDTight", &muIDTight);
    t4->Branch("muIDGlobalHighPt", &muIDGlobalHighPt);
    t4->Branch("muIDTrkHighPt", &muIDTrkHighPt);
    t4->Branch("muIDInTime", &muIDInTime);

    }
    
    if(keepTracks){
    t4->Branch("nTrk", &nTrk);
    t4->Branch("trkisTag", &trkisTag);
    t4->Branch("trkisProbe", &trkisProbe);
    t4->Branch("trkPt", &trkPt);
    t4->Branch("trkP", &trkP);
    t4->Branch("trkEta", &trkEta);
    t4->Branch("trkPhi", &trkPhi);
    t4->Branch("trkcharge", &trkcharge);
    t4->Branch("trkvx", &trkvx);
    t4->Branch("trkvy", &trkvy);
    t4->Branch("trkvz", &trkvz);
    t4->Branch("trknormchi2", &trknormchi2);
    t4->Branch("trkchi2", &trkchi2);
    t4->Branch("trkd0", &trkd0);
    t4->Branch("trkdxy", &trkdxy);
    t4->Branch("trkdz", &trkdz);
    // t4->Branch("trkdxyBS", &trkdxyBS);
    // t4->Branch("trkdzBS", &trkdzBS);
    t4->Branch("trkdxyError", &trkdxyError);
    t4->Branch("trkdzError", &trkdzError);
    t4->Branch("trkValidHits", &trkValidHits);
    t4->Branch("trkMissHits", &trkMissHits);
    t4->Branch("trkPurity", &trkPurity);
    }


   if(keepElectrons){
   t4->Branch("nEle", &nEle);  
   t4->Branch("ele_isProbe", &ele_isProbe);
   t4->Branch("ele_isTag", &ele_isTag);
   t4->Branch("eleCharge", &eleCharge);
   t4->Branch("eleChargeConsistent", &eleChargeConsistent);
   t4->Branch("eleSCPixCharge", &eleSCPixCharge);
   // t4->Branch("eleCtfCharge", &eleCtfCharge);
   t4->Branch("eleEn", &eleEn);
   t4->Branch("eleD0", &eleD0);
   t4->Branch("eleDz", &eleDz);
   t4->Branch("eleIP3D", &eleIP3D);
   t4->Branch("eleD0Err", &eleD0Err);
   t4->Branch("eleDzErr", &eleDzErr);
   t4->Branch("eleIP3DErr", &eleIP3DErr);
   t4->Branch("eleTrkPt", &eleTrkPt);
   t4->Branch("eleTrkEta", &eleTrkEta);
   t4->Branch("eleTrkPhi", &eleTrkPhi);
   t4->Branch("eleTrkCharge", &eleTrkCharge);
   t4->Branch("eleTrkPtErr", &eleTrkPtErr);
   t4->Branch("eleTrkChi2", &eleTrkChi2);
   t4->Branch("eleTrkNdof", &eleTrkNdof);
   t4->Branch("eleTrkNormalizedChi2", &eleTrkNormalizedChi2);
   t4->Branch("eleTrkValidHits", &eleTrkValidHits);
   t4->Branch("eleTrkLayers", &eleTrkLayers);

   t4->Branch("elePt", &elePt);
   t4->Branch("eleEta", &eleEta);
   t4->Branch("elePhi", &elePhi);
   t4->Branch("eleSCEn", &eleSCEn);
   t4->Branch("eleESEn", &eleESEn);
   t4->Branch("eleSCEta", &eleSCEta);
   t4->Branch("eleSCPhi", &eleSCPhi);
   t4->Branch("eleSCRawEn", &eleSCRawEn);
   t4->Branch("eleSCEtaWidth", &eleSCEtaWidth);
   t4->Branch("eleSCPhiWidth", &eleSCPhiWidth);
   t4->Branch("eleHoverE", &eleHoverE);
   t4->Branch("eleHoverEBc", &eleHoverEBc);
   t4->Branch("eleEoverP", &eleEoverP);
   t4->Branch("eleEoverPInv", &eleEoverPInv);
   t4->Branch("eleEcalE", &eleEcalE);
   t4->Branch("elePAtVtx", &elePAtVtx);
   t4->Branch("elePAtSC", &elePAtSC);
   t4->Branch("elePAtCluster", &elePAtCluster);
   t4->Branch("elePAtSeed", &elePAtSeed);
   t4->Branch("eleBrem", &eleBrem);
   t4->Branch("eledEtaAtVtx", &eledEtaAtVtx);
   t4->Branch("eledPhiAtVtx", &eledPhiAtVtx);
   t4->Branch("eledEtaSeedAtVtx", &eledEtaSeedAtVtx);
   t4->Branch("eleSigmaIEtaIEta", &eleSigmaIEtaIEta);
   t4->Branch("eleSigmaIEtaIEta_2012", &eleSigmaIEtaIEta_2012);
   t4->Branch("eleSigmaIPhiIPhi", &eleSigmaIPhiIPhi);
   t4->Branch("eleMissHits", &eleMissHits);
   }
   
   if(keepPFcands){
  t4->Branch("nPFCands", &nPFCands);
  t4->Branch("pfcand_isProbe", &pfcand_isProbe);  
  t4->Branch("pfcand_pdgId", &pfcand_pdgId);
  t4->Branch("pfcand_charge", &pfcand_charge);  
  t4->Branch("pfcand_pt", &pfcand_pt);  
  t4->Branch("pfcand_eta", &pfcand_eta);  
  t4->Branch("pfcand_phi", &pfcand_phi);    
  t4->Branch("pfcand_vtx_x", &pfcand_vtx_x); 
  t4->Branch("pfcand_vtx_y", &pfcand_vtx_y);    
  t4->Branch("pfcand_vtx_z", &pfcand_vtx_z);      
   }
  if(keepPhotons){
  t4->Branch("nPho", &nPho);
  t4->Branch("phoE", &phoE);
  t4->Branch("phoEt", &phoEt);
  t4->Branch("phoEta", &phoEta);
  t4->Branch("phoPhi", &phoPhi);

  t4->Branch("phoEcorrStdEcal", &phoEcorrStdEcal);
  t4->Branch("phoEcorrPhoEcal", &phoEcorrPhoEcal);
  t4->Branch("phoEcorrRegr1", &phoEcorrRegr1);
  t4->Branch("phoEcorrRegr2", &phoEcorrRegr2);
  t4->Branch("phoEcorrErrStdEcal", &phoEcorrErrStdEcal);
  t4->Branch("phoEcorrErrPhoEcal", &phoEcorrErrPhoEcal);
  t4->Branch("phoEcorrErrRegr1", &phoEcorrErrRegr1);
  t4->Branch("phoEcorrErrRegr2", &phoEcorrErrRegr2);

  t4->Branch("phoSCE", &phoSCE);
  t4->Branch("phoSCEt", &phoSCEt);
  t4->Branch("phoSCRawE", &phoSCRawE);
  t4->Branch("phoSCEta", &phoSCEta);
  t4->Branch("phoSCPhi", &phoSCPhi);
  t4->Branch("phoSCEtaWidth", &phoSCEtaWidth);
  t4->Branch("phoSCPhiWidth", &phoSCPhiWidth);
  t4->Branch("phoSCBrem", &phoSCBrem);
  t4->Branch("phoSCnHits", &phoSCnHits);
  t4->Branch("phoSCflags", &phoSCflags);
  t4->Branch("phoSCinClean", &phoSCinClean);
  t4->Branch("phoSCinUnClean", &phoSCinUnClean);
  t4->Branch("phoSCnBC", &phoSCnBC);
  t4->Branch("phoESEn", &phoESEn);
    
  t4->Branch("phoIsPFPhoton", &phoIsPFPhoton);
  t4->Branch("phoIsStandardPhoton", &phoIsStandardPhoton);
  t4->Branch("phoHasPixelSeed", &phoHasPixelSeed);
  t4->Branch("phoHasConversionTracks", &phoHasConversionTracks);
    // tree_->Branch("phoEleVeto", &phoEleVeto);        // TODO: not available in reco::
  t4->Branch("phoHadTowerOverEm", &phoHadTowerOverEm);
  t4->Branch("phoHoverE", &phoHoverE);
  t4->Branch("phoHoverEValid", &phoHoverEValid);
  t4->Branch("phoSigmaIEtaIEta", &phoSigmaIEtaIEta);
  t4->Branch("phoR9", &phoR9);
  }
  
  if(keepCaloTowers){
  t4->Branch("nTower", &nTower);
  t4->Branch("maxHFp", &maxHFp);
  t4->Branch("maxHFm", &maxHFm);
  t4->Branch("CaloTower_hadE", &CaloTower_hadE);
  t4->Branch("CaloTower_emE", &CaloTower_emE);
  t4->Branch("CaloTower_e", &CaloTower_e);
  t4->Branch("CaloTower_et", &CaloTower_et);
  t4->Branch("CaloTower_eta", &CaloTower_eta);
  t4->Branch("CaloTower_phi", &CaloTower_phi);
  t4->Branch("CaloTower_HF_AboveThreshold", &CaloTower_HF_AboveThreshold );
  t4->Branch("CaloTower_Had_AboveThreshold", &CaloTower_Had_AboveThreshold );
  t4->Branch("CaloTower_EM_AboveThreshold", &CaloTower_EM_AboveThreshold );
  t4->Branch("CaloTower_HF_inNoiseRegion", &CaloTower_HF_inNoiseRegion );
  t4->Branch("CaloTower_Had_inNoiseRegion", &CaloTower_Had_inNoiseRegion );
  t4->Branch("CaloTower_EM_inNoiseRegion", &CaloTower_EM_inNoiseRegion );
  }
  if(keepZDC){
  t4->Branch("ZDC_n",&ZDC_n);
  t4->Branch("ZDC_e",ZDC_e,"ZDC_e[ZDC_n]/F");
  t4->Branch("ZDC_zside",ZDC_zside,"ZDC_zside[ZDC_n]/I");
  t4->Branch("ZDC_section",ZDC_section,"ZDC_section[ZDC_n]/I");
  t4->Branch("ZDC_channel",ZDC_channel,"ZDC_channel[ZDC_n]/I");
  t4->Branch("ZDC_PM_Total_Energy",&ZDC_PM_Total_Energy );  
  t4->Branch("ZDC_P_Total_Energy",&ZDC_P_Total_Energy );
  t4->Branch("ZDC_P_ECal_Energy",&ZDC_P_ECal_Energy );
  t4->Branch("ZDC_P_HCal_Energy",&ZDC_P_HCal_Energy );          
  t4->Branch("ZDC_M_Total_Energy",&ZDC_M_Total_Energy );
  t4->Branch("ZDC_M_ECal_Energy",&ZDC_M_ECal_Energy );
  t4->Branch("ZDC_M_HCal_Energy",&ZDC_M_HCal_Energy); 
  }
     
  }
  
  
void Analysis_NtupleContent::SetTree_EfficiencyTree(TTree *mytree) { Eff_Tree = mytree; }
  void Analysis_NtupleContent::CreateBranches_EfficiencyTree(const std::vector<std::string> &HLTs,
                                   const std::vector<std::string> &selectorNames){
     
  Eff_Tree->Branch("run", &run);
  Eff_Tree->Branch("event", &event);
  Eff_Tree->Branch("ls", &ls);
  Eff_Tree->Branch("fromFullAOD", &fromFullAOD);
  Eff_Tree->Branch("genWeight", &genWeight);
  Eff_Tree->Branch("BSpot_x", &BSpot_x);
  Eff_Tree->Branch("BSpot_y", &BSpot_y);
  Eff_Tree->Branch("BSpot_z", &BSpot_z);
  Eff_Tree->Branch("pv_x", &pv_x);
  Eff_Tree->Branch("pv_y", &pv_y);
  Eff_Tree->Branch("pv_z", &pv_z);

  
  Eff_Tree->Branch("gentrk1_pt", &gentrk1_pt);
  Eff_Tree->Branch("gentrk1_eta", &gentrk1_eta);
  Eff_Tree->Branch("gentrk1_phi", &gentrk1_phi);
  Eff_Tree->Branch("gentrk1_charge", &gentrk1_charge);
  Eff_Tree->Branch("gentrk1_M", &gentrk1_M);
  Eff_Tree->Branch("gentrk1_vtx_x", &gentrk1_vtx_x);
  Eff_Tree->Branch("gentrk1_vtx_y", &gentrk1_vtx_y);
  Eff_Tree->Branch("gentrk1_vtx_z", &gentrk1_vtx_z);
  Eff_Tree->Branch("gentrk1_pdgId", &gentrk1_pdgId);  
  Eff_Tree->Branch("gentrk2_pt", &gentrk2_pt);
  Eff_Tree->Branch("gentrk2_eta", &gentrk2_eta);
  Eff_Tree->Branch("gentrk2_phi", &gentrk2_phi);
  Eff_Tree->Branch("gentrk2_charge", &gentrk2_charge);
  Eff_Tree->Branch("gentrk2_M", &gentrk2_M);
  Eff_Tree->Branch("gentrk2_vtx_x", &gentrk2_vtx_x);
  Eff_Tree->Branch("gentrk2_vtx_y", &gentrk2_vtx_y);
  Eff_Tree->Branch("gentrk2_vtx_z", &gentrk2_vtx_z);
  Eff_Tree->Branch("gentrk2_pdgId", &gentrk2_pdgId); 
  Eff_Tree->Branch("genFinalMass", &genFinalMass);
  Eff_Tree->Branch("genFinalPt", &genFinalPt);
  Eff_Tree->Branch("genFinalEta", &genFinalEta);
  Eff_Tree->Branch("genFinalPhi", &genFinalPhi);
  Eff_Tree->Branch("genFinalE", &genFinalE);
  Eff_Tree->Branch("genFinalPx", &genFinalPx);
  Eff_Tree->Branch("genFinalPy", &genFinalPy);
  Eff_Tree->Branch("genFinalPz", &genFinalPz);  
  

  Eff_Tree->Branch("gentau1_pt", &gentau1_pt);
  Eff_Tree->Branch("gentau1_eta", &gentau1_eta);
  Eff_Tree->Branch("gentau1_phi", &gentau1_phi);
  Eff_Tree->Branch("gentau1_charge", &gentau1_charge);
  Eff_Tree->Branch("gentau1_vtx_x", &gentau1_vtx_x);
  Eff_Tree->Branch("gentau1_vtx_y", &gentau1_vtx_y);
  Eff_Tree->Branch("gentau1_vtx_z", &gentau1_vtx_z);
  Eff_Tree->Branch("gentau1_pdgId", &gentau1_pdgId);  
  Eff_Tree->Branch("gentau1_decay_el", &gentau1_decay_el);
  Eff_Tree->Branch("gentau1_decay_mu", &gentau1_decay_mu); 
  Eff_Tree->Branch("gentau1_decay_1prong", &gentau1_decay_1prong); 
  Eff_Tree->Branch("gentau1_decay_3prong", &gentau1_decay_3prong); 
  Eff_Tree->Branch("gentau1_decay_other", &gentau1_decay_other); 
  Eff_Tree->Branch("gentau1_decay", &gentau1_decay); 
  Eff_Tree->Branch("gentau1_N_pi0", &gentau1_N_pi0);   
  Eff_Tree->Branch("gentau2_pt", &gentau2_pt);
  Eff_Tree->Branch("gentau2_eta", &gentau2_eta);
  Eff_Tree->Branch("gentau2_phi", &gentau2_phi);
  Eff_Tree->Branch("gentau2_charge", &gentau2_charge);
  Eff_Tree->Branch("gentau2_vtx_x", &gentau2_vtx_x);
  Eff_Tree->Branch("gentau2_vtx_y", &gentau2_vtx_y);
  Eff_Tree->Branch("gentau2_vtx_z", &gentau2_vtx_z);
  Eff_Tree->Branch("gentau2_pdgId", &gentau2_pdgId);
  Eff_Tree->Branch("gentau2_decay_el", &gentau2_decay_el);
  Eff_Tree->Branch("gentau2_decay_mu", &gentau2_decay_mu); 
  Eff_Tree->Branch("gentau2_decay_1prong", &gentau2_decay_1prong); 
  Eff_Tree->Branch("gentau2_decay_3prong", &gentau2_decay_3prong); 
  Eff_Tree->Branch("gentau2_decay_other", &gentau2_decay_other); 
  Eff_Tree->Branch("gentau2_decay", &gentau2_decay);  
  Eff_Tree->Branch("gentau2_N_pi0", &gentau2_N_pi0);   
  Eff_Tree->Branch("genDiTauMass", &genDiTauMass);
  Eff_Tree->Branch("genDiTauPt", &genDiTauPt);
  Eff_Tree->Branch("genDiTauEta", &genDiTauEta);
  Eff_Tree->Branch("genDiTauPhi", &genDiTauPhi);
  Eff_Tree->Branch("genDiTauE", &genDiTauE);
  Eff_Tree->Branch("genDiTauPx", &genDiTauPx);
  Eff_Tree->Branch("genDiTauPy", &genDiTauPy);
  Eff_Tree->Branch("genDiTauPz", &genDiTauPz);
  
  Eff_Tree->Branch("nGen", &nGen );
  Eff_Tree->Branch("gen_pdgId", &gen_pdgId );
  Eff_Tree->Branch("gen_pT", &gen_pT );
  Eff_Tree->Branch("gen_eta", &gen_eta );
  Eff_Tree->Branch("gen_phi", &gen_phi ); 
  Eff_Tree->Branch("gen_charge", &gen_charge );
  Eff_Tree->Branch("gen_mass", &gen_mass );  
  Eff_Tree->Branch("gen_vtx_x", &gen_vtx_x );
  Eff_Tree->Branch("gen_vtx_y", &gen_vtx_y );
  Eff_Tree->Branch("gen_vtx_z", &gen_vtx_z );
  Eff_Tree->Branch("gen_status", &gen_status );
  Eff_Tree->Branch("gen_E", &gen_E );
  Eff_Tree->Branch("gen_Et", &gen_Et );  
  
  
  Eff_Tree->Branch("nGen_tau1", &nGen_tau1 );
  Eff_Tree->Branch("gen_pdgId_tau1", &gen_pdgId_tau1 );
  Eff_Tree->Branch("gen_pT_tau1", &gen_pT_tau1 );
  Eff_Tree->Branch("gen_eta_tau1", &gen_eta_tau1 );
  Eff_Tree->Branch("gen_phi_tau1", &gen_phi_tau1 ); 
  Eff_Tree->Branch("gen_charge_tau1", &gen_charge_tau1 );
  Eff_Tree->Branch("gen_mass_tau1", &gen_mass_tau1 );  
  Eff_Tree->Branch("gen_vtx_x_tau1", &gen_vtx_x_tau1 );
  Eff_Tree->Branch("gen_vtx_y_tau1", &gen_vtx_y_tau1 );
  Eff_Tree->Branch("gen_vtx_z_tau1", &gen_vtx_z_tau1 );
  Eff_Tree->Branch("gen_status_tau1", &gen_status_tau1 );
  Eff_Tree->Branch("gen_E_tau1", &gen_E_tau1 );
  Eff_Tree->Branch("gen_Et_tau1", &gen_Et_tau1 );

  Eff_Tree->Branch("nGen_tau2", &nGen_tau2 );
  Eff_Tree->Branch("gen_pdgId_tau2", &gen_pdgId_tau2 );
  Eff_Tree->Branch("gen_pT_tau2", &gen_pT_tau2 );
  Eff_Tree->Branch("gen_eta_tau2", &gen_eta_tau2 );
  Eff_Tree->Branch("gen_phi_tau2", &gen_phi_tau2 ); 
  Eff_Tree->Branch("gen_charge_tau2", &gen_charge_tau2 );
  Eff_Tree->Branch("gen_mass_tau2", &gen_mass_tau2 );  
  Eff_Tree->Branch("gen_vtx_x_tau2", &gen_vtx_x_tau2 );
  Eff_Tree->Branch("gen_vtx_y_tau2", &gen_vtx_y_tau2 );
  Eff_Tree->Branch("gen_vtx_z_tau2", &gen_vtx_z_tau2 );
  Eff_Tree->Branch("gen_status_tau2", &gen_status_tau2 );
  Eff_Tree->Branch("gen_E_tau2", &gen_E_tau2 );
  Eff_Tree->Branch("gen_Et_tau2", &gen_Et_tau2 );

  // Eff_Tree->Branch("gentrk1_match_dr", &gentrk1_match_dr );
  // Eff_Tree->Branch("gentrk1_match_dphi", &gentrk1_match_dphi );
  // Eff_Tree->Branch("gentrk1_match_diff_eta", &gentrk1_match_diff_eta );
  // Eff_Tree->Branch("gentrk1_match_diff_pt", &gentrk1_match_diff_pt ); 
  // Eff_Tree->Branch("gentrk1_diff_vtx_x", &gentrk1_diff_vtx_x );
  // Eff_Tree->Branch("gentrk1_diff_vtx_y", &gentrk1_diff_vtx_y );
  // Eff_Tree->Branch("gentrk1_diff_vtx_z", & gentrk1_diff_vtx_z);
  // Eff_Tree->Branch("gentrk1_isTag", &gentrk1_isTag );
  
  // Eff_Tree->Branch("gentrk2_match_dr", &gentrk2_match_dr );
  // Eff_Tree->Branch("gentrk2_match_dphi", &gentrk2_match_dphi );
  // Eff_Tree->Branch("gentrk2_match_diff_eta", &gentrk2_match_diff_eta );
  // Eff_Tree->Branch("gentrk2_match_diff_pt", &gentrk2_match_diff_pt ); 
  // Eff_Tree->Branch("gentrk2_diff_vtx_x", &gentrk2_diff_vtx_x );
  // Eff_Tree->Branch("gentrk2_diff_vtx_y", &gentrk2_diff_vtx_y );
  // Eff_Tree->Branch("gentrk2_diff_vtx_z", &gentrk2_diff_vtx_z );  
  // Eff_Tree->Branch("gentrk2_isTag", &gentrk2_isTag );
  
  
  for (unsigned int ihlt = 0; ihlt < HLTs.size(); ihlt++)
    Eff_Tree->Branch(TString(HLTs[ihlt]), &trigger[ihlt]);


 //commented out to save space, added useful details for HIUPC collisions
    if(keepMuons){
    Eff_Tree->Branch("nMu", &nMu);
    Eff_Tree->Branch("mu_isProbe", &mu_isProbe);
    Eff_Tree->Branch("mu_isTag", &mu_isTag);    
    Eff_Tree->Branch("muPt", &muPt);
    Eff_Tree->Branch("muEta", &muEta);
    Eff_Tree->Branch("muPhi", &muPhi);
    Eff_Tree->Branch("muCharge", &muCharge);
    Eff_Tree->Branch("muType", &muType);
    Eff_Tree->Branch("muIsGood", &muIsGood);
    
    Eff_Tree->Branch("muIsGlobal", &muIsGlobal);
    Eff_Tree->Branch("muIsTracker", &muIsTracker);
    Eff_Tree->Branch("muIsPF", &muIsPF);
    Eff_Tree->Branch("muIsSTA", &muIsSTA);

    Eff_Tree->Branch("muD0", &muD0);
    Eff_Tree->Branch("muDz", &muDz);
    Eff_Tree->Branch("muIP3D", &muIP3D);
    Eff_Tree->Branch("muD0Err", &muD0Err);
    Eff_Tree->Branch("muDzErr", &muDzErr);
    Eff_Tree->Branch("muIP3DErr", &muIP3DErr);
    Eff_Tree->Branch("muChi2NDF", &muChi2NDF);
    Eff_Tree->Branch("muInnerD0", &muInnerD0);
    Eff_Tree->Branch("muInnerDz", &muInnerDz);
    
    Eff_Tree->Branch("muInnerD0Err", &muInnerD0Err);
    Eff_Tree->Branch("muInnerDzErr", &muInnerDzErr);
    Eff_Tree->Branch("muInnerPt", &muInnerPt);
    Eff_Tree->Branch("muInnerPtErr", &muInnerPtErr);
    Eff_Tree->Branch("muInnerEta", &muInnerEta);

    Eff_Tree->Branch("muTrkLayers", &muTrkLayers);
    Eff_Tree->Branch("muPixelLayers", &muPixelLayers);
    Eff_Tree->Branch("muPixelHits", &muPixelHits);
    Eff_Tree->Branch("muMuonHits", &muMuonHits);
    Eff_Tree->Branch("muTrkQuality", &muTrkQuality);
    Eff_Tree->Branch("muStations", &muStations);
    Eff_Tree->Branch("muIsoTrk", &muIsoTrk);
    Eff_Tree->Branch("muPFChIso", &muPFChIso);
    Eff_Tree->Branch("muPFPhoIso", &muPFPhoIso);
    Eff_Tree->Branch("muPFNeuIso", &muPFNeuIso);
    Eff_Tree->Branch("muPFPUIso", &muPFPUIso);
    Eff_Tree->Branch("muIDSoft", &muIDSoft);
    Eff_Tree->Branch("muIDLoose", &muIDLoose);
    Eff_Tree->Branch("muIDMedium", &muIDMedium);
    Eff_Tree->Branch("muIDMediumPrompt", &muIDMediumPrompt);
    Eff_Tree->Branch("muIDTight", &muIDTight);
    Eff_Tree->Branch("muIDGlobalHighPt", &muIDGlobalHighPt);
    Eff_Tree->Branch("muIDTrkHighPt", &muIDTrkHighPt);
    Eff_Tree->Branch("muIDInTime", &muIDInTime);

    }
    
    if(keepTracks){
    Eff_Tree->Branch("nTrk", &nTrk);
    Eff_Tree->Branch("trkisTag", &trkisTag);
    Eff_Tree->Branch("trkisProbe", &trkisProbe);
    Eff_Tree->Branch("trkPt", &trkPt);
    Eff_Tree->Branch("trkP", &trkP);
    Eff_Tree->Branch("trkEta", &trkEta);
    Eff_Tree->Branch("trkPhi", &trkPhi);
    Eff_Tree->Branch("trkcharge", &trkcharge);
    Eff_Tree->Branch("trkvx", &trkvx);
    Eff_Tree->Branch("trkvy", &trkvy);
    Eff_Tree->Branch("trkvz", &trkvz);
    Eff_Tree->Branch("trknormchi2", &trknormchi2);
    Eff_Tree->Branch("trkchi2", &trkchi2);
    Eff_Tree->Branch("trkd0", &trkd0);
    Eff_Tree->Branch("trkdxy", &trkdxy);
    Eff_Tree->Branch("trkdz", &trkdz);
    // Eff_Tree->Branch("trkdxyBS", &trkdxyBS);
    // Eff_Tree->Branch("trkdzBS", &trkdzBS);
    Eff_Tree->Branch("trkdxyError", &trkdxyError);
    Eff_Tree->Branch("trkdzError", &trkdzError);
    Eff_Tree->Branch("trkValidHits", &trkValidHits);
    Eff_Tree->Branch("trkMissHits", &trkMissHits);
    Eff_Tree->Branch("trkPurity", &trkPurity);
    }


   if(keepElectrons){
   Eff_Tree->Branch("nEle", &nEle);  
   Eff_Tree->Branch("ele_isProbe", &ele_isProbe);
   Eff_Tree->Branch("ele_isTag", &ele_isTag);
   Eff_Tree->Branch("eleCharge", &eleCharge);
   Eff_Tree->Branch("eleChargeConsistent", &eleChargeConsistent);
   Eff_Tree->Branch("eleSCPixCharge", &eleSCPixCharge);
   // Eff_Tree->Branch("eleCtfCharge", &eleCtfCharge);
   Eff_Tree->Branch("eleEn", &eleEn);
   Eff_Tree->Branch("eleD0", &eleD0);
   Eff_Tree->Branch("eleDz", &eleDz);
   Eff_Tree->Branch("eleIP3D", &eleIP3D);
   Eff_Tree->Branch("eleD0Err", &eleD0Err);
   Eff_Tree->Branch("eleDzErr", &eleDzErr);
   Eff_Tree->Branch("eleIP3DErr", &eleIP3DErr);
   Eff_Tree->Branch("eleTrkPt", &eleTrkPt);
   Eff_Tree->Branch("eleTrkEta", &eleTrkEta);
   Eff_Tree->Branch("eleTrkPhi", &eleTrkPhi);
   Eff_Tree->Branch("eleTrkCharge", &eleTrkCharge);
   Eff_Tree->Branch("eleTrkPtErr", &eleTrkPtErr);
   Eff_Tree->Branch("eleTrkChi2", &eleTrkChi2);
   Eff_Tree->Branch("eleTrkNdof", &eleTrkNdof);
   Eff_Tree->Branch("eleTrkNormalizedChi2", &eleTrkNormalizedChi2);
   Eff_Tree->Branch("eleTrkValidHits", &eleTrkValidHits);
   Eff_Tree->Branch("eleTrkLayers", &eleTrkLayers);

   Eff_Tree->Branch("elePt", &elePt);
   Eff_Tree->Branch("eleEta", &eleEta);
   Eff_Tree->Branch("elePhi", &elePhi);
   Eff_Tree->Branch("eleSCEn", &eleSCEn);
   Eff_Tree->Branch("eleESEn", &eleESEn);
   Eff_Tree->Branch("eleSCEta", &eleSCEta);
   Eff_Tree->Branch("eleSCPhi", &eleSCPhi);
   Eff_Tree->Branch("eleSCRawEn", &eleSCRawEn);
   Eff_Tree->Branch("eleSCEtaWidth", &eleSCEtaWidth);
   Eff_Tree->Branch("eleSCPhiWidth", &eleSCPhiWidth);
   Eff_Tree->Branch("eleHoverE", &eleHoverE);
   Eff_Tree->Branch("eleHoverEBc", &eleHoverEBc);
   Eff_Tree->Branch("eleEoverP", &eleEoverP);
   Eff_Tree->Branch("eleEoverPInv", &eleEoverPInv);
   Eff_Tree->Branch("eleEcalE", &eleEcalE);
   Eff_Tree->Branch("elePAtVtx", &elePAtVtx);
   Eff_Tree->Branch("elePAtSC", &elePAtSC);
   Eff_Tree->Branch("elePAtCluster", &elePAtCluster);
   Eff_Tree->Branch("elePAtSeed", &elePAtSeed);
   Eff_Tree->Branch("eleBrem", &eleBrem);
   Eff_Tree->Branch("eledEtaAtVtx", &eledEtaAtVtx);
   Eff_Tree->Branch("eledPhiAtVtx", &eledPhiAtVtx);
   Eff_Tree->Branch("eledEtaSeedAtVtx", &eledEtaSeedAtVtx);
   Eff_Tree->Branch("eleSigmaIEtaIEta", &eleSigmaIEtaIEta);
   Eff_Tree->Branch("eleSigmaIEtaIEta_2012", &eleSigmaIEtaIEta_2012);
   Eff_Tree->Branch("eleSigmaIPhiIPhi", &eleSigmaIPhiIPhi);
   Eff_Tree->Branch("eleMissHits", &eleMissHits);
   }
   
   // if(keepPFcands){
  // Eff_Tree->Branch("nPFCands", &nPFCands);
  // Eff_Tree->Branch("pfcand_isProbe", &pfcand_isProbe);  
  // Eff_Tree->Branch("pfcand_pdgId", &pfcand_pdgId);
  // Eff_Tree->Branch("pfcand_charge", &pfcand_charge);  
  // Eff_Tree->Branch("pfcand_pt", &pfcand_pt);  
  // Eff_Tree->Branch("pfcand_eta", &pfcand_eta);  
  // Eff_Tree->Branch("pfcand_phi", &pfcand_phi);    
  // Eff_Tree->Branch("pfcand_vtx_x", &pfcand_vtx_x); 
  // Eff_Tree->Branch("pfcand_vtx_y", &pfcand_vtx_y);    
  // Eff_Tree->Branch("pfcand_vtx_z", &pfcand_vtx_z);      
   // }
   
  Eff_Tree->Branch("Pos_pt", &tag_pt);
  Eff_Tree->Branch("Pos_eta", &tag_eta);
  Eff_Tree->Branch("Pos_phi", &tag_phi);
  Eff_Tree->Branch("Pos_charge", &tag_charge);
  Eff_Tree->Branch("Pos_vtx_x", &tag_vtx_x);
  Eff_Tree->Branch("Pos_vtx_y", &tag_vtx_y);
  Eff_Tree->Branch("Pos_vtx_z", &tag_vtx_z);
   
  Eff_Tree->Branch("Neg_pt", &probe_pt);
  Eff_Tree->Branch("Neg_eta", &probe_eta);
  Eff_Tree->Branch("Neg_phi", &probe_phi);
  Eff_Tree->Branch("Neg_charge", &probe_charge);
  Eff_Tree->Branch("Neg_vtx_x", &probe_vtx_x);
  Eff_Tree->Branch("Neg_vtx_y", &probe_vtx_y);
  Eff_Tree->Branch("Neg_vtx_z", &probe_vtx_z);
 
  }
  
  
  
  void Analysis_NtupleContent::CreateExtraTrgBranches_EfficiencyTree(const std::vector<std::string> &HLTs, bool isPos) {
  for (unsigned int ihlt = 0; ihlt < HLTs.size(); ihlt++) {
  if(isPos){
      Eff_Tree->Branch(TString("Pos_" + HLTs[ihlt]), &tag_trg[ihlt]);
      Eff_Tree->Branch(TString("Pos_" + HLTs[ihlt] + "_pt"), &tag_trg_pt[ihlt]);
      Eff_Tree->Branch(TString("Pos_" + HLTs[ihlt] + "_eta"), &tag_trg_eta[ihlt]);
      Eff_Tree->Branch(TString("Pos_" + HLTs[ihlt] + "_phi"), &tag_trg_phi[ihlt]);
      Eff_Tree->Branch(TString("Pos_" + HLTs[ihlt] + "_dr"), &tag_trg_dr[ihlt]);
  }if(!isPos){
      Eff_Tree->Branch(TString("Neg_" + HLTs[ihlt]), &probe_trg[ihlt]);
      Eff_Tree->Branch(TString("Neg_" + HLTs[ihlt] + "_pt"), &probe_trg_pt[ihlt]);
      Eff_Tree->Branch(TString("Neg_" + HLTs[ihlt] + "_eta"), &probe_trg_eta[ihlt]);
      Eff_Tree->Branch(TString("Neg_" + HLTs[ihlt] + "_phi"), &probe_trg_phi[ihlt]);
      Eff_Tree->Branch(TString("Neg_" + HLTs[ihlt] + "_dr"), &probe_trg_dr[ihlt]);
  }
}
  }


void Analysis_NtupleContent::ClearBranches() {
   
  CutThrough_Num = 0;
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
  hasValidVertex = false;
  hasFakeVertex = false;
  trueNumInteractions = -1.0;
  puNumInteractions = -1;
  Rho = -1;
  ntag = 0;
  ntag_electron = 0;
  npairs = 0;

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
  gentrk1_pt = 0;
  gentrk1_eta = -99;
  gentrk1_phi = -99;
  gentrk1_charge = 0;
  gentrk1_M = 0;
  gentrk1_vtx_x = -99;
  gentrk1_vtx_y = -99;
  gentrk1_vtx_z = -99;
  gentrk1_px = -99;
  gentrk1_py = -99;
  gentrk1_pz = -99;
  gentrk1_pdgId = -99;
  gentrk2_pt = 0;
  gentrk2_eta = -99;
  gentrk2_phi = -99;
  gentrk2_charge = 0;
  gentrk2_M = 0;
  gentrk2_vtx_x = -99;
  gentrk2_vtx_y = -99;
  gentrk2_vtx_z = -99;
  gentrk2_px = -99;
  gentrk2_py = -99;
  gentrk2_pz = -99;
  gentrk2_pdgId = -99;
  genFinalMass = -99;
  genFinalPt = -99;
  genFinalEta = -99;
  genFinalPhi = -99;
  genFinalE = -99;
  genFinalPx = -99;
  genFinalPy= -99;
  genFinalPz = -99; 
  
  
  indep_pt = -99;
  indep_eta = -99;
  indep_phi = -99;
  indep_mass = -99;
  gentau1_gamma_pz = -99;
  gentau2_gamma_pz = -99;  
  gentau1_pt = 0;
  gentau1_eta = -99;
  gentau1_phi = -99;
  gentau1_charge = 0;
  gentau1_vtx_x = -99;
  gentau1_vtx_y = -99;
  gentau1_vtx_z = -99;
  gentau1_px = -99;
  gentau1_py = -99;
  gentau1_pz = -99;
  gentau1_pdgId = -99;
  gentau1_decay_el = false;
  gentau1_decay_mu = false;
  gentau1_decay_1prong = false;  
  gentau1_decay_3prong = false;
  gentau1_decay_other = false;
  gentau1_decay = -99;
  gentau1_N_pi0 = 0;
  gentau2_pt = 0;
  gentau2_eta = -99;
  gentau2_phi = -99;
  gentau2_charge = 0;
  gentau2_vtx_x = -99;
  gentau2_vtx_y = -99;
  gentau2_vtx_z = -99;
  gentau2_px = -99;
  gentau2_py = -99;
  gentau2_pz = -99;
  gentau2_pdgId = -99;
  gentau2_decay_el = false;
  gentau2_decay_mu = false;
  gentau2_decay_1prong = false;  
  gentau2_decay_3prong = false;
  gentau2_decay_other = false;
  gentau2_decay = -99;
  gentau2_N_pi0 =0;
  genDiTauMass = -99;
  genDiTauPt = -99;
  genDiTauEta = -99;
  genDiTauPhi = -99;
  genDiTauE = -99;
  genDiTauPx = -99;
  genDiTauPy= -99;
  genDiTauPz = -99;   
  
  // genDiTau_pair_svprob = -1;
  // genDiTau_pair_fit_mass = -1 ;
  // genDiTau_pair_normalchi2 = -1;   
  genDiTau_vtx_x = -99;
  genDiTau_vtx_y = -99;
  genDiTau_vtx_z = -99;  
 
  // genFinal_pair_svprob = -1;
  // genFinal_pair_fit_mass = -1 ;
  // genFinal_pair_normalchi2 = -1;    
  genFinal_vtx_x = -99;
  genFinal_vtx_y = -99;
  genFinal_vtx_z = -99; 
   
  genDiTau_vtx_gentau1_IP = -99;
  genDiTau_vtx_gentau2_IP = -99;
  genDiTau_vtx_gentrk1_IP = -99;
  genDiTau_vtx_gentrk2_IP = -99;
  
  genFinal_vtx_gentau1_IP = -99;
  genFinal_vtx_gentau2_IP = -99;
  genFinal_vtx_gentrk1_IP = -99;
  genFinal_vtx_gentrk2_IP = -99;

  // genDiTau_refit_vtx_gentau1_IP = -99;
  // genDiTau_refit_vtx_gentau2_IP = -99;
  // genDiTau_refit_vtx_gentrk1_IP = -99;
  // genDiTau_refit_vtx_gentrk2_IP = -99;
  
  // genFinal_refit_vtx_gentau1_IP = -99;
  // genFinal_refit_vtx_gentau2_IP = -99;
  // genFinal_refit_vtx_gentrk1_IP = -99;
  // genFinal_refit_vtx_gentrk2_IP = -99;   


//new editions
   nGen = 0;
   gen_pdgId.clear();
   gen_pT.clear();
   gen_eta.clear();
   gen_phi.clear(); 
   gen_charge.clear();
   gen_mass.clear();  
   gen_vtx_x.clear();
   gen_vtx_y.clear();
   gen_vtx_z.clear();
   gen_status.clear();
   gen_E.clear();
   gen_Et.clear();
   
   nGen_tau1 = 0;
   gen_pdgId_tau1.clear();
   gen_pT_tau1.clear();
   gen_eta_tau1.clear();
   gen_phi_tau1.clear(); 
   gen_charge_tau1.clear();
   gen_mass_tau1.clear();  
   gen_vtx_x_tau1.clear();
   gen_vtx_y_tau1.clear();
   gen_vtx_z_tau1.clear();
   gen_status_tau1.clear();
   gen_E_tau1.clear();
   gen_Et_tau1.clear();
   
   nGen_tau2 = 0;
   gen_pdgId_tau2.clear();
   gen_pT_tau2.clear();
   gen_eta_tau2.clear();
   gen_phi_tau2.clear(); 
   gen_charge_tau2.clear();
   gen_mass_tau2.clear();  
   gen_vtx_x_tau2.clear();
   gen_vtx_y_tau2.clear();
   gen_vtx_z_tau2.clear();
   gen_status_tau2.clear();
   gen_E_tau2.clear();
   gen_Et_tau2.clear();
   
  gentrk1_match_dr = -99;
  gentrk1_match_dphi = -99;
  gentrk1_match_diff_eta= -99;
  gentrk1_match_diff_pt = -99; 
  gentrk1_diff_vtx_x = -99;
  gentrk1_diff_vtx_y = -99;
  gentrk1_diff_vtx_z = -99;  
  gentrk1_isTag = false;
  
  gentrk2_match_dr = -99;
  gentrk2_match_dphi = -99;
  gentrk2_match_diff_eta= -99;
  gentrk2_match_diff_pt = -99; 
  gentrk2_diff_vtx_x = -99;
  gentrk2_diff_vtx_y = -99;
  gentrk2_diff_vtx_z = -99;  
  gentrk2_isTag = false;


  trg_filter.clear();
  trg_pt.clear();
  trg_eta.clear();
  trg_phi.clear();
  prb_filter.clear();
  prb_pt.clear();
  prb_eta.clear();
  prb_phi.clear();

  tag_isMuon = false;
  tag_isElectron = false;
  tag_pt = 0;
  tag_eta = -99;
  tag_phi = -99;
  tag_charge = -99;
  tag_pterr = 0;
  tag_dxy = -99;
  tag_dz = -99;
  tag_vtx_x = -99;
  tag_vtx_y = -99;
  tag_vtx_z = -99;
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
  tag_pfIso04_charged = -99;
  tag_pfIso04_neutral = -99;
  tag_pfIso04_photon = -99;
  tag_pfIso04_sumPU = -99;
  tag_tuneP_pt = -99;
  tag_tuneP_pterr = -99;
  tag_nsegments = -99;
  tag_hasTrackMatch = false;
  tag_TrackMatchDR = 99;

  iprobe = 0;
  probe_pt = 0;
  probe_eta = -99;
  probe_phi = -99;
  probe_charge = -99;
  probe_vtx_x = -99;
  probe_vtx_y = -99;
  probe_vtx_z = -99;
  // probe_isLoose = false;
  // probe_isMedium = false;
  // probe_isTight = false;
  // probe_isSoft = false;
  // probe_isHighPt = false;
  // probe_isArbitratedTracker = false;
  // probe_isMuMatched = false;
  // probe_isPF = false;
  // probe_isSA = false;
  // probe_isTracker = false;
  // probe_isGlobal = false;
  // probe_isdSA = false;
  // probe_isdGlobal = false;
  probe_isCosmic = false;
  probe_ncosmic = -99;
  probe_cosmic_minDR = +99;
  // probe_isGood = false;
  probe_isHighPurity = false;
  // probe_validFraction = -99;
  // probe_trkChi2 = -99;
  // probe_positionChi2 = -99;
  // probe_trkKink = -99;
  // probe_segmentCompatibility = -99;
  // probe_trackerLayers = -99;
  // probe_pixelLayers = -99;
  // probe_muonStations = -99;
  // probe_muonHits = -99;
  // probe_DTHits = -99;
  // probe_CSCHits = -99;
  // probe_pterr = 0;
  // probe_dxy = -99;
  // probe_dz = -99;
  probe_relIso04 = -99;
  probe_miniIso = -1.;
  probe_miniIsoCharged = 0.;
  probe_miniIsoPhotons = 0.;
  probe_miniIsoNeutrals = 0.;
  probe_isMatchedGen = false;
  // probe_minDR = 0.;
  // probe_ptRel_minDR = 0.;
  // probe_iso03_sumPt = -99;
  // probe_pfIso04_charged = -99;
  // probe_pfIso04_neutral = -99;
  // probe_pfIso04_photon = -99;
  // probe_pfIso04_sumPU = -99;
  // probe_pixelHits = -99;
  // probe_matchedStations = -99;
  // probe_expectedMatchedStations = -99;
  // probe_RPCLayers = -99;
  // probe_stationMask = 0;
  // probe_nShowers = -99;
  // probe_tuneP_pt = -99;
  // probe_tuneP_pterr = -99;
  // probe_tuneP_muonHits = -99;
  // probe_nsegments = -99;
  
  probe_hasMuonMatch = false;
  probe_MuonMatchDR = 99;
  probe_hasElectronMatch = false;
  probe_ElectronMatchDR = 99;
  probe_hasPFCMatch = false ;
  probe_PFCMatchDR = 99;
  probe_pfcID = -9999;

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


  pair_pt = 0;
  pair_mass = 0;
  pair_eta = -99;
  pair_phi = -99;
  pair_fit_mass = 0;
  pair_svprob = 0;
  pair_normalchi2 = 0;
  pair_dz = -99;
  pair_dR = -99;
  pair_dphi = -99;
  pair_rank_vtx_prob = -1;
  pair_rank_dPhi_muons = -1;
  pair_rank_Mass_Mmumu = -1;

  
  
  
//new editions
  refit_vtx_x = -99;
  refit_vtx_y = -99;
  refit_vtx_z = -99;
  tag_refit_transverse_IP = -99;
  probe_refit_transverse_IP = -99;
  tag_transverse_IP = -99;
  probe_transverse_IP = -99;
  
  pair_tuneP_pt = -99;
  pair_tuneP_mass = -99;
  pair_tuneP_eta = -99;
  pair_tuneP_phi = -99;
  pair_tuneP_fit_mass = -99;
  pair_tuneP_svprob = -99;
  pair_tuneP_normalchi2 = -99;
  pair_tuneP_dz = -99;
  pair_tuneP_dR = -99;
  pair_tuneP_dphi = -99;

  refit_tuneP_vtx_x = -99;
  refit_tuneP_vtx_y = -99;
  refit_tuneP_vtx_z = -99;
  tag_refit_tuneP_transverse_IP = -99;
  probe_refit_tuneP_transverse_IP = -99;
  tag_tuneP_transverse_IP = -99;
  probe_tuneP_transverse_IP = -99;
  
    nMu = 0;
    mu_isProbe.clear();
    mu_isTag.clear();
    muPt.clear();
    muEta.clear();
    muPhi.clear();
    muCharge.clear();
    muType.clear();
    muIsGood.clear();

    muIsGlobal.clear();
    muIsTracker.clear();
    muIsPF.clear();
    muIsSTA.clear();

    muD0.clear();
    muDz.clear();
    muIP3D.clear();
    muD0Err.clear();
    muDzErr.clear();
    muIP3DErr.clear();
    muChi2NDF.clear();
    muInnerD0.clear();
    muInnerDz.clear();
    
    muInnerD0Err.clear();
    muInnerDzErr.clear();
    muInnerPt.clear();
    muInnerPtErr.clear();
    muInnerEta.clear();

    muTrkLayers.clear();
    muPixelLayers.clear();
    muPixelHits.clear();
    muMuonHits.clear();
    muTrkQuality.clear();
    muStations.clear();
    muIsoTrk.clear();
    muPFChIso.clear();
    muPFPhoIso.clear();
    muPFNeuIso.clear();
    muPFPUIso.clear();
    muIDSoft.clear();
    muIDLoose.clear();
    muIDMedium.clear();
    muIDMediumPrompt.clear();
    muIDTight.clear();
    muIDGlobalHighPt.clear();
    muIDTrkHighPt.clear();
    muIDInTime.clear();

    nTrk =0;
    trkisProbe.clear();
    trkisTag.clear();
    trkPt.clear();
    trkP.clear();
    trkEta.clear();
    trkPhi.clear();
    trkcharge.clear();
    trkvx.clear();
    trkvy.clear(); 
    trkvz.clear();
    trknormchi2.clear();                 
    trkchi2.clear();
    trkd0.clear(); 
    trkdxy.clear();
    trkdz.clear();
    // trkdxyBS.clear();
    // trkdzBS.clear();
    trkdxyError.clear();     
    trkdzError.clear();
    trkValidHits.clear();                     
    trkMissHits.clear();
    trkPurity.clear();

  
    nEle = 0;
    ele_isProbe.clear();
    ele_isTag.clear();
    eleCharge.clear();
    eleChargeConsistent.clear();
    eleSCPixCharge.clear();
    // eleCtfCharge.clear();
    eleEn.clear();
    eleD0.clear();
    eleDz.clear();
    eleIP3D.clear();
    eleD0Err.clear();
    eleDzErr.clear();
    eleIP3DErr.clear();
    eleTrkPt.clear();
    eleTrkEta.clear();
    eleTrkPhi.clear();
    eleTrkCharge.clear();
    eleTrkPtErr.clear();
    eleTrkChi2.clear();
    eleTrkNdof.clear();
    eleTrkNormalizedChi2.clear();
    eleTrkValidHits.clear();
    eleTrkLayers.clear();
    elePt.clear();
    eleEta.clear();
    elePhi.clear();
    eleSCEn.clear();
    eleESEn.clear();
    eleSCEta.clear();
    eleSCPhi.clear();
    eleSCRawEn.clear();
    eleSCEtaWidth.clear();
    eleSCPhiWidth.clear();
    eleHoverE.clear();
    eleHoverEBc.clear();
    eleEoverP.clear();
    eleEoverPInv.clear();
    eleEcalE.clear();
    elePAtVtx.clear();
    elePAtSC.clear();
    elePAtCluster.clear();
    elePAtSeed.clear();
    eleBrem.clear();
    eledEtaAtVtx.clear();
    eledPhiAtVtx.clear();
    eledEtaSeedAtVtx.clear();
    eleSigmaIEtaIEta.clear();
    eleSigmaIEtaIEta_2012.clear();
    eleSigmaIPhiIPhi.clear();
    eleMissHits.clear();
  
  nPFCands = 0;
  pfcand_isProbe.clear();
  pfcand_pdgId.clear();
  pfcand_charge.clear(); 
  pfcand_pt.clear();
  pfcand_eta.clear();
  pfcand_phi.clear(); 
  pfcand_vtx_x.clear(); 
  pfcand_vtx_y.clear(); 
  pfcand_vtx_z.clear();   
  
  //photons
  
    nPho = 0;
    phoE.clear();
    phoEt.clear();
    phoEta.clear();
    phoPhi.clear();

    phoEcorrStdEcal.clear();
    phoEcorrPhoEcal.clear();
    phoEcorrRegr1.clear();
    phoEcorrRegr2.clear();
    phoEcorrErrStdEcal.clear();
    phoEcorrErrPhoEcal.clear();
    phoEcorrErrRegr1.clear();
    phoEcorrErrRegr2.clear();

    phoSCE.clear();
    phoSCEt.clear();
    phoSCRawE.clear();
    phoSCEta.clear();
    phoSCPhi.clear();
    phoSCEtaWidth.clear();
    phoSCPhiWidth.clear();
    phoSCBrem.clear();
    phoSCnHits.clear();
    phoSCflags.clear();
    phoSCinClean.clear();
    phoSCinUnClean.clear();
    phoSCnBC.clear();
    phoESEn.clear();
    
    phoIsPFPhoton.clear();
    phoIsStandardPhoton.clear();
    phoHasPixelSeed.clear();
    phoHasConversionTracks.clear();
    // phoEleVeto.clear();  // TODO: not available in reco::
    phoHadTowerOverEm.clear();
    phoHoverE.clear();
    phoHoverEValid.clear();
    phoSigmaIEtaIEta.clear();
    phoR9.clear();
  
  nTower= 0;
  maxHFp =0;
  maxHFm =0;
  CaloTower_hadE.clear();
  CaloTower_emE.clear();
  CaloTower_e.clear();
  CaloTower_et.clear();
  CaloTower_eta.clear();
  CaloTower_phi.clear();
  CaloTower_HF_AboveThreshold.clear();
  CaloTower_Had_AboveThreshold.clear();
  CaloTower_EM_AboveThreshold.clear();
  CaloTower_HF_inNoiseRegion.clear();
  CaloTower_Had_inNoiseRegion.clear();
  CaloTower_EM_inNoiseRegion.clear();
  
  ZDC_n = 0;
  for (unsigned int i = 0; i < 18; i++) {
  ZDC_e[i] = 0;
  ZDC_zside[i] = -99;
  ZDC_section [i]= -99;
  ZDC_channel[i] = -99;
  ZDC_saturation[i] = -99;
  }
  ZDC_PM_Total_Energy = 0;  
  ZDC_P_Total_Energy = 0;
  ZDC_P_ECal_Energy = 0;
  ZDC_P_HCal_Energy = 0;          
  ZDC_M_Total_Energy = 0;
  ZDC_M_ECal_Energy = 0;
  ZDC_M_HCal_Energy = 0; 
  
  
  // probe info
  
  probe_isMuon = false ;
  probe_isElectron = false ;
  probe_isPion = false ;
  probe_isOther = false ;
  probe_typeSum = 0;
  
  
   // 3 prong things
  
  
  iprobe_3prong = 0;
  npairs_3prong = 0;
  
  tau_3prong_trk1_pt = 0;
  tau_3prong_trk1_eta = 0;
  tau_3prong_trk1_phi = 0;
  tau_3prong_trk1_charge = 0;
  tau_3prong_trk1_vtx_x = 0;
  tau_3prong_trk1_vtx_y = 0;
  tau_3prong_trk1_vtx_z = 0;
        
  tau_3prong_trk2_pt = 0;
  tau_3prong_trk2_eta = 0;
  tau_3prong_trk2_phi = 0;
   tau_3prong_trk2_charge = 0;
  tau_3prong_trk2_vtx_x = 0;
  tau_3prong_trk2_vtx_y = 0;
  tau_3prong_trk2_vtx_z = 0;
     
  tau_3prong_trk3_pt = 0;
  tau_3prong_trk3_eta = 0;
  tau_3prong_trk3_phi = 0;
   tau_3prong_trk3_charge = 0;
  tau_3prong_trk3_vtx_x = 0;
  tau_3prong_trk3_vtx_y = 0;
  tau_3prong_trk3_vtx_z  = 0;
        
  tau_3prong_total_pt = 0;
  tau_3prong_total_eta= 0;
  tau_3prong_total_phi = 0;
  tau_3prong_total_M = 0;
  tau_3prong_total_pt_sum = 0;
   tau_3prong_total_charge = 0;
   tau_3prong_total_Dz = 0;
  tau_3prong_total_vtx_x = 0;
  tau_3prong_total_vtx_y = 0;
  tau_3prong_total_vtx_z = 0;
  tau_3prong_total_vtx_prob = 0;
  tau_3prong_total_vtx_chi2 = 0;
  tau_3prong_DR_1_2 = 0;
  tau_3prong_DR_1_3 = 0;
  tau_3prong_DR_2_3 = 0;
  tau_3prong_DR_tag_1 = 0;
  tau_3prong_DR_tag_3 = 0;
  tau_3prong_DR_tag_3 = 0;
  tau_3prong_total_area = 0;
  
  
  pfcand_isFirstProbe.clear();
  pfcand_isSecondProbe.clear();
  pfcand_isThirdProbe.clear();
  trk_isFirstProbe.clear();
  trk_isSecondProbe.clear();
  trk_isThirdProbe.clear();
 
  pair_rank_vtx_prob_3prong = -1;
  pair_rank_dPhi_muons_3prong = -1;
  pair_rank_Mass_Mmumu_3prong = -1;
  

}

// Clears Vectors that are Filled both the 3prong and 1 prong loops to avoid double filling
void Analysis_NtupleContent::ClearVectors(){
   
       nMu = 0;
    mu_isProbe.clear();
    mu_isTag.clear();
    muPt.clear();
    muEta.clear();
    muPhi.clear();
    muCharge.clear();
    muType.clear();
    muIsGood.clear();

    muIsGlobal.clear();
    muIsTracker.clear();
    muIsPF.clear();
    muIsSTA.clear();

    muD0.clear();
    muDz.clear();
    muIP3D.clear();
    muD0Err.clear();
    muDzErr.clear();
    muIP3DErr.clear();
    muChi2NDF.clear();
    muInnerD0.clear();
    muInnerDz.clear();
    
    muInnerD0Err.clear();
    muInnerDzErr.clear();
    muInnerPt.clear();
    muInnerPtErr.clear();
    muInnerEta.clear();

    muTrkLayers.clear();
    muPixelLayers.clear();
    muPixelHits.clear();
    muMuonHits.clear();
    muTrkQuality.clear();
    muStations.clear();
    muIsoTrk.clear();
    muPFChIso.clear();
    muPFPhoIso.clear();
    muPFNeuIso.clear();
    muPFPUIso.clear();
    muIDSoft.clear();
    muIDLoose.clear();
    muIDMedium.clear();
    muIDMediumPrompt.clear();
    muIDTight.clear();
    muIDGlobalHighPt.clear();
    muIDTrkHighPt.clear();
    muIDInTime.clear();

    nTrk =0;
    trkisProbe.clear();
    trkisTag.clear();
    trkPt.clear();
    trkP.clear();
    trkEta.clear();
    trkPhi.clear();
    trkcharge.clear();
    trkvx.clear();
    trkvy.clear(); 
    trkvz.clear();
    trknormchi2.clear();                 
    trkchi2.clear();
    trkd0.clear(); 
    trkdxy.clear();
    trkdz.clear();
    // trkdxyBS.clear();
    // trkdzBS.clear();
    trkdxyError.clear();     
    trkdzError.clear();
    trkValidHits.clear();                     
    trkMissHits.clear();
    trkPurity.clear();

  
    nEle = 0;
    ele_isProbe.clear();
    ele_isTag.clear();
    eleCharge.clear();
    eleChargeConsistent.clear();
    eleSCPixCharge.clear();
    // eleCtfCharge.clear();
    eleEn.clear();
    eleD0.clear();
    eleDz.clear();
    eleIP3D.clear();
    eleD0Err.clear();
    eleDzErr.clear();
    eleIP3DErr.clear();
    eleTrkPt.clear();
    eleTrkEta.clear();
    eleTrkPhi.clear();
    eleTrkCharge.clear();
    eleTrkPtErr.clear();
    eleTrkChi2.clear();
    eleTrkNdof.clear();
    eleTrkNormalizedChi2.clear();
    eleTrkValidHits.clear();
    eleTrkLayers.clear();
    elePt.clear();
    eleEta.clear();
    elePhi.clear();
    eleSCEn.clear();
    eleESEn.clear();
    eleSCEta.clear();
    eleSCPhi.clear();
    eleSCRawEn.clear();
    eleSCEtaWidth.clear();
    eleSCPhiWidth.clear();
    eleHoverE.clear();
    eleHoverEBc.clear();
    eleEoverP.clear();
    eleEoverPInv.clear();
    eleEcalE.clear();
    elePAtVtx.clear();
    elePAtSC.clear();
    elePAtCluster.clear();
    elePAtSeed.clear();
    eleBrem.clear();
    eledEtaAtVtx.clear();
    eledPhiAtVtx.clear();
    eledEtaSeedAtVtx.clear();
    eleSigmaIEtaIEta.clear();
    eleSigmaIEtaIEta_2012.clear();
    eleSigmaIPhiIPhi.clear();
    eleMissHits.clear();
  
  nPFCands = 0;
  pfcand_isProbe.clear();
  pfcand_pdgId.clear();
  pfcand_charge.clear(); 
  pfcand_pt.clear();
  pfcand_eta.clear();
  pfcand_phi.clear(); 
  pfcand_vtx_x.clear(); 
  pfcand_vtx_y.clear(); 
  pfcand_vtx_z.clear();   
  
  pfcand_isFirstProbe.clear();
  pfcand_isSecondProbe.clear();
  pfcand_isThirdProbe.clear();
  trk_isFirstProbe.clear();
  trk_isSecondProbe.clear();
  trk_isThirdProbe.clear();
  
}
