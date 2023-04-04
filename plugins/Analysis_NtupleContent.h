//
// Original Author:
//                george karathanasis
//         Created:  Thu, 23 Mar 2019 17:40:23 GMT
//
// flat tree branches/var declaration

#ifndef ANALYSIS_NTUPLECONTENT_H
#define ANALYSIS_NTUPLECONTENT_H
#include <string>
#include <vector>
#include "TString.h"
#include "TTree.h"

class Analysis_NtupleContent {
public:
  Analysis_NtupleContent();
  virtual ~Analysis_NtupleContent();
  void SetTreeVariables(bool, bool,bool, bool, bool, bool, bool);
  void SetTree(TTree *t1); 
  void CreateBranches(const std::vector<std::string> &, const std::vector<std::string> &);
  void SetTree_GenVtxStudy(TTree *t2);
  void SetTree_Test(TTree *t3);  
  void CreateBranches_GenVtxStudy(); 
  void CreateExtraTrgBranches(const std::vector<std::string> &, bool);
  void CreateExtraTrgBranches_3prong(const std::vector<std::string> &, bool);
  void SetTree_3ProngStudy(TTree *t4);
  void CreateBranches_3ProngStudy(const std::vector<std::string> &, const std::vector<std::string> &);
  void SetTree_EfficiencyTree(TTree *Eff_Tree);
  void CreateBranches_EfficiencyTree(const std::vector<std::string> &, const std::vector<std::string> &);
  void CreateExtraTrgBranches_EfficiencyTree(const std::vector<std::string> &, bool);  
  void ClearBranches();
  void ClearVectors();



  // Standard stuff
  int CutThrough_Num;
  int run;
  int event;
  int ls;
  bool fromFullAOD;

  float genWeight;

  // Beamspot and vertex
  float BSpot_x;
  float BSpot_y;
  float BSpot_z;
  float pv_x;
  float pv_y;
  float pv_z;
  int nvertices;
  bool hasValidVertex;
  bool hasFakeVertex;
  

   int nGen;
   std::vector<int> gen_pdgId;
   std::vector<float> gen_pT;
   std::vector<float> gen_eta;
   std::vector<float> gen_phi; 
   std::vector<int> gen_charge;
   std::vector<float> gen_mass;  
   std::vector<float> gen_vtx_x;
   std::vector<float> gen_vtx_y;
   std::vector<float> gen_vtx_z;
   std::vector<int> gen_status;
   std::vector<float> gen_E;
   std::vector<float> gen_Et;
   
   int nGen_tau1;
   std::vector<int> gen_pdgId_tau1;
   std::vector<float> gen_pT_tau1;
   std::vector<float> gen_eta_tau1;
   std::vector<float> gen_phi_tau1; 
   std::vector<int> gen_charge_tau1;
   std::vector<float> gen_mass_tau1;  
   std::vector<float> gen_vtx_x_tau1;
   std::vector<float> gen_vtx_y_tau1;
   std::vector<float> gen_vtx_z_tau1;
   std::vector<int> gen_status_tau1;
   std::vector<float> gen_E_tau1;
   std::vector<float> gen_Et_tau1;
   
   
   int nGen_tau2;
   std::vector<int> gen_pdgId_tau2;
   std::vector<float> gen_pT_tau2;
   std::vector<float> gen_eta_tau2;
   std::vector<float> gen_phi_tau2; 
   std::vector<int> gen_charge_tau2;
   std::vector<float> gen_mass_tau2;  
   std::vector<float> gen_vtx_x_tau2;
   std::vector<float> gen_vtx_y_tau2;
   std::vector<float> gen_vtx_z_tau2;
   std::vector<int> gen_status_tau2;
   std::vector<float> gen_E_tau2;
   std::vector<float> gen_Et_tau2;

   
  // Pileup
  float trueNumInteractions;
  int puNumInteractions;
  double Rho;

  // Number of muons

  int ntag;
  int ntag_electron;
  int iprobe;
  int npairs;

  // Triggers
  static const int NTRIGGERMAX = 100;
  bool trigger[NTRIGGERMAX];
  std::vector<TString> trg_filter;
  std::vector<float> trg_pt;
  std::vector<float> trg_eta;
  std::vector<float> trg_phi;

  std::vector<TString> prb_filter;
  std::vector<float> prb_pt;
  std::vector<float> prb_eta;
  std::vector<float> prb_phi;

  // Trigger matches
  bool tag_trg[NTRIGGERMAX];
  float tag_trg_pt[NTRIGGERMAX];
  float tag_trg_eta[NTRIGGERMAX];
  float tag_trg_phi[NTRIGGERMAX];
  float tag_trg_dr[NTRIGGERMAX];
  bool probe_trg[NTRIGGERMAX];
  float probe_trg_pt[NTRIGGERMAX];
  float probe_trg_eta[NTRIGGERMAX];
  float probe_trg_phi[NTRIGGERMAX];
  float probe_trg_dr[NTRIGGERMAX];

  // Standard selectors in reco::muon::Selector
  bool probe_selectors[100];

  float l1pt;
  int l1q;
  float l1dr;
  float l1ptByQ;
  int l1qByQ;
  float l1drByQ;

  float tag_l1pt;
  int tag_l1q;
  float tag_l1dr;
  float tag_l1ptByQ;
  int tag_l1qByQ;
  float tag_l1drByQ;

  // gentrk1: mu-, gentrk2: mu+
  float gentrk1_pt;
  float gentrk1_eta;
  float gentrk1_phi;
  float gentrk1_charge;
  float gentrk1_M;
  float gentrk1_vtx_x;
  float gentrk1_vtx_y;
  float gentrk1_vtx_z;
  float gentrk1_px;
  float gentrk1_py;
  float gentrk1_pz;
  int gentrk1_pdgId;
  float gentrk2_pt;
  float gentrk2_eta;
  float gentrk2_phi;
  float gentrk2_charge;
  float gentrk2_M;
  float gentrk2_vtx_x;
  float gentrk2_vtx_y;
  float gentrk2_vtx_z;
  float gentrk2_px;
  float gentrk2_py;
  float gentrk2_pz;
  int gentrk2_pdgId;
  float genFinalMass;
  float genFinalPt;
  float genFinalEta;
  float genFinalPhi;
  float genFinalE;
  float genFinalPx;
  float genFinalPy;
  float genFinalPz; 
  
  float gentrk1_match_dr;
  float gentrk1_match_dphi;
  float gentrk1_match_diff_eta;
  float gentrk1_match_diff_pt; 
  float gentrk1_diff_vtx_x;
  float gentrk1_diff_vtx_y;
  float gentrk1_diff_vtx_z;  
  bool gentrk1_isTag;
  
  float gentrk2_match_dr;
  float gentrk2_match_dphi;
  float gentrk2_match_diff_eta;
  float gentrk2_match_diff_pt; 
  float gentrk2_diff_vtx_x;
  float gentrk2_diff_vtx_y;
  float gentrk2_diff_vtx_z;  
  bool gentrk2_isTag;
  
  
  //NEW
  float indep_pt;
  float indep_eta;
  float indep_phi;
  float indep_mass;
  float gentau1_gamma_pz;
  float gentau2_gamma_pz;
  float gentau1_pt;
  float gentau1_eta;
  float gentau1_phi;
  int gentau1_charge;
  float gentau1_vtx_x;
  float gentau1_vtx_y;
  float gentau1_vtx_z;
  float gentau1_px;
  float gentau1_py;
  float gentau1_pz;
  int gentau1_pdgId;
  bool gentau1_decay_el;
  bool gentau1_decay_mu;
  bool gentau1_decay_1prong;  
  bool gentau1_decay_3prong;
  bool gentau1_decay_other;
  int gentau1_decay;
  int gentau1_N_pi0;
  float gentau2_pt;
  float gentau2_eta;
  float gentau2_phi;
  int gentau2_charge;
  float gentau2_vtx_x;
  float gentau2_vtx_y;
  float gentau2_vtx_z;
  float gentau2_px;
  float gentau2_py;
  float gentau2_pz;
  int gentau2_pdgId;
  bool gentau2_decay_el;
  bool gentau2_decay_mu;
  bool gentau2_decay_1prong;  
  bool gentau2_decay_3prong;
  bool gentau2_decay_other;
  int gentau2_decay;
  int gentau2_N_pi0;
  float genDiTauMass;
  float genDiTauPt;
  float genDiTauEta;
  float genDiTauPhi;
  float genDiTauE;
  float genDiTauPx;
  float genDiTauPy;
  float genDiTauPz; 
  
  
  // float genDiTau_pair_svprob;
  // float genDiTau_pair_fit_mass;
  // float genDiTau_pair_normalchi2;  
  float genDiTau_vtx_x;
  float genDiTau_vtx_y;
  float genDiTau_vtx_z;  
 
  // float genFinal_pair_svprob;
  // float genFinal_pair_fit_mass;
  // float genFinal_pair_normalchi2;  
  float genFinal_vtx_x;
  float genFinal_vtx_y;
  float genFinal_vtx_z; 
   
  float genDiTau_vtx_gentau1_IP;
  float genDiTau_vtx_gentau2_IP;
  float genDiTau_vtx_gentrk1_IP;
  float genDiTau_vtx_gentrk2_IP;
  
  float genFinal_vtx_gentau1_IP;
  float genFinal_vtx_gentau2_IP;
  float genFinal_vtx_gentrk1_IP;
  float genFinal_vtx_gentrk2_IP;

  // float genDiTau_refit_vtx_gentau1_IP;
  // float genDiTau_refit_vtx_gentau2_IP;
  // float genDiTau_refit_vtx_gentrk1_IP;
  // float genDiTau_refit_vtx_gentrk2_IP;
  
  // float genFinal_refit_vtx_gentau1_IP;
  // float genFinal_refit_vtx_gentau2_IP;
  // float genFinal_refit_vtx_gentrk1_IP;
  // float genFinal_refit_vtx_gentrk2_IP;   
   
  // Tag properties
  bool tag_isMuon;
  bool tag_isElectron;
  float tag_pt;
  float tag_eta;
  float tag_phi;
  int tag_charge;
  float tag_pterr;
  float tag_dxy;
  float tag_dz;
  float tag_vtx_x;
  float tag_vtx_y;
  float tag_vtx_z;
  bool tag_isPF;
  bool tag_isSA;
  bool tag_isdSA;
  bool tag_isTracker;
  bool tag_isGlobal;
  bool tag_isLoose;
  bool tag_isMedium;
  bool tag_isTight;
  bool tag_isSoft;
  bool tag_isHighPt;
  float tag_relIso04;
  float tag_miniIso;
  float tag_miniIsoCharged;
  float tag_miniIsoPhotons;
  float tag_miniIsoNeutrals;
  bool tag_isMatchedGen;
  float tag_minDR;
  float tag_ptRel_minDR;
  float tag_iso03_sumPt;
  float tag_pfIso04_charged;
  float tag_pfIso04_neutral;
  float tag_pfIso04_photon;
  float tag_pfIso04_sumPU;
  float tag_tuneP_pt;
  float tag_tuneP_pterr;
  int tag_nsegments;
  bool tag_hasTrackMatch;
  float tag_TrackMatchDR;

  // Probe properties
  float probe_pt;
  float probe_eta;
  float probe_phi;
  int probe_charge;
  float probe_vtx_x;
  float probe_vtx_y;
  float probe_vtx_z;
  // bool probe_isLoose;
  // bool probe_isMedium;
  // bool probe_isTight;
  // bool probe_isSoft;
  // bool probe_isHighPt;
  // bool probe_isArbitratedTracker;
  // bool probe_isMuMatched;
  // bool probe_isPF;
  // bool probe_isSA;
  // bool probe_isTracker;
  // bool probe_isGlobal;
  // bool probe_isdSA;
  // bool probe_isdGlobal;
  bool probe_isCosmic;
  int probe_ncosmic;
  float probe_cosmic_minDR;
  // bool probe_isGood;
  bool probe_isHighPurity;
  float probe_relIso04;
  float probe_miniIso;
  float probe_miniIsoCharged;
  float probe_miniIsoPhotons;
  float probe_miniIsoNeutrals;
  bool probe_isMatchedGen;

  // float probe_validFraction;
  // float probe_trkChi2;
  // float probe_positionChi2;
  // float probe_trkKink;
  // float probe_segmentCompatibility;
  // float probe_trackerLayers;
  // float probe_pixelLayers;
  // float probe_muonStations;
  // float probe_muonHits;
  // float probe_DTHits;
  // float probe_CSCHits;
  // float probe_pterr;
  // float probe_dxy;
  // float probe_dz;
  // float probe_minDR;
  // float probe_ptRel_minDR;
  // float probe_iso03_sumPt;
  // float probe_pfIso04_charged;
  // float probe_pfIso04_neutral;
  // float probe_pfIso04_photon;
  // float probe_pfIso04_sumPU;
  // int probe_pixelHits;
  // int probe_matchedStations;
  // int probe_expectedMatchedStations;
  // int probe_RPCLayers;
  // unsigned int probe_stationMask;
  // int probe_nShowers;
  // float probe_tuneP_pt;
  // float probe_tuneP_pterr;
  // int probe_tuneP_muonHits;
  // int probe_nsegments;
  bool probe_hasMuonMatch;
  float probe_MuonMatchDR;
  bool probe_hasElectronMatch;
  float probe_ElectronMatchDR;
  bool probe_hasPFCMatch;
  float probe_PFCMatchDR;
  int probe_pfcID;

  // Pair properties
  float pair_pt;
  float pair_mass;
  float pair_eta;
  float pair_phi;
  float pair_fit_mass;
  float pair_svprob;
  float pair_normalchi2;
  float pair_dz;
  float pair_dR;
  float pair_dphi;
  int pair_rank_vtx_prob;
  int pair_rank_dPhi_muons;
  int pair_rank_Mass_Mmumu;

  float refit_vtx_x;
  float refit_vtx_y;
  float refit_vtx_z;
  float tag_refit_transverse_IP;
  float probe_refit_transverse_IP;
  float tag_transverse_IP;
  float probe_transverse_IP;
  
  float pair_tuneP_pt;
  float pair_tuneP_mass;
  float pair_tuneP_eta;
  float pair_tuneP_phi;
  float pair_tuneP_fit_mass;
  float pair_tuneP_svprob;
  float pair_tuneP_normalchi2;
  float pair_tuneP_dz;
  float pair_tuneP_dR;
  float pair_tuneP_dphi;

  float refit_tuneP_vtx_x;
  float refit_tuneP_vtx_y;
  float refit_tuneP_vtx_z;
  float tag_refit_tuneP_transverse_IP;
  float probe_refit_tuneP_transverse_IP;
  float tag_tuneP_transverse_IP;
  float probe_tuneP_transverse_IP;
  
   int nMu;
   std::vector<bool>  mu_isProbe;
   std::vector<bool>  mu_isTag;   
   std::vector<float>  muPt;
   std::vector<float>  muEta;
   std::vector<float>  muPhi;
   std::vector<int>    muCharge;
   std::vector<int>    muType;
   std::vector<int>    muIsGood;
   
   std::vector<int>    muIsGlobal;
   std::vector<int>    muIsTracker;
   std::vector<int>    muIsPF;
   std::vector<int>    muIsSTA;

   std::vector<float>  muD0;
   std::vector<float>  muDz;
   std::vector<float>  muIP3D;
   std::vector<float>  muD0Err;
   std::vector<float>  muDzErr;
   std::vector<float>  muIP3DErr;
   std::vector<float>  muChi2NDF;
   std::vector<float>  muInnerD0;
   std::vector<float>  muInnerDz;

   std::vector<float>  muInnerD0Err;
   std::vector<float>  muInnerDzErr;
   std::vector<float>  muInnerPt;
   std::vector<float>  muInnerPtErr;
   std::vector<float>  muInnerEta;

   std::vector<int>    muTrkLayers;
   std::vector<int>    muPixelLayers;
   std::vector<int>    muPixelHits;
   std::vector<int>    muMuonHits;
   std::vector<int>    muTrkQuality;
   std::vector<int>    muStations;
   std::vector<float>  muIsoTrk;
   std::vector<float>  muPFChIso;
   std::vector<float>  muPFPhoIso;
   std::vector<float>  muPFNeuIso;
   std::vector<float>  muPFPUIso;
   //https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_selectors_Since_9_4_X
   std::vector<int>    muIDSoft;
   std::vector<int>    muIDLoose;
   std::vector<int>    muIDMedium;
   std::vector<int>    muIDMediumPrompt;
   std::vector<int>    muIDTight;
   std::vector<int>    muIDGlobalHighPt;
   std::vector<int>    muIDTrkHighPt;
   std::vector<int>    muIDInTime;


   int nTrk;
   std::vector<bool> trkisTag;
   std::vector<bool> trkisProbe;   
   std::vector<float> trkPt;            
   std::vector<float> trkP;            
   std::vector<float> trkEta;           
   std::vector<float> trkPhi;     
   std::vector<int>   trkcharge; 
   std::vector<float> trkvx;            
   std::vector<float> trkvy;           
   std::vector<float> trkvz;     
   std::vector<float> trknormchi2;     
   std::vector<float> trkchi2;  
   std::vector<float> trkd0;                         
   std::vector<float> trkdxy;                      
   std::vector<float> trkdz;  
   // std::vector<float> trkdxyBS;
   // std::vector<float> trkdzBS;
   std::vector<float> trkdxyError;  
   std::vector<float> trkdzError;   
   std::vector<int>   trkValidHits;           
   std::vector<int>   trkMissHits;   
   std::vector<bool>   trkPurity; 

   int nEle;
   std::vector<bool>   ele_isProbe;
   std::vector<bool>   ele_isTag;
   std::vector<int>    eleCharge;
   std::vector<int>    eleChargeConsistent;
   std::vector<int>    eleSCPixCharge;
   // std::vector<int>    eleCtfCharge;
   std::vector<float>  eleEn;
   std::vector<float>  eleD0;
   std::vector<float>  eleDz;
   std::vector<float>  eleIP3D;
   std::vector<float>  eleD0Err;
   std::vector<float>  eleDzErr;
   std::vector<float>  eleIP3DErr;
   std::vector<float>  eleTrkPt;
   std::vector<float>  eleTrkEta;
   std::vector<float>  eleTrkPhi;
   std::vector<int>    eleTrkCharge;
   std::vector<float>  eleTrkPtErr;
   std::vector<float>  eleTrkChi2;
   std::vector<float>  eleTrkNdof;
   std::vector<float>  eleTrkNormalizedChi2;
   std::vector<int>    eleTrkValidHits;
   std::vector<int>    eleTrkLayers;
   std::vector<float>  elePt;
   std::vector<float>  eleEta;
   std::vector<float>  elePhi;
   std::vector<float>  eleSCEn;
   std::vector<float>  eleESEn;
   std::vector<float>  eleSCEta;
   std::vector<float>  eleSCPhi;
   std::vector<float>  eleSCRawEn;
   std::vector<float>  eleSCEtaWidth;
   std::vector<float>  eleSCPhiWidth;
   std::vector<float>  eleHoverE;
   std::vector<float>  eleHoverEBc;
   std::vector<float>  eleEoverP;
   std::vector<float>  eleEoverPInv;
   std::vector<float>  eleEcalE;
   std::vector<float>  elePAtVtx;
   std::vector<float>  elePAtSC;
   std::vector<float>  elePAtCluster;
   std::vector<float>  elePAtSeed;
   std::vector<float>  eleBrem;
   std::vector<float>  eledEtaAtVtx;
   std::vector<float>  eledPhiAtVtx;
   std::vector<float>  eledEtaSeedAtVtx;
   std::vector<float>  eleSigmaIEtaIEta;
   std::vector<float>  eleSigmaIEtaIEta_2012;
   std::vector<float>  eleSigmaIPhiIPhi;
   std::vector<int>    eleMissHits;
  
   int nPFCands;
   std::vector<bool>  pfcand_isProbe;
  // static const int maxsize = 500;
  std::vector<int> pfcand_pdgId;
  std::vector<int> pfcand_charge;
  std::vector<float> pfcand_pt;
  std::vector<float> pfcand_eta;
  std::vector<float> pfcand_phi;
  std::vector<float> pfcand_vtx_x; 
  std::vector<float> pfcand_vtx_y; 
  std::vector<float> pfcand_vtx_z;   
  
   int nTower;
   float maxHFp;
   float maxHFm;
   std::vector<float> CaloTower_hadE;
   std::vector<float> CaloTower_emE;
   std::vector<float> CaloTower_e;
   std::vector<float> CaloTower_et;
   std::vector<float> CaloTower_eta;
   std::vector<float> CaloTower_phi;
   std::vector<bool> CaloTower_HF_AboveThreshold;
   std::vector<bool> CaloTower_Had_AboveThreshold;
   std::vector<bool> CaloTower_EM_AboveThreshold;
   std::vector<bool> CaloTower_HF_inNoiseRegion;
   std::vector<bool> CaloTower_Had_inNoiseRegion;
   std::vector<bool> CaloTower_EM_inNoiseRegion;
  
  int ZDC_n;
  float  ZDC_e[18];
  int    ZDC_zside[18];
  int    ZDC_section [18];
  int    ZDC_channel[18];
  int    ZDC_saturation[18];
  float  ZDC_PM_Total_Energy;  
  float  ZDC_P_Total_Energy;
  float  ZDC_P_ECal_Energy;
  float  ZDC_P_HCal_Energy;          
  float  ZDC_M_Total_Energy;
  float  ZDC_M_ECal_Energy;
  float  ZDC_M_HCal_Energy; 
  
  //Photon info
  int nPho;
  
   std::vector<float>  phoE;
   std::vector<float>  phoEt;
   std::vector<float>  phoEta;
   std::vector<float>  phoPhi;

   std::vector<float>  phoEcorrStdEcal;
   std::vector<float>  phoEcorrPhoEcal;
   std::vector<float>  phoEcorrRegr1;
   std::vector<float>  phoEcorrRegr2;
   std::vector<float>  phoEcorrErrStdEcal;
   std::vector<float>  phoEcorrErrPhoEcal;
   std::vector<float>  phoEcorrErrRegr1;
   std::vector<float>  phoEcorrErrRegr2;

   std::vector<float>  phoSCE;
   std::vector<float>  phoSCEt;
   std::vector<float>  phoSCRawE;
   std::vector<float>  phoSCEta;
   std::vector<float>  phoSCPhi;
   std::vector<float>  phoSCEtaWidth;
   std::vector<float>  phoSCPhiWidth;
   std::vector<float>  phoSCBrem;
   std::vector<int>    phoSCnHits;
   std::vector<uint32_t> phoSCflags;
   std::vector<int>    phoSCinClean;
   std::vector<int>    phoSCinUnClean;
   std::vector<int>    phoSCnBC;
   std::vector<float>  phoESEn;

   std::vector<int>    phoIsPFPhoton;
   std::vector<int>    phoIsStandardPhoton;
   std::vector<int>    phoHasPixelSeed;
   std::vector<int>    phoHasConversionTracks;
// std::vector<int>    phoEleVeto;         // TODO: not available in reco::
   std::vector<float>  phoHadTowerOverEm;
   std::vector<float>  phoHoverE;
   std::vector<int>    phoHoverEValid;
   std::vector<float>  phoSigmaIEtaIEta;
   std::vector<float>  phoR9;
  
  //probe types
  bool probe_isMuon;
  bool probe_isElectron;
  bool probe_isPion;
  bool probe_isOther;
  int probe_typeSum;
  
  
  
  // 3 prong things
  
  
  int iprobe_3prong;
  int npairs_3prong;
  
  float  tau_3prong_trk1_pt;
  float  tau_3prong_trk1_eta;
  float  tau_3prong_trk1_phi;
  int tau_3prong_trk1_charge;
  float  tau_3prong_trk1_vtx_x;
  float  tau_3prong_trk1_vtx_y;
  float  tau_3prong_trk1_vtx_z;
        
  float  tau_3prong_trk2_pt;
  float  tau_3prong_trk2_eta;
  float  tau_3prong_trk2_phi;
  int  tau_3prong_trk2_charge;
  float  tau_3prong_trk2_vtx_x;
  float  tau_3prong_trk2_vtx_y;
  float  tau_3prong_trk2_vtx_z;
     
  float  tau_3prong_trk3_pt;
  float  tau_3prong_trk3_eta;
  float  tau_3prong_trk3_phi;
  int  tau_3prong_trk3_charge;
  float  tau_3prong_trk3_vtx_x;
  float  tau_3prong_trk3_vtx_y;
  float  tau_3prong_trk3_vtx_z;
        
  float  tau_3prong_total_pt;
  float  tau_3prong_total_eta;
  float  tau_3prong_total_phi;
  float  tau_3prong_total_M;
  float  tau_3prong_total_pt_sum;
  int    tau_3prong_total_charge;
  float  tau_3prong_total_Dz;
  float  tau_3prong_total_vtx_x;
  float  tau_3prong_total_vtx_y;
  float  tau_3prong_total_vtx_z;
  float  tau_3prong_total_vtx_prob;
  float  tau_3prong_total_vtx_chi2;
  float  tau_3prong_DR_1_2;
  float  tau_3prong_DR_1_3;
  float  tau_3prong_DR_2_3;
  float  tau_3prong_DR_tag_1;
  float  tau_3prong_DR_tag_2;
  float  tau_3prong_DR_tag_3;
  float  tau_3prong_total_area;
  
  
  std::vector<bool>  pfcand_isFirstProbe;
  std::vector<bool>  pfcand_isSecondProbe;
  std::vector<bool>  pfcand_isThirdProbe;
 
  std::vector<bool>  trk_isFirstProbe;
  std::vector<bool>  trk_isSecondProbe;
  std::vector<bool>  trk_isThirdProbe;
 
  int pair_rank_vtx_prob_3prong;
  int pair_rank_dPhi_muons_3prong;
  int pair_rank_Mass_Mmumu_3prong;


private:
  TTree *t1;
  TTree *t2;
  TTree *t3;
  TTree *t4;
  TTree *Eff_Tree;
  
  bool keepMuons= true ;
  bool keepElectrons= true ;
  bool keepTracks = true ;
  bool keepPFcands= true ;
  bool keepPhotons= true ;
  bool keepCaloTowers= true ;
  bool keepZDC= true ;
};
#endif
