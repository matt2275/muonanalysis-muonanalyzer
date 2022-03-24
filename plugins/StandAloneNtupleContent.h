//
// Original Author:
//                george karathanasis
//         Created:  Thu, 23 Mar 2019 17:40:23 GMT
//
// flat tree branches/var declaration

#ifndef STANDALONENTUPLECONTENT_H
#define STANDALONENTUPLECONTENT_H
#include <string>
#include <vector>
#include "TString.h"
#include "TTree.h"

class StandAloneNtupleContent {
public:
  StandAloneNtupleContent();
  virtual ~StandAloneNtupleContent();
  void SetTree(TTree *t1);
  void CreateBranches(const std::vector<std::string> &);
  void CreateExtraTrgBranches(const std::vector<std::string> &, bool);
  void ClearBranches();

  // Standard stuff
  // Standard stuff
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

  // Pileup
  float trueNumInteractions;
  int puNumInteractions;
  double Rho;

  // Number of muons
  int nmuons;
  int ntag;
  int iprobe;
  int npairs;

  // genmu1: mu-, genmu2: mu+
  float genmu1_pt;
  float genmu1_eta;
  float genmu1_phi;
  float genmu1_charge;
  float genmu2_pt;
  float genmu2_eta;
  float genmu2_phi;
  float genmu2_charge;
  float genMass;

  // genmuFSfromHP1: mu-, genmuFSfromHP2: mu+
  float genmuFSfromHP1_pt;
  float genmuFSfromHP1_eta;
  float genmuFSfromHP1_phi;
  float genmuFSfromHP1_charge;
  float genmuFSfromHP2_pt;
  float genmuFSfromHP2_eta;
  float genmuFSfromHP2_phi;
  float genmuFSfromHP2_charge;
  float genMassFSfromHP;

  // Tag properties
  float tag_pt;
  float tag_eta;
  float tag_phi;
  int tag_charge;
  float tag_pterr;
  float tag_dxy;
  float tag_dz;
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

  // Probe properties
  float probe_pt;
  float probe_eta;
  float probe_phi;
  int probe_charge;
  float probe_mupt;
  float probe_mueta;
  float probe_muphi;
  float probe_mucharge;
  float probe_mu_SAmu_DeltaR;
  bool probe_isLoose;
  bool probe_isMedium;
  bool probe_isTight;
  bool probe_isSoft;
  bool probe_isHighPt;
  bool probe_isArbitratedTracker;
  bool probe_isMuMatched;
  bool probe_isPF;
  bool probe_isSA;
  bool probe_isTracker;
  bool probe_isGlobal;
  bool probe_isdSA;
  bool probe_isdGlobal;
  bool probe_isCosmic;
  int probe_ncosmic;
  float probe_cosmic_minDR;
  bool probe_isGood;
  //  bool probe_isHighPurity;
  float probe_relIso04;
  float probe_miniIso;
  float probe_miniIsoCharged;
  float probe_miniIsoPhotons;
  float probe_miniIsoNeutrals;
  bool probe_isMatchedGen;

  float probe_validFraction;
  float probe_trkChi2;
  float probe_positionChi2;
  float probe_trkKink;
  float probe_segmentCompatibility;
  float probe_trackerLayers;
  float probe_pixelLayers;
  float probe_muonStations;
  float probe_muonHits;
  float probe_DTHits;
  float probe_CSCHits;
  float probe_pterr;
  float probe_dxy;
  float probe_dz;
  float probe_minDR;
  float probe_ptRel_minDR;
  float probe_iso03_sumPt;
  float probe_pfIso04_charged;
  float probe_pfIso04_neutral;
  float probe_pfIso04_photon;
  float probe_pfIso04_sumPU;
  int probe_pixelHits;
  int probe_matchedStations;
  int probe_expectedMatchedStations;
  int probe_RPCLayers;
  unsigned int probe_stationMask;
  int probe_nShowers;
  int probe_nsegments;

  // probe track information
  bool probe_isTrkMatch;
  float probe_trkPt;
  float probe_trkEta;
  float probe_trkPhi;
  int probe_trkCharge;
  float probe_trkDxy;
  float probe_trkDz;
  int probe_trkHits;
  int probe_trkStripHits;
  int probe_trkPixelHits;
  float probe_trk_SAmu_DeltaR;

  bool probeSA_isTrkMatch;
  float probeSA_trkPt;
  float probeSA_trkEta;
  float probeSA_trkPhi;
  int probeSA_trkCharge;
  float probeSA_trkDxy;
  float probeSA_trkDz;
  int probeSA_trkHits;
  int probeSA_trkStripHits;
  int probeSA_trkPixelHits;
  float probeSA_trk_SAmu_DeltaR;

  // Pair properties
  float pair_pt;
  float pair_mass;
  float pair_eta;
  float pair_phi;
  float pair_dz;
  float pair_dR;
  int pair_rank_vtx_prob;
  int pair_rank_dz_PV_SV;
  int pair_rank_dPhi_muons;
  int pair_rank_dM_Z_Mmumu;

  // sim matching information
  int tag_simType;
  int tag_simExtType;
  int tag_simFlavour;
  int tag_simHeaviestMotherFlavour;
  int tag_simPdgId;
  int tag_simMotherPdgId;
  int tag_simBX;
  float tag_simProdRho;
  float tag_simProdZ;
  float tag_simPt;
  float tag_simEta;
  float tag_simPhi;

  int probe_simType;
  int probe_simExtType;
  int probe_simFlavour;
  int probe_simHeaviestMotherFlavour;
  int probe_simPdgId;
  int probe_simMotherPdgId;
  int probe_simBX;
  float probe_simProdRho;
  float probe_simProdZ;
  float probe_simPt;
  float probe_simEta;
  float probe_simPhi;
  
//commented out to save space, added useful details for HIUPC collisions

  // int nPFCands;
  // // static const int maxsize = 500;
  // std::vector<int> pfcand_pdgId;
  // std::vector<int> pfcand_charge;
  // std::vector<float> pfcand_pt;
  // std::vector<float> pfcand_eta;
  // std::vector<float> pfcand_phi; 
  
   // int nTower;
   // std::vector<float> CaloTower_hadE;
   // std::vector<float> CaloTower_emE;
   // std::vector<float> CaloTower_e;
   // std::vector<float> CaloTower_et;
   // std::vector<float> CaloTower_eta;
   // std::vector<float> CaloTower_phi;
  
  // int ZDC_n;
  // float  ZDC_e[18];
  // int    ZDC_zside[18];
  // int    ZDC_section [18];
  // int    ZDC_channel[18];
  // int    ZDC_saturation[18];
  // float  ZDC_PM_Total_Energy;  
  // float  ZDC_P_Total_Energy;
  // float  ZDC_P_ECal_Energy;
  // float  ZDC_P_HCal_Energy;          
  // float  ZDC_M_Total_Energy;
  // float  ZDC_M_ECal_Energy;
  // float  ZDC_M_HCal_Energy; 
  
  float probe_numofassoctrks;
  
  
  static const int NTRIGGERMAX = 100;
  bool trigger[NTRIGGERMAX];

  // Trigger matches
  bool tag_trg[NTRIGGERMAX];
  float tag_trg_pt[NTRIGGERMAX];
  float tag_trg_eta[NTRIGGERMAX];
  float tag_trg_phi[NTRIGGERMAX];
  float tag_trg_dr[NTRIGGERMAX];

private:
  TTree *t1;
};
#endif
