// Package:    MuonAnalyzer for Run 3
//             version 2.0
//
/**\class

 Description: Ntuplizer class for full AOD files
*/
//
// Original Author:
//                george karathanasis
//         Created:  Thu, 20 feb 2020 17:40:23 GMT
//
// Modified:
//                Andre Frankenthal (Sept. 2020)
//
// Modified:
//                Minseok Oh (Feb. 2021)
//
//

// system include files
#include <iostream>
#include <memory>
#include <random>
#include <queue>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include <vector>
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonSimInfo.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"


//Adding particle Flow candidates 
// #include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
// #include "DataFormats/ParticleFlowReco/interface/PFSimParticle.h"
// #include "DataFormats/ParticleFlowReco/interface/PFSimParticleFwd.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

// #include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
// #include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"


// #include "DataFormats/Math/interface/AlgebraicROOTObjects.h"
// #include "KlFitter.h"
// #include "MuonBranches.h"
// #include "StandAloneMuonBranches.h"
// #include "MuonGenAnalyzer.h"
// #include "NtupleContent.h"
// #include "StandAloneNtupleContent.h"

#include "Analysis_KlFitter.h"
#include "Analysis_MuonBranches.h"
#include "Analysis_MuonGenAnalyzer.h"
#include "Analysis_CaloAnalyzer.h"
#include "Analysis_NtupleContent.h"
#include "Analysis_MuonMiniIsolation.h"

#include "helper.h"
// #include "MuonMiniIsolation.h"
#include "JetsBranches.h"

// #include "Riostream.h"
// #include "TMath.h"
// #include "Math/Vector3D.h"
// #include "TMatrixD.h"
// #include "TVectorD.h"
// #include "TGraphErrors.h"
// #include "TDecompChol.h"
// #include "TDecompSVD.h"
// #include "TF1.h"

#include <cmath>
#include <bits/stdc++.h>
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include <tuple>

using namespace std;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.
class HIUPC_Analysis_GenEfficiency : public edm::one::EDAnalyzer<> {
public:

struct PFC_3prong_data
{
   unsigned int trk_idx;
   unsigned int pfc_idx;
   double trk_pt;
   int trk_charge;
};


struct tau_cand_3prong
{
   unsigned int trk1_idx;
   unsigned int trk2_idx;
   unsigned int trk3_idx;
   unsigned int pfc1_idx;
   unsigned int pfc2_idx;
   unsigned int pfc3_idx;
   double pt_sum;
   int charge;
   TLorentzVector Tau;
};

struct tag_cand
{
   int tag_idx;
   double pt;
   string tag_type;
};


static bool compare_PFC_3prong_data_pt(PFC_3prong_data i1, PFC_3prong_data i2)
{
    return (i1.trk_pt > i2.trk_pt);
}

static bool compare_Tau_3prong_data_pt(tau_cand_3prong i1, tau_cand_3prong i2)
{
    return (i1.Tau.Pt() > i2.Tau.Pt());
}

static bool compare_tag_pt(tag_cand i1, tag_cand i2)
{
    return (i1.pt > i2.pt);
}
  

  typedef std::vector<std::pair<reco::Muon, reco::TransientTrack>> RecoTrkAndTransientTrkCollection;
  typedef std::vector<std::pair<reco::GsfElectron, reco::TransientTrack>> RecoElectronTrkAndTransientTrkCollection;
  explicit HIUPC_Analysis_GenEfficiency(const edm::ParameterSet&);
  ~HIUPC_Analysis_GenEfficiency() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  bool HLTaccept(const edm::Event&, Analysis_NtupleContent&, std::vector<std::string>&);
  void fillHLTmuon(const edm::Event&,
                   std::vector<TString>&,
                   std::vector<float>&,
                   std::vector<float>&,
                   std::vector<float>&,
                   std::vector<std::string>&,
                   const int&);
  void embedTriggerMatching(const reco::Track&,
                            std::vector<TString>&,
                            std::vector<float>&,
                            std::vector<float>&,
                            std::vector<float>&,
                            std::vector<std::string>&,
                            bool,
                            const int&);
                            
  void embedTriggerMatching_GenEfficiency(const reco::Track&,
                            std::vector<TString>&,
                            std::vector<float>&,
                            std::vector<float>&,
                            std::vector<float>&,
                            std::vector<std::string>&,
                            bool);
  // void StandAlone_embedTriggerMatching(const reco::Muon&,
                                       // std::vector<TString>&,
                                       // std::vector<float>&,
                                       // std::vector<float>&,
                                       // std::vector<float>&,
                                       // std::vector<std::string>&,
                                       // bool,
                                       // const int&);
                                       
                                       
  float Find3ProngArea( float,float,float,float,float,float);
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileupSummaryToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken_;
  edm::EDGetToken muonsToken_;
  edm::EDGetToken electronsToken_;
  edm::EDGetTokenT<edm::View<reco::Photon>> photonsToken_;
  edm::EDGetTokenT<edm::View<reco::Muon>> muonsViewToken_;
  edm::EDGetToken tracksToken_;
  edm::EDGetTokenT<edm::SortedCollection<CaloTower>> calotowerToken_;  
//  edm::EDGetToken centToken_;
// only used for J/Psi centrality 
 // edm::EDGetTokenT<reco::Centrality> CentralityTag_;
 // edm::EDGetTokenT<int> CentralityBinTag_;
 
  edm::EDGetTokenT<edm::SortedCollection<ZDCRecHit>> RecHitsToken_;
  edm::EDGetToken SAmuonsToken_;
  edm::EDGetToken dSAToken_;
  edm::EDGetToken dglToken_;
  edm::EDGetToken cosmicToken_;
  edm::EDGetTokenT<edm::TriggerResults> trgresultsToken_;
  edm::EDGetTokenT<trigger::TriggerEvent> trigobjectsToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneMatch> l1MatchesToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> l1MatchesQualityToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> l1MatchesDeltaRToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneMatch> l1MatchesByQToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> l1MatchesByQQualityToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> l1MatchesByQDeltaRToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle>> genToken_;
  edm::EDGetToken PFCands_;
  edm::EDGetTokenT<double> rhoJetsNC_;
  edm::EDGetToken jetsToken_;
  edm::EDGetToken jetCorrectorToken_;
  edm::EDGetToken genJetsToken_;
  edm::EDGetToken deepCSVProbbToken_;
  edm::EDGetToken deepCSVProbbbToken_;
  //  edm::EDGetToken deepFlavProbbToken_;
  //  edm::EDGetToken deepFlavProbbbToken_;
  edm::EDGetTokenT<edm::ValueMap<reco::MuonSimInfo>> simInfoToken_;

  std::vector<std::string> HLTPaths_;      // trigger fired
  std::vector<std::string> tagFilters_;    // tag-trigger matching
  std::vector<std::string> probeFilters_;  // probe-trigger matching
  std::vector<std::string> probeSelectorNames_;
  std::vector<unsigned> probeSelectorBits_;

  std::mt19937 m_random_generator = std::mt19937(37428479);
  const bool isMC_, includeJets_;
  const std::string era_;
  const bool keepMuons_;
  const bool keepElectrons_;
  const bool keepTracks_;
  const bool keepPFcands_;
  const bool keepPhotons_;
  const bool keepCaloTowers_;
  const bool keepZDC_;
  const int debug_;
  const string MCType_;

  edm::Service<TFileService> fs;
  // TTree* t1;
  // TTree* t2;
  TTree* t3;
  // TTree* t4;
  TTree* Eff_Tree;
  Analysis_NtupleContent nt;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HIUPC_Analysis_GenEfficiency::HIUPC_Analysis_GenEfficiency(const edm::ParameterSet& iConfig)
    :  // inputs

      genEventInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo"))),
      rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("Rho"))),
      pileupSummaryToken_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupInfo"))),
      beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
      vtxToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
      muonsToken_(consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
      electronsToken_(consumes<std::vector<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
      photonsToken_(consumes<edm::View<reco::Photon>>(iConfig.getParameter<edm::InputTag>("photons"))),
      muonsViewToken_(consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
      tracksToken_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
      calotowerToken_(consumes<edm::SortedCollection<CaloTower>>(iConfig.getParameter<edm::InputTag>("CaloTowers"))), 
//      centToken_(consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("Centrality"))),
// only used for J/Psi centrality 
     // CentralityTag_(consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("CentralitySrc"))),
       // CentralityBinTag_(consumes<int>(iConfig.getParameter<edm::InputTag>("CentralityBinSrc"))), 
       
      RecHitsToken_(consumes<edm::SortedCollection<ZDCRecHit>>(iConfig.getParameter<edm::InputTag>("RecHits"))),         
      SAmuonsToken_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("SAmuons"))),
      dSAToken_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("dSAmuons"))),
      dglToken_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("dGlmuons"))),
      cosmicToken_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("staCosmic"))),
      trgresultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
      trigobjectsToken_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
      l1MatchesToken_(consumes<pat::TriggerObjectStandAloneMatch>(iConfig.getParameter<edm::InputTag>("l1Matches"))),
      l1MatchesQualityToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("l1MatchesQuality"))),
      l1MatchesDeltaRToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("l1MatchesDeltaR"))),
      l1MatchesByQToken_(
          consumes<pat::TriggerObjectStandAloneMatch>(iConfig.getParameter<edm::InputTag>("l1MatchesByQ"))),
      l1MatchesByQQualityToken_(
          consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("l1MatchesByQQuality"))),
      l1MatchesByQDeltaRToken_(
          consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("l1MatchesByQDeltaR"))),
      genToken_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("gen"))),
      PFCands_(consumes<std::vector<reco::PFCandidate>>(iConfig.getParameter<edm::InputTag>("PFCands"))),
      rhoJetsNC_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoJetsNC"))),
      jetsToken_(consumes<std::vector<reco::PFJet>>(iConfig.getParameter<edm::InputTag>("jets"))),
      jetCorrectorToken_(consumes<reco::JetCorrector>(iConfig.getParameter<edm::InputTag>("jetCorrector"))),
      genJetsToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genJets"))),
      deepCSVProbbToken_(consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("deepCSVProbb"))),
      deepCSVProbbbToken_(consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("deepCSVProbbb"))),
      simInfoToken_(consumes<edm::ValueMap<reco::MuonSimInfo>>(iConfig.getParameter<edm::InputTag>("muonSimInfo"))),
      HLTPaths_(iConfig.getParameter<std::vector<std::string>>("triggerPaths")),
      tagFilters_(iConfig.getParameter<std::vector<std::string>>("tagFilters")),
      probeFilters_(iConfig.getParameter<std::vector<std::string>>("probeFilters")),
      probeSelectorNames_(iConfig.getParameter<std::vector<std::string>>("probeSelectorNames")),
      probeSelectorBits_(iConfig.getParameter<std::vector<unsigned>>("probeSelectorBits")),
      isMC_(iConfig.getParameter<bool>("isMC")),
      includeJets_(iConfig.getParameter<bool>("includeJets")),
      era_(iConfig.getParameter<std::string>("era")),
      keepMuons_(iConfig.getParameter<bool>("keepMuons")),
      keepElectrons_(iConfig.getParameter<bool>("keepElectrons")),
      keepTracks_(iConfig.getParameter<bool>("keepTracks")),
      keepPFcands_(iConfig.getParameter<bool>("keepPFcands")),
      keepPhotons_(iConfig.getParameter<bool>("keepPhotons")),
      keepCaloTowers_(iConfig.getParameter<bool>("keepCaloTowers")),
      keepZDC_(iConfig.getParameter<bool>("keepZDC")),
      debug_(iConfig.getParameter<int>("debug")),
      MCType_(iConfig.getParameter<std::string>("MCType")) {
  //  edm::ParameterSet
  //  runParameters=iConfig.getParameter<edm::ParameterSet>("RunParameters");

  if (probeSelectorNames_.size() != probeSelectorBits_.size()) {
    throw cms::Exception("ParameterError")
        << "length of probeSelectorNames and probeSelectorBits should be identical\n";
  }
}

HIUPC_Analysis_GenEfficiency::~HIUPC_Analysis_GenEfficiency() {
  // cout << "total " << trg_counter << " fires " << fire_counter << " l3"
  // << l3_counter << endl; do anything here that needs to be done at
  // desctruction time
}


  float HIUPC_Analysis_GenEfficiency::Find3ProngArea( float prong1_eta, float prong1_phi,
                             float prong2_eta, float prong2_phi, float prong3_eta, float prong3_phi){
                                
   float prong1_theta = 2*atan(exp(-prong1_eta));
   float prong2_theta = 2*atan(exp(-prong2_eta));
   float prong3_theta = 2*atan(exp(-prong3_eta));
   float prong1_x = cos(prong1_phi)*sin(prong1_theta);
   float prong1_y = sin(prong1_phi)*sin(prong1_theta);
   float prong1_z = cos(prong1_theta);
   float prong2_x = cos(prong2_phi)*sin(prong2_theta);
   float prong2_y = sin(prong2_phi)*sin(prong2_theta);
   float prong2_z = cos(prong2_theta);
   float prong3_x = cos(prong3_phi)*sin(prong3_theta);
   float prong3_y = sin(prong3_phi)*sin(prong3_theta);
   float prong3_z = cos(prong3_theta);
   float diff_1_2_x = prong1_x - prong2_x;
   float diff_1_2_y = prong1_y - prong2_y;
   float diff_1_2_z = prong1_z - prong2_z;
   float diff_3_2_x = prong3_x - prong2_x;
   float diff_3_2_y = prong3_y - prong2_y;
   float diff_3_2_z = prong3_z - prong2_z;
   float cross_x = diff_1_2_y*diff_3_2_z - diff_1_2_z*diff_3_2_y ;
   float cross_y = diff_1_2_z*diff_3_2_x - diff_1_2_x*diff_3_2_z ;
   float cross_z = diff_1_2_x*diff_3_2_y - diff_1_2_y*diff_3_2_x ;
   float norm_cross = sqrt( cross_x*cross_x + cross_y*cross_y + cross_z*cross_z);
   return (.5* norm_cross);
   
  }
  
  
bool HIUPC_Analysis_GenEfficiency::HLTaccept(const edm::Event& iEvent,
                                              Analysis_NtupleContent& nt,
                                              std::vector<std::string>& HLTPaths) {
  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trgresultsToken_, trigResults);
  edm::TriggerNames trigName;
  trigName = iEvent.triggerNames(*trigResults);
  bool EvtFire = false;
  unsigned int ipath = 0;
  for (auto path : HLTPaths) {
    bool TrgFire = false;
    for (unsigned int itrg = 0; itrg < trigResults->size(); ++itrg) {
      TString TrigPath = trigName.triggerName(itrg);
      if (!trigResults->accept(itrg))
        continue;
      if (!TrigPath.Contains(path))
        continue;
      EvtFire = true;
      TrgFire = true;
    }
    nt.trigger[ipath] = TrgFire;
    ipath++;
  }
  return EvtFire;
}

void HIUPC_Analysis_GenEfficiency::fillHLTmuon(const edm::Event& iEvent,
                                                std::vector<TString>& trg_filter,
                                                std::vector<float>& trg_pt,
                                                std::vector<float>& trg_eta,
                                                std::vector<float>& trg_phi,
                                                std::vector<std::string>& HLTFilters,
                                                const int& debug_) {
  edm::Handle<trigger::TriggerEvent> triggerObjects;
  iEvent.getByToken(trigobjectsToken_, triggerObjects);
  trigger::TriggerObjectCollection allTriggerObjects = triggerObjects->getObjects();
  for (auto ifilter : HLTFilters) {
    size_t filterIndex = (*triggerObjects).filterIndex(edm::InputTag(ifilter, "", "HLT"));
    if (filterIndex < (*triggerObjects).sizeFilters()) {
      const trigger::Keys& keys = (*triggerObjects).filterKeys(filterIndex);
      for (size_t j = 0; j < keys.size(); j++) {
        // L1 objects for HIUPC seem to have id =0 
        TString filter_name = TString(ifilter);
        trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
        if (!(fabs(foundObject.id()) == 11 || fabs(foundObject.id()) == 13 || ((fabs(foundObject.id()) == 0) && (filter_name.Contains("HLT_HIUPC") || filter_name.Contains("L1_Single") || filter_name.Contains("hltL1s")))))
          continue;
        trg_filter.push_back(TString(ifilter));
        trg_pt.push_back(foundObject.pt());
        trg_eta.push_back(foundObject.eta());
        trg_phi.push_back(foundObject.phi());
        if (debug_ > 0)
          std::cout << "Trg muon " << foundObject.pt() << std::endl;
      }
    }
  }
}


void HIUPC_Analysis_GenEfficiency::embedTriggerMatching(const reco::Track& mu,
                                                         std::vector<TString>& trg_filter,
                                                         std::vector<float>& trg_pt,
                                                         std::vector<float>& trg_eta,
                                                         std::vector<float>& trg_phi,
                                                         std::vector<std::string>& HLTFilters,
                                                         bool isTag,
                                                         const int& debug_ = 0) {
  
  // // Commented out because not using for GenEfficiency Study
  
  // for (const auto& thefilter : HLTFilters) {
    // TString thefilter_tstr = TString(thefilter);
    // // temporary method to tag L2 filters for dSA paths...
    // bool isL2DSA =
        // thefilter_tstr.BeginsWith("hltL2") && (thefilter_tstr.Contains("NoVtx") || thefilter_tstr.Contains("NoVertex"));

    // bool matched = false;
    // float matched_pt = -99;
    // float matched_eta = -99;
    // float matched_phi = -99;
    // float matched_dr = 99;
    // for (unsigned itrg = 0; itrg < trg_filter.size(); ++itrg) {
      // TString filter_tstr = TString(trg_filter.at(itrg));
      // if (!filter_tstr.Contains(thefilter_tstr))
        // continue;
      // float dR_tmp = deltaR(mu.eta(), mu.phi(), trg_eta.at(itrg), trg_phi.at(itrg));
      // if ((dR_tmp < matched_dr && dR_tmp < trgDRwindow_) || 
        // (isL2DSA && dR_tmp < matched_dr && dR_tmp < maxdr_trk_dsa_)) {
        // matched = true;
        // matched_pt = trg_pt.at(itrg);
        // matched_eta = trg_eta.at(itrg);
        // matched_phi = trg_phi.at(itrg);
        // matched_dr = dR_tmp;

        // if (debug_ > 0) {
          // std::cout << "embedTriggerMatching: isTag=" << isTag << "  filter=" << thefilter_tstr << "  dR=" << dR_tmp
                    // << "  matched=" << matched << std::endl;
        // }
      // }
    // }
    // if (isTag) {
      // nt.tag_trg[&thefilter - &HLTFilters[0]] = matched;
      // nt.tag_trg_pt[&thefilter - &HLTFilters[0]] = matched_pt;
      // nt.tag_trg_eta[&thefilter - &HLTFilters[0]] = matched_eta;
      // nt.tag_trg_phi[&thefilter - &HLTFilters[0]] = matched_phi;
      // nt.tag_trg_dr[&thefilter - &HLTFilters[0]] = matched_dr;
    // } else {
      // nt.probe_trg[&thefilter - &HLTFilters[0]] = matched;
      // nt.probe_trg_pt[&thefilter - &HLTFilters[0]] = matched_pt;
      // nt.probe_trg_eta[&thefilter - &HLTFilters[0]] = matched_eta;
      // nt.probe_trg_phi[&thefilter - &HLTFilters[0]] = matched_phi;
      // nt.probe_trg_dr[&thefilter - &HLTFilters[0]] = matched_dr;
    // }
  // }
 return;
}


void HIUPC_Analysis_GenEfficiency::embedTriggerMatching_GenEfficiency(const reco::Track& mu,
                                                         std::vector<TString>& trg_filter,
                                                         std::vector<float>& trg_pt,
                                                         std::vector<float>& trg_eta,
                                                         std::vector<float>& trg_phi,
                                                         std::vector<std::string>& HLTFilters,
                                                         bool isPos) {
  
  // Commented out because not using for GenEfficiency Study
  
  for (const auto& thefilter : HLTFilters) {
    TString thefilter_tstr = TString(thefilter);
    // temporary method to tag L2 filters for dSA paths...
    bool matched = false;
    float matched_pt = -99;
    float matched_eta = -99;
    float matched_phi = -99;
    float matched_dr = 99;
    for (unsigned itrg = 0; itrg < trg_filter.size(); ++itrg) {
      TString filter_tstr = TString(trg_filter.at(itrg));
      if (!filter_tstr.Contains(thefilter_tstr))
        continue;
      float dR_tmp = deltaR(mu.eta(), mu.phi(), trg_eta.at(itrg), trg_phi.at(itrg));
      if (dR_tmp < matched_dr && dR_tmp < 1) {
        matched = true;
        matched_pt = trg_pt.at(itrg);
        matched_eta = trg_eta.at(itrg);
        matched_phi = trg_phi.at(itrg);
        matched_dr = dR_tmp;

    }
    }
    if (isPos) {
      nt.tag_trg[&thefilter - &HLTFilters[0]] = matched;
      nt.tag_trg_pt[&thefilter - &HLTFilters[0]] = matched_pt;
      nt.tag_trg_eta[&thefilter - &HLTFilters[0]] = matched_eta;
      nt.tag_trg_phi[&thefilter - &HLTFilters[0]] = matched_phi;
      nt.tag_trg_dr[&thefilter - &HLTFilters[0]] = matched_dr;
    } else {
      nt.probe_trg[&thefilter - &HLTFilters[0]] = matched;
      nt.probe_trg_pt[&thefilter - &HLTFilters[0]] = matched_pt;
      nt.probe_trg_eta[&thefilter - &HLTFilters[0]] = matched_eta;
      nt.probe_trg_phi[&thefilter - &HLTFilters[0]] = matched_phi;
      nt.probe_trg_dr[&thefilter - &HLTFilters[0]] = matched_dr;
    }
  }
 return;
}
  
  void HIUPC_Analysis_GenEfficiency::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;
  t3->Fill();

  
  //clear branches
  nt.ClearBranches();
  
  // Get data
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(beamSpotToken_, theBeamSpot);
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);

  // Skip evts if there are no vertices
  // bool hasVertex = true;
  if (vertices->empty()){
   cout << " no vertex" << endl;
   // hasVertex = false;

   return;   
  }
 

  edm::Handle<std::vector<reco::Track>> SAmuons;
  iEvent.getByToken(SAmuonsToken_, SAmuons);

  edm::Handle<std::vector<reco::Muon>> muons;
  iEvent.getByToken(muonsToken_, muons);
  edm::Handle<std::vector<reco::GsfElectron>> electrons;
  iEvent.getByToken(electronsToken_, electrons);
  edm::Handle<edm::View<reco::Muon>> muonsView;
  iEvent.getByToken(muonsViewToken_, muonsView);
  edm::Handle<std::vector<reco::Track>> tracks;
  iEvent.getByToken(tracksToken_, tracks);
  edm::Handle<std::vector<reco::Track>> dSAmuons;
  iEvent.getByToken(dSAToken_, dSAmuons);
  edm::Handle<std::vector<reco::Track>> dGlmuons;
  iEvent.getByToken(dglToken_, dGlmuons);
  edm::Handle<std::vector<reco::Track>> staCosmic;
  iEvent.getByToken(cosmicToken_, staCosmic);
  // mini isolation
  edm::Handle<std::vector<reco::PFCandidate>> pfcands;
  iEvent.getByToken(PFCands_, pfcands);
  edm::Handle<double> rhoJetsNC;
  iEvent.getByToken(rhoJetsNC_, rhoJetsNC);
  // jets
  edm::Handle<std::vector<reco::PFJet>> jets;
  edm::Handle<reco::JetCorrector> jetCorrector;
  edm::Handle<reco::JetTagCollection> deepCSVProbb;
  edm::Handle<reco::JetTagCollection> deepCSVProbbb;
  JME::JetResolution resolution;
  JME::JetResolutionScaleFactor resolution_sf;
  if (includeJets_) {
    iEvent.getByToken(jetsToken_, jets);
    iEvent.getByToken(jetCorrectorToken_, jetCorrector);
    iEvent.getByToken(deepCSVProbbToken_, deepCSVProbb);
    iEvent.getByToken(deepCSVProbbbToken_, deepCSVProbbb);
    //  edm::Handle<reco::JetTagCollection> deepFlavProbb;
    //  iEvent.getByToken(deepFlavProbbToken_, deepFlavProbb);
    //  edm::Handle<reco::JetTagCollection> deepFlavProbbb;
    //  iEvent.getByToken(deepFlavProbbbToken_, deepFlavProbbb);
    resolution = JME::JetResolution::get(iSetup, "AK4PFchs_pt");
    resolution_sf = JME::JetResolutionScaleFactor::get(iSetup, "AK4PFchs");
  }
  edm::ESHandle<MagneticField> bField;
  iSetup.get<IdealMagneticFieldRecord>().get(bField);

  edm::Handle<trigger::TriggerEvent> triggerObjects;
  iEvent.getByToken(trigobjectsToken_, triggerObjects);
  edm::Handle<pat::TriggerObjectStandAloneMatch> l1Matches;
  iEvent.getByToken(l1MatchesToken_, l1Matches);
  edm::Handle<edm::ValueMap<int>> l1Qualities;
  iEvent.getByToken(l1MatchesQualityToken_, l1Qualities);
  edm::Handle<edm::ValueMap<float>> l1Drs;
  iEvent.getByToken(l1MatchesDeltaRToken_, l1Drs);
  edm::Handle<pat::TriggerObjectStandAloneMatch> l1MatchesByQ;
  iEvent.getByToken(l1MatchesByQToken_, l1MatchesByQ);
  edm::Handle<edm::ValueMap<int>> l1QualitiesByQ;
  iEvent.getByToken(l1MatchesByQQualityToken_, l1QualitiesByQ);
  edm::Handle<edm::ValueMap<float>> l1DrsByQ;
  iEvent.getByToken(l1MatchesByQDeltaRToken_, l1DrsByQ);

// // only used for J/Psi centrality 
// // added to deal with inclusive HI j/psi with a bunch of tracks

    // edm::Handle<reco::Centrality> centrality;
    // iEvent.getByToken(CentralityTag_, centrality);
    // // cout << "Cent N trk  " << centrality->Ntracks() << endl;
    
    
    // int hiBin;
    // edm::Handle<int> cbin_;
    // iEvent.getByToken(CentralityBinTag_,cbin_);
    // hiBin = *cbin_;
    // // cout << "hiBin  " << hiBin << endl;
    
    // if(hiBin < 175) return;

  // Information about run
  nt.run = iEvent.id().run();
  nt.ls = iEvent.luminosityBlock();
  nt.event = iEvent.id().event();
  nt.fromFullAOD = true;
  nt.BSpot_x = theBeamSpot->x0();
  nt.BSpot_y = theBeamSpot->y0();
  nt.BSpot_z = theBeamSpot->z0();
  nt.nvertices = vertices->size();
  
  
  Analysis_MuonGenAnalyzer genmu;
  
  // Added for gen to reco efficiency measurements
  
  
  // if(saveGenEfficiencyTree_) continue;
  
  if (!iEvent.isRealData()) {
    if(MCType_ == "TauTau") genmu.SetInputsandFillNtuple_TauTau(nt, iEvent, genToken_);
    else if(MCType_ == "MuMu") genmu.SetInputsandFillNtuple_MuMu(nt, iEvent, genToken_);
    else if(MCType_ == "EE") genmu.SetInputsandFillNtuple_EE(nt, iEvent, genToken_);
    else genmu.SetInputsandFillNtuple_Other(nt, iEvent, genToken_, .5);

  std::vector<reco::TransientTrack> test_trk_pair = {};
  Analysis_KlFitter genstudy_vtx(test_trk_pair);

  if( MCType_ == "TauTau")genstudy_vtx.genDiTau(nt);
  
  genstudy_vtx.genFinal(nt);

  if( MCType_ == "TauTau")genstudy_vtx.genDiTau_genFinal_vertex(nt); 
  
  }
  
  
    reco::TrackBase::Point vertex_point;
  bool goodVtx = false;
  reco::Vertex const* pv;
  for (const reco::Vertex& vtx : *vertices) {
     // HIUPC seem to have fake vertex here, maybe due to low particle multiplicity
    //if (vtx.isFake() || !vtx.isValid())
    if (!vtx.isValid())
      continue;
    nt.pv_x = vtx.x();
    nt.pv_y = vtx.y();
    nt.pv_z = vtx.z();
    goodVtx = true;
    nt.hasFakeVertex = vtx.isFake();
    pv = &vtx;
    break;
  }
  if (!goodVtx){
      const reco::Vertex& vtx = vertices->at(0);
      nt.pv_x = vtx.x();
      nt.pv_y = vtx.y();
      nt.pv_z = vtx.z();
     nt.hasValidVertex = false;
     nt.hasFakeVertex = vtx.isFake();
     pv = &vtx;
  cout << " no valid vertex " << endl;
  }
  
  if(goodVtx){
  nt.hasValidVertex = true;
  }

  vertex_point.SetCoordinates(nt.pv_x, nt.pv_y, nt.pv_z);


  fillHLTmuon(iEvent, nt.trg_filter, nt.trg_pt, nt.trg_eta, nt.trg_phi, tagFilters_, debug_);
  fillHLTmuon(iEvent, nt.prb_filter, nt.prb_pt, nt.prb_eta, nt.prb_phi, probeFilters_, debug_);
  
  
   bool passHLT = HLTaccept(iEvent, nt, HLTPaths_);

    bool hasPosTrk = false;
    unsigned int idx_trk_pos;
    float pt_trk_pos = 0;
    bool hasNegTrk = false;
    unsigned int idx_trk_neg;
    float pt_trk_neg = 0;
    for (const reco::Track& trk : *tracks) {
      if (trk.charge() > 0){
         if(trk.pt() < pt_trk_pos || trk.pt() > 1000) continue;
         hasPosTrk = true;
         pt_trk_pos = trk.pt();
         idx_trk_pos = &trk - &tracks->at(0);
         nt.tag_pt = trk.pt();
         nt.tag_eta = trk.eta();
         nt.tag_phi = trk.phi();
         nt.tag_charge = trk.charge();
         nt.tag_vtx_x = trk.vx();
         nt.tag_vtx_y = trk.vy();
         nt.tag_vtx_z = trk.vz();
      }
      if (trk.charge() < 0){
         if(trk.pt() < pt_trk_neg || trk.pt() > 1000) continue;
         hasNegTrk = true;
         pt_trk_neg = trk.pt();
         idx_trk_neg = &trk - &tracks->at(0);
         nt.probe_pt = trk.pt();
         nt.probe_eta = trk.eta();
         nt.probe_phi = trk.phi();
         nt.probe_charge = trk.charge();
         nt.probe_vtx_x = trk.vx();
         nt.probe_vtx_y = trk.vy();
         nt.probe_vtx_z = trk.vz();
      }
    }
   if(hasPosTrk) embedTriggerMatching_GenEfficiency(tracks->at(idx_trk_pos), nt.trg_filter, nt.trg_pt, nt.trg_eta, nt.trg_phi, tagFilters_, true);
   if(hasNegTrk) embedTriggerMatching_GenEfficiency(tracks->at(idx_trk_neg), nt.trg_filter, nt.trg_pt, nt.trg_eta, nt.trg_phi, probeFilters_, false);
   int fake_idx = -99 ;
   FillMuonBranches<reco::Muon>(*muons, nt,fake_idx, fake_idx, *pv);
   FillElectronBranches<reco::GsfElectron>(*electrons, nt, fake_idx, fake_idx, *pv);
        
   // FillPFCandBranches<reco::PFCandidate>(*pfcands, nt, -1);
   FillTrackBranches<reco::Track>(*tracks, nt, fake_idx,fake_idx); 
   
   Eff_Tree->Fill();
   
  
  
  // end of gen to reco efficiency measurements

}

// ------------ method called once each job just before starting event loop
// ------------
void HIUPC_Analysis_GenEfficiency::beginJob() {
  // t1 = fs->make<TTree>("Events", "Events");
  // t2 = fs->make<TTree>("GenVtxStudy","GenVtxStudy");
  t3 = fs->make<TTree>("test","test");
  // t4 = fs->make<TTree>("Events_3prong","Events_3prong");
  Eff_Tree = fs->make<TTree>("GenEfficiency","GenEfficiency");
  nt.SetTreeVariables(keepMuons_ ,keepElectrons_ ,keepTracks_ , keepPFcands_ , keepPhotons_ ,keepCaloTowers_ ,keepZDC_ );
  // nt.SetTree(t1);
  // if(saveCutTree_) nt.SetTree_GenVtxStudy(t2);
  nt.SetTree_Test(t3); 
  // nt.CreateBranches(HLTPaths_, probeSelectorNames_);
  // if(saveCutTree_) nt.CreateBranches_GenVtxStudy();
  // nt.SetTree_3ProngStudy(t4);
  // nt.CreateBranches_3ProngStudy(HLTPaths_, probeSelectorNames_);
  nt.SetTree_EfficiencyTree(Eff_Tree);
  nt.CreateBranches_EfficiencyTree(HLTPaths_, probeSelectorNames_);
nt.CreateExtraTrgBranches_EfficiencyTree(tagFilters_, true);
nt.CreateExtraTrgBranches_EfficiencyTree(probeFilters_, false);
  // if (!tagFilters_.empty()) {
    // nt.CreateExtraTrgBranches(tagFilters_, true);
    // nt.CreateExtraTrgBranches_3prong(tagFilters_, true);
    // // StandAlone_nt.CreateExtraTrgBranches(tagFilters_, true);
  // }
  // if (!probeFilters_.empty())
    // nt.CreateExtraTrgBranches(probeFilters_, false);
}

// ------------ method called once each job just after ending the event loop
// ------------
void HIUPC_Analysis_GenEfficiency::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the
// module  ------------
void HIUPC_Analysis_GenEfficiency::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no
  // validation
  // Please change this to state exactly what you do use, even if it is no
  // parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

///////////////////////

// define this as a plug-in
DEFINE_FWK_MODULE(HIUPC_Analysis_GenEfficiency);
