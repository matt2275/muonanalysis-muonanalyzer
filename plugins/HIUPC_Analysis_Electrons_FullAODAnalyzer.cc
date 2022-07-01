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
using namespace std;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.
class HIUPC_Analysis_Electrons_FullAODAnalyzer : public edm::one::EDAnalyzer<> {
public:
  typedef std::vector<std::pair<reco::Muon, reco::TransientTrack>> RecoTrkAndTransientTrkCollection;
  typedef std::vector<std::pair<reco::GsfElectron, reco::TransientTrack>> RecoElectronTrkAndTransientTrkCollection;
  explicit HIUPC_Analysis_Electrons_FullAODAnalyzer(const edm::ParameterSet&);
  ~HIUPC_Analysis_Electrons_FullAODAnalyzer() override;

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
  // void StandAlone_embedTriggerMatching(const reco::Muon&,
                                       // std::vector<TString>&,
                                       // std::vector<float>&,
                                       // std::vector<float>&,
                                       // std::vector<float>&,
                                       // std::vector<std::string>&,
                                       // bool,
                                       // const int&);
                                       
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileupSummaryToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken_;
  edm::EDGetToken muonsToken_;
  edm::EDGetToken electronsToken_;
  edm::EDGetTokenT<edm::View<reco::Muon>> muonsViewToken_;
  edm::EDGetToken tracksToken_;
  edm::EDGetTokenT<edm::SortedCollection<CaloTower>> calotowerToken_;  

//  edm::EDGetToken centToken_;
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
  const double trgDRwindow_;
  const unsigned int tagQual_;
  const StringCutObjectSelector<reco::Muon> tagSelection_;  // kinematic cuts for tag
  const StringCutObjectSelector<reco::GsfElectron> tagElectronSelection_;
  const bool HighPurity_;
  const StringCutObjectSelector<reco::Track> probeSelection_;    // kinematic cuts for probe
  const StringCutObjectSelector<reco::Track> probeSelectionSA_;  // kinematic cuts for probe
  const bool muonOnly_;
  const StringCutObjectSelector<reco::Muon> probeMuonSelection_;
  const double pairMassMin_;
  const double pairMassMax_;
  const double pairDz_;
  const bool RequireVtxCreation_;  // if true skip pairs that do not create
                                   // that do not have a vertex
  const double minSVtxProb_;       // min probability of a vertex to be kept. If <0 inactive
  const double maxdz_trk_mu_;
  const double maxpt_relative_dif_trk_mu_;
  const double maxdr_trk_mu_;
  const double maxdz_trk_ele_;
  const double maxpt_relative_dif_trk_ele_;
  const double maxdr_trk_ele_;
  const double maxdz_trk_pfc_;
  const double maxpt_relative_dif_trk_pfc_;
  const double maxdr_trk_pfc_;  
  const double minpt_trkSA_;
  const double maxdz_trk_SAmu_;
  const double maxpt_relative_dif_trk_SAmu_;
  const double maxdr_trk_SAmu_;
  const double maxdr_trk_dsa_;
  const double minTrkPt_;  
  const double maxTrkEta_;
  const double maxTrkNum_;
  const unsigned momPdgId_;
  const double genRecoDrMatch_;
  const bool saveCutTree_;
  const bool keepelectrontags_;
  const int debug_;
  const string MCType_;
  PropagateToMuon prop1_;

  edm::Service<TFileService> fs;
  TTree* t1;
  TTree* t2;
  TTree* t3;
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
HIUPC_Analysis_Electrons_FullAODAnalyzer::HIUPC_Analysis_Electrons_FullAODAnalyzer(const edm::ParameterSet& iConfig)
    :  // inputs

      genEventInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo"))),
      rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("Rho"))),
      pileupSummaryToken_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupInfo"))),
      beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
      vtxToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
      muonsToken_(consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
      electronsToken_(consumes<std::vector<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
      muonsViewToken_(consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
      tracksToken_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
      calotowerToken_(consumes<edm::SortedCollection<CaloTower>>(iConfig.getParameter<edm::InputTag>("CaloTowers"))), 
//      centToken_(consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("Centrality"))),  
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
      trgDRwindow_(iConfig.getParameter<double>("trgDRwindow")),
      tagQual_(iConfig.getParameter<unsigned>("tagQuality")),
      tagSelection_(iConfig.getParameter<std::string>("tagSelection")),
      tagElectronSelection_(iConfig.getParameter<std::string>("tagElectronSelection")),
      HighPurity_(iConfig.getParameter<bool>("probeHPurity")),
      probeSelection_(iConfig.getParameter<std::string>("probeSelection")),
      probeSelectionSA_(iConfig.getParameter<std::string>("probeSelectionSA")),
      muonOnly_(iConfig.getParameter<bool>("muonOnly")),
      probeMuonSelection_(iConfig.getParameter<std::string>("probeMuonSelection")),
      pairMassMin_(iConfig.getParameter<double>("pairMassMin")),
      pairMassMax_(iConfig.getParameter<double>("pairMassMax")),
      pairDz_(iConfig.getParameter<double>("pairDz")),
      RequireVtxCreation_(iConfig.getParameter<bool>("RequireVtxCreation")),
      minSVtxProb_(iConfig.getParameter<double>("minSVtxProb")),
      maxdz_trk_mu_(iConfig.getParameter<double>("maxDzProbeTrkMuon")),
      maxpt_relative_dif_trk_mu_(iConfig.getParameter<double>("maxRelPtProbeTrkMuon")),
      maxdr_trk_mu_(iConfig.getParameter<double>("maxDRProbeTrkMuon")),
      maxdz_trk_ele_(iConfig.getParameter<double>("maxDzProbeTrkElectron")),
      maxpt_relative_dif_trk_ele_(iConfig.getParameter<double>("maxRelPtProbeTrkElectron")),
      maxdr_trk_ele_(iConfig.getParameter<double>("maxDRProbeTrkElectron")),
      maxdz_trk_pfc_(iConfig.getParameter<double>("maxDzProbeTrkPFC")),
      maxpt_relative_dif_trk_pfc_(iConfig.getParameter<double>("maxRelPtProbeTrkPFC")),
      maxdr_trk_pfc_(iConfig.getParameter<double>("maxDRProbeTrkPFC")),      
      minpt_trkSA_(iConfig.getParameter<double>("minPtTrkSA")),
      maxdz_trk_SAmu_(iConfig.getParameter<double>("maxDzProbeTrkSAMuon")),
      maxpt_relative_dif_trk_SAmu_(iConfig.getParameter<double>("maxRelPtProbeTrkSAMuon")),
      maxdr_trk_SAmu_(iConfig.getParameter<double>("maxDRProbeTrkSAMuon")),
      maxdr_trk_dsa_(iConfig.getParameter<double>("maxDRProbeTrkDSA")),
      minTrkPt_(iConfig.getParameter<double>("minTrkPt")),
      maxTrkEta_(iConfig.getParameter<double>("maxTrkEta")),
      maxTrkNum_(iConfig.getParameter<int>("maxTrkNum")),      
      momPdgId_(iConfig.getParameter<unsigned>("momPdgId")),
      genRecoDrMatch_(iConfig.getParameter<double>("genRecoDrMatch")),
      saveCutTree_(iConfig.getParameter<bool>("saveCutTree")),
      keepelectrontags_(iConfig.getParameter<bool>("keepelectrontags")),
      debug_(iConfig.getParameter<int>("debug")),
      MCType_(iConfig.getParameter<std::string>("MCType")),
      prop1_(iConfig.getParameter<edm::ParameterSet>("propM1")) {
  //  edm::ParameterSet
  //  runParameters=iConfig.getParameter<edm::ParameterSet>("RunParameters");

  if (probeSelectorNames_.size() != probeSelectorBits_.size()) {
    throw cms::Exception("ParameterError")
        << "length of probeSelectorNames and probeSelectorBits should be identical\n";
  }
}

HIUPC_Analysis_Electrons_FullAODAnalyzer::~HIUPC_Analysis_Electrons_FullAODAnalyzer() {
  // cout << "total " << trg_counter << " fires " << fire_counter << " l3"
  // << l3_counter << endl; do anything here that needs to be done at
  // desctruction time
}

bool HIUPC_Analysis_Electrons_FullAODAnalyzer::HLTaccept(const edm::Event& iEvent,
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

void HIUPC_Analysis_Electrons_FullAODAnalyzer::fillHLTmuon(const edm::Event& iEvent,
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
        if (!(fabs(foundObject.id()) == 13 || ((fabs(foundObject.id()) == 0) && (filter_name.Contains("hltL1sSingleMu") || filter_name.Contains("hltL1MuOpen") || filter_name.Contains("hltL1sSingleEG")))))
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


void HIUPC_Analysis_Electrons_FullAODAnalyzer::embedTriggerMatching(const reco::Track& mu,
                                                         std::vector<TString>& trg_filter,
                                                         std::vector<float>& trg_pt,
                                                         std::vector<float>& trg_eta,
                                                         std::vector<float>& trg_phi,
                                                         std::vector<std::string>& HLTFilters,
                                                         bool isTag,
                                                         const int& debug_ = 0) {
  for (const auto& thefilter : HLTFilters) {
    TString thefilter_tstr = TString(thefilter);
    // temporary method to tag L2 filters for dSA paths...
    bool isL2DSA =
        thefilter_tstr.BeginsWith("hltL2") && (thefilter_tstr.Contains("NoVtx") || thefilter_tstr.Contains("NoVertex"));

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
      if ((dR_tmp < matched_dr && dR_tmp < trgDRwindow_) || 
        (isL2DSA && dR_tmp < matched_dr && dR_tmp < maxdr_trk_dsa_)) {
        matched = true;
        matched_pt = trg_pt.at(itrg);
        matched_eta = trg_eta.at(itrg);
        matched_phi = trg_phi.at(itrg);
        matched_dr = dR_tmp;

        if (debug_ > 0) {
          std::cout << "embedTriggerMatching: isTag=" << isTag << "  filter=" << thefilter_tstr << "  dR=" << dR_tmp
                    << "  matched=" << matched << std::endl;
        }
      }
    }
    if (isTag) {
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

void HIUPC_Analysis_Electrons_FullAODAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;
  t3->Fill();

  prop1_.init(iSetup);
  
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
   if(saveCutTree_) t2->Fill();

   return;   
  }
  nt.CutThrough_Num = __LINE__; // 1
 

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
  if(saveCutTree_){
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
  }
  // Gen weights, sim info
  // bool simInfoIsAvailalbe = false;
  edm::Handle<edm::ValueMap<reco::MuonSimInfo>> simInfo;
  if (!iEvent.isRealData()) {
    edm::Handle<GenEventInfoProduct> genEventInfoHandle;
    iEvent.getByToken(genEventInfoToken_, genEventInfoHandle);
    nt.genWeight = genEventInfoHandle->weight();

    // simInfoIsAvailalbe = iEvent.getByToken(simInfoToken_, simInfo);
  } else {  // data
    nt.genWeight = 1.;
  }

  // Pileup information
  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken_, rhoHandle);
  nt.Rho = *rhoHandle;

  float trueNumInteractions = -1;
  int puNumInteractions = -1;
  if (isMC_) {
    edm::Handle<std::vector<PileupSummaryInfo>> PupInfo;
    iEvent.getByToken(pileupSummaryToken_, PupInfo);

    for (auto PVI : *PupInfo) {
      int BX = PVI.getBunchCrossing();
      if (BX == 0) {
        trueNumInteractions = PVI.getTrueNumInteractions();
        puNumInteractions = PVI.getPU_NumInteractions();
        continue;
      }
    }
  }

  nt.trueNumInteractions = trueNumInteractions;
  nt.puNumInteractions = puNumInteractions;
  // if(nt.event == 453817906) cout << " bad event " << endl;
  // if(nt.event == 453817906){
  // Analysis_CaloAnalyzer CaloAnalyzer;
  // if(iEvent.isRealData()){
  // CaloAnalyzer.FillZDC(nt, iEvent,RecHitsToken_);
  // }     
     
  // }
  if (debug_ > 0)
    std::cout << "New Evt " << nt.run << " event " << nt.event <<  std::endl;

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
  // check if path fired, if so save hlt muons
  if (!HLTaccept(iEvent, nt, HLTPaths_)){
     if(saveCutTree_) t2->Fill();
     return;
  } 
  nt.CutThrough_Num = __LINE__; // 2
  
   if(tracks->size() < 2){
     if(saveCutTree_) t2->Fill();
     return;
  }  
  nt.CutThrough_Num = __LINE__; // 3
 
     // add cut to max number of tracks for HIUPC events deal with exclusivity
  int nTrk = 0;
  for (const auto& track : *tracks) {
     if(track.pt() > minTrkPt_ && fabs(track.eta()) < maxTrkEta_) nTrk++;
  }
 
  if(nTrk > maxTrkNum_){
     if(saveCutTree_) t2->Fill();
     return;
  } 
  nt.CutThrough_Num = __LINE__; // 4
 
 
 

  fillHLTmuon(iEvent, nt.trg_filter, nt.trg_pt, nt.trg_eta, nt.trg_phi, tagFilters_, debug_);
  fillHLTmuon(iEvent, nt.prb_filter, nt.prb_pt, nt.prb_eta, nt.prb_phi, probeFilters_, debug_);
  

  //calo event info
  Analysis_CaloAnalyzer CaloAnalyzer;
  if(iEvent.isRealData()){
  CaloAnalyzer.FillZDC(nt, iEvent,RecHitsToken_);
  }
  CaloAnalyzer.FillCaloTowers(nt, iEvent, calotowerToken_);

  
  // gen information
//  Analysis_MuonGenAnalyzer genmu;
  std::vector<unsigned> matched_muon_idx;
  std::vector<unsigned> matched_track_idx;
  std::vector<unsigned> matched_electron_idx;
  if (!iEvent.isRealData()) {
     if(!saveCutTree_){
    if(MCType_ == "TauTau") genmu.SetInputsandFillNtuple_TauTau(nt, iEvent, genToken_);
    else if(MCType_ == "MuMu") genmu.SetInputsandFillNtuple_MuMu(nt, iEvent, genToken_);
    else if(MCType_ == "EE") genmu.SetInputsandFillNtuple_EE(nt, iEvent, genToken_);
    else genmu.SetInputsandFillNtuple_Other(nt, iEvent, genToken_, 1); // mininum pt of final state particle is 1 GeV
     }
  
    auto reco_match_gentrk1 =
        MatchReco<reco::Muon>(*muons, nt.gentrk1_eta, nt.gentrk1_phi, nt.gentrk1_charge, genRecoDrMatch_);
    auto reco_match_gentrk2 =
        MatchReco<reco::Muon>(*muons, nt.gentrk2_eta, nt.gentrk2_phi, nt.gentrk2_charge, genRecoDrMatch_);

    if (reco_match_gentrk1.first)
      matched_muon_idx.push_back(reco_match_gentrk1.second);

    if (reco_match_gentrk2.first)
      matched_muon_idx.push_back(reco_match_gentrk2.second);

    reco_match_gentrk1 =
        MatchReco<reco::GsfElectron>(*electrons, nt.gentrk1_eta, nt.gentrk1_phi, nt.gentrk1_charge, genRecoDrMatch_);
    reco_match_gentrk2 =
        MatchReco<reco::GsfElectron>(*electrons, nt.gentrk2_eta, nt.gentrk2_phi, nt.gentrk2_charge, genRecoDrMatch_);

    if (reco_match_gentrk1.first)
      matched_electron_idx.push_back(reco_match_gentrk1.second);
    if (reco_match_gentrk2.first)
      matched_electron_idx.push_back(reco_match_gentrk2.second);
   
   
  }

  // match hlt with offline muon
  std::vector<unsigned> trg_idx;
  for (unsigned itrg = 0; itrg < nt.trg_pt.size(); ++itrg) {
    float minDR = 1000;
    unsigned idx = 0;
    for (auto& mu : *muons) {
      if (minDR < deltaR(nt.trg_eta[itrg], nt.trg_phi[itrg], mu.eta(), mu.phi()))
        continue;
      minDR = deltaR(nt.trg_eta[itrg], nt.trg_phi[itrg], mu.eta(), mu.phi());
      idx = &mu - &muons->at(0);
    }
    if (debug_ > 0)
      std::cout << "Trg " << itrg << ", min DR " << minDR << std::endl;
    if (minDR < trgDRwindow_)
      trg_idx.push_back(idx);
    if (minDR < trgDRwindow_ && debug_ > 0)
      std::cout << "Matched!" << std::endl;
  }
  nt.nMu = muons->size();


  // match hlt with offline electrons
  std::vector<unsigned> electron_trg_idx;
  for (unsigned itrg = 0; itrg < nt.trg_pt.size(); ++itrg) {
    float minDR = 1000;
    unsigned idx = 0;
    for (auto& ele : *electrons) {
      if (minDR < deltaR(nt.trg_eta[itrg], nt.trg_phi[itrg], ele.eta(), ele.phi()))
        continue;
      minDR = deltaR(nt.trg_eta[itrg], nt.trg_phi[itrg], ele.eta(), ele.phi());
      idx = &ele - &electrons->at(0);
    }
    if (debug_ > 0)
      std::cout << "Trg " << itrg << ", min DR " << minDR << std::endl;
    if (minDR < trgDRwindow_)
      electron_trg_idx.push_back(idx);
    if (minDR < trgDRwindow_ && debug_ > 0)
      std::cout << "Matched! Electron" << std::endl;
  }
  nt.nEle = electrons->size();
  // select muon tags
  int tag_muon_index = 0;
  std::vector<unsigned> tag_muon_map;  // idx of tag muon in muons
  RecoTrkAndTransientTrkCollection tag_trkttrk;
  std::vector<bool> genmatched_tag;
  
  
  for (const auto& mu : *muons) {
    tag_muon_index++;
    if (mu.selectors() != 0) {  // Only 9_4_X and later have selector bits
      if (!mu.passed(pow(2, tagQual_)))
        continue;
    } else {  // For 2016, assume loose ID on the tag (can be tightened at spark level)
      if (!muon::isLooseMuon(mu))
        continue;
    }
    if (!tagSelection_(mu))
      continue;
    if (std::find(trg_idx.begin(), trg_idx.end(), &mu - &muons->at(0)) == trg_idx.end())
      continue;
    tag_trkttrk.emplace_back(std::make_pair(mu, reco::TransientTrack(*mu.bestTrack(), &(*bField))));
    tag_muon_map.push_back(&mu - &muons->at(0));
    if (std::find(matched_muon_idx.begin(), matched_muon_idx.end(), &mu - &muons->at(0)) != matched_muon_idx.end())
      genmatched_tag.push_back(true);
    else
      genmatched_tag.push_back(false);
  }
  
   // select muon tags
  int tag_electron_index = 0;
  std::vector<unsigned> tag_electron_map;  // idx of tag muon in muons
  RecoElectronTrkAndTransientTrkCollection tag_electron_trkttrk;
  std::vector<bool> genmatched_electron_tag;
  
  
  for (const auto& ele : *electrons) {
    tag_electron_index++;
    // if (ele.selectors() != 0) {  // Only 9_4_X and later have selector bits
      // if (!.passed(pow(2, tagQual_)))
        // continue;
    // }
    // else {  // For 2016, assume loose ID on the tag (can be tightened at spark level)
      // if (!muon::isLooseMuon(mu))
        // continue;
    // }
    if (!tagElectronSelection_(ele))
      continue;
    if (std::find(electron_trg_idx.begin(), electron_trg_idx.end(), &ele - &electrons->at(0)) == electron_trg_idx.end())
      continue;
    tag_electron_trkttrk.emplace_back(std::make_pair(ele, reco::TransientTrack(*ele.bestTrack(), &(*bField))));
    tag_electron_map.push_back(&ele - &electrons->at(0));
    if (std::find(matched_electron_idx.begin(), matched_electron_idx.end(), &ele - &electrons->at(0)) != matched_electron_idx.end())
      genmatched_electron_tag.push_back(true);
    else
      genmatched_electron_tag.push_back(false);
  }
   
  
  
  
  nt.ntag = tag_trkttrk.size();
  nt.ntag_electron = tag_electron_trkttrk.size();
  
  
  if(nt.ntag == 0 && (nt.ntag_electron == 0 || !keepelectrontags_)){
     if(saveCutTree_) t2->Fill();
     return;
  } 
  nt.CutThrough_Num = __LINE__; // 5
  
  
  if (debug_ > 0)
    std::cout << "Tag muons " << tag_trkttrk.size() << std::endl;

  // probe track mapping with muon object
  std::pair<std::vector<unsigned>, std::vector<unsigned>> trk_muon_map;
  for (const auto& mu : *muons) {
    if (muonOnly_ && !probeMuonSelection_(mu))
      continue;
    float minDR = 1000;
    unsigned int idx_trk;
    if (debug_ > 1)
      std::cout << "New trk-muon map entry pt " << mu.pt() << " eta " << mu.eta() << " phi " << mu.phi() << std::endl;
    for (const reco::Track& trk : *tracks) {
      if (mu.charge() != trk.charge())
        continue;
      if (fabs(mu.vz() - trk.vz()) > maxdz_trk_mu_ && maxdz_trk_mu_ > 0)
        continue;
      if (fabs(mu.pt() - trk.pt()) / mu.pt() > maxpt_relative_dif_trk_mu_ && maxpt_relative_dif_trk_mu_ > 0)
        continue;
      float DR = deltaR(mu.eta(), mu.phi(), trk.eta(), trk.phi());
      if (debug_ > 1)
        std::cout << "   DR " << DR << "  " << mu.eta() << "  " << mu.phi() << "  " << trk.eta() << "  " << trk.phi()
                  << std::endl;
      if (minDR < DR)
        continue;
      minDR = DR;
      idx_trk = &trk - &tracks->at(0);
    }
    if (minDR > maxdr_trk_mu_)
      continue;
    trk_muon_map.first.push_back(idx_trk);
    trk_muon_map.second.push_back(&mu - &muons->at(0));
  }
  if (debug_ > 0)
    std::cout << "Matched trk-mu " << trk_muon_map.first.size() << std::endl;
 
 
   // probe track mapping with electron object
  std::pair<std::vector<unsigned>, std::vector<unsigned>> trk_electron_map;
  for (const auto& ele : *electrons) {

    float minDR = 1000;
    unsigned int idx_trk;
    if (debug_ > 1)
      std::cout << "New trk-electron map entry pt " << ele.pt() << " eta " << ele.eta() << " phi " << ele.phi() << std::endl;
    for (const reco::Track& trk : *tracks) {
      if (ele.charge() != trk.charge())
        continue;
      if (fabs(ele.vz() - trk.vz()) > maxdz_trk_ele_ && maxdz_trk_ele_ > 0)
        continue;
      if (fabs(ele.pt() - trk.pt()) / ele.pt() > maxpt_relative_dif_trk_ele_ && maxpt_relative_dif_trk_ele_ > 0)
        continue;
      float DR = deltaR(ele.eta(), ele.phi(), trk.eta(), trk.phi());
      if (debug_ > 1)
        std::cout << "   DR " << DR << "  " << ele.eta() << "  " << ele.phi() << "  " << trk.eta() << "  " << trk.phi()
                  << std::endl;
      if (minDR < DR)
        continue;
      minDR = DR;
      idx_trk = &trk - &tracks->at(0);
    }
    if (minDR > maxdr_trk_mu_)
      continue;
    trk_electron_map.first.push_back(idx_trk);
    trk_electron_map.second.push_back(&ele - &electrons->at(0));
  }
  if (debug_ > 0)
    std::cout << "Matched trk-electron " << trk_electron_map.first.size() << std::endl;

   // probe track mapping with pf cands
  std::pair<std::vector<unsigned>, std::vector<unsigned>> trk_pfc_map;
  for (const reco::Track& trk : *tracks) {
    float minDR = 1000;
    unsigned int idx_pfc;
    if (debug_ > 1)
      std::cout << "New trk-pf cand map entry pt " << trk.pt() << " eta " << trk.eta() << " phi " << trk.phi() << std::endl;
    for (const auto &pfc : *pfcands) {
      if (pfc.charge() != trk.charge())
        continue;
      if (fabs(pfc.vz() - trk.vz()) > maxdz_trk_pfc_ && maxdz_trk_pfc_ > 0)
        continue;
      if (fabs(pfc.pt() - trk.pt()) / pfc.pt() > maxpt_relative_dif_trk_pfc_ && maxpt_relative_dif_trk_pfc_ > 0)
        continue;
      float DR = deltaR(pfc.eta(), pfc.phi(), trk.eta(), trk.phi());
      if (debug_ > 1)
        std::cout << " PFC TRK  DR " << DR << "  " << pfc.eta() << "  " << pfc.phi() << "  " << trk.eta() << "  " << trk.phi()
                  << std::endl;
      if (minDR < DR)
        continue;
      minDR = DR;
      idx_pfc = &pfc - &pfcands->at(0);
    }
    if (minDR > maxdr_trk_pfc_)
      continue;
    trk_pfc_map.first.push_back(&trk - &tracks->at(0));
    trk_pfc_map.second.push_back(idx_pfc);
  }
  if (debug_ > 0)
    std::cout << "Matched trk-pfc " << trk_pfc_map.first.size() << std::endl; 


  using map_type = std::pair<std::vector<unsigned>, std::vector<unsigned>>;

  // Generic functions to map probe track with a second collection of tracks (e.g. dSA, cosmics, dGl, etc)

  // Simple mapping (not currently used) -- just check for the closest match between collections
  // If closest match has dR < config cone value, include match in mapping
  // This is not flexible for later use in spark (see Inclusive function below)
  auto mapTrackCollectionsSimple =
      [&](const std::vector<reco::Track>& coll_1, const std::vector<reco::Track>& coll_2, map_type& coll_map) {
        for (const reco::Track& trk_1 : coll_1) {
          unsigned idx_1 = &trk_1 - &coll_1.at(0);
          float min_dR = 1000;
          unsigned min_idx_2;
          for (const reco::Track& trk_2 : coll_2) {
            unsigned idx_2 = &trk_2 - &coll_2.at(0);
            float dR = deltaR(trk_1.eta(), trk_1.phi(), trk_2.eta(), trk_2.phi());
            if (min_dR < dR)
              continue;
            min_dR = dR;
            min_idx_2 = idx_2;
          }
          if (min_dR > maxdr_trk_dsa_)
            continue;
          coll_map.first.push_back(min_idx_2);
          coll_map.second.push_back(idx_1);
        }
      };
  (void)mapTrackCollectionsSimple;  // prevent "unused method" error

  // Inclusive mapping (currently used for all) -- check and count all possible matches within a config dR cone value
  // Also pick the pair with smallest dR, *even* if it is greater than config dR cone -- this allows to customize match criteria in spark later on
  // Also this should work even for collimated muons where nmatched may be greater than 1
  auto mapTrackCollectionsInclusive = [&]<typename T>(const vector<T>& coll_1,
                                                      const vector<reco::Track>& coll_2,
                                                      map_type& coll_map,
                                                      vector<float>& min_dRs,
                                                      vector<unsigned>& n_matched) {
    // Start with probe at first loop level and dSA/cosmic/dGl at second level instead, to ensure each probe that has a match is actually recorded
    // i.e. a dSA/cosmic/dGl track can be matched to more than one generalTrack probem, but all generalTrack probes that match should be pushed to the map
    for (const T& trk_1 : coll_1) {
      unsigned idx_1 = &trk_1 - &coll_1.at(0);
      unsigned matched = 0;
      float min_dR = 1000;
      unsigned min_idx_2 = 1000;
      for (const reco::Track& trk_2 : coll_2) {
        unsigned idx_2 = &trk_2 - &coll_2.at(0);
        float dR = deltaR(trk_1.eta(), trk_1.phi(), trk_2.eta(), trk_2.phi());
        // first add to count of matched items in coll 2 within X dR cone of track in coll 1
        if (dR < maxdr_trk_dsa_)
          matched++;
        // now check if this is the closest one so far
        if (dR < min_dR) {
          min_dR = dR;
          min_idx_2 = idx_2;
        }
      }
      if (min_idx_2 == 1000)  // corner case: coll 2 is empty
        continue;
      coll_map.first.push_back(idx_1);
      coll_map.second.push_back(min_idx_2);
      min_dRs.push_back(min_dR);
      n_matched.push_back(matched);
    }
  };

  // // Map probe and dGl
  // map_type probe_dGl_map;
  // vector<float> probe_dGl_dRs;
  // vector<unsigned> probe_dGl_segmentmatches;
  // // Don't use dR matching anymore for dGl, use segment matching below
  // //  mapTrackCollectionsInclusive(*tracks, *dGlmuons, probe_dGl_map, probe_dGl_dRs, probe_dGl_nmatched);

  // Map probe and cosmics
  map_type probe_cosmic_map;
  vector<float> probe_cosmic_dRs;
  vector<unsigned> probe_cosmic_nmatched;
  mapTrackCollectionsInclusive(*tracks, *staCosmic, probe_cosmic_map, probe_cosmic_dRs, probe_cosmic_nmatched);

  // // Map probe and dSA
  // map_type probe_dSA_map;
  // vector<float> probe_dSA_dRs;
  // vector<unsigned> probe_dSA_segmentmatches;
  // // Don't use dR matching anymore for dSA, use segment matching below
  // // mapTrackCollectionsInclusive(*tracks, *dSAmuons, probe_dSA_map, probe_dSA_dRs, probe_dSA_nmatched);

  // // Map tag and dSA
  // map_type tag_dSA_map;
  // vector<float> tag_dSA_dRs;
  // vector<unsigned> tag_dSA_segmentmatches;
  // // Don't use dR matching anymore for dSA, use segment matching below
  // // mapTrackCollectionsInclusive(*muons, *dSAmuons, tag_dSA_map, tag_dSA_dRs, tag_dSA_nmatched);

//dont use segment matching for analysis rn

  // // [Adapted from displaced dimuon analysis]
  // // Segment match: consider probe reco::Muons and dSA tracks
  // // matched if the segments used to build the dSA are the
  // // same as or a subset of segments of the reco::Muons.
  // for (const auto& track : *tracks) {
    // // can only do segment matching on probe tracks that are matched to muons
    // auto it = std::find(trk_muon_map.first.begin(), trk_muon_map.first.end(), &track - &tracks->at(0));
    // if (it == trk_muon_map.first.end())
      // continue;
    // unsigned idx_map = std::distance(trk_muon_map.first.begin(), it);
    // unsigned idx_track = trk_muon_map.first[idx_map];
    // unsigned idx_muon = trk_muon_map.second[idx_map];
    // auto muon = muons->at(idx_muon);
    // // also make sure probes are arbitrated tracker muons
    // if (!(muon.isTrackerMuon() && muon::isGoodMuon(muon, muon::TrackerMuonArbitrated)))
      // continue;

    // int max_nmatches = -1;
    // float min_dR = +99.f;
    // unsigned matched_idx;
    // for (const auto& dsa : *dSAmuons) {
      // unsigned idx_dsa = &dsa - &dSAmuons->at(0);

      // float dR = deltaR(dsa.eta(), dsa.phi(), muon.eta(), muon.phi());
      // // Don't waste time with far away dSA tracks and muons

      // // Aug. 2021: comment out for now to test dSA-dGl correspondence
      // //      if (dR > 0.7)
      // //        continue;

      // int nmatches = 0;
      // for (auto& hit : dsa.recHits()) {
        // if (!hit->isValid())
          // continue;
        // DetId id = hit->geographicalId();
        // if (id.det() != DetId::Muon)
          // continue;
        // if (id.subdetId() == MuonSubdetId::DT || id.subdetId() == MuonSubdetId::CSC) {
          // for (auto& chamber : muon.matches()) {
            // if (chamber.id.rawId() != id.rawId())
              // continue;
            // for (auto& segment : chamber.segmentMatches) {
              // if (fabs(segment.x - hit->localPosition().x()) < 1e-6 &&
                  // fabs(segment.y - hit->localPosition().y()) < 1e-6) {
                // if (debug_)
                  // std::cout << "matched segment found!!! subdet "
                            // << ((chamber.id.subdetId() == MuonSubdetId::DT) ? "DT" : "CSC")
                            // << " id = " << chamber.id.rawId() << " x = " << segment.x << " y = " << segment.y
                            // << std::endl;
                // nmatches++;
                // break;
              // }
            // }
          // }
        // }
      // }
      // if (nmatches > max_nmatches) {
        // max_nmatches = nmatches;
        // min_dR = dR;
        // matched_idx = idx_dsa;
      // } else if (nmatches == max_nmatches && dR < min_dR) {
        // min_dR = dR;
        // matched_idx = idx_dsa;
      // }
    // }
    // if (max_nmatches > -1) {
      // probe_dSA_map.first.push_back(idx_track);
      // probe_dSA_map.second.push_back(matched_idx);
      // probe_dSA_dRs.push_back(min_dR);
      // probe_dSA_segmentmatches.push_back(max_nmatches);
    // }
  // }

  // // [Adapted from displaced dimuon analysis]
  // // Segment match: consider probe reco::Muons and dGl tracks
  // // matched if the segments used to build the dGl are the
  // // same as or a subset of segments of the reco::Muons.
  // for (const auto& track : *tracks) {
    // // can only do segment matching on probe tracks that are matched to muons
    // auto it = std::find(trk_muon_map.first.begin(), trk_muon_map.first.end(), &track - &tracks->at(0));
    // if (it == trk_muon_map.first.end())
      // continue;
    // unsigned idx_map = std::distance(trk_muon_map.first.begin(), it);
    // unsigned idx_track = trk_muon_map.first[idx_map];
    // unsigned idx_muon = trk_muon_map.second[idx_map];
    // auto muon = muons->at(idx_muon);
    // // also make sure probes are arbitrated tracker muons
    // if (!(muon.isTrackerMuon() && muon::isGoodMuon(muon, muon::TrackerMuonArbitrated)))
      // continue;

    // int max_nmatches = -1;
    // float min_dR = +99.f;
    // unsigned matched_idx;
    // for (const auto& dgl : *dGlmuons) {
      // unsigned idx_dgl = &dgl - &dGlmuons->at(0);

      // float dR = deltaR(dgl.eta(), dgl.phi(), muon.eta(), muon.phi());
      // // Don't waste time with far away dGl tracks and muons
      // // Aug. 2021: comment out for now to test dSA-dGl correspondence
      // // if (dR > 0.7)
      // // continue;

      // int nmatches = 0;
      // for (auto& hit : dgl.recHits()) {
        // if (!hit->isValid())
          // continue;
        // DetId id = hit->geographicalId();
        // if (id.det() != DetId::Muon)
          // continue;
        // if (id.subdetId() == MuonSubdetId::DT || id.subdetId() == MuonSubdetId::CSC) {
          // for (auto& chamber : muon.matches()) {
            // if (chamber.id.rawId() != id.rawId())
              // continue;
            // for (auto& segment : chamber.segmentMatches) {
              // if (fabs(segment.x - hit->localPosition().x()) < 1e-6 &&
                  // fabs(segment.y - hit->localPosition().y()) < 1e-6) {
                // if (debug_)
                  // std::cout << "matched segment found!!! subdet "
                            // << ((chamber.id.subdetId() == MuonSubdetId::DT) ? "DT" : "CSC")
                            // << " id = " << chamber.id.rawId() << " x = " << segment.x << " y = " << segment.y
                            // << std::endl;
                // nmatches++;
                // break;
              // }
            // }
          // }
        // }
      // }
      // if (nmatches > max_nmatches) {
        // max_nmatches = nmatches;
        // min_dR = dR;
        // matched_idx = idx_dgl;
      // } else if (nmatches == max_nmatches && dR < min_dR) {
        // min_dR = dR;
        // matched_idx = idx_dgl;
      // }
    // }
    // if (max_nmatches > -1) {
      // probe_dGl_map.first.push_back(idx_track);
      // probe_dGl_map.second.push_back(matched_idx);
      // probe_dGl_dRs.push_back(min_dR);
      // probe_dGl_segmentmatches.push_back(max_nmatches);
    // }
  // }

  // // [Adapted from displaced dimuon analysis]

  // // Segment match: consider tag reco::Muons and dSA tracks
  // // matched if the segments used to build the dSA are the
  // // same as or a subset of segments of the reco::Muons.
  // for (const auto& muon : *muons) {
    // unsigned idx_muon = &muon - &muons->at(0);

    // int max_nmatches = -1;
    // float min_dR = +99.f;
    // unsigned matched_idx;
    // for (const auto& dsa : *dSAmuons) {
      // unsigned idx_dsa = &dsa - &dSAmuons->at(0);

      // float dR = deltaR(dsa.eta(), dsa.phi(), muon.eta(), muon.phi());
      // // Don't waste time with far away dSA tracks and muons
      // // Aug. 2021: comment out for now to test dSA-dGl correspondence
      // // if (dR > 0.7)
      // // continue;

      // int nmatches = 0;
      // for (auto& hit : dsa.recHits()) {
        // if (!hit->isValid())
          // continue;
        // DetId id = hit->geographicalId();
        // if (id.det() != DetId::Muon)
          // continue;
        // if (id.subdetId() == MuonSubdetId::DT || id.subdetId() == MuonSubdetId::CSC) {
          // for (auto& chamber : muon.matches()) {
            // if (chamber.id.rawId() != id.rawId())
              // continue;
            // for (auto& segment : chamber.segmentMatches) {
              // if (fabs(segment.x - hit->localPosition().x()) < 1e-6 &&
                  // fabs(segment.y - hit->localPosition().y()) < 1e-6) {
                // if (debug_)
                  // std::cout << "matched segment found!!! subdet "
                            // << ((chamber.id.subdetId() == MuonSubdetId::DT) ? "DT" : "CSC")
                            // << " id = " << chamber.id.rawId() << " x = " << segment.x << " y = " << segment.y
                            // << std::endl;
                // nmatches++;
                // break;
              // }
            // }
          // }
        // }
      // }
      // if (nmatches > max_nmatches) {
        // max_nmatches = nmatches;
        // min_dR = dR;
        // matched_idx = idx_dsa;
      // } else if (nmatches == max_nmatches && dR < min_dR) {
        // min_dR = dR;
        // matched_idx = idx_dsa;
      // }
    // }
    // if (max_nmatches > -1) {
      // tag_dSA_map.first.push_back(idx_muon);
      // tag_dSA_map.second.push_back(matched_idx);
      // tag_dSA_dRs.push_back(min_dR);
      // tag_dSA_segmentmatches.push_back(max_nmatches);
    // }
  // }


// dont use jets for analysis
  // // Muon collection for jet cleaning
  // std::vector<reco::Muon> muForJetCleaning;
  // for (const auto& mu : *muons) {
    // if (!muon::isLooseMuon(mu))
      // continue;
    // muForJetCleaning.push_back(mu);
  // }
  // sort(muForJetCleaning.begin(), muForJetCleaning.end(), [](const auto& l, const auto& r) { return l.pt() > r.pt(); });

  // std::vector<reco::PFJet> corrJets;
  // if (includeJets_) {
    // // Get selected jets and fill branches
    // for (size_t i = 0; i < jets->size(); ++i) {
      // reco::PFJetRef jet(jets, i);
      // if (CrossClean(*jet, muForJetCleaning))
        // continue;
      // std::unique_ptr<reco::PFJet> corrJet(jet->clone());
      // double jec = jetCorrector->correction(*jet);
      // corrJet->scaleEnergy(jec);
      // double smearFactor = 1.0;
      // if (isMC_) {
        // // Gen Jet Info
        // edm::Handle<std::vector<reco::GenJet>> genJets;
        // iEvent.getByToken(genJetsToken_, genJets);
        // for (const auto& genJet : *genJets) {
          // nt.genJets_pt.push_back(genJet.pt());
          // nt.genJets_eta.push_back(genJet.eta());
          // nt.genJets_phi.push_back(genJet.phi());
          // nt.genJets_mass.push_back(genJet.mass());
        // }

        // double jet_resolution = resolution.getResolution({{JME::Binning::JetPt, corrJet->pt()},
                                                          // {JME::Binning::JetEta, corrJet->eta()},
                                                          // {JME::Binning::Rho, *rhoHandle}});
        // double jer_sf = resolution_sf.getScaleFactor({{JME::Binning::JetPt, corrJet->pt()},
                                                      // {JME::Binning::JetEta, corrJet->eta()},
                                                      // {JME::Binning::Rho, *rhoHandle}},
                                                     // Variation::NOMINAL);
        // // gen matching
        // double min_dR = std::numeric_limits<double>::infinity();
        // const reco::GenJet* matched_genJet = nullptr;
        // for (const auto& genJet : *genJets) {
          // double dR = deltaR(genJet, *corrJet);
          // if (dR > min_dR)
            // continue;
          // if (dR >= 0.2)
            // continue;
          // min_dR = dR;
          // matched_genJet = &genJet;
        // }
        // if (matched_genJet) {
          // double dPt = corrJet->pt() - matched_genJet->pt();
          // smearFactor = 1 + (jer_sf - 1.) * dPt / corrJet->pt();
        // } else if (jer_sf > 1) {
          // double sigma = jet_resolution * std::sqrt(jer_sf * jer_sf - 1);
          // std::normal_distribution<> d(0, sigma);
          // smearFactor = 1. + d(m_random_generator);
        // }
        // if (corrJet->pt() * smearFactor < 0.01) {
          // smearFactor = 0.01 / corrJet->energy();
        // }
      // }
      // corrJet->scaleEnergy(smearFactor);
      // corrJets.push_back(*corrJet);
      // FillJetBranches(*jet, *corrJet, nt, era_);
      // if (deepCSVProbb->size() > i && deepCSVProbbb->size() > i) {
        // nt.jets_bTag_deepCSV.push_back((*deepCSVProbb)[i].second + (*deepCSVProbbb)[i].second);
      // } else
        // nt.jets_bTag_deepCSV.push_back(-9999.);
    // }
  // }



  // StandAloneFillEventInfo(nt, StandAlone_nt, HLTPaths_);
  
  
  // // run over tracks and probes once prior to filling tree to determine ordering of pairs
  // // this is necessary to use tag-probe pair with highest "quality" later on in spark_tnp.
  // // Also calculate the total number of pairs per event to include in ntuple
  using t_pair_prob = std::pair<float, std::pair<int, int>>;
  std::priority_queue<t_pair_prob> pair_vtx_probs;
  std::priority_queue<t_pair_prob> pair_dPhi_muons;
  std::priority_queue<t_pair_prob> pair_Mass_Mmumu;
  // std::priority_queue<t_pair_prob, vector<t_pair_prob>, std::greater<t_pair_prob>> pair_dz_PV_SV;    // inverse sort
  // std::priority_queue<t_pair_prob, vector<t_pair_prob>, std::greater<t_pair_prob>> pair_Mass_Mmumu;  // inverse sort

  // std::priority_queue<t_pair_prob> SA_pair_vtx_probs;
  // std::priority_queue<t_pair_prob> SA_pair_dPhi_muons;
  // std::priority_queue<t_pair_prob, vector<t_pair_prob>, std::greater<t_pair_prob>> SA_pair_dz_PV_SV;    // inverse sort
  // std::priority_queue<t_pair_prob, vector<t_pair_prob>, std::greater<t_pair_prob>> SA_pair_dM_Z_Mmumu;  // inverse sort
  // // loop over tags
  
  int TnPSize = 0;
  int tag_idx = 0;
  int tag_pt = 0;
  // int tag_electron_idx = 0;
  // int tag_electron_pt = 0;
  std::string tagtype = "none";
    for (auto& tag : tag_trkttrk) {
       if( tag.first.pt() < tag_pt) continue;
       tag_pt = tag.first.pt();
       tag_idx = &tag - &tag_trkttrk[0];
       tagtype = "muon";

    
    }
    
    if(tagtype != "muon"  && keepelectrontags_){
    for (auto& tag : tag_electron_trkttrk) {
       if( tag.first.pt() < tag_pt) continue;
       tag_pt = tag.first.pt();
       tag_idx = &tag - &tag_electron_trkttrk[0];
       tagtype = "electron";
    }
    }
    
      reco::Muon tag_muon;
      reco::GsfElectron tag_electron;
      reco::TransientTrack tag_transtrk;
      reco::Track tag_track;
      if(tagtype == "muon"){
         auto tag = tag_trkttrk[tag_idx];
         tag_muon = tag.first;
         tag_transtrk = tag.second;
         tag_track = *tag_muon.bestTrack();
         nt.tag_isElectron = false;
         nt.tag_isMuon = true;
      }
      
      if(tagtype == "electron"){
         auto tag = tag_electron_trkttrk[tag_idx];
         tag_electron = tag.first;
         tag_transtrk = tag.second;
         tag_track = *tag_electron.bestTrack();
         nt.tag_isElectron = true;
         nt.tag_isMuon = false;
      }
    if(tagtype == "none" ){
       if(saveCutTree_) t2->Fill();
       return;
    }
       nt.CutThrough_Num = __LINE__; // 6
    std::pair<reco::Track, reco::TransientTrack> tag = {tag_track, tag_transtrk};      
      
    for (const reco::Track& probe : *tracks) {
       // auto tag = tag_trkttrk[tag_idx];
      auto probe_idx = &probe - &tracks->at(0);
      if((tag.first.charge() == probe.charge()) && (deltaR(tag.first.eta(), tag.first.phi(), probe.eta(),probe.phi()) < .0001))
         continue;
      // apply cuts on probe
      if (HighPurity_ && !probe.quality(Track::highPurity))
        continue;
      if (!probeSelection_(probe))
        continue;

      // apply cuts on pairs
      if (tag.first.charge() == probe.charge())
        continue;
      if (fabs(tag.first.vz() - probe.vz()) > pairDz_ && pairDz_ > 0)
        continue;
      float mass = DimuonMass(tag.first.pt(), tag.first.eta(), tag.first.phi(), probe.pt(), probe.eta(), probe.phi());
      if (mass < pairMassMin_ || mass > pairMassMax_)
        continue;
  
      // // compute vtx
      std::vector<reco::TransientTrack> trk_pair = {tag.second, reco::TransientTrack(probe, &(*bField))};
      Analysis_KlFitter vtx(trk_pair);
      // KalmanVertexFitter fitter;
      // TransientVertex myVertex = fitter.vertex(trk_pair);
      // GlobalPoint vert(myVertex.position().x(), myVertex.position().y(), myVertex.position().z());
      // TrajectoryStateClosestToPoint  traj = trk_pair.at(0).trajectoryStateClosestToPoint(vert);
      // cout << "Test Vtx " << traj.perigeeParameters().transverseImpactParameter() << endl;
      // if (RequireVtxCreation_ && !vtx.status())
        // continue;
      // if (minSVtxProb_ > 0 && vtx.prob() < minSVtxProb_)
        // continue;

      float dPhi_muons = reco::deltaPhi(tag.first.phi(), probe.phi());
      math::PtEtaPhiMLorentzVector mu1(tag.first.pt(), tag.first.eta(), tag.first.phi(), MU_MASS);
      math::PtEtaPhiMLorentzVector mu2(probe.pt(), probe.eta(), probe.phi(), MU_MASS);
      float Mass_Mmumu = abs((mu1 + mu2).mass());

      // // save quantities to ordered heap
      auto pair_idx = std::make_pair(tag_idx, probe_idx);
      pair_vtx_probs.push(std::make_pair(vtx.prob(), pair_idx));
      // //commented out b/c crashed with HI data file
// //      pair_dz_PV_SV.push(std::make_pair(vtx.dz_PV_SV(nt.pv_z), pair_idx));
      pair_dPhi_muons.push(std::make_pair(dPhi_muons, pair_idx));
      pair_Mass_Mmumu.push(std::make_pair(Mass_Mmumu, pair_idx));
    TnPSize++;
    }
    
  nt.npairs = TnPSize;


  // assign sorted vtx indices to ranking
  map<std::pair<int, int>, int> pair_rank_vtx_prob;
  while (!pair_vtx_probs.empty()) {
    pair_rank_vtx_prob[pair_vtx_probs.top().second] = pair_rank_vtx_prob.size();  // careful: RHS evaluated first
    pair_vtx_probs.pop();
  }
  // // map<std::pair<int, int>, int> pair_rank_dz_PV_SV;
  // // while (!pair_dz_PV_SV.empty()) {
    // // pair_rank_dz_PV_SV[pair_dz_PV_SV.top().second] = pair_rank_dz_PV_SV.size();  // careful: RHS evaluated first
    // // pair_dz_PV_SV.pop();
  // // }
  map<std::pair<int, int>, int> pair_rank_dPhi_muons;
  while (!pair_dPhi_muons.empty()) {
    pair_rank_dPhi_muons[pair_dPhi_muons.top().second] = pair_rank_dPhi_muons.size();  // careful: RHS evaluated first
    pair_dPhi_muons.pop();
  }
  map<std::pair<int, int>, int> pair_rank_Mass_Mmumu;
  while (!pair_Mass_Mmumu.empty()) {
    pair_rank_Mass_Mmumu[pair_Mass_Mmumu.top().second] = pair_rank_Mass_Mmumu.size();  // careful: RHS evaluated first
    pair_Mass_Mmumu.pop();
  }

    // loop over probe tracks
      for (const reco::Track& probe : *tracks) {
      nt.CutThrough_Num = __LINE__; //7
         
      // auto tag = tag_trkttrk[tag_idx];
      if((tag.first.charge() == probe.charge()) && (deltaR(tag.first.eta(), tag.first.phi(), probe.eta(),probe.phi()) < .0001)){
         continue;
      }
      if (debug_ > 0)
      std::cout << "New tag pt " << tag.first.pt() << " eta " << tag.first.eta() << " phi " << tag.first.phi()
                << std::endl;
      if (debug_ > 1)
          std::cout << "    Probe pt " << probe.pt() << " eta " << probe.eta() << " phi " << probe.phi() << "  charge "
                    << probe.charge() << std::endl;

        // apply cuts on probe
        if (HighPurity_ && !probe.quality(Track::highPurity)){
          if(saveCutTree_) t2->Fill();
          continue;
          }
        nt.CutThrough_Num++;  //8
        if (!probeSelection_(probe)){
          if(saveCutTree_) t2->Fill();
          continue;
          }
        nt.CutThrough_Num = __LINE__; //9
        // apply cuts on pairs; selected will be saved
        if (tag.first.charge() == probe.charge()){
          if(saveCutTree_) t2->Fill();
          continue;
          }
        nt.CutThrough_Num = __LINE__; //9
        if (fabs(tag.first.vz() - probe.vz()) > pairDz_ && pairDz_ > 0){
          if(saveCutTree_) t2->Fill();
          continue;
          }
        nt.CutThrough_Num = __LINE__; //10

        float mass = DimuonMass(tag.first.pt(), tag.first.eta(), tag.first.phi(), probe.pt(), probe.eta(), probe.phi());
        if (mass < pairMassMin_ || mass > pairMassMax_){
          if(saveCutTree_) t2->Fill();
          continue;
          }
        nt.CutThrough_Num = __LINE__; // 11
        std::vector<reco::TransientTrack> trk_pair = {tag.second, reco::TransientTrack(probe, &(*bField))};
        Analysis_KlFitter vtx(trk_pair);
        // if (RequireVtxCreation_ && !vtx.status())
          // continue;
        // if (minSVtxProb_ > 0 && vtx.prob() < minSVtxProb_)
          // continue;
        int probe_trk_idx = &probe - &tracks->at(0);
        auto it = std::find(trk_muon_map.first.begin(), trk_muon_map.first.end(), probe_trk_idx);
        auto it_pfc = std::find(trk_pfc_map.first.begin(), trk_pfc_map.first.end(), probe_trk_idx);
        auto it_electron = std::find(trk_electron_map.first.begin(), trk_electron_map.first.end(), probe_trk_idx);

        if (muonOnly_ && it == trk_muon_map.first.end()){
          if(saveCutTree_) t2->Fill();
          continue;
          }
        nt.CutThrough_Num = __LINE__; // 12
       

        if(!iEvent.isRealData()){
           if(tag.first.charge() < 0 ){
                nt.gentrk1_match_dr = deltaR(tag.first.eta(), tag.first.phi(), nt.gentrk1_eta, nt.gentrk1_phi);
                nt.gentrk1_match_dphi = fabs(tag.first.phi()-nt.gentrk1_phi);
                nt.gentrk1_match_diff_eta= fabs(tag.first.eta()-nt.gentrk1_eta);
                nt.gentrk1_match_diff_pt = fabs(tag.first.pt()-nt.gentrk1_pt); 
                nt.gentrk1_diff_vtx_x = fabs(nt.pv_x-nt.gentrk1_vtx_x);
                nt.gentrk1_diff_vtx_y= fabs(nt.pv_y-nt.gentrk1_vtx_y);
                nt.gentrk1_diff_vtx_z= fabs(nt.pv_z-nt.gentrk1_vtx_z); 
                nt.gentrk1_isTag = true;
  
                nt.gentrk2_match_dr = deltaR(probe.eta(), probe.phi(), nt.gentrk2_eta, nt.gentrk2_phi);
                nt.gentrk2_match_dphi = fabs(probe.phi()-nt.gentrk2_phi);
                nt.gentrk2_match_diff_eta= fabs(probe.eta()-nt.gentrk2_eta);
                nt.gentrk2_match_diff_pt = fabs(probe.pt()-nt.gentrk2_pt); 
                nt.gentrk2_diff_vtx_x = fabs(nt.pv_x-nt.gentrk2_vtx_x);
                nt.gentrk2_diff_vtx_y= fabs(nt.pv_y-nt.gentrk2_vtx_y);
                nt.gentrk2_diff_vtx_z= fabs(nt.pv_z-nt.gentrk2_vtx_z); 
                nt.gentrk2_isTag = false;
           }
           else{
                nt.gentrk1_match_dr = deltaR(probe.eta(), probe.phi(), nt.gentrk1_eta, nt.gentrk1_phi);
                nt.gentrk1_match_dphi = fabs(probe.phi()-nt.gentrk1_phi);
                nt.gentrk1_match_diff_eta= fabs(probe.eta()-nt.gentrk1_eta);
                nt.gentrk1_match_diff_pt = fabs(probe.pt()-nt.gentrk1_pt); 
                nt.gentrk1_diff_vtx_x = fabs(nt.pv_x-nt.gentrk1_vtx_x);
                nt.gentrk1_diff_vtx_y= fabs(nt.pv_y-nt.gentrk1_vtx_y);
                nt.gentrk1_diff_vtx_z= fabs(nt.pv_z-nt.gentrk1_vtx_z); 
                nt.gentrk1_isTag = false;
  
                nt.gentrk2_match_dr = deltaR(tag.first.eta(), tag.first.phi(), nt.gentrk2_eta, nt.gentrk2_phi);
                nt.gentrk2_match_dphi = fabs(tag.first.phi()-nt.gentrk2_phi);
                nt.gentrk2_match_diff_eta= fabs(tag.first.eta()-nt.gentrk2_eta);
                nt.gentrk2_match_diff_pt = fabs(tag.first.pt()-nt.gentrk2_pt); 
                nt.gentrk2_diff_vtx_x = fabs(nt.pv_x-nt.gentrk2_vtx_x);
                nt.gentrk2_diff_vtx_y= fabs(nt.pv_y-nt.gentrk2_vtx_y);
                nt.gentrk2_diff_vtx_z= fabs(nt.pv_z-nt.gentrk2_vtx_z); 
                nt.gentrk2_isTag = true;
           }
        }

        int probe_muon_idx = -99 ;
        int tag_trk_idx = -99;
        int tag_muon_idx = -99;
        int tag_electron_idx = -99;
        if( tagtype  =="muon") tag_muon_idx = tag_muon_map[tag_idx];
        if( tagtype  =="electron") tag_electron_idx = tag_electron_map[tag_idx]; 
        
        nt.probe_vtx_x =probe.vx();
        nt.probe_vtx_y =probe.vy();
        nt.probe_vtx_z =probe.vz();
        
        if (it != trk_muon_map.first.end()){
           unsigned tmp_idx = std::distance(trk_muon_map.first.begin(), it);
           probe_muon_idx = trk_muon_map.second[tmp_idx];
           nt.probe_hasMuonMatch = true;
           nt.probe_MuonMatchDR = deltaR(probe.eta(), probe.phi(), muons->at(probe_muon_idx).eta(), muons->at(probe_muon_idx).phi());
        }
        if(tagtype == "muon"){
        auto it_tag_trk = std::find(trk_muon_map.second.begin(), trk_muon_map.second.end(), tag_muon_idx);
        if (it_tag_trk != trk_muon_map.second.end()){
           unsigned tmp_idx = std::distance(trk_muon_map.second.begin(), it_tag_trk);
           tag_trk_idx = trk_muon_map.first[tmp_idx];
           nt.tag_hasTrackMatch = true;
           nt.tag_TrackMatchDR = deltaR(tag_muon.eta(), tag_muon.phi(), tracks->at(tag_trk_idx).eta(), tracks->at(tag_trk_idx).phi());
        }
        }
        if(tagtype == "electron"){
        auto it_tag_trk = std::find(trk_electron_map.second.begin(), trk_electron_map.second.end(), tag_electron_idx);
        if (it_tag_trk != trk_electron_map.second.end()){
           unsigned tmp_idx = std::distance(trk_electron_map.second.begin(), it_tag_trk);
           tag_trk_idx = trk_electron_map.first[tmp_idx];
           nt.tag_hasTrackMatch = true;
           nt.tag_TrackMatchDR = deltaR(tag_electron.eta(), tag_electron.phi(), tracks->at(tag_trk_idx).eta(), tracks->at(tag_trk_idx).phi());
        }
        }
        int probe_pfc_idx = -99;
        if (it_pfc != trk_pfc_map.first.end()){
           unsigned tmp_idx = std::distance(trk_pfc_map.first.begin(), it_pfc);
           probe_pfc_idx = trk_pfc_map.second[tmp_idx];
           nt.probe_hasPFCMatch = true;
           nt.probe_pfcID = pfcands->at(probe_pfc_idx).pdgId();
           nt.probe_PFCMatchDR = deltaR(probe.eta(), probe.phi(), pfcands->at(probe_pfc_idx).eta(), pfcands->at(probe_pfc_idx).phi());
        }
        int probe_electron_idx = -99;
        if (it_electron != trk_electron_map.first.end()){
           unsigned tmp_idx = std::distance(trk_electron_map.first.begin(), it_electron);
           probe_electron_idx = trk_electron_map.second[tmp_idx];
           nt.probe_hasElectronMatch = true;
           nt.probe_ElectronMatchDR = deltaR(probe.eta(), probe.phi(), electrons->at(probe_electron_idx).eta(), electrons->at(probe_electron_idx).phi());
        }

        FillMuonBranches<reco::Muon>(*muons, nt,tag_muon_idx, probe_muon_idx, *pv);
        FillElectronBranches<reco::GsfElectron>(*electrons, nt, tag_electron_idx, probe_electron_idx, *pv);
        FillPFCandBranches<reco::PFCandidate>(*pfcands, nt, probe_pfc_idx);
        FillTrackBranches<reco::Track>(*tracks, nt, tag_trk_idx,probe_trk_idx);
        
        
        math::PtEtaPhiMLorentzVector mu1(tag.first.pt(), tag.first.eta(), tag.first.phi(), MU_MASS);
        math::PtEtaPhiMLorentzVector mu2(probe.pt(), probe.eta(), probe.phi(), MU_MASS);

        if(tagtype == "muon" )nt.tag_isMatchedGen = genmatched_tag[tag_idx];
        if(tagtype == "electron" )nt.tag_isMatchedGen = genmatched_electron_tag[tag_electron_idx];        
                                                                                    

        if(tagtype == "muon" ) FillTagBranches<reco::Muon, reco::Track>(tag_muon, *tracks, nt, *pv);
        if(tagtype == "electron" ) FillTagElectronBranches<reco::GsfElectron,reco::Track>(tag_electron, *tracks, nt, *pv);
        if(tagtype == "muon" ) FillMiniIso<reco::Muon, reco::PFCandidate>(*pfcands, tag_muon, *rhoJetsNC, nt, true);
        // Tag-trigger matching
        if(tagtype == "muon"){
        auto tagRef = muonsView->refAt(tag_muon_map[tag_idx]);
        pat::TriggerObjectStandAloneRef tagl1Match = (*l1Matches)[tagRef];
        if (tagl1Match.isNonnull()) {
          nt.tag_l1pt = tagl1Match->pt();
          nt.tag_l1q = (*l1Qualities)[tagRef];
          nt.tag_l1dr = (*l1Drs)[tagRef];
        } else {
          nt.tag_l1pt = -99.;
          nt.tag_l1q = -99;
          nt.tag_l1dr = 99.;
        }

        pat::TriggerObjectStandAloneRef tagl1MatchByQ = (*l1MatchesByQ)[tagRef];
        if (tagl1MatchByQ.isNonnull()) {
          nt.tag_l1ptByQ = tagl1MatchByQ->pt();
          nt.tag_l1qByQ = (*l1QualitiesByQ)[tagRef];
          nt.tag_l1drByQ = (*l1DrsByQ)[tagRef];
        } else {
          nt.tag_l1ptByQ = -99.;
          nt.tag_l1qByQ = -99;
          nt.tag_l1drByQ = 99.;
        }
        }                                
        embedTriggerMatching(tag.first, nt.trg_filter, nt.trg_pt, nt.trg_eta, nt.trg_phi, tagFilters_, true, debug_);

        // if (!simInfoIsAvailalbe) {
          // FillSimMatchingBranchesDummy(nt, true);
        // } else {
          // const auto& msi = (*simInfo)[tagRef];
          // FillSimMatchingBranchesAOD(msi, nt, true);
        // }

        // auto itdsa = std::find(probe_dSA_map.first.begin(), probe_dSA_map.first.end(), &probe - &tracks->at(0));
        // auto itdsa_tag = std::find(tag_dSA_map.first.begin(), tag_dSA_map.first.end(), tag_idx);
        // auto itdgl = std::find(probe_dGl_map.first.begin(), probe_dGl_map.first.end(), &probe - &tracks->at(0));
        auto itcosmic =
            std::find(probe_cosmic_map.first.begin(), probe_cosmic_map.first.end(), &probe - &tracks->at(0));

        if (it == trk_muon_map.first.end()) {
          if (debug_ > 0)
            std::cout << "  Unsuccessful probe " << std::endl;
          reco::Muon fakeMuon;
          fakeMuon.setP4(mu2);
          fakeMuon.setCharge(probe.charge());
          FillProbeBranches<reco::Muon, reco::Track>(fakeMuon, *tracks, nt, false, *pv);
          FillProbeBranchesSelector<reco::Muon>(fakeMuon, nt, probeSelectorBits_, false);
          FillMiniIso<reco::Muon, reco::PFCandidate>(*pfcands, fakeMuon, *rhoJetsNC, nt, false);
          // if (includeJets_)
            // FindJetProbePair<reco::PFJet, reco::Muon>(corrJets, fakeMuon, nt);

          // store dummy trigger variables if offline muon is not found
          for (const auto& path : probeFilters_) {
            nt.probe_trg[&path - &probeFilters_[0]] = false;
            nt.probe_trg_pt[&path - &probeFilters_[0]] = -99;
            nt.probe_trg_eta[&path - &probeFilters_[0]] = -99;
            nt.probe_trg_phi[&path - &probeFilters_[0]] = -99;
            nt.probe_trg_dr[&path - &probeFilters_[0]] = 99;
          }
          nt.l1pt = -99.;
          nt.l1q = -99;
          nt.l1dr = 99.;
          nt.l1ptByQ = -99.;
          nt.l1qByQ = -99;
          nt.l1drByQ = 99.;

          // FillSimMatchingBranchesDummy(nt, false);
          FillTunePPairBranchesDummy(nt);
          
        } else {
          unsigned idx = std::distance(trk_muon_map.first.begin(), it);
          if (debug_ > 0)
            std::cout << "  Successful probe pt " << muons->at(trk_muon_map.second[idx]).pt() << " eta "
                      << muons->at(trk_muon_map.second[idx]).eta() << " phi "
                      << muons->at(trk_muon_map.second[idx]).phi() << std::endl;
          FillProbeBranches<reco::Muon, reco::Track>(muons->at(trk_muon_map.second[idx]), *tracks, nt, true, *pv);
          FillProbeBranchesSelector<reco::Muon>(muons->at(trk_muon_map.second[idx]), nt, probeSelectorBits_, true);
          FillMiniIso<reco::Muon, reco::PFCandidate>(
              *pfcands, muons->at(trk_muon_map.second[idx]), *rhoJetsNC, nt, false);
          // if (includeJets_)
            // FindJetProbePair<reco::PFJet, pat::Muon>(corrJets, muons->at(trk_muon_map.second[idx]), nt);

          // Probe-trigger matching
          auto muRef = muonsView->refAt(trk_muon_map.second[idx]);
          pat::TriggerObjectStandAloneRef l1Match = (*l1Matches)[muRef];
          if (l1Match.isNonnull()) {
            nt.l1pt = l1Match->pt();
            nt.l1q = (*l1Qualities)[muRef];
            nt.l1dr = (*l1Drs)[muRef];
          } else {
            nt.l1pt = -99.;
            nt.l1q = -99;
            nt.l1dr = 99.;
          }

          pat::TriggerObjectStandAloneRef l1MatchByQ = (*l1MatchesByQ)[muRef];
          if (l1MatchByQ.isNonnull()) {
            nt.l1ptByQ = l1MatchByQ->pt();
            nt.l1qByQ = (*l1QualitiesByQ)[muRef];
            nt.l1drByQ = (*l1DrsByQ)[muRef];
          } else {
            nt.l1ptByQ = -99.;
            nt.l1qByQ = -99;
            nt.l1drByQ = 99.;
          }

          embedTriggerMatching(probe,
                               nt.prb_filter,
                               nt.prb_pt,
                               nt.prb_eta,
                               nt.prb_phi,
                               probeFilters_,
                               false,
                               debug_);

          // if (!simInfoIsAvailalbe) {
            // FillSimMatchingBranchesDummy(nt, false);
          // } else {
            // const auto& msi = (*simInfo)[muRef];
            // FillSimMatchingBranchesAOD(msi, nt, false);
          // }

          // TuneP pair branches
          if(tagtype == "muon"){
          if (tag_muon.tunePMuonBestTrack().isNonnull() &&
              muons->at(trk_muon_map.second[idx]).tunePMuonBestTrack().isNonnull()) {
            const reco::TrackRef tag_tuneP = tag_muon.tunePMuonBestTrack();
            const reco::TrackRef probe_tuneP = muons->at(trk_muon_map.second[idx]).tunePMuonBestTrack();
            FillTunePPairBranches<reco::Track, reco::Track>(*tag_tuneP, *probe_tuneP, nt);
            std::vector<reco::TransientTrack> ttrk_pair_tuneP = {reco::TransientTrack(*tag_tuneP, &(*bField)),
                                                                 reco::TransientTrack(*probe_tuneP, &(*bField))};
            Analysis_KlFitter vtx_tuneP(ttrk_pair_tuneP);
            vtx_tuneP.fillNtuple(nt, true);
          } else {
            FillTunePPairBranchesDummy(nt);
          }
        }
        }        

        // if (itdsa == probe_dSA_map.first.end()) {
          // nt.probe_dsa_segmentMatches = -1;
          // FillProbeBranchesdSA<reco::Track>(probe, nt, false);
        // } else {
          // unsigned idx = std::distance(probe_dSA_map.first.begin(), itdsa);
          // nt.probe_dsa_segmentMatches = probe_dSA_segmentmatches[idx];
          // nt.probe_dsa_minDR = probe_dSA_dRs[idx];
          // if (debug_ > 0)
            // std::cout << "Successful probe dSA " << dSAmuons->at(probe_dSA_map.second[idx]).pt() << " eta "
                      // << dSAmuons->at(probe_dSA_map.second[idx]).eta() << " phi "
                      // << dSAmuons->at(probe_dSA_map.second[idx]).phi() << std::endl;
          // FillProbeBranchesdSA<reco::Track>(dSAmuons->at(probe_dSA_map.second[idx]), nt, true);
        // }

        // if (itdsa_tag == tag_dSA_map.first.end()) {
          // nt.tag_dsa_segmentMatches = -1;
          // FillTagBranchesdSA<reco::Track>(probe, nt, false);
        // } else {
          // unsigned idx = std::distance(tag_dSA_map.first.begin(), itdsa_tag);
          // nt.tag_dsa_segmentMatches = tag_dSA_segmentmatches[idx];
          // nt.tag_dsa_minDR = tag_dSA_dRs[idx];
          // if (debug_ > 0)
            // std::cout << "Successful tag dSA " << dSAmuons->at(tag_dSA_map.second[idx]).pt() << " eta "
                      // << dSAmuons->at(tag_dSA_map.second[idx]).eta() << " phi "
                      // << dSAmuons->at(tag_dSA_map.second[idx]).phi() << std::endl;
          // FillTagBranchesdSA<reco::Track>(dSAmuons->at(tag_dSA_map.second[idx]), nt, true);
        // }

        // if (itdgl == probe_dGl_map.first.end()) {
          // nt.probe_dgl_segmentMatches = -1;
          // FillProbeBranchesdgl<reco::Track>(probe, nt, false);
        // } else {
          // unsigned idx = std::distance(probe_dGl_map.first.begin(), itdgl);
          // nt.probe_dgl_segmentMatches = probe_dGl_segmentmatches[idx];
          // nt.probe_dgl_minDR = probe_dGl_dRs[idx];
          // if (debug_ > 0)
            // std::cout << "Successful probe displaced global " << dGlmuons->at(probe_dGl_map.second[idx]).pt() << " eta "
                      // << dGlmuons->at(probe_dGl_map.second[idx]).eta() << " phi "
                      // << dGlmuons->at(probe_dGl_map.second[idx]).phi() << std::endl;
          // FillProbeBranchesdgl<reco::Track>(dGlmuons->at(probe_dGl_map.second[idx]), nt, true);
        // }

        if (itcosmic == probe_cosmic_map.first.end()) {
          nt.probe_ncosmic = 0;
          FillProbeBranchesCosmic<reco::Track>(probe, nt, false);
        } else {
          unsigned idx = std::distance(probe_cosmic_map.first.begin(), itcosmic);
          nt.probe_ncosmic = probe_cosmic_nmatched[idx];
          nt.probe_cosmic_minDR = probe_cosmic_dRs[idx];
          if (debug_ > 0)
            std::cout << "Successful probe cosmic " << staCosmic->at(probe_cosmic_map.second[idx]).pt() << " eta "
                      << staCosmic->at(probe_cosmic_map.second[idx]).eta() << " phi "
                      << staCosmic->at(probe_cosmic_map.second[idx]).phi() << std::endl;
          FillProbeBranchesCosmic<reco::Track>(staCosmic->at(probe_cosmic_map.second[idx]), nt, true);
        }

        FillPairBranches<reco::Track, reco::Track>(tag.first, probe, nt, prop1_);
        vtx.fillNtuple(nt);

        auto it_genmatch = std::find(matched_track_idx.begin(), matched_track_idx.end(), &probe - &tracks->at(0));
        nt.probe_isMatchedGen = (it_genmatch != matched_track_idx.end());

        nt.iprobe++;
        nt.pair_rank_vtx_prob = pair_rank_vtx_prob[{tag_idx, &probe - &tracks->at(0)}];
        // //commented out b/c HI data causes crash with filling this 
        // // nt.pair_rank_dz_PV_SV = pair_rank_dz_PV_SV[{tag_idx, &probe - &tracks->at(0)}];
        nt.pair_rank_dPhi_muons = pair_rank_dPhi_muons[{tag_idx, &probe - &tracks->at(0)}];
        nt.pair_rank_Mass_Mmumu = pair_rank_Mass_Mmumu[{tag_idx, &probe - &tracks->at(0)}];
        nt.probe_isHighPurity = probe.quality(Track::highPurity);

        t1->Fill();
        if(saveCutTree_) t2->Fill();
      }
    
  
}

// ------------ method called once each job just before starting event loop
// ------------
void HIUPC_Analysis_Electrons_FullAODAnalyzer::beginJob() {
  t1 = fs->make<TTree>("Events", "Events");
  t2 = fs->make<TTree>("GenVtxStudy","GenVtxStudy");
  t3 = fs->make<TTree>("test","test");
  nt.SetTree(t1);
  if(saveCutTree_) nt.SetTree_GenVtxStudy(t2);
  nt.SetTree_Test(t3); 
  nt.CreateBranches(HLTPaths_, probeSelectorNames_);
  if(saveCutTree_) nt.CreateBranches_GenVtxStudy();
  if (!tagFilters_.empty()) {
    nt.CreateExtraTrgBranches(tagFilters_, true);
    // StandAlone_nt.CreateExtraTrgBranches(tagFilters_, true);
  }
  if (!probeFilters_.empty())
    nt.CreateExtraTrgBranches(probeFilters_, false);
}

// ------------ method called once each job just after ending the event loop
// ------------
void HIUPC_Analysis_Electrons_FullAODAnalyzer::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the
// module  ------------
void HIUPC_Analysis_Electrons_FullAODAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(HIUPC_Analysis_Electrons_FullAODAnalyzer);
