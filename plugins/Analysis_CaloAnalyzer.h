#ifndef MuonAnalysis_MuonAnalyzer_plugins_Analysis_CaloAnalyzer
#define MuonAnalysis_MuonAnalyzer_plugins_Analysis_CaloAnalyzer

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"


#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"


#include "Analysis_NtupleContent.h"
#include "TLorentzVector.h"
#include "helper.h"

class Analysis_CaloAnalyzer {
public:
  Analysis_CaloAnalyzer();
  virtual ~Analysis_CaloAnalyzer();

  void FillZDC(Analysis_NtupleContent &, const edm::Event &, const edm::EDGetTokenT<edm::SortedCollection<ZDCRecHit>> &);
  void FillCaloTowers(Analysis_NtupleContent &, const edm::Event &, const edm::EDGetTokenT<edm::SortedCollection<CaloTower>> &);
  // void FillNtuple(NtupleContent &);


private:
  edm::Handle<edm::SortedCollection<CaloTower>> calotowercollection;
  edm::Handle<edm::SortedCollection<ZDCRecHit>> zdcrechits;

  const double maxEtaEB = 1.479;
  const double minEtaEE = maxEtaEB;
  const double maxEtaEE = 3.0;

  const double maxEtaHB = 1.305;
  const double minEtaHE = maxEtaHB;
  const double maxEtaHE = 3.0;

  const double minEtaHF = 2.9;
  const double maxEtaHF = 5.2;
  
  
  TLorentzVector gmuon1, gmuon2;
  unsigned int index1,index2;
  // bool success = true;
};

#endif
