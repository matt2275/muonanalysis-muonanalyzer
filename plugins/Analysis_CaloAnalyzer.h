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
#include "DataFormats/EgammaCandidates/interface/Photon.h"


#include "Analysis_NtupleContent.h"
#include "TLorentzVector.h"
#include "helper.h"
#include <sstream>
#include <map>
using namespace std;
class Analysis_CaloAnalyzer {
public:
  Analysis_CaloAnalyzer();
  virtual ~Analysis_CaloAnalyzer();
  // int  eta2ieta(double eta);
  void FillPhotons(Analysis_NtupleContent &, const edm::Event &, const edm::EDGetTokenT<edm::View<reco::Photon>> &);
  void FillZDC(Analysis_NtupleContent &, const edm::Event &, const edm::EDGetTokenT<edm::SortedCollection<ZDCRecHit>> &);
  void FillCaloTowers(Analysis_NtupleContent &, const edm::Event &, const edm::EDGetTokenT<edm::SortedCollection<CaloTower>> &);
  // void FillNtuple(NtupleContent &);


private:
  edm::Handle<edm::SortedCollection<CaloTower>> calotowercollection;
  edm::Handle<edm::SortedCollection<ZDCRecHit>> zdcrechits;
  edm::Handle<edm::View<reco::Photon> > recoPhotonsHandle;  

  const double maxEtaEB = 1.479;
  const double minEtaEE = maxEtaEB;
  const double maxEtaEE = 3.0;

  const double maxEtaHB = 1.305;
  const double minEtaHE = maxEtaHB;
  const double maxEtaHE = 3.0;

  const double minEtaHF = 2.9;
  const double maxEtaHF = 5.2;
  
  bool doNoiseEEetaDependant = true;
  
std::map<std::string, float> config_params = {


// configs_params.insert(pair<string, float>("noiseThresholdEB", 0.7));
{"noiseThresholdEB", 0.7},
{"noiseThresholdEE", 7.5},
{"noiseThresholdHB", 2.8},
{"noiseThresholdHE", 2.4},
{"noiseThresholdHFp", 7.2},
{"noiseThresholdHFm", 7.5},

{"noiseEEetaStep", 0.1},
{"noiseEEetaMin", 1.5},
{"noiseEEetaMax", 3.0},

{"noiseThresholdEE_1.5", 9.61},
{"noiseThresholdEE_1.6", 1.5},
{"noiseThresholdEE_1.7", 1.2},
{"noiseThresholdEE_1.8", 1.5},
{"noiseThresholdEE_1.9", 1.5},
{"noiseThresholdEE_2.0", 1.8},
{"noiseThresholdEE_2.1", 2.11},
{"noiseThresholdEE_2.2", 4.41},
{"noiseThresholdEE_2.3", 7.91},
{"noiseThresholdEE_2.4", 9.21},
{"noiseThresholdEE_2.5", 9.71},
{"noiseThresholdEE_2.6", 9.91},
{"noiseThresholdEE_2.7", 10.01},
{"noiseThresholdEE_2.8", 10.01},
{"noiseThresholdEE_2.9", 10.01},
{"noiseThresholdEE_3.0", 10.01}
};  
  
  TLorentzVector gmuon1, gmuon2;
  unsigned int index1,index2;
  // bool success = true;
template <typename T>
string to_string_with_precision(const T a_value, const int n = 6)
{
  ostringstream out;
  out.precision(n);
  out << fixed << a_value;
  return out.str();
}

 int ietaMax = 42;

 std::array<double, 42> etaedge = {
      {0.000, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.870, 0.957, 1.044, 1.131,
       1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.740, 1.830, 1.930, 2.043, 2.172, 2.322, 2.500, 2.650,
       2.853, 3.000, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191}};

int eta2ieta(double eta) {
  // binary search in the array of towers eta edges

  int ieta = 1;
  double xeta = fabs(eta);
  while (xeta > etaedge[ieta] && ieta < ietaMax - 1) {
    ++ieta;
  }

  if (eta < 0)
    ieta = -ieta;
  return ieta;
}


};

#endif
