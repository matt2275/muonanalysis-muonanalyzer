#ifndef MuonAnalysis_MuonAnalyzer_plugins_Analysis_MuonGenAnalyzer
#define MuonAnalysis_MuonAnalyzer_plugins_Analysis_MuonGenAnalyzer

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
#include "Analysis_NtupleContent.h"
#include "TLorentzVector.h"
#include "helper.h"

class Analysis_MuonGenAnalyzer {
public:
  Analysis_MuonGenAnalyzer();
  virtual ~Analysis_MuonGenAnalyzer();

  void SetInputsandFillNtuple_MuMu(Analysis_NtupleContent &, const edm::Event &, const edm::EDGetTokenT<edm::View<reco::GenParticle>> &);
  void SetInputsandFillNtuple_EE(Analysis_NtupleContent &, const edm::Event &, const edm::EDGetTokenT<edm::View<reco::GenParticle>> &);
  void SetInputsandFillNtuple_TauTau(Analysis_NtupleContent &, const edm::Event &, const edm::EDGetTokenT<edm::View<reco::GenParticle>> &);
  void SetInputsandFillNtuple_Other(Analysis_NtupleContent &, const edm::Event &, const edm::EDGetTokenT<edm::View<reco::GenParticle>> &, const float &);
  // void FillNtuple(NtupleContent &);

private:
  edm::Handle<edm::View<reco::GenParticle>> gens;
  TLorentzVector gmuon1, gmuon2;
  int index1,index2;
  // bool success = true;
};

#endif
