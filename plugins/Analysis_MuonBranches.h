//
// Original Author:
//                george karathanasis
//         Created:  Thu, 23 Mar 2019 17:40:23 GMT
//
// filling functions for aod and miniaod tag/probe

#ifndef MuonAnalysis_MuonAnalyzer_MuonBranches
#define MuonAnalysis_MuonAnalyzer_MuonBranches

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/MuonSimInfo.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include <type_traits>
#include "Analysis_NtupleContent.h"
#include "helper.h"

template <typename MUON, typename TRK>
inline void FillTagBranches(const MUON &muon,
                            const std::vector<TRK> &tracks,
                            Analysis_NtupleContent &nt,
                            const reco::Vertex &vertex) {
  nt.tag_pt = muon.pt();
  nt.tag_eta = muon.eta();
  nt.tag_phi = muon.phi();
  nt.tag_charge = muon.charge();
  nt.tag_pterr = muon.innerTrack()->ptError() / muon.innerTrack()->pt();
  nt.tag_dxy = muon.innerTrack()->dxy(reco::TrackBase::Point(nt.pv_x, nt.pv_y, nt.pv_z));
  nt.tag_dz = muon.innerTrack()->dz(reco::TrackBase::Point(nt.pv_x, nt.pv_y, nt.pv_z));
  nt.tag_vtx_x = muon.vx();
  nt.tag_vtx_y = muon.vy();
  nt.tag_vtx_z = muon.vz();
  nt.tag_isPF = muon.isPFMuon();
  nt.tag_isSA = muon.isStandAloneMuon();
  nt.tag_isTracker = muon.isTrackerMuon();
  nt.tag_isGlobal = muon.isGlobalMuon();
  // Use selectors instead of 'muon.passed' method which is only introduced in CMSSW_9_4_X
  nt.tag_isLoose = muon::isLooseMuon(muon);
  nt.tag_isMedium = muon::isMediumMuon(muon);
  nt.tag_isTight = muon::isTightMuon(muon, vertex);
  nt.tag_isSoft = muon::isSoftMuon(muon, vertex, false);
  nt.tag_isHighPt = muon::isHighPtMuon(muon, vertex);
  float iso04 = (TrackerEnergy04<TRK>(muon.eta(), muon.phi(), tracks) - muon.pt()) / muon.pt();
  nt.tag_relIso04 = (iso04 > 0) ? iso04 : 0;
  nt.tag_iso03_sumPt = muon.isolationR03().sumPt;
  nt.tag_pfIso04_charged = muon.pfIsolationR04().sumChargedHadronPt;
  nt.tag_pfIso04_neutral = muon.pfIsolationR04().sumNeutralHadronEt;
  nt.tag_pfIso04_photon = muon.pfIsolationR04().sumPhotonEt;
  nt.tag_pfIso04_sumPU = muon.pfIsolationR04().sumPUPt;
  if (muon.tunePMuonBestTrack().isNonnull()) {
    nt.tag_tuneP_pt = muon.tunePMuonBestTrack()->pt();
    nt.tag_tuneP_pterr = muon.tunePMuonBestTrack()->ptError();
  } else {
    nt.tag_tuneP_pt = -99.;
    nt.tag_tuneP_pterr = -99.;
  }
  int nsegments = 0;
  for (auto &chamber : muon.matches()) {
    if (chamber.id.det() != DetId::Muon)
      continue;
    if (chamber.id.subdetId() != MuonSubdetId::DT && chamber.id.subdetId() != MuonSubdetId::CSC)
      continue;
    nsegments += chamber.segmentMatches.size();
  }
  nt.tag_nsegments = nsegments;
}

template <typename MUON>
inline void FillMuonBranches(
    const std::vector<MUON> &muons, Analysis_NtupleContent &nt, int &tag_idx, int &probe_idx,const reco::Vertex &pv){

       nt.nMu = muons.size();

       for (const auto& mu : muons) {
          // nt.mu_pT.push_back(mu.pt());
          // nt.mu_eta.push_back(mu.eta());
          // nt.mu_phi.push_back(mu.phi()); 
          // nt.mu_charge.push_back(mu.charge()); 
          
          
          nt.muPt.push_back(mu.pt());
          nt.muEta.push_back(mu.eta());
          nt.muPhi.push_back(mu.phi());
          nt.muCharge.push_back(mu.charge());
          nt.muType.push_back(mu.type());
          nt.muIsGood.push_back(muon::isGoodMuon(mu, muon::selectionTypeFromString("TMOneStationTight")));
          
          nt.muIsGlobal.push_back((int)mu.isGlobalMuon());
          nt.muIsTracker.push_back((int)mu.isTrackerMuon());
          nt.muIsPF.push_back((int)mu.isPFMuon());
          nt.muIsSTA.push_back((int)mu.isStandAloneMuon());

          nt.muD0.push_back(mu.muonBestTrack()->dxy(pv.position()));
          nt.muDz.push_back(mu.muonBestTrack()->dz(pv.position()));
          nt.muD0Err.push_back(mu.muonBestTrack()->dxyError());
          nt.muDzErr.push_back(mu.muonBestTrack()->dzError());

          // //Initialize with nonphysical values
          // float muIP3D = -999;
          // float muIP3DErr = -999;
          // if (pv.isValid()) {
            // //3DImpact parameter
            // reco::TransientTrack tt = tb->build(mu.muonBestTrack().get());
            // nt.muIP3D = IPTools::absoluteImpactParameter3D(tt, pv).second.value();
            // nt.muIP3DErr = IPTools::absoluteImpactParameter3D(tt, pv).second.error();
          // }
          // nt.muIP3D.push_back(muIP3D);
          // nt.muIP3DErr.push_back(muIP3DErr);

          const reco::TrackRef glbMu = mu.globalTrack();
          const reco::TrackRef innMu = mu.innerTrack();

          if (glbMu.isNull()) {
            nt.muChi2NDF.push_back(-99);
            nt.muMuonHits.push_back(-99);
          } else {
            nt.muChi2NDF.push_back(glbMu->normalizedChi2());
            nt.muMuonHits.push_back(glbMu->hitPattern().numberOfValidMuonHits());
          }

          if (innMu.isNull()) {
            nt.muInnerD0.push_back(-99);
            nt.muInnerDz.push_back(-99);
            
            nt.muInnerD0Err.push_back(-99);
            nt.muInnerDzErr.push_back(-99);
            nt.muInnerPt.push_back(-99);
            nt.muInnerPtErr.push_back(-99);
            nt.muInnerEta.push_back(-99);
           
            
            nt.muTrkLayers.push_back(-99);
            nt.muPixelLayers.push_back(-99);
            nt.muPixelHits.push_back(-99);
            nt.muTrkQuality.push_back(-99);
          } else {
            nt.muInnerD0.push_back(innMu->dxy(pv.position()));
            nt.muInnerDz.push_back(innMu->dz(pv.position()));
            
            nt.muInnerD0Err.push_back(innMu->dxyError());
            nt.muInnerDzErr.push_back(innMu->dzError());
            nt.muInnerPt.push_back(innMu->pt());
            nt.muInnerPtErr.push_back(innMu->ptError());
            nt.muInnerEta.push_back(innMu->eta());
            
            nt.muTrkLayers.push_back(innMu->hitPattern().trackerLayersWithMeasurement());
            nt.muPixelLayers.push_back(innMu->hitPattern().pixelLayersWithMeasurement());
            nt.muPixelHits.push_back(innMu->hitPattern().numberOfValidPixelHits());
            nt.muTrkQuality.push_back(innMu->quality(reco::TrackBase::highPurity));
          }

          nt.muStations.push_back(mu.numberOfMatchedStations());
          nt.muIsoTrk.push_back(mu.isolationR03().sumPt);
          nt.muPFChIso.push_back(mu.pfIsolationR04().sumChargedHadronPt);
          nt.muPFPhoIso.push_back(mu.pfIsolationR04().sumPhotonEt);
          nt.muPFNeuIso.push_back(mu.pfIsolationR04().sumNeutralHadronEt);
          nt.muPFPUIso.push_back(mu.pfIsolationR04().sumPUPt);
          nt.muIDSoft.push_back(mu.passed(reco::Muon::SoftMvaId));
          nt.muIDLoose.push_back(mu.passed(reco::Muon::CutBasedIdLoose));
          nt.muIDMedium.push_back(mu.passed(reco::Muon::CutBasedIdMedium));
          nt.muIDMediumPrompt.push_back(mu.passed(reco::Muon::CutBasedIdMediumPrompt));
          nt.muIDTight.push_back(mu.passed(reco::Muon::CutBasedIdTight));
          nt.muIDGlobalHighPt.push_back(mu.passed(reco::Muon::CutBasedIdGlobalHighPt));
          nt.muIDTrkHighPt.push_back(mu.passed(reco::Muon::CutBasedIdTrkHighPt));
          nt.muIDInTime.push_back(mu.passed(reco::Muon::InTimeMuon)); 

          if((&mu - &muons.at(0)) == probe_idx){
          nt.mu_isProbe.push_back(true);
          nt.mu_isTag.push_back(false);
          }
          else if((&mu - &muons.at(0)) == tag_idx){
          nt.mu_isProbe.push_back(false);
          nt.mu_isTag.push_back(true);
          }
          else{
          nt.mu_isProbe.push_back(false);
          nt.mu_isTag.push_back(false);
          }           
       
    }
    }
    
template <typename TRK>
inline void FillTrackBranches(
    const std::vector<TRK> &tracks, Analysis_NtupleContent &nt, int &tag_idx, int &probe_idx){
       nt.nTrk = tracks.size();
       for (const auto& trk : tracks) { 
          nt.trkPt.push_back(trk.pt());            
          nt.trkP .push_back(trk.p());            
          nt.trkEta.push_back(trk.eta());           
          nt.trkPhi.push_back(trk.phi());     
          nt.trkcharge.push_back(trk.charge()); 
          nt.trkvx.push_back(trk.vx());       
          nt.trkvy.push_back(trk.vy());             
          nt.trkvz.push_back(trk.vz());  
          nt.trknormchi2.push_back(trk.normalizedChi2());                   
          nt.trkchi2.push_back(trk.chi2());    
          nt.trkd0.push_back(trk.d0());               
          nt.trkdxy.push_back(trk.dxy());             
          nt.trkdz.push_back(trk.dz());
          // nt.trkdxyBS.push_back(trk.dxy(beamSpotPosition));
          // nt.trkdzBS.push_back(trk.dz(beamSpotPosition));
          nt.trkdxyError.push_back(trk.dxyError());             
          nt.trkdzError.push_back(trk.dzError());   
          nt.trkValidHits.push_back(trk.numberOfValidHits());                     
          nt.trkMissHits.push_back(trk.numberOfLostHits());  
          nt.trkPurity.push_back(trk.reco::TrackBase::qualityByName("highPurity")); 
          if((&trk - &tracks.at(0)) == probe_idx){
          nt.trkisProbe.push_back(true);
          nt.trkisTag.push_back(false);
          }
          if((&trk - &tracks.at(0)) == tag_idx){
          nt.trkisProbe.push_back(false);
          nt.trkisTag.push_back(true);
          }
          else{
          nt.trkisProbe.push_back(false);
          nt.trkisTag.push_back(false);
          }           
       }

    }
    
template <typename ELE>
inline void FillElectronBranches(
    const std::vector<ELE> &electrons, Analysis_NtupleContent &nt, int &probe_idx, const reco::Vertex &pv){

       nt.nEle = electrons.size();
       for (const auto& ele : electrons) {


          nt.eleCharge.push_back(ele.charge());
          nt.eleChargeConsistent.push_back((int)ele.isGsfCtfScPixChargeConsistent());
          nt.eleSCPixCharge.push_back(ele.scPixCharge());
          // if (!(ele.closestCtfTrack().isNull())) {
            // nt.eleCtfCharge.push_back(ele.closestCtfTrack()->charge());
          // } else {
            // nt.eleCtfCharge.push_back(-99.);
          // }
          nt.eleEn.push_back(ele.energy());
          nt.eleD0.push_back(ele.gsfTrack()->dxy(pv.position()));
          nt.eleDz.push_back(ele.gsfTrack()->dz(pv.position()));
          nt.eleD0Err.push_back(ele.gsfTrack()->dxyError());
          nt.eleDzErr.push_back(ele.gsfTrack()->dzError());
          nt.eleTrkPt.push_back(ele.gsfTrack()->pt());
          nt.eleTrkEta.push_back(ele.gsfTrack()->eta());
          nt.eleTrkPhi.push_back(ele.gsfTrack()->phi());
          nt.eleTrkCharge.push_back(ele.gsfTrack()->charge());
          nt.eleTrkPtErr.push_back(ele.gsfTrack()->ptError());
          nt.eleTrkChi2.push_back(ele.gsfTrack()->chi2());
          nt.eleTrkNdof.push_back(ele.gsfTrack()->ndof());
          nt.eleTrkNormalizedChi2.push_back(ele.gsfTrack()->normalizedChi2());
          nt.eleTrkValidHits.push_back(ele.gsfTrack()->numberOfValidHits());
          nt.eleTrkLayers.push_back(ele.gsfTrack()->hitPattern().trackerLayersWithMeasurement());
          nt.elePt.push_back(ele.pt());
          nt.eleEta.push_back(ele.eta());
          nt.elePhi.push_back(ele.phi());
          nt.eleSCEn.push_back(ele.superCluster()->energy());
          nt.eleESEn.push_back(ele.superCluster()->preshowerEnergy());
          nt.eleSCEta.push_back(ele.superCluster()->eta());
          nt.eleSCPhi.push_back(ele.superCluster()->phi());
          nt.eleSCRawEn.push_back(ele.superCluster()->rawEnergy());
          nt.eleSCEtaWidth.push_back(ele.superCluster()->etaWidth());
          nt.eleSCPhiWidth.push_back(ele.superCluster()->phiWidth());
          nt.eleHoverE.push_back(ele.hcalOverEcal());
          nt.eleHoverEBc.push_back(ele.hcalOverEcalBc());
          nt.eleEoverP.push_back(ele.eSuperClusterOverP());
          nt.eleEoverPInv.push_back(1./ele.ecalEnergy() - 1./ele.trackMomentumAtVtx().R());
          nt.eleEcalE.push_back(ele.ecalEnergy());
          nt.elePAtVtx.push_back(ele.trackMomentumAtVtx().R());
          nt.elePAtSC.push_back(ele.trackMomentumAtCalo().R());
          nt.elePAtCluster.push_back(ele.trackMomentumAtEleClus().R());
          nt.elePAtSeed.push_back(ele.trackMomentumOut().R());
          nt.eleBrem.push_back(ele.fbrem());
          nt.eledEtaAtVtx.push_back(ele.deltaEtaSuperClusterTrackAtVtx());
          nt.eledPhiAtVtx.push_back(ele.deltaPhiSuperClusterTrackAtVtx());
          nt.eledEtaSeedAtVtx.push_back(ele.deltaEtaSeedClusterTrackAtVtx());
          nt.eleSigmaIEtaIEta.push_back(ele.sigmaIetaIeta());
          nt.eleSigmaIPhiIPhi.push_back(ele.sigmaIphiIphi());
          nt.eleMissHits.push_back(ele.gsfTrack()->numberOfLostHits());  
          if((&ele - &electrons.at(0)) == probe_idx) nt.ele_isProbe.push_back(true);
          else nt.ele_isProbe.push_back(false);        
       }
    }
    
    
template <typename PFC>
inline void FillPFCandBranches(
    const std::vector<PFC> &pfcands, Analysis_NtupleContent &nt, int &probe_idx){

       nt.nPFCands = pfcands.size();
    for (const auto &pfc : pfcands) {
      if(pfc.charge() == 0) continue;
      nt.pfcand_pdgId.push_back(pfc.pdgId());
      nt.pfcand_charge.push_back(pfc.charge()); 
      nt.pfcand_pt.push_back(pfc.pt());
      nt.pfcand_eta.push_back(pfc.eta());
      nt.pfcand_phi.push_back(pfc.phi()); 
      nt.pfcand_vtx_x.push_back(pfc.vx());
      nt.pfcand_vtx_y.push_back(pfc.vy());
      nt.pfcand_vtx_z.push_back(pfc.vz());
      if((&pfc - &pfcands.at(0)) == probe_idx) nt.pfcand_isProbe.push_back(true);
      else nt.pfcand_isProbe.push_back(false);
    }
    }

template <typename MUON, typename TRK>
inline void FillProbeBranches(
    const MUON &mu, const std::vector<TRK> &tracks, Analysis_NtupleContent &nt, bool success, const reco::Vertex &vertex) {
  nt.probe_pt = mu.pt();
  nt.probe_eta = mu.eta();
  nt.probe_phi = mu.phi();
  nt.probe_charge = mu.charge();
  float iso04 = (TrackerEnergy04<TRK>(mu.eta(), mu.phi(), tracks) - mu.pt()) / mu.pt();
  nt.probe_relIso04 = (iso04 > 0) ? iso04 : 0;
}

template <typename MUON>
inline void FillProbeBranchesSelector(const MUON &mu,
                                      Analysis_NtupleContent &nt,
                                      const std::vector<unsigned> selectorBits,
                                      bool success) {
  for (unsigned int ibit = 0; ibit < selectorBits.size(); ++ibit) {
    if (success && mu.selectors() != 0) {
      unsigned bit = selectorBits.at(ibit);
      nt.probe_selectors[ibit] = mu.passed(1UL << bit);
    } else {
      nt.probe_selectors[ibit] = false;
    }
  }
}



template <typename TRK>
inline void FillProbeBranchesCosmic(const TRK &trk, Analysis_NtupleContent &nt, bool passcosmic) {
  nt.probe_isCosmic = passcosmic;
}

template <typename MUO, typename TRK>
inline void FillPairBranches(const MUO &muon, const TRK &trk, Analysis_NtupleContent &nt, PropagateToMuon &prop1_) {
  math::PtEtaPhiMLorentzVector mu1(muon.pt(), muon.eta(), muon.phi(), MU_MASS);
  math::PtEtaPhiMLorentzVector mu2(trk.pt(), trk.eta(), trk.phi(), MU_MASS);
  nt.pair_pt = (mu1 + mu2).pt();
  nt.pair_mass = (mu1 + mu2).mass();
  nt.pair_eta = (mu1 + mu2).eta();
  nt.pair_phi = (mu1 + mu2).phi();
  nt.pair_dz = muon.vz() - trk.vz();
  nt.pair_dR = deltaR(muon.eta(), muon.phi(), trk.eta(), trk.phi());
  nt.pair_dphi = std::fabs(muon.phi() - trk.phi());
  //added conditional because breaks on mini aod with prop2_M1
  // trk variable which is now pf packed candidates instead of reco::Track
  // if (is_same<TRK, reco::Track>::value) {
    // TrajectoryStateOnSurface prop1_M1 = prop1_.extrapolate(muon);
    // TrajectoryStateOnSurface prop2_M1 = prop1_.extrapolate(trk);
    // if (prop1_M1.isValid() && prop2_M1.isValid()) {
      // float dphiM1 = deltaPhi<float>(prop1_M1.globalPosition().phi(), prop2_M1.globalPosition().phi());
      // nt.pair_drM1 = hypot(dphiM1, std::abs<float>(prop1_M1.globalPosition().eta() - prop2_M1.globalPosition().eta()));
    // }
  // }
}

template <typename MUO, typename TRK>
inline void FillTunePPairBranches(const MUO &muon, const TRK &trk, Analysis_NtupleContent &nt) {
  math::PtEtaPhiMLorentzVector mu1(muon.pt(), muon.eta(), muon.phi(), MU_MASS);
  math::PtEtaPhiMLorentzVector mu2(trk.pt(), trk.eta(), trk.phi(), MU_MASS);
  nt.pair_tuneP_pt = (mu1 + mu2).pt();
  nt.pair_tuneP_mass = (mu1 + mu2).mass();
  nt.pair_tuneP_eta = (mu1 + mu2).eta();
  nt.pair_tuneP_phi = (mu1 + mu2).phi();
  nt.pair_tuneP_dz = muon.vz() - trk.vz();
  nt.pair_tuneP_dR = deltaR(muon.eta(), muon.phi(), trk.eta(), trk.phi());
  nt.pair_tuneP_dphi = std::fabs(muon.phi() - trk.phi());
}

inline void FillTunePPairBranchesDummy(Analysis_NtupleContent &nt) {
  nt.pair_tuneP_pt = -99;
  nt.pair_tuneP_mass = -99;
  nt.pair_tuneP_eta = -99;
  nt.pair_tuneP_phi = -99;
  nt.pair_tuneP_dz = -99;
  nt.pair_tuneP_dR = -99;
  nt.pair_tuneP_fit_mass = -99;
  nt.pair_tuneP_svprob = -99;
  nt.pair_tuneP_normalchi2 = -99;
}


#endif
