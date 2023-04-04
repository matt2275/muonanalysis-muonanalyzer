#include "Analysis_CaloAnalyzer.h"
#include "ConfigManager.h"

Analysis_CaloAnalyzer::Analysis_CaloAnalyzer(){};

Analysis_CaloAnalyzer::~Analysis_CaloAnalyzer(){};

// static constexpr int ietaMax = 42;

// constexpr std::array<double, 42> etaedge = {
      // {0.000, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.870, 0.957, 1.044, 1.131,
       // 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.740, 1.830, 1.930, 2.043, 2.172, 2.322, 2.500, 2.650,
       // 2.853, 3.000, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191}};

// int eta2ieta(double eta) {
  // // binary search in the array of towers eta edges

  // int ieta = 1;
  // double xeta = fabs(eta);
  // while (xeta > etaedge[ieta] && ieta < ietaMax - 1) {
    // ++ieta;
  // }

  // if (eta < 0)
    // ieta = -ieta;
  // return ieta;
// }


void Analysis_CaloAnalyzer::FillCaloTowers(Analysis_NtupleContent& nt,const edm::Event& iEvent, const edm::EDGetTokenT<edm::SortedCollection<CaloTower>>& calotowerToken_) {
  iEvent.getByToken(calotowerToken_, calotowercollection);
   
   
 int nTowers = 0;
 float maxHFp =0;
 float maxHFm =0;
// string configPath = "/afs/cern.ch/user/m/mnickel/private/MUONPDG/CMSSW_10_6_18/src/MuonAnalysis/MuonAnalyzer/efficiencies.md";
// config = ConfigManager(configPath);
for (edm::SortedCollection<CaloTower>::const_iterator calo = calotowercollection->begin(); calo != calotowercollection->end(); ++calo) {
  float Tower_eta = calo->eta();
  float Tower_energy = calo->energy();
  float Tower_energyHad = calo->hadEnergy();
  float Tower_energyEM = calo->emEnergy();

  int ieta = eta2ieta(Tower_eta);
  string DetType = "";
  string DetTypeEM = "";
  bool HasTowerAboveThresholdHF = false;
  bool HasTowerAboveThresholdHad = false;
  bool HasTowerAboveThresholdEM = false;
  bool inNoiseRegionHF = false;
  bool inNoiseRegionHad = false;
  bool inNoiseRegionEM = false;
  
  if(Tower_eta > -maxEtaHF && Tower_eta < -minEtaHF) DetType = "HFm";
  if(Tower_eta >  minEtaHF && Tower_eta <  maxEtaHF) DetType = "HFp";
   if(fabs(Tower_eta) > 0 && fabs(Tower_eta) < maxEtaHB) DetType = "HB";
   if(fabs(Tower_eta) > minEtaHE && fabs(Tower_eta) < maxEtaHE) DetType = "HE";
   
  if(Tower_eta > -maxEtaHF && Tower_eta < -minEtaHF) DetTypeEM = "HFm";
  if(Tower_eta >  minEtaHF && Tower_eta <  maxEtaHF) DetTypeEM = "HFp";
  
  if(fabs(Tower_eta) > 0 && fabs(Tower_eta) < maxEtaEB) DetTypeEM = "EB";
  if(fabs(Tower_eta) > minEtaEE && fabs(Tower_eta) < maxEtaEE) DetTypeEM = "EE"; 
  
  // if(DetType == "") continue;
  
 
    if(DetType == "HFp" && Tower_energy > maxHFp) maxHFp = Tower_energy;
    if(DetType == "HFm" && Tower_energy > maxHFm) maxHFm = Tower_energy;
    if(DetType != ""){ // Check HF exclusivity
     if(DetType == "HFm" || DetType == "HFp"){
      if(fabs(ieta) == 29 || fabs(ieta) == 30) inNoiseRegionHF = true;
      if(Tower_energy > config_params["noiseThreshold"+DetType]){
        HasTowerAboveThresholdHF = true;
      }
     }
     if(DetType == "HB" || DetType == "HE"){
      if(DetType == "HE" && (ieta == 16 || ieta == -16 )) inNoiseRegionHad = true;
      if(Tower_energyHad > config_params["noiseThreshold"+DetType]){
        HasTowerAboveThresholdHad = true;
      }
     }
    }
    if(DetTypeEM != ""){
       double fabs_Tower_eta = fabs(Tower_eta);
       double threshold = -1;
    if(ieta == 27 || ieta == 28 || ieta == 29 || ieta == -27 || ieta == -28 || ieta == -29) inNoiseRegionEM = true;
    // if(DetTypeEM=="EE" && config_params["doNoiseEEetaDependant"]){
    if(DetTypeEM=="EE" && doNoiseEEetaDependant){
    double etaStep = config_params["noiseEEetaStep"];
    double etaMin = config_params["noiseEEetaMin"];
    double etaMax = config_params["noiseEEetaMax"];
    

    
    if(fabs_Tower_eta < etaMin || fabs_Tower_eta > etaMax) threshold = 99999;
    
    for(double eta=etaMin; eta<etaMax; eta += etaStep){
      if(fabs_Tower_eta> eta && fabs_Tower_eta < eta+etaStep){
        threshold = config_params["noiseThresholdEE_"+to_string_with_precision(eta, 1)];
        
      }
    }
  }
  else{
    threshold = config_params["noiseThreshold"+DetTypeEM];
  }
  if(Tower_energyEM > threshold) HasTowerAboveThresholdEM = true;
  }
    
       if(HasTowerAboveThresholdHF || HasTowerAboveThresholdHad || HasTowerAboveThresholdEM){
       nTowers++;
       nt.CaloTower_emE.push_back(calo->emEnergy());
       nt.CaloTower_hadE.push_back(calo->hadEnergy());
       nt.CaloTower_e.push_back(calo->energy());
       nt.CaloTower_et.push_back(calo->et());
       nt.CaloTower_phi.push_back(calo->phi());
       nt.CaloTower_eta.push_back(calo->eta());
       nt.CaloTower_HF_AboveThreshold.push_back(HasTowerAboveThresholdHF);
       nt.CaloTower_Had_AboveThreshold.push_back(HasTowerAboveThresholdHad);
       nt.CaloTower_EM_AboveThreshold.push_back(HasTowerAboveThresholdEM);
       nt.CaloTower_HF_inNoiseRegion.push_back(inNoiseRegionHF);
       nt.CaloTower_Had_inNoiseRegion.push_back(inNoiseRegionHad);
       nt.CaloTower_EM_inNoiseRegion.push_back(inNoiseRegionEM);
       }
}
  nt.nTower = nTowers;
  nt.maxHFp = maxHFp;  
  nt.maxHFm = maxHFm;  
  
}


void Analysis_CaloAnalyzer::FillZDC(Analysis_NtupleContent& nt,const edm::Event& iEvent, const edm::EDGetTokenT<edm::SortedCollection<ZDCRecHit>>& RecHitsToken_) {
  iEvent.getByToken(RecHitsToken_, zdcrechits);

    nt.ZDC_n = zdcrechits->size();   
    int nhits = 0;
    for (auto const& rh : *zdcrechits) {
      HcalZDCDetId zdcid = rh.id();
      if (nhits  < 18) {
       float zdc_energy = rh.energy();
       int zdc_side = zdcid.zside();
       int zdc_section = zdcid.section();
       nt.ZDC_e[nhits] = zdc_energy;
       nt.ZDC_zside[nhits] = zdc_side;
       nt.ZDC_section[nhits] = zdc_section;
       nt.ZDC_channel[nhits] = zdcid.channel();
       nt.ZDC_PM_Total_Energy+= zdc_energy;
       if(zdc_side >0 ){
       nt.ZDC_P_Total_Energy += zdc_energy;
       if(zdc_section == 1) nt.ZDC_P_ECal_Energy += zdc_energy;
       if(zdc_section == 2) nt.ZDC_P_HCal_Energy += zdc_energy;        
       }
       if(zdc_side <0 ){
       nt.ZDC_M_Total_Energy += zdc_energy;
       if(zdc_section == 1) nt.ZDC_M_ECal_Energy += zdc_energy;
       if(zdc_section == 2) nt.ZDC_M_HCal_Energy += zdc_energy;         
       }
       // StandAlone_nt.ZDC_saturation[nhits] = static_cast<int>( rh.flagField(HcalCaloFlagLabels::ADCSaturationBit) );
      }

      nhits++;
    } // end loop zdc rechits 

}


void Analysis_CaloAnalyzer::FillPhotons(Analysis_NtupleContent& nt,const edm::Event& iEvent, const edm::EDGetTokenT<edm::View<reco::Photon>>& recoPhotonsCollection) {
  edm::Handle<edm::View<reco::Photon> > recoPhotonsHandle;
  iEvent.getByToken(recoPhotonsCollection, recoPhotonsHandle);
  
  
  for (auto pho = recoPhotonsHandle->begin(); pho != recoPhotonsHandle->end(); ++pho) {

    //if(abs(pho->eta())>2.4) continue;
    //if(abs(pho->superCluster()->eta())>2.4) continue;

    nt.phoE.push_back(pho->energy());
    nt.phoEt.push_back(pho->et());
    nt.phoEta.push_back(pho->eta());
    nt.phoPhi.push_back(pho->phi());
    
        // energies from different types of corrections
    nt.phoEcorrStdEcal.push_back(pho->getCorrectedEnergy(reco::Photon::P4type::ecal_standard));
    nt.phoEcorrPhoEcal.push_back(pho->getCorrectedEnergy(reco::Photon::P4type::ecal_photons));
    nt.phoEcorrRegr1.push_back(pho->getCorrectedEnergy(reco::Photon::P4type::regression1));
    nt.phoEcorrRegr2.push_back(pho->getCorrectedEnergy(reco::Photon::P4type::regression2));
    // errors for those corrections
    nt.phoEcorrErrStdEcal.push_back(pho->getCorrectedEnergyError(reco::Photon::P4type::ecal_standard));
    nt.phoEcorrErrPhoEcal.push_back(pho->getCorrectedEnergyError(reco::Photon::P4type::ecal_photons));
    nt.phoEcorrErrRegr1.push_back(pho->getCorrectedEnergyError(reco::Photon::P4type::regression1));
    nt.phoEcorrErrRegr2.push_back(pho->getCorrectedEnergyError(reco::Photon::P4type::regression2));

    // SuperCluster info
    nt.phoSCE.push_back(pho->superCluster()->energy());
    nt.phoSCEt.push_back(pho->superCluster()->energy()/cosh(pho->superCluster()->eta()));
    nt.phoSCRawE.push_back(pho->superCluster()->rawEnergy());
    nt.phoSCEta.push_back(pho->superCluster()->eta());
    nt.phoSCPhi.push_back(pho->superCluster()->phi());
    nt.phoSCEtaWidth.push_back(pho->superCluster()->etaWidth());
    nt.phoSCPhiWidth.push_back(pho->superCluster()->phiWidth());
    nt.phoSCBrem.push_back(pho->superCluster()->phiWidth()/pho->superCluster()->etaWidth());
    nt.phoSCnHits.push_back(pho->superCluster()->size());
    nt.phoSCflags.push_back(pho->superCluster()->flags());
    nt.phoSCinClean.push_back((int)pho->superCluster()->isInClean());
    nt.phoSCinUnClean.push_back((int)pho->superCluster()->isInUnclean());
    nt.phoSCnBC.push_back((int)pho->superCluster()->clustersSize());
    nt.phoESEn.push_back(pho->superCluster()->preshowerEnergy());
    
    nt.phoIsPFPhoton.push_back((int)pho->isPFlowPhoton());
    nt.phoIsStandardPhoton.push_back((int)pho->isStandardPhoton());
    nt.phoHasPixelSeed.push_back((int)pho->hasPixelSeed());
    nt.phoHasConversionTracks.push_back((int)pho->hasConversionTracks());
    // phoEleVeto.push_back((int)pho->passElectronVeto());   // TODO: not available in reco::
    nt.phoHadTowerOverEm.push_back(pho->hadTowOverEm());
    nt.phoHoverE.push_back(pho->hadronicOverEm());
    nt.phoHoverEValid.push_back(pho->hadronicOverEmValid());
    nt.phoSigmaIEtaIEta.push_back(pho->sigmaIetaIeta());
    nt.phoR9.push_back(pho->r9());
    nt.nPho++;
  }
}