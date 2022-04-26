#include "Analysis_CaloAnalyzer.h"
#include "ConfigManager.h"

using namespace std;
Analysis_CaloAnalyzer::Analysis_CaloAnalyzer(){};

Analysis_CaloAnalyzer::~Analysis_CaloAnalyzer(){};




void Analysis_CaloAnalyzer::FillCaloTowers(Analysis_NtupleContent& nt,const edm::Event& iEvent, const edm::EDGetTokenT<edm::SortedCollection<CaloTower>>& calotowerToken_) {
  iEvent.getByToken(calotowerToken_, calotowercollection);
   
   
 int nTowers = 0;
string configPath = "/afs/cern.ch/user/m/mnickel/private/LightByLight2018/analysis/configs/applySelections.md";
config = ConfigManager(configPath);

for (edm::SortedCollection<CaloTower>::const_iterator calo = calotowercollection->begin(); calo != calotowercollection->end(); ++calo) {
  float eta = calo->eta();
  float energy = calo->energy();

  string DetType = "";
  bool HasTowerAboveThreshold = false;
  
  if(eta > -maxEtaHF && eta < -minEtaHF) DetType = "HFm";
  if(eta >  minEtaHF && eta <  maxEtaHF) DetType = "HFp";
  // if(fabs(eta) > 0 && fabs(eta) < maxEtaHB) DetType = "HB";
  // if(fabs(eta) > minEtaHE && fabs(eta) < maxEtaHE) DetType = "HE";

  if(DetType == "") continue;
  
    if(DetType != ""){ // Check HF exclusivity
      if(energy > config.params("noiseThreshold"+DetType)){
        HasTowerAboveThreshold = true;
      }
    }
    
       if(HasTowerAboveThreshold){
       nTowers++;
       nt.CaloTower_emE.push_back(calo->emEnergy());
       nt.CaloTower_hadE.push_back(calo->hadEnergy());
       nt.CaloTower_e.push_back(calo->energy());
       nt.CaloTower_et.push_back(calo->et());
       nt.CaloTower_phi.push_back(calo->phi());
       nt.CaloTower_eta.push_back(calo->eta());
       }
}
  nt.nTower = nTowers; 
}


void Analysis_CaloAnalyzer::FillZDC(Analysis_NtupleContent& nt,const edm::Event& iEvent, const edm::EDGetTokenT<edm::SortedCollection<ZDCRecHit>>& RecHitsToken_) {
  iEvent.getByToken(RecHitsToken_, zdcrechits);

    nt.ZDC_n = zdcrechits->size();   
    int nhits = 0;
    for (auto const& rh : *zdcrechits) {
      HcalZDCDetId zdcid = rh.id();
      if (nhits  < 18) {
       nt.ZDC_e[nhits] = rh.energy();
       nt.ZDC_zside[nhits] = zdcid.zside();
       nt.ZDC_section[nhits] = zdcid.section();
       nt.ZDC_channel[nhits] = zdcid.channel();
       nt.ZDC_PM_Total_Energy+= rh.energy()  ;
       if(zdcid.zside() >0 ){
       nt.ZDC_P_Total_Energy += rh.energy();
       if(zdcid.section() == 1) nt.ZDC_P_ECal_Energy += rh.energy();
       if(zdcid.section() == 2) nt.ZDC_P_HCal_Energy += rh.energy()  ;        
       }
       if(zdcid.zside() <0 ){
       nt.ZDC_M_Total_Energy += rh.energy();
       if(zdcid.section() == 1) nt.ZDC_M_ECal_Energy += rh.energy();
       if(zdcid.section() == 2) nt.ZDC_M_HCal_Energy += rh.energy() ;         
       }
       // StandAlone_nt.ZDC_saturation[nhits] = static_cast<int>( rh.flagField(HcalCaloFlagLabels::ADCSaturationBit) );
      }

      nhits++;
    } // end loop zdc rechits 

}


