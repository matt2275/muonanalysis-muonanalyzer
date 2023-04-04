#include "Analysis_MuonGenAnalyzer.h"
#include <vector>

Analysis_MuonGenAnalyzer::Analysis_MuonGenAnalyzer(){};

Analysis_MuonGenAnalyzer::~Analysis_MuonGenAnalyzer(){};

void Analysis_MuonGenAnalyzer::SetInputsandFillNtuple_AllGenParticles(Analysis_NtupleContent& nt,const edm::Event& iEvent,
                                const edm::EDGetTokenT<edm::View<reco::GenParticle>>& gens_){
                                   
    nt.nGen = gens->size();
    for (unsigned int i = 0; i < gens->size() ; i++) {
          auto gen = gens->at(i);
          nt.gen_pdgId.push_back(gen.pdgId());
          nt.gen_pT.push_back(gen.pt());
          nt.gen_eta.push_back(gen.eta());
          nt.gen_phi.push_back(gen.phi()); 
          nt.gen_charge.push_back(gen.charge());
          nt.gen_mass.push_back(gen.mass());  
          nt.gen_vtx_x.push_back(gen.vx());
          nt.gen_vtx_y.push_back(gen.vy());
          nt.gen_vtx_z.push_back(gen.vz()); 
          nt.gen_status.push_back(gen.status());
          nt.gen_E.push_back(gen.energy());
          nt.gen_Et.push_back(gen.et());    
    }       
}


void Analysis_MuonGenAnalyzer::SetInputsandFillNtuple_MuMu(Analysis_NtupleContent& nt,const edm::Event& iEvent,
                                const edm::EDGetTokenT<edm::View<reco::GenParticle>>& gens_) {
  iEvent.getByToken(gens_, gens);
  bool found1 = false;
  bool found2 = false;
  index1 = -99;
  index2 = -99;
  gmuon1.SetPtEtaPhiM(0, 0, 0, 0);
  gmuon2.SetPtEtaPhiM(0, 0, 0, 0);
  nt.nGen = gens->size();
  for (unsigned int i = 0; i < gens->size() ; i++) {

    auto gen = gens->at(i);
    
    if(gen.pdgId() == 93){
       nt.indep_pt = gen.pt();
       nt.indep_eta = gen.eta();
       nt.indep_phi = gen.phi();
       nt.indep_mass = gen.mass();       
    }
   
    if (fabs(gen.pdgId()) != 13)
      continue;

    // look for final state muons
    if (gen.status() != 1)
      continue;

    if (gen.charge() < 0.) {
      if (found1) {
        std::cout << "Warning|Analysis_MuonGenAnalyzer::SetInputs| there are more than 2 final state mu- from the same "
                     "resonance -> skip"
                  << std::endl;
      } else {
        gmuon1.SetPtEtaPhiM(gen.pt(), gen.eta(), gen.phi(), MU_MASS);
        found1 = true;
        index1 = i;
      }
    } else {
      if (found2) {
        std::cout << "Warning|Analysis_MuonGenAnalyzer::SetInputs| there are more than 2 final state mu+ from the same "
                     "resonance -> skip"
                  << std::endl;
      } else {
        gmuon2.SetPtEtaPhiM(gen.pt(), gen.eta(), gen.phi(), MU_MASS);
        found2 = true;
        index2 = i;
      }
    }
  }
  
  if (index1 >= 0 && index2 >= 0 ) {
    auto gen1 = gens->at(index1);
    auto gen2 = gens->at(index2);
    gmuon1.SetPtEtaPhiM(gen1.pt(), gen1.eta(), gen1.phi(), gen1.mass());
    gmuon2.SetPtEtaPhiM(gen2.pt(), gen2.eta(), gen2.phi(), gen2.mass());
    nt.gentrk1_pt = gmuon1.Pt();
    nt.gentrk1_eta = gmuon1.Eta();
    nt.gentrk1_phi = gmuon1.Phi();
    nt.gentrk1_vtx_x = gen1.vx();
    nt.gentrk1_vtx_y = gen1.vy();
    nt.gentrk1_vtx_z = gen1.vz();
    nt.gentrk1_px = gen1.px();
    nt.gentrk1_py = gen1.py();
    nt.gentrk1_pz = gen1.pz();
    nt.gentrk1_charge = gen1.charge();
    nt.gentrk1_pdgId = gen1.pdgId();
    nt.gentrk2_pt = gmuon2.Pt();
    nt.gentrk2_eta = gmuon2.Eta();
    nt.gentrk2_phi = gmuon2.Phi();
    nt.gentrk2_vtx_x = gen2.vx();
    nt.gentrk2_vtx_y = gen2.vy();
    nt.gentrk2_vtx_z = gen2.vz();
    nt.gentrk2_px = gen2.px();
    nt.gentrk2_py = gen2.py();
    nt.gentrk2_pz = gen2.pz();
    nt.gentrk2_charge = gen2.charge();
    nt.gentrk2_pdgId = gen2.pdgId();
    nt.genFinalMass = (gmuon1 + gmuon2).M();
    nt.genFinalPt = (gmuon1 + gmuon2).Pt();
    nt.genFinalEta = (gmuon1 + gmuon2).Eta();
    nt.genFinalPhi = (gmuon1 + gmuon2).Phi();
    nt.genFinalE = (gmuon1 + gmuon2).E();
    nt.genFinalPx = (gmuon1 + gmuon2).Px();
    nt.genFinalPy= (gmuon1 + gmuon2).Py();
    nt.genFinalPz = (gmuon1 + gmuon2).Pz(); 
  } else {
    nt.gentrk1_pt = -99;
    nt.gentrk1_eta = -99;
    nt.gentrk1_phi = -99;
    nt.gentrk1_charge = -99;
    nt.gentrk1_vtx_x = -99;
    nt.gentrk1_vtx_y = -99;
    nt.gentrk1_vtx_z = -99;
    nt.gentrk1_px = -99;
    nt.gentrk1_py = -99;
    nt.gentrk1_pz = -99;
    nt.gentrk1_pdgId = -99;
    nt.gentrk2_pt = -99;
    nt.gentrk2_eta = -99;
    nt.gentrk2_phi = -99;
    nt.gentrk2_charge = -99;
    nt.gentrk2_vtx_x = -99;
    nt.gentrk2_vtx_y = -99;
    nt.gentrk2_vtx_z = -99;
    nt.gentrk2_px = -99;
    nt.gentrk2_py = -99;
    nt.gentrk2_pz = -99;
    nt.gentrk2_pdgId = -99;
    nt.genFinalMass = -99;
    nt.genFinalPt = -99;
    nt.genFinalEta = -99;
    nt.genFinalPhi = -99;
    nt.genFinalE = -99;
    nt.genFinalPx = -99;
    nt.genFinalPy= -99;
    nt.genFinalPz = -99;
  } 
  
}

void Analysis_MuonGenAnalyzer::SetInputsandFillNtuple_EE(Analysis_NtupleContent& nt,const edm::Event& iEvent,
                                const edm::EDGetTokenT<edm::View<reco::GenParticle>>& gens_) {
  iEvent.getByToken(gens_, gens);
  bool found1 = false;
  bool found2 = false;
  index1 = -99;
  index2 = -99;
  gmuon1.SetPtEtaPhiM(0, 0, 0, 0);
  gmuon2.SetPtEtaPhiM(0, 0, 0, 0);
  nt.nGen = gens->size();
  for (unsigned int i = 0; i < gens->size() ; i++) {

    auto gen = gens->at(i);
    
    if(gen.pdgId() == 93){
       nt.indep_pt = gen.pt();
       nt.indep_eta = gen.eta();
       nt.indep_phi = gen.phi();
       nt.indep_mass = gen.mass();       
    }
     
    if (fabs(gen.pdgId()) != 11)
      continue;

    // look for final state muons
    if (gen.status() != 1)
      continue;

    if (gen.charge() < 0.) {
      if (found1) {
        std::cout << "Warning|Analysis_MuonGenAnalyzer::SetInputs| there are more than 2 final state mu- from the same "
                     "resonance -> skip"
                  << std::endl;
      } else {
        gmuon1.SetPtEtaPhiM(gen.pt(), gen.eta(), gen.phi(), MU_MASS);
        found1 = true;
        index1 = i;
      }
    } else {
      if (found2) {
        std::cout << "Warning|Analysis_MuonGenAnalyzer::SetInputs| there are more than 2 final state mu+ from the same "
                     "resonance -> skip"
                  << std::endl;
      } else {
        gmuon2.SetPtEtaPhiM(gen.pt(), gen.eta(), gen.phi(), MU_MASS);
        found2 = true;
        index2 = i;
      }
    }
  }
  
  if (index1 >= 0 && index2 >= 0 ) {
    auto gen1 = gens->at(index1);
    auto gen2 = gens->at(index2);
    gmuon1.SetPtEtaPhiM(gen1.pt(), gen1.eta(), gen1.phi(), gen1.mass());
    gmuon2.SetPtEtaPhiM(gen2.pt(), gen2.eta(), gen2.phi(), gen2.mass());
    nt.gentrk1_pt = gmuon1.Pt();
    nt.gentrk1_eta = gmuon1.Eta();
    nt.gentrk1_phi = gmuon1.Phi();
    nt.gentrk1_vtx_x = gen1.vx();
    nt.gentrk1_vtx_y = gen1.vy();
    nt.gentrk1_vtx_z = gen1.vz();
    nt.gentrk1_px = gen1.px();
    nt.gentrk1_py = gen1.py();
    nt.gentrk1_pz = gen1.pz();
    nt.gentrk1_charge = gen1.charge();
    nt.gentrk1_pdgId = gen1.pdgId();
    nt.gentrk2_pt = gmuon2.Pt();
    nt.gentrk2_eta = gmuon2.Eta();
    nt.gentrk2_phi = gmuon2.Phi();
    nt.gentrk2_vtx_x = gen2.vx();
    nt.gentrk2_vtx_y = gen2.vy();
    nt.gentrk2_vtx_z = gen2.vz();
    nt.gentrk2_px = gen2.px();
    nt.gentrk2_py = gen2.py();
    nt.gentrk2_pz = gen2.pz();
    nt.gentrk2_charge = gen2.charge();
    nt.gentrk2_pdgId = gen2.pdgId();
    nt.genFinalMass = (gmuon1 + gmuon2).M();
    nt.genFinalPt = (gmuon1 + gmuon2).Pt();
    nt.genFinalEta = (gmuon1 + gmuon2).Eta();
    nt.genFinalPhi = (gmuon1 + gmuon2).Phi();
    nt.genFinalE = (gmuon1 + gmuon2).E();
    nt.genFinalPx = (gmuon1 + gmuon2).Px();
    nt.genFinalPy= (gmuon1 + gmuon2).Py();
    nt.genFinalPz = (gmuon1 + gmuon2).Pz();    
  } else {
    nt.gentrk1_pt = -99;
    nt.gentrk1_eta = -99;
    nt.gentrk1_phi = -99;
    nt.gentrk1_charge = -99;
    nt.gentrk1_vtx_x = -99;
    nt.gentrk1_vtx_y = -99;
    nt.gentrk1_vtx_z = -99;
    nt.gentrk1_px = -99;
    nt.gentrk1_py = -99;
    nt.gentrk1_pz = -99;
    nt.gentrk1_pdgId = -99;
    nt.gentrk2_pt = -99;
    nt.gentrk2_eta = -99;
    nt.gentrk2_phi = -99;
    nt.gentrk2_charge = -99;
    nt.gentrk2_vtx_x = -99;
    nt.gentrk2_vtx_y = -99;
    nt.gentrk2_vtx_z = -99;
    nt.gentrk2_px = -99;
    nt.gentrk2_py = -99;
    nt.gentrk2_pz = -99;
    nt.gentrk2_pdgId = -99;
    nt.genFinalMass = -99;
    nt.genFinalPt = -99;
    nt.genFinalEta = -99;
    nt.genFinalPhi = -99;
    nt.genFinalE = -99;
    nt.genFinalPx = -99;
    nt.genFinalPy= -99;
    nt.genFinalPz = -99; 
  } 
  
}




void Analysis_MuonGenAnalyzer::SetInputsandFillNtuple_TauTau(Analysis_NtupleContent& nt,const edm::Event& iEvent,
                                const edm::EDGetTokenT<edm::View<reco::GenParticle>>& gens_) {
  iEvent.getByToken(gens_, gens);
  // bool found1 = false;
  // bool found2 = false;
  index1 = -99;
  index2 = -99;
  float index1_pt = 0;
  float index2_pt = 0;
  int other_index1 = -99;
  int other_index2 = -99;
  float other_index1_pt = 0;
  float other_index2_pt = 0;
  int a[] = {11, 13, 211, 111};
  int final_parts[] = {11,13,211,111,12,14,16};
  int tau_index1 = -99;
  int tau_index2 = -99;
  std::vector<int> tau_3prong_index1;
  std::vector<int> tau_3prong_index2;
  int tau_plus_mu = 0;
  int tau_plus_e = 0;
  int tau_plus_pi_plus = 0;
  int tau_plus_pi_minus = 0;
  int tau_plus_pi0 = 0;
  int tau_minus_mu = 0;
  int tau_minus_e = 0;
  int tau_minus_pi_plus = 0;
  int tau_minus_pi_minus = 0;  
  int tau_minus_pi0 = 0;
  gmuon1.SetPtEtaPhiM(0, 0, 0, 0);
  gmuon2.SetPtEtaPhiM(0, 0, 0, 0);
  nt.nGen = gens->size();
  // std::cout << " new entry " << std::endl;
  for (unsigned int i = 0; i < gens->size() ; i++) {
    auto gen = gens->at(i);
    int genID= gen.pdgId();
    int status  = gen.status();
    if(genID == 93){
       nt.indep_pt = gen.pt();
       nt.indep_eta = gen.eta();
       nt.indep_phi = gen.phi();
       nt.indep_mass = gen.mass();       
    }
    if(genID == 22 && status == 21){
       float gamma_pz = gen.pz();
       if(gamma_pz > 0) nt.gentau1_gamma_pz = gamma_pz;
       else nt.gentau2_gamma_pz = gamma_pz;
    }
       // std::cout << " pdgID " << genID << " status " << gen.status() << " charge " << gen.charge() << " pt " << gen.pt() << " eta " 
       // << gen.eta()<< "phi " << gen.phi() << " vx " <<  gen.vx() << " vy " <<  gen.vy() << " vz " <<  gen.vz() <<std::endl;

         if(abs(genID) ==15){
       const reco::Candidate* mother_tmp = gen.mother();
       if (abs(mother_tmp->pdgId()) != 15){
       
       // std::cout << " status " << gen.status() << "charge " << gen.charge() << " pt " << gen.pt() << " eta " << gen.eta()<< "phi " << gen.phi() << " vx " <<  gen.vx() << " vy " <<  gen.vy() << " vz " <<  gen.vz() <<std::endl;
       if(genID==15) tau_index1 = i;
       if(genID==-15) tau_index2 = i;
       }
    }
    // if(genID==15 && gen.status()==2) tau_index1 = i;
    // if(genID==-15 && gen.status()==2) tau_index2 = i;
    if(gen.status() != 1 ){
    if(gen.pt() > other_index1_pt && gen.charge() < 0){
       other_index1 = i;
       other_index1_pt = gen.pt();
    }
    if(gen.pt() > other_index2_pt && gen.charge() > 0){
       other_index2 = i;
       other_index2_pt = gen.pt();
    }
    }
    if(status != 1 and status != 2) continue;
        bool has_final_part = std::find(std::begin(final_parts), std::end(final_parts), fabs(genID)) != std::end(final_parts);
    if(has_final_part){
          const reco::Candidate* mother_tmp = gen.mother();
      int max_depth =0 ;
      while(mother_tmp){
      if(mother_tmp->pdgId()== 15 && mother_tmp->status() == 2){
             nt.nGen_tau1++;
             nt.gen_pdgId_tau1.push_back(gen.pdgId());
             nt.gen_pT_tau1.push_back(gen.pt());
             nt.gen_eta_tau1.push_back(gen.eta());
             nt.gen_phi_tau1.push_back(gen.phi()); 
             nt.gen_charge_tau1.push_back(gen.charge());
             nt.gen_mass_tau1.push_back(gen.mass());  
             nt.gen_vtx_x_tau1.push_back(gen.vx());
             nt.gen_vtx_y_tau1.push_back(gen.vy());
             nt.gen_vtx_z_tau1.push_back(gen.vz()); 
             nt.gen_status_tau1.push_back(gen.status());
             nt.gen_E_tau1.push_back(gen.energy());
             nt.gen_Et_tau1.push_back(gen.et()); 
             break;
                    
      }
      
      
      if(mother_tmp->pdgId()== -15 && mother_tmp->status() == 2){
             nt.nGen_tau2++;
             nt.gen_pdgId_tau2.push_back(gen.pdgId());
             nt.gen_pT_tau2.push_back(gen.pt());
             nt.gen_eta_tau2.push_back(gen.eta());
             nt.gen_phi_tau2.push_back(gen.phi()); 
             nt.gen_charge_tau2.push_back(gen.charge());
             nt.gen_mass_tau2.push_back(gen.mass());  
             nt.gen_vtx_x_tau2.push_back(gen.vx());
             nt.gen_vtx_y_tau2.push_back(gen.vy());
             nt.gen_vtx_z_tau2.push_back(gen.vz()); 
             nt.gen_status_tau2.push_back(gen.status());
             nt.gen_E_tau2.push_back(gen.energy());
             nt.gen_Et_tau2.push_back(gen.et()); 
             break;
                   
      }
        if(max_depth > 2) break;
        max_depth ++;
        mother_tmp = mother_tmp->mother(); 
     }

    }
    bool exists = std::find(std::begin(a), std::end(a), fabs(genID)) != std::end(a);
    if(!exists) continue;   
    const reco::Candidate* mother_tmp = gen.mother();
    while (mother_tmp) {
       if (fabs(mother_tmp->pdgId()) == 15){      
          if(mother_tmp->pdgId() == 15){
             if(genID == 13) tau_minus_mu++;
             if(genID == 11) tau_minus_e++;
             if(genID == 211) tau_minus_pi_plus++;
             if(genID == -211) tau_minus_pi_minus++;             
             if(genID == 111) tau_minus_pi0++;
             if(fabs(genID)==211) tau_3prong_index1.push_back(i);
             if(gen.pt() > index1_pt  && gen.charge() < 0){
             index1=i;
             index1_pt = gen.pt();
             
             
             }
          }
          if(mother_tmp->pdgId() == -15){
             if(genID == -13) tau_plus_mu++;
             if(genID == -11) tau_plus_e++;
             if(genID == 211) tau_plus_pi_plus++;
             if(genID == -211) tau_plus_pi_minus++;             
             if(genID == 111) tau_plus_pi0++;
             if(fabs(genID)==211) tau_3prong_index2.push_back(i);
             if(gen.pt() > index2_pt  && gen.charge() > 0){
             index2=i;
             index2_pt = gen.pt();
             }
          }
            
          break;      
       }
      mother_tmp = mother_tmp->mother();       
   }
 


    
  }

  if(tau_minus_e == 1 && tau_minus_mu == 0 && tau_minus_pi_minus == 0 && tau_minus_pi_plus == 0){
     nt.gentau1_decay_el =true;  
     nt.gentau1_decay = 0;
  }
  else if(tau_minus_e == 0 && tau_minus_mu == 1 && tau_minus_pi_minus == 0 && tau_minus_pi_plus == 0){
     nt.gentau1_decay_mu =true; 
     nt.gentau1_decay = 1;
  }
  else if(tau_minus_e == 0 && tau_minus_mu == 0 && tau_minus_pi_minus == 1 && tau_minus_pi_plus == 0){
     nt.gentau1_decay_1prong =true;
      nt.gentau1_decay = 2;
  }
  else if(tau_minus_e == 0 && tau_minus_mu == 0 && tau_minus_pi_minus == 2 && tau_minus_pi_plus == 1){
     nt.gentau1_decay_3prong =true;
     nt.gentau1_decay = 3;
  }
  else{
     nt.gentau1_decay_other = true;
     nt.gentau1_decay = 4;
     index1 = other_index1;
  }
  
  
  if(tau_plus_e == 1 && tau_plus_mu == 0 && tau_plus_pi_minus == 0 && tau_plus_pi_plus == 0){
     nt.gentau2_decay_el =true;  
     nt.gentau2_decay = 0;
  }
  else if(tau_plus_e == 0 && tau_plus_mu == 1 && tau_plus_pi_minus == 0 && tau_plus_pi_plus == 0){
     nt.gentau2_decay_mu =true; 
     nt.gentau2_decay = 1;
  }
  else if(tau_plus_e == 0 && tau_plus_mu == 0 && tau_plus_pi_minus == 0 && tau_plus_pi_plus == 1){
     nt.gentau2_decay_1prong =true;
      nt.gentau2_decay = 2;
  }
  else if(tau_plus_e == 0 && tau_plus_mu == 0 && tau_plus_pi_minus == 1 && tau_plus_pi_plus == 2){
     nt.gentau2_decay_3prong =true;
     nt.gentau2_decay = 3;
  }
  else{
     nt.gentau2_decay_other = true;
     nt.gentau2_decay = 4;
     index2 = other_index2;     
  }
  
  nt.gentau1_N_pi0 = tau_minus_pi0;
  nt.gentau2_N_pi0 = tau_plus_pi0;
  
  if (index1 >= 0 && index2 >= 0 ) {
    auto gen1 = gens->at(index1);
    auto gen2 = gens->at(index2);
    gmuon1.SetPtEtaPhiM(0, 0, 0, 0);
    gmuon2.SetPtEtaPhiM(0, 0, 0, 0);
    if(nt.gentau1_decay_3prong){
       for(int gen_index : tau_3prong_index1){
          TLorentzVector tmp_vec;
          auto tmp_gen = gens->at(gen_index);
          tmp_vec.SetPtEtaPhiM(tmp_gen.pt(), tmp_gen.eta(), tmp_gen.phi(), tmp_gen.mass());
          gmuon1 = gmuon1 + tmp_vec;
       }
    }
    if(!nt.gentau1_decay_3prong)    gmuon1.SetPtEtaPhiM(gen1.pt(), gen1.eta(), gen1.phi(), gen1.mass());
    if(nt.gentau2_decay_3prong){
       for(int gen_index : tau_3prong_index2){
          TLorentzVector tmp_vec;
          auto tmp_gen = gens->at(gen_index);
          tmp_vec.SetPtEtaPhiM(tmp_gen.pt(), tmp_gen.eta(), tmp_gen.phi(), tmp_gen.mass());
          gmuon2 = gmuon2 + tmp_vec;
       }
    }
    if(!nt.gentau2_decay_3prong)    gmuon2.SetPtEtaPhiM(gen2.pt(), gen2.eta(), gen2.phi(), gen2.mass());
    nt.gentrk1_pt = gmuon1.Pt();
    nt.gentrk1_eta = gmuon1.Eta();
    nt.gentrk1_phi = gmuon1.Phi();
    nt.gentrk1_vtx_x = gen1.vx();
    nt.gentrk1_vtx_y = gen1.vy();
    nt.gentrk1_vtx_z = gen1.vz();
    nt.gentrk1_px = gmuon1.Px();
    nt.gentrk1_py = gmuon1.Py();
    nt.gentrk1_pz = gmuon1.Pz();
    nt.gentrk1_charge = gen1.charge();
    nt.gentrk1_M = gmuon1.M();
    nt.gentrk1_pdgId = gen1.pdgId();
    nt.gentrk2_pt = gmuon2.Pt();
    nt.gentrk2_eta = gmuon2.Eta();
    nt.gentrk2_phi = gmuon2.Phi();
    nt.gentrk2_vtx_x = gen2.vx();
    nt.gentrk2_vtx_y = gen2.vy();
    nt.gentrk2_vtx_z = gen2.vz();
    nt.gentrk2_px = gmuon2.Px();
    nt.gentrk2_py = gmuon2.Py();
    nt.gentrk2_pz = gmuon2.Pz();
    nt.gentrk2_charge = gen2.charge();
    nt.gentrk2_M = gmuon2.M();
    nt.gentrk2_pdgId = gen2.pdgId();
    nt.genFinalMass = (gmuon1 + gmuon2).M();
    nt.genFinalPt = (gmuon1 + gmuon2).Pt();
    nt.genFinalEta = (gmuon1 + gmuon2).Eta();
    nt.genFinalPhi = (gmuon1 + gmuon2).Phi();
    nt.genFinalE = (gmuon1 + gmuon2).E();
    nt.genFinalPx = (gmuon1 + gmuon2).Px();
    nt.genFinalPy= (gmuon1 + gmuon2).Py();
    nt.genFinalPz = (gmuon1 + gmuon2).Pz(); 
  } else {
    nt.gentrk1_pt = -99;
    nt.gentrk1_eta = -99;
    nt.gentrk1_phi = -99;
    nt.gentrk1_charge = -99;
    nt.gentrk1_vtx_x = -99;
    nt.gentrk1_vtx_y = -99;
    nt.gentrk1_vtx_z = -99;
    nt.gentrk1_px = -99;
    nt.gentrk1_py = -99;
    nt.gentrk1_pz = -99;
    nt.gentrk1_pdgId = -99;
    nt.gentrk2_pt = -99;
    nt.gentrk2_eta = -99;
    nt.gentrk2_phi = -99;
    nt.gentrk2_charge = -99;
    nt.gentrk2_vtx_x = -99;
    nt.gentrk2_vtx_y = -99;
    nt.gentrk2_vtx_z = -99;
    nt.gentrk2_px = -99;
    nt.gentrk2_py = -99;
    nt.gentrk2_pz = -99;
    nt.gentrk2_pdgId = -99;
    nt.genFinalMass = -99;
    nt.genFinalPt = -99;
    nt.genFinalEta = -99;
    nt.genFinalPhi = -99;
    nt.genFinalE = -99;
    nt.genFinalPx = -99;
    nt.genFinalPy= -99;
    nt.genFinalPz = -99; 
  } 
  
  if (tau_index1 >= 0 && tau_index2 >= 0 ) {
    auto gen1 = gens->at(tau_index1);
    auto gen2 = gens->at(tau_index2);
    gmuon1.SetPtEtaPhiM(gen1.pt(), gen1.eta(), gen1.phi(), gen1.mass());
    gmuon2.SetPtEtaPhiM(gen2.pt(), gen2.eta(), gen2.phi(), gen2.mass());
    nt.gentau1_pt = gmuon1.Pt();
    nt.gentau1_eta = gmuon1.Eta();
    nt.gentau1_phi = gmuon1.Phi();
    nt.gentau1_vtx_x = gen1.vx();
    nt.gentau1_vtx_y = gen1.vy();
    nt.gentau1_vtx_z = gen1.vz();
    nt.gentau1_px = gen1.px();
    nt.gentau1_py = gen1.py();
    nt.gentau1_pz = gen1.pz();
    nt.gentau1_charge = gen1.charge();
    nt.gentau1_pdgId = gen1.pdgId();
    nt.gentau2_pt = gmuon2.Pt();
    nt.gentau2_eta = gmuon2.Eta();
    nt.gentau2_phi = gmuon2.Phi();
    nt.gentau2_vtx_x = gen2.vx();
    nt.gentau2_vtx_y = gen2.vy();
    nt.gentau2_vtx_z = gen2.vz();
    nt.gentau2_px = gen2.px();
    nt.gentau2_py = gen2.py();
    nt.gentau2_pz = gen2.pz();
    nt.gentau2_charge = gen2.charge();
    nt.gentau2_pdgId = gen2.pdgId();
    nt.genDiTauMass = (gmuon1 + gmuon2).M();
    nt.genDiTauPt = (gmuon1 + gmuon2).Pt();
    nt.genDiTauEta = (gmuon1 + gmuon2).Eta();
    nt.genDiTauPhi = (gmuon1 + gmuon2).Phi();
    nt.genDiTauE = (gmuon1 + gmuon2).E();
    nt.genDiTauPx = (gmuon1 + gmuon2).Px();
    nt.genDiTauPy= (gmuon1 + gmuon2).Py();
    nt.genDiTauPz = (gmuon1 + gmuon2).Pz(); 
  } else {
    nt.gentau1_pt = -99;
    nt.gentau1_eta = -99;
    nt.gentau1_phi = -99;
    nt.gentau1_charge = -99;
    nt.gentau1_vtx_x = -99;
    nt.gentau1_vtx_y = -99;
    nt.gentau1_vtx_z = -99;
    nt.gentau1_px = -99;
    nt.gentau1_py = -99;
    nt.gentau1_pz = -99;
    nt.gentau1_pdgId = -99;
    nt.gentau2_pt = -99;
    nt.gentau2_eta = -99;
    nt.gentau2_phi = -99;
    nt.gentau2_charge = -99;
    nt.gentau2_vtx_x = -99;
    nt.gentau2_vtx_y = -99;
    nt.gentau2_vtx_z = -99;
    nt.gentau2_px = -99;
    nt.gentau2_py = -99;
    nt.gentau2_pz = -99;
    nt.gentau2_pdgId = -99;
    nt.genDiTauMass = -99;
    nt.genDiTauPt = -99;
    nt.genDiTauEta = -99;
    nt.genDiTauPhi = -99;
    nt.genDiTauE = -99;
    nt.genDiTauPx = -99;
    nt.genDiTauPy= -99;
    nt.genDiTauPz = -99; 
  } 
  
}

void Analysis_MuonGenAnalyzer::SetInputsandFillNtuple_Other(Analysis_NtupleContent& nt,const edm::Event& iEvent, const edm::EDGetTokenT<edm::View<reco::GenParticle>>& gens_, const float& pt_cut) {
  iEvent.getByToken(gens_, gens);
  // bool found1 = false;
  // bool found2 = false;
  index1 = -99;
  index2 = -99;
  float index1_pt = 0;
  float index2_pt = 0;
  gmuon1.SetPtEtaPhiM(0, 0, 0, 0);
  gmuon2.SetPtEtaPhiM(0, 0, 0, 0);
  nt.nGen = gens->size();
  for (unsigned int i = 0; i < gens->size() ; i++) {

    auto gen = gens->at(i);
    // if (fabs(gen.pdgId()) != 11)
      // continue;

    // look for final state muons
    if (gen.status() != 1)
      continue;
    if (gen.pt() < pt_cut) 
       continue;
    if(gen.charge() < 0 ){
       if(gen.pt() > index1_pt){
          index1  = i;
          index1_pt = gen.pt();
       }
    }
    if(gen.charge() > 0 ){
       if(gen.pt() > index2_pt){
          index2  = i;
          index2_pt = gen.pt();
       }
    }
  }
  
  if (index1 >= 0 && index2 >= 0 ) {
    auto gen1 = gens->at(index1);
    auto gen2 = gens->at(index2);
    gmuon1.SetPtEtaPhiM(gen1.pt(), gen1.eta(), gen1.phi(), gen1.mass());
    gmuon2.SetPtEtaPhiM(gen2.pt(), gen2.eta(), gen2.phi(), gen2.mass());
    nt.gentrk1_pt = gmuon1.Pt();
    nt.gentrk1_eta = gmuon1.Eta();
    nt.gentrk1_phi = gmuon1.Phi();
    nt.gentrk1_vtx_x = gen1.vx();
    nt.gentrk1_vtx_y = gen1.vy();
    nt.gentrk1_vtx_z = gen1.vz();
    nt.gentrk1_px = gen1.px();
    nt.gentrk1_py = gen1.py();
    nt.gentrk1_pz = gen1.pz();
    nt.gentrk1_charge = gen1.charge();
    nt.gentrk1_pdgId = gen1.pdgId();
    nt.gentrk2_pt = gmuon2.Pt();
    nt.gentrk2_eta = gmuon2.Eta();
    nt.gentrk2_phi = gmuon2.Phi();
    nt.gentrk2_vtx_x = gen2.vx();
    nt.gentrk2_vtx_y = gen2.vy();
    nt.gentrk2_vtx_z = gen2.vz();
    nt.gentrk2_px = gen2.px();
    nt.gentrk2_py = gen2.py();
    nt.gentrk2_pz = gen2.pz();
    nt.gentrk2_charge = gen2.charge();
    nt.gentrk2_pdgId = gen2.pdgId();
    nt.genFinalMass = (gmuon1 + gmuon2).M();
    nt.genFinalPt = (gmuon1 + gmuon2).Pt();
    nt.genFinalEta = (gmuon1 + gmuon2).Eta();
    nt.genFinalPhi = (gmuon1 + gmuon2).Phi();
    nt.genFinalE = (gmuon1 + gmuon2).E();
    nt.genFinalPx = (gmuon1 + gmuon2).Px();
    nt.genFinalPy= (gmuon1 + gmuon2).Py();
    nt.genFinalPz = (gmuon1 + gmuon2).Pz(); 
  } else {
    nt.gentrk1_pt = -99;
    nt.gentrk1_eta = -99;
    nt.gentrk1_phi = -99;
    nt.gentrk1_charge = -99;
    nt.gentrk1_vtx_x = -99;
    nt.gentrk1_vtx_y = -99;
    nt.gentrk1_vtx_z = -99;
    nt.gentrk1_px = -99;
    nt.gentrk1_py = -99;
    nt.gentrk1_pz = -99;
    nt.gentrk1_pdgId = -99;
    nt.gentrk2_pt = -99;
    nt.gentrk2_eta = -99;
    nt.gentrk2_phi = -99;
    nt.gentrk2_charge = -99;
    nt.gentrk2_vtx_x = -99;
    nt.gentrk2_vtx_y = -99;
    nt.gentrk2_vtx_z = -99;
    nt.gentrk2_px = -99;
    nt.gentrk2_py = -99;
    nt.gentrk2_pz = -99;
    nt.gentrk2_pdgId = -99;
    nt.genFinalMass = -99;
    nt.genFinalPt = -99;
    nt.genFinalEta = -99;
    nt.genFinalPhi = -99;
    nt.genFinalE = -99;
    nt.genFinalPx = -99;
    nt.genFinalPy= -99;
    nt.genFinalPz = -99; 
  } 
  
}