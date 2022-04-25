#include "Analysis_KlFitter.h"

Analysis_KlFitter::Analysis_KlFitter(std::vector<reco::TransientTrack> &vecttrk) {
  if(vecttrk.size()==0){
     status_ = false;
     return;
  }
  KalmanVertexFitter vtxFitter(true);
  dimuvtx = vtxFitter.vertex(vecttrk);
  prefit.push_back(vecttrk.at(0));
  prefit.push_back(vecttrk.at(1));
  if (!dimuvtx.isValid())
    status_ = false;
  else {
    refited = dimuvtx.refittedTracks();
    prob_ = ChiSquaredProbability(dimuvtx.totalChiSquared(), dimuvtx.degreesOfFreedom());
    normalchi2_ = dimuvtx.totalChiSquared() / dimuvtx.degreesOfFreedom();
  }
}



Analysis_KlFitter::~Analysis_KlFitter(){};

// void Analysis_KlFitter::genDiTau(std::vector<reco::TransientTrack> &vecttrk) {
  // std::cout << " Genditau " << std::endl;
  // KalmanVertexFitter vtxFitter(true);
  // genDiTau_dimuvtx = vtxFitter.vertex(vecttrk);
  // genDiTau_prefit.push_back(vecttrk.at(0));
  // genDiTau_prefit.push_back(vecttrk.at(1));
  // if (!genDiTau_dimuvtx.isValid())
    // genDiTau_status_ = false;
  // else {
    // genDiTau_status_ = true;
    // genDiTau_refited = dimuvtx.refittedTracks();
    // genDiTau_prob_ = ChiSquaredProbability(dimuvtx.totalChiSquared(), dimuvtx.degreesOfFreedom());
    // genDiTau_normalchi2_ = dimuvtx.totalChiSquared() / dimuvtx.degreesOfFreedom();
  // }
// }


void Analysis_KlFitter::genDiTau(Analysis_NtupleContent & nt){
  double distance;
  gentau1_x.SetXYZ(nt.gentau1_vtx_x,nt.gentau1_vtx_y,nt.gentau1_vtx_z);
  gentau1_p.SetXYZ(nt.gentau1_px,nt.gentau1_py,nt.gentau1_pz);
 
     
  gentau2_x.SetXYZ(nt.gentau2_vtx_x,nt.gentau2_vtx_y,nt.gentau2_vtx_z);
  gentau2_p.SetXYZ(nt.gentau2_px,nt.gentau2_py,nt.gentau2_pz);
   
   XYZVector diff_x = gentau1_x - gentau2_x; 
   XYZVector sum_p = gentau1_p + gentau2_p; 
   if(isZeroVector(sum_p,.000000001)){
         gentau1_refit_vtx = gentau1_x;
         gentau2_refit_vtx = gentau2_x;
         genDiTau_vtx = .5*(gentau1_refit_vtx + gentau2_refit_vtx);
         distance = sqrt((gentau1_refit_vtx - gentau2_refit_vtx).Mag2());
         genDiTau_status_ = true;
   }else{                           
   double diff_xdp1 = diff_x.Dot(gentau1_p);                           
   double diff_xdp2 = diff_x.Dot(gentau2_p);
   double p1dp1 = gentau1_p.Dot(gentau1_p);
   double p1dp2 = gentau1_p.Dot(gentau2_p);
   double p2dp2 = gentau2_p.Dot(gentau2_p);

   double det = p1dp1*p2dp2 - p1dp2*p1dp2;
   if(det==0) genDiTau_status_ = false;
   double s = (-diff_xdp1*p2dp2 + diff_xdp2*p1dp2)/det; 
   double t = (diff_xdp1*p1dp2 - diff_xdp2*p1dp1)/det;
   
   gentau1_refit_vtx = gentau1_x + gentau1_p*s;
   gentau2_refit_vtx = gentau2_x + gentau2_p*t;
   genDiTau_vtx = .5*(gentau1_refit_vtx + gentau2_refit_vtx);
   distance = sqrt((gentau1_refit_vtx - gentau2_refit_vtx).Mag2());
   }
  
      if (genDiTau_status_ ) {
      nt.genDiTau_vtx_x = genDiTau_vtx.X();
      nt.genDiTau_vtx_y =genDiTau_vtx.Y();
      nt.genDiTau_vtx_z = genDiTau_vtx.Z();
      nt.genDiTau_vtx_gentau1_IP = distance/2;
      nt.genDiTau_vtx_gentau2_IP = distance/2; 
    } else {
      nt.genDiTau_vtx_x = -99;
      nt.genDiTau_vtx_y = -99;
      nt.genDiTau_vtx_z = -99;
      nt.genDiTau_vtx_gentau1_IP = -99;
      nt.genDiTau_vtx_gentau1_IP = -99; 
    }
                             
}

void Analysis_KlFitter::genFinal(Analysis_NtupleContent & nt){
  double distance;
  gentrk1_x.SetXYZ(nt.gentrk1_vtx_x,nt.gentrk1_vtx_y,nt.gentrk1_vtx_z);
  gentrk1_p.SetXYZ(nt.gentrk1_px,nt.gentrk1_py,nt.gentrk1_pz);
 
     
  gentrk2_x.SetXYZ(nt.gentrk2_vtx_x,nt.gentrk2_vtx_y,nt.gentrk2_vtx_z);
  gentrk2_p.SetXYZ(nt.gentrk2_px,nt.gentrk2_py,nt.gentrk2_pz);
   
   XYZVector diff_x = gentrk1_x - gentrk2_x; 
   XYZVector sum_p = gentrk1_p + gentrk2_p; 
   if(isZeroVector(sum_p,.000000001)){
         gentrk1_refit_vtx = gentrk1_x;
         gentrk2_refit_vtx = gentrk2_x;
         genFinal_vtx = .5*(gentrk1_refit_vtx + gentrk2_refit_vtx);
         distance = sqrt((gentrk1_refit_vtx - gentrk2_refit_vtx).Mag2());
         genFinal_status_ = true;
   }else{   
   double diff_xdp1 = diff_x.Dot(gentrk1_p);                           
   double diff_xdp2 = diff_x.Dot(gentrk2_p);
   double p1dp1 = gentrk1_p.Dot(gentrk1_p);
   double p1dp2 = gentrk1_p.Dot(gentrk2_p);
   double p2dp2 = gentrk2_p.Dot(gentrk2_p);

   double det = p1dp1*p2dp2 - p1dp2*p1dp2;
   if(det==0) genFinal_status_ = false;
   double s = (-diff_xdp1*p2dp2 + diff_xdp2*p1dp2)/det; 
   double t = (diff_xdp1*p1dp2 - diff_xdp2*p1dp1)/det;
   
   gentrk1_refit_vtx = gentrk1_x + gentrk1_p*s;
   gentrk2_refit_vtx = gentrk2_x + gentrk2_p*t;
   genFinal_vtx = .5*(gentrk1_refit_vtx + gentrk2_refit_vtx);
   distance = sqrt((gentrk1_refit_vtx - gentrk2_refit_vtx).Mag2());
   }
  
      if (genFinal_status_ ) {
      nt.genFinal_vtx_x = genFinal_vtx.X();
      nt.genFinal_vtx_y =genFinal_vtx.Y();
      nt.genFinal_vtx_z = genFinal_vtx.Z();
      nt.genFinal_vtx_gentrk1_IP = distance/2;
      nt.genFinal_vtx_gentrk2_IP = distance/2;    
    } else {
      nt.genFinal_vtx_x = -99;
      nt.genFinal_vtx_y = -99;
      nt.genFinal_vtx_z = -99;
      nt.genFinal_vtx_gentrk1_IP = -99;
      nt.genFinal_vtx_gentrk1_IP = -99; 
    }
                             
}

void Analysis_KlFitter::genDiTau_genFinal_vertex(Analysis_NtupleContent &nt) {
   if(genDiTau_status_ and genFinal_status_){
      XYZVector V1 = gentau1_x -genFinal_vtx;
      nt.genFinal_vtx_gentau1_IP = distancePointLine(V1,gentau1_p);
      V1 = gentau2_x -genFinal_vtx;
      nt.genFinal_vtx_gentau2_IP = distancePointLine(V1 ,gentau2_p);
      V1 = gentrk1_x -genDiTau_vtx;
      nt.genDiTau_vtx_gentrk1_IP = distancePointLine(V1,gentrk1_p);
      V1 = gentrk2_x -genDiTau_vtx;
      nt.genDiTau_vtx_gentrk2_IP = distancePointLine(V1,gentrk2_p);        
      
      
   }else{
      nt.genFinal_vtx_gentau1_IP = -99;
      nt.genFinal_vtx_gentau2_IP = -99;
      nt.genDiTau_vtx_gentrk1_IP = -99;
      nt.genDiTau_vtx_gentrk2_IP = -99;          
   }
}

// void Analysis_KlFitter::genFinal(std::vector<reco::TransientTrack> &vecttrk) {
  // std::cout << " GenFinal " << std::endl;
  // KalmanVertexFitter vtxFitter(true);
  // genFinal_dimuvtx = vtxFitter.vertex(vecttrk);
  // genFinal_prefit.push_back(vecttrk.at(0));
  // genFinal_prefit.push_back(vecttrk.at(1));
  // if (!  genFinal_dimuvtx.isValid())
      // genFinal_status_ = false;
  // else {
    // genFinal_status_ = true;
    // genFinal_refited = dimuvtx.refittedTracks();
    // genFinal_prob_ = ChiSquaredProbability(dimuvtx.totalChiSquared(), dimuvtx.degreesOfFreedom());
    // genFinal_normalchi2_ = dimuvtx.totalChiSquared() / dimuvtx.degreesOfFreedom();
  // }
// }

void Analysis_KlFitter::fillNtuple(Analysis_NtupleContent &nt, bool isTuneP) {
  if (isTuneP) {
    if (status_) {
      nt.pair_tuneP_svprob = prob_;
      nt.pair_tuneP_fit_mass = DimuonMass(refited[0].track().pt(),
                                          refited[0].track().eta(),
                                          refited[0].track().phi(),
                                          refited[1].track().pt(),
                                          refited[1].track().eta(),
                                          refited[1].track().phi());
      nt.pair_tuneP_normalchi2 = normalchi2_;
      nt.refit_tuneP_vtx_x = dimuvtx.position().x();
      nt.refit_tuneP_vtx_y =dimuvtx.position().y();
      nt.refit_tuneP_vtx_z = dimuvtx.position().z();
      GlobalPoint vert(nt.refit_vtx_x, nt.refit_vtx_y, nt.refit_vtx_z);
      TrajectoryStateClosestToPoint  traj0 = refited[0].trajectoryStateClosestToPoint(vert);
      TrajectoryStateClosestToPoint  traj1 = refited[1].trajectoryStateClosestToPoint(vert);
      nt.tag_refit_tuneP_transverse_IP = traj0.perigeeParameters().transverseImpactParameter();
      nt.probe_refit_tuneP_transverse_IP = traj1.perigeeParameters().transverseImpactParameter();
      GlobalPoint vert2(nt.pv_x, nt.pv_y, nt.pv_z);
      traj0 = refited[0].trajectoryStateClosestToPoint(vert2);
      traj1 = refited[1].trajectoryStateClosestToPoint(vert2);
      nt.tag_tuneP_transverse_IP = traj0.perigeeParameters().transverseImpactParameter();
      nt.probe_tuneP_transverse_IP = traj1.perigeeParameters().transverseImpactParameter();
    } else {
      nt.pair_tuneP_svprob = -1;
      nt.pair_tuneP_fit_mass = -1;
      nt.pair_tuneP_normalchi2 = -1;
      nt.refit_tuneP_vtx_x = -99;
      nt.refit_tuneP_vtx_y = -99;
      nt.refit_tuneP_vtx_z = -99;
      nt.tag_refit_tuneP_transverse_IP = -99;
      nt.probe_refit_tuneP_transverse_IP = -99;
      nt.tag_tuneP_transverse_IP = -99;
      nt.probe_tuneP_transverse_IP = -99;
    }
  } else {
    if (status_) {
      nt.pair_svprob = prob_;
      nt.pair_fit_mass = DimuonMass(refited[0].track().pt(),
                                    refited[0].track().eta(),
                                    refited[0].track().phi(),
                                    refited[1].track().pt(),
                                    refited[1].track().eta(),
                                    refited[1].track().phi());
      nt.pair_normalchi2 = normalchi2_;
      nt.refit_vtx_x = dimuvtx.position().x();
      nt.refit_vtx_y =dimuvtx.position().y();
      nt.refit_vtx_z = dimuvtx.position().z();
      GlobalPoint vert(nt.refit_vtx_x, nt.refit_vtx_y, nt.refit_vtx_z);
      TrajectoryStateClosestToPoint  traj0 = refited[0].trajectoryStateClosestToPoint(vert);
      TrajectoryStateClosestToPoint  traj1 = refited[1].trajectoryStateClosestToPoint(vert);
      nt.tag_refit_transverse_IP = traj0.perigeeParameters().transverseImpactParameter();
      nt.probe_refit_transverse_IP = traj1.perigeeParameters().transverseImpactParameter();
      GlobalPoint vert2(nt.pv_x, nt.pv_y, nt.pv_z);
      traj0 = refited[0].trajectoryStateClosestToPoint(vert2);
      traj1 = refited[1].trajectoryStateClosestToPoint(vert2);
      nt.tag_transverse_IP = traj0.perigeeParameters().transverseImpactParameter();
      nt.probe_transverse_IP = traj1.perigeeParameters().transverseImpactParameter();
    } else {
      nt.pair_svprob = -1;
      nt.pair_fit_mass = -1;
      nt.pair_normalchi2 = -1;
      nt.refit_vtx_x = -99;
      nt.refit_vtx_y = -99;
      nt.refit_vtx_z = -99;
      nt.tag_refit_transverse_IP = -99;
      nt.probe_refit_transverse_IP = -99;
      nt.tag_transverse_IP = -99;
      nt.probe_transverse_IP = -99;
    }
  }
}


// void Analysis_KlFitter::fillNtuple_genFinal(Analysis_NtupleContent &nt) {
  // std::cout << " genFinal Fill ntuple " << std::endl;
    // if (genFinal_status_ ) {
      // nt.genFinal_pair_svprob = genFinal_prob_;
      // nt.genFinal_pair_fit_mass = DimuonMass(genFinal_refited[0].track().pt(),
                                    // genFinal_refited[0].track().eta(),
                                    // genFinal_refited[0].track().phi(),
                                    // genFinal_refited[1].track().pt(),
                                    // genFinal_refited[1].track().eta(),
                                    // genFinal_refited[1].track().phi());
      // nt.genFinal_pair_normalchi2 = genFinal_normalchi2_;
      // nt.genFinal_vtx_x = genFinal_dimuvtx.position().x();
      // nt.genFinal_vtx_y =genFinal_dimuvtx.position().y();
      // nt.genFinal_vtx_z = genFinal_dimuvtx.position().z();
      // GlobalPoint vert(nt.genFinal_vtx_x, nt.genFinal_vtx_y, nt.genFinal_vtx_z);
      // TrajectoryStateClosestToPoint  traj0 = genFinal_refited[0].trajectoryStateClosestToPoint(vert);
      // TrajectoryStateClosestToPoint  traj1 = genFinal_refited[1].trajectoryStateClosestToPoint(vert);
      // nt.genFinal_refit_vtx_gentrk1_IP = traj0.perigeeParameters().transverseImpactParameter();
      // nt.genFinal_refit_vtx_gentrk2_IP = traj1.perigeeParameters().transverseImpactParameter();
      // // GlobalPoint vert2(nt.pv_x, nt.pv_y, nt.pv_z);
      // // traj0 = genFinal_refited[0].trajectoryStateClosestToPoint(vert2);
      // // traj1 = genFinal_refited[1].trajectoryStateClosestToPoint(vert2);
      // traj0 = genFinal_prefit[0].trajectoryStateClosestToPoint(vert);
      // traj1 = genFinal_prefit[1].trajectoryStateClosestToPoint(vert);
      // nt.genFinal_vtx_gentrk1_IP = traj0.perigeeParameters().transverseImpactParameter();
      // nt.genFinal_vtx_gentrk2_IP = traj1.perigeeParameters().transverseImpactParameter();
      // if(genDiTau_status_){
      // traj0 = genDiTau_refited[0].trajectoryStateClosestToPoint(vert);
      // traj1 = genDiTau_refited[1].trajectoryStateClosestToPoint(vert);
      // nt.genFinal_refit_vtx_gentau1_IP = traj0.perigeeParameters().transverseImpactParameter();
      // nt.genFinal_refit_vtx_gentau2_IP = traj1.perigeeParameters().transverseImpactParameter();           
      // traj0 = genDiTau_prefit[0].trajectoryStateClosestToPoint(vert);
      // traj1 = genDiTau_prefit[1].trajectoryStateClosestToPoint(vert);
      // nt.genFinal_vtx_gentau1_IP = traj0.perigeeParameters().transverseImpactParameter();
      // nt.genFinal_vtx_gentau2_IP = traj1.perigeeParameters().transverseImpactParameter();
      // }else{
      // nt.genFinal_refit_vtx_gentau1_IP = -99;
      // nt.genFinal_refit_vtx_gentau2_IP = -99;
      // nt.genFinal_vtx_gentau1_IP = -99;
      // nt.genFinal_vtx_gentau2_IP = -99;               
      // }      
    // } else {
      // nt.genFinal_pair_svprob = -1;
      // nt.genFinal_pair_fit_mass = -1;
      // nt.genFinal_pair_normalchi2 = -1;
      // nt.genFinal_vtx_x = -99;
      // nt.genFinal_vtx_y = -99;
      // nt.genFinal_vtx_z = -99;
      // nt.genFinal_refit_vtx_gentrk1_IP = -99;
      // nt.genFinal_refit_vtx_gentrk2_IP = -99;
      // nt.genFinal_vtx_gentrk1_IP = -99;
      // nt.genFinal_vtx_gentrk1_IP = -99;
      // nt.genFinal_refit_vtx_gentau1_IP = -99;
      // nt.genFinal_refit_vtx_gentau2_IP = -99;
      // nt.genFinal_vtx_gentau1_IP = -99;
      // nt.genFinal_vtx_gentau2_IP = -99;  
    // }
// }


// void Analysis_KlFitter::fillNtuple_genDiTau(Analysis_NtupleContent &nt) {
  // std::cout << " Genditau fill ntuple " << std::endl;
    // if (genDiTau_status_ ) {
      // nt.genDiTau_pair_svprob = genDiTau_prob_;
      // nt.genDiTau_pair_fit_mass = DimuonMass(genDiTau_refited[0].track().pt(),
                                    // genDiTau_refited[0].track().eta(),
                                    // genDiTau_refited[0].track().phi(),
                                    // genDiTau_refited[1].track().pt(),
                                    // genDiTau_refited[1].track().eta(),
                                    // genDiTau_refited[1].track().phi());
      // nt.genDiTau_pair_normalchi2 = genDiTau_normalchi2_;
      // nt.genDiTau_vtx_x = genDiTau_dimuvtx.position().x();
      // nt.genDiTau_vtx_y =genDiTau_dimuvtx.position().y();
      // nt.genDiTau_vtx_z = genDiTau_dimuvtx.position().z();
      // GlobalPoint vert(nt.genDiTau_vtx_x, nt.genDiTau_vtx_y, nt.genDiTau_vtx_z);
      // TrajectoryStateClosestToPoint  traj0 = genDiTau_refited[0].trajectoryStateClosestToPoint(vert);
      // TrajectoryStateClosestToPoint  traj1 = genDiTau_refited[1].trajectoryStateClosestToPoint(vert);
      // nt.genDiTau_refit_vtx_gentau1_IP = traj0.perigeeParameters().transverseImpactParameter();
      // nt.genDiTau_refit_vtx_gentau2_IP = traj1.perigeeParameters().transverseImpactParameter();
      // // GlobalPoint vert2(nt.pv_x, nt.pv_y, nt.pv_z);
      // // traj0 = genDiTau_refited[0].trajectoryStateClosestToPoint(vert2);
      // // traj1 = genDiTau_refited[1].trajectoryStateClosestToPoint(vert2);
      // traj0 = genDiTau_prefit[0].trajectoryStateClosestToPoint(vert);
      // traj1 = genDiTau_prefit[1].trajectoryStateClosestToPoint(vert);
      // nt.genDiTau_vtx_gentau1_IP = traj0.perigeeParameters().transverseImpactParameter();
      // nt.genDiTau_vtx_gentau2_IP = traj1.perigeeParameters().transverseImpactParameter();
      // if(genFinal_status_){
      // traj0 = genFinal_refited[0].trajectoryStateClosestToPoint(vert);
      // traj1 = genFinal_refited[1].trajectoryStateClosestToPoint(vert);
      // nt.genFinal_refit_vtx_gentrk1_IP = traj0.perigeeParameters().transverseImpactParameter();
      // nt.genFinal_refit_vtx_gentrk2_IP = traj1.perigeeParameters().transverseImpactParameter();           
      // traj0 = genFinal_prefit[0].trajectoryStateClosestToPoint(vert);
      // traj1 = genFinal_prefit[1].trajectoryStateClosestToPoint(vert);
      // nt.genFinal_vtx_gentrk1_IP = traj0.perigeeParameters().transverseImpactParameter();
      // nt.genFinal_vtx_gentrk2_IP = traj1.perigeeParameters().transverseImpactParameter();
      // }else{
      // nt.genFinal_refit_vtx_gentrk1_IP = -99;
      // nt.genFinal_refit_vtx_gentrk2_IP = -99;
      // nt.genFinal_vtx_gentrk1_IP = -99;
      // nt.genFinal_vtx_gentrk2_IP = -99;               
      // }      
    // } else {
      // nt.genDiTau_pair_svprob = -1;
      // nt.genDiTau_pair_fit_mass = -1;
      // nt.genDiTau_pair_normalchi2 = -1;
      // nt.genDiTau_vtx_x = -99;
      // nt.genDiTau_vtx_y = -99;
      // nt.genDiTau_vtx_z = -99;
      // nt.genDiTau_refit_vtx_gentrk1_IP = -99;
      // nt.genDiTau_refit_vtx_gentrk2_IP = -99;
      // nt.genDiTau_vtx_gentrk1_IP = -99;
      // nt.genDiTau_vtx_gentrk1_IP = -99;
      // nt.genDiTau_refit_vtx_gentau1_IP = -99;
      // nt.genDiTau_refit_vtx_gentau2_IP = -99;
      // nt.genDiTau_vtx_gentau1_IP = -99;
      // nt.genDiTau_vtx_gentau1_IP = -99;  
    // }
// }


