//
// Original Author:
//                george karathanasis
//         Created:  Thu, 23 Mar 2019 17:40:23 GMT
//
// dimuon vertex fitter class

#ifndef ANALYSIS_KLFITTER_H
#define ANALYSIS_KLFITTER_H

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "Analysis_NtupleContent.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "helper.h"

#include "TMath.h"
#include "Math/Vector3D.h"
#include "TMatrixD.h"
#include "TVectorD.h"



class Analysis_KlFitter {
public:
  typedef ROOT::Math::XYZVector XYZVector;
  Analysis_KlFitter(std::vector<reco::TransientTrack> &);

  ~Analysis_KlFitter();

  void  genFinal(Analysis_NtupleContent &nt);
  void  genDiTau(Analysis_NtupleContent &nt);
  double distancePointLine(XYZVector &V1, XYZVector &V2){
     double v1_mag2 = V1.Mag2();
     double v2_mag2 = V2.Mag2();
     double v1dv2 = V1.Dot(V2);
     return std::sqrt(v1_mag2 - v1dv2*v1dv2/v2_mag2);
     
  } 
  bool isZeroVector(XYZVector &V1, double eps){
     bool XisZero = V1.X() < eps;
     bool YisZero = V1.Y() < eps;
     bool ZisZero = V1.Z() < eps;
     if(XisZero && YisZero && ZisZero) return true;
     else return false;
  }
  
  void fillNtuple(Analysis_NtupleContent &nt, bool isTuneP = false);
  bool status() { return status_; }
  float prob() { return prob_; }
  float normalchi2() { return normalchi2_; }
  float dz_PV_SV(float pv_z) { return abs(dimuvtx.position().z() - pv_z); }


  void genDiTau_genFinal_vertex(Analysis_NtupleContent &nt);
  // void fillNtuple_genFinal(Analysis_NtupleContent &nt);
  bool genFinal_status() { return genFinal_status_; }
  
  // void fillNtuple_genDiTau(Analysis_NtupleContent &nt);
  bool genDiTau_status() { return genDiTau_status_; }
private:
  TransientVertex dimuvtx;
  std::vector<reco::TransientTrack> prefit;
  std::vector<reco::TransientTrack> refited;
  bool status_ = true;
  float prob_ = 0;
  float normalchi2_ = 0;
  
  
  XYZVector gentrk1_x;
  XYZVector gentrk1_p;
  XYZVector gentrk2_x;
  XYZVector gentrk2_p; 
  XYZVector gentrk1_refit_vtx;
  XYZVector gentrk2_refit_vtx;
  XYZVector genFinal_vtx;
  bool genFinal_status_ = true;
  
  XYZVector gentau1_x;
  XYZVector gentau1_p;
  XYZVector gentau2_x;
  XYZVector gentau2_p; 
  XYZVector gentau1_refit_vtx;
  XYZVector gentau2_refit_vtx;
  XYZVector genDiTau_vtx;
  bool genDiTau_status_ = true;
};

#endif
