#include <iostream>
#include "TVectorD.h"
#include "TMatrix.h"
#include "TMath.h"
#include "TVector3.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TF2.h"
#include "Math/Vector3D.h"
#include <algorithm>
#include <vector>
using namespace ROOT::Math;

Double_t BetheBlochGeant(float bg,
         float kp0,
         float kp1,
         float kp2,
         float kp3,
         float kp4) {
  //
  // This is the parameterization of the Bethe-Bloch formula inspired by Geant.
  //
  // bg  - beta*gamma
  // kp0 - density [g/cm^3]
  // kp1 - density effect first junction point
  // kp2 - density effect second junction point
  // kp3 - mean excitation energy [GeV]
  // kp4 - mean Z/A
  
  //
  // The default values for the kp* parameters are for silicon. 
  // The returned value is in [GeV/(g/cm^2)].
  // 

  const float mK  = 0.307075e-3; // [GeV*cm^2/g]
  const float me  = 0.511e-3;    // [GeV/c^2]
  const float rho = kp0;
  const float x0  = kp1*2.303;
  const float x1  = kp2*2.303;
  const float mI  = kp3;
  const float mZA = kp4;
  const float bg2 = bg*bg;
  const float maxT= 2*me*bg2;    // neglecting the electron mass
  
  //*** Density effect
  float d2=0.; 
  const float x=TMath::Log(bg);
  const float lhwI=TMath::Log(28.816*1e-9*TMath::Sqrt(rho*mZA)/mI);
  if (x > x1) {
    d2 = lhwI + x - 0.5;
  } else if (x > x0) {
    const float r=(x1-x)/(x1-x0);
    d2 = lhwI + x - 0.5 + (0.5 - lhwI - x0)*r*r*r;
  }

  return mK*mZA*(1+bg2)/bg2*
         (0.5*TMath::Log(2*me*bg2*maxT/(mI*mI)) - bg2/(1+bg2) - d2);
}

Double_t dPdxEulerStep(double p, double mass,  Double_t xTimesRho, double step,float bg,
         float kp0,
         float kp1,
         float kp2,
         float kp3,
         float kp4){
    const Double_t kBGStop = 0.02;
   
    if (bg<kBGStop) return 0;
    Double_t dPdx=TMath::Abs(BetheBlochGeant(bg,kp0,kp1,kp2,kp3,kp4))*TMath::Sqrt(1.+1./(bg*bg));
    bg=p/mass;
    if (bg<kBGStop) return 0;
    Double_t dPdx2=TMath::Abs(BetheBlochGeant(bg,kp0,kp1,kp2,kp3,kp4))*TMath::Sqrt(1.+1./(bg*bg));
    //
    Int_t nSteps=1+(TMath::Abs(dPdx*xTimesRho)+TMath::Abs(dPdx2*xTimesRho))/step;
    if (nSteps==1) return 0.5*(dPdx+dPdx2)*xTimesRho;
    Float_t xTimesRhoS=xTimesRho/nSteps;
    Float_t sumP=0;
    for (Int_t i=0; i<nSteps;i++){
      p+=dPdx*xTimesRhoS;
      sumP+=dPdx*xTimesRhoS;
      bg=p/mass;
      if (bg<kBGStop) return 0;
      dPdx=TMath::Abs(BetheBlochGeant(bg,kp0,kp1,kp2,kp3,kp4))*TMath::Sqrt(1.+1./(bg*bg));
    }
    return sumP;
};

Bool_t CorrectForMeanMaterial(Double_t xTimesRho, Double_t mass, Float_t stepFraction, Double_t p, float bg,
         float kp0,
         float kp1,
         float kp2,
         float kp3,
         float kp4,
         TVectorD &parvec,
         TMatrixD &P){
  const Double_t kBGStop=0.02;
  Double_t mass2=mass*mass;
  //p*=q;
  if ((p/mass)<kBGStop) return kFALSE;
  Double_t p2=p*p;
  Double_t Ein=TMath::Sqrt(p2+mass2);
  Double_t dP= dPdxEulerStep(p,mass,xTimesRho,stepFraction,bg,kp0,kp1,kp2,kp3,kp4);
  if (dP==0) return kFALSE;
  Double_t pOut=p+dP;
  if ((pOut/mass)<kBGStop) return kFALSE;
  Double_t Eout=TMath::Sqrt(pOut*pOut+mass2);
  p=(p+pOut)*0.5;
  // Use mean values for p2, p and beta2
  p2=p*p;
  Double_t beta2=p2/(p2+mass2);
  //
  
  //
  //Calculating the multiple scattering corrections******************
  Double_t cC22 = 0.;
  Double_t cC33 = 0.;
  Double_t cC43 = 0.;
  Double_t cC44 = 0.;

  /*
  if (xOverX0 != 0) {
    //Double_t theta2=1.0259e-6*14*14/28/(beta2*p2)*TMath::Abs(d)*9.36*2.33;
    Double_t theta2=0.0136*0.0136/(beta2*p2)*TMath::Abs(xOverX0);
    if (GetUseLogTermMS()) {
      double lt = 1+0.038*TMath::Log(TMath::Abs(xOverX0));
      if (lt>0) theta2 *= lt*lt;
    }
    theta2 *= q*q;    // q=2 particle
    if(theta2>TMath::Pi()*TMath::Pi()) return kFALSE;
    cC22 = theta2*((1.-fP2)*(1.+fP2))*(1. + fP3*fP3);
    cC33 = theta2*(1. + fP3*fP3)*(1. + fP3*fP3);
    cC43 = theta2*fP3*fP4*(1. + fP3*fP3);
    cC44 = theta2*fP3*fP4*fP3*fP4;
  }
  */

  //Calculating the energy loss corrections************************
  Double_t cP4=1.;
  if ((xTimesRho != 0.) && (beta2 < 1.)) {
    Double_t dE=Eout-Ein;
    if ( (1.+ dE/p2*(dE + 2*Ein)) < 0. ) return kFALSE;
    cP4 = 1./TMath::Sqrt(1.+ dE/p2*(dE + 2*Ein));  //A precise formula by Ruben !
    //if (TMath::Abs(fP4*cP4)>100.) return kFALSE; //Do not track below 10 MeV/c -dsiable controlled by the BG cut
    // Approximate energy loss fluctuation (M.Ivanov)
    const Double_t knst=0.07; // To be tuned.
    Double_t sigmadE=knst*TMath::Sqrt(TMath::Abs(dE));
    cC44 += ((sigmadE*Ein/p2*parvec[4])*(sigmadE*Ein/p2*parvec[4]));
  }

  //Applying the corrections*****************************
  P[2][2] += cC22;
  P[3][3] += cC33;
  P[4][3] += cC43;
  P[4][4] += cC44;
  parvec[4]  *= cP4;
  //CheckCovariance();
  return kTRUE;
}


Bool_t CorrectForMeanMaterial(Double_t xTimesRho, Double_t mass, Double_t stepFraction, Double_t p, float bg,
         float kp0,
         float kp1,
         float kp2,
         float kp3,
         float kp4,
         Double_t &invpT){
  const Double_t kBGStop=0.02;
  Double_t mass2=mass*mass;
  //p*=q;
  if ((p/mass)<kBGStop) return kFALSE;
  Double_t p2=p*p;
  Double_t Ein=TMath::Sqrt(p2+mass2);
  Double_t dP= dPdxEulerStep(p,mass,xTimesRho,stepFraction,bg,kp0,kp1,kp2,kp3,kp4);
  if (dP==0) return kFALSE;
  Double_t pOut=p+dP;
  if ((pOut/mass)<kBGStop) return kFALSE;
  Double_t Eout=TMath::Sqrt(pOut*pOut+mass2);
  p=(p+pOut)*0.5;
  // Use mean values for p2, p and beta2
  p2=p*p;
  Double_t beta2=p2/(p2+mass2);
  //
  
  //
  //Calculating the multiple scattering corrections******************


  /*
  if (xOverX0 != 0) {
    //Double_t theta2=1.0259e-6*14*14/28/(beta2*p2)*TMath::Abs(d)*9.36*2.33;
    Double_t theta2=0.0136*0.0136/(beta2*p2)*TMath::Abs(xOverX0);
    if (GetUseLogTermMS()) {
      double lt = 1+0.038*TMath::Log(TMath::Abs(xOverX0));
      if (lt>0) theta2 *= lt*lt;
    }
    theta2 *= q*q;    // q=2 particle
    if(theta2>TMath::Pi()*TMath::Pi()) return kFALSE;
    cC22 = theta2*((1.-fP2)*(1.+fP2))*(1. + fP3*fP3);
    cC33 = theta2*(1. + fP3*fP3)*(1. + fP3*fP3);
    cC43 = theta2*fP3*fP4*(1. + fP3*fP3);
    cC44 = theta2*fP3*fP4*fP3*fP4;
  }
  */

  //Calculating the energy loss corrections************************
  Double_t cP4=1.;
  if ((xTimesRho != 0.) && (beta2 < 1.)) {
    Double_t dE=Eout-Ein;
    if ( (1.+ dE/p2*(dE + 2*Ein)) < 0. ) return kFALSE;
    cP4 = 1./TMath::Sqrt(1.+ dE/p2*(dE + 2*Ein));  //A precise formula by Ruben !
    //if (TMath::Abs(fP4*cP4)>100.) return kFALSE; //Do not track below 10 MeV/c -dsiable controlled by the BG cut
    // Approximate energy loss fluctuation (M.Ivanov)
    const Double_t knst=0.07; // To be tuned.
    Double_t sigmadE=knst*TMath::Sqrt(TMath::Abs(dE));
  }

  //Applying the corrections*****************************
  
  invpT  *= cP4;
  //CheckCovariance();
  return kTRUE;
}



template<typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}

template<typename Iter>
Iter select_randomly(Iter start, Iter end) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
}