#pragma once
#ifndef CORRECTMEANMATERIAL_H_
#define CORRECTMEANMATERIAL_H_
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
#include "material_part_utils.h"
using namespace ROOT::Math;
using namespace utils;


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
         TMatrixD &P,
         Double_t &dErec,
         int dir,
         std::string MS,
         Double_t xOverX0){

  const Double_t kBGStop=0.02;
  Double_t mass2=mass*mass;
  //p*=q;
  if ((p/mass)<kBGStop) return kFALSE;
  Double_t p2=p*p;
  Double_t Ein=TMath::Sqrt(p2+mass2);
  Double_t dP= dPdxEulerStep(p,mass,xTimesRho,stepFraction,bg,kp0,kp1,kp2,kp3,kp4)*dir;
  if (dP==0) return kFALSE;
  Double_t pOut=p+dP;
  //if(dir<0) std::cout<<"dir:"<<dir<<" dP:"<<dP<<std::endl;
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

  if (xOverX0 != 0 && MS=="addMS_Smearing_Corr") {
    //Double_t theta2=1.0259e-6*14*14/28/(beta2*p2)*TMath::Abs(d)*9.36*2.33;
    Double_t theta2=0.0136*0.0136/(beta2*p2)*TMath::Abs(xOverX0);
    
    double lt = 1+0.038*TMath::Log(TMath::Abs(xOverX0));
    if (lt>0) theta2 *= lt*lt;
    
    //theta2 *= q*q;    // q=2 particle
    if(theta2>TMath::Pi()*TMath::Pi()) return kFALSE;
    cC22 = theta2*((1.-parvec[2])*(1.+parvec[2]))*(1. + parvec[3]*parvec[3]);
    cC33 = theta2*(1. + parvec[3]*parvec[3])*(1. + parvec[3]*parvec[3]);
    cC43 = theta2*parvec[3]*parvec[4]*(1. + parvec[3]*parvec[3]);
    cC44 = theta2*parvec[3]*parvec[3]*parvec[4]*parvec[4];
  }

  //Calculating the energy loss corrections************************
  Double_t cP4=1.;
  if ((xTimesRho != 0.) && (beta2 < 1.)) {
    Double_t dE=Eout-Ein;
    dErec=dE;
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
  P[3][4] += cC43;
  P[4][4] += cC44;
  parvec[4]  *= cP4;
  //std::cout<<"dir="<<dir<<" cP4="<<cP4<<" cC44= "<<cC44<<std::endl;
  //CheckCovariance();
  return kTRUE;
}


Bool_t CorrectForMeanMaterial(Double_t xTimesRho, Double_t mass, Double_t stepFraction, Double_t p, float bg,
         float kp0,
         float kp1,
         float kp2,
         float kp3,
         float kp4,
         Double_t &sinphi,
         Double_t &tanlambda,
         Double_t &invpT,
         Double_t &dErec,
         std::string Energy_smear,
         std::string MS,
         Double_t xOverX0){
           
  const Double_t kBGStop=0.02;
  Double_t mass2=mass*mass;
  //p*=q;
  if ((p/mass)<kBGStop) return kFALSE;
  Double_t p2=p*p;
  Double_t Ein=TMath::Sqrt(p2+mass2);
  Double_t dP= dPdxEulerStep(p,mass,xTimesRho,stepFraction,bg,kp0,kp1,kp2,kp3,kp4);
  if(Energy_smear=="landau")
  {
    //std::cout<<"dP: "<<dP;
    Double_t sign = TMath::Sign(1,dP);
    dP=gRandom->Landau(abs(dP),0.15*abs(dP));
    if (dP<0.0001) dP=0.0001;
    dP*=sign;
    //std::cout<<" dPsmear: "<<dP<<std::endl;
  }
  if(Energy_smear=="gauss")
  {
    //std::cout<<"dP: "<<dP;
    Double_t sign = TMath::Sign(1,dP);
    dP=gRandom->Gaus(abs(dP),0.15*abs(dP));
    if (dP<0.0001) dP=0.0001;
    dP*=sign;
    //std::cout<<" dPsmear: "<<dP<<std::endl;
  }
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
  Double_t cC44 = 0.;
  if (xOverX0 != 0) {
    //Double_t theta2=1.0259e-6*14*14/28/(beta2*p2)*TMath::Abs(d)*9.36*2.33;
    Double_t theta2=0.0136*0.0136/(beta2*p2)*TMath::Abs(xOverX0);
    
    double lt = 1+0.038*TMath::Log(TMath::Abs(xOverX0));
    if (lt>0) theta2 *= lt*lt;
    
    //theta2 *= q*q;    // q=2 particle
    if(theta2>TMath::Pi()*TMath::Pi()) return kFALSE;
    cC22 = theta2*((1.-sinphi)*(1.+sinphi))*(1. + tanlambda*tanlambda);
    cC33 = theta2*(1. + tanlambda*tanlambda)*(1. + tanlambda*tanlambda);
    cC44 = theta2*tanlambda*tanlambda*invpT*invpT;
  }
  

  //Calculating the energy loss corrections************************
  Double_t cP4=1.;
  if ((xTimesRho != 0.) && (beta2 < 1.)) {
    Double_t dE=Eout-Ein;
    dErec=dE;
    if ( (1.+ dE/p2*(dE + 2*Ein)) < 0. ) return kFALSE;
    cP4 = 1./TMath::Sqrt(1.+ dE/p2*(dE + 2*Ein));  //A precise formula by Ruben !
    //if (TMath::Abs(fP4*cP4)>100.) return kFALSE; //Do not track below 10 MeV/c -dsiable controlled by the BG cut
    // Approximate energy loss fluctuation (M.Ivanov)
    const Double_t knst=0.07; // To be tuned.
    Double_t sigmadE=knst*TMath::Sqrt(TMath::Abs(dE));
    cC44 += ((sigmadE*Ein/p2*invpT)*(sigmadE*Ein/p2*invpT));
  }

  //Applying the corrections*****************************
  if (MS=="addMS_Smearing" || MS=="addMS_Smearing_Corr"){
    const float kMaxP3=0.5;
    const float kMaxP4=0.3;
    if (TMath::Sqrt(cC44)>kMaxP4*TMath::Abs(invpT)) return kFALSE;
    Float_t p2New=sinphi+gRandom->Gaus(0,TMath::Sqrt(cC22));
    Float_t dp3New=gRandom->Gaus(0,TMath::Sqrt(cC33));
    if (TMath::Abs(p2New)>1.) return kFALSE;
    if (TMath::Abs(dp3New)>kMaxP3) return kFALSE;
    sinphi=p2New;
    tanlambda+=dp3New;
    invpT+=gRandom->Gaus(0,TMath::Sqrt(cC44));
  }
  
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


Bool_t Propagate(double dz, TMatrixD &PPred, TMatrix P, TVectorD & predstep, TVectorD parvec, 
                 Double_t kAlmost0, Double_t kAlmost1, int fPrintLevel, int dir, Bool_t InPlane, 
                 Double_t& dErec ,Double_t &dxyzrec, Bool_t Energy_loss, std::string CorrTime,
                 std::string MS,Double_t xx0)
{
            ////Prepare all the parameters for the prediction
          /*
          double ZA =0.54141;
          double Ipar=64.7e-9;      ///GeV
          double rho= 1.032;   //g/cm^3
          double X1=2.49;
          double X0=0.1469;
          double muon_mass=0.1056583755; //GeV/c^2
          double a=0.1610;
          double m=3.24;
          double hw=21.54e-9;
          double Z=0.085+6*0.915;
          */

          double curvature = (B*0.3e-2)*parvec[4];
          double sinphi = parvec[2];
          double tanlambda = parvec[3];

          Double_t z2r = curvature*dz;
          //std::cout<<"check0"<<std::endl;
          Double_t f1=sinphi, f2=f1 + z2r;
          //std::cout<<"check 0.1"<<std::endl;
          if (TMath::Abs(f1) >= kAlmost1) 
          {
            std::cout<<"f1= "<<f1<<std::endl;
            return false;
          }
          //std::cout<<"check 0.2"<<std::endl;
          if (TMath::Abs(f2) >= kAlmost1) 
          {
            std::cout<<"f2= "<<f2<<std::endl;
            return false;
          }
          if (TMath::Abs(tanlambda)< kAlmost0) 
          {
            std::cout<<"tanlambda= "<<tanlambda<<std::endl;
            return false;
          }

          //std::cout<<"check1"<<std::endl;
          
          
          Double_t r1=TMath::Sqrt((1.-f1)*(1.+f1)), r2=TMath::Sqrt((1.-f2)*(1.+f2));
          if (TMath::Abs(r1)<kAlmost0)  
          {
            std::cout<<"r1= "<<r1<<std::endl;
            return false;
          }
          if (TMath::Abs(r2)<kAlmost0)  
          {
            std::cout<<"r2= "<<r1<<std::endl;
            return false;
          }

          //std::cout<<"check2"<<std::endl;
 
          double dy2dz = (f1+f2)/(r1+r2);
          double rot = TMath::ASin(r1*f2 - r2*f1); 
            if (f1*f1+f2*f2>1 && f1*f2<0) {          // special cases of large rotations or large abs angles
              if (f2>0) rot =  TMath::Pi() - rot;    //
              else      rot = -TMath::Pi() - rot;
            }

          //std::cout<<"check3"<<std::endl;

          

          // predicted step
          if (fPrintLevel > 1)
            {
              std::cout << "P Matrix: " << std::endl;
              P.Print();
            }
          
          //std::cout<<"check4"<<std::endl;
          predstep = parvec;

          if(fPrintLevel>0) std::cout<<"z2r: "<<z2r<<std::endl;
          
          predstep[0] += dz*dy2dz;  // update y
          predstep[2] += z2r;        // update sinphi
          predstep[1] +=tanlambda/curvature*rot; //update x
                                 // update tree values
          
          
          
          
          
          // equations from the extended Kalman filter
          //f = F - 1
          /*
          Double_t f02=    dx/(r1*r1*r1);            Double_t cc=crv/fP4;
          Double_t f04=0.5*dx*dx/(r1*r1*r1);         f04*=cc;
          Double_t f12=    dx*fP3*f1/(r1*r1*r1);
          Double_t f14=0.5*dx*dx*fP3*f1/(r1*r1*r1);  f14*=cc;
          Double_t f13=    dx/r1;
          Double_t f24=    dx;                       f24*=cc;
          */
          Double_t rinv = 1./r1;
          Double_t r3inv = rinv*rinv*rinv;
          Double_t f24=    z2r/parvec[4];
          Double_t f02=    dz*r3inv;
          Double_t f04=0.5*f24*f02;
          Double_t f12=    f02*parvec[3]*f1;
          Double_t f14=0.5*f24*f02*parvec[3]*f1;
          Double_t f13=    dz*rinv;


          //b = C*ft
          Double_t b00=f02*P[2][0] + f04*P[4][0], b01=f12*P[2][0] + f14*P[4][0] + f13*P[3][0];
          Double_t b02=f24*P[4][0];
          Double_t b10=f02*P[2][1] + f04*P[4][1], b11=f12*P[2][1] + f14*P[4][1] + f13*P[3][1];
          Double_t b12=f24*P[4][1];
          Double_t b20=f02*P[2][2] + f04*P[4][2], b21=f12*P[2][2] + f14*P[4][2] + f13*P[3][2];
          Double_t b22=f24*P[4][2];
          Double_t b40=f02*P[4][2] + f04*P[4][4], b41=f12*P[4][2] + f14*P[4][4] + f13*P[4][3];
          Double_t b42=f24*P[4][4];
          Double_t b30=f02*P[3][2] + f04*P[4][3], b31=f12*P[3][2] + f14*P[4][3] + f13*P[3][3];
          Double_t b32=f24*P[4][3];
          

          //a = f*b = f*C*ft
          Double_t a00=f02*b20+f04*b40,a01=f02*b21+f04*b41,a02=f02*b22+f04*b42;
          Double_t a11=f12*b21+f14*b41+f13*b31,a12=f12*b22+f14*b42+f13*b32;
          Double_t a22=f24*b42;

          
          PPred=P;
          //F*C*Ft = C + (b + bt + a)
          PPred[0][0] += b00 + b00 + a00;
          PPred[1][0] += b10 + b01 + a01; 
          PPred[0][1] += b10 + b01 + a01;
          PPred[2][0] += b20 + b02 + a02;
          PPred[0][2] += b20 + b02 + a02;
          PPred[3][0] += b30;
          PPred[0][3] += b30;
          PPred[4][0] += b40;
          PPred[0][4] += b40;
          PPred[1][1] += b11 + b11 + a11;
          PPred[2][1] += b21 + b12 + a12;
          PPred[1][2] += b21 + b12 + a12;
          PPred[3][1] += b31;
          PPred[1][3] += b31; 
          PPred[4][1] += b41;
          PPred[1][4] += b41;
          PPred[2][2] += b22 + b22 + a22;
          PPred[3][2] += b32;
          PPred[2][3] += b32;
          PPred[4][2] += b42;
          PPred[2][4] += b42;

          
            float tanPhi2 = predstep[2]*predstep[2];
            tanPhi2/=(1-tanPhi2);
            double deltaxyz=dz*TMath::Sqrt(1.+tanPhi2+predstep[3]*predstep[3]);           
          
          dxyzrec=deltaxyz;
          //std::cout<<"deltaxyz dist"<<deltaxyz<<std::endl;

          //double deltaxyz = rot*dir /abs( curvature * cos(atan(tanlambda)));
          
          //std::cout<<"deltaxyz curv"<<deltaxyz<<std::endl;

          double p = sqrt(pow(parvec[3]/parvec[4],2)+pow(1/parvec[4],2));

          
          
          if(InPlane && Energy_loss && CorrTime!="after") 
          {
          //std::cout<<"deltaxyz Kalman:"<<deltaxyz<<std::endl;
          Bool_t checkstatus= CorrectForMeanMaterial(-dir*deltaxyz*rho,muon_mass,0.0005,p,(p/muon_mass),rho,X0,X1,Ipar,ZA,predstep,PPred,dErec,dir,MS,deltaxyz/xx0);
          if (fPrintLevel >0 )
          {
              std::cout << " Predstep: y " << predstep[0] << " x " << predstep[1] << " sinphi " << predstep[2] << " tanlambda " << predstep[3] << " 1/pT " << predstep[4] << std::endl;
          }
          //std::cout<<"I'm correcting"<<std::endl;
          return checkstatus;
          }
          else
          {
            if (fPrintLevel >0 )
          {
              std::cout << " Predstep: y " << predstep[0] << " x " << predstep[1] << " sinphi " << predstep[2] << " tanlambda " << predstep[3] << " 1/pT " << predstep[4] << std::endl;
          }

          }


          

          return 1.;
}

#endif