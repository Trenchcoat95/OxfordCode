#pragma once
#ifndef HELIXNEW_H_
#define HELIXNEW_H_
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
#include "correctMeanmaterial.h"
#include <algorithm>
#include <vector>
#include "material_part_utils.h"
using namespace ROOT::Math;
using namespace utils;

Double_t makeC(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3){
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature
  //-----------------------------------------------------------------
  x3 -=x1;
  x2 -=x1;
  y3 -=y1;
  y2 -=y1;
  //  
  Double_t det = x3*y2-x2*y3;
  if (TMath::Abs(det)<1e-10){
    return 100;
  }
  //
  Double_t u = 0.5* (x2*(x2-x3)+y2*(y2-y3))/det;
  Double_t x0 = x3*0.5-y3*u;
  Double_t y0 = y3*0.5+x3*u;
  Double_t c2 = 1/TMath::Sqrt(x0*x0+y0*y0);
  if (det<0) c2*=-1;
  return c2;
}


Double_t makeSnp(Double_t x1,Double_t y1, Double_t x2,Double_t y2, Double_t x3,Double_t y3){
  //-----------------------------------------------------------------
  // Initial approximation of the track snp at position x1
  //-----------------------------------------------------------------
  x3 -=x1;
  x2 -=x1;
  y3 -=y1;
  y2 -=y1;
  //  
  Double_t det = x3*y2-x2*y3;
  if (TMath::Abs(det)<1e-10) {
    return 100;
  }
  //
  Double_t u = 0.5* (x2*(x2-x3)+y2*(y2-y3))/det;
  Double_t x0 = x3*0.5-y3*u; 
  Double_t y0 = y3*0.5+x3*u;
  Double_t c2 = 1/TMath::Sqrt(x0*x0+y0*y0);
  if (det<0) c2*=-1;
  x0*=c2;  
  return x0;
}

//_____________________________________________________________________________
Double_t makeTgln(Double_t x1,Double_t y1, Double_t x2,Double_t y2,Double_t z1,Double_t z2,Double_t c){
  //-----------------------------------------------------------------
  // Initial approximation of the tangent of the track dip angle
  //-----------------------------------------------------------------
  Double_t d  =  TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
  if (TMath::Abs(d*c*0.5)>1) return 0;
  Double_t   angle2    = asin(d*c*0.5);

  angle2  = (z1-z2)*c/(angle2*2.);    //dz /(R*dPhi)
  return angle2;
  //return (z1 - z2)/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

Bool_t SeedMaterialCorrection(Double_t xTimesRho, Double_t mass, Float_t stepFraction, Double_t p, float bg,
         float kp0,
         float kp1,
         float kp2,
         float kp3,
         float kp4,
         double &curvature_init,
         double &tanlambda_init,
         double &sinphi_init,
         TMatrixD &P,
         int dir,
         std::string Helix_Corr,
         Double_t xOverX0)
              
{
 const Double_t kBGStop=0.02;
  Double_t mass2=mass*mass;
  //p*=q;
  if ((p/mass)<kBGStop) return kFALSE;
  Double_t p2=p*p;
  Double_t Ein=TMath::Sqrt(p2+mass2);


  Double_t dP= dPdxEulerStep(p,mass,xTimesRho,stepFraction,bg,kp0,kp1,kp2,kp3,kp4);

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
  double invpTinit= curvature_init/(0.5*0.3e-2);
  //
  //Calculating the multiple scattering corrections******************

  Double_t cC22 = 0.;
  Double_t cC33 = 0.;
  //Double_t cC43 = 0.;
  Double_t cC44 = 0.;



  if (xOverX0 != 0 && Helix_Corr=="Eloss_MS") {
    //Double_t theta2=1.0259e-6*14*14/28/(beta2*p2)*TMath::Abs(d)*9.36*2.33;
    Double_t theta2=0.0136*0.0136/(beta2*p2)*TMath::Abs(xOverX0);
    
    double lt = 1+0.038*TMath::Log(TMath::Abs(xOverX0));
    if (lt>0) theta2 *= lt*lt;
    
    //theta2 *= q*q;    // q=2 particle
    if(theta2>TMath::Pi()*TMath::Pi()) return kFALSE;
    cC22 = theta2*((1.-sinphi_init)*(1.+sinphi_init))*(1. + tanlambda_init*tanlambda_init);
    cC33 = theta2*(1. + tanlambda_init*tanlambda_init)*(1. + tanlambda_init*tanlambda_init);
    //cC43 = theta2*parvec[3]*parvec[4]*(1. + tanlambda_init*tanlambda_init);
    cC44 = theta2*tanlambda_init*tanlambda_init*invpTinit*invpTinit;
  }


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
    cC44 += ((sigmadE*Ein/p2*invpTinit)*(sigmadE*Ein/p2*invpTinit));
    
  }



  //Applying the corrections*****************************


  P[2][2] += cC22;
  P[3][3] += cC33;
  //P[4][3] += cC43;
  P[4][4] += cC44;
  curvature_init  = invpTinit*cP4*0.5*0.3e-2;
  //if(dir>0)std::cout<<"dir="<<dir<<" cP4="<<cP4<<" cC44= "<<cC44<<std::endl;
  //CheckCovariance();

  return kTRUE;
}

double CalculatePath(const std::vector<XYZVector>  TPCClusters,
              double &curvature_init,
              double &tanlambda_init,
              double &sinphi_init
              )
{

  const int nplanes = 6;
  double Planes_Z[]={1179,1183,1199,1203,1219,1223,1239,1243,1339,1343,1539,1543};
  double Planes_Z_bkw[]={1539,1543,1339,1343,1239,1243,1219,1223,1199,1203,1179,1183};

  ////In development

  double xpos = TPCClusters[0].X();
  double ypos = TPCClusters[0].Y();
  double zpos = TPCClusters[0].Z();
  double sinphipos = sinphi_init;
  double zh;
  double dtot = 0;
  int fPrintLevel=0;
  if (fPrintLevel > 0)
            {
              std::cout << std::endl;
              std::cout << "Starting Position: x:" << TPCClusters[0].X() << " y: " << TPCClusters[0].Y() << " z: " << TPCClusters[0].Z() << std::endl;
            }

  for (size_t iTPCCluster=1; iTPCCluster<TPCClusters.size(); ++iTPCCluster)
        {

          if (fPrintLevel > 0)
            {
              std::cout << std::endl;
              std::cout << "Adding a new TPCCluster: x:" << TPCClusters[iTPCCluster].X() << " y: " << TPCClusters[iTPCCluster].Y() << " z: " << TPCClusters[iTPCCluster].Z() << std::endl;
            }

          zh=TPCClusters[iTPCCluster].Z();
          double dz=zh-zpos;
          for(size_t p=0;p<nplanes;p++)
            {
                if(abs(zh-zpos)<=Plane_thick&&abs(zh-zpos)>=0&&((dz>0 && zpos>=Planes_Z[p*2] && zpos<=Planes_Z[p*2+1]) || (dz<0 && zpos>=Planes_Z_bkw[p*2] && zpos<=Planes_Z_bkw[p*2+1])))
                  {
                   
                    double dzmid=zh-zpos;
                    
                   
                    
                    Double_t z2r = curvature_init*dzmid;
                    Double_t f1=sinphi_init;
                    Double_t f2=f1 + z2r;
                    Double_t r1=TMath::Sqrt((1.-f1)*(1.+f1)), r2=TMath::Sqrt((1.-f2)*(1.+f2));

                    Double_t dy2dz = (f1+f2)/(r1+r2);
                    Double_t rot = TMath::ASin(r1*f2 - r2*f1);
                    if (f1*f1+f2*f2>1 && f1*f2<0) {          // special cases of large rotations or large abs angles
                      if (f2>0) rot =  TMath::Pi() - rot;    //
                      else      rot = -TMath::Pi() - rot;
                    }

                    xpos+=tanlambda_init/curvature_init*rot;
                    ypos+=dzmid*dy2dz;
                    zpos+=dzmid;       
                    sinphipos+=z2r; 

                    if (fPrintLevel>0) std::cout<<zpos<<std::endl;

                    if (TMath::Abs(f1) >= kAlmost1 || TMath::Abs(f2) >= kAlmost1 || abs(tanlambda_init)<kAlmost0 || TMath::Abs(r1)<kAlmost0 || TMath::Abs(r2)<kAlmost0) 
                    {
                      std::cout<<TMath::Abs(f1)<<" "<<TMath::Abs(f2)<<" "<<tanlambda_init<<" "<<TMath::Abs(r1)<<" "<<TMath::Abs(r2)<<std::endl;
                      return -1;
                      
                    }

                    float tanPhi2 = sinphipos*sinphipos;
                    tanPhi2/=(1-tanPhi2);
                    dtot+=abs(dzmid*TMath::Sqrt(1.+tanPhi2+tanlambda_init*tanlambda_init)); 
                    if (fPrintLevel>0) std::cout<<"dtot: "<<dtot<<std::endl;


                                   
                    break;
                    
                    
                      
                  }
                

               else if((dz>0 && zpos>=Planes_Z[p*2] && zpos<=Planes_Z[p*2+1]) || (dz<0 && zpos>=Planes_Z_bkw[p*2] && zpos<=Planes_Z_bkw[p*2+1]))
                  {

                    double dzmid;
                    
                   
                    
                    if (dz>0) dzmid=Planes_Z[p*2+1]-zpos;
                    else dzmid=Planes_Z_bkw[p*2]-zpos;

                    Double_t z2r, f1, f2,r1, r2, dy2dz,rot = 0;
                    
                    if(dzmid!=0)
                    {
                    
                    z2r = curvature_init*dzmid;
                    f1=sinphi_init;
                    f2=f1 + z2r;
                    r1=TMath::Sqrt((1.-f1)*(1.+f1)), r2=TMath::Sqrt((1.-f2)*(1.+f2));

                    dy2dz = (f1+f2)/(r1+r2);
                    rot = TMath::ASin(r1*f2 - r2*f1);
                    if (f1*f1+f2*f2>1 && f1*f2<0) {          // special cases of large rotations or large abs angles
                      if (f2>0) rot =  TMath::Pi() - rot;    //
                      else      rot = -TMath::Pi() - rot;
                    }

                    xpos+=tanlambda_init/curvature_init*rot;
                    ypos+=dzmid*dy2dz;
                    zpos+=dzmid;       
                    sinphipos+=z2r; 
                    if (fPrintLevel>0) std::cout<<zpos<<std::endl;

                    if (TMath::Abs(f1) >= kAlmost1 || TMath::Abs(f2) >= kAlmost1 || abs(tanlambda_init)<kAlmost0 || TMath::Abs(r1)<kAlmost0 || TMath::Abs(r2)<kAlmost0) 
                    {
                      return -1;
                    }

                    float tanPhi2 = sinphipos*sinphipos;
                    tanPhi2/=(1-tanPhi2);
                    dtot+=abs(dzmid*TMath::Sqrt(1.+tanPhi2+tanlambda_init*tanlambda_init)); 
                    if (fPrintLevel>0) std::cout<<"dtot: "<<dtot<<std::endl;
                    }
                    
                    
                    

                    
                    if (dz>0) dzmid=Planes_Z[p*2+2]-Planes_Z[p*2+1];
                    else dzmid=Planes_Z_bkw[p*2+3]-Planes_Z_bkw[p*2];

                    z2r = curvature_init*dzmid;
                    f1=sinphi_init;
                    f2=f1 + z2r;
                    r1=TMath::Sqrt((1.-f1)*(1.+f1)), r2=TMath::Sqrt((1.-f2)*(1.+f2));

                    dy2dz = (f1+f2)/(r1+r2);
                    rot = TMath::ASin(r1*f2 - r2*f1);
                    if (f1*f1+f2*f2>1 && f1*f2<0) {          // special cases of large rotations or large abs angles
                      if (f2>0) rot =  TMath::Pi() - rot;    //
                      else      rot = -TMath::Pi() - rot;
                    }

                    xpos+=tanlambda_init/curvature_init*rot;
                    ypos+=dzmid*dy2dz;
                    zpos+=dzmid;       
                    sinphipos+=z2r; 
                    if (fPrintLevel>0) std::cout<<zpos<<std::endl;

                    if (TMath::Abs(f1) >= kAlmost1 || TMath::Abs(f2) >= kAlmost1 || abs(tanlambda_init)<kAlmost0 || TMath::Abs(r1)<kAlmost0 || TMath::Abs(r2)<kAlmost0) 
                    {
                      return -1;
                    }
                    


                    continue;
                    
                  }
            }
          
         
          
        }
    return dtot; 

}

void makeSeed(const std::vector<XYZVector>  TPCClusters,
              XYZVector  &TPCClustersSeed,
              double &curvature_init,
              double &tanlambda_init,
              double &sinphi_init,
              double dir,
              int printlevel,
              TMatrixD &P,
              double sxy,
              std::string Helix_Corr,
              size_t nCrossedPlanes)
{

  //std::cout<<"I'm in function"<<std::endl;
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
  double xx0 = 42.54;  //radiation length in cm (different from Alice's code where it's xx0=1/X0 in cm^-1)
  double Plane_thick = 4;

  size_t InitialTPNTPCClusters = 100;
  size_t nTPCClusters = TPCClusters.size();
  size_t firstTPCCluster = 0;
  size_t farTPCCluster = TMath::Min(nTPCClusters-1, InitialTPNTPCClusters);
  size_t intTPCCluster = farTPCCluster/2;
  size_t lastTPCCluster = nTPCClusters-1;

  double xyz0[3] = {TPCClusters.at(firstTPCCluster).X(),
                       TPCClusters.at(firstTPCCluster).Y(),
                       TPCClusters.at(firstTPCCluster).Z()};

  double xyz1[3] = {TPCClusters.at(intTPCCluster).X(),
                  TPCClusters.at(intTPCCluster).Y(),
                  TPCClusters.at(intTPCCluster).Z()};

  double xyz2[3] = {TPCClusters.at(farTPCCluster).X(),
                  TPCClusters.at(farTPCCluster).Y(),
                  TPCClusters.at(farTPCCluster).Z()};

  if (printlevel>1)
    {
      std::cout << "TPCCluster Dump in initial_trackpar_estimate: " << std::endl;
      for (size_t i=0;i<nTPCClusters;++i)
        {
          size_t ihf = i;
          std::cout << i << " : " <<
            TPCClusters.at(ihf).X() << " " <<
            TPCClusters.at(ihf).Y() << " " <<
            TPCClusters.at(ihf).Z() << std::endl;
        }
    }
  if (printlevel>0)
    {
      std::cout << "First TPCCluster x, y, z: " << xyz0[0] << " " << xyz0[1] << " " << xyz0[2] << std::endl;
      std::cout << "Inter TPCCluster x, y, z: " << xyz1[0] << " " << xyz1[1] << " " << xyz1[2] << std::endl;
      std::cout << "Far   TPCCluster x, y, z: " << xyz2[0] << " " << xyz2[1] << " " << xyz2[2] << std::endl;
    }

  
  Double_t sxy2=sxy*sxy;
  
  // calculate initial param
  TPCClustersSeed.SetXYZ(xyz0[0],xyz0[1],xyz0[2]);              
  sinphi_init=dir*makeSnp(xyz0[2],xyz0[1],xyz1[2],xyz1[1],xyz2[2],xyz2[1]);
  curvature_init=dir*makeC(xyz2[2],xyz2[1],xyz1[2],xyz1[1],xyz0[2],xyz0[1]);
  tanlambda_init=dir*makeTgln(xyz2[2],xyz2[1],xyz1[2],xyz1[1],xyz2[0],xyz1[0],curvature_init);
  //
  Double_t f40=(dir*makeC(xyz2[2],xyz2[1]+sxy,xyz1[2],xyz1[1],xyz0[2],xyz0[1])-curvature_init)/sxy;
  Double_t f42=(dir*makeC(xyz2[2],xyz2[1],xyz1[2],xyz1[1]+sxy,xyz0[2],xyz0[1])-curvature_init)/sxy;
  Double_t f43=(dir*makeC(xyz2[2],xyz2[1],xyz1[2],xyz1[1],xyz0[2],xyz0[1]+sxy)-curvature_init)/sxy;
  //
  Double_t f20=(dir*makeSnp(xyz0[2],xyz0[1]+sxy,xyz1[2],xyz1[1],xyz2[2],xyz2[1])-sinphi_init)/sxy;
  Double_t f22=(dir*makeSnp(xyz0[2],xyz0[1],xyz1[2],xyz1[1]+sxy,xyz2[2],xyz2[1])-sinphi_init)/sxy;
  Double_t f23=(dir*makeSnp(xyz0[2],xyz0[1],xyz1[2],xyz1[1],xyz2[2],xyz2[1]+sxy)-sinphi_init)/sxy;
  //
  Double_t f30=(dir*makeTgln(xyz2[2],xyz2[1]+sxy,xyz1[2],xyz1[1],xyz2[0],xyz1[0],curvature_init)-tanlambda_init)/sxy;
  Double_t f31=(dir*makeTgln(xyz2[2],xyz2[1],xyz1[2],xyz1[1],xyz2[0]+sxy,xyz1[0],curvature_init)-tanlambda_init)/sxy;
  Double_t f32=(dir*makeTgln(xyz2[2],xyz2[1],xyz1[2],xyz1[1]+sxy,xyz2[0],xyz1[0],curvature_init)-tanlambda_init)/sxy;
  Double_t f34=(dir*makeTgln(xyz2[2],xyz2[1],xyz1[2],xyz1[1],xyz2[0],xyz1[0]+sxy,curvature_init)-tanlambda_init)/sxy;
  P.Zero();

  /*
  P[0][0]=sxy2;        P[0][1]=0.;         P[0][2]=f20*sxy2;                                 P[0][3]=f30*sxy2;                                              P[0][4]=f40*sxy2;
  P[1][0]=0.;          P[1][1]=sxy2;       P[1][2]=0.;                                       P[1][3]=f31*sxy2;                                              P[1][4]=0.;
  P[2][0]=f20*sxy2;    P[2][1]=0.;         P[2][2]=f20*sxy2*f20+f22*sxy2*f22+f23*sxy2*f23;   P[2][3]=f30*sxy2*f20+f32*sxy2*f22;                             P[2][4]=f40*sxy2*f20+f42*sxy2*f22+f43*sxy2*f23;
  P[3][0]=f30*sxy2;    P[3][1]=f31*sxy2;   P[3][2]=f30*sxy2*f20+f32*sxy2*f22;                P[3][3]=f30*sxy2*f30+f31*sxy2*f31+f32*sxy2*f32+f34*sxy2*f34;   P[3][4]=f30*sxy2*f40+f32*sxy2*f42;
  P[4][0]=f40*sxy2;    P[4][1]=0.;         P[4][2]=f40*sxy2*f20+f42*sxy2*f22+f43*sxy2*f23;   P[4][3]=f30*sxy2*f40+f32*sxy2*f42;                             P[4][4]=f40*sxy2*f40+f42*sxy2*f42+f43*sxy2*f43;
  */
  //For now only consider diagonal elements
  P[0][0]=sxy2;
  P[1][1]=sxy2;
  P[2][2]=f20*sxy2*f20+f22*sxy2*f22+f23*sxy2*f23;
  P[3][3]=f30*sxy2*f30+f31*sxy2*f31+f32*sxy2*f32+f34*sxy2*f34; 
  P[4][4]=f40*sxy2*f40+f42*sxy2*f42+f43*sxy2*f43;

  P[4][4]/=(0.5*0.3e-2)*(0.5*0.3e-2); // transform to 1/pt



  //P[4][2]/=(0.5*0.3e-2);
  //P[4][0]/=(0.5*0.3e-2);

  if (printlevel>0)
    {
      std::cout << "phi calc: dz, dy " << xyz2[2]-xyz0[2] << " " <<  xyz2[1]-xyz0[1] << std::endl;
      std::cout << "initial curvature, phi, lambda: " << curvature_init << " " << asin(sinphi_init) << " " << atan(tanlambda_init) << std::endl;
    }

  double invpT = curvature_init/(0.5*0.3e-2);

  double p = sqrt(pow(tanlambda_init/(curvature_init/(0.5*0.3e-2)),2)+pow((0.5*0.3e-2)/curvature_init,2));


  double crossLength = CalculatePath(TPCClusters,curvature_init,tanlambda_init,sinphi_init);
  if(crossLength==-1) crossLength=nCrossedPlanes*Plane_thick;
  
  
  //std::cout<<"crossLength: "<<crossLength<<" dtot: "<<dtot<<std::endl;

  //float tanPhi2 = sinphi_init*sinphi_init;
  //tanPhi2/=(1-tanPhi2);
  //crossLength*=TMath::Sqrt(1.+tanPhi2+tanlambda_init*tanlambda_init);

  if (Helix_Corr == "Eloss_MS" || Helix_Corr == "Eloss") {
    SeedMaterialCorrection(nCrossedPlanes*Plane_thick*rho,muon_mass,0.0005,p,(p/muon_mass),rho,X0,X1,Ipar,ZA,
                                                                                curvature_init,tanlambda_init,sinphi_init,P,dir, Helix_Corr,crossLength/xx0);
  }
  
}


#endif