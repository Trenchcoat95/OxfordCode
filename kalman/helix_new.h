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

void makeSeed(const std::vector<XYZVector>  TPCClusters,
              XYZVector  &TPCClustersSeed,
              double &curvature_init,
              double &tanlambda_init,
              double &sinphi_init,
              double dir,
              int printlevel,
              TMatrixD &P,
              double sxy)
{

  //std::cout<<"I'm in function"<<std::endl;

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

  //P[0][0]=sxy2;        P[0][1]=0.;         P[0][2]=f20*sxy2;                                 P[0][3]=f30*sxy2;                                              P[0][4]=f40*sxy2;
  //P[1][0]=0.;          P[1][1]=sxy2;       P[1][2]=0.;                                       P[1][3]=f31*sxy2;                                              P[1][4]=0.;
  //P[2][0]=f20*sxy2;    P[2][1]=0.;         P[2][2]=f20*sxy2*f20+f22*sxy2*f22+f23*sxy2*f23;   P[2][3]=f30*sxy2*f20+f32*sxy2*f22;                             P[2][4]=f40*sxy2*f20+f42*sxy2*f22+f43*sxy2*f23;
  //P[3][0]=f30*sxy2;    P[3][1]=f31*sxy2;   P[3][2]=f30*sxy2*f20+f32*sxy2*f22;                P[3][3]=f30*sxy2*f30+f31*sxy2*f31+f32*sxy2*f32+f34*sxy2*f34;   P[3][4]=f30*sxy2*f40+f32*sxy2*f42;
  //P[4][0]=f40*sxy2;    P[4][1]=0.;         P[4][2]=f40*sxy2*f20+f42*sxy2*f22+f43*sxy2*f23;   P[4][3]=f30*sxy2*f40+f32*sxy2*f42;                             P[4][4]=f40*sxy2*f40+f42*sxy2*f42+f43*sxy2*f43;
  
  //For now only consider diagonal elements
  P[0][0]=sxy2;
  P[1][1]=sxy2;
  P[2][2]=f20*sxy2*f20+f22*sxy2*f22+f23*sxy2*f23;
  P[3][3]=f30*sxy2*f30+f31*sxy2*f31+f32*sxy2*f32+f34*sxy2*f34; 
  P[4][4]=f40*sxy2*f40+f42*sxy2*f42+f43*sxy2*f43;

  P[4][4]/=(0.5*0.299792458e-2)*(0.5*0.299792458e-2); // transform to 1/pt

  //P[4][2]/=(0.5*0.299792458e-2);
  //P[4][0]/=(0.5*0.299792458e-2);

  if (printlevel>0)
    {
      std::cout << "phi calc: dz, dy " << xyz2[2]-xyz0[2] << " " <<  xyz2[1]-xyz0[1] << std::endl;
      std::cout << "initial curvature, sinphi, tanlambda: " << curvature_init << " " << sinphi_init << " " << tanlambda_init << std::endl;
    }
 
  
}