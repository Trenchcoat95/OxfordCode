#pragma once
#ifndef KALMAN_H_
#define KALMAN_H_
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
#include <ctime>
#include "material_part_utils.h"
using namespace ROOT::Math;
using namespace utils;

    // KalmanFit does a forwards or backwards Kalman fit using the sorted TPCCluster list
    // variables:  z is the independent variable
    // 0: y
    // 1: x
    // 2: sin(phi)
    // 3: tan(lambda)
    // 4: 1/pT (1/(GeV/c))

void KalmanFit(  
                   std::vector<double> &xht, 
                   std::vector<double> &yht, 
                   std::vector<double> &zht, 
                   std::vector<TVectorD> &parvect,
                   std::vector<TVectorD> &predstept,
                   std::vector<TMatrixD> &Pt,
                   std::vector<TMatrixD> &PPredt,
                   std::vector<TMatrixD> &Rt,
                   std::vector<double> &zpost, 
                   std::vector<XYZVector> xyz_plane,
                   double Ry,
                   double Rx,
                   double Ryx,
                   XYZVector  xyz_seed,
                   double  tanlambda_seed,
                   double curvature_seed,
                   double sinphi_seed,
                   double dir,
                   bool &status,
                   int fPrintLevel,
                   std::vector<std::vector<Double_t>> &dEreco,
                   std::vector<std::vector<Double_t>> &dxreco,
                   Bool_t Energy_loss,
                   std::string CorrTime,
                   Bool_t Fixed_Cov,
                   Bool_t Smear,
                   double xy_smear,
                   TMatrixD &P_seed,
                   std::string Seedtype,
                   std::string MS


                   

                                )
    {
     
      const int nplanes = 6;
      double Planes_Z[]={1179,1183,1199,1203,1219,1223,1239,1243,1339,1343,1539,1543};
      double Planes_Z_bkw[]={1539,1543,1339,1343,1239,1243,1219,1223,1199,1203,1179,1183};
      /*
      if (dir<0) 
      {
        Planes_Z[0]=1539;  Planes_Z[1]=1543;
        Planes_Z[2]=1339;  Planes_Z[3]=1343;
        Planes_Z[4]=1239;  Planes_Z[5]=1243;
        Planes_Z[6]=1219;  Planes_Z[7]=1223;
        Planes_Z[8]=1199;  Planes_Z[9]=1203;
        Planes_Z[10]=1179; Planes_Z[11]=1183;
      }
      */
      
        
      //for(size_t i=0;i<10;i++) std::cout<<Planes_Z[i]<<" ";
      //std::cout<<"\n";

      
      
      
      
      
      
      if(fPrintLevel>0)
      {
        std::cout << " Seed: y " << xyz_seed.Y() << " x " << xyz_seed.X() <<"z:"<<xyz_seed.Z()<< " sinphi " << sinphi_seed << " tanlambda " << tanlambda_seed << " 1/pT " << curvature_seed/(B*0.3e-2) << " p: " << (1/TMath::Cos(TMath::ATan(tanlambda_seed)))*(B*0.3e-2)/curvature_seed<<std::endl;
      }
      

      //std::cout<<"kalman: |p| "<<(1/TMath::Cos(TMath::ATan(tanlambda_seed)))*(0.5*0.3e-2)/curvature_seed<<" dxyz: not given "<< " 1/pT " << curvature_seed/(0.5*0.3e-2) << std::endl;

      //std::cout<<i<<std::endl;

      // Kalman fitter variables

      double zpos = xyz_seed.Z();

      TMatrixD P(5,5);  // covariance matrix of parameters
      // fill in initial guesses -- generous uncertainties on first value.
      P.Zero();
      if (!Smear)
      {
        P[0][0] = TMath::Sq(0.0000001);   // initial position uncertainties -- y 0.1
        P[1][1] = TMath::Sq(0.0000001);   // and x 0.1
        P[2][2] = TMath::Sq(0.0032);  // sinphi uncertainty 3
        P[3][3] = TMath::Sq(0.000034);  // tanlambda uncertainty 3
        P[4][4] = TMath::Sq(0.073);  // q/pT uncertainty 75
      }
      else if (xy_smear == 0.1)
      {
        P[0][0] = TMath::Sq(0.15);   // initial position uncertainties -- y 1
        P[1][1] = TMath::Sq(0.15);   // and x 1
        P[2][2] = TMath::Sq(0.0033);  // sinphi uncertainty 0.5
        P[3][3] = TMath::Sq(0.001);  // tanlambda uncertainty 0.5
        P[4][4] = TMath::Sq(0.073);  // q/pT uncertainty 0.5
      }
      else if (xy_smear == 0.5)
      {
        P[0][0] = TMath::Sq(0.5);   // initial position uncertainties -- y 1
        P[1][1] = TMath::Sq(0.5);   // and x 1
        P[2][2] = TMath::Sq(0.008);  // sinphi uncertainty 0.5
        P[3][3] = TMath::Sq(0.004);  // tanlambda uncertainty 0.5
        P[4][4] = TMath::Sq(0.078);  // q/pT uncertainty 0.5
      }
      if(Seedtype=="alice")
      {
        P = P_seed;
      }
      if(dir<0 &&  !Fixed_Cov) P=Pt.at(Pt.size()-1);

      if (fPrintLevel > 1)
            {
              std::cout << "P0 guess Matrix: " << std::endl;
              P.Print();
            }

      TMatrixD PPred(5,5);
      PPred.Zero();

      // per-step additions to the covariance matrix 
      //need to check this
      TMatrixD Q(5,5);
      Q.Zero();
      //Q[4][4] = 0.00001;     // allow for some curvature uncertainty between points
      //Q[2][2] = 1e-09;      // phi
      //Q[3][3] = 0.0001;   // lambda

      // Noise covariance on the measured points.
      // 16 cm2 initially, might reasonably be lowered to typicalResidual near line 552-67
      TMatrixD R(2,2);
      R.Zero();
      R[0][0] = Ry;  // in cm^2 usually 4 1.1921e-07
      R[1][1] = Rx;//TMath::Sq(0.01); //TMath::Sq(0);  // in cm^2 usually 4
      R[1][0] = Ryx;
      R[0][1] = Ryx;
      //Rt.push_back(R);
      // add the TPCClusters and update the track parameters and uncertainties.  Put in additional terms to keep uncertainties from shrinking when
      // scattering and energy loss can change the track parameters along the way.

      // F = partial(updatefunc)/partial(parvec).  Update functions are in the comments below.

      
      
      TVectorD parvec(5);
      parvec[0] = xyz_seed.Y();
      parvec[1] = xyz_seed.X();
      parvec[2] = sinphi_seed;
      parvec[3] = tanlambda_seed;
      parvec[4] = curvature_seed/(B*0.3e-2);

      TVectorD parvecprev = parvec;

      
      //std::cout << " Parvec: y " << parvec[0] << " z " << parvec[1] << " c " << parvec[2] << " phi " << parvec[3] << " lambda " << parvec[4] << std::endl;
      TVectorD predstep(5);
      predstep.Zero();
      
      
      

      TMatrixD H(2,5);   // partial(obs)/partial(params)
      H.Zero();
      H[0][0] = 1;  // y
      H[1][1] = 1;  // x
      TMatrixD HT(5,2);

      TVectorD zv(2);
      TVectorD ytilde(2);
      TVectorD hx(2);
      TMatrixD S(2,2);
      TMatrixD K(5,2);

      TMatrixD I(5,5);
      I.Zero();
      for (int i=0;i<5;++i) I[i][i] = 1;


      
      //Create a parametrized helix as a substitute for the measurement
      //t1s.GetEntry(0);
      double xh = xyz_plane.at(0).X(); // -23.3112;  ///need to check what this is
      double yh = xyz_plane.at(0).Y(); // -337.045;
      double zh = xyz_plane.at(0).Z(); // v1634.73;      
      
      size_t planecounterk = 10;
      std::vector<double> dEvecreck;
      std::vector<double> dxvecreck;

      for (size_t iTPCCluster=1; iTPCCluster<xyz_plane.size(); ++iTPCCluster)
        {

          
          
          
          //std::cout<<std::endl<<"Real parameters: "<<xyz_plane.at(iTPCCluster).X()<<" "<<xyz_plane.at(iTPCCluster).Y()<<" "<<xyz_plane.at(iTPCCluster).Z()<<" "<<" "<<phi_plane.at(iTPCCluster)<<std::endl<<std::endl;

          xh=xyz_plane.at(iTPCCluster).X(); 
          yh=xyz_plane.at(iTPCCluster).Y();
          zh=xyz_plane.at(iTPCCluster).Z();
          //std::cout<<"zpos: "<<zpos<<" zh: "<<zh<<std::endl;
          
          
          //std::cout<<"iTPCCluster "<<iTPCCluster<<std::endl;
          //std::cout << "xh,yh,zh: "<< xh <<" "<<yh<<" "<<zh<<std::endl;
          //std::cout << " Parvec: y " << parvec[0] << " z " << parvec[1] << " c " << parvec[2] << " phi " << parvec[3] << " lambda " << parvec[4] << std::endl;
          
          
          xht.push_back(xh);  //add measured position to ttree
          yht.push_back(yh);
          zht.push_back(zh);


          ///rotate coordinate system

          
          
          if (fPrintLevel > 0)
            {
              std::cout << std::endl;
              std::cout << "Current zpos: "<<zpos<<std::endl;
              std::cout << "Adding a new TPCCluster: x:" << xh << " y: " << yh << " z: " << zh << std::endl;
            }


          
          
          double dz;
          dz = zh - zpos;
          if (dz == 0) dz = 1E-3;
          
          ////Prepare all the parameters for the prediction

          if (fPrintLevel>0) std::cout<<"dz: "<<dz<<std::endl;

         

          //std::cout<<"dz: "<<dz<<std::endl;
          Bool_t InPlane;
          Double_t dxyzrecord = 0;
          Double_t dErecord = 0;


           PPred=P;
           predstep=parvec;
          
            for(size_t p=0;p<nplanes;p++)
            {
                if(abs(zh-zpos)<=Plane_thick&&abs(zh-zpos)>=0&&((dz>0 && zpos>=Planes_Z[p*2] && zpos<=Planes_Z[p*2+1]) || (dz<0 && zpos>=Planes_Z_bkw[p*2] && zpos<=Planes_Z_bkw[p*2+1])))
                  {

                    InPlane = 1;
                    /// do the prediction
                    if(planecounterk!=p && !dEvecreck.empty() && !dxvecreck.empty() && dir>0)
                    {
                      if (fPrintLevel>0) std::cout<<"Plane "<<p<<" dxreco: ";
                      if (fPrintLevel>0) for(int j =0; j<dxvecreck.size();j++) std::cout<<dxvecreck.at(j)<<" ";
                      if (fPrintLevel>0) std::cout<<"\n";
                      dEreco.push_back(dEvecreck);
                      dxreco.push_back(dxvecreck);
                      dEvecreck.clear();
                      dxvecreck.clear();
                    }
                    if(planecounterk!=p && dir>0) {planecounterk=p;}

                    
                    Double_t dxyzrec, dErec;
                    double dzmid=zh-zpos;
                    if(dzmid!=0)
                    {
                    zpos+=dzmid;
                    if (fPrintLevel>0) std::cout<<zpos<<std::endl;
                    Bool_t checkstatus = Propagate(dzmid,PPred, PPred, predstep,  predstep, kAlmost0, kAlmost1, fPrintLevel, dir, InPlane, dErec, dxyzrec, Energy_loss, CorrTime, MS, xx0);
                      if (!checkstatus)
                      { 
                        status=checkstatus;
                        break;
                      }


                    //if(dir>0) std::cout<<dErec<<" "<<dxyzrec<<std::endl;
                    if(dir>0) dEvecreck.push_back(dErec);
                    if(dir>0) dxvecreck.push_back(dxyzrec);

                    dxyzrecord = dxyzrec;
                    dErecord = dErec;
                    break;
                    }
                    else break;
                    //std::cout<<dEvecreck.size()<<std::endl;
                      
                  }
                

                else if((dz>0 && zpos>=Planes_Z[p*2] && zpos<=Planes_Z[p*2+1]) || (dz<0 && zpos>=Planes_Z_bkw[p*2] && zpos<=Planes_Z_bkw[p*2+1]))
                  {
                    if(planecounterk!=p && !dEvecreck.empty() && !dxvecreck.empty()  && dir>0)
                    {
                      
                      dEreco.push_back(dEvecreck);
                      dxreco.push_back(dxvecreck);
                      dEvecreck.clear();
                      dxvecreck.clear();
                    }
                    if(planecounterk!=p  && dir>0) {planecounterk=p;}

                    double dzmid;
                    Double_t dxyzrec, dErec;
                    Double_t dxyzrectot=0;
                    Double_t dErectot=0;
                    
                    if (dz>0) dzmid=Planes_Z[p*2+1]-zpos;
                    else 
                    {
                    dzmid=Planes_Z_bkw[p*2]-zpos;
                    //std::cout<<"hell yeah"<<std::endl;
                    }

                    Bool_t checkstatus=1;

                    if(dzmid!=0)
                    {
                      zpos+=dzmid;
                      
                      if (fPrintLevel>0) std::cout<<zpos<<std::endl;
                      //if(dir<0)std::cout<<"dz from point to plane: "<< dzmid <<std::endl;
                      InPlane=1;
                      
                      checkstatus = Propagate(dzmid,PPred, PPred, predstep,  predstep, kAlmost0, kAlmost1, fPrintLevel, dir, InPlane, dErec, dxyzrec, Energy_loss, CorrTime, MS, xx0);
                      dErectot+=dErec;
                      dxyzrectot+=dxyzrec;
                      if (!checkstatus)
                      { 
                        status=checkstatus;
                        break;
                      }
                    }                    
                    //if(dir>0)std::cout<<"dxyzmid: "<<dxyzrec;
                    

                    InPlane = false;
                    if (dz>0) dzmid=Planes_Z[p*2+2]-Planes_Z[p*2+1];
                    else 
                    {
                    dzmid=Planes_Z_bkw[p*2+3]-Planes_Z_bkw[p*2];
                    //std::cout<<"hell yeah 2"<<std::endl;
                    }
                    zpos+=dzmid;
                    if (fPrintLevel>0) std::cout<<zpos<<std::endl;
                    //TVectorD predmidstep = predstep;
                    //TMatrixD Predmidstep = PPred;
                    //if(dir<0) std::cout<<"dz between planes: "<< dzmid <<std::endl;
                    checkstatus = Propagate(dzmid,PPred, PPred, predstep,  predstep, kAlmost0, kAlmost1, fPrintLevel, dir, InPlane, dErec, dxyzrec, Energy_loss, CorrTime, MS, xx0);
                    if (!checkstatus)
                    { 
                      status=checkstatus;
                      break;
                    }

                    continue;
                  }
            }

          

          PPred=PPred+Q;
          
          predstept.push_back(predstep);
          PPredt.push_back(PPred);                           // update tree covariance


          
          if (fPrintLevel > 1)
            {
              std::cout << "PPred Matrix: " << std::endl;
              PPred.Print();
            }

          ytilde[0] = yh - predstep[0];
          ytilde[1] = xh - predstep[1];
          

          
          if (fPrintLevel > 1)
            {
              std::cout << "ytilde (residuals): " << std::endl;
              ytilde.Print();
            }
          if (fPrintLevel > 1)
            {
              std::cout << "H Matrix: " << std::endl;
              H.Print();
            }

          HT.Transpose(H);
          S = H*PPred*HT + R;

          if (fPrintLevel > 1)
            {
              std::cout << "S Matrix: " << std::endl;
              S.Print();
            }

          S.Invert();
          if (fPrintLevel > 1)
            {
              std::cout << "Inverted S Matrix: " << std::endl;
              S.Print();
            }

          K = PPred*HT*S;
          if (fPrintLevel > 1)
            {
              std::cout << "K Matrix: " << std::endl;
              K.Print();
            }
          

          parvec = predstep + K*ytilde;
          
          P = (I-K*H)*PPred;

          // Correcting for energy loss after a posteriori
          if(InPlane && Energy_loss && CorrTime == "after") 
          {
          //std::cout<<"deltaxyz Kalman:"<<deltaxyz<<std::endl;
          double p = sqrt(pow(parvec[3]/parvec[4],2)+pow(1/parvec[4],2));
          Bool_t checkstatus= CorrectForMeanMaterial(-dxyzrecord*rho,muon_mass,0.005,p,(p/muon_mass),rho,X0,X1,Ipar,ZA,parvec,PPred,dErecord,dir, MS,dxyzrecord/xx0);
          //std::cout<<"I'm correcting"<<std::endl;
          if (!checkstatus)
                    { 
                      status=checkstatus;
                      break;
                    }
          
          }

          parvect.push_back(parvec);               //add best estimate to ttree
          Pt.push_back(P);                           //add covariance best estimate to tree
       
                               //add xpos to tree

          

          parvecprev = parvec;
          //zpos = zpos + dz;
          
          zpost.push_back(zpos);

          if (fPrintLevel > 0)
            {
              std::cout << "z: " << zpos << " dz: " << dz <<  std::endl;
              std::cout << " Parvec:   y " << parvec[0] << " x " << parvec[1] << " sinphi " << parvec[2] << " tanlambda " << parvec[3] << " 1/pT " << parvec[4] << " p: " << (1/TMath::Cos(TMath::ATan(parvec[3])))/parvec[4] <<std::endl;
            }
          
          //std::cout << " Updated xpos: " << xpos << " " << dx << std::endl;
          
          
        //std::cout<<"dir="<<dir<<" invpTplane="<<invpTplane.at(iTPCCluster)<<" predstep[4]="<<predstep[4]<<" parvec[4]"<<parvec[4]<<" invpTseed="<<curvature_seed/(0.5*0.299792458e-2)<<" p00="<<P[1][1]<<std::endl;
        //std::cout<<"dir="<<dir<<" ytilde: "<<ytilde[0]<<" "<<ytilde[1]<<std::endl;
        }
        if(dir>0)dEreco.push_back(dEvecreck);
        if(dir>0)dxreco.push_back(dxvecreck);
      
    }


#endif