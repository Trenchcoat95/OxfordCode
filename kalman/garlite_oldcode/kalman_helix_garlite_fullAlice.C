



// ROOT includes

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







//--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

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
                   std::vector<double> sinphi_plane,
                   double Ry,
                   double Rx,
                   double Rxy,
                   XYZVector  xyz_seed,
                   double  tanlambda_seed,
                   double curvature_seed,
                   double sinphi_seed,
                   double dir
                   

                                )
    {

      

      
      double fPrintLevel=1;
      const Double_t kAlmost1=1. - Double_t(FLT_EPSILON);
      const Double_t kAlmost0=Double_t(FLT_MIN);

      
      
      
      
      
      if(fPrintLevel>0)
      {
        std::cout << " Seed: y " << xyz_seed.Y() << " x " << xyz_seed.X() <<"z:"<<xyz_seed.Z()<< " sinphi " << sinphi_seed << " tanlambda " << tanlambda_seed << " 1/pT " << curvature_seed/(0.5*0.299792458e-2) << " p: " << (1/TMath::Cos(TMath::ATan(tanlambda_seed)))*(0.5*0.299792458e-2)/curvature_seed<<std::endl;
      }
      

      

      //std::cout<<i<<std::endl;

      // Kalman fitter variables

      double zpos = xyz_seed.Z();

      TMatrixD P(5,5);  // covariance matrix of parameters
      // fill in initial guesses -- generous uncertainties on first value.
      P.Zero();
      P[0][0] = TMath::Sq(1);   // initial position uncertainties -- y
      P[1][1] = TMath::Sq(1);   // and z
      P[2][2] = TMath::Sq(.5);  // curvature of zero gets us to infinite momentum, and curvature of 2 is curled up tighter than the pads
      P[3][3] = TMath::Sq(.5);  // phi uncertainty
      P[4][4] = TMath::Sq(.5);  // lambda uncertainty
      if(dir<0) P=Pt.at(Pt.size()-1);
      //TMatrixD PPred(5,5);
      //PPred.Zero();

      // per-step additions to the covariance matrix 
      //need to check this
      /*
      TMatrixD Q(5,5);
      Q.Zero();
      Q[4][4] = 2;     // allow for some curvature uncertainty between points
      Q[2][2] = 1e-09;      // phi
      Q[3][3] = 0.0001;   // lambda
      */

      // Noise covariance on the measured points.
      // 16 cm2 initially, might reasonably be lowered to typicalResidual near line 552-67
      /*
      TMatrixD R(2,2);
      R.Zero();
      R[0][0] = Ry;  // in cm^2 usually 4 1.1921e-07
      R[1][1] = Rz;//TMath::Sq(0.01); //TMath::Sq(0);  // in cm^2 usually 4
      R[1][0] = Ryz;
      R[0][1] = Ryz;
      */
      //Rt.push_back(R);
      // add the TPCClusters and update the track parameters and uncertainties.  Put in additional terms to keep uncertainties from shrinking when
      // scattering and energy loss can change the track parameters along the way.

      // F = partial(updatefunc)/partial(parvec).  Update functions are in the comments below.

      
      
      TVectorD parvec(5);
      parvec[0] = xyz_seed.Y();
      parvec[1] = xyz_seed.X();
      parvec[2] = sinphi_seed;
      parvec[3] = tanlambda_seed;
      parvec[4] = curvature_seed/(0.5*0.299792458e-2);

      
      //std::cout << " Parvec: y " << parvec[0] << " z " << parvec[1] << " c " << parvec[2] << " phi " << parvec[3] << " lambda " << parvec[4] << std::endl;
      //TVectorD predstep(5);
      //predstep.Zero();
      
      
      
      /*
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
      */

      TVectorD ytilde(2);
      
      //Create a parametrized helix as a substitute for the measurement
      //t1s.GetEntry(0);
      double xh = xyz_plane.at(0).X(); // -23.3112;  ///need to check what this is
      double yh = xyz_plane.at(0).Y(); // -337.045;
      double zh = xyz_plane.at(0).Z(); // v1634.73;      
      

      for (size_t iTPCCluster=1; iTPCCluster<xyz_plane.size(); ++iTPCCluster)
        {
          
          
          //std::cout<<std::endl<<"Real parameters: "<<xyz_plane.at(iTPCCluster).X()<<" "<<xyz_plane.at(iTPCCluster).Y()<<" "<<xyz_plane.at(iTPCCluster).Z()<<" "<<" "<<phi_plane.at(iTPCCluster)<<std::endl<<std::endl;

          xh=xyz_plane.at(iTPCCluster).X(); 
          yh=xyz_plane.at(iTPCCluster).Y();
          zh=xyz_plane.at(iTPCCluster).Z();

          
          
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
              std::cout << "Adding a new TPCCluster: x:" << xh << " y: " << yh << " z: " << zh << " sinphi: " << sinphi_plane.at(iTPCCluster) << std::endl;
            }

          // for readability
          double curvature = (0.5*0.299792458e-2)*parvec[4];
          double sinphi = parvec[2];
          double tanlambda = parvec[3];

         

          // update prediction to the plane containing x.  Maybe we need to find
          // the closest point on the helix to the TPCCluster we are adding,
          // and not necessarily force it to be at this x

          // y = yold + slope*dx*Sin(phi).   F[0][i] = dy/dtrackpar[i], where f is the update function slope*dx*Sin(phi)

          
          
          double dz;
          dz = zh - zpos;
          if (dz == 0) dz = 1E-3;
          
          ////Prepare all the parameters for the prediction

          if (fPrintLevel>0) std::cout<<"dz: "<<dz<<std::endl;

          Double_t z2r = curvature*dz;
          Double_t f1=sinphi, f2=f1 + z2r;
          if (TMath::Abs(f1) >= kAlmost1) break;
          if (TMath::Abs(f2) >= kAlmost1) break;
          if (TMath::Abs(tanlambda)< kAlmost0) break;

          
          
          Double_t r1=TMath::Sqrt((1.-f1)*(1.+f1)), r2=TMath::Sqrt((1.-f2)*(1.+f2));
          if (TMath::Abs(r1)<kAlmost0)  break;
          if (TMath::Abs(r2)<kAlmost0)  break;
 
          double dy2dz = (f1+f2)/(r1+r2);
          double rot = TMath::ASin(r1*f2 - r2*f1); 
            if (f1*f1+f2*f2>1 && f1*f2<0) {          // special cases of large rotations or large abs angles
              if (f2>0) rot =  TMath::Pi() - rot;    //
              else      rot = -TMath::Pi() - rot;
            }

          // predicted step
          if (fPrintLevel > 1)
            {
              std::cout << "P Matrix: " << std::endl;
              P.Print();
            }
          

          

          if(fPrintLevel>0) std::cout<<"z2r: "<<z2r<<std::endl;
          
          parvec[0] += dz*dy2dz;  // update y
          parvec[2] += z2r;        // update sinphi
          parvec[1] +=tanlambda/curvature*rot; //update x
          predstept.push_back(parvec);                       // update tree values
          
          if (fPrintLevel >0 )
            {
              std::cout << " Predstep: y " << parvec[0] << " x " << parvec[1] << " sinphi " << parvec[2] << " tanlambda " << parvec[3] << " 1/pT " << parvec[4] << std::endl;
            }
          
          
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
          
          
          //F*C*Ft = C + (b + bt + a)
          P[0][0] += b00 + b00 + a00;
          P[1][0] += b10 + b01 + a01; 
          P[2][0] += b20 + b02 + a02;
          P[3][0] += b30;
          P[4][0] += b40;
          P[1][1] += b11 + b11 + a11;
          P[2][1] += b21 + b12 + a12;
          P[3][1] += b31; 
          P[4][1] += b41;
          P[2][2] += b22 + b22 + a22;
          P[3][1] += b32;
          P[4][1] += b42;

          PPredt.push_back(P);                           // update tree pred covariance

          if (fPrintLevel > 1)
            {
              std::cout << "PPred Matrix: " << std::endl;
              P.Print();
            }


          ///////Update step

          Double_t r00=Ry, r01=Rxy, r11=Rx;
          r00+=P[0][0]; r01+=P[1][0]; r11+=P[1][1];
          Double_t det=r00*r11 - r01*r01;

          if (TMath::Abs(det) < kAlmost0) break;

          Double_t tmp=r00; r00=r11/det; r11=tmp/det; r01=-r01/det;

          ///////Kalman Gain
 
          Double_t k00=P[0][0]*r00+P[1][0]*r01, k01=P[0][0]*r01+P[1][0]*r11;
          Double_t k10=P[1][0]*r00+P[1][1]*r01, k11=P[1][0]*r01+P[1][1]*r11;
          Double_t k20=P[2][0]*r00+P[2][1]*r01, k21=P[2][0]*r01+P[2][1]*r11;
          Double_t k30=P[3][0]*r00+P[3][1]*r01, k31=P[3][0]*r01+P[3][1]*r11;
          Double_t k40=P[4][0]*r00+P[4][1]*r01, k41=P[4][0]*r01+P[4][1]*r11;
          

          /////Residuals
          
          ytilde[0] = yh - parvec[0];
          ytilde[1] = xh - parvec[1]; 
          
          if (fPrintLevel > 0)
            {
              std::cout << "ytilde (residuals): " << std::endl;
              ytilde.Print();
            }

          ////Update parameter vector
          
          Double_t sf=parvec[2] + k20*ytilde[0] + k21*ytilde[1];
          if (TMath::Abs(sf) > kAlmost1) break; 

          parvec[0] +=  k00*ytilde[0] + k01*ytilde[1];
          parvec[1] +=  k10*ytilde[0] + k11*ytilde[1];
          parvec[2]  =  sf;
          parvec[3] +=  k30*ytilde[0] + k31*ytilde[1];
          parvec[4] +=  k40*ytilde[0] + k41*ytilde[1];
                    
          parvect.push_back(parvec);               //add best estimate to ttree


          /////Update covariance matrix

          Double_t c01=P[1][0], c02=P[2][0], c03=P[3][0], c04=P[4][0];
          Double_t c12=P[2][1], c13=P[3][1], c14=P[4][1];

          P[0][0]-=k00*P[0][0]+k01*P[1][0]; P[1][0]-=k00*c01+k01*P[1][1];
          P[2][0]-=k00*c02+k01*c12;   P[3][0]-=k00*c03+k01*c13;
          P[4][0]-=k00*c04+k01*c14; 

          P[1][1]-=k10*c01+k11*P[1][1];
          P[2][1]-=k10*c02+k11*c12;   P[3][1]-=k10*c03+k11*c13;
          P[4][1]-=k10*c04+k11*c14; 

          P[2][2]-=k20*c02+k21*c12;   P[3][2]-=k20*c03+k21*c13;
          P[4][2]-=k20*c04+k21*c14; 

          P[3][3]-=k30*c03+k31*c13;
          P[4][3]-=k30*c04+k31*c14; 
          
          P[4][4]-=k40*c04+k41*c14; 
          
          Pt.push_back(P);                           //add covariance best estimate to tree
       
          zpos = zpos + dz;
          zpost.push_back(zpos);                     //add xpos to tree

          if (fPrintLevel > 0)
            {
              std::cout << "z: " << zpos << " dz: " << dz <<  std::endl;
              std::cout << " Parvec:   y " << parvec[0] << " x " << parvec[1] << " sinphi " << parvec[2] << " tanlambda " << parvec[3] << " 1/pT " << parvec[4] << "p: " << (1/TMath::Cos(TMath::ATan(parvec[3])))/parvec[4] <<std::endl;
            }
          
          //std::cout << " Updated xpos: " << xpos << " " << dx << std::endl;
          
          
           
          
        }
        
      
    }
    
    void kalman_helix_garlite_fullAlice()

    {
      const Double_t kAlmost1=1. - Double_t(FLT_EPSILON);
      const Double_t kAlmost0=Double_t(FLT_MIN);
      // variables:  z is the independent variable
      // 0: y
      // 1: x
      // 2: sin(phi)
      // 3: tan(lambda)
      // 4: 1/pT (1/(GeV/c))
      ///Prepare a simple tree with just the coordinates

      

      //Right now with perfect helix
      std::cout << std::setprecision(17);
      

      TFile fs("garlitetest.root","recreate");
      TTree t1s("t1s","helix simple tree");
      double  Plane_XY[4] = {-300, +300, -400, 100}; //all in cm
      double  Plane_XY_fid[4] = {-200, +200, -300, 0};
      double  Planes_Z[10] =   {1244,1248,1334,1338,1484,1488,1634,1638,1724,1728};
      XYZVector GArCenter(0,-150.473,1486); 
      double GAr_r = 349.9;
      double GAr_L = 669.6;
      double Plane_thick = 4;
      double  B=0.5; //in Tesla
      XYZVector Z_axis(0,0,1);
      XYZVector Y_axis(0,1,0);
     
      
      
      ////MC

      std::vector<XYZVector> xyz;
      std::vector<XYZVector> pxyz;
      std::vector<XYZVector> xyz_plane;
      std::vector<XYZVector> xyz_plane_sm;
      std::vector<XYZVector> xyz_plane_mean;
      std::vector<double> sinphi;
      std::vector<double> sinphi_plane;
      std::vector<double> tanlambda;
      double  charge;
      int nPoints;
      std::vector<int> nHits_perPlane;
      std::vector<XYZVector> xyz_firstinPlane;
      std::vector<XYZVector> pxyz_firstinPlane;
      std::vector<double> sinphi_firstinPlane;
      std::vector<double> tanlambda_firstinPlane;
      std::vector<double> curvature_firstinPlane;
      

      ////Kalman

      double Ry = TMath::Sq(0.2);//TMath::Sq(2);//TMath::Sq(0.01);//TMath::Sq(3);
      double Rz = TMath::Sq(0.2);//TMath::Sq(0.01);//TMath::Sq(0.01); //1.1921e-07
      double Ryz = TMath::Sq(0);

      ////FWD
      XYZVector xyz_seed;
      double sinphi_seed,tanlambda_seed,curvature_seed;
      std::vector<double> xht,yht,zht,zpost;
      std::vector<TVectorD> parvect; //(5)
      std::vector<TVectorD> predstept; //(5)
      std::vector<TMatrixD> Pt;//(5,5);
      std::vector<TMatrixD> PPredt;//(5,5);
      std::vector<TMatrixD> Rt;//(2,2);

      ///BKW
      XYZVector xyz_seed_bkw;
      double sinphi_seed_bkw,tanlambda_seed_bkw,curvature_seed_bkw;
      std::vector<double> xht_bkw,yht_bkw,zht_bkw,zpost_bkw;
      std::vector<TVectorD> parvect_bkw; //(5)
      std::vector<TVectorD> predstept_bkw; //(5)
      std::vector<TMatrixD> Pt_bkw;//(5,5);
      std::vector<TMatrixD> PPredt_bkw;//(5,5);
      std::vector<TMatrixD> Rt_bkw;//(2,2);


      ///////MC
      t1s.Branch("xyz",&xyz);
      t1s.Branch("pxyz",&pxyz);
      t1s.Branch("sinphi",&sinphi);
      t1s.Branch("sinphi_plane",&sinphi_plane);
      t1s.Branch("tanlambda",&tanlambda);
      t1s.Branch("xyz_plane",&xyz_plane);
      t1s.Branch("xyz_plane_sm",&xyz_plane_sm);
      t1s.Branch("xyz_plane_mean",&xyz_plane_mean);
      t1s.Branch("charge",&charge);
      t1s.Branch("nPoints",&nPoints);
      t1s.Branch("nHits_perPlane",&nHits_perPlane);
      t1s.Branch("xyz_firstinPlane",&xyz_firstinPlane);
      t1s.Branch("sinphi_firstinPlane",&sinphi_firstinPlane);
      t1s.Branch("tanlambda_firstinPlane",&tanlambda_firstinPlane);
      t1s.Branch("pxyz_firstinPlane",&pxyz_firstinPlane);
      t1s.Branch("curvature_firstinPlane",&curvature_firstinPlane);
      

      ////Kalman
      t1s.Branch("xht",&xht);
      t1s.Branch("yht",&yht);
      t1s.Branch("zht",&zht);
      t1s.Branch("zpost",&zpost);
      t1s.Branch("parvect",&parvect);
      t1s.Branch("predstept",&predstept);
      t1s.Branch("Pt",&Pt);
      t1s.Branch("PPredt",&PPredt);
      t1s.Branch("Rt",&Rt);
      t1s.Branch("xyz_seed",&xyz_seed);
      t1s.Branch("sinphi_seed",&sinphi_seed);
      t1s.Branch("tanlambda_seed",&tanlambda_seed);
      t1s.Branch("curvature_seed",&curvature_seed);

      t1s.Branch("xht_bkw",&xht_bkw);
      t1s.Branch("yht_bkw",&yht_bkw);
      t1s.Branch("zht_bkw",&zht_bkw);
      t1s.Branch("zpost_bkw",&zpost_bkw);
      t1s.Branch("parvect_bkw",&parvect_bkw);
      t1s.Branch("predstept_bkw",&predstept_bkw);
      t1s.Branch("Pt_bkw",&Pt_bkw);
      t1s.Branch("PPredt_bkw",&PPredt_bkw);
      t1s.Branch("Rt_bkw",&Rt_bkw);
      t1s.Branch("xyz_seed_bkw",&xyz_seed_bkw);
      t1s.Branch("sinphi_seed_bkw",&sinphi_seed_bkw);
      t1s.Branch("tanlambda_seed_bkw",&tanlambda_seed_bkw);
      t1s.Branch("curvature_seed_bkw",&curvature_seed_bkw);
      
     
      
      size_t nevents=1;
      
      for(size_t i=0;i<nevents;i++)
      {
      XYZVector xyztemp(gRandom->Rndm()*(Plane_XY_fid[1]-Plane_XY_fid[0])+Plane_XY_fid[0],gRandom->Rndm()*(Plane_XY_fid[3]-Plane_XY_fid[2])+Plane_XY_fid[2],gRandom->Rndm()*10+(Planes_Z[0]-10));
      xyz.push_back(xyztemp);
      
      charge = (gRandom->Rndm()<0.5) ? -1 : 1;
      
      //std::cout<<i<<std::endl;
      
      
      double ptotal=0.1+gRandom->Rndm()*4; //in GeV/c2
      double px=-0.5+gRandom->Rndm();
      double py=-0.5+gRandom->Rndm();
      double pz=1.0+gRandom->Rndm();
      double norm=ptotal/sqrt(px*px+py*py+pz*pz);
      
      XYZVector pxyztemp(px*norm,py*norm,pz*norm);
      pxyz.push_back(pxyztemp);

      
      XYZVector pprojyz(0,pxyztemp.Y(),pxyztemp.Z()); ///pT
      //XYZVector pprojx(pxyztemp.X(),0,0);
      
      double sinphitemp = 0;
      if(pprojyz.Y()>0) sinphitemp=TMath::Sin(TMath::ACos(pprojyz.Unit().Dot(Z_axis.Unit())));
      else sinphitemp=-TMath::Sin(TMath::ACos(pprojyz.Unit().Dot(Z_axis.Unit())));
      //std::cout<<pxyztemp<<" "<<sinphitemp<<std::endl;
      double curvaturetemp =charge*((0.299792458e-2)*B)/(sqrt(pxyztemp.Y()*pxyztemp.Y()+pxyztemp.Z()*pxyztemp.Z())); //in cm^-1
      double tanlambdatemp =0;
      if(pxyztemp.X()>0) tanlambdatemp =TMath::Tan(TMath::ACos(pxyztemp.Unit().Dot(pprojyz.Unit())));
      else tanlambdatemp =-TMath::Tan(TMath::ACos(pxyztemp.Unit().Dot(pprojyz.Unit())));

      double pyzmodule=sqrt(pxyztemp.Y()*pxyztemp.Y()+pxyztemp.Z()*pxyztemp.Z());
     
      //if(pxyztemp.X()!= (sqrt(pxyztemp.Y()*pxyztemp.Y()+pxyztemp.Z()*pxyztemp.Z())*tanlambdatemp))std::cout<<"px= "<<pxyztemp.X()<<" "<<(sqrt(pxyztemp.Y()*pxyztemp.Y()+pxyztemp.Z()*pxyztemp.Z())*tanlambdatemp)<<std::endl;
      //std::cout<<"py= "<<pxyztemp.Y()<<" pz: "<<pxyztemp.Z()<<std::endl;
      //std::cout<<"ptotal: "<<ptotal<<" "<< (sqrt(pxyztemp.Y()*pxyztemp.Y()+pxyztemp.Z()*pxyztemp.Z()+pxyztemp.X()*pxyztemp.X()))<<std::endl;
      //std::cout<<"CosLambda: "<<pxyztemp.Unit().Dot(pprojyz.Unit())<<" "<<pprojyz.Dot(pxyztemp)/(sqrt(pxyztemp.Y()*pxyztemp.Y()+pxyztemp.Z()*pxyztemp.Z()+pxyztemp.X()*pxyztemp.X())*sqrt(pprojyz.Y()*pprojyz.Y()+pprojyz.Z()*pprojyz.Z()+pprojyz.X()*pprojyz.X()))<<std::endl;
      
      //TF2 *f2 = new TF2("f2","exp(-0.5*((x)/3)**2)*exp(-0.5*((y)/3)**2)",-20,20,-20,20);
      
      //TCanvas *mccanvaspT = new TCanvas("mccanvaspT","",1000,800);
      //f2->Draw();

      Bool_t InGAr=true;
      Bool_t InPlane[5]={0,0,0,0,0};
      double dz=0.50000000000000;

      //std::cout << " Initial coordinates: y " << xyztemp.Y() << " x " << xyztemp.X() << " sinphi " << sinphitemp << " tanlambda " << tanlambdatemp << " 1/pT " << charge/(sqrt(pprojyz.Y()*pprojyz.Y()+pprojyz.Z()*pprojyz.Z())) << std::endl;
     

      while(InGAr)
        {
            
             //Randomized x
            //double dx=0.04;

            Double_t z2r = curvaturetemp*dz;
            Double_t f1=sinphitemp, f2=f1 + z2r;
            
            if (TMath::Abs(f1) >= kAlmost1) break;
            if (TMath::Abs(f2) >= kAlmost1) break;
            if (TMath::Abs(tanlambdatemp)< kAlmost0) break;

            Double_t r1=TMath::Sqrt((1.-f1)*(1.+f1)), r2=TMath::Sqrt((1.-f2)*(1.+f2));
            
            if (TMath::Abs(r1)<kAlmost0)  break;
            if (TMath::Abs(r2)<kAlmost0)  break;

            double dy2dz = (f1+f2)/(r1+r2);
            
            xyztemp.SetY(xyztemp.Y()+dz*dy2dz);
            xyztemp.SetZ(xyztemp.Z()+dz);
            sinphitemp+=z2r;
            pxyztemp.SetY(pyzmodule*sinphitemp);
            pxyztemp.SetZ(pyzmodule*sqrt(1-sinphitemp*sinphitemp));
            //std::cout<<sqrt(pxyztemp.Y()*pxyztemp.Y()+pxyztemp.Z()*pxyztemp.Z())<<" "<<sinphitemp<<std::endl;



            double rot = TMath::ASin(r1*f2 - r2*f1); 
            if (f1*f1+f2*f2>1 && f1*f2<0) {          // special cases of large rotations or large abs angles
              if (f2>0) rot =  TMath::Pi() - rot;    //
              else      rot = -TMath::Pi() - rot;
            }
             

            xyztemp.SetX(xyztemp.X()+tanlambdatemp/curvaturetemp*rot);

            xyz.push_back(xyztemp);
            pxyz.push_back(pxyztemp);
            sinphi.push_back(sinphitemp);
            tanlambda.push_back(tanlambdatemp);

            //std::cout<<"XYZphi MC: "<< xyztemp.X()<< " " << xyztemp.Y()<< " " << xyztemp.Z()<< " " <<phitemp<<std::endl;
            //std::cout<<"XYZ MC:"<<x<<std::endl;
          
            
            
            for(size_t p=0;p<5;p++)
            {
              
              if(xyztemp.X()>Plane_XY[0] && xyztemp.X()<Plane_XY[1] && xyztemp.Y()>Plane_XY[2] && xyztemp.Y()<Plane_XY[3])
              {
                if(xyztemp.Z()>Planes_Z[p*2] && xyztemp.Z()<Planes_Z[p*2+1] && InPlane[p]==0)
                  {
                    InPlane[p]=1;
                    int hitsinplane=1+gRandom->Gaus(10000,1);
                    //std::cout<<hitsinplane<<std::endl;
                    nHits_perPlane.push_back(hitsinplane);
                    xyz_firstinPlane.push_back(xyztemp);
                    pxyz_firstinPlane.push_back(pxyztemp);
                    tanlambda_firstinPlane.push_back(tanlambdatemp);
                    sinphi_firstinPlane.push_back(sinphitemp);
                    curvature_firstinPlane.push_back(curvaturetemp);

                    XYZVector xyz_plane_temp=xyztemp;
                    XYZVector xyz_plane_temp_sm=xyztemp;
                    XYZVector pxyz_plane_temp=pxyztemp;
               
                    std::vector<XYZVector> xyz_plane_m;
                    double sinphi_plane_temp=sinphitemp;
                    double tanlamda_plane_temp=tanlambdatemp;
                    double curvature_plane_temp=curvaturetemp;

                    //std::cout<<"XYZ MC: "<< xyztemp.X()<< " " << xyztemp.Y()<< " " << xyztemp.Z()<< std::endl;
                    for(int i=0;i<hitsinplane;i++)
                      {
                        //XYZVector hit(gRandom->Gaus(xyztemp.X(),1),gRandom->Gaus(xyztemp.Y(),1),Planes_Z[p*2]+Plane_thick*gRandom->Rndm());
                        double dzp=Plane_thick/(hitsinplane+1);

                        Double_t z2rp = curvature_plane_temp*dzp;
                        Double_t f1p=sinphi_plane_temp, f2p=f1p + z2rp;
                        Double_t r1p=TMath::Sqrt((1.-f1p)*(1.+f1p)), r2p=TMath::Sqrt((1.-f2p)*(1.+f2p));
                        double dy2dzp = (f1p+f2p)/(r1p+r2p);
                        
                        
                        xyz_plane_temp.SetY(xyz_plane_temp.Y()+dzp*dy2dzp);
                        xyz_plane_temp.SetZ(xyz_plane_temp.Z()+dzp);
                        sinphi_plane_temp+=z2rp;
                        pxyz_plane_temp.SetY(sqrt(pxyz_plane_temp.Y()*pxyz_plane_temp.Y()+pxyz_plane_temp.Z()*pxyz_plane_temp.Z())*sinphi_plane_temp);
                        pxyz_plane_temp.SetZ(sqrt(pxyz_plane_temp.Y()*pxyz_plane_temp.Y()+pxyz_plane_temp.Z()*pxyz_plane_temp.Z())*sqrt(1-sinphitemp*sinphitemp));


                        double rotp = TMath::ASin(r1p*f2p - r2p*f1p); 
                        if (f1p*f1p+f2p*f2p>1 && f1p*f2p<0) {          // special cases of large rotations or large abs angles
                          if (f2p>0) rotp =  TMath::Pi() - rotp;    //
                          else      rotp = -TMath::Pi() - rotp;
                        }
                        

                        xyz_plane_temp.SetX(xyz_plane_temp.X()+tanlamda_plane_temp/curvature_plane_temp*rotp);

                        //std::cout<<"XYZ MC: "<< xyz_plane_temp.X()<< " " << xyz_plane_temp.Y()<< " " << xyz_plane_temp.Z()<<" "<<phi_plane_temp<< std::endl;

                        //XYZVector hit(gRandom->Gaus(xyz_plane_temp.X(),0.5),gRandom->Gaus(xyz_plane_temp.Y(),0.5),gRandom->Gaus(xyz_plane_temp.Z(),0.5));
                        
                        xyz_plane.push_back( xyz_plane_temp);
                        sinphi_plane.push_back(sinphi_plane_temp);
                        xyz_plane_m.push_back( xyz_plane_temp);
                        
                        //std::cout<<"XYZ hit: "<< hit.X()<< " " << hit.Y()<< " " << hit.Z()<< std::endl;
                      }
                    //xyz_plane.push_back(xyztemp);
                    std::sort(xyz_plane.begin(),xyz_plane.end(),
                    [](const XYZVector &a, const XYZVector &b)->bool
                    { return a.Z() < b.Z(); } );
                    //std::cout<< " "<<std::endl;
                    XYZVector mean(0,0,0);
                    for(int k=0; k<xyz_plane_m.size(); k++) mean+=xyz_plane_m.at(k);
                    xyz_plane_mean.push_back(mean/xyz_plane_m.size());
                    
                  }
              }
            }
            

            double r=sqrt(TMath::Power((xyztemp.Y()-GArCenter.Y()),2)+TMath::Power((xyztemp.Z()-GArCenter.Z()),2));

            if(xyztemp.X()>0.5*GAr_L || xyztemp.X()<-0.5*GAr_L || r>GAr_r) InGAr=false;
            

        }

      nPoints=xyz.size();
      //std::cout<<"nPoints: "<<nPoints<<std::endl;
      //nHits= xyz_plane.size();

      //std::cout<<"Recorded hit z: ";
      //for(int i=0;i<xyz_plane->size();i++) std::cout<<xyz_plane->at(i).Z()<<" ";
      //std::cout<<""<<std::endl;
      /*
      if(xyz_plane.size()>0)
      {
      for(int k=0; k<xyz_plane.size(); k++) std::cout<<xyz_plane.at(k).X()<<" "<<xyz_plane.at(k).Y()<<" "<<xyz_plane.at(k).Z()<<" "<<std::endl;
      std::cout<<" "<<std::endl;
      }
      */
      if(xyz_plane.size()>0 && nHits_perPlane.size()>=3) 
      {
        
        //////FWD

        ////Smear the plane points

        
        
        xyz_seed.SetX(gRandom->Gaus(xyz_firstinPlane.at(0).X(),2));
        xyz_seed.SetY(gRandom->Gaus(xyz_firstinPlane.at(0).Y(),2));
        xyz_seed.SetZ(gRandom->Gaus(xyz_firstinPlane.at(0).Z(),2));
        curvature_seed=gRandom->Gaus(curvature_firstinPlane.at(0),curvature_firstinPlane.at(0)/10);
        sinphi_seed=gRandom->Gaus(sinphi_firstinPlane.at(0),sinphi_firstinPlane.at(0)/10);
        tanlambda_seed=gRandom->Gaus(tanlambda_firstinPlane.at(0),tanlambda_firstinPlane.at(0)/10);
        double forward=1.;
        
        /*
        xyz_seed=xyz_firstinPlane.at(0);
        curvature_seed=curvature_firstinPlane.at(0);
        sinphi_seed=sinphi_firstinPlane.at(0);
        tanlambda_seed=tanlambda_firstinPlane.at(0);
        double forward=1.;
        */

        std::cout << " Real Values: y " << xyz_firstinPlane.at(0).Y() << " x " << xyz_firstinPlane.at(0).X() << " sinphi " << sinphi_firstinPlane.at(0) << " tanlambda " << tanlambda_firstinPlane.at(0) << " 1/pT " << curvature_firstinPlane.at(0)/(0.5*0.299792458e-2) << " p: " <<sqrt(pxyz_firstinPlane.at(0).Mag2())<< std::endl;
        KalmanFit(xht,yht,zht,parvect,predstept,Pt,PPredt,Rt,zpost,xyz_plane,sinphi_plane,Ry,Rz,Ryz,xyz_seed,tanlambda_seed,curvature_seed,sinphi_seed,forward);
        
        //////BKW
        //std::cout<<xht.size()<<" "<<xyz_plane.size()<<" "<<nHits_perPlane.size()<<std::endl;
        xyz_seed_bkw.SetX(parvect.at(parvect.size()-1)[1]);
        xyz_seed_bkw.SetY(parvect.at(parvect.size()-1)[0]);
        xyz_seed_bkw.SetZ(zht.at(xht.size()-1));
        curvature_seed_bkw=parvect.at(parvect.size()-1)[4]*(0.5*0.299792458e-2);
        sinphi_seed_bkw=parvect.at(parvect.size()-1)[2];
        tanlambda_seed_bkw=parvect.at(parvect.size()-1)[3];
        std::vector<XYZVector> xyz_plane_bkw=xyz_plane;
        std::reverse(xyz_plane_bkw.begin(),xyz_plane_bkw.end());
        std::vector<double> sinphi_plane_bkw=sinphi_plane;
        std::reverse(sinphi_plane_bkw.begin(),sinphi_plane_bkw.end());
        Pt_bkw=Pt;
        //std::cout<<xyz_seed_bkw.X()<<" "<<xyz_seed_bkw.X()<<" "<<xyz_seed_bkw.X()<<" "<<curvature_seed<<" "<<phi_seed_bkw<<" "<<lambda_seed_bkw<<std::endl;
        std::cout<<" "<<std::endl<<std::endl;
        double backwards=-1.;
        KalmanFit(xht_bkw,yht_bkw,zht_bkw,parvect_bkw,predstept_bkw,Pt_bkw,PPredt_bkw,Rt_bkw,zpost_bkw,xyz_plane_bkw,sinphi_plane_bkw,Ry,Rz,Ryz,xyz_seed_bkw,tanlambda_seed_bkw,curvature_seed_bkw,sinphi_seed_bkw,backwards);
        std::cout<<" "<<std::endl<<std::endl;
      }
      t1s.Fill();

      //std::cout<<"nHits"<<xyz_plane.size()<<std::endl;


      xyz.clear();
      pxyz.clear();
      xyz_plane.clear();
      xht.clear();
      yht.clear();
      zht.clear();
      zpost.clear();
      zpost_bkw.clear();
      parvect.clear();
      predstept.clear();
      Pt.clear();//(5,5);
      PPredt.clear();//(5,5);
      Rt.clear();//(2,2);
      sinphi.clear();
      tanlambda.clear();
      
      xyz_firstinPlane.clear();
      sinphi_firstinPlane.clear();
      tanlambda_firstinPlane.clear();
      curvature_firstinPlane.clear();
      pxyz_firstinPlane.clear();
      nHits_perPlane.clear();
      xyz_plane_mean.clear();
      
    }
    
    t1s.Write();
          
      
  }

    

    