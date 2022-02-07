////////////////////////////////////////////////////////////////////////
// Class:       tpctrackfit2
// Plugin Type: producer (art v3_00_00)
// File:        tpctrackfit2_module.cc
//
// Generated at Tue Feb  5 11:34:54 2019 by Thomas Junk using cetskelgen
// from cetlib version v3_04_00.
////////////////////////////////////////////////////////////////////////



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
    // variables:  x is the independent variable
    // 0: x
    // 1: y
    // 2: curvature
    // 3: phi
    // 4: lambda

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
                   std::vector<double> phi_plane,
                   double Ry,
                   double Rz,
                   double Ryz,
                   XYZVector  xyz_seed,
                   double  lambda_seed,
                   double curvature_seed,
                   double phi_seed,
                   double dir
                   

                                )
    {

      

      
      double fPrintLevel=0;

      
      
      
      
      
      if(fPrintLevel>0)
      {
        std::cout<<std::endl<<"Seed parameters: "<<xyz_seed.X()<<" "<<xyz_seed.Y()<<" "<<xyz_seed.Z()<<" "<<phi_seed<<" "<<lambda_seed<<std::endl<<std::endl;
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
      TMatrixD PPred(5,5);
      PPred.Zero();

      // per-step additions to the covariance matrix 
      //need to check this
      TMatrixD Q(5,5);
      Q.Zero();
      Q[2][2] = 1e-09;     // allow for some curvature uncertainty between points
      Q[3][3] = 1e-09;      // phi
      Q[4][4] = 0.0001;   // lambda

      // Noise covariance on the measured points.
      // 16 cm2 initially, might reasonably be lowered to typicalResidual near line 552-67
      TMatrixD R(2,2);
      R.Zero();
      R[0][0] = Ry;  // in cm^2 usually 4 1.1921e-07
      R[1][1] = Rz;//TMath::Sq(0.01); //TMath::Sq(0);  // in cm^2 usually 4
      R[1][0] = Ryz;
      R[0][1] = Ryz;
      //Rt.push_back(R);
      // add the TPCClusters and update the track parameters and uncertainties.  Put in additional terms to keep uncertainties from shrinking when
      // scattering and energy loss can change the track parameters along the way.

      // F = partial(updatefunc)/partial(parvec).  Update functions are in the comments below.

      
      TMatrixD F(5,5);
      TMatrixD FT(5,5);
      TVectorD parvec(5);
      parvec[0] = xyz_seed.X();
      parvec[1] = xyz_seed.Y();
      parvec[2] = curvature_seed;
      parvec[3] = phi_seed;
      parvec[4] = lambda_seed;

      
      //std::cout << " Parvec: y " << parvec[0] << " z " << parvec[1] << " c " << parvec[2] << " phi " << parvec[3] << " lambda " << parvec[4] << std::endl;
      TVectorD predstep(5);
      predstep.Zero();
      
      
      

      TMatrixD H(2,5);   // partial(obs)/partial(params)
      H.Zero();
      H[0][0] = 1;  // y
      H[1][1] = 1;  // z
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
      double fTPCClusterResolYZ=1;//1;
      

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
          
          if (fPrintLevel > 0)
            {
              std::cout << std::endl;
              std::cout << "Adding a new TPCCluster: " << xh << " " << yh << " " << zh << std::endl;
            }

          // for readability
          double curvature = parvec[2];
          double phi = parvec[3];
          double lambda = parvec[4];

         

          // update prediction to the plane containing x.  Maybe we need to find
          // the closest point on the helix to the TPCCluster we are adding,
          // and not necessarily force it to be at this x

          F.Zero();

          // y = yold + slope*dx*Sin(phi).   F[0][i] = dy/dtrackpar[i], where f is the update function slope*dx*Sin(phi)

          double slope = TMath::Tan(lambda);
          if (slope != 0)
            {
              slope = 1.0/slope;
            }
          else
            {
              slope = 1E9;
            }

          
          double dz;

          
          dz = zh - zpos;
          if (dz == 0) dz = 1E-3;
          
          
            
          
          //TODO check this -- are these the derivatives?
          // x = xold + dz/(slope*cos(phi))
          // slope = cot(lambda), so dslope/dlambda = -csc^2(lambda) = -1 - slope^2
          F[0][0] = 1.;
          F[0][3] = dz*TMath::Tan(phi)/(slope*TMath::Cos(phi));
          F[0][4] = dz*(1+slope*slope)/(slope*slope*TMath::Cos(phi));

          // y = yold + dz*tan(phi)
          F[1][1] = 1.;
          F[1][3] = dz/(TMath::Cos(phi)*TMath::Cos(phi));
          

          // curvature = old curvature -- doesn't change but put in an uncertainty
          F[2][2] = 1.;

          // phi = old phi + curvature*slope*dx
          // need to take the derivative of a product here
          F[3][2] = dz/(TMath::Cos(phi));
          F[3][3] = 1.+dz*TMath::Tan(phi)/(curvature*TMath::Cos(phi));


          // lambda -- same -- but put in an uncertainty in case it changes
          F[4][4] = 1.;

          // predicted step
          if (fPrintLevel > 1)
            {
              std::cout << "F Matrix: " << std::endl;
              F.Print();
              std::cout << "P Matrix: " << std::endl;
              P.Print();
            }
          

          predstep = parvec;
          
          predstep[0] += dz/(slope*TMath::Cos(phi));  // update x
          predstep[1] += dz*TMath::Tan(phi);  // update y
          predstep[3] += dz*curvature/(TMath::Cos(phi));        // update phi
          predstept.push_back(predstep);                       // update tree values
          
          if (fPrintLevel >0 )
            {
              std::cout << " Predstep: x " << predstep[0] << " y " << predstep[1] << " c " << predstep[2] << " phi " << predstep[3] << " lambda " << predstep[4] << std::endl;
            }
          std::cout << " Predstep: x " << predstep[0] << " y " << predstep[1] << " c " << predstep[2] << " phi " << predstep[3] << " lambda " << predstep[4] << std::endl;
          
          // equations from the extended Kalman filter
          FT.Transpose(F);
          PPred = F*P*FT + Q;
          PPredt.push_back(PPred);                           // update tree covariance

          if (fPrintLevel > 1)
            {
              std::cout << "PPred Matrix: " << std::endl;
              PPred.Print();
            }

          ytilde[0] = xh - predstep[0];
          ytilde[1] = yh - predstep[1];
          

          
          if (fPrintLevel > 0)
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
          
          double xprev = parvec[0];
          double yprev = parvec[1];
          parvec = predstep + K*ytilde;
          parvect.push_back(parvec);               //add best estimate to ttree
          P = (I-K*H)*PPred;
          Pt.push_back(P);                           //add covariance best estimate to tree
       
          zpos = zpos + dz;
          zpost.push_back(zpos);                     //add xpos to tree

          if (fPrintLevel > 0)
            {
              std::cout << "z: " << zpos << " dz: " << dz <<  std::endl;
              std::cout << " Parvec:   x " << parvec[0] << " y " << parvec[1] << " c " << parvec[2] << " phi " << parvec[3] << " lambda " << parvec[4] << std::endl;
            }
          
          //std::cout << " Updated xpos: " << xpos << " " << dx << std::endl;
          
          

          
        }
        
      
    }
    
    void kalman_helix_garlitenewparam()

    {
      // variables:  z is the independent variable
      // 0: x
      // 1: y
      // 2: curvature
      // 3: phi
      // 4: lambda = angle from the cathode plane
      // 5: z   /// added on to the end

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
      std::vector<XYZVector> xyz_plane_mean;
      std::vector<double> phi;
      std::vector<double> phi_plane;
      std::vector<double> lambda;
      std::vector<double> slope;
      double  charge;
      int nPoints;
      std::vector<int> nHits_perPlane;
      std::vector<XYZVector> xyz_firstinPlane;
      std::vector<double> phi_firstinPlane;
      std::vector<double> lambda_firstinPlane;
      std::vector<double> curvature_firstinPlane;

      ////Kalman

      double Ry = TMath::Sq(4);//TMath::Sq(2);//TMath::Sq(0.01);//TMath::Sq(3);
      double Rz = TMath::Sq(4);//TMath::Sq(0.01);//TMath::Sq(0.01); //1.1921e-07
      double Ryz = TMath::Sq(0);

      ////FWD
      XYZVector xyz_seed;
      double phi_seed,lambda_seed,curvature_seed;
      std::vector<double> xht,yht,zht,zpost;
      std::vector<TVectorD> parvect; //(5)
      std::vector<TVectorD> predstept; //(5)
      std::vector<TMatrixD> Pt;//(5,5);
      std::vector<TMatrixD> PPredt;//(5,5);
      std::vector<TMatrixD> Rt;//(2,2);

      ///BKW
      XYZVector xyz_seed_bkw;
      double phi_seed_bkw,lambda_seed_bkw,curvature_seed_bkw;
      std::vector<double> xht_bkw,yht_bkw,zht_bkw,zpost_bkw;
      std::vector<TVectorD> parvect_bkw; //(5)
      std::vector<TVectorD> predstept_bkw; //(5)
      std::vector<TMatrixD> Pt_bkw;//(5,5);
      std::vector<TMatrixD> PPredt_bkw;//(5,5);
      std::vector<TMatrixD> Rt_bkw;//(2,2);


      ///////MC
      t1s.Branch("xyz",&xyz);
      t1s.Branch("pxyz",&pxyz);
      t1s.Branch("phi",&phi);
      t1s.Branch("phi_plane",&phi_plane);
      t1s.Branch("lambda",&lambda);
      t1s.Branch("slope",&slope);
      t1s.Branch("xyz_plane",&xyz_plane);
      t1s.Branch("xyz_plane_mean",&xyz_plane_mean);
      t1s.Branch("charge",&charge);
      t1s.Branch("nPoints",&nPoints);
      t1s.Branch("nHits_perPlane",&nHits_perPlane);
      t1s.Branch("xyz_firstinPlane",&xyz_firstinPlane);
      t1s.Branch("phi_firstinPlane",&phi_firstinPlane);
      t1s.Branch("lambda_firstinPlane",&lambda_firstinPlane);
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
      t1s.Branch("phi_seed",&phi_seed);
      t1s.Branch("lambda_seed",&lambda_seed);
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
      t1s.Branch("phi_seed_bkw",&phi_seed_bkw);
      t1s.Branch("lambda_seed_bkw",&lambda_seed_bkw);
      t1s.Branch("curvature_seed_bkw",&curvature_seed_bkw);
      
     
      
      size_t nevents=1;
      
      for(size_t i=0;i<nevents;i++)
      {
      XYZVector xyztemp(gRandom->Rndm()*(Plane_XY_fid[1]-Plane_XY_fid[0])+Plane_XY_fid[0],gRandom->Rndm()*(Plane_XY_fid[3]-Plane_XY_fid[2])+Plane_XY_fid[2],gRandom->Rndm()*10+(Planes_Z[0]-10));
      xyz.push_back(xyztemp);
      
      charge = (gRandom->Rndm()<0.5) ? -1 : 1;
      
      //std::cout<<charge<<std::endl;
      
      
      double ptotal=0.1+gRandom->Rndm()*4; //in GeV/c2
      double px=-0.5+gRandom->Rndm();
      double py=-0.5+gRandom->Rndm();
      double pz=1.0+gRandom->Rndm();
      double norm=ptotal/sqrt(px*px+py*py+pz*pz);
      
      XYZVector pxyztemp(px*norm,py*norm,pz*norm);
      pxyz.push_back(pxyztemp);

      XYZVector pyz(0,pxyztemp.Y(),pxyztemp.Z()); ///pT
      XYZVector pxy(pxyztemp.X(),pxyztemp.Y(),0);
      
      double phitemp = TMath::ACos(pyz.Unit().Dot(Z_axis.Unit()));
      double curvaturetemp =charge*0.01*(0.3*B)/(sqrt(pyz.Mag2())); //in cm^-1
      double lambdatemp =TMath::ACos(pxy.Unit().Dot(Y_axis.Unit()));
      double slopetemp = TMath::Tan(lambdatemp);
      double xd =0;
      double yd =0;
      double zd =0;
      if (slopetemp != 0)
            {
              slopetemp = 1.0/slopetemp;
            }
          else
            {
              slopetemp = 1E9;
            }
      
      //TF2 *f2 = new TF2("f2","exp(-0.5*((x)/3)**2)*exp(-0.5*((y)/3)**2)",-20,20,-20,20);
      
      //TCanvas *mccanvaspT = new TCanvas("mccanvaspT","",1000,800);
      //f2->Draw();

      Bool_t InGAr=true;
      Bool_t InPlane[5]={0,0,0,0,0};
      double dz=0.020000000000000;
     

      while(InGAr)
        {
            
             //Randomized x
            //double dx=0.04;
            
            xyztemp.SetX(xyztemp.X()+dz/(slopetemp*TMath::Cos(phitemp)));
            xyztemp.SetY(xyztemp.Y()+dz*TMath::Sin(phitemp)/TMath::Cos(phitemp));
            xyztemp.SetZ(xyztemp.Z()+dz);
            phitemp+=dz*curvaturetemp/TMath::Cos(phitemp);
            pxyztemp.SetY(TMath::Sqrt(pyz.Mag2())*TMath::Sin(phitemp));
            pxyztemp.SetZ(TMath::Sqrt(pyz.Mag2())*TMath::Cos(phitemp));

            xyz.push_back(xyztemp);
            pxyz.push_back(pxyztemp);
            phi.push_back(phitemp);
            slope.push_back(slopetemp);

            //std::cout<<"XYZphi MC: "<< xyztemp.X()<< " " << xyztemp.Y()<< " " << xyztemp.Z()<< " " <<phitemp<<std::endl;
            //std::cout<<"XYZ MC:"<<x<<std::endl;
          
            
            
            for(size_t p=0;p<5;p++)
            {
              
              if(xyztemp.X()>Plane_XY[0] && xyztemp.X()<Plane_XY[1] && xyztemp.Y()>Plane_XY[2] && xyztemp.Y()<Plane_XY[3])
              {
                if(xyztemp.Z()>Planes_Z[p*2] && xyztemp.Z()<Planes_Z[p*2+1] && InPlane[p]==0)
                  {
                    InPlane[p]=1;
                    int hitsinplane=1+gRandom->Gaus(2,1);
                    //std::cout<<hitsinplane<<std::endl;
                    nHits_perPlane.push_back(hitsinplane);
                    xyz_firstinPlane.push_back(xyztemp);
                    curvature_firstinPlane.push_back(curvaturetemp);
                    lambda_firstinPlane.push_back(lambdatemp);
                    phi_firstinPlane.push_back(phitemp);

                    XYZVector xyz_plane_temp=xyztemp;
                    std::vector<XYZVector> xyz_plane_m;
                    double phi_plane_temp=phitemp;
                    double lamda_plane_temp=lambdatemp;
                    double slope_plane_temp=slopetemp;
                    double curvature_plane_temp=curvaturetemp;

                    //std::cout<<"XYZ MC: "<< xyztemp.X()<< " " << xyztemp.Y()<< " " << xyztemp.Z()<< std::endl;
                    for(int i=0;i<hitsinplane;i++)
                      {
                        //XYZVector hit(gRandom->Gaus(xyztemp.X(),1),gRandom->Gaus(xyztemp.Y(),1),Planes_Z[p*2]+Plane_thick*gRandom->Rndm());
                        double dzp=Plane_thick/(hitsinplane+1);

                        xyz_plane_temp.SetX(xyz_plane_temp.X()+dzp/(slope_plane_temp*TMath::Cos(phi_plane_temp)));
                        xyz_plane_temp.SetY(xyz_plane_temp.Y()+dzp*TMath::Sin(phi_plane_temp)/TMath::Cos(phi_plane_temp));
                        xyz_plane_temp.SetZ(xyz_plane_temp.Z()+dzp);
                        phi_plane_temp+=dzp*curvature_plane_temp/TMath::Cos(phi_plane_temp);

                        //std::cout<<"XYZ MC: "<< xyz_plane_temp.X()<< " " << xyz_plane_temp.Y()<< " " << xyz_plane_temp.Z()<<" "<<phi_plane_temp<< std::endl;

                        //XYZVector hit(gRandom->Gaus(xyz_plane_temp.X(),0.5),gRandom->Gaus(xyz_plane_temp.Y(),0.5),gRandom->Gaus(xyz_plane_temp.Z(),0.5));
                        
                        xyz_plane.push_back( xyz_plane_temp);
                        phi_plane.push_back(phi_plane_temp);
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
        /*
        xyz_seed.SetX(gRandom->Gaus(xyz_firstinPlane.at(0).X(),1));
        xyz_seed.SetY(gRandom->Gaus(xyz_firstinPlane.at(0).Y(),1));
        xyz_seed.SetZ(gRandom->Gaus(xyz_firstinPlane.at(0).Z(),1));
        curvature_seed=gRandom->Gaus(curvature_firstinPlane.at(0),curvature_firstinPlane.at(0)/20);
        phi_seed=gRandom->Gaus(phi_firstinPlane.at(0),phi_firstinPlane.at(0)/20);
        lambda_seed=gRandom->Gaus(lambda_firstinPlane.at(0),lambda_firstinPlane.at(0)/20);
        */
        xyz_seed=xyz_firstinPlane.at(0);
        curvature_seed=curvature_firstinPlane.at(0);
        phi_seed=phi_firstinPlane.at(0);
        lambda_seed=lambda_firstinPlane.at(0);
        double forward=1.;
        //std::cout<<std::endl<<"Real Values: "<<xyz_firstinPlane.at(0).X()<<" "<<xyz_firstinPlane.at(0).Y()<<" "<<xyz_firstinPlane.at(0).Z()<<" "<<curvature_firstinPlane.at(0)<<" "<<phi_firstinPlane.at(0)<<" "<<lambda_firstinPlane.at(0)<<std::endl;
        KalmanFit(xht,yht,zht,parvect,predstept,Pt,PPredt,Rt,zpost,xyz_plane,phi_plane,Ry,Rz,Ryz,xyz_seed,lambda_seed,curvature_seed,phi_seed,forward);

        //////BKW
        //std::cout<<xht.size()<<" "<<xyz_plane.size()<<" "<<nHits_perPlane.size()<<std::endl;
        xyz_seed_bkw.SetX(xht.at(xht.size()-1));
        xyz_seed_bkw.SetY(parvect.at(parvect.size()-1)[0]);
        xyz_seed_bkw.SetZ(parvect.at(parvect.size()-1)[1]);
        curvature_seed_bkw=parvect.at(parvect.size()-1)[2];
        phi_seed_bkw=parvect.at(parvect.size()-1)[3];
        lambda_seed_bkw=parvect.at(parvect.size()-1)[4];
        std::vector<XYZVector> xyz_plane_bkw=xyz_plane;
        std::reverse(xyz_plane_bkw.begin(),xyz_plane_bkw.end());
        //std::cout<<xyz_seed_bkw.X()<<" "<<xyz_seed_bkw.X()<<" "<<xyz_seed_bkw.X()<<" "<<curvature_seed<<" "<<phi_seed_bkw<<" "<<lambda_seed_bkw<<std::endl;
        double backwards=-1.;
        //KalmanFit(xht_bkw,yht_bkw,zht_bkw,parvect_bkw,predstept_bkw,Pt_bkw,PPredt_bkw,Rt_bkw,zpost_bkw,xyz_plane_bkw,Ry,Rz,Ryz,xyz_seed_bkw,lambda_seed_bkw,curvature_seed_bkw,phi_seed_bkw,backwards);
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
      parvect.clear();
      predstept.clear();
      Pt.clear();//(5,5);
      PPredt.clear();//(5,5);
      Rt.clear();//(2,2);
      phi.clear();
      lambda.clear();
      slope.clear();
      xyz_firstinPlane.clear();
      phi_firstinPlane.clear();
      lambda_firstinPlane.clear();
      curvature_firstinPlane.clear();
      nHits_perPlane.clear();
      xyz_plane_mean.clear();
      
    }
    
    t1s.Write();
          
      
  }

    

    