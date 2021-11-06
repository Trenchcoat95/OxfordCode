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
#include "TVectorF.h"
#include "TMatrix.h"
#include "TMath.h"
#include "TVector3.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TF2.h"
#include "Math/Vector3D.h"
using namespace ROOT::Math;







//--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    // KalmanFit does a forwards or backwards Kalman fit using the sorted TPCCluster list
    // variables:  x is the independent variable
    // 0: y
    // 1: z
    // 2: curvature
    // 3: phi
    // 4: lambda

    void KalmanFit( TTree &t, 
                   std::vector<float> &xht, 
                   std::vector<float> &yht, 
                   std::vector<float> &zht, 
                   std::vector<TVectorF> &parvect,
                   std::vector<TVectorF> &predstept,
                   std::vector<TMatrixF> &Pt,
                   std::vector<TMatrixF> &PPredt,
                   std::vector<TMatrixF> &Rt,
                   std::vector<float> &xpost, 
                   TTree &t1s,
                   std::vector<XYZVector> * & xyz_plane,
                   float Ry,
                   float Rz,
                   float Ryz,
                   std::vector<XYZVector> * & xyz_firstinPlane,
                   std::vector<float> * & lambda_firstinPlane,
                   std::vector<float> * & curvature_firstinPlane,
                   std::vector<float> * & phi_firstinPlane,
                   XYZVector  & xyz_seed,
                   float  & lambda_seed,
                   float & curvature_seed,
                   float  & phi_seed
                   

                                )
    {

      

      size_t nentries = t1s.GetEntries();
      float fPrintLevel=0;

      bool showprog = true;  
      if(showprog==true) std::cout<<"Progress:  "<<std::endl;
      
      
      for (size_t i=0;i<nentries;i++)
      {
      // estimate curvature, lambda, phi, xpos from the initial track parameters
      /*
      float curvature_init= 0.001336709363386035;         ///need to check what this is
      float phi_init = 2.415579080581665;
      float lambda_init = 0.037335269153118134;
      float xpos_init= -23.076112747192383;
      float ypos_init=-340.07403564453125;
      float zpos_init=1638.2122802734375;
      */

      int prog = 100*i/nentries;
      std::string strprog = std::to_string(prog);
      if(showprog==true) std::cout<<strprog<<"%";
      

      

      t1s.GetEntry(i);

      //std::cout<<i<<std::endl;

      float curvature_init= gRandom->Gaus(curvature_firstinPlane->at(0),curvature_firstinPlane->at(0)/20);         ///need to check what this is
      float phi_init = gRandom->Gaus(phi_firstinPlane->at(0),phi_firstinPlane->at(0)/20);
      float lambda_init = gRandom->Gaus(lambda_firstinPlane->at(0),lambda_firstinPlane->at(0)/20);
      float xpos_init= gRandom->Gaus(xyz_firstinPlane->at(0).X(),xyz_firstinPlane->at(0).X()/20);
      float ypos_init= gRandom->Gaus(xyz_firstinPlane->at(0).Y(),xyz_firstinPlane->at(0).Y()/20);
      float zpos_init= gRandom->Gaus(xyz_firstinPlane->at(0).Z(),xyz_firstinPlane->at(0).Z()/20);

      xyz_seed.SetXYZ(xpos_init,ypos_init,zpos_init);
      lambda_seed=lambda_init;
      curvature_seed=curvature_init;
      phi_seed=phi_init;

      // Kalman fitter variables

      float xpos = xpos_init;

      TMatrixF P(5,5);  // covariance matrix of parameters
      // fill in initial guesses -- generous uncertainties on first value.
      P.Zero();
      P[0][0] = TMath::Sq(1);   // initial position uncertainties -- y
      P[1][1] = TMath::Sq(1);   // and z
      P[2][2] = TMath::Sq(.5);  // curvature of zero gets us to infinite momentum, and curvature of 2 is curled up tighter than the pads
      P[3][3] = TMath::Sq(.5);  // phi uncertainty
      P[4][4] = TMath::Sq(.5);  // lambda uncertainty
      TMatrixF PPred(5,5);
      PPred.Zero();

      // per-step additions to the covariance matrix 
      //need to check this
      TMatrixF Q(5,5);
      Q.Zero();
      Q[2][2] = 1e-09;     // allow for some curvature uncertainty between points
      Q[3][3] = 1e-09;      // phi
      Q[4][4] = 0.0001;   // lambda

      // Noise covariance on the measured points.
      // 16 cm2 initially, might reasonably be lowered to typicalResidual near line 552-67
      TMatrixF R(2,2);
      R.Zero();
      R[0][0] = Ry;  // in cm^2 usually 4 1.1921e-07
      R[1][1] = Rz;//TMath::Sq(0.01); //TMath::Sq(0);  // in cm^2 usually 4
      R[1][0] = Ryz;
      R[0][1] = Ryz;
      //Rt.push_back(R);
      // add the TPCClusters and update the track parameters and uncertainties.  Put in additional terms to keep uncertainties from shrinking when
      // scattering and energy loss can change the track parameters along the way.

      // F = partial(updatefunc)/partial(parvec).  Update functions are in the comments below.
      TMatrixF F(5,5);
      TMatrixF FT(5,5);
      TVectorF parvec(5);
      parvec[0] = ypos_init;
      parvec[1] = zpos_init;
      parvec[2] = curvature_init;
      parvec[3] = phi_init;
      parvec[4] = lambda_init;
      //std::cout << " Parvec: y " << parvec[0] << " z " << parvec[1] << " c " << parvec[2] << " phi " << parvec[3] << " lambda " << parvec[4] << std::endl;
      TVectorF predstep(5);
      predstep.Zero();
      /*
      //Filling the TTree for step zero
      Pt.push_back(P);
      PPredt.push_back(PPred);
      parvect.push_back(parvec);
      predstept.push_back(predstep);
      xht.push_back(0);
      yht.push_back(0);
      zht.push_back(0);
      Rt.push_back(R);
      xpost.push_back(xpos);
      t.Fill();
      */
      

      TMatrixF H(2,5);   // partial(obs)/partial(params)
      H.Zero();
      H[0][0] = 1;  // y
      H[1][1] = 1;  // z
      TMatrixF HT(5,2);

      TVectorF zv(2);
      TVectorF ytilde(2);
      TVectorF hx(2);
      TMatrixF S(2,2);
      TMatrixF K(5,2);

      TMatrixF I(5,5);
      I.Zero();
      for (int i=0;i<5;++i) I[i][i] = 1;


      
      //Create a parametrized helix as a substitute for the measurement
      //t1s.GetEntry(0);
      float xh = xyz_plane->at(0).X(); // -23.3112;  ///need to check what this is
      float yh = xyz_plane->at(0).Y(); // -337.045;
      float zh = xyz_plane->at(0).Z(); // v1634.73;
      float fTPCClusterResolYZ=1;
      float fTPCClusterResolX=0.5;

      for (size_t iTPCCluster=1; iTPCCluster<xyz_plane->size(); ++iTPCCluster)
        {
          
          if (iTPCCluster>1)
          {
            
            xh=xyz_plane->at(iTPCCluster).X(); 
            yh=xyz_plane->at(iTPCCluster).Y();
            zh=xyz_plane->at(iTPCCluster).Z();

          }
          
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
          float curvature = parvec[2];
          float phi = parvec[3];
          float lambda = parvec[4];

          // update prediction to the plane containing x.  Maybe we need to find
          // the closest point on the helix to the TPCCluster we are adding,
          // and not necessarily force it to be at this x

          F.Zero();

          // y = yold + slope*dx*Sin(phi).   F[0][i] = dy/dtrackpar[i], where f is the update function slope*dx*Sin(phi)

          float slope = TMath::Tan(lambda);
          if (slope != 0)
            {
              slope = 1.0/slope;
            }
          else
            {
              slope = 1E9;
            }

          // relocate dx to be the location along the helix which minimizes
		      // [ (Xhit -Xhelix)/sigmaX ]^2 + [ (Yhit -Yhelix)/sigmaY ]^2 + [ (Zhit -Zhelix)/sigmaZ ]^2
          // Linearize for now near xpos:
          //        x = xpos + dx
          //        y = parvec[0] + slope * dx * sin(phi)
          //        z = parvec[1] + slope * dx * cos(phi)
          // parvec is updated as the fit progresses so the 'zero point' where y_0, z_0, phi_0
          // are defined is at the end of the fit, not at the place where the fit begins.
          //
          // old calc was just based on TPCCluster position in x:
          // float dx = xh - xpos;
          float dx;

          
          dx = xh - xpos;
          if (dx == 0) dx = 1E-3;
          
          
            
          /*
            float dxnum = (slope/(fTPCClusterResolYZ*fTPCClusterResolYZ))*( (yh - parvec[0])*TMath::Sin(phi) + (zh - parvec[1])*TMath::Cos(phi) )
            + (xh - xpos)/(fTPCClusterResolX*fTPCClusterResolX);
            float dxdenom = slope*slope/(fTPCClusterResolYZ*fTPCClusterResolYZ) + 1.0/(fTPCClusterResolX*fTPCClusterResolX);
            dx = dxnum/dxdenom;
            if (dx == 0) dx = 1E-3;
          */
            //int precision = std::numeric_limits<double>::max_digits10;
            //std::cout << std::setprecision(precision);
            //std::cout<<"dxnum "<<dxnum<<std::endl;
            //std::cout<<"dxdenom "<<dxdenom<<std::endl;
            //std::cout<<"dx "<<dx<<std::endl;
            //std::cout<<slope<<" "<<fTPCClusterResolYZ<<" "<<yh<<" "<<parvec[0]<<" "<<TMath::Sin(phi)<<" "<<zh<<" "<<parvec[1]<<" "<<TMath::Cos(phi)<<" "<<xh<<" "<<xpos<<" "<<fTPCClusterResolX<<" "<<phi<<std::endl;
          //std::cout<<"dyzpart "<<((dxnum-(xh - xpos)/(fTPCClusterResolX*fTPCClusterResolX))/dxdenom)<<" dxpart "<<((xh - xpos)/(fTPCClusterResolX*fTPCClusterResolX))/dxdenom<<std::endl;
          //std::cout << "dxdenom, dxnum: " << dxdenom << " " << dxnum << std::endl;
          //std::cout << "Track pos: " << xpos << " " << parvec[0] << " " << parvec[1] << " " << " TPCCluster pos: " << xh << " " << yh << " " << zh << std::endl;
          //std::cout << "dx old and new: " << xh - xpos << " " << dx << std::endl<<std::endl;
          

          //TODO check this -- are these the derivatives?
          // y = yold + dx*slope*TMath::Sin(phi)
          // slope = cot(lambda), so dslope/dlambda = -csc^2(lambda) = -1 - slope^2
          F[0][0] = 1.;
          F[0][3] = dx*slope*TMath::Cos(phi);
          F[0][4] = dx*TMath::Sin(phi)*(-1.0-slope*slope);

          // z = zold + slope*dx*Cos(phi)
          F[1][1] = 1.;
          F[1][3] = -dx*slope*TMath::Sin(phi);
          F[1][4] = dx*TMath::Cos(phi)*(-1.0-slope*slope);

          // curvature = old curvature -- doesn't change but put in an uncertainty
          F[2][2] = 1.;

          // phi = old phi + curvature*slope*dx
          // need to take the derivative of a product here
          F[3][2] = dx*slope;
          F[3][3] = 1.;
          F[3][4] = dx*curvature*(-1.0-slope*slope);

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
          if (fPrintLevel > 0)
            {
              std::cout << "x: " << xpos << " dx: " << dx <<  std::endl;
              std::cout << " Parvec:   y " << parvec[0] << " z " << parvec[1] << " c " << parvec[2] << " phi " << parvec[3] << " lambda " << parvec[4] << std::endl;
            }

          predstep = parvec;
          
          predstep[0] += slope*dx*TMath::Sin(phi);  // update y
          predstep[1] += slope*dx*TMath::Cos(phi);  // update z
          predstep[3] += slope*dx*curvature;        // update phi
          predstept.push_back(predstep);                       // update tree values
          
          if (fPrintLevel > 1)
            {
              std::cout << " Predstep: y " << predstep[0] << " z " << predstep[1] << " c " << predstep[2] << " phi " << predstep[3] << " lambda " << predstep[4] << std::endl;
            }
          // equations from the extended Kalman filter
          FT.Transpose(F);
          PPred = F*P*FT + Q;
          PPredt.push_back(PPred);                           // update tree covariance

          if (fPrintLevel > 1)
            {
              std::cout << "PPred Matrix: " << std::endl;
              PPred.Print();
            }

          ytilde[0] = yh - predstep[0];
          ytilde[1] = zh - predstep[1];
          

          
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
          
          float yprev = parvec[0];
          float zprev = parvec[1];
          parvec = predstep + K*ytilde;
          parvect.push_back(parvec);               //add best estimate to ttree
          P = (I-K*H)*PPred;
          Pt.push_back(P);                           //add covariance best estimate to tree
       
          xpos = xpos + dx;
          xpost.push_back(xpos);                     //add xpos to tree

          
          
          //std::cout << " Updated xpos: " << xpos << " " << dx << std::endl;
          
          

          
        }
        t.Fill();
        if(showprog==true) std::cout << std::string(strprog.length(),'\b')<<"\b";
      }
      
    }
    
    void kalman_helix_garlite()

    {
      // variables:  x is the independent variable
      // 0: y
      // 1: z
      // 2: curvature
      // 3: phi
      // 4: lambda = angle from the cathode plane
      // 5: x   /// added on to the end

      ///Prepare a simple tree with just the coordinates

      

      //Right now with perfect helix

      TFile fs("garlitetest.root","recreate");
      TTree t1s("t1s","helix simple tree");
      float  Plane_XY[4] = {-300, +300, -400, 100}; //all in cm
      float  Plane_XY_fid[4] = {-200, +200, -300, 0};
      float  Planes_Z[10] =   {1244,1248,1334,1338,1484,1488,1634,1638,1724,1728};
      XYZVector GArCenter(0,-150.473,1486); 
      float GAr_r = 349.9;
      float GAr_L = 669.6;
      float Plane_thick = 4;
      float  B=0.5; //in Tesla
      XYZVector Z_axis(0,0,1);
      XYZVector Y_axis(0,1,0);
      
      
      

      std::vector<XYZVector> *xyz=0;
      std::vector<XYZVector> *pxyz=0;
      std::vector<XYZVector> *xyz_plane=0;
      std::vector<float> *phi=0;
      std::vector<float> *lambda=0;
      std::vector<float> *slope=0;
      float  charge;
      int nPoints;
      int nHits;

      std::vector<XYZVector> *xyz_firstinPlane=0;
      std::vector<float> *phi_firstinPlane=0;
      std::vector<float> *lambda_firstinPlane=0;
      std::vector<float> *curvature_firstinPlane=0;
      //TBranch     *b_xyz;
      //TBranch     *b_pxyz;

      t1s.Branch("xyz",&xyz);
      t1s.Branch("pxyz",&pxyz);
      t1s.Branch("phi",&phi);
      t1s.Branch("lambda",&lambda);
      t1s.Branch("slope",&slope);
      t1s.Branch("xyz_plane",&xyz_plane);
      t1s.Branch("charge",&charge);
      t1s.Branch("nPoints",&nPoints);
      t1s.Branch("nHits",&nHits);

      t1s.Branch("xyz_firstinPlane",&xyz_firstinPlane);
      t1s.Branch("phi_firstinPlane",&phi_firstinPlane);
      t1s.Branch("lambda_firstinPlane",&lambda_firstinPlane);
      t1s.Branch("curvature_firstinPlane",&curvature_firstinPlane);
      
      size_t nevents=100;
      
      for(size_t i=0;i<nevents;i++)
      {
      XYZVector xyztemp(gRandom->Rndm()*(Plane_XY_fid[1]-Plane_XY_fid[0])+Plane_XY_fid[0],gRandom->Rndm()*(Plane_XY_fid[3]-Plane_XY_fid[2])+Plane_XY_fid[2],gRandom->Rndm()*10+(Planes_Z[0]-10));
      xyz->push_back(xyztemp);
      
      charge = (gRandom->Rndm()<0.5) ? -1 : 1;
      
      //std::cout<<charge<<std::endl;
      
      
      float ptotal=0.1+gRandom->Rndm()*4; //in GeV/c2
      float px=-0.5+gRandom->Rndm();
      float py=-0.5+gRandom->Rndm();
      float pz=1.0+gRandom->Rndm();
      float norm=ptotal/sqrt(px*px+py*py+pz*pz);
      
      XYZVector pxyztemp(px*norm,py*norm,pz*norm);
      pxyz->push_back(pxyztemp);

      XYZVector pyz(0,pxyztemp.Y(),pxyztemp.Z());
      XYZVector pxy(pxyztemp.X(),pxyztemp.Y(),0);
      
      float phitemp = TMath::ACos(pyz.Unit().Dot(Z_axis.Unit()));
      float curvaturetemp =charge*0.01*(0.3*B)/ptotal; //in cm^-1
      float lambdatemp =TMath::ACos(pxy.Unit().Dot(Y_axis.Unit()));
      float slopetemp = TMath::Tan(lambdatemp);
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
      

      while(InGAr)
        {
          //if (iTPCCluster>1)
          //{
            float dx=0.02+0.04*gRandom->Rndm(); //Randomized x
            //float dx=0.04;
            xyztemp.SetX(xyztemp.X()+dx);
            xyztemp.SetY(xyztemp.Y()+slopetemp*dx*TMath::Sin(phitemp));
            xyztemp.SetZ(xyztemp.Z()+slopetemp*dx*TMath::Cos(phitemp));
            //z=rnd->Gaus(ztemp,3);    //Smeared z 
            //f2->GetRandom2(yd,zd);
            //x=xd+xtemp;
            //y=yd+ytemp;
            //z=zd+ztemp;
            phitemp+=slopetemp*dx*curvaturetemp;
            pxyztemp.SetY(TMath::Sqrt(pyz.Mag2())*TMath::Sin(phitemp));
            pxyztemp.SetZ(TMath::Sqrt(pyz.Mag2())*TMath::Cos(phitemp));

            xyz->push_back(xyztemp);
            pxyz->push_back(pxyztemp);
            phi->push_back(phitemp);
            slope->push_back(slopetemp);

          //}
          
            
            
            for(size_t p=0;p<5;p++)
            {
              
              if(xyztemp.X()>Plane_XY[0] && xyztemp.X()<Plane_XY[1] && xyztemp.Y()>Plane_XY[2] && xyztemp.Y()<Plane_XY[3])
              {
                if(xyztemp.Z()>Planes_Z[p*2] && xyztemp.Z()<Planes_Z[p*2+1] && InPlane[p]==0)
                  {
                    InPlane[p]=1;
                    float hitsinplane=gRandom->Gaus(2,1);
                    xyz_firstinPlane->push_back(xyztemp);
                    curvature_firstinPlane->push_back(slopetemp);
                    lambda_firstinPlane->push_back(lambdatemp);
                    phi_firstinPlane->push_back(phitemp);

                    //std::cout<<"XYZ MC: "<< xyztemp.X()<< " " << xyztemp.Y()<< " " << xyztemp.Z()<< std::endl;
                    for(int i=0;i<hitsinplane;i++)
                      {
                        XYZVector hit(gRandom->Gaus(xyztemp.X(),1),gRandom->Gaus(xyztemp.Y(),1),Planes_Z[p*2]+Plane_thick*gRandom->Rndm());
                        xyz_plane->push_back(hit);
                        //std::cout<<"XYZ hit: "<< hit.X()<< " " << hit.Y()<< " " << hit.Z()<< std::endl;
                      }
                    //std::cout<< " "<<std::endl;
                    
                  }
              }
            }
            

            float r=sqrt(TMath::Power((xyztemp.Y()-GArCenter.Y()),2)+TMath::Power((xyztemp.Z()-GArCenter.Z()),2));

            if(xyztemp.X()>0.5*GAr_L || xyztemp.X()<-0.5*GAr_L || r>GAr_r) InGAr=false;

        }

      nPoints=xyz->size();
      nHits= xyz_plane->size();

      //std::cout<<"Recorded hit z: ";
      //for(int i=0;i<xyz_plane->size();i++) std::cout<<xyz_plane->at(i).Z()<<" ";
      //std::cout<<""<<std::endl;
      
      t1s.Fill();
      xyz->clear();
      pxyz->clear();
      xyz_plane->clear();
    }
      t1s.Write();
      
      
      float Ry = TMath::Sq(3);//TMath::Sq(2);//TMath::Sq(0.01);//TMath::Sq(3);
      float Rz = TMath::Sq(2);//TMath::Sq(0.01);//TMath::Sq(0.01); //1.1921e-07
      float Ryz = TMath::Sq(0);

      TFile f("garlite_kalman.root","recreate");
      TTree t1_FWD("t1_FWD","Forward fitter tree");
      
      XYZVector xyz_seed;
      float phi_seed,lambda_seed,curvature_seed;
      std::vector<float> xht,yht,zht,xpost;
      std::vector<TVectorF> parvect(5);
      std::vector<TVectorF> predstept(5);
      std::vector<TMatrixF> Pt;//(5,5);
      std::vector<TMatrixF> PPredt;//(5,5);
      std::vector<TMatrixF> Rt;//(2,2);
      
      t1_FWD.Branch("xht",&xht);
      t1_FWD.Branch("yht",&yht);
      t1_FWD.Branch("zht",&zht);
      t1_FWD.Branch("xpost",&xpost);
      t1_FWD.Branch("parvect",&parvect);
      t1_FWD.Branch("predstept",&predstept);
      t1_FWD.Branch("Pt",&Pt);
      t1_FWD.Branch("PPredt",&PPredt);
      t1_FWD.Branch("Rt",&Rt);

      t1_FWD.Branch("xyz_seed",&xyz_seed);
      t1_FWD.Branch("phi_seed",&phi_seed);
      t1_FWD.Branch("lambda_seed",&lambda_seed);
      t1_FWD.Branch("curvature_seed",&curvature_seed);
      
      KalmanFit(t1_FWD,xht,yht,zht,parvect,predstept,Pt,PPredt,Rt,xpost,t1s,xyz,Ry,Rz,Ryz,xyz_firstinPlane,lambda_firstinPlane,curvature_firstinPlane,phi_firstinPlane,xyz_seed,lambda_seed,curvature_seed,phi_seed);
      
      
      t1_FWD.Write();
      
      
    }

    

    