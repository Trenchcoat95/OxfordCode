////////////////////////////////////////////////////////////////////////
// Class:       tpctrackfit2
// Plugin Type: producer (art v3_00_00)
// File:        tpctrackfit2_module.cc
//
// Generated at Tue Feb  5 11:34:54 2019 by Thomas Junk using cetskelgen
// from cetlib version v3_04_00.
////////////////////////////////////////////////////////////////////////



// ROOT includes

#include "TVectorF.h"
#include "TMatrix.h"
#include "TMath.h"
#include "TVector3.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "kalman_helix_algs.h"


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
                   float &xht, 
                   float &yht, 
                   float &zht, 
                   TVectorF &parvect,
                   TVectorF &predstept,
                   TMatrixF &Pt,
                   TMatrixF &PPredt,
                   TMatrixF &Rt,
                   float &xpost,
                   size_t n, 
                   TTree &t1s,
                   float &x,
                   float &y,
                   float &z
                                )
    {

      
      
      float fPrintLevel=0;
      
      // estimate curvature, lambda, phi, xpos from the initial track parameters
      /*
      float curvature_init= 0.001336709363386035;         ///need to check what this is
      float phi_init = 2.415579080581665;
      float lambda_init = 0.037335269153118134;
      float xpos_init= -23.076112747192383;
      float ypos_init=-340.07403564453125;
      float zpos_init=1638.2122802734375;
      */

      float curvature_init=0.1;
      float phi_init = 0;
      float lambda_init = 0;
      float xpos_init=0;
      float ypos_init=0;
      float zpos_init=0;

      initial_trackpar_estimate(t1s,x,y,z,n,curvature_init,phi_init,lambda_init,xpos_init,ypos_init,zpos_init,fPrintLevel);
      
      std::cout<<"Initial estimate: "<<curvature_init<<" "<<phi_init<<" "<<lambda_init<<" "<<xpos_init<<" "<<ypos_init<<" "<<zpos_init<<std::endl;

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
      R[0][0] = TMath::Sq(0.001);  // in cm^2 usually 4
      R[1][1] = TMath::Sq(3); //TMath::Sq(0);  // in cm^2 usually 4
      Rt=R;
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

      //Filling the TTree for step zero
      Pt=P;
      PPredt.Zero();
      parvect=parvec;
      predstept.Zero();
      xht=0;
      yht=0;
      zht=0;
      Rt=R;
      xpost=xpos;
      t.Fill();

      

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
      int i=0;
      t1s.GetEntry(i);
      float xh = x; // -23.3112;  ///need to check what this is
      float yh = y; // -337.045;
      float zh = z; // v1634.73;
      float fTPCClusterResolYZ=1;
      float fTPCClusterResolX=0.5;
      //int i=0;
      for (size_t iTPCCluster=1; iTPCCluster<n; ++iTPCCluster)
        {
          i+=1;
          if (true)
          {
            t1s.GetEntry(i);
            xh=x; 
            yh=y;
            zh=z;

          }
          
          //std::cout<<"iTPCCluster "<<iTPCCluster<<std::endl;
          //std::cout << "xh,yh,zh: "<< xh <<" "<<yh<<" "<<zh<<std::endl;
          //std::cout << " Parvec: y " << parvec[0] << " z " << parvec[1] << " c " << parvec[2] << " phi " << parvec[3] << " lambda " << parvec[4] << std::endl;
          
          
          xht=xh;  //add measured position to ttree
          yht=yh;
          zht=zh;
          
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
          
            int precision = std::numeric_limits<double>::max_digits10;
            std::cout << std::setprecision(precision);
            std::cout<<"dxnum "<<dxnum<<std::endl;
            std::cout<<"dxdenom "<<dxdenom<<std::endl;
            std::cout<<"dx "<<dx<<std::endl;
            std::cout<<slope<<" "<<fTPCClusterResolYZ<<" "<<yh<<" "<<parvec[0]<<" "<<TMath::Sin(phi)<<" "<<zh<<" "<<parvec[1]<<" "<<TMath::Cos(phi)<<" "<<xh<<" "<<xpos<<" "<<fTPCClusterResolX<<" "<<phi<<std::endl;
          //std::cout<<"dyzpart "<<((dxnum-(xh - xpos)/(fTPCClusterResolX*fTPCClusterResolX))/dxdenom)<<" dxpart "<<((xh - xpos)/(fTPCClusterResolX*fTPCClusterResolX))/dxdenom<<std::endl;
          //std::cout << "dxdenom, dxnum: " << dxdenom << " " << dxnum << std::endl;
          //std::cout << "Track pos: " << xpos << " " << parvec[0] << " " << parvec[1] << " " << " TPCCluster pos: " << xh << " " << yh << " " << zh << std::endl;
          //std::cout << "dx old and new: " << xh - xpos << " " << dx << std::endl<<std::endl;
          */

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
          predstept=predstep;                       // update tree values
          
          if (fPrintLevel > 1)
            {
              std::cout << " Predstep: y " << predstep[0] << " z " << predstep[1] << " c " << predstep[2] << " phi " << predstep[3] << " lambda " << predstep[4] << std::endl;
            }
          // equations from the extended Kalman filter
          FT.Transpose(F);
          PPred = F*P*FT + Q;
          PPredt = PPred;                           // update tree covariance

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
          parvect = parvec;               //add best estimate to ttree
          P = (I-K*H)*PPred;
          Pt=P;                           //add covariance best estimate to tree
       
          xpos = xpos + dx;
          xpost=xpos;                     //add xpos to tree

          
          t.Fill();
          //std::cout << " Updated xpos: " << xpos << " " << dx << std::endl;
          
          

          
        }
      
    }
    
    void kalman_helix_test()

    {
      // variables:  x is the independent variable
      // 0: y
      // 1: z
      // 2: curvature
      // 3: phi
      // 4: lambda = angle from the cathode plane
      // 5: x   /// added on to the end

      ///Prepare a simple tree with just the coordinates
      size_t n = 452;

      //Right now with perfect helix

      TFile fs("m_perfect_helix_simple_rndx_smz3_stdK.root","recreate");
      TTree t1s("t1s","helix simple tree");
      float x,y,z;
      t1s.Branch("x",&x,"x/F");
      t1s.Branch("y",&y,"y/F");
      t1s.Branch("z",&z,"z/F");
      TRandom *rnd = new TRandom();
      x = -23.311233520507812;  
      //float xtemp=x;
      y = -337.04501342773438;
      float ytemp=y;
      z = 1634.73486328125;
      float ztemp=z;
      float phi = 6;
      float curvature =-0.014;
      float lambda =-0.05;
      float slope = TMath::Tan(lambda);
      if (slope != 0)
            {
              slope = 1.0/slope;
            }
          else
            {
              slope = 1E9;
            }
      int precision = std::numeric_limits<double>::max_digits10;
      std::cout << std::setprecision(precision);
      std::cout<<"Initial values: "<<curvature<<" "<<phi<<" "<<lambda<<" "<<x<<" "<<y<<" "<<z<<std::endl;
      for (size_t iTPCCluster=1; iTPCCluster<n; ++iTPCCluster)
        {
          if (iTPCCluster>1)
          {
            float dx=0.02+0.04*rnd->Rndm(); //Randomized x
            //float dx=0.04;
            x+=dx;
            //x=rnd->Gaus(xtemp,0.04);   //Smeared x
            y+=slope*dx*TMath::Sin(phi);
            //y=rnd->Gaus(ytemp,3);
            ztemp+=slope*dx*TMath::Cos(phi);
            z=rnd->Gaus(ztemp,3);    //Smeared z

            phi+=slope*dx*curvature;
          }
          t1s.Fill();

        }
      t1s.Write();

      

      

      TFile f("m_perfect_helix_rndx_smz3_R_0001_3_stdK.root","recreate");
      TTree t1_FWD("t1_FWD","Forward fitter tree");
      float xht,yht,zht,xpost;
      TVectorF parvect(5);
      TVectorF predstept(5);
      TMatrixF Pt(5,5);
      TMatrixF PPredt(5,5);
      TMatrixF Rt(2,2);
     
      t1_FWD.Branch("xht",&xht,"xht/F");
      t1_FWD.Branch("yht",&yht,"yht/F");
      t1_FWD.Branch("zht",&zht,"zht/F");
      t1_FWD.Branch("xpost",&xpost,"xpost/F");
      t1_FWD.Branch("parvect",&parvect);
      t1_FWD.Branch("predstept",&predstept);
      t1_FWD.Branch("Pt",&Pt);
      t1_FWD.Branch("PPredt",&PPredt);
      t1_FWD.Branch("Rt",&Rt);

      KalmanFit(t1_FWD,xht,yht,zht,parvect,predstept,Pt,PPredt,Rt,xpost,n,t1s,x,y,z);
      
     
      t1_FWD.Write();

      
    }

    

    

    

  

    

    