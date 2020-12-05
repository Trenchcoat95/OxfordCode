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
                   int &ev, 
                   TVectorF &parvect,
                   TVectorF &predstept,
                   TMatrixF &Pt,
                   TMatrixF &PPredt,
                   TMatrixF &Rt,
                   float &xpost
                                )
    {

      
      size_t n = 452;
      float fPrintLevel=0;
      
      // estimate curvature, lambda, phi, xpos from the initial track parameters
      float curvature_init= 0.00133671;         ///need to check what this is
      float phi_init = 2.41558;
      float lambda_init = 0.0373353;
      float xpos_init= -23.0761;
      float ypos_init=-340.074;
      float zpos_init=1638.21;
      
      

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
      R[0][0] = TMath::Sq(4);  // in cm^2
      R[1][1] = TMath::Sq(4);  // in cm^2
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
      ev=0;
      Rt=R;
      xpost=xpos;
      t.Fill();

      

      TMatrixF H(2,5);   // partial(obs)/partial(params)
      H.Zero();
      H[0][0] = 1;  // y
      H[1][1] = 1;  // z
      TMatrixF HT(5,2);

      TVectorF z(2);
      TVectorF ytilde(2);
      TVectorF hx(2);
      TMatrixF S(2,2);
      TMatrixF K(5,2);

      TMatrixF I(5,5);
      I.Zero();
      for (int i=0;i<5;++i) I[i][i] = 1;


      
      //Create a parametrized helix as a substitute for the measurement
      float xh = -23.3112;  ///need to check what this is
      float yh = -337.045;
      float zh = 1634.73;
      float phih = 6;
      float curvatureh =-0.014;
      float lambdah =-0.05;
      float slopeh = TMath::Tan(lambdah);
      if (slopeh != 0)
            {
              slopeh = 1.0/slopeh;
            }
          else
            {
              slopeh = 1E9;
            }
      
      for (size_t iTPCCluster=1; iTPCCluster<n; ++iTPCCluster)
        {
          ev=iTPCCluster;
          //std::cout<<"iTPCCluster "<<iTPCCluster<<std::endl;
          if (iTPCCluster>1)
          {
            xh+=0.04;
            yh+=slopeh*0.04*TMath::Sin(phih);
            zh+=slopeh*0.04*TMath::Cos(phih);
            phih+=slopeh*0.04*curvatureh;

          }
          
          
          
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

          
          //dx = xh - xpos;
          //if (dx == 0) dx = 1E-3;
          
          
            float fTPCClusterResolYZ=1;
            float fTPCClusterResolX=0.5;

            float dxnum = (slope/(fTPCClusterResolYZ*fTPCClusterResolYZ))*( (yh - parvec[0])*TMath::Sin(phi) + (zh - parvec[1])*TMath::Cos(phi) )
            + (xh - xpos)/(fTPCClusterResolX*fTPCClusterResolX);
            float dxdenom = slope*slope/(fTPCClusterResolYZ*fTPCClusterResolYZ) + 1.0/(fTPCClusterResolX*fTPCClusterResolX);
            dx = dxnum/dxdenom;
            if (dx == 0) dx = 1E-3;
          
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
    
    void kalman_helix()

    {
      // variables:  x is the independent variable
      // 0: y
      // 1: z
      // 2: curvature
      // 3: phi
      // 4: lambda = angle from the cathode plane
      // 5: x   /// added on to the end


      

      TFile f("rootmacro_perfect_helix.root","recreate");
      TTree t1_FWD("t1_FWD","Forward fitter tree");
      float xht,yht,zht,xpost;
      int ev;
      TVectorF parvect(5);
      TVectorF predstept(5);
      TMatrixF Pt(5,5);
      TMatrixF PPredt(5,5);
      TMatrixF Rt(2,2);
     
      t1_FWD.Branch("xht",&xht,"xht/F");
      t1_FWD.Branch("yht",&yht,"yht/F");
      t1_FWD.Branch("zht",&zht,"zht/F");
      t1_FWD.Branch("ev",&ev,"ev/I");
      t1_FWD.Branch("xpost",&xpost,"xpost/F");
      t1_FWD.Branch("parvect",&parvect);
      t1_FWD.Branch("predstept",&predstept);
      t1_FWD.Branch("Pt",&Pt);
      t1_FWD.Branch("PPredt",&PPredt);
      t1_FWD.Branch("Rt",&Rt);

      KalmanFit(t1_FWD,xht,yht,zht,ev,parvect,predstept,Pt,PPredt,Rt,xpost);
      
     
      t1_FWD.Write();

      
    }

    

    