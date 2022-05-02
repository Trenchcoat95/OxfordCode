



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
#include "correctMeanmaterial.h"
#include "helix_new.h"
#include <ctime>

using namespace ROOT::Math;




double capprox(double x1,double y1,
                        double x2,double y2,
                        double x3,double y3,
                        double &xc, double &yc)
{
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature -- copied from ALICE
  // here x is y and y is z for us
  //-----------------------------------------------------------------
  x3 -=x1;
  x2 -=x1;
  y3 -=y1;
  y2 -=y1;
  //
  double det = x3*y2-x2*y3;
  if (TMath::Abs(det)<1e-10){
    return 100;
  }
  //
  double u = 0.5* (x2*(x2-x3)+y2*(y2-y3))/det;
  double x0 = x3*0.5-y3*u;
  double y0 = y3*0.5+x3*u;
  double c2 = 1/TMath::Sqrt(x0*x0+y0*y0);
  xc = x0 + x1;
  yc = y0 + y1;
  if (det<0) c2*=-1;
  return c2;
}

void ouchef(double *x, double *y, int np, double &xcirccent, double &ycirccent,
                      double &rcirc, double &chisq, int &ifail)
{

  // converted from the original Fortran by T. Junk and de-OPALified -- apologies for the gotos

  //      SUBROUTINE OUCHEF(X,Y,NP,XCIRCCENT,YCIRCCENT,RCIRC,CHISQ,IFAIL)
  // *
  // *...OUCHEF  circle fit by linear regression method
  // *.
  // *.            OUCHEF(X,Y,NP,XCIRCCENT,YCIRCCENT,RCIRC,CHISQ,IFAIL)
  // *.
  // *.  INPUT
  // *.       X(.)      x coordinates of points
  // *.       Y(.)      y coordinates of points
  // *.       NP        number of points
  // *.
  // *.  OUTPUT
  // *.       XMED, YMED   center of circle
  // *.       R            radius of circle
  // *.       CHISQ     residual sum of squares : to be a chi**2 ,it should
  // *.                 be divided by (NP-3)*(mean residual)**2
  // *.
  // *.       IFAIL     error flag : 0   OK
  // *.                              1   less than 3 points
  // *.                              2   no convergence
  // *.                              5   number of points exceeding
  // *.                                  arrays dimensions -> fit performed
  // *.                                  with allowed maximum (NPTMAX)
  // *.
  // *. AUTHOR :  N.Chernov, G.Ososkov
  // *.
  // *.     written by N. Chernov and G. Ososkov, Dubna,
  // *.     reference:  computer physics communications vol 33, p329
  // *.     adapted slightly by F. James, March, 1986
  // *.     further adapted slightly by M.Hansroul May,1989
  // *.
  // *. CALLED   : any
  // *.
  // *. CREATED  : 31 May 1989
  // *. LAST MOD : 11 January 1996
  // *.
  // *. Modification Log.:
  // *. 11-jan-96  M.Hansroul  preset chisq to 9999.
  // *. 12-oct-95  M.Hansroul  clean up
  // *. 12-jul-93  M.Hansroul  protection against huge d0
  // *.                        put OUATG2 in line
  // *. 31-May-89  M.Hansroul  use OPAL output parameters
  // *.*********************************************************************

  const double eps = 1.0e-10;
  const int itmax=20;
  const int nptmax=1000;
  //const double twopi=8.0*atan(1.0);
  //const double piby2=2.0*atan(1.0);

  int nptlim=0;
  int npf=0;
  int n=0;
  int ihalf=0;
  int i=0;
  int iter=0;
  int iswap=0;
  //double sfi0=0;
  //double cfi0=0;
  double xmid=0;
  double ymid=0;
  double *xf;
  double *yf;
  double vlf[2];
  double vmf[2];
  //double vlm[2];

  double alf=0;
  double alm=0;
  double a0=0;
  double a1=0;
  double a2=0;
  double a22=0;
  double bem=0;
  double bet=0;
  double cur=0;
  double cu2=0;
  double dd=0;
  double den=0;
  double det=0;
  double dy=0;
  double d2=0;
  double f=0;
  double fact=0;
  double fg=0;
  double f1=0;
  double g=0;
  double gam=0;
  double gam0=0;
  double gmm=0;
  double g1=0;
  double h=0;
  double h2=0;
  double p2=0;
  double q2=0;
  double rm=0;
  double rn=0;
  double xa=0;
  double xb=0;
  double xd=0;
  double xi=0;
  double xm=0;
  double xx=0;
  double xy=0;
  double x1=0;
  double x2=0;
  double ya=0;
  double yb=0;
  double yd=0;
  double yi=0;
  double ym=0;
  double yy=0;
  double y1=0;
  double y2=0;
  double xcd=0;
  double ycd=0;
  //double xcycsq=0;
  //double rcd;

  //*     -----------------------------------------------------------------
  //*     xmid,ymid       = centre of fitted circle of radius radfit
  //*     chisq           = residual sum of squares per point
  //*     alf,bet,gam,cur = internal parameters of the circle

  ifail = 0;
  chisq = 9999.;

  nptlim = np;
  if (nptlim>nptmax)
    {
      nptlim = nptmax;
      ifail  = 5;
      return;
    }

  npf = nptlim + 1;
 label_12:
  npf = npf - 1;
  if (npf < 2)
    {
      ifail = 1;
      return;
    }
  vlf[0] = x[npf-1] - x[0];
  vlf[1] = y[npf-1] - y[0];

  //*                get a point near the middle

  ihalf  = (npf+1)/2;

  vmf[0] = x[ihalf-1] - x[0];
  vmf[1] = y[ihalf-1] - y[0];
  //vlm[0] = x[npf-1] - x[ihalf-1];
  //vlm[1] = y[npf-1] - y[ihalf-1];

  //*            check for more than half a cicle of points

  if (vlf[0]*vmf[0] + vlf[1]*vmf[1] < 0.) goto label_12;

  //*            if the track is nearer to the y axis
  //*            interchange x and y

  if (fabs(y[npf-1] - y[0]) > fabs(x[npf-1]-x[0]))
    {
      yf = x;
      xf = y;
      iswap = 1;
    }
  else
    {
      xf = x;
      yf = y;
      iswap = 0;
    }

  n  = npf;
  rn = 1./((double) n);
  xm = 0.;
  ym = 0.;
  for (i=0;i<n;++i)
    {
      xm = xm + xf[i];
      ym = ym + yf[i];
    }

  xm = xm * rn;
  ym = ym * rn;
  x2 = 0.;
  y2 = 0.;
  xy = 0.;
  xd = 0.;
  yd = 0.;
  d2 = 0.;
  for (i=0; i<n; ++i)
    {
      xi = xf[i] - xm;
      yi = yf[i] - ym;
      xx = xi*xi;
      yy = yi*yi;
      x2 = x2 + xx;
      y2 = y2 + yy;
      xy = xy + xi*yi;
      dd = xx + yy;
      xd = xd + xi*dd;
      yd = yd + yi*dd;
      d2 = d2 + dd*dd;
    }

  x2   = x2*rn;
  y2   = y2*rn;
  xy   = xy*rn;
  d2   = d2*rn;
  xd   = xd*rn;
  yd   = yd*rn;
  f    = 3.*x2 + y2;
  g    = 3.*y2 + x2;
  fg   = f*g;
  h    = xy + xy;
  h2   = h*h;
  p2   = xd*xd;
  q2   = yd*yd;
  gam0 = x2 + y2;
  fact = gam0*gam0;
  a2   = (fg-h2-d2)/fact;
  fact = fact*gam0;
  a1   = (d2*(f+g) - 2.*(p2+q2))/fact;
  fact = fact*gam0;
  a0   = (d2*(h2-fg) + 2.*(p2*g + q2*f) - 4.*xd*yd*h)/fact;
  a22  = a2 + a2;
  yb   = 1.0e30;
  iter = 0;
  xa   = 1.;
  //*                main iteration
 label_3:
  ya = a0 + xa*(a1 + xa*(a2 + xa*(xa-4.0)));
  if (fabs(ya) > fabs(yb))  goto label_4;
  if (iter >= itmax)  goto label_4;
  dy = a1 + xa*(a22 + xa*(4.*xa - 12.));
  xb = xa - ya/dy;
  if (fabs(xa-xb) <= eps)  goto label_5;
  xa = xb;
  yb = ya;
  iter = iter + 1;
  goto label_3;

 label_4:
  ifail = 2;
  //*      print 99, iter
  //*   99 format ('  circle - no convergence after',i4,' iterations.')
  xb = 1.;

 label_5:
  gam   = gam0*xb;
  f1    = f - gam;
  g1    = g - gam;
  x1    = xd*g1 - yd*h;
  y1    = yd*f1 - xd*h;
  det   = f1*g1 - h2;
  den   = 1./sqrt(x1*x1 + y1*y1 + gam*det*det);
  cur   = det*den;
  alf   = -(xm*det + x1)*den;
  bet   = -(ym*det + y1)*den;
  rm    = xm*xm + ym*ym;
  gam   = ((rm-gam)*det + 2.*(xm*x1 + ym*y1))*den*0.5;
  alm   = alf + cur*xm;
  bem   = bet + cur*ym;

  gmm   = gam + alf*xm + bet*ym + 0.5*cur*rm;
  chisq = ((0.5*cur)*(0.5*cur)*d2 + cur*(alm*xd + bem*yd + gmm*gam0)
           + alm*alm*x2 + bem*bem*y2 + 2.*alm*bem*xy + gmm*gmm) /rn;
  cu2   = 0.;
  if (cur != 0.)
    {
      cu2    =  1./cur;
      xcd    = -alf*cu2;
      ycd    = -bet*cu2;
      //rcd    =  fabs(cu2);
      xmid   =  xcd;
      ymid   =  ycd;
      xcirccent = xmid;
      ycirccent = ymid;
      if (iswap == 1)
        {
          xcirccent = ymid;
          ycirccent = xmid;
        }
      rcirc = cu2;
    }
  else
    {
      rcirc = 1e6;
    }

  return;
}


int Helix_Fit(const std::vector<XYZVector>  TPCClusters,
              XYZVector  &TPCClustersSeed,
              double &curvature_init,
              double &lambda_init,
              double &phi_init,
              double dir,
              int printlevel)
{

  size_t InitialTPNTPCClusters = 100;
  size_t nTPCClusters = TPCClusters.size();
  size_t firstTPCCluster = 0;
  size_t farTPCCluster = TMath::Min(nTPCClusters-1, InitialTPNTPCClusters);
  size_t intTPCCluster = farTPCCluster/2;
  size_t lastTPCCluster = nTPCClusters-1;
  

  double trackbeg[3] = {TPCClusters.at(firstTPCCluster).X(),
                       TPCClusters.at(firstTPCCluster).Y(),
                       TPCClusters.at(firstTPCCluster).Z()};

  double tp1[3] = {TPCClusters.at(intTPCCluster).X(),
                  TPCClusters.at(intTPCCluster).Y(),
                  TPCClusters.at(intTPCCluster).Z()};

  double tp2[3] = {TPCClusters.at(farTPCCluster).X(),
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
      std::cout << "First TPCCluster x, y, z: " << trackbeg[0] << " " << trackbeg[1] << " " << trackbeg[2] << std::endl;
      std::cout << "Inter TPCCluster x, y, z: " << tp1[0] << " " << tp1[1] << " " << tp1[2] << std::endl;
      std::cout << "Far   TPCCluster x, y, z: " << tp2[0] << " " << tp2[1] << " " << tp2[2] << std::endl;
    }

  TPCClustersSeed.SetXYZ(trackbeg[0],trackbeg[1],trackbeg[2]);
  double x_other_end = TPCClusters.at(lastTPCCluster).X();


  double ycc=0;
  double zcc=0;
  curvature_init = capprox(trackbeg[1],trackbeg[2],tp1[1],tp1[2],tp2[1],tp2[2],ycc,zcc);
  //std::cout << " inputs to trackpar circ fit (y,z): " << trackbeg[1] << " " << trackbeg[2] << " : "
  //        << tp1[1] << " " << tp1[2] << " : " << tp2[1] << " " << tp2[2] << std::endl;
  //std::cout << "curvature output: " << curvature_init << std::endl;

  // redo calc with a circle fit to all points, not just three

  double dycc=0;
  double dzcc=0;
  double drcc=0;
  double dchisq=0;
  int ifail=0;
  std::vector<double> dtpcclusy;
  std::vector<double> dtpcclusz;
  for (size_t i=0;i<nTPCClusters; ++i)
    {
      dtpcclusy.push_back(TPCClusters.at(i).Y());
      dtpcclusz.push_back(TPCClusters.at(i).Z());
    }
  int ict = nTPCClusters;
  ouchef(dtpcclusy.data(),dtpcclusz.data(),ict,dycc,dzcc,drcc,dchisq,ifail);
  if (ifail == 0 && drcc != 0)
    {
      ycc = dycc;
      zcc = dzcc;
      if (curvature_init < 0) drcc = -std::abs(drcc);
      //if (ict > 3)
      //  {
      //    std::cout << "updating curvature " << ict << " " << curvature_init << " with " << 1.0/drcc << std::endl;
      //  }
      curvature_init = 1.0/drcc;
    }


  phi_init = TMath::ATan2( trackbeg[2] - zcc, ycc - trackbeg[1] );
  double phi2 = phi_init;
  if (curvature_init<0) phi_init += TMath::Pi();
  double radius_init = 10000;
  if (curvature_init != 0) radius_init = 1.0/curvature_init;

  double dx1 = tp2[0] - TPCClustersSeed.X();
  if (dx1 != 0)
    {
      double dphi2 = TMath::ATan2(tp2[2]-zcc,ycc-tp2[1])-phi2;
      if (dphi2 > TMath::Pi()) dphi2 -= 2.0*TMath::Pi();
      if (dphi2 < -TMath::Pi()) dphi2 += 2.0*TMath::Pi();
      lambda_init = TMath::ATan(1.0/((radius_init/dx1)*dphi2));
    }
  else
    {
      //std::cout << "initial track par estimate failure" << std::endl;
      lambda_init = 0;
      return 1;
    } // got fMinNumTPCClusters all at exactly the same value of x (they were sorted).  Reject track.

  if (printlevel>0)
    {
      std::cout << "phi calc: dz, dy " << tp2[2]-trackbeg[2] << " " <<  tp2[1]-trackbeg[1] << std::endl;
      std::cout << "initial curvature, phi, lambda: " << curvature_init << " " << phi_init << " " << lambda_init << std::endl;
    }
  return 0;

}



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
                   std::vector<double> invpTplane,
                   std::string CorrTime,
                   Bool_t Fixed_Cov,
                   Bool_t Smear,
                   double xy_smear,
                   TMatrixD &P_seed,
                   std::string Seedtype,
                   std::string MS


                   

                                )
    {

      

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
      
      
      const Double_t kAlmost1=1. - Double_t(FLT_EPSILON);
      const Double_t kAlmost0=Double_t(FLT_MIN);
      const int nplanes = 6;
      double Planes_Z[]={1179,1183,1199,1203,1219,1223,1239,1243,1339,1343,1539,1543};
      if (dir<0) 
      {
        Planes_Z[0]=1539;  Planes_Z[1]=1543;
        Planes_Z[2]=1339;  Planes_Z[3]=1343;
        Planes_Z[4]=1239;  Planes_Z[5]=1243;
        Planes_Z[6]=1219;  Planes_Z[7]=1223;
        Planes_Z[8]=1199;  Planes_Z[9]=1203;
        Planes_Z[10]=1179; Planes_Z[11]=1183;
      }
        
      //for(size_t i=0;i<10;i++) std::cout<<Planes_Z[i]<<" ";
      //std::cout<<"\n";

      
      
      
      
      
      
      if(fPrintLevel>0)
      {
        std::cout << " Seed: y " << xyz_seed.Y() << " x " << xyz_seed.X() <<"z:"<<xyz_seed.Z()<< " sinphi " << sinphi_seed << " tanlambda " << tanlambda_seed << " 1/pT " << curvature_seed/(0.5*0.299792458e-2) << " p: " << (1/TMath::Cos(TMath::ATan(tanlambda_seed)))*(0.5*0.299792458e-2)/curvature_seed<<std::endl;
      }
      

      //std::cout<<"kalman: |p| "<<(1/TMath::Cos(TMath::ATan(tanlambda_seed)))*(0.5*0.299792458e-2)/curvature_seed<<" dxyz: not given "<< " 1/pT " << curvature_seed/(0.5*0.299792458e-2) << std::endl;

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
      //Q[4][4] = 0.0001;     // allow for some curvature uncertainty between points
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
      parvec[4] = curvature_seed/(0.5*0.299792458e-2);

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
              std::cout << "Adding a new TPCCluster: x:" << xh << " y: " << yh << " z: " << zh << " sinphi: " << sinphi_plane.at(iTPCCluster) << std::endl;
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
          if (abs(dz)<10)
          {
            for(size_t p=0;p<nplanes;p++)
            {
              if((zh>Planes_Z[p*2] && zh<Planes_Z[p*2+1]))
                  {

                    InPlane = 1;
                    /// do the prediction
                    if(planecounterk!=p && !dEvecreck.empty() && !dxvecreck.empty() && dir>0)
                    {
                      std::cout<<"Plane "<<p<<" dxreco: ";
                      for(int j =0; j<dxvecreck.size();j++) std::cout<<dxvecreck.at(j)<<" ";
                      std::cout<<"\n";
                      dEreco.push_back(dEvecreck);
                      dxreco.push_back(dxvecreck);
                      dEvecreck.clear();
                      dxvecreck.clear();
                    }
                    if(planecounterk!=p && dir>0) {planecounterk=p;}

                    
                    Double_t dxyzrec, dErec;
                    Bool_t checkstatus = Propagate(dz,PPred, P, predstep,  parvec, kAlmost0, kAlmost1, fPrintLevel, dir, InPlane, dErec, dxyzrec, Energy_loss, CorrTime, MS, xx0);
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
                    //std::cout<<dEvecreck.size()<<std::endl;
                      
                  }
            }
          }
          else
          {
            for(size_t p=0;p<nplanes;p++)
            {
              if((zh>Planes_Z[p*2] && zh<Planes_Z[p*2+1]))
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
                    if (dir>0) dzmid=Planes_Z[p*2-1]-zpos;
                    else dzmid=Planes_Z[p*2-2]-zpos;
                    //if(dir<0)std::cout<<"dz from point to plane: "<< dzmid <<std::endl;
                    InPlane=1;
                    Bool_t checkstatus = Propagate(dzmid,PPred, P, predstep,  parvec, kAlmost0, kAlmost1, fPrintLevel, dir, InPlane, dErec, dxyzrec, Energy_loss, CorrTime, MS, xx0);
                    dErectot+=dErec;
                    dxyzrectot+=dxyzrec;
                    if (!checkstatus)
                    { 
                      status=checkstatus;
                      break;
                    }
                    //if(dir>0)std::cout<<"dxyzmid: "<<dxyzrec;

                    InPlane = false;
                    if (dir>0) dzmid=Planes_Z[p*2]-Planes_Z[p*2-1];
                    else dzmid=Planes_Z[p*2+1]-Planes_Z[p*2-2];
                    TVectorD predmidstep = predstep;
                    TMatrixD Predmidstep = PPred;
                    //if(dir<0) std::cout<<"dz between planes: "<< dzmid <<std::endl;
                    checkstatus = Propagate(dzmid,PPred, Predmidstep, predstep,  predmidstep, kAlmost0, kAlmost1, fPrintLevel, dir, InPlane, dErec, dxyzrec, Energy_loss, CorrTime, MS, xx0);
                    if (!checkstatus)
                    { 
                      status=checkstatus;
                      break;
                    }

                    InPlane = true;
                    if (dir>0) dzmid=zh-Planes_Z[p*2];
                    else dzmid=zh-Planes_Z[2*p+1];
                    predmidstep = predstep;
                    Predmidstep = PPred;
                    //if(dir<0) std::cout<<"dz from plane to new point: "<< dzmid <<std::endl;
                    checkstatus = Propagate(dzmid,PPred, Predmidstep, predstep,  predmidstep, kAlmost0, kAlmost1, fPrintLevel, dir, InPlane, dErec, dxyzrec, Energy_loss,CorrTime, MS, xx0);
                    
                    if (!checkstatus)
                    { 
                      status=checkstatus;
                      break;
                    }
                    //if(dir>0)std::cout<<" "<<dxyzrec;
                    dErectot+=dErec;
                    dxyzrectot+=dxyzrec;
                    //if(dir>0)std::cout<<" dxyztot: "<<dxyzrectot<<std::endl;
                    
                    if(dir>0) dEvecreck.push_back(dErectot);
                    if(dir>0) dxvecreck.push_back(dxyzrectot);

                    dxyzrecord = dxyzrectot;
                    dErecord = dErectot;

                  }
            }

          }

          
          
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
          zpos = zpos + dz;
          
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
    
void kalman_helix_garlite_6planes(size_t nevents)

    {
      //srand ( time(NULL) );
      const Double_t kAlmost1=1. - Double_t(FLT_EPSILON);
      const Double_t kAlmost0=Double_t(FLT_MIN);
      // variables:  z is the independent variable
      // 0: y
      // 1: x
      // 2: sin(phi)
      // 3: tan(lambda)
      // 4: q/pT (1/(GeV/c))
      ///Prepare a simple tree with just the coordinates

      Double_t dEdx = 0.01;

      //Right now with perfect helix
      //std::cout << std::setprecision(17);
      

      TFile fs("./toygarlite/6planes/smearx05y05_aliceseed_elosslandaucorr_MScorr_fixedp0.7/garlitetest_smearx05y05_aliceseed_elosslandaucorr_MScorr_r05_05_fixedp0.7.root","recreate");
      //TFile fs("test.root","recreate");
      TTree t1s("t1s","helix simple tree");

      ///Parameters regulating simulation
      std::string Seedtype = "alice"; //perfect, real or alice
      std::string Helix_Corr = ""; //"Eloss" or "Eloss_MS"
      std::string Energy_Smear = "landau"; //gauss or landau
      std::string CorrTime = ""; //select if apply energy loss correction before or after a posteriori step. Either "after" or anything else for "before"
      Bool_t Backward_separate = false; // apply Kalman filter backwards reusing the Helix fit and not the final point in the forward Kalman
      Bool_t Fixed_Cov = true; // apply Kalman filter backwards using fixed guess values for the covariance matrix
      Bool_t Energy_loss_corr = true;
      Bool_t Energy_loss = true;
      Bool_t Smear = true;
      Bool_t Seed_only = false;
      std::string MS = "addMS_Smearing_Corr"; //use "addMS_Smearing" for just the multiple scattering smearing or "addMS_Smearing_Corr" to also have the correction
      double xy_smear = 0.5; //smear due to plane precision
      double  fixedp = 0.7; ///GeV/c Set to 0 or negative value if you don't want it fixed
      ////Kalman
      double Ry = TMath::Sq(0.5);//TMath::Sq(2);//TMath::Sq(0.01);//TMath::Sq(3);
      double Rx = TMath::Sq(0.5);//TMath::Sq(0.01);//TMath::Sq(0.01); //1.1921e-07
      double Ryx = TMath::Sq(0);

      
      double  Plane_XY_fid[4] = {-200, +200, -300, 0};

      const int  nplanes = 6 ;
      double  Planes_X[] =   {-300,300,-300,300,-300,300,-300,300,-300,300,-300,300};
      double  Planes_Y[] =   {-283.5,-16.5,-322.5,22.5,-350,50,-375,75,-400,100,-400,100};
      double  Planes_Z[] =   {1179,1183,1199,1203,1219,1223,1239,1243,1339,1343,1539,1543};
      XYZVector GArCenter(0,-150.473,1486); 
      double GAr_r = 349.9;
      double GAr_L = 669.6;
      double Plane_thick = 4;
      double  B=0.5; //in Tesla
      XYZVector Z_axis(0,0,1);
      XYZVector Y_axis(0,1,0);
      double ZA =0.54141;
      double I=64.7e-9;      ///GeV
      double rho= 1.032;   //g/cm^3
      double X1=2.49;
      double X0=0.1469;
      double muon_mass=0.1056583755; //GeV/c^2   ///material properties are for typical plastic polymere scintillator: polyvinyltoluene
      double a=0.1610;
      double m=3.24;
      double hw=21.54e-9;
      double Z=0.085+6*0.915;
      double xx0 = 42.54;  //radiation length in cm (different from Alice's code where it's xx0=1/X0 in cm^-1)

      int printlevelKalman = 0;
      int printlevelHelix = 0;
     
      
      
      ////MC

      

      std::vector<XYZVector> xyz;
      std::vector<XYZVector> pxyz;
      std::vector<XYZVector> xyz_plane;
      std::vector<XYZVector> pxyz_plane;
      std::vector<XYZVector> xyz_plane_sm;
      //std::vector<XYZVector> xyz_plane_mean;
      std::vector<double> sinphi;
      std::vector<double> sinphi_plane;
      std::vector<double> tanlambda;
      std::vector<double> tanlambda_plane;
      std::vector<double> invpT;  //this is q/pt
      std::vector<double> invpT_plane; //this is q/pt
      double  charge;
      int nPoints;
      Bool_t status;
      std::vector<int> nHits_perPlane;
      //std::vector<XYZVector> xyz_firstinPlane;
      //std::vector<XYZVector> pxyz_firstinPlane;
      //std::vector<double> sinphi_firstinPlane;
      //std::vector<double> tanlambda_firstinPlane;
      //std::vector<double> curvature_firstinPlane;
      

      

      ////FWD
      XYZVector xyz_seed;
      double sinphi_seed,tanlambda_seed,curvature_seed, phi_seed, lambdaGar_seed;
      TMatrixD P_seed(5,5);
      std::vector<double> xht,yht,zht,zpost;
      std::vector<TVectorD> parvect; //(5)
      std::vector<TVectorD> predstept; //(5)
      std::vector<TMatrixD> Pt;//(5,5);
      std::vector<TMatrixD> PPredt;//(5,5);
      std::vector<TMatrixD> Rt;//(2,2);

      ///BKW
      XYZVector xyz_seed_bkw;
      double sinphi_seed_bkw,tanlambda_seed_bkw,curvature_seed_bkw;
      TMatrixD P_seed_bkw(5,5);
      std::vector<double> xht_bkw,yht_bkw,zht_bkw,zpost_bkw;
      std::vector<TVectorD> parvect_bkw; //(5)
      std::vector<TVectorD> predstept_bkw; //(5)
      std::vector<TMatrixD> Pt_bkw;//(5,5);
      std::vector<TMatrixD> PPredt_bkw;//(5,5);
      std::vector<TMatrixD> Rt_bkw;//(2,2);

      std::vector<std::vector<double>> dEsim;
      std::vector<std::vector<double>> dxsim;
      std::vector<std::vector<double>> dEreco;
      std::vector<std::vector<double>> dxreco;
      std::vector<std::vector<double>> hitid;


      ///////MC
      t1s.Branch("status",&status);
      t1s.Branch("xyz",&xyz);
      t1s.Branch("pxyz",&pxyz);
      t1s.Branch("sinphi",&sinphi);
      t1s.Branch("sinphi_plane",&sinphi_plane);
      t1s.Branch("tanlambda",&tanlambda);
      t1s.Branch("tanlambda_plane",&tanlambda_plane);
      t1s.Branch("invpT",&invpT);
      t1s.Branch("invpT_plane",&invpT_plane);
      t1s.Branch("xyz_plane",&xyz_plane);
      t1s.Branch("pxyz_plane",&pxyz_plane);
      t1s.Branch("xyz_plane_sm",&xyz_plane_sm);
      //t1s.Branch("xyz_plane_mean",&xyz_plane_mean);
      t1s.Branch("charge",&charge);
      t1s.Branch("nPoints",&nPoints);
      t1s.Branch("nHits_perPlane",&nHits_perPlane);
      t1s.Branch("dEsim",&dEsim);
      t1s.Branch("dxsim",&dxsim);
      //t1s.Branch("xyz_firstinPlane",&xyz_firstinPlane);
      //t1s.Branch("sinphi_firstinPlane",&sinphi_firstinPlane);
      //t1s.Branch("tanlambda_firstinPlane",&tanlambda_firstinPlane);
      //t1s.Branch("pxyz_firstinPlane",&pxyz_firstinPlane);
      //t1s.Branch("curvature_firstinPlane",&curvature_firstinPlane);
      

      ////Kalman
      if(!Seed_only)
      {
        t1s.Branch("xht",&xht);
        t1s.Branch("yht",&yht);
        t1s.Branch("zht",&zht);
        t1s.Branch("zpost",&zpost);
        t1s.Branch("parvect",&parvect);
        t1s.Branch("predstept",&predstept);
        t1s.Branch("Pt",&Pt);
        t1s.Branch("PPredt",&PPredt);
        t1s.Branch("Rt",&Rt);
      }

      t1s.Branch("xyz_seed",&xyz_seed);
      t1s.Branch("sinphi_seed",&sinphi_seed);
      t1s.Branch("phi_seed",&phi_seed);
      t1s.Branch("lambdaGar_seed",&lambdaGar_seed);
      t1s.Branch("tanlambda_seed",&tanlambda_seed);
      t1s.Branch("curvature_seed",&curvature_seed);
      t1s.Branch("P_seed",&P_seed);
      t1s.Branch("dEreco",&dEreco);
      t1s.Branch("dxreco",&dxreco);
      t1s.Branch("hitid",&hitid);
      

      if(!Seed_only)
      {
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
        t1s.Branch("P_seed_bkw",&P_seed_bkw);
      }
      
     
      
      //size_t nevents=1;
      
      for(size_t i=0;i<nevents;i++)
      {
      XYZVector xyztemp(gRandom->Rndm()*(Plane_XY_fid[1]-Plane_XY_fid[0])+Plane_XY_fid[0],gRandom->Rndm()*(Plane_XY_fid[3]-Plane_XY_fid[2])+Plane_XY_fid[2],gRandom->Rndm()*10+(Planes_Z[0]-10));
      xyz.push_back(xyztemp);
      std::cout<<i<<std::endl;
      
      charge = (gRandom->Rndm()<0.5) ? -1 : 1;
      
      //std::cout<<i<<std::endl;
      
      //std::cout<<"check0"<<std::endl;

      double ptotal=0.1+gRandom->Rndm()*4; //in GeV/c
      if (fixedp>0) ptotal=fixedp;
      double px=-0.5+gRandom->Rndm();
      double py=-0.5+gRandom->Rndm();
      double pz=1.0+gRandom->Rndm();
      double norm=ptotal/sqrt(px*px+py*py+pz*pz);
      
      XYZVector pxyztemp(px*norm,py*norm,pz*norm);
      pxyz.push_back(pxyztemp);

      //std::cout<<"check0.1"<<std::endl;
      
      XYZVector pprojyz(0,pxyztemp.Y(),pxyztemp.Z()); ///pT
      //XYZVector pprojx(pxyztemp.X(),0,0);
      
      double sinphitemp = 0;
      if(pprojyz.Y()>0) sinphitemp=TMath::Sin(TMath::ACos(pprojyz.Unit().Dot(Z_axis.Unit())));
      else sinphitemp=-TMath::Sin(TMath::ACos(pprojyz.Unit().Dot(Z_axis.Unit())));
      //std::cout<<pxyztemp<<" "<<sinphitemp<<std::endl;
      double invpTtemp=charge/(sqrt(pxyztemp.Y()*pxyztemp.Y()+pxyztemp.Z()*pxyztemp.Z()));
      double curvaturetemp =((0.299792458e-2)*B)*invpTtemp; //in cm^-1
      
      double tanlambdatemp =0;
      if(pxyztemp.X()>0) tanlambdatemp =TMath::Tan(TMath::ACos(pxyztemp.Unit().Dot(pprojyz.Unit())));
      else tanlambdatemp =-TMath::Tan(TMath::ACos(pxyztemp.Unit().Dot(pprojyz.Unit())));

      double pyzmodule=sqrt(pxyztemp.Y()*pxyztemp.Y()+pxyztemp.Z()*pxyztemp.Z());
     
      //std::cout<<"check0.2"<<std::endl;

      Bool_t InGAr=true;
      Bool_t InPlane[nplanes]={0,0,0,0,0};
      double dz=0.500000000000;
      status=true;

      std::vector<XYZVector> xyzplanecontainer[nplanes];
      std::vector<XYZVector> pxyzplanecontainer[nplanes];
      std::vector<double> sinphiplanecontainer[nplanes];
      std::vector<double> tanlambdaplanecontainer[nplanes];
      std::vector<double> invpTplanecontainer[nplanes];

      //std::cout << " Initial coordinates: y " << xyztemp.Y() << " x " << xyztemp.X() << " sinphi " << sinphitemp << " tanlambda " << tanlambdatemp << " 1/pT " << charge/(sqrt(pprojyz.Y()*pprojyz.Y()+pprojyz.Z()*pprojyz.Z())) << std::endl;
     


      size_t planecounter = 10;
      std::vector<double> dEvecrec;
      std::vector<double> dxvecrec;

      while(InGAr)
        {

            
            
             //Randomized x
            //double dx=0.04;

            XYZVector xyzprev=xyztemp;
            double tanlambdaprev=tanlambdatemp;
            double curvatureprev=curvaturetemp;
            double sinphiprev=sinphitemp;

            Double_t z2r = curvaturetemp*dz;

            
            Double_t f1=sinphitemp;
            Double_t f2=f1 + z2r;
            
            
            if (TMath::Abs(f1) >= kAlmost1) 
            {
              std::cout<<"f1: "<<f1<<"y: "<<xyztemp.Y()<<"z: "<<xyztemp.Z()<<std::endl;
              status= false;
              break;
            }
            
            if (TMath::Abs(f2) >= kAlmost1) 
            {
              std::cout<<"InGAr f2: "<<f2<<std::endl;
              status= false;
              break;
            }
            if (TMath::Abs(tanlambdatemp)< kAlmost0) 
            {
              std::cout<<"tanlambdatemp: "<<tanlambdatemp<<std::endl;
              status= false;
              break;
            }

            Double_t r1=TMath::Sqrt((1.-f1)*(1.+f1)), r2=TMath::Sqrt((1.-f2)*(1.+f2));
            
            if (TMath::Abs(r1)<kAlmost0)  
            {
              std::cout<<"r1: "<<r1<<std::endl;
              status= false;
              break;
            }
            if (TMath::Abs(r2)<kAlmost0)  
            {
              std::cout<<"r2: "<<r1<<std::endl;
              status= false;
              break;
            }

            double dy2dz = (f1+f2)/(r1+r2);
            
            
            xyztemp.SetY(xyztemp.Y()+dz*dy2dz);
            xyztemp.SetZ(xyztemp.Z()+dz);

            
            sinphitemp+=z2r;
            
            //std::cout<<sqrt(pxyztemp.Y()*pxyztemp.Y()+pxyztemp.Z()*pxyztemp.Z())<<" "<<sinphitemp<<std::endl;

            //std::cout<<"checkGAr0"<<std::endl;

            double rot = TMath::ASin(r1*f2 - r2*f1); 
            if (f1*f1+f2*f2>1 && f1*f2<0) {          // special cases of large rotations or large abs angles
              if (f2>0) rot =  TMath::Pi() - rot;    //
              else      rot = -TMath::Pi() - rot;
            }
             

            //std::cout<<"checkGAr1"<<std::endl;
            xyztemp.SetX(xyztemp.X()+tanlambdatemp/curvaturetemp*rot);

            /////apply correction for energy loss
            
            for(size_t p=0;p<nplanes;p++)
            {
              
              if((xyztemp.X()>Planes_X[p*2] && xyztemp.X()<Planes_X[p*2+1] && xyztemp.Y()>Planes_Y[p*2] && xyztemp.Y()<Planes_Y[p*2+1]) ||
                 (xyzprev.X()>Planes_X[p*2] && xyzprev.X()<Planes_X[p*2+1] && xyzprev.Y()>Planes_Y[p*2] && xyzprev.Y()<Planes_Y[p*2+1])    )              
              {
                if((xyztemp.Z()>Planes_Z[p*2] && xyztemp.Z()<Planes_Z[p*2+1])||
                   (xyzprev.Z()>Planes_Z[p*2] && xyzprev.Z()<Planes_Z[p*2+1]))
                  {
                    //std::cout<<"In Plane "<<p<<std::endl;
                    
                    if(planecounter!=p && !dEvecrec.empty() && !dxvecrec.empty())
                    {
                      dEsim.push_back(dEvecrec);
                      dxsim.push_back(dxvecrec);
                      dEvecrec.clear();
                      dxvecrec.clear();
                    }
                    if(planecounter!=p) {planecounter=p;}
                    
                    
                    double deltaxyz;
                    XYZVector xyztemptemp=xyztemp;
                    if(xyzprev.Z()<Planes_Z[p*2])
                    {
                      double dzfromplane=Planes_Z[p*2]-xyzprev.Z();
                      xyzprev.SetZ(Planes_Z[p*2]);
                      Double_t z2rin = curvatureprev*dzfromplane;            
                      Double_t f1in=sinphiprev;
                      Double_t f2in=f1in + z2rin;      
                      Double_t r1in=TMath::Sqrt((1.-f1in)*(1.+f1in)), r2in=TMath::Sqrt((1.-f2in)*(1.+f2in));                
                      if (TMath::Abs(f1in) >= kAlmost1 || TMath::Abs(f2in) >= kAlmost1 || TMath::Abs(tanlambdaprev)< kAlmost0 || TMath::Abs(r1)<kAlmost0 || TMath::Abs(r2)<kAlmost0) 
                      {
                        std::cout<<"f1in: "<<f1in<<" f2in: "<<f2in<<" r1: "<<r1<<" r2: "<<r2<<" tantambdaprev: "<<tanlambdaprev<<std::endl;
                        status= false;
                        break;
                      }
                      double dy2dzin = (f1in+f2in)/(r1in+r2in);
                      double rotin = TMath::ASin(r1in*f2in - r2in*f1in); 
                      if (f1in*f1in+f2in*f2in>1 && f1in*f2in<0) {          // special cases of large rotations or large abs angles
                        if (f2in>0) rotin =  TMath::Pi() - rotin;    //
                        else      rotin = -TMath::Pi() - rotin;
                      }
                      xyzprev.SetY(xyzprev.Y()+dzfromplane*dy2dzin);
                      xyzprev.SetX(xyzprev.X()+tanlambdaprev/curvatureprev*rotin);
            
                    }
                    
                    if(xyztemptemp.Z()>Planes_Z[p*2+1])
                    {
                      double dzfromplane=Planes_Z[p*2+1]-xyzprev.Z();
                      //std::cout<<"dzfromplane: "<<dzfromplane<<" dz: "<<xyztemptemp.Z()-xyzprev.Z()<<std::endl;
                      xyztemptemp.SetZ(Planes_Z[p*2+1]);
                      Double_t z2rin = curvatureprev*dzfromplane;            
                      Double_t f1in=sinphiprev;
                      Double_t f2in=f1in + z2rin;      
                      Double_t r1in=TMath::Sqrt((1.-f1in)*(1.+f1in)), r2in=TMath::Sqrt((1.-f2in)*(1.+f2in));                
                      if (TMath::Abs(f1in) >= kAlmost1 || TMath::Abs(f2in) >= kAlmost1 || TMath::Abs(tanlambdaprev)< kAlmost0 || TMath::Abs(r1)<kAlmost0 || TMath::Abs(r2)<kAlmost0) 
                      {
                        std::cout<<"f1in: "<<f1in<<" f2in: "<<f2in<<" r1: "<<r1<<" r2: "<<r2<<" tantambdaprev: "<<tanlambdaprev<<std::endl;
                        status= false;
                        break;
                      }
                      double dy2dzin = (f1in+f2in)/(r1in+r2in);
                      double rotin = TMath::ASin(r1in*f2in - r2in*f1in); 
                      if (f1in*f1in+f2in*f2in>1 && f1in*f2in<0) {          // special cases of large rotations or large abs angles
                        if (f2in>0) rotin =  TMath::Pi() - rotin;    //
                        else      rotin = -TMath::Pi() - rotin;
                      }
                      xyztemptemp.SetY(xyzprev.Y()+dzfromplane*dy2dzin);
                      xyztemptemp.SetX(xyzprev.X()+tanlambdaprev/curvatureprev*rotin);
            
                    }
                    XYZVector Delta=(xyztemptemp-xyzprev);
                    deltaxyz=sqrt(Delta.Mag2());
                    //std::cout<<"deltaxyz propagation:"<<deltaxyz<<std::endl;
                    
                    //std::cout<<deltaxyz<<std::endl;
                    Bool_t checkstatus=1;
                    Double_t dErec=0;
                    if(Energy_loss) checkstatus= CorrectForMeanMaterial(-deltaxyz*rho,muon_mass,0.005,sqrt(pxyztemp.Mag2()),(sqrt(pxyztemp.Mag2())/muon_mass),
                                                                        rho,X0,X1,I,ZA,sinphitemp,tanlambdatemp,invpTtemp,dErec,Energy_Smear,MS,deltaxyz/xx0);

                   
                    dEvecrec.push_back(dErec);
                    dxvecrec.push_back(deltaxyz);
                   
                    
                    if (!checkstatus) status=checkstatus;

                    curvaturetemp =((0.299792458e-2)*B)*invpTtemp;

                    //std::cout<<"propagation: |p| "<<sqrt(pxyztemp.Mag2())<<" dxyz "<< deltaxyz << " 1/pT " << invpTtemp << " curvature " << std::endl;
                    pxyztemp.SetY(abs(1/invpTtemp)*sinphitemp);
                    pxyztemp.SetZ(abs(1/invpTtemp)*sqrt(1-sinphitemp*sinphitemp));
                    pxyztemp.SetX(abs(1/invpTtemp)*tanlambdatemp);
                    
                    
                    if (xyztemp.Z()>Planes_Z[p*2] && xyztemp.Z()<Planes_Z[p*2+1])
                    {
                      xyzplanecontainer[p].push_back(xyztemp);
                      pxyzplanecontainer[p].push_back(pxyztemp);
                      sinphiplanecontainer[p].push_back(sinphitemp);
                      tanlambdaplanecontainer[p].push_back(tanlambdatemp);
                      invpTplanecontainer[p].push_back(invpTtemp);

                      //std::cout<<abs((1/np.cos(np.arctan(t.tanlambda_plane.at(0))))*1/t.invpT_plane.at(0))
                    }
                    //double deltap = dEdx*
                    //std::cout<<status<<std::endl;
                    //std::cout<<"q/pT "<<invpTtemp<< " pxyz: "<<  abs(1/invpTtemp)*sinphitemp<<" "<<abs(1/invpTtemp)*sqrt(1-sinphitemp*sinphitemp)<<" "<<abs(1/invpTtemp)*tanlambdatemp<<std::endl;
                  }
              }
            }
            

            
            

            

            pxyztemp.SetY(abs(1/invpTtemp)*sinphitemp);
            pxyztemp.SetZ(abs(1/invpTtemp)*sqrt(1-sinphitemp*sinphitemp));
            pxyztemp.SetX(abs(1/invpTtemp)*tanlambdatemp);

            

            invpT.push_back(invpTtemp);
            xyz.push_back(xyztemp);
            pxyz.push_back(pxyztemp);
            sinphi.push_back(sinphitemp);
            tanlambda.push_back(tanlambdatemp);

            
            //std::cout<<"checkGAr3"<<std::endl;
            double r=sqrt(TMath::Power((xyztemp.Y()-GArCenter.Y()),2)+TMath::Power((xyztemp.Z()-GArCenter.Z()),2));
            double ry= abs(xyztemp.Y()-GArCenter.Y());
            double rz= abs(xyztemp.Z()-GArCenter.Z());
            //std::cout<<"checkGAr4"<<std::endl;
            if(xyztemp.X()>0.5*GAr_L || xyztemp.X()<-0.5*GAr_L || r>GAr_r || ry>GAr_r || rz>GAr_r) InGAr=false;
            //std::cout<<"checkGAr4.1"<<std::endl;

        }
        dEsim.push_back(dEvecrec);
        dxsim.push_back(dxvecrec);

      ////generate and smear hits in planes

      planecounter=10;
      std::vector<double> hitvecrec;
      
      for(size_t p=0;p<nplanes;p++)
      {

        //std::cout<<"check1"<<std::endl;
        //int rndm = gRandom->Rndm();
        int hitsinplane=0;
        
        //std::cout<<"xyzplanecontainer[p].size() "<<xyzplanecontainer[p].size()<<std::endl;
        if(!xyzplanecontainer[p].empty()) 
        {
        if(xyzplanecontainer[p].size()>1) {
          hitsinplane=2+gRandom->Rndm()*(xyzplanecontainer[p].size()-1);
          //std::cout<<hitsinplane<<" "<<xyzplanecontainer[p].size()<<std::endl;
          //hitsinplane= hitsinplane % (xyzplanecontainer[p].size()-1) + 2 ;
        }
        else hitsinplane=xyzplanecontainer[p].size();
        //hitsinplane=xyzplanecontainer[p].size();
        nHits_perPlane.push_back(hitsinplane);
        }
        else continue;

        if(planecounter!=p && !hitvecrec.empty() ) {
          hitid.push_back(hitvecrec);
          hitvecrec.clear();
        }
        if(planecounter!=p) planecounter=p;
        
        
        //std::cout<<"hitsinplane "<<hitsinplane<<std::endl;
        std::vector<size_t> rec;
        for(size_t l=0;l<xyzplanecontainer[p].size();l++) {
          rec.push_back(l);
          //std::cout<<rec.at(l)<<std::endl;
        }
        //std::random_shuffle(rec.begin(), rec.end());  
        //std::cout<<rec.at(0)<<std::endl;
        std::vector<size_t> recsub;           
        for(size_t l=0;l<hitsinplane;l++) {
          size_t position = gRandom->Rndm()*(rec.size()-1);
          recsub.push_back(rec.at(position));
          //std::cout<<recsub.at(l)<<" "<<rec.at(position)<<" "<<rec.size()<<std::endl;
          rec.erase(rec.begin()+position);
          //std::cout<<rec.size()<<std::endl<<std::endl;
        }
        std::sort(recsub.begin(),recsub.end());
        //std::cout<<"plane: "<<p<<" hitsinplane "<<hitsinplane<<" recsub.size "<<recsub.size()<<std::endl;
        //for(size_t l=0;l<hitsinplane;l++) std::cout<<recsub.at(l)<<" ";
        //std::cout<<std::endl;
        for(size_t l=0;l<hitsinplane;l++)
        {
          //std::cout<<recsub.at(l)<<" "<<std::endl;
          hitvecrec.push_back(recsub.at(l));
          xyz_plane.push_back(xyzplanecontainer[p].at(recsub.at(l)));
          pxyz_plane.push_back(pxyzplanecontainer[p].at(recsub.at(l)));
          sinphi_plane.push_back(sinphiplanecontainer[p].at(recsub.at(l)));
          tanlambda_plane.push_back(tanlambdaplanecontainer[p].at(recsub.at(l)));
          invpT_plane.push_back(invpTplanecontainer[p].at(recsub.at(l)));
          

          XYZVector xyz_plane_temp_sm;

          xyz_plane_temp_sm.SetX(gRandom->Gaus(xyzplanecontainer[p].at(recsub.at(l)).X(),xy_smear));
          xyz_plane_temp_sm.SetY(gRandom->Gaus(xyzplanecontainer[p].at(recsub.at(l)).Y(),xy_smear));
          xyz_plane_temp_sm.SetZ(xyzplanecontainer[p].at(recsub.at(l)).Z());

          //std::cout<<"X: "<<xyz_plane_temp_sm.X()<<" Y: "<<xyz_plane_temp_sm.Y()<<" Z: "<<xyz_plane_temp_sm.Z()<<std::endl;

          xyz_plane_sm.push_back(xyz_plane_temp_sm);
        }
      }
      //std::cout<<"\n";

      hitid.push_back(hitvecrec);
      
      //std::cout<<sqrt(pxyz_plane.at(0).Y()*pxyz_plane.at(0).Y()+pxyz_plane.at(0).Z()*pxyz_plane.at(0).Z())<<" "<<abs(1/invpT_plane.at(0))<<std::endl;
      //std::cout<<"check2"<<std::endl;
      
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

      size_t nCrossedPlanes = nHits_perPlane.size();

      
      if(xyz_plane.size()>0 && nHits_perPlane.size()>=3 && status==true) 
      {
        //std::cout<<"checkKal1"<<std::endl;
        parvect_bkw.clear();
        //////FWD

        ////Smear the plane points

        double forward=1.;
        
        if(Seedtype=="perfect")
          {
            //Perfect seed
            
            xyz_seed=xyz_plane.at(0);
            curvature_seed=invpT_plane.at(0)*((0.299792458e-2)*B);
            //std::cout<<"calculated q/pT: "<<curvature_seed/((0.299792458e-2)*B)<<" actual value "<<invpT_plane.at(0)<<std::endl; 
            sinphi_seed=sinphi_plane.at(0);
            tanlambda_seed=tanlambda_plane.at(0);  
          }
        else if(Seedtype=="real")
          {
            ///Helix Fit Seed

            if(Smear==kFALSE) Helix_Fit(xyz_plane,xyz_seed,curvature_seed,lambdaGar_seed,phi_seed,forward,printlevelHelix);
            else Helix_Fit(xyz_plane_sm,xyz_seed,curvature_seed,lambdaGar_seed,phi_seed,forward,printlevelHelix);
            sinphi_seed = TMath::Sin(phi_seed);
            tanlambda_seed = TMath::Tan(lambdaGar_seed);
            //std::cout<<"checkKal3"<<std::endl;
          }
        else if(Seedtype=="alice")
          {
            ///Helix Fit Seed as done in Alice

            //std::cout << "initial MC curvature, sinphi, tanlambda: " << invpT_plane.at(0)*((0.299792458e-2)*B) << " " << sinphi_plane.at(0) << " " << tanlambda_plane.at(0) << std::endl;
            if(Smear==kFALSE) makeSeed(xyz_plane,xyz_seed,curvature_seed,tanlambda_seed,sinphi_seed,forward,printlevelHelix,P_seed,0.0000001,Helix_Corr,nCrossedPlanes);
            else makeSeed(xyz_plane_sm,xyz_seed,curvature_seed,tanlambda_seed,sinphi_seed,forward,printlevelHelix,P_seed,xy_smear,Helix_Corr,nCrossedPlanes);
            //std::cout << "seed curvature, sinphi, tanlambda: " << curvature_seed << " " << sinphi_seed << " " << tanlambda_seed << std::endl;
            //std::cout << "P0 guess Matrix: " << std::endl;
            //P_seed.Print();
          }
        else
          {
            //No Seed
            xyz_seed.SetXYZ(0.,0.,0.);
            curvature_seed=0.0001;
            sinphi_seed=0;
            tanlambda_seed=0.01;
          }
        
        if(!Seed_only)
        {

          if(printlevelKalman>0) std::cout << " Real Values: y " << xyz_plane.at(0).Y() << " x " << xyz_plane.at(0).X() << " sinphi " << sinphi_plane.at(0) << " tanlambda " << tanlambda_plane.at(0) << " 1/pT " << invpT_plane.at(0) << " p: " <<sqrt(pxyz_plane.at(0).Mag2())<< std::endl;
          if(Smear==kFALSE) KalmanFit(xht,yht,zht,parvect,predstept,Pt,PPredt,Rt,zpost,xyz_plane,sinphi_plane,Ry,Rx,Ryx,xyz_seed,tanlambda_seed,curvature_seed,sinphi_seed,forward,status,printlevelKalman,dEreco,dxreco,Energy_loss_corr,invpT_plane,CorrTime,Fixed_Cov,Smear,xy_smear,P_seed,Seedtype,MS);
          else KalmanFit(xht,yht,zht,parvect,predstept,Pt,PPredt,Rt,zpost,xyz_plane_sm,sinphi_plane,Ry,Rx,Ryx,xyz_seed,tanlambda_seed,curvature_seed,sinphi_seed,forward,status,printlevelKalman,dEreco,dxreco,Energy_loss_corr,invpT_plane,CorrTime,Fixed_Cov,Smear,xy_smear,P_seed,Seedtype,MS);
          //std::cout<<"status= "<<status<<std::endl;
          //std::cout<<"checkKal4"<<std::endl;
          
          //////BKW
          //std::cout<<xht.size()<<" "<<xyz_plane.size()<<" "<<nHits_perPlane.size()<<std::endl;
          //std::cout<<"checkKal5"<<std::endl;
          if(!parvect.empty())
          {
            
            double backwards=-1.;
            std::vector<XYZVector> xyz_plane_bkw;
            if(Smear==kFALSE) xyz_plane_bkw=xyz_plane;
            else xyz_plane_bkw=xyz_plane_sm;

            std::reverse(xyz_plane_bkw.begin(),xyz_plane_bkw.end());

            if(Seedtype=="alice")
              {
                ///Helix Fit Seed as done in Alice

                //std::cout<<"initial curvature, sinphi, tanlambda Kalman seed backwards: "<<curvature_seed_bkw<<" "<<sinphi_seed_bkw<<" "<<tanlambda_seed_bkw<<std::endl;
                //std::cout << "initial MC curvature, sinphi, tanlambda: " << invpT_plane.at(0)*((0.299792458e-2)*B) << " " << sinphi_plane.at(0) << " " << tanlambda_plane.at(0) << std::endl;
                if(Smear==kFALSE) makeSeed(xyz_plane_bkw,xyz_seed_bkw,curvature_seed_bkw,tanlambda_seed_bkw,sinphi_seed_bkw,backwards,printlevelHelix,P_seed_bkw,0.0000001,Helix_Corr,nCrossedPlanes);
                else makeSeed(xyz_plane_bkw,xyz_seed_bkw,curvature_seed_bkw,tanlambda_seed_bkw,sinphi_seed_bkw,backwards,printlevelHelix,P_seed_bkw,xy_smear,Helix_Corr,nCrossedPlanes); 
              }
            xyz_seed_bkw.SetX(parvect.at(parvect.size()-1)[1]);
            xyz_seed_bkw.SetY(parvect.at(parvect.size()-1)[0]);
            //xyz_seed_bkw.SetZ(zht.at(xht.size()-1));
            xyz_seed_bkw.SetZ(xyz_plane.at(xyz_plane.size()-1).Z());
            curvature_seed_bkw=parvect.at(parvect.size()-1)[4]*(0.5*0.299792458e-2);
            sinphi_seed_bkw=parvect.at(parvect.size()-1)[2];
            tanlambda_seed_bkw=parvect.at(parvect.size()-1)[3];

            
            
            std::vector<double> invpT_plane_bkw;

            

            invpT_plane_bkw=invpT_plane;

            
            std::reverse(invpT_plane_bkw.begin(),invpT_plane_bkw.end());
            std::vector<double> sinphi_plane_bkw=sinphi_plane;
            std::reverse(sinphi_plane_bkw.begin(),sinphi_plane_bkw.end());
            Pt_bkw=Pt;
          
            

            if (Backward_separate)
            {
              //std::cout<<"initial xyz, curvature, sinphi, tanlambda Kalman seed backwards: "<<xyz_seed_bkw.X()<<" "<<xyz_seed_bkw.Y()<<" "<<xyz_seed_bkw.Z()<<" "<<curvature_seed_bkw<<" "<<sinphi_seed_bkw<<" "<<tanlambda_seed_bkw<<std::endl;
              //std::cout<<"initial curvature, phi, lambda Helix seed forward: "<<curvature_seed<<" "<<phi_seed<<" "<<lambdaGar_seed<<std::endl;
              Helix_Fit(xyz_plane_bkw,xyz_seed_bkw,curvature_seed_bkw,tanlambda_seed_bkw,sinphi_seed_bkw,backwards,printlevelHelix);
              curvature_seed_bkw=backwards*curvature_seed_bkw;
              sinphi_seed_bkw=backwards*TMath::Sin(sinphi_seed_bkw);
              tanlambda_seed_bkw=backwards*TMath::Tan(tanlambda_seed_bkw);
              //std::cout<<"initial xyz, curvature, sinphi, tanlambda Helix seed backwards: "<<xyz_seed_bkw.X()<<" "<<xyz_seed_bkw.Y()<<" "<<xyz_seed_bkw.Z()<<" "<<curvature_seed_bkw<<" "<<sinphi_seed_bkw<<" "<<tanlambda_seed_bkw<<std::endl;
            }
                  //std::cout<<xyz_seed_bkw.X()<<" "<<xyz_seed_bkw.X()<<" "<<xyz_seed_bkw.X()<<" "<<curvature_seed<<" "<<phi_seed_bkw<<" "<<lambda_seed_bkw<<std::endl;
                  //std::cout<<" "<<std::endl<<std::endl;
          
          

          KalmanFit(xht_bkw,yht_bkw,zht_bkw,parvect_bkw,predstept_bkw,Pt_bkw,PPredt_bkw,Rt_bkw,zpost_bkw,xyz_plane_bkw,sinphi_plane_bkw,Ry,Rx,Ryx,xyz_seed_bkw,tanlambda_seed_bkw,curvature_seed_bkw,sinphi_seed_bkw,backwards,status,printlevelKalman,dEreco,dxreco,Energy_loss_corr,invpT_plane_bkw,CorrTime,Fixed_Cov,Smear,xy_smear,P_seed_bkw,Seedtype,MS);
          //std::cout<<"status= "<<status<<std::endl;
          }
        }
        //std::cout<<"checkKal6"<<std::endl;
        
        //std::cout<<" "<<std::endl<<std::endl;
      }
      //std::cout<<"checkKal7"<<std::endl;

      //std::cout<<status<<std::endl;

      if((parvect.empty() || parvect_bkw.empty()) && !Seed_only) status=false;
      //std::cout<<parvect.size()<<std::endl;
      //std::cout<<"checkKal8"<<std::endl;
      //std::cout<<"status "<<status<<std::endl;
      
      if(status==true) t1s.Fill();

      /*
      //std::cout<<"nHits"<<xyz_plane.size()<<std::endl;
      if(status == true && abs((sinphi_plane.at(0)-parvect_bkw.at(parvect_bkw.size()-1)[2])/sinphi_plane.at(0))>1) {
        std::cout<<"sinMC: "<<sinphi_plane.at(0)<<" sinpredbkw: "<<parvect_bkw.at(parvect_bkw.size()-1)[2]<< "status:" <<status<<std::endl;
        //std::cout<<"sinMC: "<<sinphi_plane.at(0)<<" sinpred: "<<parvect.at(0)[2]<< "status:" <<status<<std::endl;
      }
      */
      //std::cout<<"check1"<<std::endl;

      //std::cout<<nHits_perPlane.size()<<" "<<dEsim.size()<<" "<<dEreco.size()<<" "<<hitid.size()<<" "<<xyz_plane.size()<<std::endl;
      
      if(!xyz.empty())xyz.clear();
      if(!pxyz.empty())pxyz.clear();
      if(!xyz_plane.empty())xyz_plane.clear();
      if(!pxyz_plane.empty())pxyz_plane.clear();
      if(!xyz_plane_sm.empty())xyz_plane_sm.clear();
      
      if(!xht.empty())xht.clear();
      if(!yht.empty())yht.clear();
      if(!zht.empty())zht.clear();
      if(!zpost.empty())zpost.clear();
      if(!parvect.empty())parvect.clear();
      if(!predstept.empty())predstept.clear();
      if(!Pt.empty())Pt.clear();//(5,5);
      if(!PPredt.empty())PPredt.clear();//(5,5);
      if(!Rt.empty())Rt.clear();//(2,2);
      
      if(!sinphi.empty())sinphi.clear();
      if(!tanlambda.empty())tanlambda.clear();
      if(!sinphi_plane.empty())sinphi_plane.clear();
      if(!tanlambda_plane.empty())tanlambda_plane.clear();
      if(!invpT_plane.empty())invpT_plane.clear();
      if(!invpT.empty())invpT.clear();

      

      for(int i =0; i<dEsim.size();i++)
      {
        //std::cout<<"Plane "<<i<<" dxsim: ";
        //for(int j =0; j<dEsim.at(i).size();j++) std::cout<<dEsim.at(i).at(j)<<" ";
        //std::cout<<"\n";
        //if(!dEsim.at(i).empty())dEsim.at(i).clear();
        //if(!dxsim.at(i).empty())dxsim.at(i).clear();

      }
      if(!dEsim.empty()) dEsim.clear();
      if(!dxsim.empty()) dxsim.clear();
      for(int i =0; i<dEreco.size();i++)
      {
        //std::cout<<"Plane "<<i<<" dxreco: ";
        //for(int j =0; j<dEreco.at(i).size();j++) std::cout<<dEreco.at(i).at(j)<<" ";
        //std::cout<<"\n";
        //if(!dEreco.at(i).empty())dEreco.at(i).clear();
        //if(!dxreco.at(i).empty())dxreco.at(i).clear();
      }
      if(!dEreco.empty()) dEreco.clear();
      if(!dxreco.empty()) dxreco.clear();
      
      for(int i =0; i<hitid.size();i++)
      {
        //std::cout<<"Plane "<<i<<" hitid: ";
        //for(int j =0; j<hitid.at(i).size();j++) std::cout<<hitid.at(i).at(j)<<" ";
        //std::cout<<"\n";
        //if(!hitid.at(i).empty())hitid.at(i).clear();
      }
      if(!hitid.empty()) hitid.clear();
      if(!nHits_perPlane.empty())nHits_perPlane.clear();

      //parvect_bkw.clear();

      

      
    }
    
    //std::cout<<t1s.GetEntries()<<std::endl;
    t1s.Write();
          
      
  }

    

    