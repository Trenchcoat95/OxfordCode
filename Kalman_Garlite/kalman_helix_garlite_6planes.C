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
#include "helix_old.h"
#include "kalman.h"
#include <ctime>
#include "material_part_utils.h"
using namespace ROOT::Math;
using namespace utils;


    
void kalman_helix_garlite_6planes(size_t nevents)

    {
      TFile fs("./toygarlite/6planes/smearx05y05_aliceseed_elossgaussnocorr_MSnocorr_fixedp0.7/garlitetest_smearx05y05_aliceseed_elossgaussnocorr_MSnocorr_r05_05_fixedp0.7.root","recreate");
      TTree t1s("t1s","helix simple tree");

      ///Parameters regulating simulation

      
      std::string Seedtype = "alice";       //perfect, real or alice
      std::string Helix_Corr = "";          //"Eloss" or "Eloss_MS"
      std::string Energy_Smear = "gauss";        //gauss or landau
      std::string CorrTime = "";            //select if apply energy loss correction before or after a posteriori step. 
                                            //Either "after" or anything else for "before"

      Bool_t Backward_separate = false;     // apply Kalman filter backwards reusing the Helix fit and not the final point in the forward Kalman
      Bool_t Fixed_Cov = true;              // apply Kalman filter backwards using fixed guess values for the covariance matrix
      Bool_t Energy_loss_corr = false;
      Bool_t Energy_loss = true;
      Bool_t Smear = true;
      Bool_t Seed_only = false;
      std::string MS = "addMS_Smearing";    //use "addMS_Smearing" for just the multiple scattering smearing 
                                            //or "addMS_Smearing_Corr" to also have the correction

      double xy_smear = 0.5;                //smear due to plane precision
      double  fixedp = 0.7;                 //GeV/c Set to 0 or negative value if you don't want it fixed
      ////Kalman
      double Ry = TMath::Sq(0.5);           //R matrix of Kalman Filter
      double Rx = TMath::Sq(0.5); 
      double Ryx = TMath::Sq(0);

      
      double  Plane_XY_fid[4] = {-200, +200, -300, 0};

      const int  nplanes = 6 ;
      double  Planes_X[] =   {-300,300,-300,300,-300,300,-300,300,-300,300,-300,300};
      double  Planes_Y[] =   {-283.5,-16.5,-322.5,22.5,-350,50,-375,75,-400,100,-400,100};
      double  Planes_Z[] =   {1179,1183,1199,1203,1219,1223,1239,1243,1339,1343,1539,1543};
   

      int printlevelKalman = 0;
      int printlevelHelix = 0;
     
      
      
      ////MC
      // variables:  z is the independent variable
      // 0: y
      // 1: x
      // 2: sin(phi)
      // 3: tan(lambda)
      // 4: q/pT (1/(GeV/c))
      ///Prepare a simple tree with just the coordinates
      

      std::vector<XYZVector> xyz;
      std::vector<XYZVector> pxyz;
      std::vector<XYZVector> xyz_plane;
      std::vector<XYZVector> pxyz_plane;
      std::vector<XYZVector> xyz_plane_sm;
      std::vector<double> sinphi;
      std::vector<double> sinphi_plane;
      std::vector<double> tanlambda;
      std::vector<double> tanlambda_plane;
      std::vector<double> invpT;                //this is q/pt
      std::vector<double> invpT_plane;          //this is q/pt
      double  charge;
      int nPoints;
      Bool_t status;
      std::vector<int> nHits_perPlane;
   
      

      

      ////FWD
      XYZVector xyz_seed;
      double sinphi_seed,tanlambda_seed,curvature_seed, phi_seed, lambdaGar_seed;
      TMatrixD P_seed(5,5);
      std::vector<double> xht,yht,zht,zpost;
      std::vector<TVectorD> parvect;            //(5)
      std::vector<TVectorD> predstept;          //(5)
      std::vector<TMatrixD> Pt;                 //(5,5);
      std::vector<TMatrixD> PPredt;             //(5,5);
      std::vector<TMatrixD> Rt;                 //(2,2);

      ///BKW
      XYZVector xyz_seed_bkw;
      double sinphi_seed_bkw,tanlambda_seed_bkw,curvature_seed_bkw;
      TMatrixD P_seed_bkw(5,5);
      std::vector<double> xht_bkw,yht_bkw,zht_bkw,zpost_bkw;
      std::vector<TVectorD> parvect_bkw;        //(5)
      std::vector<TVectorD> predstept_bkw;      //(5)
      std::vector<TMatrixD> Pt_bkw;             //(5,5);
      std::vector<TMatrixD> PPredt_bkw;         //(5,5);
      std::vector<TMatrixD> Rt_bkw;             //(2,2);

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
      t1s.Branch("charge",&charge);
      t1s.Branch("nPoints",&nPoints);
      t1s.Branch("nHits_perPlane",&nHits_perPlane);
      t1s.Branch("dEsim",&dEsim);
      t1s.Branch("dxsim",&dxsim);


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
      
     

      
      for(size_t i=0;i<nevents;i++)
      {

      std::cout<<i<<std::endl;
      
      ////Randomize initial position  

      XYZVector xyztemp(gRandom->Rndm()*(Plane_XY_fid[1]-Plane_XY_fid[0])+Plane_XY_fid[0]
                                        ,gRandom->Rndm()*(Plane_XY_fid[3]-Plane_XY_fid[2])+Plane_XY_fid[2]
                                        ,gRandom->Rndm()*10+(Planes_Z[0]-10));
      xyz.push_back(xyztemp);
     
      ////Randomize charge  

      charge = (gRandom->Rndm()<0.5) ? -1 : 1;


      ////Randomize total momentum between 0.1 and 4 GeV/c

      double ptotal=0.1+gRandom->Rndm()*4; //in GeV/c
      if (fixedp>0) ptotal=fixedp;

      ////Randomize angle

      double px=-0.5+gRandom->Rndm();
      double py=-0.5+gRandom->Rndm();
      double pz=1.0+gRandom->Rndm();
      double norm=ptotal/sqrt(px*px+py*py+pz*pz);      
      XYZVector pxyztemp(px*norm,py*norm,pz*norm);
      pxyz.push_back(pxyztemp);

      
      ///Get helix variable accordingly

      XYZVector pprojyz(0,pxyztemp.Y(),pxyztemp.Z()); ///pT     
      double sinphitemp = 0;
      if(pprojyz.Y()>0) sinphitemp=TMath::Sin(TMath::ACos(pprojyz.Unit().Dot(Z_axis.Unit())));
      else sinphitemp=-TMath::Sin(TMath::ACos(pprojyz.Unit().Dot(Z_axis.Unit())));
      double invpTtemp=charge/(sqrt(pxyztemp.Y()*pxyztemp.Y()+pxyztemp.Z()*pxyztemp.Z()));
      double curvaturetemp =((0.299792458e-2)*B)*invpTtemp; //in cm^-1     
      double tanlambdatemp =0;
      if(pxyztemp.X()>0) tanlambdatemp =TMath::Tan(TMath::ACos(pxyztemp.Unit().Dot(pprojyz.Unit())));
      else tanlambdatemp =-TMath::Tan(TMath::ACos(pxyztemp.Unit().Dot(pprojyz.Unit())));
      double pyzmodule=sqrt(pxyztemp.Y()*pxyztemp.Y()+pxyztemp.Z()*pxyztemp.Z());
     

      Bool_t InGAr=true;
      Bool_t InPlane[nplanes]={0,0,0,0,0};
      double dz=0.500000000000;
      status=true;

      std::vector<XYZVector> xyzplanecontainer[nplanes];
      std::vector<XYZVector> pxyzplanecontainer[nplanes];
      std::vector<double> sinphiplanecontainer[nplanes];
      std::vector<double> tanlambdaplanecontainer[nplanes];
      std::vector<double> invpTplanecontainer[nplanes];

      size_t planecounter = 10;
      std::vector<double> dEvecrec;
      std::vector<double> dxvecrec;

      ///Propagate particle until it leaves tracking area

      while(InGAr)
        {

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
            
            double rot = TMath::ASin(r1*f2 - r2*f1); 
            if (f1*f1+f2*f2>1 && f1*f2<0) {          // special cases of large rotations or large abs angles
              if (f2>0) rot =  TMath::Pi() - rot;    //
              else      rot = -TMath::Pi() - rot;
            }
             

            xyztemp.SetX(xyztemp.X()+tanlambdatemp/curvaturetemp*rot);

            /////apply correction for energy loss and MS
            
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

                    Bool_t checkstatus=1;
                    Double_t dErec=0;
                    if(Energy_loss) checkstatus= CorrectForMeanMaterial(-deltaxyz*rho,muon_mass,0.005,sqrt(pxyztemp.Mag2()),(sqrt(pxyztemp.Mag2())/muon_mass),
                                                                        rho,X0,X1,Ipar,ZA,sinphitemp,tanlambdatemp,invpTtemp,dErec,Energy_Smear,MS,deltaxyz/xx0);

                   
                    dEvecrec.push_back(dErec);
                    dxvecrec.push_back(deltaxyz);
                   
                    
                    if (!checkstatus) status=checkstatus;

                    curvaturetemp =((0.299792458e-2)*B)*invpTtemp;

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
                    }
    
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

            

            double r=sqrt(TMath::Power((xyztemp.Y()-GArCenter.Y()),2)+TMath::Power((xyztemp.Z()-GArCenter.Z()),2));
            double ry= abs(xyztemp.Y()-GArCenter.Y());
            double rz= abs(xyztemp.Z()-GArCenter.Z());
            if(xyztemp.X()>0.5*GAr_L || xyztemp.X()<-0.5*GAr_L || r>GAr_r || ry>GAr_r || rz>GAr_r) InGAr=false;

        }
        dEsim.push_back(dEvecrec);
        dxsim.push_back(dxvecrec);

      ////generate and smear hits in planes

      planecounter=10;
      std::vector<double> hitvecrec;
      
      for(size_t p=0;p<nplanes;p++)
      {
        int hitsinplane=0;
        
        if(!xyzplanecontainer[p].empty()) 
        {
        if(xyzplanecontainer[p].size()>1) {
          hitsinplane=2+gRandom->Rndm()*(xyzplanecontainer[p].size()-1);
        }
        else hitsinplane=xyzplanecontainer[p].size();
        nHits_perPlane.push_back(hitsinplane);
        }
        else continue;

        if(planecounter!=p && !hitvecrec.empty() ) {
          hitid.push_back(hitvecrec);
          hitvecrec.clear();
        }
        if(planecounter!=p) planecounter=p;
        
        
        std::vector<size_t> rec;
        for(size_t l=0;l<xyzplanecontainer[p].size();l++) {
          rec.push_back(l);
        }

        std::vector<size_t> recsub;           
        for(size_t l=0;l<hitsinplane;l++) {
          size_t position = gRandom->Rndm()*(rec.size()-1);
          recsub.push_back(rec.at(position));
          rec.erase(rec.begin()+position);
        }
        std::sort(recsub.begin(),recsub.end());

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



          xyz_plane_sm.push_back(xyz_plane_temp_sm);
        }
      }


      hitid.push_back(hitvecrec);
      
      nPoints=xyz.size();

      size_t nCrossedPlanes = nHits_perPlane.size();

      
      if(xyz_plane.size()>0 && nHits_perPlane.size()>=3 && status==true) 
      {
        parvect_bkw.clear();
        //////FWD


        double forward=1.;
        
        if(Seedtype=="perfect")
          {
            //Perfect seed          
            xyz_seed=xyz_plane.at(0);
            curvature_seed=invpT_plane.at(0)*((0.299792458e-2)*B);
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
          }
        else if(Seedtype=="alice")
          {
            ///Helix Fit Seed as done in Alice
            if(Smear==kFALSE) makeSeed(xyz_plane,xyz_seed,curvature_seed,tanlambda_seed,sinphi_seed,forward,printlevelHelix,P_seed,0.0000001,Helix_Corr,nCrossedPlanes);
            else makeSeed(xyz_plane_sm,xyz_seed,curvature_seed,tanlambda_seed,sinphi_seed,forward,printlevelHelix,P_seed,xy_smear,Helix_Corr,nCrossedPlanes);
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
          if(Smear==kFALSE) KalmanFit(xht,yht,zht,parvect,predstept,Pt,PPredt,Rt,zpost,xyz_plane,Ry,Rx,Ryx,xyz_seed,tanlambda_seed,curvature_seed,sinphi_seed,forward,status,printlevelKalman,dEreco,dxreco,Energy_loss_corr,CorrTime,Fixed_Cov,Smear,xy_smear,P_seed,Seedtype,MS);
          else KalmanFit(xht,yht,zht,parvect,predstept,Pt,PPredt,Rt,zpost,xyz_plane_sm,Ry,Rx,Ryx,xyz_seed,tanlambda_seed,curvature_seed,sinphi_seed,forward,status,printlevelKalman,dEreco,dxreco,Energy_loss_corr,CorrTime,Fixed_Cov,Smear,xy_smear,P_seed,Seedtype,MS);
          
          //////BKW

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

                if(Smear==kFALSE) makeSeed(xyz_plane_bkw,xyz_seed_bkw,curvature_seed_bkw,tanlambda_seed_bkw,sinphi_seed_bkw,backwards,printlevelHelix,P_seed_bkw,0.0000001,Helix_Corr,nCrossedPlanes);
                else makeSeed(xyz_plane_bkw,xyz_seed_bkw,curvature_seed_bkw,tanlambda_seed_bkw,sinphi_seed_bkw,backwards,printlevelHelix,P_seed_bkw,xy_smear,Helix_Corr,nCrossedPlanes); 
              }

            xyz_seed_bkw.SetX(parvect.at(parvect.size()-1)[1]);
            xyz_seed_bkw.SetY(parvect.at(parvect.size()-1)[0]);
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
              Helix_Fit(xyz_plane_bkw,xyz_seed_bkw,curvature_seed_bkw,tanlambda_seed_bkw,sinphi_seed_bkw,backwards,printlevelHelix);
              curvature_seed_bkw=backwards*curvature_seed_bkw;
              sinphi_seed_bkw=backwards*TMath::Sin(sinphi_seed_bkw);
              tanlambda_seed_bkw=backwards*TMath::Tan(tanlambda_seed_bkw);
            }
          
          KalmanFit(xht_bkw,yht_bkw,zht_bkw,parvect_bkw,predstept_bkw,Pt_bkw,PPredt_bkw,Rt_bkw,zpost_bkw,xyz_plane_bkw,Ry,Rx,Ryx,xyz_seed_bkw,tanlambda_seed_bkw,curvature_seed_bkw,sinphi_seed_bkw,backwards,status,printlevelKalman,dEreco,dxreco,Energy_loss_corr,CorrTime,Fixed_Cov,Smear,xy_smear,P_seed_bkw,Seedtype,MS);
          }
        }
      }

      if((parvect.empty() || parvect_bkw.empty()) && !Seed_only) status=false;
      if(status==true) t1s.Fill();

      
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

      
      if(!dEsim.empty()) dEsim.clear();
      if(!dxsim.empty()) dxsim.clear();
      if(!dEreco.empty()) dEreco.clear();
      if(!dxreco.empty()) dxreco.clear();
      if(!hitid.empty()) hitid.clear();
      if(!nHits_perPlane.empty())nHits_perPlane.clear();
    }
    
    t1s.Write();
          
      
  }  