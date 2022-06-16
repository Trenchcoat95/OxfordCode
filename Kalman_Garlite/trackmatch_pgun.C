#define garana_cxx
#include <iostream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TString.h>
#include <TVector3.h>
#include <TH2F.h>
#include <TMath.h>
#include <TROOT.h>
#include <TPad.h>
#include <TEfficiency.h>
#include <TTree.h>
#include <TVectorF.h>
#include <TMatrix.h>
#include <TChain.h>
#include <TF1.h>
#include <TF1NormSum.h>
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
#include "kalman_new.h"
#include <ctime>
#include "material_part_utils.h"
#include "garana.h"
using namespace ROOT::Math;
using namespace utils;


//   In a ROOT session, you can do:
  //      root> .L l2g_trackmatch_t.C
  //      root> garana t
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.l2g_Trackmatch();       // Use the trackmatching function
//int t=0;
void perm(std::vector<vector<XYZVector>> &v,std::vector<XYZVector> vtot,std::vector<int> a, int n, int k, int i)
{

    
    if(i == 0)
    {
        //cout<< "permutation "<<t<<" "<<std::endl;
        //for(int j=n; j<n+k; j++) cout << a[j] << " ";
        //cout << endl;
        std::vector<XYZVector> temp;
        for(int j=n; j<n+k; j++) temp.push_back(vtot[a[j]]) ;
        //if(a[n]==1&&a[n+1]==0&&a[n+2]==2) std::cout<<"Gotcha! "<<t<<std::endl;
        v.push_back(temp);
        
        //t+=1;
        return;
    }

    for(int j=0; j<n; j++)
    {
        std::swap(a[j], a[n-1]);
        perm(v,vtot,a, n-1, k, i-1);
        std::swap(a[j], a[n-1]);
    }

}


void trackmatch_pgun(size_t start_entry, size_t n_entries)
{
    ////Geometry parameters
    TVector3 GArCenter(0,-150.473,1486); 
    float GAr_r = 349.9;
    float GAr_L = 669.6;
    const int  nplanes = 6 ;
    double  Planes_X[] =   {-300,300,-300,300,-300,300,-300,300,-300,300,-300,300};
    double  Planes_Y[] =   {-283.5,-16.5,-322.5,22.5,-350,50,-375,75,-400,100,-400,100};
    double  Planes_Z[] =   {1179,1183,1199,1203,1219,1223,1239,1243,1339,1343,1539,1543};
    std::cout << std::setprecision(10);

    /////Kalman parameters
    std::string Seedtype = "real";       //perfect, real or alice
    Bool_t OldCov = true;
    std::string Helix_Corr = "Eloss_MS";          //"Eloss" or "Eloss_MS"
    std::string Energy_Smear = "";        //gauss or landau
    std::string CorrTime = "";            //select if apply energy loss correction before or after a posteriori step. 
                                        //Either "after" or anything else for "before"
    Bool_t Backward_separate = false;     // apply Kalman filter backwards reusing the Helix fit and not the final point in the forward Kalman
    Bool_t Fixed_Cov = true;              // apply Kalman filter backwards using fixed guess values for the covariance matrix
    Bool_t Energy_loss_corr = true;
    Bool_t Smear = true;
    Bool_t Seed_only = false;
    std::string MS = "addMS_Smearing_Corr";    //use "addMS_Smearing" for just the multiple scattering smearing 

    double Reuse_mat=1.0;
    double Reuse_mat_cycle=0.0;
    int    ncycles=1;
                                            //or "addMS_Smearing_Corr" to also have the correction   
    double xy_smear = 0.3;                //smear due to plane precision                                    
    double Ry = TMath::Sq(0.3);           //R matrix of Kalman Filter
    double Rx = TMath::Sq(0.3); 
    double Ryx = TMath::Sq(0);

    int printlevelKalman = 0;
    int printlevelHelix = 0;

    //////Quantities to be stored
    TFile fs("./MCgarlite/6planes/muon_test_ana_kalman.root","recreate");
    TTree t1s("t1s","helix simple tree");

    //////////////////////////////////Read from MC and gar reco
    Bool_t status=1;
    XYZVector xyz_MC;
    XYZVector pxyz_MC;
    double charge_MC;
    double sinphi_MC,tanlambda_MC,curvature_MC,invpT_MC;

    XYZVector xyz_seed_old;
    XYZVector pxyz_seed_old;
    double charge_seed_old;
    double sinphi_seed_old,tanlambda_seed_old,curvature_seed_old,invpT_seed_old;

    XYZVector xyz_seed_old_bkw;
    XYZVector pxyz_seed_old_bkw;
    double charge_seed_old_bkw;
    double sinphi_seed_old_bkw,tanlambda_seed_old_bkw,curvature_seed_old_bkw,invpT_seed_old_bkw;

    std::vector<XYZVector> xyz_plane;
    std::vector<XYZVector> xyz_plane_as_is;

    //////////////////////////////////Produced with custom Reco

    XYZVector xyz_seed;
    XYZVector pxyz_seed;
    double charge_seed;
    double sinphi_seed,tanlambda_seed,curvature_seed,invpT_seed;
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

    std::vector<std::vector<double>> dEreco;
    std::vector<std::vector<double>> dxreco;

    


    t1s.Branch("xyz_MC",&xyz_MC);
    t1s.Branch("pxyz_MC",&pxyz_MC);
    t1s.Branch("charge_MC",&charge_MC);
    t1s.Branch("sinphi_MC",&sinphi_MC);
    t1s.Branch("tanlambda_MC",&tanlambda_MC);
    t1s.Branch("curvature_MC",&curvature_MC);
    t1s.Branch("invpT_MC",&invpT_MC);

    t1s.Branch("xyz_seed_old",&xyz_seed_old);
    t1s.Branch("pxyz_seed_old",&pxyz_seed_old);
    t1s.Branch("charge_seed_old",&charge_seed_old);
    t1s.Branch("sinphi_seed_old",&sinphi_seed_old);
    t1s.Branch("tanlambda_seed_old",&tanlambda_seed_old);
    t1s.Branch("curvature_seed_old",&curvature_seed_old);
    t1s.Branch("invpT_seed_old",&invpT_seed_old);

    t1s.Branch("xyz_seed_old_bkw",&xyz_seed_old_bkw);
    t1s.Branch("pxyz_seed_old_bkw",&pxyz_seed_old_bkw);
    t1s.Branch("charge_seed_old_bkw",&charge_seed_old_bkw);
    t1s.Branch("sinphi_seed_old_bkw",&sinphi_seed_old_bkw);
    t1s.Branch("tanlambda_seed_old_bkw",&tanlambda_seed_old_bkw);
    t1s.Branch("curvature_seed_old_bkw",&curvature_seed_old_bkw);
    t1s.Branch("invpT_seed_old_bkw",&invpT_seed_old_bkw);

    t1s.Branch("xyz_plane",&xyz_plane);

    t1s.Branch("xyz_seed",&xyz_seed);
    t1s.Branch("pxyz_seed",&pxyz_seed);
    t1s.Branch("charge_seed",&charge_seed);
    t1s.Branch("sinphi_seed",&sinphi_seed);
    t1s.Branch("tanlambda_seed",&tanlambda_seed);
    t1s.Branch("curvature_seed",&curvature_seed);
    t1s.Branch("invpT_seed",&invpT_seed);
    t1s.Branch("P_seed",&P_seed);

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

    t1s.Branch("status",&status);


    ///Read the MC data
    garana t;
    Int_t nentries = t.fChain->GetEntries();

    ////Cycle over the data

    for (size_t i=start_entry; i<start_entry+n_entries; i++) 
    {
        
        t.fChain->GetEntry(i);    
     
        if (t.TrackStartZ->size()==0) continue;
        
        
        /////Get the GArSoft Reco

        bool naive_direction;

        if(abs(t.TrackStartZ->at(0)-t.MCPStartZ->at(0))<100)
        {
            xyz_seed_old.SetXYZ(t.TrackStartX->at(0),t.TrackStartY->at(0),t.TrackStartZ->at(0));
            pxyz_seed_old.SetXYZ(t.TrackStartPX->at(0),t.TrackStartPY->at(0),t.TrackStartPZ->at(0));
            charge_seed_old=t.TrackStartQ->at(0);

            xyz_seed_old_bkw.SetXYZ(t.TrackEndX->at(0),t.TrackEndY->at(0),t.TrackEndZ->at(0));
            pxyz_seed_old_bkw.SetXYZ(t.TrackEndPX->at(0),t.TrackEndPY->at(0),t.TrackEndPZ->at(0));
            charge_seed_old_bkw=t.TrackEndQ->at(0);
            naive_direction = 1;
        }
        else
        {

            xyz_seed_old.SetXYZ(t.TrackEndX->at(0),t.TrackEndY->at(0),t.TrackEndZ->at(0));
            pxyz_seed_old.SetXYZ(t.TrackEndPX->at(0),t.TrackEndPY->at(0),t.TrackEndPZ->at(0));
            charge_seed_old=-1*t.TrackEndQ->at(0);

            xyz_seed_old_bkw.SetXYZ(t.TrackStartX->at(0),t.TrackStartY->at(0),t.TrackStartZ->at(0));
            pxyz_seed_old_bkw.SetXYZ(t.TrackStartPX->at(0),t.TrackStartPY->at(0),t.TrackStartPZ->at(0));
            charge_seed_old_bkw=-1*t.TrackStartQ->at(0);
            naive_direction = 0;
        }

        /////Get the hit clusters and order them
        bool Crossed[] = {0,0,0,0,0,0};
        size_t nCrossedPlanes = 0;
        for (Int_t j=0; j<t.TPCClusterX->size(); j++) 
        {           
            
            if (t.TrackIDNumber->at(0)==t.TPCClusterTrkIDNumber->at(j)) 
            {
                //if(i==949) std::cout<<t.TrackIDNumber->at(0)<<" "<<t.TPCClusterTrkIDNumber->at(j)<<std::endl;
                XYZVector Cluster_temp;
                Cluster_temp.SetXYZ(t.TPCClusterX->at(j),t.TPCClusterY->at(j),t.TPCClusterZ->at(j));
                xyz_plane_as_is.push_back(Cluster_temp); 

                ///Get number of crossed planes    
                for (Int_t k=0; k<6; k++)  
                {
                    if (t.TPCClusterZ->at(j)>Planes_Z[2*k] && t.TPCClusterZ->at(j)<Planes_Z[2*k+1] && Crossed[k]==0)
                    {
                        Crossed[k]+=1;
                        nCrossedPlanes+=1;
                    }
                }        
            }
        } 

         

        float lengthforwards = 0;
        std::vector<int> hlf ;
        float lengthbackwards = 0;
        std::vector<int> hlb ;
        
        sort_TPCClusters_along_track(xyz_plane_as_is,hlf,hlb,printlevelHelix,lengthforwards,lengthbackwards,1.0,2.0);

        
        if(naive_direction==0)
        {
          for(int ku=0;ku<hlb.size();ku++) xyz_plane.push_back(xyz_plane_as_is[hlb[ku]]);
        } 
        else{
            for(int ku=0;ku<hlb.size();ku++) xyz_plane.push_back(xyz_plane_as_is[hlf[ku]]);           
        }

       

        //std::vector<XYZVector> xyz_plane_bkw;
        //xyz_plane_bkw=xyz_plane;
        //std::reverse(xyz_plane_bkw.begin(),xyz_plane_bkw.end());
            
        
        ///Get the MC quantities
        xyz_MC.SetXYZ(t.MCPStartX->at(0),t.MCPStartY->at(0),t.MCPStartZ->at(0));
        pxyz_MC.SetXYZ(t.MCPStartPX->at(0),t.MCPStartPY->at(0),t.MCPStartPZ->at(0));
        charge_MC = t.PDG->at(0)>0 ? 1:-1;


        XYZVector dir_MC = pxyz_MC/(sqrt(pxyz_MC.Mag2()));
        float tanfactor_MC=1/TMath::Sqrt(dir_MC.Y()*dir_MC.Y()+dir_MC.Z()*dir_MC.Z()); //sqrt(1+tanlambda**2)
        float sign_MC=TMath::Sin(TMath::ATan2(dir_MC.Y(),dir_MC.Z()))/(dir_MC.Y()*tanfactor_MC);
        tanlambda_MC=dir_MC.X()*sign_MC*tanfactor_MC;
        sinphi_MC=dir_MC.Y()*sign_MC*tanfactor_MC;
        invpT_MC=charge_MC/(sqrt(pxyz_MC.Y()*pxyz_MC.Y()+pxyz_MC.Z()*pxyz_MC.Z()));
        curvature_MC =((0.3e-2)*B)*invpT_MC;
        

        Double_t dz_MC = xyz_plane.at(0).Z()-xyz_MC.Z();
        Double_t z2r_MC = curvature_MC*dz_MC;
        Double_t f1_MC=sinphi_MC;
        Double_t f2_MC=f1_MC + z2r_MC;
        Double_t r1_MC=TMath::Sqrt((1.-f1_MC)*(1.+f1_MC)), r2_MC=TMath::Sqrt((1.-f2_MC)*(1.+f2_MC));

        Double_t dy2dz_MC = (f1_MC+f2_MC)/(r1_MC+r2_MC);
        Double_t rot = TMath::ASin(r1_MC*f2_MC - r2_MC*f1_MC); 
                        
        xyz_MC.SetX(xyz_MC.X()+tanlambda_MC/curvature_MC*rot);
        xyz_MC.SetY(xyz_MC.Y()+dz_MC*dy2dz_MC);
        xyz_MC.SetZ(xyz_MC.Z()+dz_MC);       
        sinphi_MC+=z2r_MC;    

        

        XYZVector dir_seed_old = pxyz_seed_old/(sqrt(pxyz_seed_old.Mag2()));
        float tanfactor=1/TMath::Sqrt(dir_seed_old.Y()*dir_seed_old.Y()+dir_seed_old.Z()*dir_seed_old.Z()); //sqrt(1+tanlambda**2)
        float sign=TMath::Sin(TMath::ATan2(dir_seed_old.Y(),dir_seed_old.Z()))/(dir_seed_old.Y()*tanfactor);
        tanlambda_seed_old=dir_seed_old.X()*sign*tanfactor;
        sinphi_seed_old=dir_seed_old.Y()*sign*tanfactor;
        invpT_seed_old=charge_seed_old/(TMath::Sqrt(pxyz_seed_old.Y()*pxyz_seed_old.Y()+pxyz_seed_old.Z()*pxyz_seed_old.Z()));
        curvature_seed_old =((0.3e-2)*B)*invpT_seed_old;

        
        double forward=1.;
        double backwards=-1.;
        double lambda_seed,phi_seed;
        std::cout<<"Entry: "<<i;
        if(naive_direction==0)std::cout<<" backwards"<<std::endl;
        if(naive_direction==1)std::cout<<" forward"<<std::endl;
        
        
        ////Apply Alice Helix Fit
          
        makeSeed(xyz_plane,xyz_seed,curvature_seed,tanlambda_seed,sinphi_seed,forward,0,P_seed,xy_smear,Helix_Corr,nCrossedPlanes);
        
        invpT_seed=curvature_seed/((0.3e-2)*(B));
        if(naive_direction==1) 
        {
          invpT_seed*=-1;
          curvature_seed*=-1;
        }
        
        if(((xyz_plane[xyz_plane.size()-1].X()<xyz_plane[0].X()&&tanlambda_seed>0) || (xyz_plane[xyz_plane.size()-1].X()>xyz_plane[0].X()&&tanlambda_seed<0)))           
        {
          tanlambda_seed*=-1;
          sinphi_seed*=-1;
          invpT_seed*=-1;
          curvature_seed*=-1;
        }
        
        pxyz_seed.SetXYZ(tanlambda_seed/abs(invpT_seed),sinphi_seed/abs(invpT_seed),cos(asin(sinphi_seed))/abs(invpT_seed));
          
          
        
        if (Seedtype=="real")
        {
        ////Apply Helix Fit 
        //std::cout<<"Entry: "<<i;
        //if(naive_direction==0)std::cout<<" backwards"<<std::endl;
        //if(naive_direction==1)std::cout<<" forward"<<std::endl;
        status=Helix_Fit(xyz_plane,xyz_seed,curvature_seed,tanlambda_seed,sinphi_seed,forward,0);
        //std::cout<<std::endl;
            lambda_seed=tanlambda_seed;
            phi_seed=sinphi_seed;
            //
            if(OldCov) CalculateCovariance(xyz_plane,P_seed,curvature_seed,lambda_seed,phi_seed,1,0,0.3);
            tanlambda_seed=TMath::Tan(tanlambda_seed);
            sinphi_seed=TMath::Sin(sinphi_seed);
            if (Helix_Corr == "Eloss_MS" || Helix_Corr == "Eloss") {
              double crossLength = CalculatePath(xyz_plane,curvature_seed,tanlambda_seed,sinphi_seed);
              crossLength=nCrossedPlanes*Plane_thick;
              SeedMaterialCorrection(crossLength*rho,muon_mass,0.05,sqrt(pxyz_seed.Mag2()),(sqrt(pxyz_seed.Mag2())/muon_mass),rho,X0,X1,Ipar,ZA,curvature_seed,tanlambda_seed,sinphi_seed,P_seed,+1, Helix_Corr,nCrossedPlanes*Plane_thick/xx0);
            }
            invpT_seed=curvature_seed/((0.3e-2)*B);
            if(naive_direction==1) 
            {
              invpT_seed*=-1;
              curvature_seed*=-1;
            }
            
            if(((xyz_plane[xyz_plane.size()-1].X()<xyz_plane[0].X()&&tanlambda_seed>0) || (xyz_plane[xyz_plane.size()-1].X()>xyz_plane[0].X()&&tanlambda_seed<0)))           
            {
              tanlambda_seed*=-1;
              sinphi_seed*=-1;
              invpT_seed*=-1;
              curvature_seed*=-1;
            }

        
        
          pxyz_seed.SetXYZ(tanlambda_seed/abs(invpT_seed),sinphi_seed/abs(invpT_seed),cos(asin(sinphi_seed))/abs(invpT_seed));
          //std::cout<<"bubbolo"<<std::endl;
          
        }        
  

        
        
       //printlevelKalman=1;
       ////Apply Kalman Filter
      
      TMatrix P_seed_static=P_seed;
      std::vector<XYZVector> xyz_plane_bkw;
      xyz_plane_bkw=xyz_plane;
      std::reverse(xyz_plane_bkw.begin(),xyz_plane_bkw.end());
    
      if(Reuse_mat==0)
      {
        ///Helix Fit Seed as done in Alice
        makeSeed(xyz_plane_bkw,xyz_seed_bkw,curvature_seed_bkw,tanlambda_seed_bkw,sinphi_seed_bkw,backwards,printlevelHelix,P_seed_bkw,xy_smear,Helix_Corr,nCrossedPlanes);                
      }

       for(int nc=0;nc<ncycles;nc++)
       {
        
          if(nc>0)
          {
            xyz_seed.SetX(parvect_bkw.at(parvect_bkw.size()-1)[1]);
            xyz_seed.SetY(parvect_bkw.at(parvect_bkw.size()-1)[0]);
            xyz_seed.SetZ(xyz_plane_bkw.at(xyz_plane_bkw.size()-1).Z());
            curvature_seed=parvect_bkw.at(parvect_bkw.size()-1)[4]*B*0.3e-2;
            sinphi_seed=parvect_bkw.at(parvect_bkw.size()-1)[2];
            tanlambda_seed=parvect_bkw.at(parvect_bkw.size()-1)[3];
           
            if(Reuse_mat_cycle>0)
              {
                P_seed.Zero();
                P_seed[0][0]=Reuse_mat*Pt[Pt.size()-1][0][0];
                P_seed[1][1]=Reuse_mat*Pt[Pt.size()-1][1][1];
                P_seed[2][2]=Reuse_mat*Pt[Pt.size()-1][2][2];
                P_seed[3][3]=Reuse_mat*Pt[Pt.size()-1][3][3];
                P_seed[4][4]=Reuse_mat*Pt[Pt.size()-1][4][4];
              }
            if(!xht.empty())xht.clear();
            if(!yht.empty())yht.clear();
            if(!zht.empty())zht.clear();
            if(!zpost.empty())zpost.clear();
            if(!parvect.empty())parvect.clear();
            if(!predstept.empty())predstept.clear();
            if(!Pt.empty())Pt.clear();//(5,5);
            if(!PPredt.empty())PPredt.clear();//(5,5);
            if(!Rt.empty())Rt.clear();//(2,2);         
            if(!xht_bkw.empty())xht_bkw.clear();
            if(!yht_bkw.empty())yht_bkw.clear();
            if(!zht_bkw.empty())zht_bkw.clear();
            if(!zpost_bkw.empty())zpost_bkw.clear();
            if(!parvect_bkw.empty())parvect_bkw.clear();
            if(!predstept_bkw.empty())predstept_bkw.clear();
            if(!Pt_bkw.empty())Pt.clear();//(5,5);
            if(!PPredt_bkw.empty())PPredt_bkw.clear();//(5,5);
            if(!Rt_bkw.empty())Rt_bkw.clear();//(2,2);
          }

          KalmanFit(xht,yht,zht,parvect,predstept,Pt,PPredt,Rt,zpost,xyz_plane,Ry,Rx,Ryx,xyz_seed,tanlambda_seed,curvature_seed,sinphi_seed,forward,status,printlevelKalman,dEreco,dxreco,Energy_loss_corr,CorrTime,Fixed_Cov,Smear,xy_smear,P_seed,Seedtype,MS);
          
          
          if(!parvect.empty() && status==1)
            {    
              //std::cout<<status<<std::endl;      
              double backwards=-1.;
              
              if(Reuse_mat>0)
                {
                  P_seed_bkw.Zero();
                  P_seed_bkw[0][0]=Reuse_mat*Pt[Pt.size()-1][0][0];
                  P_seed_bkw[1][1]=Reuse_mat*Pt[Pt.size()-1][1][1];
                  P_seed_bkw[2][2]=Reuse_mat*Pt[Pt.size()-1][2][2];
                  P_seed_bkw[3][3]=Reuse_mat*Pt[Pt.size()-1][3][3];
                  P_seed_bkw[4][4]=Reuse_mat*Pt[Pt.size()-1][4][4];
                }

              xyz_seed_bkw.SetX(parvect.at(parvect.size()-1)[1]);
              xyz_seed_bkw.SetY(parvect.at(parvect.size()-1)[0]);
              xyz_seed_bkw.SetZ(xyz_plane.at(xyz_plane.size()-1).Z());
              curvature_seed_bkw=parvect.at(parvect.size()-1)[4]*B*0.3e-2;
              sinphi_seed_bkw=parvect.at(parvect.size()-1)[2];
              tanlambda_seed_bkw=parvect.at(parvect.size()-1)[3];
            
              Pt_bkw=Pt;
              //printlevelKalman=1;
                    
              KalmanFit(xht_bkw,yht_bkw,zht_bkw,parvect_bkw,predstept_bkw,Pt_bkw,PPredt_bkw,Rt_bkw,zpost_bkw,xyz_plane_bkw,Ry,Rx,Ryx,xyz_seed_bkw,tanlambda_seed_bkw,curvature_seed_bkw,sinphi_seed_bkw,backwards,status,printlevelKalman,dEreco,dxreco,Energy_loss_corr,CorrTime,Fixed_Cov,Smear,xy_smear,P_seed_bkw,Seedtype,MS);
              
              //printlevelKalman=0;
              //for(int nu=0;nu<parvect_bkw.size();nu++)
              //{
              //std::cout<<"parvect_bkw: "<<parvect_bkw[nu][0]<<" "<<parvect_bkw[nu][1]<<" "<<parvect_bkw[nu][2]<<" "<<parvect_bkw[nu][3]<<" "<<parvect_bkw[nu][4]<<std::endl;
              //}
              
            }
          else break;
       }

       if(status==1)
       {
        if((parvect_bkw[parvect_bkw.size()-1][4]-invpT_MC)/invpT_MC>1)
        {
        std::cout<<"Seed     "<<xyz_seed.X()<<" "<<xyz_seed.Y()<<" "<<xyz_seed.Z()<<" "<<sinphi_seed<<" "<<tanlambda_seed<<" "<<invpT_seed<<std::endl;
        std::cout<<"MC       "<<xyz_MC.X()<<" "<<xyz_MC.Y()<<" "<<xyz_MC.Z()<<" "<<sinphi_MC<<" "<<tanlambda_MC<<" "<<invpT_MC<<std::endl;
        std::cout<<"Seed Old "<<xyz_seed_old.X()<<" "<<xyz_seed_old.Y()<<" "<<xyz_seed_old.Z()<<" "<<sinphi_seed_old<<" "<<tanlambda_seed_old<<" "<<invpT_seed_old<<std::endl;
        if(status==1) std::cout<<"Kalman   "<<xyz_plane_bkw[xyz_plane_bkw.size()-1].X()<<" "<<xyz_plane_bkw[xyz_plane_bkw.size()-1].Y()<<" "<<xyz_plane_bkw[xyz_plane_bkw.size()-1].Z()<<" "<<parvect_bkw[parvect_bkw.size()-1][2]<<" "<<parvect_bkw[parvect_bkw.size()-1][3]<<" "<<parvect_bkw[parvect_bkw.size()-1][4]<<std::endl;
        }
       }

    
        
       
       
        ////Fill The tree
        if(status!=1)std::cout<<"status: "<<status<<std::endl;
        t1s.Fill();
        
       

        ////Empty the containers
        if(!xyz_plane.empty())xyz_plane.clear();
        //if(!xyz_plane_bkw.empty())xyz_plane_bkw.clear();
        if(!xyz_plane_as_is.empty())xyz_plane_as_is.clear();
        if(!xht.empty())xht.clear();
        if(!yht.empty())yht.clear();
        if(!zht.empty())zht.clear();
        if(!zpost.empty())zpost.clear();
        if(!parvect.empty())parvect.clear();
        if(!predstept.empty())predstept.clear();
        if(!Pt.empty())Pt.clear();//(5,5);

        if(!PPredt.empty())PPredt.clear();//(5,5);
        if(!Rt.empty())Rt.clear();//(2,2);

        
        if(!xht_bkw.empty())xht_bkw.clear();
        if(!yht_bkw.empty())yht_bkw.clear();
        if(!zht_bkw.empty())zht_bkw.clear();
        if(!zpost_bkw.empty())zpost_bkw.clear();
        if(!parvect_bkw.empty())parvect_bkw.clear();
        if(!predstept_bkw.empty())predstept_bkw.clear();
        if(!Pt_bkw.empty())Pt.clear();//(5,5);
        if(!PPredt_bkw.empty())PPredt_bkw.clear();//(5,5);
        if(!Rt_bkw.empty())Rt_bkw.clear();//(2,2);
        

        status=true;

        
        
  }

    ///////Write the tree to File
  t1s.Write();
    
}