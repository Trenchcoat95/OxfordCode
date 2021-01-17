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

void l2g_nu()
{
    TChain chain("/anatree/GArAnaTree");
    chain.Add("/pnfs/dune/persistent/users/ebrianne/ProductionSamples/ND-LAr/nd_hall_dayone_lar_SPY_v2_wMuID/Anatree/neutrino/neutrino.nd_hall_dayone_lar_SPY_v2_wMuID.volArgonCubeActive.Ev973000.Ev973999.2037.anatree.root");

    //neutrino.nd_hall_dayone_lar_SPY_v2_wMuID.volArgonCubeActive.Ev973000.Ev973999.2037.anatree.root

    vector<int>     *PDG=0;
    vector<int>     *PDGMother=0;
    vector<float>   *MCNuPx=0;
    vector<float>   *MCNuPy=0;
    vector<float>   *MCNuPz=0;
    vector<float>   *MCPStartPX=0;
    vector<float>   *MCPStartPY=0;
    vector<float>   *MCPStartPZ=0;


   
    TBranch *b_PDG=0;
    TBranch *b_PDGMother=0;
    TBranch *b_MCNuPx=0;
    TBranch *b_MCNuPy=0;
    TBranch *b_MCNuPz=0;
    TBranch *b_MCPStartPX=0;
    TBranch *b_MCPStartPY=0;
    TBranch *b_MCPStartPZ=0;




    
    chain.SetBranchAddress("PDG", &PDG, &b_PDG);
    chain.SetBranchAddress("PDGMother", &PDGMother, &b_PDGMother);
    chain.SetBranchAddress("MCNuPx", &MCNuPx, &b_MCNuPx);
    chain.SetBranchAddress("MCNuPy", &MCNuPy, &b_MCNuPy);
    chain.SetBranchAddress("MCNuPz", &MCNuPz, &b_MCNuPz);
    chain.SetBranchAddress("MCPStartPX", &MCPStartPX, &b_MCPStartPX);
    chain.SetBranchAddress("MCPStartPY", &MCPStartPY, &b_MCPStartPY);
    chain.SetBranchAddress("MCPStartPZ", &MCPStartPZ, &b_MCPStartPZ);



    chain.SetBranchStatus("*",0);
    chain.SetBranchStatus("PDG",1);
    chain.SetBranchStatus("PDGMother",1);
    chain.SetBranchStatus("MCNuPx",1);
    chain.SetBranchStatus("MCNuPy",1);
    chain.SetBranchStatus("MCNuPz",1);
    chain.SetBranchStatus("MCPStartPX",1);
    chain.SetBranchStatus("MCPStartPY",1);
    chain.SetBranchStatus("MCPStartPZ",1);


    

    Int_t nentries = (Int_t)chain.GetEntries();

    TH3F* h_p_theta_Enu_prim_mu = new TH3F("h_p_theta_Enu_prim_mu", "(p,#theta,E(#nu)) prim mu", 100, 0, 20, 100, 0, 5,100,0,20);
    TH3F* h_p_theta_Enu_nonprim_mu = new TH3F("h_p_theta_Enu_nonprim_mu", "(p,#theta,E(#nu)) non prim mu", 100, 0, 20, 100, 0, 5,100,0,20);
    TH3F* h_p_theta_Enu_prim_antimu = new TH3F("h_p_theta_Enu_prim_antimu", "(p,#theta,E(#nu)) prim antimu", 100, 0, 20, 100, 0, 5,100,0,20);
    TH3F* h_p_theta_Enu_nonprim_antimu = new TH3F("h_p_theta_Enu_nonprim_antimu", "(p,#theta,E(#nu)) nonprim antimu", 100, 0, 20, 100, 0, 5,100,0,20);
    
    

    for (Int_t i=0; i<nentries; i++) 
    {
        chain.GetEntry(i);
        float pxNu= MCNuPx->at(0);
        float pyNu= MCNuPy->at(0);
        float pzNu= MCNuPz->at(0);
        float Enu = sqrt(pxNu*pxNu+pyNu*pyNu+pzNu*pzNu);

        for (Int_t j=0; j<PDG->size(); j++) 
        {
           if(PDG->at(j)==13 || PDG->at(j)==-13)
           {
               float px = MCPStartPX->at(j);
               float py = MCPStartPY->at(j);
               float pz = MCPStartPZ->at(j);
               float p = sqrt(px*px+py*py+pz*pz);              
               float Theta = acos(pz/p);

               if(PDG->at(j)==13 && PDGMother->at(j)==0){h_p_theta_Enu_prim_mu->Fill(p,Theta,Enu);}
               if(PDG->at(j)==13 && PDGMother->at(j)!=0){h_p_theta_Enu_nonprim_mu->Fill(p,Theta,Enu);}
               if(PDG->at(j)==-13 && PDGMother->at(j)==0){h_p_theta_Enu_prim_antimu->Fill(p,Theta,Enu);}
               if(PDG->at(j)==-13 && PDGMother->at(j)!=0){h_p_theta_Enu_nonprim_antimu->Fill(p,Theta,Enu);}
               

               
           }
        }
        
        

    }
 

    
    TCanvas *mccanvasp_theta_Enu_prim_mu = new TCanvas("mccanvasp_theta_Enu_prim_mu","",1000,800);
    h_p_theta_Enu_prim_mu->SetTitle("(p,#theta,E_#nu) prim mu;p [GeV/c];#theta [degree];E_#nu [Gev]");
    h_p_theta_Enu_prim_mu->Draw("ISO");
    mccanvasp_theta_Enu_prim_mu->Print("16_p_theta_Enu_prim_mu.png");

    TCanvas *mccanvasp_theta_Enu_nonprim_mu = new TCanvas("mccanvasp_theta_Enu_nonprim_mu","",1000,800);
    h_p_theta_Enu_nonprim_mu->SetTitle("(p,#theta,E_#nu) non prim mu;p [GeV/c];#theta [degree];E_#nu [Gev]");
    h_p_theta_Enu_nonprim_mu->Draw("ISO");
    mccanvasp_theta_Enu_nonprim_mu->Print("16_p_theta_Enu_nonprim_mu.png");

    TCanvas *mccanvasp_theta_Enu_prim_antimu = new TCanvas("mccanvasp_theta_Enu_prim_antimu","",1000,800);
    h_p_theta_Enu_prim_antimu->SetTitle("(p,#theta,E_#nu) prim antimu;p [GeV/c];#theta [degree];E_#nu [Gev]");
    h_p_theta_Enu_prim_antimu->Draw("ISO");
    mccanvasp_theta_Enu_prim_antimu->Print("17_p_theta_Enu_prim_antimu.png");

    TCanvas *mccanvasp_theta_Enu_nonprim_antimu = new TCanvas("mccanvasp_theta_Enu_nonprim_antimu","",1000,800);
    h_p_theta_Enu_nonprim_antimu->SetTitle("(p,#theta,E_#nu) non prim antimu;p [GeV/c];#theta [degree];E_#nu [Gev]");
    h_p_theta_Enu_nonprim_antimu->Draw("ISO");
    mccanvasp_theta_Enu_nonprim_antimu->Print("18_p_theta_Enu_nonprim_antimu.png");
    

}