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
    chain.Add("/pnfs/dune/persistent/users/ebrianne/ProductionSamples/ND-LAr/nd_hall_dayone_lar_SPY_v2_wMuID/Anatree/neutrino/*");

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

    TGraph2D* g_p_theta_Enu_prim_mu=new TGraph2D();
    TGraph2D* g_p_theta_Enu_nonprim_mu=new TGraph2D();
    TGraph2D* g_p_theta_Enu_prim_antimu=new TGraph2D();
    TGraph2D* g_p_theta_Enu_nonprim_antimu=new TGraph2D();

    //TH3F* h_p_theta_Enu_prim_mu = new TH3F("h_p_theta_Enu_prim_mu", "(p,#theta,E(#nu)) prim mu", 100, 0, 20, 100, 0, 5,100,0,20);
    //TH3F* h_p_theta_Enu_nonprim_mu = new TH3F("h_p_theta_Enu_nonprim_mu", "(p,#theta,E(#nu)) non prim mu", 100, 0, 20, 100, 0, 5,100,0,20);
    //TH3F* h_p_theta_Enu_prim_antimu = new TH3F("h_p_theta_Enu_prim_antimu", "(p,#theta,E(#nu)) prim antimu", 100, 0, 20, 100, 0, 5,100,0,20);
    //TH3F* h_p_theta_Enu_nonprim_antimu = new TH3F("h_p_theta_Enu_nonprim_antimu", "(p,#theta,E(#nu)) nonprim antimu", 100, 0, 20, 100, 0, 5,100,0,20);
    
    
    int g1=0;
    int g2=0;
    int g3=0;
    int g4=0;
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

               if(PDG->at(j)==13 && PDGMother->at(j)==0 && p<=20 && Theta<=5){
                   //h_p_theta_Enu_prim_mu->Fill(p,Theta,Enu);
                   g_p_theta_Enu_prim_mu->SetPoint(g1,p,Theta,Enu);
                   g1++;
                   }
               if(PDG->at(j)==13 && PDGMother->at(j)!=0 && p<=20 && Theta<=5){
                   //h_p_theta_Enu_nonprim_mu->Fill(p,Theta,Enu);
                   g_p_theta_Enu_nonprim_mu->SetPoint(g2,p,Theta,Enu);
                   g2++;
                   }
               if(PDG->at(j)==-13 && PDGMother->at(j)==0 && p<=20 && Theta<=5){
                   //h_p_theta_Enu_prim_antimu->Fill(p,Theta,Enu);       
                   g_p_theta_Enu_prim_antimu->SetPoint(g3,p,Theta,Enu);
                   g3++;            
                   }
               if(PDG->at(j)==-13 && PDGMother->at(j)!=0 && p<=20 && Theta<=5){
                   //h_p_theta_Enu_nonprim_antimu->Fill(p,Theta,Enu);
                   g_p_theta_Enu_nonprim_antimu->SetPoint(g4,p,Theta,Enu);
                   g4++;
                   }
               

               
           }
        }
        
        

    }
 
    
    //TH2D *hProjection = (TH2D*) h_p_theta_Enu_prim_mu->Project3D("xy");
    
    TCanvas *mccanvasgraph = new TCanvas("mccanvasgraph","",1000,800);

    const Int_t NRGBs = 5;
    const Int_t NCont = 30;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 }; 
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);   
    gStyle->SetNumberContours(NCont);
    
    TH1 *frame = mccanvasgraph->DrawFrame(0,0,20,5);
    frame->GetXaxis()->SetTitle("p [GeV/c]");
    frame->GetYaxis()->SetTitle("#theta [degree]"); 
    frame->SetTitle("Primary #mu^{-} E_{#nu} [GeV] contour");

    g_p_theta_Enu_prim_mu->Draw("CONT z same");
    mccanvasgraph->Print("16g_p_theta_Enu_prim_mu.png");




    TCanvas *mccanvasgraph2 = new TCanvas("mccanvasgraph2","",1000,800);

    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);   
    gStyle->SetNumberContours(NCont);
    
    TH1 *frame2 = mccanvasgraph2->DrawFrame(0,0,20,5);
    frame2->GetXaxis()->SetTitle("p [GeV/c]");
    frame2->GetYaxis()->SetTitle("#theta [degree]"); 
    frame2->SetTitle("Non primary #mu^{-} E_{#nu} [GeV] contour");

    g_p_theta_Enu_nonprim_mu->Draw("CONT z same");
    mccanvasgraph2->Print("17g_p_theta_Enu_nonprim_mu.png");




    TCanvas *mccanvasgraph3 = new TCanvas("mccanvasgraph3","",1000,800);

    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);   
    gStyle->SetNumberContours(NCont);
    
    TH1 *frame3 = mccanvasgraph3->DrawFrame(0,0,20,5);
    frame3->GetXaxis()->SetTitle("p [GeV/c]");
    frame3->GetYaxis()->SetTitle("#theta [degree]"); 
    frame3->SetTitle("Primary #mu^{+} E_{#nu} [GeV] contour");

    g_p_theta_Enu_prim_antimu->Draw("CONT z same");
    mccanvasgraph3->Print("18g_p_theta_Enu_prim_antimu.png");




    TCanvas *mccanvasgraph4 = new TCanvas("mccanvasgraph4","",1000,800);

    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);   
    gStyle->SetNumberContours(NCont);
    
    TH1 *frame4 = mccanvasgraph4->DrawFrame(0,0,20,5);
    frame4->GetXaxis()->SetTitle("p [GeV/c]");
    frame4->GetYaxis()->SetTitle("#theta [degree]"); 
    frame4->SetTitle("Non primary #mu^{+} E_{#nu} [GeV] contour");

    g_p_theta_Enu_nonprim_antimu->Draw("CONT z same");
    mccanvasgraph4->Print("19g_p_theta_Enu_nonprim_antimu.png");
    /*
    TCanvas *mccanvasproj = new TCanvas("mccanvasproj","",1000,800);
    hProjection->SetTitle("(p,#theta,E_#nu) prim mu;p [GeV/c];#theta [degree]");
    hProjection->Draw("COLZ");
    mccanvasproj->Print("projection.png");

    
    
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
    */

}