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

void l2gread()
{
    TChain chain("/anatree/GArAnaTree");
    chain.Add("/mnt/c/Users/fedeb/Documents/Universita/Federico_2020-2021/OxfordCode/l2g/data/*.root");

    //neutrino.nd_hall_dayone_lar_SPY_v2_wMuID.volArgonCubeActive.Ev973000.Ev973999.2037.anatree.root

    vector<int>     *PDG=0;
    vector<int>     *PDGMother=0;
    vector<float>   *MCPStartPX=0;
    vector<float>   *MCPStartPY=0;
    vector<float>   *MCPStartPZ=0;
    vector<float>   *MCPStartX=0;
    vector<float>   *MCPStartY=0;
    vector<float>   *MCPStartZ=0;


   
    TBranch *b_PDG=0;
    TBranch *b_PDGMother=0;
    TBranch *b_MCPStartPX=0;
    TBranch *b_MCPStartPY=0;
    TBranch *b_MCPStartPZ=0;
    TBranch *b_MCPStartX=0;
    TBranch *b_MCPStartY=0;
    TBranch *b_MCPStartZ=0;



    
    chain.SetBranchAddress("PDG", &PDG, &b_PDG);
    chain.SetBranchAddress("PDGMother", &PDGMother, &b_PDGMother);
    chain.SetBranchAddress("MCPStartPX", &MCPStartPX, &b_MCPStartPX);
    chain.SetBranchAddress("MCPStartPY", &MCPStartPY, &b_MCPStartPY);
    chain.SetBranchAddress("MCPStartPZ", &MCPStartPZ, &b_MCPStartPZ);
    chain.SetBranchAddress("MCPStartX", &MCPStartX, &b_MCPStartX);
    chain.SetBranchAddress("MCPStartY", &MCPStartY, &b_MCPStartY);
    chain.SetBranchAddress("MCPStartZ", &MCPStartZ, &b_MCPStartZ);




    chain.SetBranchStatus("*",0);
    chain.SetBranchStatus("PDG",1);
    chain.SetBranchStatus("PDGMother",1);
    chain.SetBranchStatus("MCPStartPX",1);
    chain.SetBranchStatus("MCPStartPY",1);
    chain.SetBranchStatus("MCPStartPZ",1);
    chain.SetBranchStatus("MCPStartX",1);
    chain.SetBranchStatus("MCPStartY",1);
    chain.SetBranchStatus("MCPStartZ",1);

    

    Int_t nentries = (Int_t)chain.GetEntries();


    
    TH1F* hp = new TH1F("hp", "Muon Momentum spectrum", 100, 0.0, 20.0);
    TH2F* hxy = new TH2F("hxy", "Muon production point distribution", 100, -500, 500, 100, -300, 300);
    TH2F* hzy = new TH2F("hzy", "Muon production point distribution", 100, 0, 2000, 100, -300, 300);

    TH2F* h_theta_p = new TH2F("h_theta_p", "Theta VS p", 100, 0, 20, 100, 0, 100);
    TH2F* h_p_theta = new TH2F("h_p_theta", "p VS theta", 100, 0, 100, 100, 0, 20);
    TH2F* h_thetaProb_p = new TH2F("h_thetaProb_p", "P(#theta|p)", 100, 0, 20, 100, 0, 100);
    TH2F* h_theta_pProb = new TH2F("h_theta_pProb", "P(p|#theta)", 100, 0, 20, 100, 0, 100);
    TH2F* h_theta_pProbswap = new TH2F("h_theta_pProbswap", "P(p|#theta)", 100, 0, 100, 100, 0, 20);

    TH2F* h_px_p = new TH2F("h_px_p", "px VS p", 100, 0, 20, 100, 0, 5);
    TH2F* h_pT_p = new TH2F("h_pT_p", "pT VS p", 100, 0, 20, 100, 0, 20);
    TH2F* h_P_px_p = new TH2F("h_P_px_p", "P(px|p)", 100, 0, 20, 100, 0, 5);
    TH2F* h_P_pT_p = new TH2F("h_P_pT_p", "P(pT|p)", 100, 0, 20, 100, 0, 20);

    TH2F* h_frac_pxp_p = new TH2F("h_frac_pxp_p", "px/p VS p", 100, 0, 20, 100, 0, 1);
    TH2F* h_frac_pTp_p = new TH2F("h_frac_pTp_p", "pT/p VS p", 100, 0, 20, 100, 0, 1);
    TH2F* h_P_frac_pxp_p = new TH2F("h_P_frac_pxp_p", "P(px/p|p)", 100, 0, 20, 100, 0, 1);
    TH2F* h_P_frac_pTp_p = new TH2F("h_P_frac_pTp_p", "P(pT/p|p)", 100, 0, 20, 100, 0, 1);
    
    
    int doublemuon=0;
    int pp=0;
    int pnp=0;
    int npnp=0;
 
    for (Int_t i=0; i<nentries; i++) 
    {
        chain.GetEntry(i);
        int ID=0;
        int muons=0;
        int prim=0;
    

        for (Int_t j=0; j<PDG->size(); j++) 
        {
           if(PDG->at(j)==13)// || PDG->at(j)==-13)
           {
               float px = MCPStartPX->at(j);
               float py = MCPStartPY->at(j);
               float pz = MCPStartPZ->at(j);
               float p = sqrt(px*px+py*py+pz*pz);
               float pT = sqrt(pz*pz+py*py);
               

               float x = MCPStartX->at(j);
               float y = MCPStartY->at(j);
               float z = MCPStartZ->at(j);
               
               float Theta = acos(pz/p)* (180.0/3.141592653589793238463);

               hp->Fill(p);
               if(PDGMother->at(j)==0)hxy->Fill(x,y);
               if(PDGMother->at(j)==0)hzy->Fill(z,y);
               h_theta_p->Fill(p,Theta);
               h_p_theta->Fill(Theta,p);
               h_px_p->Fill(p,px);
               h_pT_p->Fill(p,pT);
               h_frac_pxp_p->Fill(p,px/p);
               h_frac_pTp_p->Fill(p,pT/p);
  
               muons++;
               //std::cout<<"PDGMother "<<PDGMother->at(j)<<std::endl;
               if (PDGMother->at(j)==0)
               {
                   prim++;
                   //std::cout<<"prim "<<prim<<std::endl;
                }
               //if (muons>1){std::cout<<MCPProc->at(j)<<std::endl;}

               

               
           }
        }
        
        if (muons>1){doublemuon++;}
        if (muons>1 && prim==0){npnp++;}
        if (muons>1 && prim==1){pnp++;}
        if (muons>1 && prim==2){pp++;}
        
 

    }
    std::cout<<"doubelmuon="<<doublemuon<<" tot="<<nentries<<std::endl;
    std::cout<<"pp="<<pp<<"pnp="<<pnp<<"npnp="<<npnp<<std::endl;
 
    for (Int_t x=1; x<=100; x++)
    {
        float sum=0;
        float sum2=0;
        float sumpx=0;
        float sumpT=0;
        float sumpx_p=0;
        float sumpT_p=0;
        float sumcheck=0;
        /*
        for (Int_t y=1; y<=100; y++)
        {
            sum+=h_theta_p->GetBinContent(h_theta_p->GetBin(x,y));
            sumpx+=h_px_p->GetBinContent(h_px_p->GetBin(x,y));
            sumpT+=h_pT_p->GetBinContent(h_pT_p->GetBin(x,y));
            sumpx_p+=h_frac_pxp_p->GetBinContent(h_frac_pxp_p->GetBin(x,y));
            sumpT_p+=h_frac_pTp_p->GetBinContent(h_frac_pTp_p->GetBin(x,y));
        }
        */
        for (Int_t y=1; y<=100; y++)
        {
            if(sum<h_theta_p->GetBinContent(h_theta_p->GetBin(x,y))) {sum=h_theta_p->GetBinContent(h_theta_p->GetBin(x,y));}
            if(sum2<h_p_theta->GetBinContent(h_p_theta->GetBin(x,y))) {sum2=h_p_theta->GetBinContent(h_p_theta->GetBin(x,y));}
            if(sumpx<h_px_p->GetBinContent(h_px_p->GetBin(x,y)))  {sumpx=h_px_p->GetBinContent(h_px_p->GetBin(x,y));}
            if(sumpT<h_pT_p->GetBinContent(h_pT_p->GetBin(x,y)))  {sumpT=h_pT_p->GetBinContent(h_pT_p->GetBin(x,y));}
            if(sumpx_p<h_frac_pxp_p->GetBinContent(h_frac_pxp_p->GetBin(x,y)))  {sumpx_p=h_frac_pxp_p->GetBinContent(h_frac_pxp_p->GetBin(x,y));}
            if(sumpT_p<h_frac_pTp_p->GetBinContent(h_frac_pTp_p->GetBin(x,y)))  {sumpT_p=h_frac_pTp_p->GetBinContent(h_frac_pTp_p->GetBin(x,y));}
        }
        
        for (Int_t y=1; y<=100; y++)
        {
            float bin=0;
            float bin2=0;
            float binpx=0;
            float binpT=0;
            float binpx_p=0;
            float binpT_p=0;

            if (sum!=0){ bin=h_theta_p->GetBinContent(h_theta_p->GetBin(x,y))/sum;}
            if (sum2!=0){ bin2=h_p_theta->GetBinContent(h_p_theta->GetBin(x,y))/sum2;}
            if (sumpx!=0){binpx=h_px_p->GetBinContent(h_px_p->GetBin(x,y))/sumpx;}
            if (sumpT!=0){binpT=h_pT_p->GetBinContent(h_pT_p->GetBin(x,y))/sumpT;}
            if (sumpx_p!=0){binpx_p=h_frac_pxp_p->GetBinContent(h_frac_pxp_p->GetBin(x,y))/sumpx_p;}
            if (sumpT_p!=0){binpT_p=h_frac_pTp_p->GetBinContent(h_frac_pTp_p->GetBin(x,y))/sumpT_p;}
            
            h_theta_pProb->SetBinContent(x,y,bin);
            h_theta_pProbswap->SetBinContent(x,y,bin2);
            h_P_px_p->SetBinContent(x,y,binpx);
            h_P_pT_p->SetBinContent(x,y,binpT);
            h_P_frac_pxp_p->SetBinContent(x,y,binpx_p);
            h_P_frac_pTp_p->SetBinContent(x,y,binpT_p);
            sumcheck+=bin;
        }
        //std::cout<<"sumcheck: "<<sumcheck<<std::endl;
    }

    for (Int_t y=1; y<=100; y++)
    {
        float sum=0;
        float sumcheck=0;
        /*
        for (Int_t x=1; x<=100; x++)
        {
            sum+=h_theta_p->GetBinContent(h_theta_p->GetBin(x,y));
        }
        */
        for (Int_t x=1; x<=100; x++)
        {
            if (sum<h_theta_p->GetBinContent(h_theta_p->GetBin(x,y)))  {sum=h_theta_p->GetBinContent(h_theta_p->GetBin(x,y));}
        }
        
        for (Int_t x=1; x<=100; x++)
        {
            float bin=0;

            if (sum!=0){bin=h_theta_p->GetBinContent(h_theta_p->GetBin(x,y))/sum;}     

            h_thetaProb_p->SetBinContent(x,y,bin);
            sumcheck+=bin;
        }
        //std::cout<<"sum: "<<sum<<std::endl;
        //std::cout<<"sumcheck: "<<sumcheck<<std::endl;
    }
    
    TCanvas *mccanvasp = new TCanvas("mccanvasp","",1000,800);
    hp->SetTitle("Muon Momentum spectrum;p (GeV/#it{c});Muons");
    hp->Draw();
    mccanvasp->Print("0_Muon_Spectrum.png");

    TCanvas *mccanvasxy = new TCanvas("mccanvasxy","",1000,800);
    gStyle->SetOptStat(0);
    hxy->SetTitle("Muon production vertex distribution;x (cm);y (cm)");
    hxy->Draw("COLZ");
    mccanvasxy->Print("1_Muon_Startxy_prim.png");

    TCanvas *mccanvaszy = new TCanvas("mccanvaszy","",1000,800);
    gStyle->SetOptStat(0);
    hzy->SetTitle("Muon production vertex distribution;z (cm);y (cm)");
    hzy->Draw("COLZ");
    mccanvaszy->Print("2_Muon_Startzy_prim.png");
    



    TCanvas *mccanvasThetap = new TCanvas("mccanvasThetap","",1000,800);
    gStyle->SetOptStat(0);
    h_theta_p->SetTitle("#theta VS p;p (GeV/#it{c});#theta (degree)");
    h_theta_p->Draw("COLZ");
    mccanvasThetap->Print("3_ThetaVSp.png");
    
    TCanvas *mccanvasThetapProb = new TCanvas("mccanvasThetapProb","",1000,800);
    gStyle->SetOptStat(0);
    h_theta_pProb->SetTitle("P(#theta|p);p (GeV/#it{c});#theta (degree)");
    h_theta_pProb->Draw("COLZ");
    mccanvasThetapProb->Print("4_ThetaProbVSp.png");
    
    TCanvas *mccanvasThetaProbp = new TCanvas("mccanvasThetaProbp","",1000,800);
    gStyle->SetOptStat(0);
    h_theta_pProbswap->SetTitle("P(p|#theta);#theta (degree);p (GeV/#it{c})");
    h_theta_pProbswap->Draw("COLZ");
    mccanvasThetaProbp->Print("5_ThetaVSpProb.png");

    


    TCanvas *mccanvaspxp = new TCanvas("mccanvaspxp","",1000,800);
    gStyle->SetOptStat(0);
    h_px_p->SetTitle("px VS p;p (GeV/#it{c});px (GeV/#it{c})");
    h_px_p->Draw("COLZ");
    mccanvaspxp->Print("6_pxVSp.png");
    
    TCanvas *mccanvaspTp = new TCanvas("mccanvaspTp","",1000,800);
    gStyle->SetOptStat(0);
    h_pT_p->SetTitle("pT VS p;p (GeV/#it{c});pT (GeV/#it{c})");
    h_pT_p->Draw("COLZ");
    mccanvaspTp->Print("7_pTVSp.png");

    TCanvas *mccanvasPpxp = new TCanvas("mccanvasPpxp","",1000,800);
    gStyle->SetOptStat(0);
    h_P_px_p->SetTitle("P(px|p);p (GeV/#it{c});px (GeV/#it{c})");
    h_P_px_p->Draw("COLZ");
    mccanvasPpxp->Print("8_P_ppx.png");
    
    TCanvas *mccanvasPpTp = new TCanvas("mccanvasPpTp","",1000,800);
    gStyle->SetOptStat(0);
    h_P_pT_p->SetTitle("P(pT|p);p (GeV/#it{c});pT (GeV/#it{c})");
    h_P_pT_p->Draw("COLZ");
    mccanvasPpTp->Print("9_P_ppT.png");
    


    TCanvas *mccanvasfracpxp = new TCanvas("mccanvasfracpxp","",1000,800);
    gStyle->SetOptStat(0);
    h_frac_pxp_p->SetTitle("px/p VS p;p;px/p (GeV/#it{c})");
    h_frac_pxp_p->Draw("COLZ");
    mccanvasfracpxp->Print("10_frac_pxpVSp.png");
    
    TCanvas *mccanvasfracpTp = new TCanvas("mccanvasfracpTp","",1000,800);
    gStyle->SetOptStat(0);
    h_frac_pTp_p->SetTitle("pT/p VS p;pT/p;p (GeV/c)");
    h_frac_pTp_p->Draw("COLZ");
    mccanvasfracpTp->Print("11_frac_pTpVSp.png");

    TCanvas *mccanvasfracPpxp = new TCanvas("mccanvasfracPpxp","",1000,800);
    gStyle->SetOptStat(0);
    h_P_frac_pxp_p->SetTitle("P(px/p|p);p (GeV/#it{c});px/p");
    h_P_frac_pxp_p->Draw("COLZ");
    mccanvasfracPpxp->Print("12_P_fracppx.png");
    
    TCanvas *mccanvasfracPpTp = new TCanvas("mccanvasfracPpTp","",1000,800);
    gStyle->SetOptStat(0);
    h_P_frac_pTp_p->SetTitle("P(pT/p|p);p (GeV/#it{c});pT/p");
    h_P_frac_pTp_p->Draw("COLZ");
    mccanvasfracPpTp->Print("13_P_fracppT.png");
    
}