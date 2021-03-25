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



void l2g_multimu_spill()
{
    std::cout<<"Reading data"<<std::endl;
    TChain chain("/anatree/GArAnaTree");
    chain.Add("/mnt/c/Users/fedeb/Documents/Universita/Federico_2020-2021/OxfordCode/l2g/data_single_spill/*.root");

    TVector3 GArCenter(0,-68.287,1486);
    float GAr_r = 349.9;
    float GAr_L = 669.6;
    
     
    TVector3 LArCenter(0,0,666);
    TVector3 LArDim(714.7,301.02,509.1);
    TVector3 LArUpLimit = LArCenter+0.5*LArDim;
    TVector3 LArLowLimit = LArCenter-0.5*LArDim;

    ///pnfs/dune/persistent/users/ebrianne/ProductionSamples/ND-LAr/nd_hall_dayone_lar_SPY_v2_wMuID/Anatree/neutrino/neutrino.nd_hall_dayone_lar_SPY_v2_wMuID.volArgonCubeActive.Ev973000.Ev973999.2037.anatree.root
    //neutrino.nd_hall_dayone_lar_SPY_v2_wMuID.volArgonCubeActive.Spill96000.Spill96999.28100.overlay.anatree.root

    vector<int>     *PDG=0;
    vector<int>     *PDGMother=0;
    vector<float>   *MCPStartPX=0;
    vector<float>   *MCPStartPY=0;
    vector<float>   *MCPStartPZ=0;
    vector<float>   *MCPStartX=0;
    vector<float>   *MCPStartY=0;
    vector<float>   *MCPStartZ=0;
    vector<float>   *TrajMCPX =0;
    vector<float>   *TrajMCPY =0;
    vector<float>   *TrajMCPZ =0;
    vector<float>   *TrajMCPT =0;
    vector<int>     *MCTrkID=0;
    vector<int>     *TrajMCPTrackID=0;


   
    TBranch *b_PDG=0;
    TBranch *b_PDGMother=0;
    TBranch *b_MCPStartPX=0;
    TBranch *b_MCPStartPY=0;
    TBranch *b_MCPStartPZ=0;
    TBranch *b_MCPStartX=0;
    TBranch *b_MCPStartY=0;
    TBranch *b_MCPStartZ=0;
    TBranch *b_TrajMCPX=0;
    TBranch *b_TrajMCPY=0;
    TBranch *b_TrajMCPZ=0;
    TBranch *b_TrajMCPT=0;
    TBranch *b_MCTrkID=0;
    TBranch *b_TrajMCPTrackID=0;



    
    chain.SetBranchAddress("PDG", &PDG, &b_PDG);
    chain.SetBranchAddress("PDGMother", &PDGMother, &b_PDGMother);
    chain.SetBranchAddress("MCPStartPX", &MCPStartPX, &b_MCPStartPX);
    chain.SetBranchAddress("MCPStartPY", &MCPStartPY, &b_MCPStartPY);
    chain.SetBranchAddress("MCPStartPZ", &MCPStartPZ, &b_MCPStartPZ);
    chain.SetBranchAddress("MCPStartX", &MCPStartX, &b_MCPStartX);
    chain.SetBranchAddress("MCPStartY", &MCPStartY, &b_MCPStartY);
    chain.SetBranchAddress("MCPStartZ", &MCPStartZ, &b_MCPStartZ);
    chain.SetBranchAddress("TrajMCPX", &TrajMCPX, &b_TrajMCPX);
    chain.SetBranchAddress("TrajMCPY", &TrajMCPY, &b_TrajMCPY);
    chain.SetBranchAddress("TrajMCPZ", &TrajMCPZ, &b_TrajMCPZ);
    chain.SetBranchAddress("TrajMCPT", &TrajMCPT, &b_TrajMCPT);
    chain.SetBranchAddress("MCTrkID", &MCTrkID, &b_MCTrkID);
    chain.SetBranchAddress("TrajMCPTrackID", &TrajMCPTrackID, &b_TrajMCPTrackID);




    chain.SetBranchStatus("*",0);
    chain.SetBranchStatus("PDG",1);
    chain.SetBranchStatus("PDGMother",1);
    chain.SetBranchStatus("MCPStartPX",1);
    chain.SetBranchStatus("MCPStartPY",1);
    chain.SetBranchStatus("MCPStartPZ",1);
    chain.SetBranchStatus("MCPStartX",1);
    chain.SetBranchStatus("MCPStartY",1);
    chain.SetBranchStatus("MCPStartZ",1);
    chain.SetBranchStatus("TrajMCPX",1);
    chain.SetBranchStatus("TrajMCPY",1);
    chain.SetBranchStatus("TrajMCPZ",1);
    chain.SetBranchStatus("TrajMCPT",1);
    chain.SetBranchStatus("TrajMCPTrackID",1);
    chain.SetBranchStatus("MCTrkID",1);

    

    Int_t nentries = (Int_t)chain.GetEntries();

    TGraph2D* g_p_theta_nmu_InLAr_mu=new TGraph2D();
    TGraph2D* g_p_theta_nmu_InLAr_mu_cut=new TGraph2D();

    TH2F* h_p_theta_nmu_InLAr_mu_q25 = new TH2F("h_p_theta_nmu_InLAr_mu_q25", "Theta VS p", 100, 0, 5, 100, 0, 180);
    TH2F* h_p_theta_nmu_InLAr_mu_q50 = new TH2F("h_p_theta_nmu_InLAr_mu_q50", "Theta VS p", 100, 0, 5, 100, 0, 180);
    TH2F* h_p_theta_nmu_InLAr_mu_q75 = new TH2F("h_p_theta_nmu_InLAr_mu_q75", "Theta VS p", 100, 0, 5, 100, 0, 180);

    
    TH2F* h_p_theta_nmu_InLAr_mu_q25_cut = new TH2F("h_p_theta_nmu_InLAr_mu_q25_cut", "Theta VS p", 100, 0, 5, 100, 0, 180);
    TH2F* h_p_theta_nmu_InLAr_mu_q50_cut = new TH2F("h_p_theta_nmu_InLAr_mu_q50_cut", "Theta VS p", 100, 0, 5, 100, 0, 180);
    TH2F* h_p_theta_nmu_InLAr_mu_q75_cut = new TH2F("h_p_theta_nmu_InLAr_mu_q75_cut", "Theta VS p", 100, 0, 5, 100, 0, 180);

    
    TH1F* hnmuons = new TH1F("hnmuons", "Muon number", 13, 1.0, 14);
    TH2F* hpnp = new TH2F("hpnp", "Portions of primamaries vs non primary muons per event", 5, 0.0, 5, 17, 0.0, 17);
    TH2F* N2_DTSpread = new TH2F("N2_DTSpread", "DTSpread at LAr Exit", 100, 0, 10000, 100, 0, 800);
    TH2F* N3_DTSpread = new TH2F("N3_DTSpread", "DTSpread at LAr Exit", 100, 0, 10000, 100, 0, 800);
    TH2F* N4_DTSpread = new TH2F("N4_DTSpread", "DTSpread at LAr Exit", 100, 0, 10000, 100, 0, 800);

    TH2F* PN2_DTSpread = new TH2F("PN2_DTSpread", "DTSpread at LAr Exit", 100, 0, 10000, 100, 0, 800);
    TH2F* PN3_DTSpread = new TH2F("PN3_DTSpread", "DTSpread at LAr Exit", 100, 0, 10000, 100, 0, 800);
    TH2F* PN4_DTSpread = new TH2F("PN4_DTSpread", "DTSpread at LAr Exit", 100, 0, 10000, 100, 0, 800);

    TH1F* hnmuons_InLAr = new TH1F("hnmuons_InLar", "Number of muons in the liquid Argon", 13, 1.0, 14);
    TH2F* hpnp_InLAr = new TH2F("hpnp_InLar", "Portions of primamaries vs non primary muons per event in the liquid Argon", 5, 0.0, 5, 17, 0.0, 17);
    
    
    
    int doublemuon=0;
    int singlemuon=0;
    int pp=0;
    int pnp=0;
    int npnp=0;
    int g=0;
    int g_cut=0;

    bool showprog = true;

    if(showprog==true) std::cout<<"Progress:  "<<std::endl;
 
    for (Int_t i=0; i<nentries; i++) 
    {
        chain.GetEntry(i);
        int muons=0;
        int muonsinLAr=0;
        int prim=0;
        int priminLAr=0;
        int close=0;
        vector<int> muonswithhits;

        int prog = 100*i/nentries;
        string strprog = std::to_string(prog);
        if(showprog==true) std::cout<<strprog<<"%";
    

        for (Int_t j=0; j<PDG->size(); j++) 
        {
           
           if(PDG->at(j)==13 || PDG->at(j)==-13)
           {
               
               
               int ID = MCTrkID->at(j);
               

               float x = MCPStartX->at(j);
               float y = MCPStartY->at(j);
               float z = MCPStartZ->at(j);
               
               

  
               muons++;
               //std::cout<<"PDGMother "<<PDGMother->at(j)<<std::endl;
               if (PDGMother->at(j)==0)
               {
                   prim++;
                   //std::cout<<"prim "<<prim<<std::endl;
                }
               //if (muons>1){std::cout<<MCPProc->at(j)<<std::endl;}

               float hit=0;

               for (Int_t k=0; k<TrajMCPX->size(); k++) 
                {
                    if (TrajMCPTrackID->at(k)==ID)
                    {
                       if(TrajMCPX->at(k)<LArUpLimit.X() && TrajMCPX->at(k)>LArLowLimit.X() && TrajMCPY->at(k)<LArUpLimit.Y() && TrajMCPY->at(k)>LArLowLimit.Y() && TrajMCPZ->at(k)<LArUpLimit.Z() && TrajMCPZ->at(k)>LArLowLimit.Z())
                       {
                           hit++;
                       }
                    }
                }

                if(hit>0)
                {
                   muonsinLAr++; 
                   muonswithhits.push_back(ID);
                } 
               
                if(hit>0 && PDGMother->at(j)==0) priminLAr++;

               
           }
           
        }
        if(muons>0) hnmuons->Fill(muons);
        if(muons>0) hpnp->Fill(prim,muons-prim);

        if(muonsinLAr>0) hnmuons_InLAr->Fill(muonsinLAr);
        if(muonsinLAr>0) hpnp_InLAr->Fill(priminLAr,muonsinLAr-priminLAr);
        
        vector<float> ydm;
        vector<float> xdm;
        vector<float> tdm;

        //std::cout<<"muonsInLAr= "<<muonsinLAr<<" muonswithhits= "<<muonswithhits.size()<<std::endl;

        for(int v=0;v<muonswithhits.size();v++)
        {
            if(muonsinLAr>1 && muonswithhits.size()>1)
            {
                int foundpoint=0;
                
                for (Int_t k=0; k<TrajMCPZ->size(); k++)
                {
                    if (TrajMCPTrackID->at(k)==muonswithhits.at(v) && TrajMCPZ->at(k)>=LArUpLimit.Z() && foundpoint==0)
                    {
                        xdm.push_back(TrajMCPX->at(k));
                        ydm.push_back(TrajMCPY->at(k)); 
                        tdm.push_back(TrajMCPT->at(k)/1e+9);
                        foundpoint++;
                    }

                }
                
            }

            for (Int_t j=0; j<PDG->size(); j++) 
            {
                if (MCTrkID->at(j)==muonswithhits.at(v))
                {
                    float px = MCPStartPX->at(j);
                    float py = MCPStartPY->at(j);
                    float pz = MCPStartPZ->at(j);
                    float p = sqrt(px*px+py*py+pz*pz);
                    float Theta = acos(pz/p)* (180.0/3.141592653589793238463);
                    g_p_theta_nmu_InLAr_mu->SetPoint(g,p,Theta,muonsinLAr);
                    g++;
                    
                }
            }

        }
        
        //std::cout<<"xdm="<<xdm.size()<<std::endl;

        vector<float> dx;
        vector<float> dy;
        vector<float> dt;
        
        for (Int_t j=0; j<xdm.size(); j++)
        {
            for (Int_t k=0; k<j; k++)
            {
              dx.push_back(abs(xdm.at(k)-xdm.at(j)));
              dy.push_back(abs(ydm.at(k)-ydm.at(j)));
              dt.push_back(abs(tdm.at(k)-tdm.at(j)));
              
              //std::cout<<"j,k,size "<<j<<" "<<k<<" "<<xdm.size()<<std::endl;
            }
        } 

        
        
        
        for (Int_t j=0; j<dx.size(); j++)
        {
            //std::cout<<dx.at(j)<<std::endl;
            if (muonsinLAr==2 )N2_DTSpread->Fill(abs(dt.at(j)),sqrt(dx.at(j)*dx.at(j)+dy.at(j)*dy.at(j)));
            if (muonsinLAr==3 )N3_DTSpread->Fill(abs(dt.at(j)),sqrt(dx.at(j)*dx.at(j)+dy.at(j)*dy.at(j)));
            if (muonsinLAr==4 )N4_DTSpread->Fill(abs(dt.at(j)),sqrt(dx.at(j)*dx.at(j)+dy.at(j)*dy.at(j)));

            if (muonsinLAr==2 )PN2_DTSpread->Fill(abs(dt.at(j)),sqrt(dx.at(j)*dx.at(j)+dy.at(j)*dy.at(j)));
            if (muonsinLAr==3 )PN3_DTSpread->Fill(abs(dt.at(j)),sqrt(dx.at(j)*dx.at(j)+dy.at(j)*dy.at(j)));
            if (muonsinLAr==4 )PN4_DTSpread->Fill(abs(dt.at(j)),sqrt(dx.at(j)*dx.at(j)+dy.at(j)*dy.at(j)));

            if(dt.at(j)<5000 && sqrt(dx.at(j)*dx.at(j)+dy.at(j)*dy.at(j))<500) close++;
            
        }
        
        //std::cout<<"close: "<<close<<std::endl;
        if (close!=0)
        {
            for(int v=0;v<muonswithhits.size();v++)
            {
                for (Int_t j=0; j<PDG->size(); j++) 
                {
                    if (MCTrkID->at(j)==muonswithhits.at(v))
                    {
                        float px = MCPStartPX->at(j);
                        float py = MCPStartPY->at(j);
                        float pz = MCPStartPZ->at(j);
                        float p = sqrt(px*px+py*py+pz*pz);
                        float Theta = acos(pz/p)* (180.0/3.141592653589793238463);
                        g_p_theta_nmu_InLAr_mu_cut->SetPoint(g_cut,p,Theta,muonsinLAr);
                        g_cut++;
                        
                    }
                }
            }
        }

        if(showprog==true) std::cout << string(strprog.length(),'\b')<<"\b";
        
 

    }

    //std::cout<<"g: "<<g<<" g_cut: "<<g_cut<<std::endl;
    ///Quantile histograms
    TH1F* h_nmu_temp_InLAr = new TH1F("h_Enu_temp_InLAr", "", 200, 0, 200);
    TH1F* h_nmu_temp_InLAr_cut = new TH1F("h_Enu_temp_InLAr_cut", "", 200, 0, 200);
    for(Int_t a=1; a<=100; a++)
    {
        for(Int_t b=1; b<=100; b++)
        {
            int n=0;
            int n_cut=0;
            for (Int_t i=0; i<g; i++) 
            {
                    Double_t gp,gTheta,gE;
                    g_p_theta_nmu_InLAr_mu->GetPoint(i,gp,gTheta,gE);
                    
                    if(gp<0.05*(a) && gp>0.05*(a-1) && gTheta<1.8*(b) && gTheta>1.8*(b-1) )
                    {
                        h_nmu_temp_InLAr->Fill(gE);
                        n++;
                    } 
            }

            for (Int_t i=0; i<g_cut; i++) 
            {
                    Double_t gp_cut,gTheta_cut,gE_cut;
                    g_p_theta_nmu_InLAr_mu_cut->GetPoint(i,gp_cut,gTheta_cut,gE_cut);
                    
                    if(gp_cut<0.05*(a) && gp_cut>0.05*(a-1) && gTheta_cut<1.8*(b) && gTheta_cut>1.8*(b-1) )
                    {
                        h_nmu_temp_InLAr_cut->Fill(gE_cut);
                        n_cut++;
                    } 
            }

            Double_t xq[3] = {0.25,0.5,0.75};
            Double_t yq_InLAr[3] = {0,0,0};
            Double_t xq_cut[3] = {0.25,0.5,0.75};
            Double_t yq_InLAr_cut[3] = {0,0,0};

            if(n!=0) h_nmu_temp_InLAr->GetQuantiles(3,yq_InLAr,xq);
            if(n_cut!=0) h_nmu_temp_InLAr_cut->GetQuantiles(3,yq_InLAr_cut,xq_cut);
                
            h_p_theta_nmu_InLAr_mu_q25->SetBinContent(a,b,yq_InLAr[0]);
            h_p_theta_nmu_InLAr_mu_q50->SetBinContent(a,b,yq_InLAr[1]);
            h_p_theta_nmu_InLAr_mu_q75->SetBinContent(a,b,yq_InLAr[2]);

            h_p_theta_nmu_InLAr_mu_q25_cut->SetBinContent(a,b,yq_InLAr_cut[0]);
            h_p_theta_nmu_InLAr_mu_q50_cut->SetBinContent(a,b,yq_InLAr_cut[1]);
            h_p_theta_nmu_InLAr_mu_q75_cut->SetBinContent(a,b,yq_InLAr_cut[2]);

            if(n!=0) h_nmu_temp_InLAr->Reset();
            if(n_cut!=0) h_nmu_temp_InLAr_cut->Reset();
        }
    }

    ////Probability histograms
    for (Int_t x=1; x<=100; x++)
    {
        float max2=0;
        float max3=0;
        float max4=0;
    
        for (Int_t y=1; y<=100; y++)
        {
            if (max2<PN2_DTSpread->GetBinContent(x,y))  max2=PN2_DTSpread->GetBinContent(x,y);
            if (max3<PN3_DTSpread->GetBinContent(x,y))  max3=PN3_DTSpread->GetBinContent(x,y);
            if (max4<PN4_DTSpread->GetBinContent(x,y))  max4=PN4_DTSpread->GetBinContent(x,y);
        }
        
        for (Int_t y=1; y<=100; y++)
        {
            float bin2=0;
            float bin3=0;
            float bin4=0;

            if (max2!=0) bin2=PN2_DTSpread->GetBinContent(x,y)/max2;
            if (max3!=0) bin3=PN3_DTSpread->GetBinContent(x,y)/max3;
            if (max4!=0) bin4=PN4_DTSpread->GetBinContent(x,y)/max4;     

            PN2_DTSpread->SetBinContent(x,y,bin2);
            PN3_DTSpread->SetBinContent(x,y,bin3);
            PN4_DTSpread->SetBinContent(x,y,bin4);
        }
        //std::cout<<"sum: "<<sum<<std::endl;
        //std::cout<<"sumcheck: "<<sumcheck<<std::endl;
    }

    TCanvas *mccanvasThetap_inLAr_nmu_q25 = new TCanvas("mccanvasThetap_inLAr_nmu_q25","",1000,800);
    gStyle->SetOptStat(0);
    h_p_theta_nmu_InLAr_mu_q25->SetTitle("Muon number in LAr 25% quantile contour;#it{p} (GeV/#it{c});#theta (deg)");
    //h_p_theta_nmu_InLAr_mu_q25->SetMaximum(4);
    h_p_theta_nmu_InLAr_mu_q25->Draw("COL Z");
    mccanvasThetap_inLAr_nmu_q25->Print("h_p_theta_nmu_InLAr_mu_q25_spill.png");

    TCanvas *mccanvasThetap_inLAr_nmu_q50 = new TCanvas("mccanvasThetap_inLAr_nmu_q50","",1000,800);
    gStyle->SetOptStat(0);
    h_p_theta_nmu_InLAr_mu_q50->SetTitle("Muon number in LAr 50% quantile contour;#it{p} (GeV/#it{c});#theta (deg)");
    //h_p_theta_nmu_InLAr_mu_q50->SetMaximum(4);
    h_p_theta_nmu_InLAr_mu_q50->Draw("COL Z");
    mccanvasThetap_inLAr_nmu_q50->Print("h_p_theta_nmu_InLAr_mu_q50_spill.png");

    TCanvas *mccanvasThetap_inLAr_nmu_q75 = new TCanvas("mccanvasThetap_inLAr_nmu_q75","",1000,800);
    gStyle->SetOptStat(0);
    h_p_theta_nmu_InLAr_mu_q75->SetTitle("Muon number in LAr 75% quantile contour;#it{p} (GeV/#it{c});#theta (deg)");
    //h_p_theta_nmu_InLAr_mu_q75->SetMaximum(4);
    h_p_theta_nmu_InLAr_mu_q75->Draw("COL Z");
    mccanvasThetap_inLAr_nmu_q75->Print("h_p_theta_nmu_InLAr_mu_q75_spill.png");




    TCanvas *mccanvasThetap_inLAr_nmu_q25_cut = new TCanvas("mccanvasThetap_inLAr_nmu_q25_cut","",1000,800);
    gStyle->SetOptStat(0);
    h_p_theta_nmu_InLAr_mu_q25_cut->SetTitle("Muon number in LAr 25% quantile contour for events having close exiting muons;#it{p} (GeV/#it{c});#theta (deg)");
    //h_p_theta_nmu_InLAr_mu_q25->SetMaximum(4);
    h_p_theta_nmu_InLAr_mu_q25_cut->Draw("COL Z");
    mccanvasThetap_inLAr_nmu_q25_cut->Print("h_p_theta_nmu_InLAr_mu_q25_spill_cut.png");

    TCanvas *mccanvasThetap_inLAr_nmu_q50_cut = new TCanvas("mccanvasThetap_inLAr_nmu_q50_cut","",1000,800);
    gStyle->SetOptStat(0);
    h_p_theta_nmu_InLAr_mu_q50_cut->SetTitle("Muon number in LAr 50% quantile contour for events having close exiting muons;#it{p} (GeV/#it{c});#theta (deg)");
    //h_p_theta_nmu_InLAr_mu_q50->SetMaximum(4);
    h_p_theta_nmu_InLAr_mu_q50_cut->Draw("COL Z");
    mccanvasThetap_inLAr_nmu_q50_cut->Print("h_p_theta_nmu_InLAr_mu_q50_spill_cut.png");

    TCanvas *mccanvasThetap_inLAr_nmu_q75_cut = new TCanvas("mccanvasThetap_inLAr_nmu_q75_cut","",1000,800);
    gStyle->SetOptStat(0);
    h_p_theta_nmu_InLAr_mu_q75_cut->SetTitle("Muon number in LAr 75% quantile contour for events having close exiting muons;#it{p} (GeV/#it{c});#theta (deg)");
    //h_p_theta_nmu_InLAr_mu_q75->SetMaximum(4);
    h_p_theta_nmu_InLAr_mu_q75_cut->Draw("COL Z");
    mccanvasThetap_inLAr_nmu_q75_cut->Print("h_p_theta_nmu_InLAr_mu_q75_spill_cut.png");




    
    TCanvas *mccanvas_DT_N2 = new TCanvas("mccanvas_DT_N2","",1000,800);
    gStyle->SetOptStat(0);
    N2_DTSpread->SetTitle("#Deltad (cm) VS #Deltat(s) at first traj point after LAr exit for multimuon events N=2;#Deltat(ns);#Deltad(cm)");
    N2_DTSpread->Draw("COL Z");
    mccanvas_DT_N2->Print("N2_DTSpread_spill.png");

    TCanvas *mccanvas_DT_N3 = new TCanvas("mccanvas_DT_N3","",1000,800);
    gStyle->SetOptStat(0);
    N3_DTSpread->SetTitle("#Deltad (cm) VS #Deltat(s) at first traj point after LAr exit for multimuon events N=3;#Deltat(ns);#Deltad(cm)");
    N3_DTSpread->Draw("COL Z");
    mccanvas_DT_N3->Print("N3_DTSpread_spill.png");

    TCanvas *mccanvas_DT_N4 = new TCanvas("mccanvas_DT_N4","",1000,800);
    gStyle->SetOptStat(0);
    N4_DTSpread->SetTitle("#Deltad (cm) VS #Deltat(s) at first traj point after LAr exit for multimuon events N=4;#Deltat(ns);#Deltad(cm)");
    N4_DTSpread->Draw("COL Z");
    mccanvas_DT_N4->Print("N4_DTSpread_spill.png");




    TCanvas *mccanvas_DT_N2P = new TCanvas("mccanvas_DT_N2P","",1000,800);
    gStyle->SetOptStat(0);
    PN2_DTSpread->SetTitle("P(#Deltad|#Deltat) at first traj point after LAr exit for multimuon events N=2;#Deltat(ns);#Deltad(cm)");
    PN2_DTSpread->Draw("COL Z");
    mccanvas_DT_N2P->Print("PN2_DTSpread_spill.png");

    TCanvas *mccanvas_DT_N3P = new TCanvas("mccanvas_DT_N3P","",1000,800);
    gStyle->SetOptStat(0);
    PN3_DTSpread->SetTitle("P(#Deltad|#Deltat) at first traj point after LAr exit for multimuon events N=3;#Deltat(ns);#Deltad(cm)");
    PN3_DTSpread->Draw("COL Z");
    mccanvas_DT_N3P->Print("PN3_DTSpread_spill.png");

    TCanvas *mccanvas_DT_N4P = new TCanvas("mccanvas_DT_N4P","",1000,800);
    gStyle->SetOptStat(0);
    PN4_DTSpread->SetTitle("P(#Deltad|#Deltat) at first traj point after LAr exit for multimuon events N=4;#Deltat(ns);#Deltad(cm)");
    PN4_DTSpread->Draw("COL Z");
    mccanvas_DT_N4P->Print("PN4_DTSpread_spill.png");




 
    gStyle->SetPaintTextFormat("3.5f");
    
    TCanvas *mccanvas_nmuons = new TCanvas("mccanvas_nmuons","",1000,800);
    hnmuons->SetTitle("Number of #mu^{+}+#mu^{-} per event;n muons per event;Probability");
    hnmuons->Scale(1/hnmuons->GetEntries());
    hnmuons->Draw();
    hnmuons->Draw("HIST TEXT0 SAME");
    mccanvas_nmuons->Print("Muons_per_event_atleast1_spill.png");

    TCanvas *mccanvas_nmuons_InLAr = new TCanvas("mccanvas_nmuons_InLAr","",1000,800);
    hnmuons_InLAr->SetTitle("Number of #mu^{+}+#mu^{-} in LAr per event;n muons per event;Probability");
    hnmuons_InLAr->Scale(1/hnmuons_InLAr->GetEntries());
    hnmuons_InLAr->Draw();
    hnmuons_InLAr->Draw("HIST TEXT0 SAME");
    mccanvas_nmuons_InLAr->Print("Muons_per_event_InLAr_atleast1_spill.png");

    TCanvas *mccanvas_pnp = new TCanvas("mccanvas_pnp","",1000,800);
    gStyle->SetOptStat(0);
    hpnp->SetTitle("Frequencies of ratios of primamaries vs non primary muons per event;Primaries;Non-Primaries");
    hpnp->Scale(1/hpnp->GetEntries());
    hpnp->Draw("COL Z");
    hpnp->Draw("TEXT0 SAME");
    mccanvas_pnp->Print("PNP_atleast1_spill.png");

    TCanvas *mccanvas_pnp_InLAr = new TCanvas("mccanvas_pnp_InLAr","",1000,800);
    gStyle->SetOptStat(0);
    hpnp_InLAr->SetTitle("Frequencies of ratios of primaries vs non primary muons in LAr per event;Primaries;Non-Primaries");
    hpnp_InLAr->Scale(1/hpnp_InLAr->GetEntries());
    hpnp_InLAr->Draw("COL Z");
    hpnp_InLAr->Draw("TEXT0 SAME");
    mccanvas_pnp_InLAr->Print("PNP_InLAr_atleast1_spill.png");

    
    
}