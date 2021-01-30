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

void l2g_traj()
{
    std::cout<<"Reading data"<<std::endl;

    TChain chain("/anatree/GArAnaTree");
    chain.Add("/mnt/c/Users/fedeb/Documents/Universita/Federico_2020-2021/OxfordCode/l2g/data/*.root");

    TVector3 GArCenter(0,-68.287,1486);
    float GAr_r = 349.9;
    float GAr_L = 669.6;
    
     
    TVector3 LArCenter(0,0,666);
    TVector3 LArDim(714.7,301.02,509.1);
    TVector3 LArUpLimit = LArCenter+0.5*LArDim;
    TVector3 LArLowLimit = LArCenter-0.5*LArDim;

    //neutrino.nd_hall_dayone_lar_SPY_v2_wMuID.volArgonCubeActive.Ev973000.Ev973999.2037.anatree.root

    vector<int>     *PDG=0;
    vector<int>     *PDGMother=0;
    vector<int>     *MCTrkID=0;
    vector<int>     *TrajMCPTrackID=0;
    vector<float>   *TrajMCPX =0;
    vector<float>   *TrajMCPY =0;
    vector<float>   *TrajMCPZ =0;
    vector<float>   *MCPStartPX=0;
    vector<float>   *MCPStartPY=0;
    vector<float>   *MCPStartPZ=0;
    vector<float>   *MCNuPx=0;
    vector<float>   *MCNuPy=0;
    vector<float>   *MCNuPz=0;
    vector<float>   *GPartE=0;

   
    TBranch *b_PDG=0;
    TBranch *b_PDGMother=0;
    TBranch *b_MCTrkID=0;
    TBranch *b_TrajMCPTrackID=0;
    TBranch *b_TrajMCPX=0;
    TBranch *b_TrajMCPY=0;
    TBranch *b_TrajMCPZ=0;
    TBranch *b_MCPStartPX=0;
    TBranch *b_MCPStartPY=0;
    TBranch *b_MCPStartPZ=0;
    TBranch *b_MCNuPx=0;
    TBranch *b_MCNuPy=0;
    TBranch *b_MCNuPz=0;
    TBranch *b_GPartE=0;


    
    chain.SetBranchAddress("PDG", &PDG, &b_PDG);
    chain.SetBranchAddress("PDGMother", &PDGMother, &b_PDGMother);
    chain.SetBranchAddress("MCTrkID", &MCTrkID, &b_MCTrkID);
    chain.SetBranchAddress("TrajMCPTrackID", &TrajMCPTrackID, &b_TrajMCPTrackID);
    //chain.SetBranchAddress("MCPProc", &MCPProc, &b_MCPProc);
    chain.SetBranchAddress("TrajMCPX", &TrajMCPX, &b_TrajMCPX);
    chain.SetBranchAddress("TrajMCPY", &TrajMCPY, &b_TrajMCPY);
    chain.SetBranchAddress("TrajMCPZ", &TrajMCPZ, &b_TrajMCPZ);
    chain.SetBranchAddress("MCPStartPX", &MCPStartPX, &b_MCPStartPX);
    chain.SetBranchAddress("MCPStartPY", &MCPStartPY, &b_MCPStartPY);
    chain.SetBranchAddress("MCPStartPZ", &MCPStartPZ, &b_MCPStartPZ);
    chain.SetBranchAddress("MCNuPx", &MCNuPx, &b_MCNuPx);
    chain.SetBranchAddress("MCNuPy", &MCNuPy, &b_MCNuPy);
    chain.SetBranchAddress("MCNuPz", &MCNuPz, &b_MCNuPz);
    chain.SetBranchAddress("GPartE", &GPartE, &b_GPartE);



    chain.SetBranchStatus("*",0);
    chain.SetBranchStatus("PDG",1);
    chain.SetBranchStatus("TrajMCPTrackID",1);
    chain.SetBranchStatus("PDGMother",1);
    chain.SetBranchStatus("MCTrkID",1);
    chain.SetBranchStatus("TrajMCPX",1);
    chain.SetBranchStatus("TrajMCPY",1);
    chain.SetBranchStatus("TrajMCPZ",1);
    chain.SetBranchStatus("MCPStartPX",1);
    chain.SetBranchStatus("MCPStartPY",1);
    chain.SetBranchStatus("MCPStartPZ",1);
    chain.SetBranchStatus("MCNuPx",1);
    chain.SetBranchStatus("MCNuPy",1);
    chain.SetBranchStatus("MCNuPz",1);
    chain.SetBranchStatus("GPartE",1);
    

    Int_t nentries = (Int_t)chain.GetEntries();  

 
    
    TH1F* hEnu = new TH1F("hEnu", "Enu", 100, 0, 100);
    TH2F* hEnu_Theta_InLAr = new TH2F("hEnu_Theta_InLAr ", "Enu VS Theta", 30, 0, 180, 30, 0, 8);
    TH2F* hP_Enu_Theta_InLAr = new TH2F("hP_Enu_Theta_InLAr", "P(Enu|Theta)", 30, 0, 180, 30, 0, 8);

    TH2F* hxyLArEnd = new TH2F("hxyLArEnd", "XY LArEnd", 100, -1000, 1000, 100, -1000, 1000);
    TH2F* hxyGArStart = new TH2F("hxyGArStart", "XY GArStart", 100, -1000, 1000, 100, -1000, 1000);

    TGraph2D* g_p_theta_Enu_InLAr_prim_mu=new TGraph2D();
    TGraph2D* g_p_theta_Enu_Crossing_prim_mu=new TGraph2D();

    TH2F* h_p_theta_Enu_InLAr_prim_mu = new TH2F("h_p_theta_Enu_InLAr_prim_mu", "Theta VS p", 100, 0, 5, 100, 0, 180);
    TH2F* h_p_theta_Enu_InLAr_prim_mu_q25 = new TH2F("h_p_theta_Enu_InLAr_prim_mu_q25", "Theta VS p", 100, 0, 5, 100, 0, 180);
    TH2F* h_p_theta_Enu_InLAr_prim_mu_q50 = new TH2F("h_p_theta_Enu_InLAr_prim_mu_q50", "Theta VS p", 100, 0, 5, 100, 0, 180);
    TH2F* h_p_theta_Enu_InLAr_prim_mu_q75 = new TH2F("h_p_theta_Enu_InLAr_prim_mu_q75", "Theta VS p", 100, 0, 5, 100, 0, 180);

    TH2F* h_p_theta_Enu_Crossing_prim_mu = new TH2F("h_p_theta_Enu_Crossing_prim_mu", "Theta VS p", 100, 0, 8, 100, 0, 180);
    TH2F* h_p_theta_Enu_Crossing_prim_mu_q25 = new TH2F("h_p_theta_Enu_Crossing_prim_mu_q25", "Theta VS p", 100, 0, 8, 100, 0, 180);
    TH2F* h_p_theta_Enu_Crossing_prim_mu_q50 = new TH2F("h_p_theta_Enu_Crossing_prim_mu_q50", "Theta VS p", 100, 0, 8, 100, 0, 180);
    TH2F* h_p_theta_Enu_Crossing_prim_mu_q75 = new TH2F("h_p_theta_Enu_Crossing_prim_mu_q75", "Theta VS p", 100, 0, 8, 100, 0, 180);

    TH1F* hbin = new TH1F("hbin", "Enu", 100, 0, 25);
    TH1* hbin_CDF;


    TH2F* h_theta_p_InLAr_prim_mu = new TH2F("h_theta_p_InLAr_prim_mu", "Theta VS p", 100, 0, 5, 100, 0, 180);
    TH2F* h_theta_p_Crossing_prim_mu = new TH2F("h_theta_p_Crossing_prim_mu", "Theta VS p", 100, 0, 8, 100, 0, 180);
    TH2F* h_theta_p_All_prim_mu = new TH2F("h_theta_p_All_prim_mu", "Theta VS p", 100, 0, 8, 100, 0, 180);
    
    
    
    int g1=0;
    int g2=0;
    bool showprog = true;

    if(showprog==true) std::cout<<"Progress:  "<<std::endl;

    for (Int_t i=0; i<nentries; i++) 
    {
        chain.GetEntry(i);
        float pxNu= MCNuPx->at(0);
        float pyNu= MCNuPy->at(0);
        float pzNu= MCNuPz->at(0);
        float Enu = sqrt(pxNu*pxNu+pyNu*pyNu+pzNu*pzNu);
        float Enu_check= GPartE->at(0);
        float ThetaNu = acos(pzNu/Enu);
        float ThetaNu_deg = ThetaNu * (180.0/3.141592653589793238463);

        //std::cout<<"pxNu= "<<pxNu<<std::endl;
        //std::cout<<"pyNu= "<<pyNu<<std::endl;
        //std::cout<<"pzNu= "<<pzNu<<std::endl;
        int prog = 100*i/nentries;
        string strprog = std::to_string(prog);
        if(showprog==true) std::cout<<strprog<<"%"; 
        hEnu->Fill(Enu);
        int prim=0;
        

          
 

        for (Int_t j=0; j<PDG->size(); j++) 
        {
           if(PDG->at(j)==13 || PDG->at(j)==-13)
           {
               int ID = MCTrkID->at(j);
               int PDGMom = PDGMother->at(j);
               float px = MCPStartPX->at(j);
               float py = MCPStartPY->at(j);
               float pz = MCPStartPZ->at(j);
               float p = sqrt(px*px+py*py+pz*pz);
               float Theta = acos(pz/p);
               float Theta_deg = Theta * (180.0/3.141592653589793238463);
               //if(PDGMother->at(j)==0)prim++;
               
               

               int no1=0;
               int no2=0;
               int OutLAr=0;
               int InLAr=0;
               int InGAr=0;

               for (Int_t k=0; k<TrajMCPX->size(); k++) 
                {
                    if (TrajMCPTrackID->at(k)==ID)
                    {
                        //Find actual limits LAr and GAr
                        if(TrajMCPX->at(k)<LArUpLimit.X() && TrajMCPX->at(k)>LArLowLimit.X() && TrajMCPY->at(k)<LArUpLimit.Y() && TrajMCPY->at(k)>LArLowLimit.Y() && TrajMCPZ->at(k)<LArUpLimit.Z() && TrajMCPZ->at(k)>LArLowLimit.Z())
                        {
                          InLAr++;
                        }

                        if(TrajMCPX->at(k)<LArLowLimit.X() || TrajMCPX->at(k)>LArUpLimit.X() || TrajMCPY->at(k)<LArLowLimit.Y() || TrajMCPY->at(k)>LArUpLimit.Y() || TrajMCPZ->at(k)<LArLowLimit.Z() || TrajMCPZ->at(k)>LArUpLimit.Z())
                        {
                          OutLAr++;
                        }
                        
                        float r = sqrt((TrajMCPY->at(k)-GArCenter.Y())*(TrajMCPY->at(k)-GArCenter.Y())+(TrajMCPZ->at(k)-GArCenter.Z())*(TrajMCPZ->at(k)-GArCenter.Z()));

                        if(TrajMCPX->at(k)<0.5*GAr_L && TrajMCPX->at(k)>-0.5*GAr_L && r<GAr_r)
                        {
                          InGAr++;
                        }
                        
                        //Propagation
                        if(TrajMCPZ->at(k)<950 && TrajMCPZ->at(k)>850  && no1==0 && PDGMother->at(j)==0  && PDG->at(j)==13)
                        {
                            hxyLArEnd->Fill(TrajMCPX->at(k),TrajMCPY->at(k));
                            no1++;
                        }

                        if(TrajMCPZ->at(k)<1200 && TrajMCPZ->at(k)>1100  && no2==0 && PDGMother->at(j)==0  && PDG->at(j)==13)
                        {
                            hxyGArStart->Fill(TrajMCPX->at(k),TrajMCPY->at(k));
                            no2++;
                        }
                    }
                }

                if(OutLAr==0 && InLAr!=0 && PDG->at(j)==13 && PDGMother->at(j)==0)
                    {
                        h_theta_p_InLAr_prim_mu->Fill(p,Theta_deg);
                        if(p<5 && Theta_deg<180 ){
                        g_p_theta_Enu_InLAr_prim_mu->SetPoint(g1,p,Theta_deg,Enu);
                        hEnu_Theta_InLAr->Fill(Theta_deg,Enu);
                        hP_Enu_Theta_InLAr->Fill(Theta_deg,Enu);
                        g1++;
                        }
                    }
                if(OutLAr!=0 && InGAr!=0 && InLAr!=0 && PDG->at(j)==13 && PDGMother->at(j)==0)
                    {
                        h_theta_p_Crossing_prim_mu->Fill(p,Theta_deg);
                        if(p<8 && Theta_deg<180){
                            g_p_theta_Enu_Crossing_prim_mu->SetPoint(g2,p,Theta_deg,Enu);
                            g2++;
                        }
                    }
                if(PDG->at(j)==13 && PDGMother->at(j)==0 )
                    {
                         h_theta_p_All_prim_mu->Fill(p,Theta_deg);
                    }

                
            }
   
        }
        
        
        if(showprog==true) std::cout << string(strprog.length(),'\b')<<"\b";   

    }

    showprog = true;
    if(showprog==true) std::cout<<std::endl<<"Computing probability plot  "<<std::endl;

    for (Int_t x=1; x<=100; x++)
    {
        float max=0;
    
        for (Int_t y=1; y<=100; y++)
        {
            if (max<hP_Enu_Theta_InLAr->GetBinContent(x,y))  max=hP_Enu_Theta_InLAr->GetBinContent(x,y);
        }
        
        for (Int_t y=1; y<=100; y++)
        {
            float bin=0;

            if (max!=0) bin=hP_Enu_Theta_InLAr->GetBinContent(x,y)/max;     

            hP_Enu_Theta_InLAr->SetBinContent(x,y,bin);
        }
        //std::cout<<"sum: "<<sum<<std::endl;
        //std::cout<<"sumcheck: "<<sumcheck<<std::endl;
    }
    
    

    showprog = true;
    if(showprog==true) std::cout<<std::endl<<"Binning histogram progress:  "<<std::endl;

    TH1F* h_Enu_temp_InLAr = new TH1F("h_Enu_temp_InLAr", "", 200, 0, 200);
    TH1F* h_Enu_temp_Crossing = new TH1F("h_Enu_temp_Crossing", "", 200, 0, 200);
    
    for(Int_t a=1; a<=100; a++)
    {
        int prog = a;
        string strprog = std::to_string(prog);
        if(showprog==true) std::cout<<strprog<<"%";

        for(Int_t b=1; b<=100; b++)
        {
                float sum=0.;
                float sum_c=0.;
                int n=0;
                int n_c=0;

                for (Int_t i=0; i<g1; i++) 
                {
                    Double_t gp,gTheta,gE;
                    g_p_theta_Enu_InLAr_prim_mu->GetPoint(i,gp,gTheta,gE);
                    
                    if(gp<0.05*(a) && gp>0.05*(a-1) && gTheta<1.8*(b) && gTheta>1.8*(b-1) )
                    {
                        sum+=gE;
                        n++;
                        h_Enu_temp_InLAr->Fill(gE);
                    } 
                  
                }

                for (Int_t i=0; i<g2; i++) 
                {
                    Double_t gp,gTheta,gE;
                    g_p_theta_Enu_Crossing_prim_mu->GetPoint(i,gp,gTheta,gE);
                    if(gp<0.08*(a) && gp>0.08*(a-1) && gTheta<1.8*(b) && gTheta>1.8*(b-1))
                    {
                        sum_c+=gE;
                        n_c++;
                        h_Enu_temp_Crossing->Fill(gE);
                        //if(a==20 && b==20) hbin->Fill(gE);
                    } 
                  
                }
                if(n!=0){sum/=n;}
                if(n_c!=0){sum_c/=n_c;}
                //Int_t bin1=h_p_theta_Enu_InLAr_prim_mu->GetBin(a,b);
                //Int_t bin2=h_p_theta_Enu_Crossing_prim_mu->GetBin(a,b);
                h_p_theta_Enu_InLAr_prim_mu->SetBinContent(a,b,sum);
                h_p_theta_Enu_Crossing_prim_mu->SetBinContent(a,b,sum_c);

                Double_t xq[3] = {0.25,0.5,0.75};
                Double_t yq_InLAr[3] = {0,0,0};
                Double_t yq_Crossing[3] = {0,0,0};

                if(n!=0 && sum!=0)
                {
                    h_Enu_temp_InLAr->GetQuantiles(3,yq_InLAr,xq);
                }
                //bin1=h_p_theta_Enu_InLAr_prim_mu_q25->GetBin(a,b);
                h_p_theta_Enu_InLAr_prim_mu_q25->SetBinContent(a,b,yq_InLAr[0]);
                //bin1=h_p_theta_Enu_InLAr_prim_mu_q50->GetBin(a,b);
                h_p_theta_Enu_InLAr_prim_mu_q50->SetBinContent(a,b,yq_InLAr[1]);
                //bin1=h_p_theta_Enu_InLAr_prim_mu_q75->GetBin(a,b);
                h_p_theta_Enu_InLAr_prim_mu_q75->SetBinContent(a,b,yq_InLAr[2]);
                //std::cout<<"debug1"<<std::endl;
                
                
                if(n_c!=0 && sum_c!=0)
                {
                    h_Enu_temp_Crossing->GetQuantiles(3,yq_Crossing,xq);
                }

                    //bin2=h_p_theta_Enu_Crossing_prim_mu_q25->GetBin(a,b);
                    h_p_theta_Enu_Crossing_prim_mu_q25->SetBinContent(a,b,yq_Crossing[0]); 
                    //bin2=h_p_theta_Enu_Crossing_prim_mu_q50->GetBin(a,b);
                    h_p_theta_Enu_Crossing_prim_mu_q50->SetBinContent(a,b,yq_Crossing[1]);              
                    //bin2=h_p_theta_Enu_Crossing_prim_mu_q75->GetBin(a,b);
                    h_p_theta_Enu_Crossing_prim_mu_q75->SetBinContent(a,b,yq_Crossing[2]);
                //std::cout<<"debug2"<<std::endl;
                

                //if(sum!=0 || sum_c!=0) std::cout<<"sum= "<<sum<<" quantile= "<<yq_InLAr[2]<<std::endl;
                //if(sum!=0 || sum_c!=0) std::cout<<"sum_c= "<<sum_c<<" quantile= "<<yq_Crossing[2]<<std::endl<<std::endl;

                
                
                if(n!=0 && sum!=0) h_Enu_temp_InLAr->Reset();
                if(n_c!=0 && sum_c!=0) h_Enu_temp_Crossing->Reset();
                //std::cout<<"debug3"<<std::endl;


                // std::cout<<"sum= "<<sum<<" sum_c= "<<sum_c<<std::endl;
                //if(sum!=0 || sum_c!=0) std::cout<<"bin1= "<<h_p_theta_Enu_InLAr_prim_mu->GetBinContent(bin1)<<" bin2= "<<h_p_theta_Enu_Crossing_prim_mu->GetBinContent(bin2)<<std::endl;
               
        }

        if(showprog==true) std::cout << string(strprog.length(),'\b')<<"\b";
    }
    

    
    
     
    
    //gStyle->SetPalette(104);
    
    TCanvas *mccanvasEnu = new TCanvas("mccanvasEnu","",1000,800);
    hEnu->SetTitle("Neutrino energy spectrum ;E_{#nu};neutrinos");
    hEnu->Draw();
    mccanvasEnu->Print("Enu.png");

    /*
    TCanvas *mccanvasbin = new TCanvas("mccanvasbin","",1000,800);
    hbin->SetTitle("Neutrino energy spectrum in one (p,#theta) bin;E_{#nu};neutrinos");
    hbin->Draw();
    mccanvasbin->Print("Enu_bin.png");

    hbin->Scale(1/hbin->Integral());
    hbin_CDF=hbin->GetCumulative();
    

    TCanvas *mccanvasbin_CDF = new TCanvas("mccanvasbin_CDF","",1000,800);
    hbin_CDF->SetTitle("Neutrino energy spectrum CDF in one (p,#theta) bin;E_{#nu};neutrinos");
    hbin_CDF->Draw();
    mccanvasbin_CDF->Print("Enu_bin_CDF.png");
    */
    TCanvas *mccanvasThetaVSEnu_InLAr = new TCanvas("mccanvasThetaVSEnu_InLAr","",1000,800);
    gStyle->SetOptStat(0);
    hEnu_Theta_InLAr->SetTitle("neutrino energy E_{#nu} VS Muon angle #theta;#theta (deg);E_{#nu} (GeV)");
    hEnu_Theta_InLAr->Draw("COLZ");
    mccanvasThetaVSEnu_InLAr->Print("ThetaVSEnu_InLAr.png");

    TCanvas *mccanvasEnuVSThetanu_InLAr = new TCanvas("mccanvasEnuVSThetanu_InLAr","",1000,800);
    gStyle->SetOptStat(0);
    hP_Enu_Theta_InLAr->SetTitle("P(E_{#nu}|#theta);#theta (deg);E_{#nu} (GeV)");
    hP_Enu_Theta_InLAr->Draw("COLZ");
    mccanvasEnuVSThetanu_InLAr->Print("P_EnuVSTheta_InLAr.png");
    

    
    
    TCanvas *mccanvasThetap_inLAr_Enu = new TCanvas("mccanvasThetap_inLAr_Enu","",1000,800);
    gStyle->SetOptStat(0);
    h_p_theta_Enu_InLAr_prim_mu->SetTitle("Primary #mu^{-} contained in LAr alone E_{#nu} [GeV] contour;#it{p} (GeV/#it{c});#theta (deg)");
    h_p_theta_Enu_InLAr_prim_mu->SetMaximum(4);
    h_p_theta_Enu_InLAr_prim_mu->Draw("COL Z");
    mccanvasThetap_inLAr_Enu->Print("h_p_theta_Enu_InLAr_prim_mu_max.png");

    TCanvas *mccanvasThetap_inLAr_Enu_q25 = new TCanvas("mccanvasThetap_inLAr_Enu_q25","",1000,800);
    gStyle->SetOptStat(0);
    h_p_theta_Enu_InLAr_prim_mu_q25->SetTitle("Primary #mu^{-} contained in LAr alone E_{#nu} 25% quantile [GeV] contour;#it{p} (GeV/#it{c});#theta (deg)");
    h_p_theta_Enu_InLAr_prim_mu_q25->SetMaximum(4);
    h_p_theta_Enu_InLAr_prim_mu_q25->Draw("COL Z");
    mccanvasThetap_inLAr_Enu_q25->Print("h_p_theta_Enu_InLAr_prim_mu_max_q25.png");

    TCanvas *mccanvasThetap_inLAr_Enu_q50 = new TCanvas("mccanvasThetap_inLAr_Enu_q50","",1000,800);
    gStyle->SetOptStat(0);
    h_p_theta_Enu_InLAr_prim_mu_q50->SetTitle("Primary #mu^{-} contained in LAr alone E_{#nu} 50% quantile [GeV] contour;#it{p} (GeV/#it{c});#theta (deg)");
    h_p_theta_Enu_InLAr_prim_mu_q50->SetMaximum(4);
    h_p_theta_Enu_InLAr_prim_mu_q50->Draw("COL Z");
    mccanvasThetap_inLAr_Enu_q50->Print("h_p_theta_Enu_InLAr_prim_mu_max_q50.png");

    TCanvas *mccanvasThetap_inLAr_Enu_q75 = new TCanvas("mccanvasThetap_inLAr_Enu_q75","",1000,800);
    gStyle->SetOptStat(0);
    h_p_theta_Enu_InLAr_prim_mu_q75->SetTitle("Primary #mu^{-} contained in LAr alone E_{#nu} 75% quantile [GeV] contour;#it{p} (GeV/#it{c});#theta (deg)");
    h_p_theta_Enu_InLAr_prim_mu_q75->SetMaximum(4);
    h_p_theta_Enu_InLAr_prim_mu_q75->Draw("COL Z");
    mccanvasThetap_inLAr_Enu_q75->Print("h_p_theta_Enu_InLAr_prim_mu_max_q75.png");


    
    TCanvas *mccanvasThetap_Crossing_Enu = new TCanvas("mccanvasThetap_Crossing_Enu","",1000,800);
    gStyle->SetOptStat(0);
    h_p_theta_Enu_Crossing_prim_mu->SetTitle("Primary crossing LAr+TMS #mu^{-} E_{#nu} [GeV] contour;#it{p} (GeV/#it{c});#theta (deg)");
    h_p_theta_Enu_Crossing_prim_mu->SetMaximum(5);
    h_p_theta_Enu_Crossing_prim_mu->Draw("COL Z");
    mccanvasThetap_Crossing_Enu->Print("h_p_theta_Enu_Crossing_prim_mu_max.png");

    TCanvas *mccanvasThetap_Crossing_Enu_q25 = new TCanvas("mccanvasThetap_Crossing_Enu_q25","",1000,800);
    gStyle->SetOptStat(0);
    h_p_theta_Enu_Crossing_prim_mu_q25->SetTitle("Primary crossing LAr+TMS #mu^{-} E_{#nu} 25% quantile [GeV] contour;#it{p} (GeV/#it{c});#theta (deg)");
    h_p_theta_Enu_Crossing_prim_mu_q25->SetMaximum(5);
    h_p_theta_Enu_Crossing_prim_mu_q25->Draw("COL Z");
    mccanvasThetap_Crossing_Enu_q25->Print("h_p_theta_Enu_Crossing_prim_mu_max_q25.png");

    TCanvas *mccanvasThetap_Crossing_Enu_q50 = new TCanvas("mccanvasThetap_Crossing_Enu_q50","",1000,800);
    gStyle->SetOptStat(0);
    h_p_theta_Enu_Crossing_prim_mu_q50->SetTitle("Primary crossing LAr+TMS #mu^{-} E_{#nu} 50% quantile [GeV] contour;#it{p} (GeV/#it{c});#theta (deg)");
    h_p_theta_Enu_Crossing_prim_mu_q50->SetMaximum(5);
    h_p_theta_Enu_Crossing_prim_mu_q50->Draw("COL Z");
    mccanvasThetap_Crossing_Enu_q50->Print("h_p_theta_Enu_Crossing_prim_mu_max_q50.png");

    TCanvas *mccanvasThetap_Crossing_Enu_q75 = new TCanvas("mccanvasThetap_Crossing_Enu_q75","",1000,800);
    gStyle->SetOptStat(0);
    h_p_theta_Enu_Crossing_prim_mu_q75->SetTitle("Primary crossing LAr+TMS #mu^{-} E_{#nu} 75% quantile [GeV] contour;#it{p} (GeV/#it{c});#theta (deg)");
    h_p_theta_Enu_Crossing_prim_mu_q75->SetMaximum(5);
    h_p_theta_Enu_Crossing_prim_mu_q75->Draw("COL Z");
    mccanvasThetap_Crossing_Enu_q75->Print("h_p_theta_Enu_Crossing_prim_mu_max_q75.png");


    






    TCanvas *mccanvasThetap_inLAr = new TCanvas("mccanvasThetap_inLAr","",1000,800);
    gStyle->SetOptStat(1);
    h_theta_p_InLAr_prim_mu->SetTitle("Primary muons contained in LAr alone #theta VS #it{p};#it{p} (GeV/#it{c});#theta (deg)");
    h_theta_p_InLAr_prim_mu->Draw("COLZ");
    mccanvasThetap_inLAr->Print("ThetaVSp_inLAr_prim_mu.png");

    TCanvas *mccanvasThetap_Crossing = new TCanvas("mccanvasThetap_Crossing","",1000,800);
    gStyle->SetOptStat(1);
    h_theta_p_Crossing_prim_mu->SetTitle("Primary muons crossing LAr+TMS #theta VS #it{p} ;#it{p} (GeV/#it{c});#theta (deg)");
    h_theta_p_Crossing_prim_mu->Draw("COLZ");
    mccanvasThetap_Crossing->Print("ThetaVSp_inGAr_prim_mu.png");

    TCanvas *mccanvasThetap_All = new TCanvas("mccanvasThetap_All","",1000,800);
    gStyle->SetOptStat(1);
    h_theta_p_All_prim_mu->SetTitle("All primary muons #theta VS #it{p} ;#it{p} (GeV/#it{c});#theta (deg)");
    h_theta_p_All_prim_mu->Draw("COLZ");
    mccanvasThetap_All->Print("ThetaVSp_All_prim_mu.png");

    

    /*
    TH2F *hxyLArGArRatio = (TH2F*) hxyGArStart->Clone("hxyLArGArRatio");
    hxyLArGArRatio->Divide(hxyLArEnd);

    TCanvas *mccanvasRatio = new TCanvas("mccanvasRatio","",1000,800);
    gStyle->SetOptStat(0);
    hxyLArGArRatio->SetTitle("Muon ratio GAr entrance/LAr exit;x[cm];y[cm]");
    hxyLArGArRatio->Draw("COLZ");
    mccanvasRatio->Print("14_LArGArRatio_prim.png");
    */

    
 
}