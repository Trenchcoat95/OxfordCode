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

Double_t CauchyDens(Double_t *x, Double_t *par)
{
   Double_t pi   = TMath::Pi();
   Double_t mean = par[0];
   Double_t fwhm = par[1];

   Double_t arg = x[0]-mean;
   Double_t top = 0.5*fwhm;
   Double_t bot = pi*(arg*arg+top*top);

   Double_t func = top/bot;
   return func;
}

Double_t CauchyPeak(Double_t *x, Double_t *par)
{
   Double_t height = par[2];
   Double_t func = height*CauchyDens(x,par);
   return func;
}

void l2g_trackmatch_test()
{
    TChain chain("/anatree/GArAnaTree");
    chain.Add("/home/federico/Documents/Universita/Federico_2020-2021/OxfordCode/l2g/data/neutrino.nd_hall_dayone_lar_SPY_v2_wMuID.volArgonCubeActive.Ev973000.Ev973999.2037.anatree.root");

    TVector3 GArCenter(0,-68.287,1486);
    float GAr_r = 349.9;
    float GAr_L = 669.6;
    //neutrino.nd_hall_dayone_lar_SPY_v2_wMuID.volArgonCubeActive.Ev973000.Ev973999.2037.anatree.root

    vector<int>     *PDG=0;
    vector<int>     *PDGMother=0;
    vector<float>   *MCNuPx=0;
    vector<float>   *MCNuPy=0;
    vector<float>   *MCNuPz=0;
    vector<int>     *MCTrkID=0;
    vector<int>     *TrajMCPTrackID=0;
    vector<float>   *MCPStartPX=0;
    vector<float>   *MCPStartPY=0;
    vector<float>   *MCPStartPZ=0;
    vector<float>   *MCPStartX=0;
    vector<float>   *MCPStartY=0;
    vector<float>   *MCPStartZ=0;
    vector<float>   *TrajMCPX =0;
    vector<float>   *TrajMCPY =0;
    vector<float>   *TrajMCPZ =0;
    vector<float>   *TrajMCPPX =0;
    vector<float>   *TrajMCPPY =0;
    vector<float>   *TrajMCPPZ =0;
    vector<int>     *NTPCClustersOnTrack=0;
    vector<float>     *TrackStartX = 0;
    vector<float>     *TrackStartY = 0;
    vector<float>     *TrackStartZ = 0;
    vector<float>     *TrackStartPX = 0;
    vector<float>     *TrackStartPY = 0;
    vector<float>     *TrackStartPZ = 0;
    vector<float>     *TrackEndX = 0;
    vector<float>     *TrackEndY = 0;
    vector<float>     *TrackEndZ = 0;
    vector<float>     *TrackEndPX = 0;
    vector<float>     *TrackEndPY = 0;
    vector<float>     *TrackEndPZ = 0;
    vector<int>     *TrackStartQ = 0;
    vector<int>     *TrackEndQ = 0;
    //vector<ULong_64t> *TrackIDNumber=0;
    //vector<string>  *MCPProc =0;
   
    TBranch *b_PDG=0;
    TBranch *b_PDGMother=0;
    TBranch *b_MCNuPx=0;
    TBranch *b_MCNuPy=0;
    TBranch *b_MCNuPz=0;
    TBranch *b_MCTrkID=0;
    TBranch *b_TrajMCPTrackID=0;
    TBranch *b_MCPStartPX=0;
    TBranch *b_MCPStartPY=0;
    TBranch *b_MCPStartPZ=0;
    TBranch *b_MCPStartX=0;
    TBranch *b_MCPStartY=0;
    TBranch *b_MCPStartZ=0;
    TBranch *b_TrajMCPPX=0;
    TBranch *b_TrajMCPPY=0;
    TBranch *b_TrajMCPPZ=0;
    TBranch *b_TrajMCPX=0;
    TBranch *b_TrajMCPY=0;
    TBranch *b_TrajMCPZ=0;
    TBranch *b_NTPCClustersOnTrack=0;
    TBranch *b_TrackStartX=0;
    TBranch *b_TrackStartY=0;
    TBranch *b_TrackStartZ=0;
    TBranch *b_TrackStartPX=0;
    TBranch *b_TrackStartPY=0;
    TBranch *b_TrackStartPZ=0;
    TBranch *b_TrackEndX=0;
    TBranch *b_TrackEndY=0;
    TBranch *b_TrackEndZ=0;
    TBranch *b_TrackEndPX=0;
    TBranch *b_TrackEndPY=0;
    TBranch *b_TrackEndPZ=0;
    TBranch *b_TrackStartQ=0;
    TBranch *b_TrackEndQ=0;
    //TBranch *b_TrackIDNumber =0;
    //TBranch *b_MCPProc=0;

    
    chain.SetBranchAddress("PDG", &PDG, &b_PDG);
    chain.SetBranchAddress("PDGMother", &PDGMother, &b_PDGMother);
    chain.SetBranchAddress("MCNuPx", &MCNuPx, &b_MCNuPx);
    chain.SetBranchAddress("MCNuPy", &MCNuPy, &b_MCNuPy);
    chain.SetBranchAddress("MCNuPz", &MCNuPz, &b_MCNuPz);
    chain.SetBranchAddress("MCTrkID", &MCTrkID, &b_MCTrkID);
    chain.SetBranchAddress("TrajMCPTrackID", &TrajMCPTrackID, &b_TrajMCPTrackID);
    //chain.SetBranchAddress("MCPProc", &MCPProc, &b_MCPProc);
    chain.SetBranchAddress("MCPStartPX", &MCPStartPX, &b_MCPStartPX);
    chain.SetBranchAddress("MCPStartPY", &MCPStartPY, &b_MCPStartPY);
    chain.SetBranchAddress("MCPStartPZ", &MCPStartPZ, &b_MCPStartPZ);
    chain.SetBranchAddress("MCPStartX", &MCPStartX, &b_MCPStartX);
    chain.SetBranchAddress("MCPStartY", &MCPStartY, &b_MCPStartY);
    chain.SetBranchAddress("MCPStartZ", &MCPStartZ, &b_MCPStartZ);
    chain.SetBranchAddress("TrajMCPX", &TrajMCPX, &b_TrajMCPX);
    chain.SetBranchAddress("TrajMCPY", &TrajMCPY, &b_TrajMCPY);
    chain.SetBranchAddress("TrajMCPZ", &TrajMCPZ, &b_TrajMCPZ);
    chain.SetBranchAddress("TrajMCPPX", &TrajMCPPX, &b_TrajMCPPX);
    chain.SetBranchAddress("TrajMCPPY", &TrajMCPPY, &b_TrajMCPPY);
    chain.SetBranchAddress("TrajMCPPZ", &TrajMCPPZ, &b_TrajMCPPZ);
    chain.SetBranchAddress("NTPCClustersOnTrack", &NTPCClustersOnTrack, &b_NTPCClustersOnTrack);
    chain.SetBranchAddress("TrackStartX", &TrackStartX, &b_TrackStartX);
    chain.SetBranchAddress("TrackStartY", &TrackStartY, &b_TrackStartY);
    chain.SetBranchAddress("TrackStartZ", &TrackStartZ, &b_TrackStartZ);
    chain.SetBranchAddress("TrackStartPX", &TrackStartPX, &b_TrackStartPX);
    chain.SetBranchAddress("TrackStartPY", &TrackStartPY, &b_TrackStartPY);
    chain.SetBranchAddress("TrackStartPZ", &TrackStartPZ, &b_TrackStartPZ); 
    chain.SetBranchAddress("TrackEndX", &TrackEndX, &b_TrackEndX);
    chain.SetBranchAddress("TrackEndY", &TrackEndY, &b_TrackEndY);
    chain.SetBranchAddress("TrackEndZ", &TrackEndZ, &b_TrackEndZ);
    chain.SetBranchAddress("TrackEndPX", &TrackEndPX, &b_TrackEndPX);
    chain.SetBranchAddress("TrackEndPY", &TrackEndPY, &b_TrackEndPY);
    chain.SetBranchAddress("TrackEndPZ", &TrackEndPZ, &b_TrackEndPZ);   
    chain.SetBranchAddress("TrackStartQ", &TrackStartQ, &b_TrackEndQ);
    chain.SetBranchAddress("TrackEndQ", &TrackEndQ, &b_TrackEndQ);
    //chain.SetBranchAddress("TrackIDNumber", &TrackIDNumber, &b_TrackIDNumber);


    chain.SetBranchStatus("*",0);
    chain.SetBranchStatus("PDG",1);
    chain.SetBranchStatus("PDGMother",1);
    chain.SetBranchStatus("MCNuPx",1);
    chain.SetBranchStatus("MCNuPy",1);
    chain.SetBranchStatus("MCNuPz",1);
    chain.SetBranchStatus("MCTrkID",1);
    chain.SetBranchStatus("TrajMCPTrackID",1);
    chain.SetBranchStatus("MCPStartPX",1);
    chain.SetBranchStatus("MCPStartPY",1);
    chain.SetBranchStatus("MCPStartPZ",1);
    chain.SetBranchStatus("MCPStartX",1);
    chain.SetBranchStatus("MCPStartY",1);
    chain.SetBranchStatus("MCPStartZ",1);
    chain.SetBranchStatus("TrajMCPX",1);
    chain.SetBranchStatus("TrajMCPY",1);
    chain.SetBranchStatus("TrajMCPZ",1);
    chain.SetBranchStatus("TrajMCPPX",1);
    chain.SetBranchStatus("TrajMCPPY",1);
    chain.SetBranchStatus("TrajMCPPZ",1);
    chain.SetBranchStatus("NTPCClustersOnTrack",1);
    chain.SetBranchStatus("TrackStartX",1);
    chain.SetBranchStatus("TrackStartY",1);
    chain.SetBranchStatus("TrackStartZ",1);
    chain.SetBranchStatus("TrackStartPX",1);
    chain.SetBranchStatus("TrackStartPY",1);
    chain.SetBranchStatus("TrackStartPZ",1);
    chain.SetBranchStatus("TrackEndX",1);
    chain.SetBranchStatus("TrackEndY",1);
    chain.SetBranchStatus("TrackEndZ",1);
    chain.SetBranchStatus("TrackEndPX",1);
    chain.SetBranchStatus("TrackEndPY",1);
    chain.SetBranchStatus("TrackEndPZ",1);
    chain.SetBranchStatus("TrackStartQ",1);
    chain.SetBranchStatus("TrackEndQ",1);
    //chain.SetBranchStatus("TrackIDNumber",1);

    Int_t nentries = (Int_t)chain.GetEntries();
    TH1F* hnhits = new TH1F("hnhits", "hnhits", 20, 0, 20);
    TH1F* hnhitsmatched = new TH1F("hnhitsmatched", "hnhitsmatched", 20, 0, 20);
    TH1F* hz = new TH1F("hz", "hz", 50, 0, 20);
    TH1F* hcos = new TH1F("hcos", "hcos", 50, 0.96, 1);
    TH1F* hdelX = new TH1F("hdelX", "hdelX", 50, 0, 50);

    TH1F* hmu_good_recocharge_p_start = new TH1F("hmu_good_recocharge_p_start", "Reconstructed charge resolution", 50, 0, 8.0);
    TH1F* hmu_all_recocharge_p_start = new TH1F("hmu_all_recocharge_p_start", "Reconstructed charge resolution", 50, 0, 8.0);

    TH1F* hmu_good_reco_thetanu = new TH1F("hmu_good_reco_thetanu", "Reconstruction efficiency", 180, 0, 180.0);
    TH1F* hmu_all_thetanu = new TH1F("hmu_all_thetanu", "Reconstruction efficiency", 180, 0, 180.0);
    TH1F* hmu_all_thetanu_InGAr = new TH1F("hmu_all_thetanu_InGAr", "Reconstruction efficiency", 180, 0, 180.0);
    
    TH1F* hmu_good_reco_pmu = new TH1F("hmu_good_reco_pmu", "Reconstruction efficiency", 100, 0, 20.0);
    TH1F* hmu_all_pmu = new TH1F("hmu_all_pmu", "Reconstruction efficiency", 100, 0, 20.0);
    TH1F* hmu_all_pmu_InGAr = new TH1F("hmu_all_pmu_InGAr", "Reconstruction efficiency", 100, 0, 20.0);
    
    

    TH1F* frac_resid = new TH1F("frac_resid", "Fractional residuals", 150, -0.4, 0.4);
    TH1F* frac_resid_c = new TH1F("frac_resid_c", "Fractional residuals", 150, -0.4, 0.4);
    TH1F* ptrue = new TH1F("ptrue", "ptrue", 50, 0, 8);
    TH2F* pVSfrac_resid = new TH2F("pVSfrac_resid", "Fractional residuals",50,0,8, 50, -1, 1);
    TH2F* PpVSfrac_resid = new TH2F("PpVSfrac_resid", "P(Residual|p)",50,0,8, 50, -1, 1);
    

    bool showprog = true;

    if(showprog==true) std::cout<<"Progress:  "<<std::endl;

    for (Int_t i=0; i<nentries; i++) 
    {
        chain.GetEntry(i);
        int prog = 100*i/nentries;
        string strprog = std::to_string(prog);
        if(showprog==true) std::cout<<strprog<<"%";
        
        float pxNu= MCNuPx->at(0);
        float pyNu= MCNuPy->at(0);
        float pzNu= MCNuPz->at(0);
        TVector3 MCNu(pxNu,pyNu,pzNu);
        int numucc=0;
        int numuccInGAr=0;
        int recoInGAr=0;
        //float pNu = sqrt(pxNu*pxNu+pyNu*pyNu+pzNu*pzNu);
        //float ThetaNu = acos(pzNu/pNu);
        //float ThetaNu_deg = ThetaNu * (180.0/3.141592653589793238463);

        for (Int_t j=0; j<PDG->size(); j++) 
        {
           if(PDG->at(j)==13 || PDG->at(j)==-13)
           {
               float px = 0;
               float py = 0;
               float pz = 0;
               float p=0;
               float x = 0;
               float y = 0;
               float z = 0;
               float pxSt = MCPStartPX->at(j);
               float pySt = MCPStartPY->at(j);
               float pzSt = MCPStartPZ->at(j);
               float pSt = sqrt(pxSt*pxSt+pySt*pySt+pzSt*pzSt);
               int ID = MCTrkID->at(j);
               TVector3 MCmu(MCPStartPX->at(j),MCPStartPY->at(j),MCPStartPZ->at(j));
               float ThetaNuMu = MCNu.Angle(MCmu);

               if (PDGMother->at(j)==0 && numucc==0)
               {
                   hmu_all_thetanu->Fill(ThetaNuMu * (180.0/3.141592653589793238463));
                   hmu_all_pmu->Fill(pSt);
                   numucc++;
               }
    

               int InGAr=0;

               //Find first point in GAr
               for (Int_t k=0; k<TrajMCPX->size(); k++) 
                {
                    if (TrajMCPTrackID->at(k)==ID)
                    {
                        float r = sqrt((TrajMCPY->at(k)-GArCenter.Y())*(TrajMCPY->at(k)-GArCenter.Y())+(TrajMCPZ->at(k)-GArCenter.Z())*(TrajMCPZ->at(k)-GArCenter.Z()));

                        if(TrajMCPX->at(k)<0.5*GAr_L && TrajMCPX->at(k)>-0.5*GAr_L && r<GAr_r && InGAr==0)
                        {
                          px=TrajMCPPX->at(k);
                          py=TrajMCPPY->at(k);
                          pz=TrajMCPPZ->at(k);
                          p = sqrt(px*px+py*py+pz*pz);
                          x=TrajMCPX->at(k);
                          y=TrajMCPY->at(k);
                          z=TrajMCPZ->at(k);
                          InGAr++;
                        }
                    }
                }

               //TRACK MAtching
               
               int nTracks = TrackStartX->size();
               TVector3 MCpart(px,py,pz);
               TVector3 MCpartHat = MCpart.Unit();
               double cosTmax = -1000;
               int iRECOpart = -1;
               TVector3 RECO_P;
               int whichEnd = -1;

               if(InGAr!=0)
               {
                    //hmu_all_thetanu->Fill(ThetaNu_deg);
                    if (PDGMother->at(j)==0 && numuccInGAr==0)
                    {
                        hmu_all_thetanu_InGAr->Fill(ThetaNuMu * (180.0/3.141592653589793238463));
                        hmu_all_pmu_InGAr->Fill(pSt);
                        numuccInGAr++;
                    }
                    
                    for (int iTrack=0; iTrack<nTracks; ++iTrack) 
                    {
                            TVector3 RECO_P_beg(TrackStartPX->at(iTrack), TrackStartPY->at(iTrack),TrackStartPZ->at(iTrack));
                            TVector3 RECO_P_end(TrackEndPX->at(iTrack), TrackEndPY->at(iTrack),  TrackEndPZ->at(iTrack));

                            // Direction matching
                            double cosTbeg = MCpartHat.Dot(RECO_P_beg)/RECO_P_beg.Mag();
                            double cosTend = MCpartHat.Dot(RECO_P_end)/RECO_P_end.Mag();

                            hnhits->Fill(NTPCClustersOnTrack->at(iTrack));

                            if (cosTbeg > cosTend) 
                            {
                                if (cosTbeg > cosTmax) 
                                {
                                    cosTmax   = cosTbeg;
                                    iRECOpart = iTrack;
                                    whichEnd  = 0;
                                    RECO_P	  = RECO_P_beg;
                                }
                            } 
                            else 
                            {
                                if (cosTend > cosTmax) 
                                {
                                    cosTmax   = cosTend;
                                    iRECOpart = iTrack;
                                    whichEnd  = 1;
                                    RECO_P	  = RECO_P_beg;
                                }
                            }   
                    } // End loop over all Tracks

                        //std::cout<<"cosTmax"<<cosTmax<<std::endl;
                        //std::cout<<"iRecoPart"<<iRECOpart<<std::endl;

                    if (iRECOpart == -1) continue;

                    // Plot offset of reco track

                    TVector3 vecX;
                    float mindist = 1000;
                    
                    for (Int_t k=0; k<TrajMCPX->size(); k++) 
                    {
                        if (TrajMCPTrackID->at(k)==ID)
                        {
                            if ( whichEnd == 0 ) 
                            {
                                float dist = sqrt((TrackStartX->at(iRECOpart) -x)*(TrackStartX->at(iRECOpart) -x)+(TrackStartY->at(iRECOpart) -y)*(TrackStartY->at(iRECOpart) -y)+(TrackStartZ->at(iRECOpart) -z)*(TrackStartZ->at(iRECOpart) -z));
                                if(dist<mindist)
                                {
                                vecX.SetXYZ(TrackStartX->at(iRECOpart) -TrajMCPX->at(k),
                                            TrackStartY->at(iRECOpart) -TrajMCPY->at(k),
                                            TrackStartZ->at(iRECOpart) -TrajMCPY->at(k));
                                MCpart.SetXYZ(TrajMCPX->at(k),TrajMCPY->at(k),TrajMCPZ->at(k));
                                MCpartHat=MCpart.Unit();
                                TVector3 RECO_P_beg(TrackStartPX->at(iRECOpart), TrackStartPY->at(iRECOpart),TrackStartPZ->at(iRECOpart));
                                cosTmax=MCpartHat.Dot(RECO_P_beg)/RECO_P_beg.Mag();
                                }
                            } 
                            else if (whichEnd==1)
                            {
                                float dist = sqrt((TrackEndX->at(iRECOpart) -x)*(TrackEndX->at(iRECOpart) -x)+(TrackEndY->at(iRECOpart) -y)*(TrackEndY->at(iRECOpart) -y)+(TrackEndZ->at(iRECOpart) -z)*(TrackEndZ->at(iRECOpart) -z));
                                if(dist<mindist)
                                {
                                vecX.SetXYZ(TrackEndX->at(iRECOpart) -TrajMCPX->at(k),
                                            TrackEndY->at(iRECOpart) -TrajMCPY->at(k),
                                            TrackEndZ->at(iRECOpart) -TrajMCPY->at(k));
                                MCpart.SetXYZ(TrajMCPX->at(k),TrajMCPY->at(k),TrajMCPZ->at(k));
                                MCpartHat=MCpart.Unit();
                                TVector3 RECO_P_beg(TrackEndPX->at(iRECOpart), TrackEndPY->at(iRECOpart),TrackEndPZ->at(iRECOpart));
                                cosTmax=MCpartHat.Dot(RECO_P_beg)/RECO_P_beg.Mag();
                                }
                            }
                        }
                    }
                    Float_t delX = vecX.Cross(MCpartHat).Mag();
                    hz->Fill(vecX.Z());
                    hcos->Fill(cosTmax);
                    hdelX->Fill(delX);

                    if ( cosTmax <= 0.997 ) continue;
                    if (  delX   >= 3.0   ) continue;

                    hnhitsmatched->Fill(NTPCClustersOnTrack->at(iRECOpart));

                    if(cosTmax > 0.997 && delX   < 3.0 && numucc!=0 && recoInGAr==0 && PDGMother->at(j)==0) 
                    {
                      hmu_good_reco_thetanu->Fill(ThetaNuMu* (180.0/3.141592653589793238463));
                      hmu_good_reco_pmu->Fill(pSt);
                      recoInGAr++;
                    }
                    

                    float preco = sqrt(TrackStartPX->at(iRECOpart)*TrackStartPX->at(iRECOpart)+TrackStartPY->at(iRECOpart)*TrackStartPY->at(iRECOpart)+TrackStartPZ->at(iRECOpart)*TrackStartPZ->at(iRECOpart));
                    
                    if (whichEnd==0)
                    {
                        
                        //if(p<=2.5 && p>=2.3) 
                        if(p<=2.5 && p>=2.3&&NTPCClustersOnTrack->at(iRECOpart)>3&&NTPCClustersOnTrack->at(iRECOpart)<11)frac_resid->Fill((preco-p)/p);
                        //if(p<=2.5 && p>=2.3) 
                        if(p<=2.5 && p>=2.3&&NTPCClustersOnTrack->at(iRECOpart)>3&&NTPCClustersOnTrack->at(iRECOpart)<11)frac_resid_c->Fill((preco-p)/p);
                        ptrue->Fill(p);
                        pVSfrac_resid->Fill(p,(preco-p)/p);
                        PpVSfrac_resid->Fill(p,(preco-p)/p);
                        if((PDG->at(j)==13 && TrackStartQ->at(iRECOpart)==-1) || (PDG->at(j)==-13 && TrackStartQ->at(iRECOpart)==1))
                        {
                            hmu_good_recocharge_p_start->Fill(p);
                            
                        }
                    }
                    else if (whichEnd==1)
                    {
                        //if(p<=2.5 && p>=2.3) 
                        if(p<=2.5 && p>=2.3&&NTPCClustersOnTrack->at(iRECOpart)>3&&NTPCClustersOnTrack->at(iRECOpart)<11)frac_resid->Fill((preco-p)/p);
                        //if(p<=2.5 && p>=2.3) 
                        if(p<=2.5 && p>=2.3&&NTPCClustersOnTrack->at(iRECOpart)>3&&NTPCClustersOnTrack->at(iRECOpart)<11)frac_resid_c->Fill((preco-p)/p);
                        ptrue->Fill(p);
                        pVSfrac_resid->Fill(p,(preco-p)/p);
                        PpVSfrac_resid->Fill(p,(preco-p)/p);
                        if((PDG->at(j)==13 && TrackEndQ->at(iRECOpart)==-1) || (PDG->at(j)==-13 && TrackEndQ->at(iRECOpart)==1))
                        {
                            hmu_good_recocharge_p_start->Fill(p);
                            
                        }   
                    }

                hmu_all_recocharge_p_start->Fill(p);
                }

               
           }
        }

        if(showprog==true) std::cout << string(strprog.length(),'\b')<<"\b";
        

    }


    for (Int_t x=1; x<=50; x++)
    {
        float max=0;
        
    
        for (Int_t y=1; y<=50; y++)
        {
            if (max<PpVSfrac_resid->GetBinContent(x,y))  max=PpVSfrac_resid->GetBinContent(x,y);
        }
        
        for (Int_t y=1; y<=100; y++)
        {
            float bin=0;

            if (max!=0) bin=PpVSfrac_resid->GetBinContent(x,y)/max;

            PpVSfrac_resid->SetBinContent(x,y,bin);
       
        }
        //std::cout<<"sum: "<<sum<<std::endl;
        //std::cout<<"sumcheck: "<<sumcheck<<std::endl;
    }

    TH2F *hmu_good_reco_thetanu_2 = (TH2F*)hmu_good_reco_thetanu->Clone("hmu_good_reco_thetanu_2");
    TH2F *hmu_good_reco_pmu_2 = (TH2F*)hmu_good_reco_pmu->Clone("hmu_good_reco_pmu_2");

    TCanvas *mccanvasnhits = new TCanvas("mccanvasnhits","",1000,800);
    hnhits->SetTitle("n hits per track;nhits;n");
    hnhits->Draw();
    mccanvasnhits->Print("test/nhits.png");

    TCanvas *mccanvasnhitsmatched = new TCanvas("mccanvasnhitsmatched","",1000,800);
    hnhitsmatched->SetTitle("n hits per matched track;nhits;n");
    hnhitsmatched->Draw();
    mccanvasnhitsmatched->Print("test/nhitsmatched.png");

    TCanvas *mccanvasdz = new TCanvas("mccanvasndz","",1000,800);
    hz->SetTitle("dz;dz(cm);n");
    hz->Draw();
    mccanvasdz->Print("test/dz.png");

    TCanvas *mccanvascos = new TCanvas("mccanvascos","",1000,800);
    //gStyle->SetOptStat(0);
    hcos->SetTitle("cos#theta;cos#theta;n");
    //hcos->Divide(hmu_all_thetanu);
    //hcos->SetMinimum(0);
    hcos->Draw();
    mccanvascos->Print("test/cos.png");

    TCanvas *mccanvasdelX = new TCanvas("mccanvasdelX","",1000,800);
    //gStyle->SetOptStat(0);
    hdelX->SetTitle("delX;delX(cm);n");
    //hcos->Divide(hmu_all_thetanu);
    //hcos->SetMinimum(0);
    hdelX->Draw();
    mccanvasdelX->Print("test/delX.png");

    TCanvas *mccanvasReco = new TCanvas("mccanvasReco","",1000,800);
    gStyle->SetOptStat(0);
    hmu_good_reco_thetanu->SetTitle("Reconstruction efficiency as function of angle;#theta_{#nu#mu} [deg];n(reco)/n(primary #mu in LAr)");
    hmu_good_reco_thetanu->Divide(hmu_all_thetanu);
    hmu_good_reco_thetanu->SetMinimum(0);
    hmu_good_reco_thetanu->Draw();
    mccanvasReco->Print("test/AngleRecoEfficiency.png");

    TCanvas *mccanvasReco2 = new TCanvas("mccanvasReco2","",1000,800);
    gStyle->SetOptStat(0);
    hmu_good_reco_thetanu_2->SetTitle("Reconstruction efficiency as function of angle;#theta_{#nu#mu} [deg];n(reco)/n(primary #mu in LAr+GAr)");
    hmu_good_reco_thetanu_2->Divide(hmu_all_thetanu_InGAr);
    hmu_good_reco_thetanu_2->SetMinimum(0);
    hmu_good_reco_thetanu_2->Draw();
    mccanvasReco2->Print("test/AngleRecoEfficiency_InGAr.png");

    TCanvas *mccanvasRecop = new TCanvas("mccanvasRecop","",1000,800);
    gStyle->SetOptStat(0);
    hmu_good_reco_pmu->SetTitle("Reconstruction efficiency as function of momentum;p_{#mu}^{start} [GeV/c];n(reco)/n(primary #mu in LAr)");
    hmu_good_reco_pmu->Divide(hmu_all_pmu);
    hmu_good_reco_pmu->SetMinimum(0);
    hmu_good_reco_pmu->Draw();
    mccanvasRecop->Print("test/PRecoEfficiency.png");

    TCanvas *mccanvasRecop2 = new TCanvas("mccanvasRecop2","",1000,800);
    gStyle->SetOptStat(0);
    hmu_good_reco_pmu_2->SetTitle("Reconstruction efficiency as function of momentum;p_{#mu}^{start} [GeV/c];n(reco)/n(primary #mu in LAr+GAr)");
    hmu_good_reco_pmu_2->Divide(hmu_all_pmu_InGAr);
    hmu_good_reco_pmu_2->SetMinimum(0);
    hmu_good_reco_pmu_2->Draw();
    mccanvasRecop2->Print("test/PRecoEfficiency_InGAr.png");
    
    TCanvas *mccanvascharge = new TCanvas("mccanvascharge","",1000,800);
    gStyle->SetOptStat(0);
    hmu_good_recocharge_p_start->SetTitle("Charge reconstruction resolution;p_{true}[GeV/c];n(reco ok)/n(all muons in GAr)");
    hmu_good_recocharge_p_start->Divide(hmu_all_recocharge_p_start);
    hmu_good_recocharge_p_start->SetMinimum(0.5);
    hmu_good_recocharge_p_start->Draw();
    mccanvascharge->Print("test/ChargeRes.png");

    TCanvas *mccanvas_fracresid = new TCanvas("mccanvas_fracresid","",1000,800);
    gStyle->SetOptStat(1);
    frac_resid->SetTitle("Momentum fractional residuals GeV/c (Gauss Fit);(p_{reco}-p_{true})/p_{true};n(muons in GAr)");
    frac_resid->Fit("gaus");
    gStyle->SetOptFit(1);
    frac_resid->Draw();
    mccanvas_fracresid->Print("test/Frac_resid_B3M10.png");

    TCanvas *mccanvas_fracresid_cauchy = new TCanvas("mccanvas_fracresid_cauchy","",1000,800);
    gStyle->SetOptStat(1);
    TF1 *fitfunc = new TF1("fitfunc",CauchyPeak,-10,10,3);
    Double_t fitpar[3];
    fitpar[0] = 0.05;
    fitpar[1] = 0.05;
    fitpar[2] = frac_resid_c->GetEntries();
    fitfunc->SetParameters(fitpar);
    frac_resid_c->Draw("same");
    frac_resid_c->SetTitle("Momentum fractional residuals GeV/c (Cauchy Fit);(p_{reco}-p_{true})/p_{true};n(muons in GAr)");
    frac_resid_c->Fit("fitfunc","","",-0.09,0.09);
    gStyle->SetOptFit(1);
    mccanvas_fracresid_cauchy->Print("test/Frac_resid_cauchy_B3M10.png");

    TCanvas *mccanvas_truep = new TCanvas("mccanvas_truep","",1000,800);
    gStyle->SetOptStat(1);
    ptrue->SetTitle("True momentum of muons in GAr;p_{true}(GeV/c);n(muons in GAr)");
    ptrue->Draw();
    mccanvas_truep->Print("test/Truep_InGAr.png");

    TCanvas *mccanvas_pVSfracresid = new TCanvas("mccanvas_pVSfracresid","",1000,800);
    gStyle->SetOptStat(0);
    pVSfrac_resid->SetTitle("Momentum fractional residuals VS p_{true};p_{true}(GeV/c);(p_{reco}-p_{true})/p_{true}");
    pVSfrac_resid->Draw("COLZ");
    mccanvas_pVSfracresid->Print("test/pVSFrac_resid.png");

    TCanvas *mccanvas_pVSfracresidP = new TCanvas("mccanvas_pVSfracresidP","",1000,800);
    gStyle->SetOptStat(0);
    PpVSfrac_resid->SetTitle("P(Residual|p_{true});p_{true}(GeV/c);(p_{reco}-p_{true})/p_{true}");
    PpVSfrac_resid->Draw("COLZ");
    mccanvas_pVSfracresidP->Print("test/PpVSFrac_resid.png");

    
}