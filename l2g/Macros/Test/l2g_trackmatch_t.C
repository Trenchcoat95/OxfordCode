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
#include "garana.h"

//   In a ROOT session, you can do:
  //      root> .L garana.C
  //      root> garana t
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.Loop();       // Loop on all entries

void GetProbabilityPlot(TH2F * HSTD, TH2F * &HPROB)
{
    Int_t nbinsx=HSTD->GetNbinsX();
    Int_t nbinsy=HSTD->GetNbinsY();

    for (Int_t x=1; x<=nbinsx; x++)
    {
        float max=0;
        for (Int_t y=1; y<=nbinsy; y++)
        {
            if (max<HSTD->GetBinContent(x,y))  max=HSTD->GetBinContent(x,y);
        }
        
        for (Int_t y=1; y<=nbinsy; y++)
        {
            float bin=0;
            if (max!=0) bin=HSTD->GetBinContent(x,y)/max;
            HPROB->SetBinContent(x,y,bin);       
        }
    }
}


void garana::l2g_Trackmatch()
{

    //TVector3 GArCenter(0,-150.473,1486); //Correct values
    //bool edge=false;
    TVector3 GArCenter(0,-68.287,1486);  //Edge sample values
    bool edge=true;
    float GAr_r = 349.9;
    float GAr_L = 669.6;
    Int_t nentries = fChain->GetEntries();


    #pragma region "Plot Declaration"

        /////////////////////////////////////1D plots of the three main variables
        TH1F* hnhits_RecoTracks_Sample = new TH1F("hnhits_RecoTracks_Sample", "nhits distribution for all muon tracks in GAr", 20, 0, 20);
        TH1F* hptrueStart_InGAr_Sample = new TH1F("hptrueStart_InGAr_Sample", "ptrue distribution for the sample of muons traversing GAr", 50, 0, 8);
        TH1F* htheta_InGAr_Sample = new TH1F("hmu_theta_InGAr_Sample", "Thetatrue distribution for the sample of muons traversing GAr", 180, 0, 180.0);

        /////////////////////////////////////1D Distribution of matched tracks
        TH1F* hnhits_RecoTracks_Matched = new TH1F("hnhits_RecoTracks_Matched", "nhits distribution for all matched muon tracks in GAr", 20, 0, 20);
        TH1F* hptrueStart_InGAr_Matched = new TH1F("hptrueStart_InGAr_Matched", "ptrue distribution for all matched muon tracks", 50, 0, 8);
        TH1F* htheta_InGAr_Matched = new TH1F("htheta_InGAr_Matched", "Thetatrue distribution for all matched muon tracks", 180, 0, 180.0);

        /////////////////////////////////////Charge resolution
        TH1F* hmu_good_recocharge_p_start = new TH1F("hmu_good_recocharge_p_start", "Reconstructed charge resolution", 50, 0, 8.0);
        TH1F* hmu_all_recocharge_p_start = new TH1F("hmu_all_recocharge_p_start", "Reconstructed charge resolution", 50, 0, 8.0);
        
        /////////////////////////////////////1D Fractional Residual distribution for the reconstructed muons
        TH1F* frac_resid = new TH1F("frac_resid", "Fractional residuals", 150, -0.4, 0.4);

        /////////////////////////////////////Fractional Residuals 2D   
        TH2F* pVSfrac_resid = new TH2F("pVSfrac_resid", "Fractional residuals",50,0,8, 50, -1, 1);
        TH2F* PpVSfrac_resid = new TH2F("PpVSfrac_resid", "P(Residual|p)",50,0,8, 50, -1, 1);
        TH2F* nhitsVSfrac_resid = new TH2F("nhitsVSfrac_resid", "Fractional residuals",20,0,20, 50, -1, 1);
        TH2F* PnhitsVSfrac_resid = new TH2F("PnhitsVSfrac_resid", "P(Residual|Theta)",20,0,20, 50, -1, 1);
        TH2F* ThetaVSfrac_resid = new TH2F("ThetaVSfrac_resid", "Fractional residuals",80,0,80, 50, -1, 1);
        TH2F* PThetaVSfrac_resid = new TH2F("PThetaVSfrac_resid", "P(Residual|Theta)",80,0,80, 50, -1, 1);

    #pragma endregion    

    bool showprog = true;  
    if(showprog==true) std::cout<<"Progress:  "<<std::endl;

    for (Int_t i=0; i<nentries; i++) 
    {
        fChain->GetEntry(i);
        int prog = 100*i/nentries;
        std::string strprog = std::to_string(prog);
        if(showprog==true) std::cout<<strprog<<"%";
        
        float pxNu= MCNuPx->at(0);
        float pyNu= MCNuPy->at(0);
        float pzNu= MCNuPz->at(0);
        TVector3 MCNu(pxNu,pyNu,pzNu);

        //Cycle over MC Tracks
        for (Int_t j=0; j<PDG->size(); j++) 
        {
           if((PDG->at(j)==13 || PDG->at(j)==-13) && PDGMother->at(j)==0) //consider only primary muons
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

               int InGAr=0;

               //Cycle over all trajectory points, Find the ones corresponding to the current MC Track and Find first point in GAr
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
                          break;
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
               int RECO_CHARGE;
               int whichEnd = -1;

               if(InGAr!=0)
               {
                    //Fill the MC track in GAr distributions  
                    htheta_InGAr_Sample->Fill(ThetaNuMu * (180.0/3.141592653589793238463));
                    hptrueStart_InGAr_Sample->Fill(pSt);
                    
                    //looping over all reconstructed tracks
                    for (int iTrack=0; iTrack<nTracks; ++iTrack) 
                    {
                            hnhits_RecoTracks_Sample->Fill(NTPCClustersOnTrack->at(iTrack));

                            TVector3 RECO_P_beg(TrackStartPX->at(iTrack), TrackStartPY->at(iTrack),TrackStartPZ->at(iTrack));
                            TVector3 RECO_P_end(TrackEndPX->at(iTrack), TrackEndPY->at(iTrack),  TrackEndPZ->at(iTrack));
                            int RECO_CHARGE_beg=TrackStartQ->at(iTrack);
                            int RECO_CHARGE_end=TrackEndQ->at(iTrack);

                            // Direction matching
                            double cosTbeg = MCpartHat.Dot(RECO_P_beg)/RECO_P_beg.Mag();
                            double cosTend = MCpartHat.Dot(RECO_P_end)/RECO_P_end.Mag();
            
                            if (cosTbeg > cosTend) 
                            {
                                if (cosTbeg > cosTmax) 
                                {
                                    cosTmax   = cosTbeg;
                                    iRECOpart = iTrack;
                                    whichEnd  = 0;
                                    RECO_P	  = RECO_P_beg;
                                    RECO_CHARGE = RECO_CHARGE_beg;
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
                                    RECO_CHARGE = RECO_CHARGE_end;
                                }
                            }   
                    } // End loop over all Tracks

                    if (iRECOpart == -1) continue;

                    TVector3 vecX;
                    if ( whichEnd == 0 ) 
                    {
                        vecX.SetXYZ(TrackStartX->at(iRECOpart) -x,
                                    TrackStartY->at(iRECOpart) -y,
                                    TrackStartZ->at(iRECOpart) -z);
                    } 
                    else if (whichEnd==1)
                    {
                        vecX.SetXYZ(TrackEndX->at(iRECOpart)   -x,
                                    TrackEndY->at(iRECOpart)   -y,
                                    TrackEndZ->at(iRECOpart)   -z);
                    }
                    Float_t delX = vecX.Cross(MCpartHat).Mag();
                    if ( cosTmax <= 0.997 ) continue;
                    if (  delX   >= 3.0   ) continue;

                    //If we are here the track was matched

                    
                    //Fill distributions for matched tracks
                    htheta_InGAr_Matched->Fill(ThetaNuMu* (180.0/3.141592653589793238463));
                    hptrueStart_InGAr_Matched->Fill(pSt);
                    hnhits_RecoTracks_Matched->Fill(NTPCClustersOnTrack->at(iRECOpart));
                         
                    //Interesting Region if(p<=2.5 && p>=2.3)
                    //Fill fractional residual plots
                    double preco = RECO_P.Mag();
                    frac_resid->Fill((preco-p)/p);
                    pVSfrac_resid->Fill(p,(preco-p)/p);
                    nhitsVSfrac_resid->Fill(NTPCClustersOnTrack->at(iRECOpart),(preco-p)/p);
                    ThetaVSfrac_resid->Fill(ThetaNuMu* (180.0/3.141592653589793238463),(preco-p)/p);

                    //Fill charge reconstruction resolution
                    if((PDG->at(j)==13 && RECO_CHARGE==-1) || (PDG->at(j)==-13 && RECO_CHARGE==1)) hmu_good_recocharge_p_start->Fill(p);
                    hmu_all_recocharge_p_start->Fill(p);
                }

               
           }
        }

        if(showprog==true) std::cout << std::string(strprog.length(),'\b')<<"\b";
        

    }


    GetProbabilityPlot(pVSfrac_resid,PpVSfrac_resid);
    GetProbabilityPlot(nhitsVSfrac_resid,PnhitsVSfrac_resid);
    GetProbabilityPlot(ThetaVSfrac_resid,PThetaVSfrac_resid);
    
    #pragma region "Plotting"

        ///////////////////////////////////////////////1D distributions for three main variables
        TCanvas *mccanvasnhits = new TCanvas("mccanvasnhits","",1000,800);
        hnhits_RecoTracks_Sample->SetTitle("number of hits distribution for all muon tracks in GAr;n_{hits};n");
        hnhits_RecoTracks_Sample->Draw();
        if (edge) mccanvasnhits->Print("Edge/nhits1D_edge.png");
        else mccanvasnhits->Print("Standard/nhits1D.png");

        TCanvas *mccanvasp = new TCanvas("mccanvasp","",1000,800);
        hptrueStart_InGAr_Sample->SetTitle("p_{true} distribution for the sample of muons traversing GAr;p_{true}(GeV/c);n");
        hptrueStart_InGAr_Sample->Draw();
        if (edge) mccanvasp->Print("Edge/ptrue1D_edge.png");
        else mccanvasp->Print("Standard/ptrue1D.png");

        TCanvas *mccanvasTheta = new TCanvas("mccanvasTheta","",1000,800);
        htheta_InGAr_Sample->SetTitle("#theta_{true} distribution for the sample of muons traversing GAr;#theta_{true}(deg);n");
        htheta_InGAr_Sample->Draw();
        if (edge) mccanvasTheta->Print("Edge/Thetatrue1D_edge.png");
        else mccanvasTheta->Print("Standard/Thetatrue1D.png");
    
        ///////////////////////////////////////////////1D distributions of matched tracksfor three main variables
        TCanvas *mccanvasnhitsm = new TCanvas("mccanvasnhitsm","",1000,800);
        hnhits_RecoTracks_Matched->SetTitle("nhits distribution of matched muon tracks in GAr;nhits;n");
        hnhits_RecoTracks_Matched->Draw();
        if (edge) mccanvasnhitsm->Print("Edge/nhits1D_matched_edge.png");
        else mccanvasnhitsm->Print("Standard/nhits1D_matched.png");

        TCanvas *mccanvaspm = new TCanvas("mccanvaspm","",1000,800);
        hptrueStart_InGAr_Matched->SetTitle("p_{true} distribution for the sample of matched muon tracks;p_{true}(GeV/c);n");
        hptrueStart_InGAr_Matched->Draw();
        if(edge) mccanvaspm->Print("Edge/ptrue1D_matched_edge.png");
        else mccanvaspm->Print("Standard/ptrue1D_matched.png");

        TCanvas *mccanvasThetam = new TCanvas("mccanvasThetam","",1000,800);
        htheta_InGAr_Matched->SetTitle("#theta_{true} distribution for the sample of matched muon tracks;#theta_{true}(deg);n");
        htheta_InGAr_Matched->Draw();
        if(edge) mccanvasThetam->Print("Edge/Thetatrue1D_matched_edge.png");
        else mccanvasThetam->Print("Standard/Thetatrue1D_matched.png");

        ///////////////////////////////////////////////Efficiency plots
        TCanvas *mccanvasHitsEff = new TCanvas("mccanvasHitsEff","",1000,800);
        gStyle->SetOptStat(0);
        hnhits_RecoTracks_Matched->SetTitle("Reconstruction (track matching) efficiency as function of number of hits;p_{true} [GeV/c];n(matched)/n(Reco Tracks)");
        hnhits_RecoTracks_Matched->Divide(hnhits_RecoTracks_Sample);
        hnhits_RecoTracks_Matched->SetMinimum(0);
        hnhits_RecoTracks_Matched->Draw();
        if (edge) mccanvasHitsEff->Print("Edge/nhits1D_Efficiency_edge.png");
        else mccanvasHitsEff->Print("Standard/nhits1D_Efficiency.png");

        TCanvas *mccanvasPEff = new TCanvas("mccanvasPEff","",1000,800);
        gStyle->SetOptStat(0);
        hptrueStart_InGAr_Matched->SetTitle("Reconstruction (track matching) efficiency as function of initial momentum;p_{true} [GeV/c];n(matched)/n(MC Tracks in GAr)");
        hptrueStart_InGAr_Matched->Divide(hptrueStart_InGAr_Sample);
        hptrueStart_InGAr_Matched->SetMinimum(0);
        hptrueStart_InGAr_Matched->Draw();
        if(edge) mccanvasPEff->Print("Edge/ptrue1D_Efficiency_edge.png");
        else mccanvasPEff->Print("Standard/ptrue1D_Efficiency.png");

        TCanvas *mccanvasAngleEff = new TCanvas("mccanvasAngleEff","",1000,800);
        gStyle->SetOptStat(0);
        htheta_InGAr_Matched->SetTitle("Reconstruction (track matching) efficiency as function of initial momentum;#theta_{#nu#mu} [deg];n(matched)/n(MC Tracks in GAr)");
        htheta_InGAr_Matched->Divide(htheta_InGAr_Sample);
        htheta_InGAr_Matched->SetMinimum(0);
        htheta_InGAr_Matched->Draw();
        if(edge) mccanvasAngleEff->Print("Edge/Thetatrue1D_Efficiency_edge.png");
        else mccanvasAngleEff->Print("Standard/Thetatrue1D_Efficiency.png");

        /////////////////////////////////////////////////////Charge Resolution
        TCanvas *mccanvascharge = new TCanvas("mccanvascharge","",1000,800);
        gStyle->SetOptStat(0);
        hmu_good_recocharge_p_start->SetTitle("Charge reconstruction resolution;p_{true}[GeV/c];n(correct)/n(matched)");
        hmu_good_recocharge_p_start->Sumw2();
        hmu_good_recocharge_p_start->Divide(hmu_good_recocharge_p_start,hmu_all_recocharge_p_start,1,1,"B");
        hmu_good_recocharge_p_start->SetMinimum(0.5);
        hmu_good_recocharge_p_start->Draw();
        if (edge) mccanvascharge->Print("Edge/ChargeRes_ErrorBars_Edge.png");
        else mccanvascharge->Print("Standard/ChargeRes_ErrorBars_Edge.png");

        ////////////////////////////////////////////////////Fractional Residual 1D plot
        TCanvas *mccanvas_fracresid = new TCanvas("mccanvas_fracresid","",1000,800);
        gStyle->SetOptStat(1);
        frac_resid->SetTitle("Momentum fractional residuals (Double Gauss Fit);(p_{reco}-p_{true})/p_{true};n(#mu in GAr-Lite)");
        TF1 *double_gauss = new TF1("double_gauss","[0]*(exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-([1]+abs([4])))/[5])^2))",-0.4,0.4);
        double_gauss->SetParameter(0,8000);
        double_gauss->SetParameter(1,-0.01);
        double_gauss->SetParameter(2,0.05);
        double_gauss->SetParameter(3,0.01);
        double_gauss->SetParameter(4,0);
        double_gauss->SetParameter(5,0.1);
        frac_resid->Fit("double_gauss");
        gStyle->SetOptFit(1);
        frac_resid->Draw();
        if (edge) mccanvas_fracresid->Print("Edge/Frac_resid_two_gauss_fullspectrum_edge.png");
        else mccanvas_fracresid->Print("Standard/Frac_resid_two_gauss_fullspectrum.png");

        ///////////////////////////////////////////////////////Fractional ResidualVSptrue
        TCanvas *mccanvas_pVSfracresid = new TCanvas("mccanvas_pVSfracresid","",1000,800);
        gStyle->SetOptStat(0);
        pVSfrac_resid->SetTitle("Momentum fractional residuals VS p_{true};p_{true}(GeV/c);(p_{reco}-p_{true})/p_{true}");
        pVSfrac_resid->Draw("COLZ");
        if (edge) mccanvas_pVSfracresid->Print("Edge/pVSFrac_resid_edge.png");
        else  mccanvas_pVSfracresid->Print("Standard/pVSFrac_resid.png");

        TCanvas *mccanvas_pVSfracresidP = new TCanvas("mccanvas_pVSfracresidP","",1000,800);
        gStyle->SetOptStat(0);
        PpVSfrac_resid->SetTitle("P(Residual|p_{true});p_{true}(GeV/c);(p_{reco}-p_{true})/p_{true}");
        PpVSfrac_resid->Draw("COLZ");
        if (edge) mccanvas_pVSfracresidP->Print("Edge/PpVSFrac_resid_edge.png");
        else  mccanvas_pVSfracresidP->Print("Standard/PpVSFrac_resid.png");

        ///////////////////////////////////////////////////////Fractional ResidualVSnhits
        TCanvas *mccanvas_hitsVSfracresid = new TCanvas("mccanvas_hitsVSfracresid","",1000,800);
        gStyle->SetOptStat(0);
        nhitsVSfrac_resid->SetTitle("Momentum fractional residuals VS n_{hits};n_{hits};(p_{reco}-p_{true})/p_{true}");
        nhitsVSfrac_resid->Draw("COLZ");
        if (edge) mccanvas_hitsVSfracresid->Print("Edge/nhitsVSFrac_resid_edge.png");
        else  mccanvas_hitsVSfracresid->Print("Standard/nhitsVSFrac_resid.png");

        TCanvas *mccanvas_hitsVSfracresidP = new TCanvas("mccanvas_hitsVSfracresidP","",1000,800);
        gStyle->SetOptStat(0);
        PnhitsVSfrac_resid->SetTitle("P(Residual|n_{hits});n_{hits};(p_{reco}-p_{true})/p_{true}");
        PnhitsVSfrac_resid->Draw("COLZ");
        if (edge) mccanvas_hitsVSfracresidP->Print("Edge/PnhitsVSFrac_resid_edge.png");
        else  mccanvas_hitsVSfracresidP->Print("Standard/PnhitsVSFrac_resid.png");

        ///////////////////////////////////////////////////////Fractional ResidualVSTheta
        TCanvas *mccanvas_ThetaVSfracresid = new TCanvas("mccanvas_ThetaVSfracresid","",1000,800);
        gStyle->SetOptStat(0);
        ThetaVSfrac_resid->SetTitle("Momentum fractional residuals VS #theta;#theta;(p_{reco}-p_{true})/p_{true}");
        ThetaVSfrac_resid->Draw("COLZ");
        if (edge) mccanvas_ThetaVSfracresid->Print("Edge/ThetaVSFrac_resid_edge.png");
        else  mccanvas_ThetaVSfracresid->Print("Standard/ThetaVSFrac_resid.png");

        TCanvas *mccanvas_ThetaVSfracresidP = new TCanvas("mccanvas_ThetaVSfracresidP","",1000,800);
        gStyle->SetOptStat(0);
        PThetaVSfrac_resid->SetTitle("P(Residual|#theta);#theta;(p_{reco}-p_{true})/p_{true}");
        PThetaVSfrac_resid->Draw("COLZ");
        if (edge) mccanvas_ThetaVSfracresidP->Print("Edge/PThetaVSFrac_resid_edge.png");
        else  mccanvas_ThetaVSfracresidP->Print("Standard/PThetaVSFrac_resid.png");

    #pragma endregion
    /*
    TH2F *hmu_good_reco_thetanu_2 = (TH2F*)hmu_good_reco_thetanu->Clone("hmu_good_reco_thetanu_2");
    TH2F *hmu_good_reco_pmu_2 = (TH2F*)hmu_good_reco_pmu->Clone("hmu_good_reco_pmu_2");

    TCanvas *mccanvasnhits = new TCanvas("mccanvasnhits","",1000,800);
    hnhits->SetTitle("n hits per track;nhits;n");
    hnhits->Draw();
    //mccanvasnhits->Print("nhits.png");

    TCanvas *mccanvasdz = new TCanvas("mccanvasndz","",1000,800);
    hz->SetTitle("dz;dz(cm);n");
    hz->Draw();
    //mccanvasdz->Print("dz.png");

    TCanvas *mccanvascos = new TCanvas("mccanvascos","",1000,800);
    //gStyle->SetOptStat(0);
    hcos->SetTitle("cos#theta;cos#theta;n");
    //hcos->Divide(hmu_all_thetanu);
    //hcos->SetMinimum(0);
    hcos->Draw();
    //mccanvascos->Print("cos.png");

    TCanvas *mccanvasdelX = new TCanvas("mccanvasdelX","",1000,800);
    //gStyle->SetOptStat(0);
    hdelX->SetTitle("delX;delX(cm);n");
    //hcos->Divide(hmu_all_thetanu);
    //hcos->SetMinimum(0);
    hdelX->Draw();
    //mccanvasdelX->Print("delX.png");

    TCanvas *mccanvasReco = new TCanvas("mccanvasReco","",1000,800);
    gStyle->SetOptStat(0);
    hmu_good_reco_thetanu->SetTitle("Reconstruction efficiency as function of angle;#theta_{#nu#mu} [deg];n(reco)/n(primary #mu in LAr)");
    hmu_good_reco_thetanu->Divide(hmu_all_thetanu);
    hmu_good_reco_thetanu->SetMinimum(0);
    hmu_good_reco_thetanu->Draw();
    //mccanvasReco->Print("Plots/New/AngleRecoEfficiency.png");

    TCanvas *mccanvasReco2 = new TCanvas("mccanvasReco2","",1000,800);
    gStyle->SetOptStat(0);
    hmu_good_reco_thetanu_2->SetTitle("Reconstruction efficiency as function of angle;#theta_{#nu#mu} [deg];n(reco)/n(primary #mu in LAr+GAr)");
    hmu_good_reco_thetanu_2->Divide(hmu_all_thetanu_InGAr);
    hmu_good_reco_thetanu_2->SetMinimum(0);
    hmu_good_reco_thetanu_2->Draw();
    //mccanvasReco2->Print("Plots/New/AngleRecoEfficiency_InGAr.png");

    TCanvas *mccanvasRecop = new TCanvas("mccanvasRecop","",1000,800);
    gStyle->SetOptStat(0);
    hmu_good_reco_pmu->SetTitle("Reconstruction efficiency as function of momentum;p_{#mu}^{start} [GeV/c];n(reco)/n(primary #mu in LAr)");
    hmu_good_reco_pmu->Divide(hmu_all_pmu);
    hmu_good_reco_pmu->SetMinimum(0);
    hmu_good_reco_pmu->Draw();
    //mccanvasRecop->Print("Plots/New/PRecoEfficiency.png");

    TCanvas *mccanvasRecop2 = new TCanvas("mccanvasRecop2","",1000,800);
    gStyle->SetOptStat(0);
    hmu_good_reco_pmu_2->SetTitle("Reconstruction efficiency as function of momentum;p_{#mu}^{start} [GeV/c];n(reco)/n(primary #mu in LAr+GAr)");
    hmu_good_reco_pmu_2->Divide(hmu_all_pmu_InGAr);
    hmu_good_reco_pmu_2->SetMinimum(0);
    hmu_good_reco_pmu_2->Draw();
    //mccanvasRecop2->Print("Plots/New/PRecoEfficiency_InGAr.png");
    
    TCanvas *mccanvascharge = new TCanvas("mccanvascharge","",1000,800);
    gStyle->SetOptStat(0);
    hmu_good_recocharge_p_start->SetTitle("Charge reconstruction resolution;p_{true}[GeV/c];n(reco ok)/n(muons matched in ND-GAr-Lite)");
    hmu_good_recocharge_p_start->Sumw2();
    hmu_good_recocharge_p_start->Divide(hmu_good_recocharge_p_start,hmu_all_recocharge_p_start,1,1,"B");
    hmu_good_recocharge_p_start->SetMinimum(0.5);
    hmu_good_recocharge_p_start->Draw();
    //mccanvascharge->Print("Plots/New/ChargeRes_ErrorBars.png");

    TCanvas *mccanvas_fracresid = new TCanvas("mccanvas_fracresid","",1000,800);
    gStyle->SetOptStat(1);
    frac_resid->SetTitle("Momentum fractional residuals with 2.3<p<2.5 GeV (Gauss Fit);(p_{reco}-p_{true})/p_{true};n(#mu in GAr-Lite)");
    TF1 *double_gauss = new TF1("double_gauss","[0]*(exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2))",-0.4,0.4);
    //double_gauss->SetParameter(0,8000);
    double_gauss->SetParameter(0,8000);
    double_gauss->SetParameter(1,-0.01);
    double_gauss->SetParameter(2,0.05);
    double_gauss->SetParameter(3,0.01);
    double_gauss->SetParameter(4,0);
    double_gauss->SetParameter(5,0.1);
    frac_resid->Fit("double_gauss");
    //frac_resid->Fit("gaus","","",-0.1,0.1);
    gStyle->SetOptFit(1);
    frac_resid->Draw();
    //mccanvas_fracresid->Print("Plots/New/Frac_resid_two_gauss_fullspectrum.png");

    TCanvas *mccanvas_fracresidone = new TCanvas("mccanvas_fracresidone","",1000,800);
    gStyle->SetOptStat(1);
    frac_resid_g->SetTitle("Momentum fractional residuals with 2.3<p<2.5 GeV (Gauss Fit);(p_{reco}-p_{true})/p_{true};n(#mu in GAr-Lite)");
    //TF1 *double_gauss = new TF1("double_gauss","[0]*(exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2))",-0.4,0.4);
    //double_gauss->SetParameter(0,11000);
    //double_gauss->SetParameter(1,-0.01);
    //double_gauss->SetParameter(2,0.02);
    //double_gauss->SetParameter(3,0.05);
    //double_gauss->SetParameter(4,0);
    //double_gauss->SetParameter(5,0.08);
    //frac_resid->Fit("double_gauss");
    frac_resid_g->Fit("gaus","","",-0.1,0.1);
    gStyle->SetOptFit(1);
    frac_resid_g->Draw();
    //mccanvas_fracresidone->Print("Plots/New/Frac_resid_one_gauss_fullspectrum.png");
    //_fullspectrum

    TCanvas *mccanvas_fracresid_cauchy = new TCanvas("mccanvas_fracresid_cauchy","",1000,800);
    gStyle->SetOptStat(1);
    TF1 *fitfunc = new TF1("fitfunc",CauchyPeak,-10,10,3);
    Double_t fitpar[3];
    fitpar[0] = 0.05;
    fitpar[1] = 0.05;
    fitpar[2] = 100;
    fitfunc->SetParameters(fitpar);
    frac_resid_c->Draw("same");
    frac_resid_c->SetTitle("Momentum fractional residuals with 2.3<p<2.5 GeV (Cauchy Fit);(p_{reco}-p_{true})/p_{true};n(#mu in GAr-Lite)");
    frac_resid_c->Fit("fitfunc","","",-0.09,0.09);
    gStyle->SetOptFit(1);
    //mccanvas_fracresid_cauchy->Print("Plots/New/Frac_resid_cauchy_fullspectrum.png");

    TCanvas *mccanvas_truep = new TCanvas("mccanvas_truep","",1000,800);
    gStyle->SetOptStat(1);
    ptrue->SetTitle("True momentum of muons in GAr;p_{true}(GeV/c);n(muons in GAr)");
    ptrue->Draw();
    //mccanvas_truep->Print("Plots/New/Truep_InGAr.png");

    TCanvas *mccanvas_pVSfracresid = new TCanvas("mccanvas_pVSfracresid","",1000,800);
    gStyle->SetOptStat(0);
    pVSfrac_resid->SetTitle("Momentum fractional residuals VS p_{true};p_{true}(GeV/c);(p_{reco}-p_{true})/p_{true}");
    pVSfrac_resid->Draw("COLZ");
    //mccanvas_pVSfracresid->Print("Plots/New/pVSFrac_resid.png");

    TCanvas *mccanvas_pVSfracresidP = new TCanvas("mccanvas_pVSfracresidP","",1000,800);
    gStyle->SetOptStat(0);
    PpVSfrac_resid->SetTitle("P(Residual|p_{true});p_{true}(GeV/c);(p_{reco}-p_{true})/p_{true}");
    PpVSfrac_resid->Draw("COLZ");
    //mccanvas_pVSfracresidP->Print("Plots/New/PpVSFrac_resid.png");

    TCanvas *mccanvas_nhitsVSfracresid = new TCanvas("mccanvas_nhitsVSfracresid","",1000,800);
    gStyle->SetOptStat(0);
    nhitsVSfrac_resid->SetTitle("Momentum fractional residuals VS n(hits) ;n(hits);(p_{reco}-p_{true})/p_{true}");
    nhitsVSfrac_resid->Draw("COLZ");
    //mccanvas_nhitsVSfracresid->Print("Plots/New/nhitsVSFrac_resid_prange.png");

    TCanvas *mccanvas_nhitsVSfracresidP = new TCanvas("mccanvas_nhitsVSfracresidP","",1000,800);
    gStyle->SetOptStat(0);
    PnhitsVSfrac_resid->SetTitle("P(Residual|n(hits));n(hits);(p_{reco}-p_{true})/p_{true}");
    PnhitsVSfrac_resid->Draw("COLZ");
    //mccanvas_nhitsVSfracresidP->Print("Plots/New/PnhitsVSFrac_resid_prange.png");
    */
    
}