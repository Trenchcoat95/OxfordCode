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
#include "garana.h"

//   In a ROOT session, you can do:
  //      root> .L garana.C
  //      root> garana t
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.l2g_Trackmatch();       // Use the trackmatching function

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

void FitSlicesYCustom(TH2F * HSTD, TF1* Func, TObjArray &HPROB, std::string type, bool edge)
{
    TH1D* h0= new TH1D("h0", "", HSTD->GetNbinsX(), HSTD->GetXaxis()->GetXmin(), HSTD->GetXaxis()->GetXmax());
    TH1D* h1= new TH1D("h1", "", HSTD->GetNbinsX(), HSTD->GetXaxis()->GetXmin(), HSTD->GetXaxis()->GetXmax());
    TH1D* h2= new TH1D("h2", "", HSTD->GetNbinsX(), HSTD->GetXaxis()->GetXmin(), HSTD->GetXaxis()->GetXmax());
    TH1D* h3= new TH1D("h3", "", HSTD->GetNbinsX(), HSTD->GetXaxis()->GetXmin(), HSTD->GetXaxis()->GetXmax());
    TH1D* h4= new TH1D("h4", "", HSTD->GetNbinsX(), HSTD->GetXaxis()->GetXmin(), HSTD->GetXaxis()->GetXmax());
    TH1D* h5= new TH1D("h5", "", HSTD->GetNbinsX(), HSTD->GetXaxis()->GetXmin(), HSTD->GetXaxis()->GetXmax());
    TH1D* hchi= new TH1D("hchi", "", HSTD->GetNbinsX(), HSTD->GetXaxis()->GetXmin(), HSTD->GetXaxis()->GetXmax());
    TCanvas *mccanvas_Slices = new TCanvas("mccanvas_Slices","",1000,800);
    gStyle->SetOptFit(1);

    for (Int_t x=HSTD->GetNbinsX(); x>=1; x--)
    {
        TH1D * projY = HSTD->ProjectionY("projY",x-1,x,"");
        if (projY->GetEntries()>1000)
        {
            if(x==HSTD->GetNbinsX()) Func->SetParameters(projY->GetEntries(),projY->GetMean(),projY->GetRMS(),0.5,projY->GetRMS(),projY->GetRMS());
            //if(x==HSTD->GetNbinsX()) Func->SetParameters(projY->GetEntries(),projY->GetMean(),projY->GetRMS(),0.5,projY->GetMean(),projY->GetRMS());
            else Func->SetParameter(0,projY->GetEntries());
            projY->Fit(Func->GetName(),"Q");
            //std::cout<<"Bin: "<<x<<" Entries: "<<projY->GetEntries()<<" Param 4: "<<Func->GetParameter(4)<<std::endl;
        }
    }

    for (Int_t x=1; x<=HSTD->GetNbinsX(); x++)
    {
        TH1D * projY = HSTD->ProjectionY("projY",x-1,x,"");
        if (projY->GetEntries()>1000)
        {
            //if(x==1) Func->SetParameters(projY->GetEntries(),projY->GetMean(),projY->GetRMS(),0.5,projY->GetRMS(),projY->GetRMS());
            if (x!=1) Func->SetParameter(0,projY->GetEntries());
            projY->Fit(Func->GetName(),"Q");

            projY->Draw();
            std::string Formula0 = "0.39894228040143*"+std::to_string(projY->GetBinWidth(0))+"*([0]/[2])*(exp(-0.5*((x-[1])/[2])^2))";
            TF1 *gauss0 = new TF1("gauss1",Formula0.c_str(),-0.4,0.4);
            gauss0->SetParameters(Func->GetParameter(0),Func->GetParameter(1),Func->GetParameter(2));
            gauss0->SetLineColor(kBlue);
            gauss0->SetLineStyle(9);
            gauss0->Draw("SAME");
            TF1 *gauss2 = new TF1("gauss2",Formula0.c_str(),-0.4,0.4);
            //gauss2->SetParameters(double_gauss->GetParameter(0)*double_gauss->GetParameter(3),double_gauss->GetParameter(1)+abs(double_gauss->GetParameter(4)),double_gauss->GetParameter(5));
            gauss2->SetParameters(Func->GetParameter(0)*Func->GetParameter(3),Func->GetParameter(4),Func->GetParameter(5));
            gauss2->SetLineColor(kBlue);
            gauss2->SetLineStyle(9);
            gauss2->Draw("SAME");

            std::string xstring;
            if(edge)  xstring ="Edge/Slices/"+type+"_SlicesFit_"+std::to_string(x)+"_edge.png";
            else  xstring ="Standard/Slices/"+type+"_SlicesFit_"+std::to_string(x)+".png";
            //std::cout<<xstring.c_str()<<std::endl;
            mccanvas_Slices->Print(xstring.c_str());
            if(Func->GetChisquare()/Func->GetNDF()<100)
            {
                h0->SetBinContent(x,(Func->GetParameter(0)));
                h0->SetBinError(x,Func->GetParError(0));
                h1->SetBinContent(x,Func->GetParameter(1));
                if(Func->GetParError(1)<100) h1->SetBinError(x,Func->GetParError(1));
                h2->SetBinContent(x,(Func->GetParameter(2)));
                if(Func->GetParError(2)<100) h2->SetBinError(x,Func->GetParError(2));
                h3->SetBinContent(x,(Func->GetParameter(3)));
                if(Func->GetParError(3)<100) h3->SetBinError(x,Func->GetParError(3));
                h4->SetBinContent(x,Func->GetParameter(4));
                if(Func->GetParError(4)<100) h4->SetBinError(x,Func->GetParError(4));
                h5->SetBinContent(x,(Func->GetParameter(5)));
                if(Func->GetParError(5)<100) h5->SetBinError(x,Func->GetParError(5));
                hchi->SetBinContent(x,Func->GetChisquare()/Func->GetNDF());
            }
        }
    }
    HPROB.Add(h0);
    HPROB.Add(h1);
    HPROB.Add(h2);
    HPROB.Add(h3);
    HPROB.Add(h4);
    HPROB.Add(h5);
    HPROB.Add(hchi);
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

    gROOT->SetBatch(); ///Run in batch mode


    #pragma region "Plot Declaration"

        /////////////////////////////////////1D plots of the three main variables and 2D correlations
        TH1F* hnhits_RecoTracks_Sample = new TH1F("hnhits_RecoTracks_Sample", "nhits distribution for all muon tracks in GAr", 20, 0, 20);
        TH1F* hptrueStart_InGAr_Sample = new TH1F("hptrueStart_InGAr_Sample", "ptrue distribution for the sample of muons traversing GAr", 50, 0, 8);
        TH1F* htheta_InGAr_Sample = new TH1F("hmu_theta_InGAr_Sample", "Thetatrue distribution for the sample of muons traversing GAr",80, 0, 80.0);

        TH2F* hnhitsVSptrueStart_InGAr_Sample = new TH2F("hptrueStartVSnhits_InGAr_Sample ", "ptrue distribution for the sample of muons traversing GAr",50, 0, 8, 20, 0, 20);
        TH2F* hnhitsVSTheta_InGAr_Sample = new TH2F("hthetaVSnhits_InGAr_Sample ", "ptrue distribution for the sample of muons traversing GAr",80, 0, 80.0, 20, 0, 20);
        TH2F* PnhitsVSptrueStart_InGAr_Sample = new TH2F("PptrueStartVSnhits_InGAr_Sample ", "ptrue distribution for the sample of muons traversing GAr", 50, 0, 8, 20, 0, 20);
        TH2F* PnhitsVSTheta_InGAr_Sample = new TH2F("PthetaVSnhits_InGAr_Sample ", "ptrue distribution for the sample of muons traversing GAr",80, 0, 80.0, 20, 0, 20);

        /////////////////////////////////////1D Distribution of matched tracks
        TH1F* hnhits_RecoTracks_Matched = new TH1F("hnhits_RecoTracks_Matched", "nhits distribution for all matched muon tracks in GAr", 20, 0, 20);
        TH1F* hptrueStart_InGAr_Matched = new TH1F("hptrueStart_InGAr_Matched", "ptrue distribution for all matched muon tracks", 50, 0, 8);
        TH1F* htheta_InGAr_Matched = new TH1F("htheta_InGAr_Matched", "Thetatrue distribution for all matched muon tracks", 80, 0, 80.0);

        /////////////////////////////////////Charge resolution
        TH1F* hmu_good_recocharge_p_start = new TH1F("hmu_good_recocharge_p_start", "Reconstructed charge resolution", 50, 0, 8.0);
        TH1F* hmu_all_recocharge_p_start = new TH1F("hmu_all_recocharge_p_start", "Reconstructed charge resolution", 50, 0, 8.0);
        TH1F* hmu_good_recocharge_nhits = new TH1F("hmu_good_recocharge_nhits", "Reconstructed charge resolution", 16, 2, 18);
        TH1F* hmu_all_recocharge_nhits = new TH1F("hmu_all_recocharge_nhits", "Reconstructed charge resolution", 16, 2, 18);
        TH1F* hmu_good_recocharge_Theta = new TH1F("hmu_good_recocharge_Theta", "Reconstructed charge resolution", 55, 0, 55);
        TH1F* hmu_all_recocharge_Theta = new TH1F("hmu_all_recocharge_Theta", "Reconstructed charge resolution", 55, 0, 55);
        
        /////////////////////////////////////1D Fractional Residual distribution for the reconstructed muons
        TH1F* frac_resid = new TH1F("frac_resid", "Fractional residuals", 150, -0.4, 0.4);

        /////////////////////////////////////Fractional Residuals 2D   
        TH2F* pVSfrac_resid = new TH2F("pVSfrac_resid", "Fractional residuals",50,0,8, 50, -0.3, 0.3);
        TH2F* PpVSfrac_resid = new TH2F("PpVSfrac_resid", "P(Residual|p)",50,0,8, 50, -0.3, 0.3);
        TH2F* nhitsVSfrac_resid = new TH2F("nhitsVSfrac_resid", "Fractional residuals",16,2,18, 50, -0.3, 0.3);
        TH2F* PnhitsVSfrac_resid = new TH2F("PnhitsVSfrac_resid", "P(Residual|Theta)",16,2,18, 50, -0.3, 0.3);
        TH2F* ThetaVSfrac_resid = new TH2F("ThetaVSfrac_resid", "Fractional residuals",50,0,50, 50, -0.3, 0.3);
        TH2F* PThetaVSfrac_resid = new TH2F("PThetaVSfrac_resid", "P(Residual|Theta)",50,0,50, 50, -0.3, 0.3);

        /////////////////////////////////////Slice fit plots
        TH2F* pVSfrac_resid_SL = new TH2F("pVSfrac_residSL", "Fractional residuals",7,1,8, 150, -0.4, 0.4);
        TH2F* nhitsVSfrac_resid_SL = new TH2F("nhitsVSfrac_residSL", "Fractional residuals",14,2,16, 150, -0.4, 0.4);
        TH2F* ThetaVSfrac_resid_SL = new TH2F("ThetaVSfrac_residSL", "Fractional residuals",8,0,32, 150, -0.4, 0.4);
        TObjArray pSlices;
        TObjArray nhitsSlices;
        TObjArray ThetaSlices;

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

        int nTracks = TrackStartX->size();
        
        //Cycle over MC Particle
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

               //Cycle over all trajectory points, Find the ones corresponding to the current MC Particle and Find first point in GAr
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

                //Fill the MC Particle in GAr distributions and Fill the nhits per track in GAr distribution
                if(InGAr!=0) 
                {
                    htheta_InGAr_Sample->Fill(ThetaNuMu * (180.0/3.141592653589793238463));
                    hptrueStart_InGAr_Sample->Fill(pSt);
                    for (int iTrack=0; iTrack<nTracks; ++iTrack) 
                    {
                        hnhits_RecoTracks_Sample->Fill(NTPCClustersOnTrack->at(iTrack));
                        hnhitsVSptrueStart_InGAr_Sample->Fill(pSt,NTPCClustersOnTrack->at(iTrack));
                        hnhitsVSTheta_InGAr_Sample->Fill(ThetaNuMu * (180.0/3.141592653589793238463),NTPCClustersOnTrack->at(iTrack));
                    }
                }

               //TRACK MAtching
               
               TVector3 MCpart(px,py,pz);
               TVector3 MCpartHat = MCpart.Unit();
               double cosTmax = -1000;
               int iRECOpart = -1;
               TVector3 RECO_P;
               int RECO_CHARGE;
               int whichEnd = -1;

               if(InGAr!=0)
               {
                    
                    //looping over all reconstructed tracks
                    for (int iTrack=0; iTrack<nTracks; ++iTrack) 
                    {

                            //hnhits_RecoTracks_Sample->Fill(NTPCClustersOnTrack->at(iTrack));
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
                    pVSfrac_resid->Fill(pSt,(preco-p)/p);
                    nhitsVSfrac_resid->Fill(NTPCClustersOnTrack->at(iRECOpart),(preco-p)/p);
                    ThetaVSfrac_resid->Fill(ThetaNuMu* (180.0/3.141592653589793238463),(preco-p)/p);
                    pVSfrac_resid_SL->Fill(pSt,(preco-p)/p);
                    nhitsVSfrac_resid_SL->Fill(NTPCClustersOnTrack->at(iRECOpart),(preco-p)/p);
                    ThetaVSfrac_resid_SL->Fill(ThetaNuMu* (180.0/3.141592653589793238463),(preco-p)/p);

                    //Fill charge reconstruction resolution
                    hmu_all_recocharge_p_start->Fill(pSt);
                    hmu_all_recocharge_nhits->Fill(NTPCClustersOnTrack->at(iRECOpart));
                    hmu_all_recocharge_Theta->Fill(ThetaNuMu* (180.0/3.141592653589793238463));
                    if((PDG->at(j)==13 && RECO_CHARGE==-1) || (PDG->at(j)==-13 && RECO_CHARGE==1)) 
                    {
                        hmu_good_recocharge_p_start->Fill(pSt);
                        hmu_good_recocharge_nhits->Fill(NTPCClustersOnTrack->at(iRECOpart));
                        hmu_good_recocharge_Theta->Fill(ThetaNuMu* (180.0/3.141592653589793238463));
                    }
                    
                    
                }

               
           }
        }

        if(showprog==true) std::cout << std::string(strprog.length(),'\b')<<"\b";
        

    }


    GetProbabilityPlot(pVSfrac_resid,PpVSfrac_resid);
    GetProbabilityPlot(nhitsVSfrac_resid,PnhitsVSfrac_resid);
    GetProbabilityPlot(ThetaVSfrac_resid,PThetaVSfrac_resid);
    GetProbabilityPlot(hnhitsVSptrueStart_InGAr_Sample,PnhitsVSptrueStart_InGAr_Sample);
    GetProbabilityPlot(hnhitsVSTheta_InGAr_Sample,PnhitsVSTheta_InGAr_Sample);

    /////Define fitting function for Fractional Residual Plots
    //std::string Formula = "0.39894228040143*"+std::to_string(frac_resid->GetBinWidth(0))+"*([0]/[2])*(exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-([1]+abs([4])))/[5])^2)*([2]/[5]))";
    std::string Formula = "0.39894228040143*"+std::to_string(frac_resid->GetBinWidth(0))+"*([0]/[2])*(exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-([1]+[4]))/[5])^2)*([2]/[5]))";
    TF1 *double_gauss = new TF1("double_gauss",Formula.c_str(),-0.4,0.4);

    /////Construct fit evolution plots
    FitSlicesYCustom(pVSfrac_resid_SL,double_gauss,pSlices,"p",edge);
    FitSlicesYCustom(nhitsVSfrac_resid_SL,double_gauss,nhitsSlices,"nhits",edge);
    FitSlicesYCustom(ThetaVSfrac_resid_SL,double_gauss,ThetaSlices,"Theta",edge);
    
    #pragma region "Plotting"

        ///////////////////////////////////////////////1D distributions for three main variables
        TCanvas *mccanvasnhits = new TCanvas("mccanvasnhits","",1000,800);
        hnhits_RecoTracks_Sample->SetTitle("number of hits distribution for all muon tracks in GAr;n_{hits};n");
        hnhits_RecoTracks_Sample->SetMaximum(48000);
        hnhits_RecoTracks_Sample->Draw();
        if (edge) mccanvasnhits->Print("Edge/nhits1D_edge.png");
        else mccanvasnhits->Print("Standard/nhits1D.png");

        TCanvas *mccanvasp = new TCanvas("mccanvasp","",1000,800);
        hptrueStart_InGAr_Sample->SetTitle("p_{true} distribution for the sample of muons traversing GAr;p_{true}^{St}(GeV/c);n");
        hptrueStart_InGAr_Sample->SetMaximum(23000);
        hptrueStart_InGAr_Sample->Draw();
        if (edge) mccanvasp->Print("Edge/ptrue1D_edge.png");
        else mccanvasp->Print("Standard/ptrue1D.png");

        TCanvas *mccanvasTheta = new TCanvas("mccanvasTheta","",1000,800);
        htheta_InGAr_Sample->SetTitle("#theta_{true} distribution for the sample of muons traversing GAr;#theta_{true}^{St}(deg);n");
        htheta_InGAr_Sample->Draw();
        htheta_InGAr_Sample->SetMaximum(22000);
        if (edge) mccanvasTheta->Print("Edge/Thetatrue1D_edge.png");
        else mccanvasTheta->Print("Standard/Thetatrue1D.png");

        ///////////////////////////////////////////////2D distributions for three main variables
        TCanvas *mccanvaspVSnhits = new TCanvas("mccanvaspVSnhits","",1000,800);
        hnhitsVSptrueStart_InGAr_Sample->SetTitle("n_{hits} VS p_{true}^{St};p_{true}^{St}(GeV/c);n_{hits}");
        hnhitsVSptrueStart_InGAr_Sample->Draw("COLZ");
        if (edge) mccanvaspVSnhits->Print("Edge/nhitsVSp_edge.png");
        else mccanvaspVSnhits->Print("Standard/nhitsVSp.png");

        TCanvas *mccanvasPpVSnhits = new TCanvas("mccanvasPpVSnhits","",1000,800);
        PnhitsVSptrueStart_InGAr_Sample->SetTitle("P(n_{hits}|p_{true}^{St});p_{true}^{St}(GeV/c);n_{hits}");
        PnhitsVSptrueStart_InGAr_Sample->Draw("COLZ");
        if (edge) mccanvasPpVSnhits->Print("Edge/PnhitsVSp_edge.png");
        else mccanvasPpVSnhits->Print("Standard/PnhitsVSp.png");

        TCanvas *mccanvasThetaVSnhits = new TCanvas("mccanvasThetaVSnhits","",1000,800);
        hnhitsVSTheta_InGAr_Sample->SetTitle("n_{hits} VS #theta_{true}^{St} ;#theta_{true}^{St}(deg);n_{hits}");
        hnhitsVSTheta_InGAr_Sample->Draw("COLZ");
        if (edge) mccanvasThetaVSnhits->Print("Edge/nhitsVSTheta_edge.png");
        else mccanvasThetaVSnhits->Print("Standard/nhitsVSTheta.png");

        TCanvas *mccanvasPThetaVSnhits = new TCanvas("mccanvasPThetaVSnhits","",1000,800);
        PnhitsVSTheta_InGAr_Sample->SetTitle("P(n_{hits}|#theta_{true}^{St}) ;#theta_{true}^{St}(deg);n_{hits}");
        PnhitsVSTheta_InGAr_Sample->Draw("COLZ");
        if (edge) mccanvasPThetaVSnhits->Print("Edge/PnhitsVSTheta_edge.png");
        else mccanvasPThetaVSnhits->Print("Standard/PnhitsVSTheta.png");

        TCanvas *mccanvas_pVSnhitsProj = new TCanvas("mccanvas_pVSnhitsProj","",1000,800);
        TH1D * pVSnhits_projY = hnhitsVSptrueStart_InGAr_Sample->ProjectionY("hnhitsVSptrueStart_InGAr_Sample",10,10,"");
        pVSnhits_projY->Draw();
        if (edge) mccanvas_pVSnhitsProj->Print("Edge/pVSpVSnhits_projY_edge.png");
        else  mccanvas_pVSnhitsProj->Print("Standard/pVSnhits_projY.png");
    
        ///////////////////////////////////////////////1D distributions of matched tracks for three main variables
        TCanvas *mccanvasnhitsm = new TCanvas("mccanvasnhitsm","",1000,800);
        hnhits_RecoTracks_Matched->SetTitle("nhits distribution of matched muon tracks in GAr;nhits;n");
        hnhits_RecoTracks_Matched->Draw();
        if (edge) mccanvasnhitsm->Print("Edge/nhits1D_matched_edge.png");
        else mccanvasnhitsm->Print("Standard/nhits1D_matched.png");

        TCanvas *mccanvaspm = new TCanvas("mccanvaspm","",1000,800);
        hptrueStart_InGAr_Matched->SetTitle("p_{true} distribution for the sample of matched muon tracks;p_{true}^{St}(GeV/c);n");
        hptrueStart_InGAr_Matched->Draw();
        hptrueStart_InGAr_Matched->SetMaximum(23000);
        if(edge) mccanvaspm->Print("Edge/ptrue1D_matched_edge.png");
        else mccanvaspm->Print("Standard/ptrue1D_matched.png");

        TCanvas *mccanvasThetam = new TCanvas("mccanvasThetam","",1000,800);
        htheta_InGAr_Matched->SetTitle("#theta_{true} distribution for the sample of matched muon tracks;#theta_{true}^{St}(deg);n");
        htheta_InGAr_Matched->Draw();
        htheta_InGAr_Matched->SetMaximum(22000);
        if(edge) mccanvasThetam->Print("Edge/Thetatrue1D_matched_edge.png");
        else mccanvasThetam->Print("Standard/Thetatrue1D_matched.png");

        ///////////////////////////////////////////////Efficiency plots
        TCanvas *mccanvasHitsEff = new TCanvas("mccanvasHitsEff","",1000,800);
        gStyle->SetOptStat(0);
        hnhits_RecoTracks_Matched->SetTitle("Reconstruction (track matching) efficiency as function of number of hits;n_{hits};n(matched)/n(Reco Tracks)");
        hnhits_RecoTracks_Matched->Divide(hnhits_RecoTracks_Sample);
        hnhits_RecoTracks_Matched->SetMinimum(0);
        hnhits_RecoTracks_Matched->SetMaximum(0.8);
        hnhits_RecoTracks_Matched->Draw();
        if (edge) mccanvasHitsEff->Print("Edge/nhits1D_Efficiency_edge.png");
        else mccanvasHitsEff->Print("Standard/nhits1D_Efficiency.png");

        TCanvas *mccanvasPEff = new TCanvas("mccanvasPEff","",1000,800);
        gStyle->SetOptStat(0);
        hptrueStart_InGAr_Matched->SetTitle("Reconstruction (track matching) efficiency as function of initial momentum;p_{true}^{St} [GeV/c];n(matched)/n(MC Tracks in GAr)");
        hptrueStart_InGAr_Matched->Divide(hptrueStart_InGAr_Sample);
        hptrueStart_InGAr_Matched->SetMinimum(0);
        hptrueStart_InGAr_Matched->SetMaximum(0.7);
        hptrueStart_InGAr_Matched->Draw();
        if(edge) mccanvasPEff->Print("Edge/ptrue1D_Efficiency_edge.png");
        else mccanvasPEff->Print("Standard/ptrue1D_Efficiency.png");

        TCanvas *mccanvasAngleEff = new TCanvas("mccanvasAngleEff","",1000,800);
        gStyle->SetOptStat(0);
        htheta_InGAr_Matched->SetTitle("Reconstruction (track matching) efficiency as function of initial momentum;#theta_{true}^{St} [deg];n(matched)/n(MC Tracks in GAr)");
        htheta_InGAr_Matched->Divide(htheta_InGAr_Sample);
        htheta_InGAr_Matched->SetMinimum(0);
        htheta_InGAr_Matched->SetMaximum(0.7);
        htheta_InGAr_Matched->Draw();
        if(edge) mccanvasAngleEff->Print("Edge/Thetatrue1D_Efficiency_edge.png");
        else mccanvasAngleEff->Print("Standard/Thetatrue1D_Efficiency.png");

        /////////////////////////////////////////////////////Charge Resolution
        TCanvas *mccanvascharge = new TCanvas("mccanvascharge","",1000,800);
        gStyle->SetOptStat(0);
        hmu_good_recocharge_p_start->SetTitle("Charge reconstruction resolution;p_{true}^{St}[GeV/c];n(correct)/n(matched)");
        hmu_good_recocharge_p_start->Sumw2();
        hmu_good_recocharge_p_start->Divide(hmu_good_recocharge_p_start,hmu_all_recocharge_p_start,1,1,"B");
        hmu_good_recocharge_p_start->SetMinimum(0.5);
        hmu_good_recocharge_p_start->Draw();
        if (edge) mccanvascharge->Print("Edge/ChargeRes_ErrorBars_p_edge.png");
        else mccanvascharge->Print("Standard/ChargeRes_ErrorBars_p.png");

        TCanvas *mccanvaschargenhits = new TCanvas("mccanvaschargenhits","",1000,800);
        gStyle->SetOptStat(0);
        hmu_good_recocharge_nhits->SetTitle("Charge reconstruction resolution;n_{hits};n(correct)/n(matched)");
        hmu_good_recocharge_nhits->Sumw2();
        hmu_good_recocharge_nhits->Divide(hmu_good_recocharge_nhits,hmu_all_recocharge_nhits,1,1,"B");
        hmu_good_recocharge_nhits->SetMinimum(0.5);
        hmu_good_recocharge_nhits->Draw();
        if (edge) mccanvaschargenhits->Print("Edge/ChargeRes_ErrorBars_nhits_edge.png");
        else mccanvaschargenhits->Print("Standard/ChargeRes_ErrorBars.png");

        TCanvas *mccanvaschargeTheta = new TCanvas("mccanvaschargeTheta","",1000,800);
        gStyle->SetOptStat(0);
        hmu_good_recocharge_Theta->SetTitle("Charge reconstruction resolution;#theta_{true}^{St}[deg];n(correct)/n(matched)");
        hmu_good_recocharge_Theta->Sumw2();
        hmu_good_recocharge_Theta->Divide(hmu_good_recocharge_Theta,hmu_all_recocharge_Theta,1,1,"B");
        hmu_good_recocharge_Theta->SetMinimum(0.5);
        hmu_good_recocharge_Theta->Draw();
        if (edge) mccanvaschargeTheta->Print("Edge/ChargeRes_ErrorBars_Theta_edge.png");
        else mccanvaschargeTheta->Print("Standard/ChargeRes_ErrorBars_Theta.png");

        ////////////////////////////////////////////////////Fractional Residual 1D plot
        TCanvas *mccanvas_fracresid = new TCanvas("mccanvas_fracresid","",1000,800);
        gStyle->SetOptStat(1);
        frac_resid->SetTitle("Momentum fractional residuals (Double Gauss Fit);(p_{reco}-p_{true})/p_{true};n(#mu in GAr-Lite)");
        double_gauss->SetParameters(frac_resid->GetEntries(),frac_resid->GetMean(),frac_resid->GetRMS(),0.5,frac_resid->GetRMS(),frac_resid->GetRMS());
        //double_gauss->SetParameters(frac_resid->GetEntries(),frac_resid->GetMean(),frac_resid->GetRMS(),0.5,frac_resid->GetMean(),frac_resid->GetRMS());
        frac_resid->Fit("double_gauss");
        gStyle->SetOptFit(1);
        frac_resid->Draw();
        std::string Formula1 = "0.39894228040143*"+std::to_string(frac_resid->GetBinWidth(0))+"*([0]/[2])*(exp(-0.5*((x-[1])/[2])^2))";
        TF1 *gauss1 = new TF1("gauss1",Formula1.c_str(),-0.4,0.4);
        gauss1->SetParameters(double_gauss->GetParameter(0),double_gauss->GetParameter(1),double_gauss->GetParameter(2));
        gauss1->SetLineColor(kBlue);
        gauss1->SetLineStyle(9);
        gauss1->Draw("SAME");
        TF1 *gauss2 = new TF1("gauss2",Formula1.c_str(),-0.4,0.4);
        //gauss2->SetParameters(double_gauss->GetParameter(0)*double_gauss->GetParameter(3),double_gauss->GetParameter(1)+abs(double_gauss->GetParameter(4)),double_gauss->GetParameter(5));
        gauss2->SetParameters(double_gauss->GetParameter(0)*double_gauss->GetParameter(3),double_gauss->GetParameter(4),double_gauss->GetParameter(5));
        gauss2->SetLineColor(kBlue);
        gauss2->SetLineStyle(9);
        gauss2->Draw("SAME");
        if (edge) mccanvas_fracresid->Print("Edge/Frac_resid_two_gauss_fullspectrum_edge.png");
        else mccanvas_fracresid->Print("Standard/Frac_resid_two_gauss_fullspectrum.png");

        ///////////////////////////////////////////////////////Fractional ResidualVSptrue
        TCanvas *mccanvas_pVSfracresid = new TCanvas("mccanvas_pVSfracresid","",1000,800);
        gStyle->SetOptStat(0);
        pVSfrac_resid->SetTitle("Momentum fractional residuals VS p_{true}^{St};p_{true}^{St}(GeV/c);(p_{reco}-p_{true})/p_{true}");
        pVSfrac_resid->Draw("COLZ");
        if (edge) mccanvas_pVSfracresid->Print("Edge/pVSFrac_resid_edge.png");
        else  mccanvas_pVSfracresid->Print("Standard/pVSFrac_resid.png");

        TCanvas *mccanvas_pVSfracresidP = new TCanvas("mccanvas_pVSfracresidP","",1000,800);
        gStyle->SetOptStat(0);
        PpVSfrac_resid->SetTitle("P(Residual|p_{true}^{St});p_{true}^{St}(GeV/c);(p_{reco}-p_{true})/p_{true}");
        PpVSfrac_resid->Draw("COLZ");
        if (edge) mccanvas_pVSfracresidP->Print("Edge/PpVSFrac_resid_edge.png");
        else  mccanvas_pVSfracresidP->Print("Standard/PpVSFrac_resid.png");

        TCanvas *mccanvas_pVSfracresidProj = new TCanvas("mccanvas_pVSfracresidProj","",1000,800);
        TH1D * pVSfrac_resid_projY = pVSfrac_resid->ProjectionY("pVSfrac_resid_projY",10,11,"");
        pVSfrac_resid_projY->Draw();
        if (edge) mccanvas_pVSfracresidProj->Print("Edge/pVSFrac_resid_proj_edge.png");
        else  mccanvas_pVSfracresidProj->Print("Standard/pVSFrac_resid_proj.png");

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
        ThetaVSfrac_resid->SetTitle("Momentum fractional residuals VS #theta_{true}^{St};#theta_{true}^{St};(p_{reco}-p_{true})/p_{true}");
        ThetaVSfrac_resid->Draw("COLZ");
        if (edge) mccanvas_ThetaVSfracresid->Print("Edge/ThetaVSFrac_resid_edge.png");
        else  mccanvas_ThetaVSfracresid->Print("Standard/ThetaVSFrac_resid.png");

        TCanvas *mccanvas_ThetaVSfracresidP = new TCanvas("mccanvas_ThetaVSfracresidP","",1000,800);
        gStyle->SetOptStat(0);
        PThetaVSfrac_resid->SetTitle("P(Residual|#theta_{true}^{St});#theta_{true}^{St};(p_{reco}-p_{true})/p_{true}");
        PThetaVSfrac_resid->Draw("COLZ");
        if (edge) mccanvas_ThetaVSfracresidP->Print("Edge/PThetaVSFrac_resid_edge.png");
        else  mccanvas_ThetaVSfracresidP->Print("Standard/PThetaVSFrac_resid.png");

        /////////////////////////////////////////////////////////Fit evolution Fractional Residual in p slices
        
        TCanvas *mccanvas_param0p = new TCanvas("mccanvas_param0p","",1000,800);
        TH1D* hpSlices0=(TH1D*) pSlices.At(0);
        hpSlices0->SetMinimum(0);
        hpSlices0->SetMaximum(65000);
        hpSlices0->Draw();
        if (edge) mccanvas_param0p->Print("Edge/pslices_0_edge.png");
        else  mccanvas_param0p->Print("Standard/pslices_0.png");

        TCanvas *mccanvas_param1p = new TCanvas("mccanvas_param1p","",1000,800);
        TH1D* hpSlices1=(TH1D*) pSlices.At(1);
        hpSlices1->SetMinimum(-0.03);
        hpSlices1->SetMaximum(0.0);
        hpSlices1->Draw();
        if (edge) mccanvas_param1p->Print("Edge/pslices_1_edge.png");
        else  mccanvas_param1p->Print("Standard/pslices_1.png");

        TCanvas *mccanvas_param2p = new TCanvas("mccanvas_param2p","",1000,800);
        TH1D* hpSlices2=(TH1D*) pSlices.At(2);
        hpSlices2->SetMinimum(0.014);
        hpSlices2->SetMaximum(0.034);
        hpSlices2->Draw();
        if (edge) mccanvas_param2p->Print("Edge/pslices_2_edge.png");
        else  mccanvas_param2p->Print("Standard/pslices_2.png");

        TCanvas *mccanvas_param3p = new TCanvas("mccanvas_param3p","",1000,800);
        TH1D* hpSlices3=(TH1D*) pSlices.At(3);
        hpSlices3->SetMinimum(0.12);
        hpSlices3->SetMaximum(0.86);
        hpSlices3->Draw();
        if (edge) mccanvas_param3p->Print("Edge/pslices_3_edge.png");
        else  mccanvas_param3p->Print("Standard/pslices_3.png");

        TCanvas *mccanvas_param4p = new TCanvas("mccanvas_param4p","",1000,800);
        TH1D* hpSlices4=(TH1D*) pSlices.At(4);
        hpSlices4->SetMinimum(-0.17);
        hpSlices4->SetMaximum(0.012);
        hpSlices4->Draw();
        if (edge) mccanvas_param4p->Print("Edge/pslices_4_edge.png");
        else  mccanvas_param4p->Print("Standard/pslices_4.png");

        TCanvas *mccanvas_param5p = new TCanvas("mccanvas_param5p","",1000,800);
        TH1D* hpSlices5=(TH1D*) pSlices.At(5);
        hpSlices5->SetMinimum(0.03);
        hpSlices5->SetMaximum(0.13);
        hpSlices5->Draw();
        if (edge) mccanvas_param5p->Print("Edge/pslices_5_edge.png");
        else  mccanvas_param5p->Print("Standard/pslices_5.png");

        TCanvas *mccanvas_paramchip = new TCanvas("mccanvas_paramchip","",1000,800);
        TH1D* hpSlices6=(TH1D*) pSlices.At(6);
        hpSlices6->SetMinimum(0);
        hpSlices6->SetMaximum(10);
        hpSlices6->Draw();
        if (edge) mccanvas_paramchip->Print("Edge/pslices_chi_edge.png");
        else  mccanvas_paramchip->Print("Standard/pslices_chi.png");
        
        /////////////////////////////////////////////////////////Fit evolution Fractional Residual in nhits slices
        TF1 *InverseSqrt = new TF1("InverseSqrt","[0]*1/(sqrt(x))",6,11);

        TCanvas *mccanvas_param0nhits = new TCanvas("mccanvas_paramnhits","",1000,800);
        TH1D* hnhitsSlices0=(TH1D*) nhitsSlices.At(0);
        hnhitsSlices0->SetMinimum(0);
        hnhitsSlices0->SetMaximum(48000);
        hnhitsSlices0->Draw();
        if (edge) mccanvas_param0nhits->Print("Edge/nhitsslices_0_edge.png");
        else  mccanvas_param0nhits->Print("Standard/nhitsslices_0.png");

        TCanvas *mccanvas_param1nhits = new TCanvas("mccanvas_param1nhits","",1000,800);
        TH1D* hnhitsSlices1=(TH1D*) nhitsSlices.At(1);
        hnhitsSlices1->SetMinimum(-0.02);
        hnhitsSlices1->SetMaximum(-0.005);
        hnhitsSlices1->Draw();
        if (edge) mccanvas_param1nhits->Print("Edge/nhitsslices_1_edge.png");
        else  mccanvas_param1nhits->Print("Standard/nhitsslices_1.png");

        TCanvas *mccanvas_param2nhits = new TCanvas("mccanvas_param2nhits","",1000,800);
        TH1D* hnhitsSlices2=(TH1D*) nhitsSlices.At(2);
        gStyle->SetOptFit(1);
        hnhitsSlices2->Fit("InverseSqrt","Q","",6,11);    
        hnhitsSlices2->SetMinimum(0.01);
        hnhitsSlices2->SetMaximum(0.05);
        hnhitsSlices2->Draw();
        if (edge) mccanvas_param2nhits->Print("Edge/nhitsslices_2_edge.png");
        else  mccanvas_param2nhits->Print("Standard/nhitsslices_2.png");

        TCanvas *mccanvas_param3nhits = new TCanvas("mccanvas_param3nhits","",1000,800);
        TH1D* hnhitsSlices3=(TH1D*) nhitsSlices.At(3);
        hnhitsSlices3->SetMinimum(0);
        hnhitsSlices3->SetMaximum(1.6);
        hnhitsSlices3->Draw();
        if (edge) mccanvas_param3nhits->Print("Edge/nhitsslices_3_edge.png");
        else  mccanvas_param3nhits->Print("Standard/nhitsslices_3.png");

        TCanvas *mccanvas_param4nhits = new TCanvas("mccanvas_param4nhits","",1000,800);
        TH1D* hnhitsSlices4=(TH1D*) nhitsSlices.At(4);
        hnhitsSlices4->SetMinimum(-0.185);
        hnhitsSlices4->SetMaximum(0.015);
        hnhitsSlices4->Draw();
        if (edge) mccanvas_param4nhits->Print("Edge/nhitsslices_4_edge.png");
        else  mccanvas_param4nhits->Print("Standard/nhitsslices_4.png");

        TCanvas *mccanvas_param5nhits = new TCanvas("mccanvas_param5nhits","",1000,800);
        TH1D* hnhitsSlices5=(TH1D*) nhitsSlices.At(5);
        hnhitsSlices5->SetMinimum(0.05);
        hnhitsSlices5->SetMaximum(0.3);
        hnhitsSlices5->Fit("InverseSqrt","Q","",6,11);
        hnhitsSlices5->Draw();
        if (edge) mccanvas_param5nhits->Print("Edge/nhitsslices_5_edge.png");
        else  mccanvas_param5nhits->Print("Standard/nhitsslices_5.png");

        TCanvas *mccanvas_paramchinhits = new TCanvas("mccanvas_paramchinhits","",1000,800);
        TH1D* hnhitsSlices6=(TH1D*) nhitsSlices.At(6);
        hnhitsSlices6->SetMinimum(0);
        hnhitsSlices6->SetMaximum(10);
        hnhitsSlices6->Draw();
        if (edge) mccanvas_paramchinhits->Print("Edge/nhitsslices_chi_edge.png");
        else  mccanvas_paramchinhits->Print("Standard/nhitsslices_chi.png");

        /////////////////////////////////////////////////////////Fit evolution Fractional Residual in Theta slices
        
        TCanvas *mccanvas_param0Theta = new TCanvas("mccanvas_param0Theta","",1000,800);
        TH1D* hThetaSlices0=(TH1D*) ThetaSlices.At(0);
        hThetaSlices0->SetMinimum(0);
        hThetaSlices0->SetMaximum(60000);
        hThetaSlices0->Draw();
        if (edge) mccanvas_param0Theta->Print("Edge/Thetaslices_0_edge.png");
        else  mccanvas_param0Theta->Print("Standard/Thetaslices_0.png");

        TCanvas *mccanvas_param1Theta = new TCanvas("mccanvas_param1Theta","",1000,800);
        TH1D* hThetaSlices1=(TH1D*) ThetaSlices.At(1);
        hThetaSlices1->SetMinimum(-0.035);
        hThetaSlices1->SetMaximum(-0.005);
        hThetaSlices1->Draw();
        if (edge) mccanvas_param1Theta->Print("Edge/Thetaslices_1_edge.png");
        else  mccanvas_param1Theta->Print("Standard/Thetaslices_1.png");

        TCanvas *mccanvas_param2Theta = new TCanvas("mccanvas_param2Theta","",1000,800);
        TH1D* hThetaSlices2=(TH1D*) ThetaSlices.At(2);
        hThetaSlices2->SetMinimum(0.015);
        hThetaSlices2->SetMaximum(0.035);
        hThetaSlices2->Draw();
        if (edge) mccanvas_param2Theta->Print("Edge/Thetaslices_2_edge.png");
        else  mccanvas_param2Theta->Print("Standard/Thetaslices_2.png");

        TCanvas *mccanvas_param3Theta = new TCanvas("mccanvas_param3Theta","",1000,800);
        TH1D* hThetaSlices3=(TH1D*) ThetaSlices.At(3);
        hThetaSlices3->SetMinimum(0);
        hThetaSlices3->SetMaximum(1);
        hThetaSlices3->Draw();
        if (edge) mccanvas_param3Theta->Print("Edge/Thetaslices_3_edge.png");
        else  mccanvas_param3Theta->Print("Standard/Thetaslices_3.png");

        TCanvas *mccanvas_param4Theta = new TCanvas("mccanvas_param4Theta","",1000,800);
        TH1D* hThetaSlices4=(TH1D*) ThetaSlices.At(4);
        hThetaSlices4->SetMinimum(-0.17);
        hThetaSlices4->SetMaximum(0.02);
        hThetaSlices4->Draw();
        if (edge) mccanvas_param4Theta->Print("Edge/Thetaslices_4_edge.png");
        else  mccanvas_param4Theta->Print("Standard/Thetaslices_4.png");

        TCanvas *mccanvas_param5Theta = new TCanvas("mccanvas_param5Theta","",1000,800);
        TH1D* hThetaSlices5=(TH1D*) ThetaSlices.At(5);
        hThetaSlices5->SetMinimum(0);
        hThetaSlices5->SetMaximum(0.2);
        hThetaSlices5->Draw();
        if (edge) mccanvas_param5Theta->Print("Edge/Thetaslices_5_edge.png");
        else  mccanvas_param5Theta->Print("Standard/Thetaslices_5.png");

        TCanvas *mccanvas_paramchiTheta = new TCanvas("mccanvas_paramchiTheta","",1000,800);
        TH1D* hThetaSlices6=(TH1D*) ThetaSlices.At(6);
        hThetaSlices6->SetMinimum(0);
        hThetaSlices6->SetMaximum(13);
        hThetaSlices6->Draw();
        if (edge) mccanvas_paramchiTheta->Print("Edge/Thetaslices_chi_edge.png");
        else  mccanvas_paramchiTheta->Print("Standard/Thetaslices_chi.png");

        

    #pragma endregion

}