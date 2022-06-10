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
  //      root> .L l2g_trackmatch_t.C
  //      root> garana t
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.l2g_Trackmatch();       // Use the trackmatching function




void l2g_trackmatch_t2()
{

    TVector3 GArCenter(0,-150.473,1486); //Correct values
    bool edge=false;
    //TVector3 GArCenter(0,-68.287,1486);  //Edge sample values
    //bool edge=true;
    float GAr_r = 349.9;
    float GAr_L = 669.6;
    const int  nplanes = 6 ;
    double  Planes_X[] =   {-300,300,-300,300,-300,300,-300,300,-300,300,-300,300};
    double  Planes_Y[] =   {-283.5,-16.5,-322.5,22.5,-350,50,-375,75,-400,100,-400,100};
    double  Planes_Z[] =   {1179,1183,1199,1203,1219,1223,1239,1243,1339,1343,1539,1543};

    garana t;
    Int_t nentries = t.fChain->GetEntries();

    gROOT->SetBatch(); ///Run in batch mode


    

    bool showprog = false;  
    if(showprog==true) std::cout<<"Progress:  "<<std::endl;

    for (Int_t i=0; i<nentries; i++) 
    {
        t.fChain->GetEntry(i);
        int prog = 100*i/nentries;
        std::string strprog = std::to_string(prog);
        if(showprog==true) std::cout<<strprog<<"%";
        
        

        int nTracks = t.TrackStartX->size();
        int nMCTraj = t.PDG->size();
        
        //Cycle over MC Particle
        for (Int_t j=0; j<nMCTraj; j++) 
        {
           
               float px = 0;
               float py = 0;
               float pz = 0;
               float p=0;
               float x = 0;
               float y = 0;
               float z = 0;
               float pxSt = t.MCPStartPX->at(j);
               float pySt = t.MCPStartPY->at(j);
               float pzSt = t.MCPStartPZ->at(j);
               float pSt = sqrt(pxSt*pxSt+pySt*pySt+pzSt*pzSt);
               int ID = t.MCTrkID->at(j);
               TVector3 MCPpart(t.MCPStartPX->at(j),t.MCPStartPY->at(j),t.MCPStartPZ->at(j));
               
               int InGAr=0;
               
               //t.TrajMCPX->at(k)>Planes_X[p*2] && t.TrajMCPX->at(k)<Planes_X[p*2+1] && 
            //                   t.TrajMCPY->at(k)>Planes_Y[p*2] && t.TrajMCPY->at(k)<Planes_Y[p*2+1] && 

               //Cycle over all trajectory points, Find the ones corresponding to the current MC Particle and Find first point in GAr
               for (Int_t k=0; k<t.TrajMCPX->size(); k++) 
                {
                    if (t.TrajMCPTrackID->at(k)==ID)
                    {
                        float r = sqrt((t.TrajMCPY->at(k)-GArCenter.Y())*(t.TrajMCPY->at(k)-GArCenter.Y())+(t.TrajMCPZ->at(k)-GArCenter.Z())*(t.TrajMCPZ->at(k)-GArCenter.Z()));
                    
                        if(t.TrajMCPX->at(k)<0.5*GAr_L && t.TrajMCPX->at(k)>-0.5*GAr_L && r<GAr_r && InGAr==0)              
                        {
                        
                                px=t.TrajMCPPX->at(k);
                                py=t.TrajMCPPY->at(k);
                                pz=t.TrajMCPPZ->at(k);
                                p = sqrt(px*px+py*py+pz*pz);
                                x=t.TrajMCPX->at(k);
                                y=t.TrajMCPY->at(k);
                                z=t.TrajMCPZ->at(k);
                                InGAr++;
                                break;
                            
                        }
                            
                        
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
                            TVector3 RECO_P_beg(t.TrackStartPX->at(iTrack), t.TrackStartPY->at(iTrack),t.TrackStartPZ->at(iTrack));
                            TVector3 RECO_P_end(t.TrackEndPX->at(iTrack), t.TrackEndPY->at(iTrack),  t.TrackEndPZ->at(iTrack));
                            int RECO_CHARGE_beg=t.TrackStartQ->at(iTrack);
                            int RECO_CHARGE_end=t.TrackEndQ->at(iTrack);

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
                                    RECO_P	  = RECO_P_end;
                                    RECO_CHARGE = RECO_CHARGE_end;
                                }
                            }   
                    } // End loop over all Tracks

                    if (iRECOpart == -1) continue;

                    TVector3 vecX;
                    if ( whichEnd == 0 ) 
                    {
                        vecX.SetXYZ(t.TrackStartX->at(iRECOpart) -x,
                                    t.TrackStartY->at(iRECOpart) -y,
                                    t.TrackStartZ->at(iRECOpart) -z);
                    } 
                    else if (whichEnd==1)
                    {
                        vecX.SetXYZ(t.TrackEndX->at(iRECOpart)   -x,
                                    t.TrackEndY->at(iRECOpart)   -y,
                                    t.TrackEndZ->at(iRECOpart)   -z);
                    }
                    Float_t delX = vecX.Cross(MCpartHat).Mag();
                    if ( cosTmax <= 0.997 ) continue;
                    if (  delX   >= 3.0   ) continue;

                    //If we are here the track was matched

                std::cout<<"Track matched! Reco p: ("<<RECO_P.X()<<" ,"<<RECO_P.Y()<<" ,"<<RECO_P.Z()<<" ) | MC p: ("<<px<<" ,"<<py<<" ,"<<pz<<" )"<<std::endl;    
                if(whichEnd==0) std::cout<<"               Reco xyz: ("<<t.TrackStartX->at(iRECOpart)<<" ,"<<t.TrackStartY->at(iRECOpart)<<" ,"<<t.TrackStartZ->at(iRECOpart)<<" ) | MC p: ("<<x<<" ,"<<y<<" ,"<<z<<" )"<<std::endl; 
                if(whichEnd==1) std::cout<<"               Reco xyz: ("<<t.TrackEndX->at(iRECOpart)<<" ,"<<t.TrackEndY->at(iRECOpart)<<" ,"<<t.TrackEndZ->at(iRECOpart)<<" ) | MC xyz: ("<<x<<" ,"<<y<<" ,"<<z<<" )"<<std::endl;   
                    
                    
                }

         
           
        }

        if(showprog==true) std::cout << std::string(strprog.length(),'\b')<<"\b";
        
    }
    


    

}