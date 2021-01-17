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

void l2g_trackmatch()
{
    TChain chain("/anatree/GArAnaTree");
    chain.Add("/pnfs/dune/persistent/users/ebrianne/ProductionSamples/ND-LAr/nd_hall_dayone_lar_SPY_v2_wMuID/Anatree/neutrino/neutrino.nd_hall_dayone_lar_SPY_v2_wMuID.volArgonCubeActive.Ev973000.Ev973999.2037.anatree.root");

    //neutrino.nd_hall_dayone_lar_SPY_v2_wMuID.volArgonCubeActive.Ev973000.Ev973999.2037.anatree.root

    vector<int>     *PDG=0;
    vector<int>     *PDGMother=0;
    vector<int>     *MCTrkID=0;
    vector<int>     *TrajMCPTrackID=0;
    vector<float>   *MCPStartPX=0;
    vector<float>   *MCPStartPY=0;
    vector<float>   *MCPStartPZ=0;
    vector<float>   *MCPStartX=0;
    vector<float>   *MCPStartY=0;
    vector<float>   *MCPStartZ=0;
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
    TBranch *b_MCTrkID=0;
    TBranch *b_TrajMCPTrackID=0;
    TBranch *b_MCPStartPX=0;
    TBranch *b_MCPStartPY=0;
    TBranch *b_MCPStartPZ=0;
    TBranch *b_MCPStartX=0;
    TBranch *b_MCPStartY=0;
    TBranch *b_MCPStartZ=0;
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
    chain.SetBranchAddress("MCTrkID", &MCTrkID, &b_MCTrkID);
    chain.SetBranchAddress("TrajMCPTrackID", &TrajMCPTrackID, &b_TrajMCPTrackID);
    //chain.SetBranchAddress("MCPProc", &MCPProc, &b_MCPProc);
    chain.SetBranchAddress("MCPStartPX", &MCPStartPX, &b_MCPStartPX);
    chain.SetBranchAddress("MCPStartPY", &MCPStartPY, &b_MCPStartPY);
    chain.SetBranchAddress("MCPStartPZ", &MCPStartPZ, &b_MCPStartPZ);
    chain.SetBranchAddress("MCPStartX", &MCPStartX, &b_MCPStartX);
    chain.SetBranchAddress("MCPStartY", &MCPStartY, &b_MCPStartY);
    chain.SetBranchAddress("MCPStartZ", &MCPStartZ, &b_MCPStartZ);
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
    chain.SetBranchStatus("MCTrkID",1);
    chain.SetBranchStatus("TrajMCPTrackID",1);
    chain.SetBranchStatus("MCPStartPX",1);
    chain.SetBranchStatus("MCPStartPY",1);
    chain.SetBranchStatus("MCPStartPZ",1);
    chain.SetBranchStatus("MCPStartX",1);
    chain.SetBranchStatus("MCPStartY",1);
    chain.SetBranchStatus("MCPStartZ",1);

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

    TH1F* hmu_good_recocharge_p_start = new TH1F("hmu_good_recocharge_p_start", "Reconstructed charge resolution", 100, 0, 20.0);
    TH1F* hmu_all_recocharge_p_start = new TH1F("hmu_all_recocharge_p_start", "Reconstructed charge resolution", 100, 0, 20.0);
    

    for (Int_t i=0; i<nentries; i++) 
    {
        chain.GetEntry(i);
        

        for (Int_t j=0; j<PDG->size(); j++) 
        {
           if(PDG->at(j)==13 || PDG->at(j)==-13)
           {
               float px = MCPStartPX->at(j);
               float py = MCPStartPY->at(j);
               float pz = MCPStartPZ->at(j);
               float p = sqrt(px*px+py*py+pz*pz);
    
               

               float x = MCPStartX->at(j);
               float y = MCPStartY->at(j);
               float z = MCPStartZ->at(j);
               
               float Theta = acos(pz/p);

               //TRACK MAtching
               int nTracks = TrackStartX->size();
               TVector3 MCpart(MCPStartX->at(j),MCPStartY->at(j),MCPStartZ->at(j));
               TVector3 MCpartHat = MCpart.Unit();
               double cosTmax = -1000;
               int iRECOpart = -1;
               TVector3 RECO_P;
               int whichEnd = -1;
               for (int iTrack=0; iTrack<nTracks; ++iTrack) 
               {
				    TVector3 RECO_P_beg(TrackStartPX->at(iTrack),
									TrackStartPY->at(iTrack),TrackStartPZ->at(iTrack));
				    TVector3 RECO_P_end(TrackEndPX->at(iTrack),
									TrackEndPY->at(iTrack),  TrackEndPZ->at(iTrack));

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
			    if ( whichEnd == 0 ) {
				vecX.SetXYZ(TrackStartX->at(iRECOpart) -x,
							TrackStartY->at(iRECOpart) -y,
							TrackStartZ->at(iRECOpart) -z);
			    } else {
				vecX.SetXYZ(TrackEndX->at(iRECOpart)   -x,
							TrackEndY->at(iRECOpart)   -y,
							TrackEndZ->at(iRECOpart)   -z);
			    }
			    Float_t delX = vecX.Cross(MCpartHat).Mag();
                std::cout<<"delX= "<<delX<<" cm"<<std::endl;
                if ( cosTmax <= 0.96 ) continue;
			    //if (  delX   >= 3.0   ) continue;
                

                if((PDG->at(j)==13 && TrackStartQ->at(iRECOpart)==-1) || (PDG->at(j)==-13 && TrackStartQ->at(iRECOpart)==1))
                {
                    hmu_good_recocharge_p_start->Fill(p);
                    
                }
                hmu_all_recocharge_p_start->Fill(p);

               
           }
        }
        

    }
    
    TCanvas *mccanvascharge = new TCanvas("mccanvascharge","",1000,800);
    hmu_good_recocharge_p_start->SetTitle("Muon ratio GAr entrance/LAr exit;p[GeV/c];Q res");
    hmu_good_recocharge_p_start->Divide(hmu_all_recocharge_p_start);
    hmu_good_recocharge_p_start->Draw();
    mccanvascharge->Print("15_ChargeRes.png");

    
}