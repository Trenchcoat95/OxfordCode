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
    TChain chain("/anatree/GArAnaTree");
    chain.Add("/pnfs/dune/persistent/users/ebrianne/ProductionSamples/ND-LAr/nd_hall_dayone_lar_SPY_v2_wMuID/Anatree/neutrino/*");

    //neutrino.nd_hall_dayone_lar_SPY_v2_wMuID.volArgonCubeActive.Ev973000.Ev973999.2037.anatree.root

    vector<int>     *PDG=0;
    vector<int>     *PDGMother=0;
    vector<int>     *MCTrkID=0;
    vector<int>     *TrajMCPTrackID=0;
    vector<float>   *TrajMCPX =0;
    vector<float>   *TrajMCPY =0;
    vector<float>   *TrajMCPZ =0;

   
    TBranch *b_PDG=0;
    TBranch *b_PDGMother=0;
    TBranch *b_MCTrkID=0;
    TBranch *b_TrajMCPTrackID=0;
    TBranch *b_TrajMCPX=0;
    TBranch *b_TrajMCPY=0;
    TBranch *b_TrajMCPZ=0;


    
    chain.SetBranchAddress("PDG", &PDG, &b_PDG);
    chain.SetBranchAddress("PDGMother", &PDGMother, &b_PDGMother);
    chain.SetBranchAddress("MCTrkID", &MCTrkID, &b_MCTrkID);
    chain.SetBranchAddress("TrajMCPTrackID", &TrajMCPTrackID, &b_TrajMCPTrackID);
    //chain.SetBranchAddress("MCPProc", &MCPProc, &b_MCPProc);
    chain.SetBranchAddress("TrajMCPX", &TrajMCPX, &b_TrajMCPX);
    chain.SetBranchAddress("TrajMCPY", &TrajMCPY, &b_TrajMCPY);
    chain.SetBranchAddress("TrajMCPZ", &TrajMCPZ, &b_TrajMCPZ);



    chain.SetBranchStatus("*",0);
    chain.SetBranchStatus("PDG",1);
    chain.SetBranchStatus("TrajMCPTrackID",1);
    chain.SetBranchStatus("PDGMother",1);
    chain.SetBranchStatus("MCTrkID",1);
    chain.SetBranchStatus("TrajMCPX",1);
    chain.SetBranchStatus("TrajMCPY",1);
    chain.SetBranchStatus("TrajMCPZ",1);
    

    Int_t nentries = (Int_t)chain.GetEntries();  

 
    

    TH2F* hxyLArEnd = new TH2F("hxyLArEnd", "XY LArEnd", 100, -1000, 1000, 100, -1000, 1000);
    TH2F* hxyGArStart = new TH2F("hxyGArStart", "XY GArStart", 100, -1000, 1000, 100, -1000, 1000);
    
    
    


    for (Int_t i=0; i<nentries; i++) 
    {
        chain.GetEntry(i);
        int ID=0;

 

        for (Int_t j=0; j<PDG->size(); j++) 
        {
           if(PDG->at(j)==13 || PDG->at(j)==-13)
           {
               ID=MCTrkID->at(j);
           }
   
        }
        
        int no1=0;
        int no2=0;

        for (Int_t k=0; k<TrajMCPX->size(); k++) 
        {
            if (TrajMCPTrackID->at(k)==ID)
            {
                if(TrajMCPZ->at(k)<950 && TrajMCPZ->at(k)>850  && no1==0)
                {
                   hxyLArEnd->Fill(TrajMCPX->at(k),TrajMCPY->at(k));
                   no1++;
                }

                if(TrajMCPZ->at(k)<1200 && TrajMCPZ->at(k)>1100  && no2==0)
                {
                   hxyGArStart->Fill(TrajMCPX->at(k),TrajMCPY->at(k));
                   no2++;
                }
            }
        }
        

    }

    
    TH2F *hxyLArGArRatio = (TH2F*) hxyGArStart->Clone("hxyLArGArRatio");
    hxyLArGArRatio->Divide(hxyLArEnd);

    TCanvas *mccanvasRatio = new TCanvas("mccanvasRatio","",1000,800);
    gStyle->SetOptStat(0);
    hxyLArGArRatio->SetTitle("Muon ratio GAr entrance/LAr exit;x[cm];y[cm]");
    hxyLArGArRatio->Draw("COLZ");
    mccanvasRatio->Print("14_LArGArRatio.png");
    

    
 
}