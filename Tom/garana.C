#define garana_cxx
#include "garana.h"
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

void garana::Loop()
{
  //   In a ROOT session, you can do:
  //      root> .L garana.C
  //      root> garana t
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.Loop();       // Loop on all entries
  //

  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;

  gStyle->SetOptFit(1);

  //TVector3 detcent(0,218.196,1124.02);
  TVector3 detcent(0,0,0);
  TVector3 xhat(1,0,0);
  float detR = 250;

  TH1F *truemuonpx = new TH1F("truemuonpx",";Muon Px (GeV);events",100,-10,10);
  TH1F *truemuonpy = new TH1F("truemuonpy",";Muon Py (GeV);events",100,-10,10);
  TH1F *truemuonpz = new TH1F("truemuonpz",";Muon Pz (GeV);events",100,0,30);
  TH1F *truemuonp  = new TH1F("truemuonp ",";Muon Ptot (GeV);events",100,0,30);

  TH1F *recotrackpx = new TH1F("recotrackpx",";Reco Muon Px (GeV);events",100,-10,10);
  TH1F *recotrackpy = new TH1F("recotrackpy",";Reco Muon Py (GeV);events",100,-10,10);
  TH1F *recotrackpz = new TH1F("recotrackpz",";Reco Muon Pz (GeV);events",100,0,30);
  TH1F *recotrackp  = new TH1F("recotrackp",";Reco Muon Ptot (GeV);events",100,0,30);

  //TH1F *fracresid  = new TH1F("fracresid","1<P<3 GeV, No ECAL;(reco P - true P)/(true P);muons",100,-1,1);
  TH1F *fracresid  = new TH1F("fracresid","3<P<5 GeV, No ECAL;(reco P - true P)/(true P);muons",100,-1,1);


  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    // find the first muon in the sample
    bool found=false;
    TVector3 muonp(0,0,0);
    for (size_t imcp=0; imcp<PDG->size(); ++imcp)
      {
	if (std::abs(PDG->at(imcp)) == 13)
	  {
	    muonp.SetX(MCPStartPX->at(imcp));
	    muonp.SetY(MCPStartPY->at(imcp));
	    muonp.SetZ(MCPStartPZ->at(imcp));
	    found = true;
	    break;
	  }
      }
    if (!found) continue;  // no muon found, go to next event
    truemuonpx->Fill(muonp.X());
    truemuonpy->Fill(muonp.Y());
    truemuonpz->Fill(muonp.Z());
    float mpvar = muonp.Mag();
    truemuonp->Fill(mpvar);

    // look for the longest track and pick the direction that points along +z
    found = false;
    TVector3 trackp(0,0,0);
    float lenmax = 0;
    for (size_t itrack=0; itrack<TrackStartX->size(); ++itrack)
      {
	float tlen = (TrackLenF->at(itrack) + TrackLenB->at(itrack))/2.0;
	if (tlen > lenmax)
	  {
	    lenmax = tlen;
	    found = true;
	    if (TrackStartPZ->at(itrack)>0)
	      {
		trackp.SetX(TrackStartPX->at(itrack));
		trackp.SetY(TrackStartPY->at(itrack));
		trackp.SetZ(TrackStartPZ->at(itrack));
	      }
	    else
	      {
		trackp.SetX(TrackEndPX->at(itrack));
		trackp.SetY(TrackEndPY->at(itrack));
		trackp.SetZ(TrackEndPZ->at(itrack));
	      }
	  }
      }

    if (!found) continue;  // no muon found, go to next event
    recotrackpx->Fill(trackp.X());
    recotrackpy->Fill(trackp.Y());
    recotrackpz->Fill(trackp.Z());
    float tpvar = trackp.Mag();
    recotrackp->Fill(tpvar);

    float dpp = (tpvar - mpvar)/mpvar;
    fracresid->Fill(dpp);

  } // end loop over tree entries

  TCanvas *mccanvas1 = new TCanvas("mccanvas1","",1000,800);
  //mccanvas1->Divide(3,3);
  //mccanvas1->cd(1);
  fracresid->Fit("gaus");
  gStyle->SetOptFit(1);
  fracresid->Draw("e0");
  mccanvas1->Print("fracresid.png");

}
