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
#include <math.h>

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
  float muon_mass = 0.1057;

  TH1F *truemuonpx = new TH1F("truemuonpx",";Muon Px (GeV);events",100,-10,10);
  TH1F *truemuonpy = new TH1F("truemuonpy",";Muon Py (GeV);events",100,-10,10);
  TH1F *truemuonpz = new TH1F("truemuonpz",";Muon Pz (GeV);events",100,0,30);
  TH1F *truemuonp  = new TH1F("truemuonp ",";Muon Ptot (GeV);events",100,0,30);

  TH1F *truemuonpxend = new TH1F("truemuonpxend",";Muon Px end (GeV);events",100,-10,10);
  TH1F *truemuonpyend = new TH1F("truemuonpyend",";Muon Py end (GeV);events",100,-10,10);
  TH1F *truemuonpzend = new TH1F("truemuonpzend",";Muon Pz end (GeV);events",100,0,30);
  TH1F *truemuonpend  = new TH1F("truemuonpend",";Muon Ptot end (GeV);events",100,0,30);

  TH1F *recotrackpx = new TH1F("recotrackpx",";Reco Muon Px (GeV);events",100,-10,10);
  TH1F *recotrackpy = new TH1F("recotrackpy",";Reco Muon Py (GeV);events",100,-10,10);
  TH1F *recotrackpz = new TH1F("recotrackpz",";Reco Muon Pz (GeV);events",100,0,30);
  TH1F *recotrackp  = new TH1F("recotrackp",";Reco Muon Ptot (GeV);events",100,0,30);

  
  TH1F *fracresid  = new TH1F("fracresid","3<P<5 GeV, Standard ECAL;(reco P - true P)/(true P);muons",100,-1,1); 
  TH1F *fracresid2  = new TH1F("fracresid2","3<P<5 GeV, Standard ECAL;#Delta E (GeV);muons",100,0,2);
  TH2F *resY  = new TH2F("resy","3<P<5 GeV, Standard ECAL;y;res",100,-200,0,100,-1,0);
  TH2F *resY2  = new TH2F("resy2","3<P<5 GeV, Standard ECAL;y;#Delta E",100,-200,0,100,0,2);
  /*
  TH1F *fracresid  = new TH1F("fracresid","1<P<3 GeV, Standard ECAL;(reco P - true P)/(true P);muons",100,-1,1);
  TH1F *fracresid2  = new TH1F("fracresid2","1<P<3 GeV, Standard ECAL;(true P end - true P)/(true P);muons",100,0,2);
  TH2F *resY  = new TH2F("resy","1<P<3 GeV, Standard ECAL;y;res",100,-200,0,100,-1,0);
  TH2F *resY2  = new TH2F("resy2","1<P<3 GeV, Standard ECAL;y;#Delta E",100,-200,0,100,0,2);
  */


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
    TVector3 muonpend(0,0,0);
    float muony;

    for (size_t imcp=0; imcp<PDG->size(); ++imcp)
      {
	    if (std::abs(PDG->at(imcp)) == 13)
        {
          muonp.SetX(MCPStartPX->at(imcp));
          muonp.SetY(MCPStartPY->at(imcp));
          muonp.SetZ(MCPStartPZ->at(imcp));
          muonpend.SetX(MCPEndPX->at(imcp));
          muonpend.SetY(MCPEndPY->at(imcp));
          muonpend.SetZ(MCPEndPZ->at(imcp));
          muony = MCPStartY->at(imcp);
          found = true;
          break;
        }
      }
    if (!found) continue;  // no muon found, go to next event
    truemuonpx->Fill(muonp.X());
    truemuonpy->Fill(muonp.Y());
    truemuonpz->Fill(muonp.Z());
    truemuonpxend->Fill(muonpend.X());
    truemuonpyend->Fill(muonpend.Y());
    truemuonpzend->Fill(muonpend.Z());
    float mpvar = muonp.Mag();
    float mpvarend = muonpend.Mag();
    truemuonp->Fill(mpvar);
    truemuonpend->Fill(mpvarend);

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
    float dpp2 = - sqrt(mpvarend*mpvarend+muon_mass*muon_mass) + sqrt(mpvar*mpvar+muon_mass*muon_mass);
    fracresid->Fill(dpp);
    fracresid2->Fill(dpp2);
    resY->Fill(muony,dpp);
    resY2->Fill(muony,dpp2);

  } // end loop over tree entries

  TCanvas *mccanvas1 = new TCanvas("mccanvas1","",1000,800);
  //mccanvas1->Divide(3,3);
  //mccanvas1->cd(1);
  fracresid->Fit("gaus");
  gStyle->SetOptFit(1);
  fracresid->Draw("e0");
  mccanvas1->Print("fracresidV3profhigh.png");

  TCanvas *mccanvas2 = new TCanvas("mccanvas2","",1000,800);
  //mccanvas1->Divide(3,3);
  //mccanvas1->cd(1);
  fracresid2->Fit("landau");
  gStyle->SetOptFit(1);
  fracresid2->Draw("e0");
  mccanvas2->Print("fracresid2V3profhigh.png");

  TCanvas *mccanvas3 = new TCanvas("mccanvas3","",1000,800);
  //mccanvas1->Divide(3,3);
  //mccanvas1->cd(1);
  //fracresid2->Fit("gaus");
  TF1 *f1 = new TF1("f1","[0]+[1]*x",0,1);
  f1->SetParameters(0.,1.);
  f1->SetLineColor(kRed);
  resY->Fit(f1);
  gStyle->SetOptFit(0);
  resY->Draw("col");
  f1->Draw("same");
  mccanvas3->Print("resYprofhigh.png");
  
  TCanvas *mccanvas4 = new TCanvas("mccanvas4","",1000,800);
  //mccanvas1->Divide(3,3);
  //mccanvas1->cd(1);
  //fracresid2->Fit("gaus");
  gStyle->SetOptFit(0);
  resY2->Fit(f1);
  resY2->Draw("col");
  f1->Draw("same");
  mccanvas4->Print("resY2profhigh.png");

}
