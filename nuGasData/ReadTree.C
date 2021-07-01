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


void ReadTree()
{
    TChain* t = new TChain("tree");
    t->Add("outAna10_DUNEGiBUU__EXCL3a10nuC_Filelist_DUNE_nu_T0_Carbon_800.root");

    Double_t muonmomentum=0;
    Double_t pionmomentum=0;
    Double_t protonmomentum=0;
    Double_t muontheta=0;
    Double_t piontheta=0;
    Double_t protontheta=0;
    Double_t perweight=0;

    TBranch* b_muonmomentum=0;
    TBranch* b_pionmomentum=0;
    TBranch* b_protonmomentum=0;
    TBranch* b_muontheta=0;
    TBranch* b_piontheta=0;
    TBranch* b_protontheta=0;
    TBranch* b_perweight=0;

    t->SetBranchAddress("muonmomentum", &muonmomentum, &b_muonmomentum);
    t->SetBranchAddress("pionmomentum", &pionmomentum, &b_pionmomentum);
    t->SetBranchAddress("protonmomentum", &protonmomentum, &b_protonmomentum);
    t->SetBranchAddress("muontheta", &muontheta, &b_muontheta);
    t->SetBranchAddress("piontheta", &piontheta, &b_piontheta);
    t->SetBranchAddress("protontheta", &protontheta, &b_protontheta);
    t->SetBranchAddress("perweight", &perweight, &b_perweight);

    TH1D* hmuonmomentum = new TH1D("hmuonmomentum", "muon momentum", 500, 0.02, 10);
    TH1D* hpionmomentum = new TH1D("hpionmomentum", "pion momentum", 500, 0.02, 10);
    TH1D* hprotonmomentum = new TH1D("hprotonmomentum", "proton momentum", 500, 0.02, 10);
    TH1D* hmuontheta = new TH1D("hmuontheta", "muon theta", 180, 0, 180);
    TH1D* hpiontheta = new TH1D("hpiontheta", "pion theta", 180, 0, 180);
    TH1D* hprotontheta = new TH1D("hprotontheta", "proton theta", 180, 0, 180);

    
    Int_t nentries = t->GetEntries();

    bool showprog = true;  
    if(showprog==true) std::cout<<"Progress:  "<<std::endl;

    for (Int_t i=0; i<nentries; i++) 
    {
        t->GetEntry(i);

        int prog = 100*i/nentries;
        std::string strprog = std::to_string(prog);
        if(showprog==true) std::cout<<strprog<<"%";

        hmuonmomentum->Fill(muonmomentum,perweight);
        hpionmomentum->Fill(pionmomentum,perweight);
        hprotonmomentum->Fill(protonmomentum,perweight);
        hmuontheta->Fill(muontheta,perweight);
        hpiontheta->Fill(piontheta,perweight);
        hprotontheta->Fill(protontheta,perweight);

        if(showprog==true) std::cout << std::string(strprog.length(),'\b')<<"\b";
    }

    
    TCanvas *mccanvasmomentum = new TCanvas("mccanvasmomentum","",1000,800);
    gPad->SetLogx();
    hpionmomentum->SetTitle("momentum;p(GeV/c);d#sigma");
    hmuonmomentum->SetMarkerColor(kBlue);
    hpionmomentum->SetMarkerColor(kGreen);
    hprotonmomentum->SetMarkerColor(kRed);
    hmuonmomentum->SetLineColor(kBlue);
    hpionmomentum->SetLineColor(kGreen);
    hprotonmomentum->SetLineColor(kRed);
    hmuonmomentum->SetLineWidth(2);
    hpionmomentum->SetLineWidth(2);
    hprotonmomentum->SetLineWidth(2);
    hpionmomentum->Draw("HIST C");
    hmuonmomentum->Draw("HIST C SAME");
    hprotonmomentum->Draw("HIST C SAME");
    mccanvasmomentum->Print("momentum.png");
    
    TCanvas *mccanvastheta = new TCanvas("mccanvastheta","",1000,800);
    gPad->SetLogx(0);
    hmuontheta->SetTitle("#theta;#theta;d#sigma");
    hmuontheta->SetMarkerColor(kBlue);
    hpiontheta->SetMarkerColor(kGreen);
    hprotontheta->SetMarkerColor(kRed);
    hmuontheta->SetLineColor(kBlue);
    hpiontheta->SetLineColor(kGreen);
    hprotontheta->SetLineColor(kRed);
    hmuontheta->SetLineWidth(2);
    hpiontheta->SetLineWidth(2);
    hprotontheta->SetLineWidth(2);
    hmuontheta->Draw("HIST C");
    hpiontheta->Draw("HIST C SAME");
    hprotontheta->Draw("HIST C SAME");
    mccanvastheta->Print("theta.png");
}