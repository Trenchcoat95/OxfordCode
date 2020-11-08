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
#include "kalana.h"
#include "kalLoop.h"
#include "garana.h"


void kalana_mac()
{
    kalana kal = kalana("t1_FWD;1","smeared_helix.root");
    kal.Loop();
    garana g;

    const double B = 0.4;
    
    TTree* t= kal.t;
    float xpost=0;
    TVectorF *parvect=0;
    TMatrixF *Pt=0;
    TBranch *b_xpost=0;
    TBranch *b_parvect=0;
    TBranch *b_Pt=0;
    t->SetBranchAddress("xpost",&xpost,&b_xpost);
    t->SetBranchAddress("parvect",&parvect,&b_parvect);
    t->SetBranchAddress("Pt",&Pt,&b_Pt);

    TTree* gt=g.fChain;
    vector<float>   *TrajMCPX=0;
    vector<float>   *TrajMCPY=0;
    vector<float>   *TrajMCPZ=0;
    vector<float>   *TrajMCPPX=0;
    vector<float>   *TrajMCPPY=0;
    vector<float>   *TrajMCPPZ=0;
    vector<float>   *TrajMCPIndex=0;
    TBranch *b_TrajMCPX=0;
    TBranch *b_TrajMCPY=0;
    TBranch *b_TrajMCPZ=0;
    TBranch *b_TrajMCPPX=0;
    TBranch *b_TrajMCPPY=0;
    TBranch *b_TrajMCPPZ=0;
    TBranch *b_TrajMCPIndex=0;
    gt->SetBranchAddress("TrajMCPX", &TrajMCPX, &b_TrajMCPX);
    gt->SetBranchAddress("TrajMCPY", &TrajMCPY, &b_TrajMCPY);
    gt->SetBranchAddress("TrajMCPZ", &TrajMCPZ, &b_TrajMCPZ);
    gt->SetBranchAddress("TrajMCPPX", &TrajMCPPX, &b_TrajMCPPX);
    gt->SetBranchAddress("TrajMCPPY", &TrajMCPPY, &b_TrajMCPPY);
    gt->SetBranchAddress("TrajMCPPZ", &TrajMCPPZ, &b_TrajMCPPZ);
    gt->SetBranchAddress("TrajMCPIndex", &TrajMCPIndex, &b_TrajMCPIndex);
    
    Int_t nentries = (Int_t)t->GetEntries();
    TGraphErrors *Ypar =new  TGraphErrors(nentries-1);
    TGraphErrors *Zpar =new  TGraphErrors(nentries-1);
    TGraphErrors *ppar =new  TGraphErrors(nentries-1);
    TGraph *pTpar =new  TGraph(nentries-1);
    TGraph *YTrue =new  TGraph(nentries-1);
    TGraph *ZTrue =new  TGraph(nentries-1);
    TGraph *pTrue =new  TGraph(nentries-1);
    TGraph *pTTrue =new  TGraph(nentries-1);

    for (Int_t i=1; i<nentries; i++) {
      t->GetEntry(i);
      Ypar->SetPoint(i-1,xpost,(*parvect)[0]);
      Ypar->SetPointError(i-1,xpost,TMath::Sqrt((*Pt)[0][0]));
      Zpar->SetPoint(i-1,xpost,(*parvect)[1]);
      Zpar->SetPointError(i-1,xpost,TMath::Sqrt((*Pt)[1][1]));
      double pT=0.01*0.3*B/((*parvect)[2]);
      double px=pT*tan((*parvect)[4]);
      //std::cout<<pT<<std::endl;
      pTpar->SetPoint(i-1,xpost,abs(pT));
      double x,y;
      //pTpar->GetPoint(i-1,x,y);
      //std::cout<<y<<std::endl;
      ppar->SetPoint(i-1,xpost,sqrt(pT*pT+px*px));
    }
    //t->GetLeaf("xpost")->GetValue()
    gt->GetEntry(0);
    //std::cout<<TrajMCPX->size()<<std::endl;
    for (Int_t i=1; i<TrajMCPX->size(); i++) 
    {
        if(TrajMCPIndex->at(i)==0)
        {
            YTrue->SetPoint(i-1,TrajMCPX->at(i),TrajMCPY->at(i));
            ZTrue->SetPoint(i-1,TrajMCPX->at(i),TrajMCPZ->at(i));
            pTTrue->SetPoint(i-1,TrajMCPX->at(i),sqrt(TrajMCPPZ->at(i)*TrajMCPPZ->at(i)+TrajMCPPY->at(i)*TrajMCPPY->at(i)));
            pTrue->SetPoint(i-1,TrajMCPX->at(i),sqrt(TrajMCPPZ->at(i)*TrajMCPPZ->at(i)+TrajMCPPX->at(i)*TrajMCPPX->at(i)+TrajMCPPY->at(i)*TrajMCPPY->at(i)));
        }
    }
   
    /*
    TCanvas *mccanvaspT = new TCanvas("mccanvaspT","",1000,800);
    pTpar->SetTitle("pT;x(cm);p(GeV/c)");
    pTpar->SetLineColor(kGreen);
    pTpar->Draw();
    pTTrue->SetLineColor(kRed);
    pTTrue->SetMarkerColor(kRed);
    pTTrue->Draw("same");
    auto legendpT = new TLegend(0.1,0.75,0.45,0.9);
    legendpT->AddEntry(pTpar,"Kalman Fitter best estimate","lep");
    legendpT->AddEntry(pTTrue,"Montecarlo Truth","lep");
    legendpT->Draw();
    mccanvaspT->Print("helix_pTTrue.png");
    
    TCanvas *mccanvasp = new TCanvas("mccanvasp","",1000,800);
    ppar->SetTitle("p;x(cm);p(GeV/c)");
    ppar->SetLineColor(kGreen);
    ppar->Draw();
    pTrue->SetLineColor(kRed);
    pTrue->SetMarkerColor(kRed);
    pTrue->Draw("same");
    auto legendp = new TLegend(0.1,0.75,0.45,0.9);
    legendp->AddEntry(ppar,"Kalman Fitter best estimate","lep");
    legendp->AddEntry(pTrue,"Montecarlo Truth","lep");
    legendp->Draw();
    mccanvasp->Print("helix_pTrue.png");
    
    TCanvas *mccanvasy = new TCanvas("mccanvasy","",1000,800);
    YTrue->SetTitle("Y;x(cm);y(cm)");
    YTrue->SetLineColor(kRed);
    YTrue->SetMarkerColor(kRed);
    YTrue->Draw();
    Ypar->SetLineColor(kGreen);
    Ypar->Draw("same");
    auto legendy = new TLegend(0.1,0.75,0.45,0.9);
    legendy->AddEntry(Ypar,"Kalman Fitter best estimate","lep");
    legendy->AddEntry(YTrue,"Montecarlo Truth","lep");
    legendy->Draw();
    mccanvasy->Print("helix_YTrue.png");

    TCanvas *mccanvasz = new TCanvas("mccanvasz","",1000,800);
    ZTrue->SetTitle("Z;x(cm);z(cm)");
    ZTrue->SetLineColor(kRed);
    ZTrue->SetMarkerColor(kRed);
    ZTrue->Draw();
    Zpar->SetLineColor(kGreen);
    Zpar->Draw("same");
    auto legendz = new TLegend(0.55,0.75,0.9,0.9);
    legendz->AddEntry(Zpar,"Kalman Fitter best estimate","lep");
    legendz->AddEntry(ZTrue,"Montecarlo Truth","lep");
    legendz->Draw();
    mccanvasz->Print("helix_ZTrue.png");
    */
}

