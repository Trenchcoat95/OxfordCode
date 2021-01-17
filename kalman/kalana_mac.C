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
    std::string folder= "helix_rndx_sm/no_sm/";
    std::string  sm= "x01z3";
    std::string  R= "3_2";
    std::string  filename = "m_perfect_helix_rndx_sm"+sm+"_R_"+R+"_stdK.root";
    kalana kal = kalana("t1_FWD;1",filename.c_str());
    kal.Loop(folder,sm,R);
    //garana g;

    const double B = 0.4;
    
    TTree* t= kal.t;
    float xpost=0;
    TVectorF *parvect=0;
    TMatrixF *Pt=0;
    float xht=0;
    float yht=0;
    float zht=0;
    TBranch *b_xpost=0;
    TBranch *b_parvect=0;
    TBranch *b_Pt=0;
    TBranch *b_xht=0;
    TBranch *b_yht=0;
    TBranch *b_zht=0;
    t->SetBranchAddress("xht",&xht,&b_xht);
    t->SetBranchAddress("yht",&yht,&b_yht);
    t->SetBranchAddress("zht",&zht,&b_zht);
    t->SetBranchAddress("xpost",&xpost,&b_xpost);
    t->SetBranchAddress("parvect",&parvect,&b_parvect);
    t->SetBranchAddress("Pt",&Pt,&b_Pt);

    /*
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
    */ 

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
      pTpar->SetPoint(i-1,xpost,abs(pT));
      //double x,y;
      ppar->SetPoint(i-1,xpost,sqrt(pT*pT+px*px));
      
      
      /////Truth for Toy Montecarlo
      YTrue->SetPoint(i-1,xht,yht);
      ZTrue->SetPoint(i-1,xht,zht);
      double pTT=0.01*0.3*B/(-0.014);
      double pxTrue=pTT*tan(-0.05);
      pTTrue->SetPoint(i-1,xht,abs(pTT));
      pTrue->SetPoint(i-1,xht,sqrt(pTT*pTT+pxTrue*pxTrue));
      
      
    }

    /*
    ///Truth for proper Montecarlo
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
    */
    
    TCanvas *mccanvaspT = new TCanvas("mccanvaspT","",1000,800);
    pTpar->SetTitle("pT;x(cm);p(GeV/c)");
    pTpar->SetLineColor(kGreen);
    pTpar->Draw();
    pTTrue->SetLineColor(kRed);
    pTTrue->SetMarkerColor(kRed);
    pTTrue->Draw("same");
    auto legendpT = new TLegend(0.55,0.75,0.9,0.9);
    legendpT->AddEntry(pTpar,"Kalman Fitter best estimate","lep");
    legendpT->AddEntry(pTTrue,"Montecarlo Truth","lep");
    legendpT->Draw();
    std::string st = folder + "helix_pTTrue_rndx_sm" + sm + "_R_" + R + "_stdK.png";
    mccanvaspT->Print(st.c_str());
    
    TCanvas *mccanvasp = new TCanvas("mccanvasp","",1000,800);
    ppar->SetTitle("p;x(cm);p(GeV/c)");
    ppar->SetLineColor(kGreen);
    ppar->Draw();
    pTrue->SetLineColor(kRed);
    pTrue->SetMarkerColor(kRed);
    pTrue->Draw("same");
    auto legendp = new TLegend(0.55,0.75,0.9,0.9);
    legendp->AddEntry(ppar,"Kalman Fitter best estimate","lep");
    legendp->AddEntry(pTrue,"Montecarlo Truth","lep");
    legendp->Draw();
    st = folder + "helix_pTrue_rndx_sm" + sm + "_R_" + R + "_stdK.png";
    mccanvasp->Print(st.c_str());
    
    TCanvas *mccanvasy = new TCanvas("mccanvasy","",1000,800);
    YTrue->SetTitle("Y;x(cm);y(cm)");
    YTrue->SetLineColor(kRed);
    YTrue->SetMarkerColor(kRed);
    YTrue->Draw();
    Ypar->SetLineColor(kGreen);
    Ypar->Draw("same");
    auto legendy = new TLegend(0.55,0.75,0.9,0.9);
    legendy->AddEntry(Ypar,"Kalman Fitter best estimate","lep");
    legendy->AddEntry(YTrue,"Montecarlo Truth","lep");
    legendy->Draw();
    st = folder + "helix_YTrue_rndx_sm" + sm + "_R_" + R + "_stdK.png";
    mccanvasy->Print(st.c_str());

    TCanvas *mccanvasz = new TCanvas("mccanvasz","",1000,800);
    ZTrue->SetTitle("Z;x(cm);z(cm)");
    ZTrue->SetLineColor(kRed);
    ZTrue->SetMarkerColor(kRed);
    ZTrue->Draw();
    Zpar->SetLineColor(kGreen);
    Zpar->Draw("same");
    auto legendz = new TLegend(0.325,0.75,0.675,0.9);
    legendz->AddEntry(Zpar,"Kalman Fitter best estimate","lep");
    legendz->AddEntry(ZTrue,"Montecarlo Truth","lep");
    legendz->Draw();
    st = folder + "helix_ZTrue_rndx_sm" + sm + "_R_" + R + "_stdK.png";
    mccanvasz->Print(st.c_str());
    
}

