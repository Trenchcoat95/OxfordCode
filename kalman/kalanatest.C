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


void kalanatest()
{
    kalana kal = kalana("t1_FWD;1","demo.root");
    kal.Loop();
    garana g;
    
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
    TBranch *b_TrajMCPX=0;
    TBranch *b_TrajMCPY=0;
    TBranch *b_TrajMCPZ=0;
    gt->SetBranchAddress("TrajMCPX", &TrajMCPX, &b_TrajMCPX);
    gt->SetBranchAddress("TrajMCPY", &TrajMCPY, &b_TrajMCPY);
    gt->SetBranchAddress("TrajMCPZ", &TrajMCPZ, &b_TrajMCPZ);
    
    Int_t nentries = (Int_t)t->GetEntries();
    TGraphErrors *Ypar =new  TGraphErrors(nentries-1);
    TGraphErrors *Zpar =new  TGraphErrors(nentries-1);
    TGraph *YTrue =new  TGraphErrors(nentries-1);
    TGraph *ZTrue =new  TGraphErrors(nentries-1);

    for (Int_t i=1; i<nentries; i++) {
      t->GetEntry(i);
      Ypar->SetPoint(i-1,t->GetLeaf("xpost")->GetValue(),(*parvect)[0]);
      Ypar->SetPointError(i-1,xpost,TMath::Sqrt((*Pt)[0][0]));
    }
    
    gt->GetEntry(0);
    std::cout<<TrajMCPX->size()<<std::endl;
    for (Int_t i=1; i<TrajMCPX->size(); i++) {
      std::cout<<(*TrajMCPX)[i]<<std::endl;
    }
   //float xpost;
   //t->SetBranchAddress("xht",&xht);
   //t->GetEntry(2);
   //std::cout<<xht<<std::endl;

}

