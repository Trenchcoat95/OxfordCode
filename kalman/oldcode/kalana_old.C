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


void kalana()
{

TFile *f = new TFile("demo.root");
TTree *t1 = (TTree*)f->Get("t1_FWD;1");

float xht = 0;
float yht = 0;
float zht = 0;
float xpost = 0;
int ev = 0;
TVectorF *parvect = 0;
TVectorF *predstept = 0;
TMatrixF *Pt = 0;
TMatrixF *PPredt = 0;
TMatrixF *Rt = 0;

TBranch *b_xht = 0;
TBranch *b_yht = 0;
TBranch *b_zht = 0;
TBranch *b_xpost = 0;
TBranch *b_ev = 0;
TBranch *b_parvect = 0;
TBranch *b_predstept = 0;
TBranch *b_Pt = 0;
TBranch *b_PPredt = 0;
TBranch *b_Rt = 0;

t1->SetBranchAddress("xht",&xht,&b_xht);
t1->SetBranchAddress("yht",&yht,&b_yht);
t1->SetBranchAddress("zht",&zht,&b_zht);
t1->SetBranchAddress("ev",&ev,&b_ev);
t1->SetBranchAddress("xpost",&xpost,&b_xpost);
t1->SetBranchAddress("parvect",&parvect,&b_parvect);
t1->SetBranchAddress("predstept",&predstept,&b_predstept);
t1->SetBranchAddress("Pt",&Pt,&b_Pt);
t1->SetBranchAddress("PPredt",&PPredt,&b_PPredt);
t1->SetBranchAddress("Rt",&Rt,&b_Rt);

Int_t nentries = (Int_t)t1->GetEntries();

TGraphErrors *Yh =new  TGraphErrors(nentries-1);
TGraphErrors *Ypred =new  TGraphErrors(nentries-1);
TGraphErrors *Ypar =new  TGraphErrors(nentries-1);

TGraphErrors *Zh =new  TGraphErrors(nentries-1);
TGraphErrors *Zpred =new  TGraphErrors(nentries-1);
TGraphErrors *Zpar =new  TGraphErrors(nentries-1);

TGraphErrors *Curvpred =new  TGraphErrors(nentries-1);
TGraphErrors *Curvpar =new  TGraphErrors(nentries-1);

TGraphErrors *Phipred =new  TGraphErrors(nentries-1);
TGraphErrors *Phipar =new  TGraphErrors(nentries-1);

TGraphErrors *Lambdapred =new  TGraphErrors(nentries-1);
TGraphErrors *Lambdapar =new  TGraphErrors(nentries-1);


for (Int_t i=1; i<nentries; i++) {

      t1->GetEntry(i);

      Yh->SetPoint(i-1,xht,yht);
      Yh->SetPointError(i-1,xht,TMath::Sqrt((*Rt)[0][0]));
      Ypred->SetPoint(i-1,xpost,(*predstept)[0]);
      Ypred->SetPointError(i-1,xpost,TMath::Sqrt((*PPredt)[0][0]));
      Ypar->SetPoint(i-1,xpost,(*parvect)[0]);
      Ypar->SetPointError(i-1,xpost,TMath::Sqrt((*Pt)[0][0]));

      Zh->SetPoint(i-1,xht,zht);
      Zh->SetPointError(i-1,0,TMath::Sqrt((*Rt)[1][1]));
      Zpred->SetPoint(i-1,xpost,(*predstept)[1]);
      Zpred->SetPointError(i-1,0,TMath::Sqrt((*PPredt)[1][1]));
      Zpar->SetPoint(i-1,xpost,(*parvect)[1]);
      Zpar->SetPointError(i-1,0,TMath::Sqrt((*Pt)[1][1]));

      Curvpred->SetPoint(i-1,xpost,(*predstept)[2]);
      Curvpred->SetPointError(i-1,0,TMath::Sqrt((*PPredt)[2][2]));
      Curvpar->SetPoint(i-1,xpost,(*parvect)[2]);
      Curvpar->SetPointError(i-1,0,TMath::Sqrt((*Pt)[2][2]));

      Phipred->SetPoint(i-1,xpost,(*predstept)[3]);
      Phipred->SetPointError(i-1,0,TMath::Sqrt((*PPredt)[3][3]));
      Phipar->SetPoint(i-1,xpost,(*parvect)[3]);
      Phipar->SetPointError(i-1,0,TMath::Sqrt((*Pt)[3][3]));

      Lambdapred->SetPoint(i-1,xpost,(*predstept)[4]);
      Lambdapred->SetPointError(i-1,0,TMath::Sqrt((*PPredt)[4][4]));
      Lambdapar->SetPoint(i-1,xpost,(*parvect)[4]);
      Lambdapar->SetPointError(i-1,0,TMath::Sqrt((*Pt)[4][4]));
      }

TCanvas *mccanvas1 = new TCanvas("mccanvas1","",1000,800);
Yh->SetTitle("Y");
Yh->SetLineColor(kRed);
Yh->Draw();
Ypred->SetLineColor(kBlue);
Ypred->Draw("same");
Ypar->SetLineColor(kGreen);
Ypar->Draw("same");
auto legend = new TLegend(0.1,0.75,0.45,0.9);
legend->AddEntry(Yh,"TPC Cluster measured value","lep");
legend->AddEntry(Ypred,"A priori prediction","lep");
legend->AddEntry(Ypar,"A posteriori prediction","lep");
legend->Draw();
mccanvas1->Print("Y.png");



TCanvas *mccanvas2 = new TCanvas("mccanvas2","",1000,800);
Zh->SetTitle("Z");
Zh->SetLineColor(kRed);
Zh->Draw();
Zpred->SetLineColor(kBlue);
Zpred->Draw("same");
Zpar->SetLineColor(kGreen);
Zpar->Draw("same");
auto legend2 = new TLegend(0.55,0.75,0.9,0.9);
legend2->AddEntry(Zh,"TPC Cluster measured value","lep");
legend2->AddEntry(Zpred,"A priori prediction","lep");
legend2->AddEntry(Zpar,"A posteriori prediction","lep");
legend2->Draw();
mccanvas2->Print("Z.png");



TCanvas *mccanvas3 = new TCanvas("mccanvas3","",1000,800);
Curvpred->SetTitle("Curv");
Curvpred->SetLineColor(kBlue);
Curvpred->Draw();
Curvpar->SetLineColor(kGreen);
Curvpar->Draw("same");
auto legend3 = new TLegend(0.55,0.75,0.9,0.9);
legend3->AddEntry(Curvpred,"A priori prediction","lep");
legend3->AddEntry(Curvpar,"A posteriori prediction","lep");
legend3->Draw();
mccanvas3->Print("Curv.png");

TCanvas *mccanvas4 = new TCanvas("mccanvas4","",1000,800);
Phipred->SetTitle("Phi");
Phipred->SetLineColor(kBlue);
Phipred->Draw();
Phipar->SetLineColor(kGreen);
Phipar->Draw("same");
auto legend4 = new TLegend(0.55,0.75,0.9,0.9);
legend4->AddEntry(Phipred,"A priori prediction","lep");
legend4->AddEntry(Phipar,"A posteriori prediction","lep");
legend4->Draw();
mccanvas4->Print("Phi.png");


TCanvas *mccanvas5 = new TCanvas("mccanvas5","",1000,800);
Lambdapred->SetTitle("Lambda");
Lambdapred->SetLineColor(kBlue);
Lambdapred->Draw();
Lambdapar->SetLineColor(kGreen);
Lambdapar->Draw("same");
auto legend5 = new TLegend(0.55,0.75,0.9,0.9);
legend5->AddEntry(Lambdapred,"A priori prediction","lep");
legend5->AddEntry(Lambdapar,"A posteriori prediction","lep");
legend5->Draw();
mccanvas2->Print("Lambda.png");






}

