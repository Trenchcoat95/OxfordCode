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
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>

void Snippet_simple()
{
    gStyle->SetOptFit(1);
   auto c0 = new TCanvas("c0","A Simple Graph with error bars",200,10,700,500);
   const Int_t n = 4;
   Double_t x[n]  = {2,3,4,5};
   Double_t y[n]  = {3.18,3.93,4.7,5.16};
   Double_t ex[n] = {0,0,0,0};
   Double_t ey[n] = {0.09,0.5,0.2,0.16};
   auto gr = new TGraphErrors(n,x,y,ex,ey);
   gr->SetTitle("Evolution of p_{0};n_{seed};p_{0}");
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->Draw("AP");
   TF1 *fa2 = new TF1("fa2","x+[c]*(x)+[b]",0,100);
   gr->Fit(fa2);
   gStyle->SetStatX(0.5);
   gStyle->SetStatY(0.9);
   c0->Print("./Plots/log_Evolution_p0.png");

   auto c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
   const Int_t n1 = 4;
   Double_t x1[n1]  = {2,3,4,5};
   Double_t y1[n1]  = {-13.2,-14.6,-16.4,-17.9};
   Double_t ex1[n1] = {0,0,0,0};
   Double_t ey1[n1] = {0.2,0.4,0.6,0.5};
   auto gr1 = new TGraphErrors(n1,x1,y1,ex1,ey1);
   gr1->SetTitle("Evolution of p_{1};n_{seed};p_{1}");
   gr1->SetMarkerColor(4);
   gr1->SetMarkerStyle(21);
   TF1 *fa1 = new TF1("fa1","[a]*x+[b]",0,100);
   gr1->Fit(fa1);
   gStyle->SetStatX(0.85);
   gStyle->SetStatY(0.9);
   gr1->Draw("AP");
   c1->Print("./Plots/log_Evolution_p1.png");

   gStyle->SetOptFit(1);
   auto c02 = new TCanvas("c02","A Simple Graph with error bars",200,10,700,500);
   const Int_t n2 = 4;
   Double_t x2[n2]  = {2,3,4,5};
   Double_t y2[n2]  = {1.18,0.93,0.7,0.16};
   Double_t ex2[n2] = {0,0,0,0};
   Double_t ey2[n2] = {0.09,0.5,0.2,0.16};
   auto gr2 = new TGraphErrors(n2,x2,y2,ex2,ey2);
   gr2->SetTitle("Evolution of p_{0};n_{seed};p_{0}");
   gr2->SetMarkerColor(4);
   gr2->SetMarkerStyle(21);
   gr2->Draw("AP");
   TF1 *fa22 = new TF1("fa22","[a]*x+[b]",0,100);
   //fa22->SetParameters(-2,1);
   gr2->Fit(fa22);
   gStyle->SetStatX(0.5);
   gStyle->SetStatY(0.9);
   c02->Print("./Plots/log_Testb.png");

   auto c03 = new TCanvas("c03","A Simple Graph with error bars",200,10,700,500);
   const Int_t n3 = 3;
   Double_t x3[n3]  = {1.31414,6.65808,0.863751};
   Double_t y3[n3]  = {1.3475,6.67582,0.87728};
   auto gr3 = new TGraph(n3,x3,y3);
   gr3->SetTitle("Total time evolution;#sumt_{bunch}[s];t[s]");
   gr3->SetMarkerColor(4);
   gr3->SetMarkerStyle(21);
   gr3->Draw("AP");
   TF1 *fa33 = new TF1("fa33","[a]*x+[b]",0,100);
   //fa22->SetParameters(-2,1);
   gr3->Fit(fa33);
   gStyle->SetStatX(0.5);
   gStyle->SetStatY(0.9);
   c03->Print("./Plots/tVSsumtbunch.png");


}