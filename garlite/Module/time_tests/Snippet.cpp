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

void Snippet(int n=3, int logx=0, int logy=0)
{
//TGraph *MyGraph = new TGraph("time.txt");
TString h;
if(n==1) h="time5cm.txt";
if(n==2) h="timeold.txt";
if(n==3) h="time.txt";
if(n==4) h="time_quadruplet.txt";
if(n==5) h="time_quintuplet.txt";

TGraph *MyGraph = new TGraph(h);


TF1 *fa2 = new TF1("fa2","abs([0])*x**[1]+(abs([0])*(x**[1]))**2",0,100);

if(n==5) fa2->SetParameters(1e-9,5); //quintuplets
if(n==4) fa2->SetParameters(1e-8,5); //quadruplets
if(n<=3) fa2->SetParameters(1e-6,3); //triplets

MyGraph->Fit("fa2","","",0,100);

gStyle->SetOptFit(1);

if(n==1) MyGraph->SetTitle("Process Time (Triplet 5cm);nTPCpoints;time [s]");
if(n==2) MyGraph->SetTitle("Process Time (Triplet Old);nTPCpoints;time [s]");
if(n==3) MyGraph->SetTitle("Process Time (Triplet);nTPCpoints;time [s]");
if(n==4) MyGraph->SetTitle("Process Time (Quadruplet);nTPCpoints;time [s]");
if(n==5) MyGraph->SetTitle("Process Time (Quintuplet);nTPCpoints;time [s]");

gStyle->SetStatX(0.5);
gStyle->SetStatY(0.9);

TCanvas *c2 = new TCanvas("c2","ac2", 200,10,700,500);

MyGraph->Draw("A*");

if (logy>0) gPad->SetLogy(1);
if (logx>0) gPad->SetLogx(1);

if (logy==0) gPad->SetLogy(0);
if (logx==0) gPad->SetLogx(0);

TString save;

if (n==1) save="./Plots/pol_timetriplet5cm";
if (n==2) save="./Plots/pol_timetripletold";
if (n==3) save="./Plots/pol_timetriplet";
if (n==4) save="./Plots/pol_timequadruplet";
if (n==5) save="./Plots/pol_timequintuplet";

if (logy>0) save+="logy";
if (logx>0) save+="logx";

save+=".png";

c2->Print(save);

}