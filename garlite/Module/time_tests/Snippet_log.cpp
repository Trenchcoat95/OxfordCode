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

void Snippet_log(int n=3)
{
//TGraph *MyGraph = new TGraph("time.txt");
TString h;
if(n==1) h="time5cm.txt";
if(n==2) h="timeold.txt";
if(n==3) h="time.txt";
if(n==4) h="time_quadruplet.txt";
if(n==5) h="time_quintuplet.txt";
if(n==6) h="time_doublet.txt";
if(n==7) h="time_timedivision_planes_ev3.txt";
if(n==8) h="time_timedivision_planes0.txt";

TGraph *MyGraph = new TGraph(h);
TGraph *MyGraphlog = new TGraph;

for (Int_t i=0; i<MyGraph->GetN() ;i++)
{
    if(MyGraph->GetPointX(i)>0 && MyGraph->GetPointY(i)>0.00001 && MyGraph->GetPointX(i)!=0) MyGraphlog->SetPoint(i,log(MyGraph->GetPointX(i)),log(MyGraph->GetPointY(i)));
}

//TF1 *fa2 = new TF1("fa2","abs([0])*x**[1]+(abs([0])*(x**[1]))**2",0,100);

TF1 *fa2 = new TF1("fa2","[0]*x+[1]",0,100);

if(n==5) fa2->SetParameters(5,-15); //quintuplets
if(n==4) fa2->SetParameters(5,-15); //quadruplets
if(n<=3 || n==7) fa2->SetParameters(5,-15); //triplets

MyGraphlog->Fit("fa2","","",0.5,100);

gStyle->SetOptFit(1);

if(n==1) MyGraphlog->SetTitle("Process Time (Triplet 5cm);log(nTPCpoints);log(time) [s]");
if(n==2) MyGraphlog->SetTitle("Process Time (Triplet Old);log(nTPCpoints);log(time) [s]");
if(n==3) MyGraphlog->SetTitle("Process Time (Triplet);log(nTPCpoints);log(time) [s]");
if(n==4) MyGraphlog->SetTitle("Process Time (Quadruplet);log(nTPCpoints);log(time) [s]");
if(n==5) MyGraphlog->SetTitle("Process Time (Quintuplet);log(nTPCpoints);log(time) [s]");
if(n==6) MyGraphlog->SetTitle("Process Time (Doublets);log(nTPCpoints);log(time) [s]");
if(n==7) MyGraphlog->SetTitle("Process Time (Time bunches);log(nTPCpoints(bunch));log(time) [s]");
if(n==8) MyGraphlog->SetTitle("Process Time (Time bunches tot);log(nTPCpoints(tot));log(time) [s]");

gStyle->SetStatX(0.5);
gStyle->SetStatY(0.9);

TCanvas *c2 = new TCanvas("c2","ac2", 200,10,700,500);

MyGraphlog->Draw("A*");


TString save;

if (n==1) save="./Plots/log_timetriplet5cm";
if (n==2) save="./Plots/log_timetripletold";
if (n==3) save="./Plots/log_timetriplet";
if (n==4) save="./Plots/log_timequadruplet";
if (n==5) save="./Plots/log_timequintuplet";
if (n==6) save="./Plots/log_timedoublet";
if (n==7) save="./Plots/log_timebunches";
if (n==8) save="./Plots/log_timebunches_tot";

if(n==7)
{
    Float_t sum=0;
    for (Int_t i=0; i<MyGraph->GetN() ;i++)
    {
      sum+=MyGraph->GetPointY(i);
    }
    std::cout<<"sum: "<<sum<<std::endl;

}


save+=".png";

c2->Print(save);

}