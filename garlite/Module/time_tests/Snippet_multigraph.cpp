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

void Snippet_multigraph()
{
TGraph *MyGraph = new TGraph("time.txt");
TGraph *MyGraph2 = new TGraph("timeold.txt");


//TF1 *fa2 = new TF1("fa2","abs([0])*x**[1]+(abs([0])*(x**[1]))**2",0,100);

TF1 *fa2 = new TF1("fa2","exp([0])*x**[1]",0,100);
TF1 *fa3 = new TF1("fa3","exp([0])*x**[1]",0,100);

gStyle->SetOptFit(0);

fa2->SetParameters(-14.6,3.94); //triplets
//MyGraph->Fit("fa2","","",0,100);
fa3->SetLineColor(kGreen);
//MyGraph2->Fit("fa2","","",0,100);
fa3->SetParameters(-14.8,4.1);



MyGraph2->SetMarkerColor(kBlue);

TMultiGraph *mg = new TMultiGraph();
mg->Add(MyGraph);
mg->Add(MyGraph2);
TCanvas *c2 = new TCanvas("c2","Temperature", 200,10,700,500);
gPad->SetLogx(1);
gPad->SetLogy(1);

mg->SetTitle("Process Time (Old VS New);nTPCPoints;time [s]");
mg->Draw("A*");
fa2->Draw("same");
fa3->Draw("same");
c2->Print("./Plots/log_timetripletoldVSnew_log.png");

}