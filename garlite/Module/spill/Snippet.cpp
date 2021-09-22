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

void Snippet()
{
TH1F *hev1 = new TH1F("hev1", "Spill 1 time distribution", 1200, 0, 12000);
TH1F *hev2 = new TH1F("hev2", "SPill 2 time distribution", 1200, 0, 12000);
TH1F *hev3 = new TH1F("hev3", "Spill 3 time distribution", 1200, 0, 12000);
TH2F *hev1z = new TH2F("hev1z", "Spill 1 time distribution", 1200, 1200, 1800,1000,0,12000);
TH2F *hev2z = new TH2F("hev2z", "SPill 2 time distribution", 1200, 1200, 1800,1000,0,12000);
TH2F *hev3z = new TH2F("hev3z", "Spill 3 time distribution", 1200, 1200, 1800,1000,0,12000);
TH1F *h1 = new TH1F("h1", "Time distance", 200, 0, 2000);

TTree *MyTree3 = new TTree("MyTree3", "MyTree3");
MyTree3->ReadFile("text_files/time_spill_std_ev3_z_acuddg.txt", "Time:z");
TTree *MyTree2 = new TTree("MyTree2", "MyTree2");
MyTree2->ReadFile("text_files/time_spill_std_ev2_z_acuddg.txt", "Time:z");
TTree *MyTree1 = new TTree("MyTree1", "MyTree1");
MyTree1->ReadFile("text_files/time_spill_std_ev1_z_acuddg.txt", "Time:z");


Float_t Time;
Float_t z;
std::vector<Float_t> Times;

MyTree3->SetBranchAddress("Time",&Time);
MyTree3->SetBranchAddress("z",&z);
for(int i=0; i<MyTree3->GetEntries(); i++){
    MyTree3->GetEntry(i);
    Times.push_back(Time);
    hev3z->Fill(z,Time);
  }
sort(Times.begin(),Times.end());
for (int i=0; i<Times.size()-1; i++)
{
    h1->Fill(Times.at(i+1)-Times.at(i));
    hev3->Fill(Times.at(i));
}

Float_t Time2;
Float_t z2;
std::vector<Float_t> Times2;

MyTree2->SetBranchAddress("Time",&Time2);
MyTree2->SetBranchAddress("z",&z2);
for(int i=0; i<MyTree2->GetEntries(); i++){
    MyTree2->GetEntry(i);
    Times2.push_back(Time2);
    hev2z->Fill(z2,Time2);
  }
sort(Times2.begin(),Times2.end());
for (int i=0; i<Times2.size()-1; i++)
{
    h1->Fill(Times2.at(i+1)-Times2.at(i));
    hev2->Fill(Times2.at(i));
    std::cout<<Times2.at(i)<<std::endl;
}

Times.clear();

Float_t Time1;
Float_t z1;
std::vector<Float_t> Times1;
MyTree1->SetBranchAddress("Time",&Time1);
MyTree1->SetBranchAddress("z",&z1);
for(int i=0; i<MyTree1->GetEntries(); i++){
    MyTree1->GetEntry(i);
    Times1.push_back(Time1);
    hev1z->Fill(z1,Time1);
  }
sort(Times1.begin(),Times1.end());
for (int i=0; i<Times1.size()-1; i++)
{
    h1->Fill(Times1.at(i+1)-Times1.at(i));
    hev1->Fill(Times1.at(i));
}

TCanvas *cev1 = new TCanvas("cev1","acev1", 200,10,700,500);
hev1->SetTitle("Spill 1 time distribution;nTicks;nTPCClusters");
hev1->Draw();
cev1->Print("hev1_LAronly_adrewcg.png");

TCanvas *cev1z = new TCanvas("cev1z","acev1z", 200,10,700,500);
hev1z->SetTitle("Spill 1 time distribution;z;t");
hev1z->SetMarkerStyle(3);
hev1z->Draw();
cev1z->Print("hev1z_LAronly_adrewcg.png");

TCanvas *cev2 = new TCanvas("cev2","acev2", 200,10,700,500);
hev2->SetTitle("Spill 2 time distribution;nTicks;nTPCClusters");
hev2->Draw();
cev2->Print("hev2_LAronly_adrewcg.png");

TCanvas *cev2z = new TCanvas("cev2z","acev2z", 200,10,700,500);
hev2z->SetTitle("Spill 2 time distribution;z;t");
hev2z->SetMarkerStyle(3);
hev2z->Draw();
cev2z->Print("hev2z_LAronly_adrewcg.png");

TCanvas *cev3 = new TCanvas("cev3","acev3", 200,10,700,500);
hev3->SetTitle("Spill 3 time distribution;nTicks;nTPCClusters");
hev3->Draw();
cev3->Print("hev3_LAronly_adrewcg.png");

TCanvas *cev3z = new TCanvas("cev3z","acev3z", 200,10,700,500);
hev3z->SetTitle("Spill 3 time distribution;z;t");
hev3z->SetMarkerStyle(3);
hev3z->Draw();
cev3z->Print("hev3z_LAronly_adrewcg.png");

TCanvas *cdiff= new TCanvas("cdiff","acdiff", 200,10,700,500);
h1->SetTitle("Time difference distribution;nTicks;nTPCClusters");
h1->Draw();
cdiff->Print("hdiff_LAronly_adrewcg.png");
     
}
