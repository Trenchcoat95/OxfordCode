#include "TROOT.h"
#include "TRandom.h"
#include "iostream"
#include "TMath.h"

void makehisto()
{
  TH1F *histo  = new TH1F("histo","0<#theta<25 , Standard ECAL;#theta (deg);#Delta E (GeV)",5,0,25);
  histo->SetBinContent(1,0.7605);
  histo->SetBinContent(2,0.762);
  histo->SetBinContent(3,0.7723);
  histo->SetBinContent(4,0.7948);
  histo->SetBinContent(5,0.8168);
  //histo->SetBinContent(6,1.051);
  //histo->SetBinContent(7,0.9271);
  //histo->SetBinContent(8,0.8898);
  

  histo->SetBinError(1,0.0015);
  histo->SetBinError(2,0.002);
  histo->SetBinError(3,0.0018);
  histo->SetBinError(4,0.0023);
  histo->SetBinError(5,0.0020);
  //histo->SetBinError(6,0.003);
  //histo->SetBinError(7,0.0023);
  //histo->SetBinError(8,0.0026);
  
  gStyle->SetOptStat(0);
  TCanvas *mccanvas2 = new TCanvas("mccanvas2","",1000,800);
  //mccanvas1->Divide(3,3);
  //mccanvas1->cd(1);
  //gStyle->SetOptFit(1);
  histo->SetFillColor(590);
  histo->Draw("E2");
  histo->Draw("E1same");
  mccanvas2->Print("anglehisto.png");
}
