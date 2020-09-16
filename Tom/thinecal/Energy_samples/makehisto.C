#include "TROOT.h"
#include "TRandom.h"
#include "iostream"
#include "TMath.h"

void makehisto()
{
  TH1F *histo  = new TH1F("histo","0.5<P<5 GeV, Standard ECAL;p start (GeV/c);#Delta E (GeV)",9,0.5,5);
  histo->SetBinContent(1,0.5389);
  histo->SetBinContent(2,0.6875);
  histo->SetBinContent(3,0.7381);
  histo->SetBinContent(4,0.7345);
  histo->SetBinContent(5,0.7592);
  histo->SetBinContent(6,0.7598);
  histo->SetBinContent(7,0.7621);
  histo->SetBinContent(8,0.7616);
  histo->SetBinContent(9,0.7657);

  histo->SetBinError(1,0.0076);
  histo->SetBinError(2,0.0010);
  histo->SetBinError(3,0.0021);
  histo->SetBinError(4,0.0025);
  histo->SetBinError(5,0.0015);
  histo->SetBinError(6,0.0015);
  histo->SetBinError(7,0.0018);
  histo->SetBinError(8,0.0019);
  histo->SetBinError(9,0.0018); 

  TCanvas *mccanvas2 = new TCanvas("mccanvas2","",1000,800);
  //mccanvas1->Divide(3,3);
  //mccanvas1->cd(1);
  gStyle->SetOptStat(0);
  histo->SetFillColor(590);
  histo->Draw("E2");
  histo->Draw("E1same");
  mccanvas2->Print("energyhisto.png");
}
