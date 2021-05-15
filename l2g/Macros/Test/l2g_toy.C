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


void l2g_toy()
{   
    TF1 *fgaus1= new TF1("fgaus2","gaus",-0.4,0.4);
    TF1 *fgaus2= new TF1("fgaus2","gaus",-0.4,0.4);

    fgaus1->SetParameters(1,-0.02,0.01);
    fgaus2->SetParameters(1,0,0.04);

    TH1F* frac_resid = new TH1F("frac_resid", "Fractional residuals", 150, -0.4, 0.4);


    for (int i=0; i<100000; i++) frac_resid->Fill(fgaus1->GetRandom());
    for (int i=0; i<30000; i++) frac_resid->Fill(fgaus2->GetRandom());
    
    TF1 *double_gauss = new TF1("double_gauss","0.39894228040143*0.005*([0]/[2])*(exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-([1]+abs([4])))/[5])^2)*([2]/[5]))",-0.4,0.4);
    
    ////////////////////////////////////////////////////Fractional Residual 1D plot
    TCanvas *mccanvas_fracresid = new TCanvas("mccanvas_fracresid","",1000,800);
    gStyle->SetOptStat(1);
    frac_resid->SetTitle("Momentum fractional residuals (Double Gauss Fit);(p_{reco}-p_{true})/p_{true};n(#mu in GAr-Lite)");
    //double_gauss->SetParameters(100000,-0.05,0.01,0.3,0,0.04);
    double_gauss->SetParameters(100000,-0.01,0.02,0.4,0.05,0.05);
    frac_resid->Fit("double_gauss");
    gStyle->SetOptFit(1);
    frac_resid->Draw();
    TF1 *gauss1 = new TF1("gauss1","0.39894228040143*0.005*([0]/[2])*(exp(-0.5*((x-[1])/[2])^2))",-0.4,0.4);
    gauss1->SetParameters(double_gauss->GetParameter(0),double_gauss->GetParameter(1),double_gauss->GetParameter(2));
    gauss1->SetLineColor(kBlue);
    gauss1->SetLineStyle(9);
    gauss1->Draw("SAME");
    TF1 *gauss2 = new TF1("gauss2","0.39894228040143*0.005*([0]/[2])*(exp(-0.5*((x-[1])/[2])^2))",-0.4,0.4);
    gauss2->SetParameters(double_gauss->GetParameter(0)*double_gauss->GetParameter(3),double_gauss->GetParameter(1)+abs(double_gauss->GetParameter(4)),double_gauss->GetParameter(5));
    gauss2->SetLineColor(kBlue);
    gauss2->SetLineStyle(9);
    gauss2->Draw("SAME");   

}