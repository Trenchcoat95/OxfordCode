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
#include <TGraph2D.h>

void kalana::Loop(std::string folder,std::string  sm,std::string  R)
{
Int_t nentries = (Int_t)t->GetEntries();

TGraph2D *XYh =new  TGraph2D(nentries-1);

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

      t->GetEntry(i);
      
      XYh->SetPoint(i-1,xht,yht,zht);
      //std::cout<<xht<<" "<<yht<<" "<<zht<<std::endl;

      Yh->SetPoint(i-1,xht,yht);
      //Yh->SetPointError(i-1,0,TMath::Sqrt((*Rt)[0][0]));
      Ypred->SetPoint(i-1,xpost,(*predstept)[0]);
      Ypred->SetPointError(i-1,0,TMath::Sqrt((*PPredt)[0][0]));
      Ypar->SetPoint(i-1,xpost,(*parvect)[0]);
      Ypar->SetPointError(i-1,0,TMath::Sqrt((*Pt)[0][0]));

      Zh->SetPoint(i-1,xht,zht);
      //Zh->SetPointError(i-1,0,TMath::Sqrt((*Rt)[1][1]));
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


TCanvas *mccanvas = new TCanvas("mccanvas","",1000,800);
XYh->SetTitle("XY;x(cm);y(cm);z(cm)");
XYh->SetLineColor(kRed);
XYh->SetMarkerStyle(3);
XYh->Draw();
std::string st = folder + "helix_XY_rndx_sm" + sm + "_R_" + R + "_stdK.png";
mccanvas->Print(st.c_str());


TCanvas *mccanvas1 = new TCanvas("mccanvas1","",1000,800);
Yh->SetTitle("Y;x(cm);y(cm)");
Yh->SetLineColor(kRed);
Yh->Draw();
Ypred->SetLineColor(kBlue);
Ypred->Draw("same");
Ypar->SetLineColor(kGreen);
Ypar->Draw("same");
auto legend = new TLegend(0.55,0.75,0.9,0.9);
legend->AddEntry(Yh,"TPC Cluster measured value","lep");
legend->AddEntry(Ypred,"A priori prediction","lep");
legend->AddEntry(Ypar,"A posteriori prediction","lep");
legend->Draw();
st = folder + "helix_Y_rndx_sm" + sm + "_R_" + R + "_stdK.png";
mccanvas1->Print(st.c_str());



TCanvas *mccanvas2 = new TCanvas("mccanvas2","",1000,800);
Zh->SetTitle("Z;x(cm);z(cm)");
Zh->SetLineColor(kRed);
Zh->Draw();
Zpred->SetLineColor(kBlue);
Zpred->Draw("same");
Zpar->SetLineColor(kGreen);
Zpar->Draw("same");
auto legend2 = new TLegend(0.325,0.75,0.675,0.9);
legend2->AddEntry(Zh,"TPC Cluster measured value","lep");
legend2->AddEntry(Zpred,"A priori prediction","lep");
legend2->AddEntry(Zpar,"A posteriori prediction","lep");
legend2->Draw();
st = folder + "helix_Z_rndx_sm" + sm + "_R_" + R + "_stdK.png";
mccanvas2->Print(st.c_str());



TCanvas *mccanvas3 = new TCanvas("mccanvas3","",1000,800);
Curvpred->SetTitle("Curv;x(cm);1/r(cm^{-1})");
Curvpred->SetLineColor(kBlue);
Curvpred->Draw();
Curvpar->SetLineColor(kGreen);
Curvpar->Draw("same");
auto legend3 = new TLegend(0.55,0.75,0.9,0.9);
legend3->AddEntry(Curvpred,"A priori prediction","lep");
legend3->AddEntry(Curvpar,"A posteriori prediction","lep");
legend3->Draw();
st = folder + "helix_Curv_rndx_sm" + sm + "_R_" + R + "_stdK.png";
mccanvas3->Print(st.c_str());

TCanvas *mccanvas4 = new TCanvas("mccanvas4","",1000,800);
Phipred->SetTitle("Phi;x(cm);#phi");
Phipred->SetLineColor(kBlue);
Phipred->Draw();
Phipar->SetLineColor(kGreen);
Phipar->Draw("same");
auto legend4 = new TLegend(0.325,0.75,0.675,0.9);
legend4->AddEntry(Phipred,"A priori prediction","lep");
legend4->AddEntry(Phipar,"A posteriori prediction","lep");
legend4->Draw();
st = folder + "helix_Phi_rndx_sm" + sm + "_R_" + R + "_stdK.png";
mccanvas4->Print(st.c_str());


TCanvas *mccanvas5 = new TCanvas("mccanvas5","",1000,800);
Lambdapred->SetTitle("Lambda;x(cm);#lambda");
Lambdapred->SetLineColor(kBlue);
Lambdapred->Draw();
Lambdapar->SetLineColor(kGreen);
Lambdapar->Draw("same");
auto legend5 = new TLegend(0.55,0.75,0.9,0.9);
legend5->AddEntry(Lambdapred,"A priori prediction","lep");
legend5->AddEntry(Lambdapar,"A posteriori prediction","lep");
legend5->Draw();
st = folder + "helix_Lambda_rndx_sm" + sm + "_R_" + R + "_stdK.png";
mccanvas5->Print(st.c_str());
}