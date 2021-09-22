#include "GetVar.h"


bool GetVar::CreateOrganicCompound(gasMixer & gQuencher, int nC, int nH, int nF, gasMixer gC, gasMixer gH, gasMixer gF)
{
    if(nC==1) gQuencher=gQuencher+gC;
    else if (nC>1) gQuencher=gQuencher+gC*nC;
    if(nH==1) gQuencher=gQuencher+gH;
    else if (nH>1) gQuencher=gQuencher+gH*nH;
    if(nF==1) gQuencher=gQuencher+gF;
    else if (nF>1) gQuencher=gQuencher+gF*nF;
    if(nH==0&&nC==0&&nF==0) 
    {
        std::cout<<"Can't have a null quencher!"<<std::endl;
        return 0;
    }
    return 1;
}

float GetVar::GetSigmaT(int fquencher, int p, int nC, int nH, int Zgas, int Agas, int E, int B)
{ 
    //Setup input file and tree
    TFile * fin=new TFile("all_merged.root");
    TTree * gastree=dynamic_cast<TTree*>(fin->Get("gastree"));

    //Define the cut to get the desired value
    TString cut;
    cut = Form("fquencher==%d && p==%d && nC==%d && nH==%d && Zgas==%d && Agas==%d && E==%d && B==%d",fquencher,p,nC,nH,Zgas,Agas,E,B);

    //obtain the ttree entry correspondent to our parameters
    gastree->Draw("Entry$>>hist(Entries$,0,Entries$)",cut,"goff");
    TH1I *hist = (TH1I*)gDirectory->Get("hist");
    Long64_t iEntry = hist->GetBinLowEdge(hist->FindFirstBinAbove(0));
    //std::cout<<"nEntries: "<<hist->GetEntries()<<"; Entry: "<<iEntry<<std::endl;

    if (hist->GetEntries()==0)
    {
        std::cout<<"No entry matches your criteria"<<std::endl;
        return 0;
    }
    if (hist->GetEntries()==0)
    {
        std::cout<<"More than one entry matches your criteria"<<std::endl;
        return 0;
    }
    
    
    //Setup the sigma diffusion branch
    float dt;
    gastree->SetBranchAddress("dt",&dt);
    gastree->GetEntry(iEntry);

    return dt;
}

float GetVar::GetSigmaTErr(int fquencher, int p, int nC, int nH, int Zgas, int Agas, int E, int B)
{ 
    //Setup input file and tree
    TFile * fin=new TFile("all_merged.root");
    TTree * gastree=dynamic_cast<TTree*>(fin->Get("gastree"));

    //Define the cut to get the desired value
    TString cut;
    cut = Form("fquencher==%d && p==%d && nC==%d && nH==%d && Zgas==%d && Agas==%d && E==%d && B==%d",fquencher,p,nC,nH,Zgas,Agas,E,B);

    //obtain the ttree entry correspondent to our parameters
    gastree->Draw("Entry$>>hist(Entries$,0,Entries$)",cut,"goff");
    TH1I *hist = (TH1I*)gDirectory->Get("hist");
    Long64_t iEntry = hist->GetBinLowEdge(hist->FindFirstBinAbove(0));
    //std::cout<<"nEntries: "<<hist->GetEntries()<<"; Entry: "<<iEntry<<std::endl;

    if (hist->GetEntries()==0)
    {
        std::cout<<"No entry matches your criteria"<<std::endl;
        return 0;
    }
    if (hist->GetEntries()==0)
    {
        std::cout<<"More than one entry matches your criteria"<<std::endl;
        return 0;
    }
    
    
    //Setup the sigma diffusion branch
    float dterr;
    gastree->SetBranchAddress("dterr",&dterr);
    gastree->GetEntry(iEntry);

    return dterr;
}


float GetVar::GetX0(int fquencher, int p , int nC, int nH, int Zgas, int Agas) //fquencher in percentage, p in mBar2
{
    gSystem->Load("libgasElement.so");
    gSystem->Load("libgasMixer.so");
    gSystem->Load("libGasUtils.so");

    //Define building blocks
    const gasMixer gC(6, 12), gH(1, 1), gF(9,19);
    
    //Define main gas
    gasMixer gMain(Zgas,Agas);
    
    //Define quencher
    gasMixer gQuencher;
    CreateOrganicCompound(gQuencher,nC,nH,0,gC,gH,gF);
    
    //Define mixture
    const gasMixer gTotal(gMain, gQuencher, fquencher/100., p/1000.);

    //gMain.Print();
    //gQuencher.Print();
    //gTotal.Print();
    
    return gTotal.getX0();    
}

float GetVar::GetX0F(int fquencher, int p , int nC, int nH, int nF, int Zgas, int Agas, int mC, int mH, int mF) //fquencher in percentage, p in mBar2
{
    gSystem->Load("libgasElement.so");
    gSystem->Load("libgasMixer.so");
    gSystem->Load("libGasUtils.so");

    //Define building blocks
    const gasMixer gC(6, 12), gH(1, 1), gF(9,19);
    
    //Define main gas
    gasMixer gMain(Zgas,Agas);
    if (mC!=0 || mH!=0 || mF!=0) CreateOrganicCompound(gMain,mC,mH,mF,gC,gH,gF);
    
    //Define quencher
    gasMixer gQuencher;
    CreateOrganicCompound(gQuencher,nC,nH,nF,gC,gH,gF);
    
    //Define mixture
    const gasMixer gTotal(gMain, gQuencher, fquencher/100., p/1000.);

    //gMain.Print();
    //gQuencher.Print();
    //gTotal.Print();
    
    return gTotal.getX0();    
}

