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



class kalana {
public :
    TTree  *t;
    TFile  *f;

    float xht;
    float yht;
    float zht;
    float xpost;
    int ev;
    TVectorF *parvect;
    TVectorF *predstept;
    TMatrixF *Pt;
    TMatrixF *PPredt;
    TMatrixF *Rt;

    TBranch *b_xht;
    TBranch *b_yht;
    TBranch *b_zht;
    TBranch *b_xpost;
    TBranch *b_ev;
    TBranch *b_parvect;
    TBranch *b_predstept;
    TBranch *b_Pt;
    TBranch *b_PPredt;
    TBranch *b_Rt;

    kalana(const char* treename, const char *filename);
    virtual ~kalana();
    virtual void  Init(TTree *tree);
    //virtual void  Loop(TGraphErrors *h, TGraphErrors *pred, TGraphErrors *par, Int_t nentries);
    //virtual void DoStuff();
    virtual void Loop();
    //virtual void Loop(float xh, float yh, TMatrixF *R, float xpos, TVectorF *predstep, TMatrixF *PPred, TVectorF *parvec, TMatrixF *P);
};

kalana::kalana(const char* treename, const char *filename)
{
    f = new TFile(filename);
    t = (TTree*)f->Get(treename);
    Init(t);
}

kalana::~kalana()
{
   if (!t) return;
   delete t->GetCurrentFile();
}

void kalana::Init(TTree *tree)
{
    xht = 0;
    yht = 0;
    zht = 0;
    xpost = 0;
    ev = 0;
    parvect = 0;
    predstept = 0;
    Pt = 0;
    PPredt = 0;
    Rt = 0;

    b_xht = 0;
    b_yht = 0;
    b_zht = 0;
    b_xpost = 0;
    b_ev = 0;
    b_parvect = 0;
    b_predstept = 0;
    b_Pt = 0;
    b_PPredt = 0;
    b_Rt = 0;

    tree->SetBranchAddress("xht",&xht,&b_xht);
    tree->SetBranchAddress("yht",&yht,&b_yht);
    tree->SetBranchAddress("zht",&zht,&b_zht);
    tree->SetBranchAddress("ev",&ev,&b_ev);
    tree->SetBranchAddress("xpost",&xpost,&b_xpost);
    tree->SetBranchAddress("parvect",&parvect,&b_parvect);
    tree->SetBranchAddress("predstept",&predstept,&b_predstept);
    tree->SetBranchAddress("Pt",&Pt,&b_Pt);
    tree->SetBranchAddress("PPredt",&PPredt,&b_PPredt);
    tree->SetBranchAddress("Rt",&Rt,&b_Rt);   
}
