//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 22 17:38:47 2018 by ROOT version 6.12/06
// from TTree GArAnaTree/GArAnaTree
// found on file: anatree.root
//////////////////////////////////////////////////////////

#ifndef garana_h
#define garana_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class garana {
public :
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Event;
   Int_t           SubRun;
   Int_t           Run;
   std::vector<int>     *NType;
   std::vector<int>     *CCNC;
   std::vector<int>     *PDG;
   std::vector<int>     *PDGMother;
   std::vector<float>   *MCVertX;
   std::vector<float>   *MCVertY;
   std::vector<float>   *MCVertZ;
   std::vector<float>   *MCNuPx;
   std::vector<float>   *MCNuPy;
   std::vector<float>   *MCNuPz;
   std::vector<float>   *MCPStartX;
   std::vector<float>   *MCPStartY;
   std::vector<float>   *MCPStartZ;
   std::vector<float>   *MCPStartPX;
   std::vector<float>   *MCPStartPY;
   std::vector<float>   *MCPStartPZ;
   std::vector<int>     *MCTrkID;
   std::vector<float>   *TrackLenF;
   std::vector<float>   *TrackLenB;
   std::vector<float>   *TrackStartX;
   std::vector<float>   *TrackStartY;
   std::vector<float>   *TrackStartZ;
   std::vector<float>   *TrackStartPX;
   std::vector<float>   *TrackStartPY;
   std::vector<float>   *TrackStartPZ;
   std::vector<int>     *TrackStartQ;
   std::vector<float>   *TrackEndX;
   std::vector<float>   *TrackEndY;
   std::vector<float>   *TrackEndZ;
   std::vector<float>   *TrackEndPX;
   std::vector<float>   *TrackEndPY;
   std::vector<float>   *TrackEndPZ;
   std::vector<int>     *TrackEndQ;
   std::vector<int>     *NTPCClustersOnTrack;
   std::vector<int>     *TrajMCPTrackID;
   std::vector<float>   *TrajMCPX;
   std::vector<float>   *TrajMCPY;
   std::vector<float>   *TrajMCPZ;
   std::vector<float>   *TrajMCPPX;
   std::vector<float>   *TrajMCPPY;
   std::vector<float>   *TrajMCPPZ;
    

   // List of branches
   TBranch        *b_Event;   //!
   TBranch        *b_SubRun;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_NType;   //!
   TBranch        *b_CCNC;   //!
   TBranch        *b_PDG;   //!
   TBranch        *b_PDGMother;
   TBranch        *b_TrackLenF;   //!
   TBranch        *b_TrackLenB;   //!
   TBranch        *b_MCVertX;
   TBranch        *b_MCVertY;
   TBranch        *b_MCVertZ;
   TBranch        *b_MCNuPx;
   TBranch        *b_MCNuPy;
   TBranch        *b_MCNuPz;
   TBranch        *b_MCPStartX;   //!
   TBranch        *b_MCPStartY;   //!
   TBranch        *b_MCPStartZ;   //!
   TBranch        *b_MCPStartPX;   //!
   TBranch        *b_MCPStartPY;   //!
   TBranch        *b_MCPStartPZ;   //!
   TBranch        *b_MCTrkID;
   TBranch        *b_TrackStartX;   //!
   TBranch        *b_TrackStartY;   //!
   TBranch        *b_TrackStartZ;   //!
   TBranch        *b_TrackStartPX;   //!
   TBranch        *b_TrackStartPY;   //!
   TBranch        *b_TrackStartPZ;   //!
   TBranch        *b_TrackStartQ;
   TBranch        *b_TrackEndX;   //!
   TBranch        *b_TrackEndY;   //!
   TBranch        *b_TrackEndZ;   //!
   TBranch        *b_TrackEndPX;   //!
   TBranch        *b_TrackEndPY;   //!
   TBranch        *b_TrackEndPZ;   //!
   TBranch        *b_TrackEndQ;
   TBranch        *b_NTPCClustersOnTrack;
   TBranch        *b_TrajMCPTrackID;
   TBranch        *b_TrajMCPPX;
   TBranch        *b_TrajMCPPY;
   TBranch        *b_TrajMCPPZ;
   TBranch        *b_TrajMCPX;
   TBranch        *b_TrajMCPY;
   TBranch        *b_TrajMCPZ;


   garana(TChain *tree=0);
   virtual ~garana();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Int_t    GetEntries();
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TChain *tree);
   //virtual void     Loop();
   virtual void     l2g_Trackmatch();
   virtual void     Plane_Eloss();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef garana_cxx
garana::garana(TChain *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("hemanatree.root");
      //if (!f || !f->IsOpen()) {
      //   f = new TFile("hemanatree.root");
      //}
      //TDirectory * dir = (TDirectory*)f->Get("hemanatree.root:/anatree");
      //dir->GetObject("GArAnaTree",tree);
      //TChain* point = new TChain("/anatree/GArAnaTree");
      //point->Add("/home/federico/Documents/Universita/Federico_2020-2021/OxfordCode/l2g/data/*.root");
      tree = new TChain("/anatree/GArAnaTree");
      //tree->Add("/home/federico/Documents/Universita/Federico_2020-2021/OxfordCode/l2g/data/*.root");
      //tree->Add("/home/federico/Documents/Universita/Federico_2020-2021/OxfordCode/l2g/Macros/Current/Standard/data/*.root");
      tree->Add("/home/federico/Documents/Universita/Federico_2020-2021/OxfordCode/l2g/Macros/Current/Standard/data/*.root");
      //neutrino.nd_hall_dayone_lar_SPY_v2_wMuID.volArgonCubeActive.Ev973000.Ev973999.2037.anatree.root

   }
   Init(tree);
}

garana::~garana()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t garana::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Int_t garana::GetEntries()
{
// Return Number of entries.
   if (!fChain) return 0;
   return fChain->GetEntries();
}

Long64_t garana::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void garana::Init(TChain *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   NType = 0;
   CCNC = 0;
   PDG = 0;
   PDGMother = 0;
   TrackLenF = 0;
   TrackLenB = 0;
   MCVertX = 0;
   MCVertY = 0;
   MCVertZ = 0;
   MCNuPx = 0;
   MCNuPy = 0;
   MCNuPz = 0;
   MCPStartX = 0;
   MCPStartY = 0;
   MCPStartZ = 0;
   MCPStartPX = 0;
   MCPStartPY = 0;
   MCPStartPZ = 0;
   MCTrkID = 0;
   TrackStartX = 0;
   TrackStartY = 0;
   TrackStartZ = 0;
   TrackStartPX = 0;
   TrackStartPY = 0;
   TrackStartPZ = 0;
   TrackStartQ = 0;
   TrackEndX = 0;
   TrackEndY = 0;
   TrackEndZ = 0;
   TrackEndPX = 0;
   TrackEndPY = 0;
   TrackEndPZ = 0;
   TrackEndQ = 0;
   NTPCClustersOnTrack = 0;
   TrajMCPTrackID = 0;
   TrajMCPPX = 0;
   TrajMCPPY = 0;
   TrajMCPPZ = 0;
   TrajMCPX = 0;
   TrajMCPY = 0;
   TrajMCPZ = 0;
   
   
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   
   fCurrent = -1;
   fChain->SetMakeClass(1);

   //fChain->SetBranchAddress("Event", &Event, &b_Event);
   //fChain->SetBranchAddress("SubRun", &SubRun, &b_SubRun);
   //fChain->SetBranchAddress("Run", &Run, &b_Run);
   //fChain->SetBranchAddress("NType", &NType, &b_NType);
   //fChain->SetBranchAddress("CCNC", &CCNC, &b_CCNC);
   fChain->SetBranchAddress("PDG", &PDG, &b_PDG);
   fChain->SetBranchAddress("PDGMother", &PDGMother, &b_PDGMother);
   //fChain->SetBranchAddress("TrackLenF", &TrackLenF, &b_TrackLenF);
   //fChain->SetBranchAddress("TrackLenB", &TrackLenB, &b_TrackLenB);
   fChain->SetBranchAddress("MCNuPx", &MCNuPx, &b_MCNuPx);
   fChain->SetBranchAddress("MCNuPy", &MCNuPy, &b_MCNuPy);
   fChain->SetBranchAddress("MCNuPz", &MCNuPz, &b_MCNuPz);
   fChain->SetBranchAddress("MCVertX", &MCVertX, &b_MCVertX);
   fChain->SetBranchAddress("MCVertY", &MCVertY, &b_MCVertY);
   fChain->SetBranchAddress("MCVertZ", &MCVertZ, &b_MCVertZ);
   fChain->SetBranchAddress("MCPStartX", &MCPStartX, &b_MCPStartX);
   fChain->SetBranchAddress("MCPStartY", &MCPStartY, &b_MCPStartY);
   fChain->SetBranchAddress("MCPStartZ", &MCPStartZ, &b_MCPStartZ);
   fChain->SetBranchAddress("MCPStartPX", &MCPStartPX, &b_MCPStartPX);
   fChain->SetBranchAddress("MCPStartPY", &MCPStartPY, &b_MCPStartPY);
   fChain->SetBranchAddress("MCPStartPZ", &MCPStartPZ, &b_MCPStartPZ);
   fChain->SetBranchAddress("MCTrkID", &MCTrkID, &b_MCTrkID);
   fChain->SetBranchAddress("TrackStartX", &TrackStartX, &b_TrackStartX);
   fChain->SetBranchAddress("TrackStartY", &TrackStartY, &b_TrackStartY);
   fChain->SetBranchAddress("TrackStartZ", &TrackStartZ, &b_TrackStartZ);
   fChain->SetBranchAddress("TrackStartPX", &TrackStartPX, &b_TrackStartPX);
   fChain->SetBranchAddress("TrackStartPY", &TrackStartPY, &b_TrackStartPY);
   fChain->SetBranchAddress("TrackStartPZ", &TrackStartPZ, &b_TrackStartPZ);
   fChain->SetBranchAddress("TrackStartQ", &TrackStartQ, &b_TrackStartQ);
   fChain->SetBranchAddress("TrackEndX", &TrackEndX, &b_TrackEndX);
   fChain->SetBranchAddress("TrackEndY", &TrackEndY, &b_TrackEndY);
   fChain->SetBranchAddress("TrackEndZ", &TrackEndZ, &b_TrackEndZ);
   fChain->SetBranchAddress("TrackEndPX", &TrackEndPX, &b_TrackEndPX);
   fChain->SetBranchAddress("TrackEndPY", &TrackEndPY, &b_TrackEndPY);
   fChain->SetBranchAddress("TrackEndPZ", &TrackEndPZ, &b_TrackEndPZ);
   fChain->SetBranchAddress("TrackEndQ", &TrackEndQ, &b_TrackEndQ);
   fChain->SetBranchAddress("NTPCClustersOnTrack", &NTPCClustersOnTrack, &b_NTPCClustersOnTrack);
   fChain->SetBranchAddress("TrajMCPTrackID", &TrajMCPTrackID, &b_TrajMCPTrackID);
   fChain->SetBranchAddress("TrajMCPX", &TrajMCPX, &b_TrajMCPX);
   fChain->SetBranchAddress("TrajMCPY", &TrajMCPY, &b_TrajMCPY);
   fChain->SetBranchAddress("TrajMCPZ", &TrajMCPZ, &b_TrajMCPZ);
   fChain->SetBranchAddress("TrajMCPPX", &TrajMCPPX, &b_TrajMCPPX);
   fChain->SetBranchAddress("TrajMCPPY", &TrajMCPPY, &b_TrajMCPPY);
   fChain->SetBranchAddress("TrajMCPPZ", &TrajMCPPZ, &b_TrajMCPPZ);
   Notify();
}

Bool_t garana::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void garana::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t garana::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef garana_cxx
