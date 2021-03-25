#include <TGeoManager.h>
#include <TString.h>
#include <TGeoNode.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TGeoTrd2.h>
#include <TGeoTube.h>
#include <TChain.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TSystem.h>
#include "TG4Event.h"
#include "TG4HitSegment.h"

#include <vector>
#include <map>
#include <iostream>



int Select(const char* finname, const char* foutname)
{
	  //TChain* t = new TChain("EDepSimEvents","EDepSimEvents");
	  //t->Add(finname);
	  //TFile f(t->GetListOfFiles()->At(0)->GetTitle());
    TFile f(finname,"READ");
	if (! f.IsOpen() )
	{
		std::cerr << "File " << finname << " not found!" << std::endl;
		return -1;
	}
    TTree* t = (TTree*) f.Get("EDepSimEvents");
	if (!t)
	{
		std::cerr << "Tree EDepSimEvents not found" << std::endl;
		return -2;
	}
    TGeoManager* geo = (TGeoManager*) f.Get("EDepSimGeometry"); 
	if (!geo)
	{
		std::cerr << "TGeoManager EDepSimGeometry not found" << std::endl;
		return -3;
	}
    TTree* gRooTracker = (TTree*) f.Get("DetSimPassThru/gRooTracker");
    TTree* InputKinem = (TTree*) f.Get("DetSimPassThru/InputKinem");
    TTree* InputFiles = (TTree*) f.Get("DetSimPassThru/InputFiles");
    //std::cout << "checkpoint #0: Digitize" << std::endl;
    //int retCode = init(geo);
    //std::cout << "checkpoint #1: Digitize" << std::endl;
    TG4Event* ev = new TG4Event;
    t->SetBranchAddress("Event",&ev);
  
    
    TFile fout(foutname,"RECREATE");
    TTree tout("tSelect","Selection");
    //tout.Branch("cell","std::vector<cell>",&vec_cell);
    
    
    const int nev = t->GetEntries();
    
    std::cout << "Events: " << nev << " [";
    std::cout << std::setw(3) << int(0) << "%]" << std::flush;

    for(int i = 0; i < nev; i++)
    {
      t->GetEntry(i);
    
      std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i)/nev*100) << "%]" << std::flush;
         
      //tout.Fill();
    }
    std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
    std::cout << std::endl;
    
    fout.cd();
    //tout.Write();
    //geo->Write();
    //t->CloneTree()->Write();
    //gRooTracker->CloneTree()->Write();
    //InputKinem->CloneTree()->Write();
    //InputFiles->CloneTree()->Write();
    fout.Close();
    
    f.Close();
	return 0;
}

void help_select()
{
  std::cout << "l2g_Select <input file> <output file>" << std::endl;
  std::cout << "input file name could contain wild card" << std::endl;
} 

int main(int argc, char* argv[])
{

  if(argc != 3)
    help_select();
  else
    Select(argv[1], argv[2]);
}
