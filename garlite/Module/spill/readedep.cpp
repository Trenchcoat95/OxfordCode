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
#include <TApplication.h> 
#include <TChain.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TDirectoryFile.h>

#include "/home/federico/Documents/edep-sim/edep-gcc-9-x86_64-linux-gnu/include/EDepSim/TG4Event.h"
#include "/home/federico/Documents/edep-sim/edep-gcc-9-x86_64-linux-gnu/include/EDepSim/TG4HitSegment.h"


#include <vector>
#include <map>
#include <iostream>
#include <iomanip>




void readedep(const char* fIn)    ///////(First file, List of file in directory using wildcard *, output file)
{
	
	//////////////////LOADING THE LIBRARIES AND TREES///////////////////////////////
    gSystem->Load("/home/federico/Documents/edep-sim/edep-gcc-9-x86_64-linux-gnu/lib/libedepsim_io.so");
    gSystem->Load("/home/federico/Documents/edep-sim/edep-gcc-9-x86_64-linux-gnu/lib/libedepsim.so");

    TH2F *htz = new TH2F("htz", "htz", 250, 12000, 18000,250,1e+12,1.2e+13);

    TFile f(fIn,"READ");
    TTree* t = (TTree*) f.Get("EDepSimEvents");
    TGeoManager* geo = (TGeoManager*) f.Get("EDepSimGeometry");

    
    TG4Event* ev = new TG4Event;
    t->SetBranchAddress("Event",&ev);
    const int nev = t->GetEntries();
    
    std::cout << "Events: " << nev << "\n";

    for(int i = 0; i < nev; i++)
    {
      t->GetEntry(i);
      for (std::map<std::string,std::vector<TG4HitSegment> >::iterator it=ev->SegmentDetectors.begin();
            it!=ev->SegmentDetectors.end(); ++it)
      {
          std::cout<<it->first<<std::endl;
          if(it->first=="Tracker_vol")
          {
              for(unsigned int j = 0; j < it->second.size(); j++)
              {
                  TG4HitSegment i=it->second.at(j);
                  //std::cout<<"X: "<<(i.GetStart().Z()+i.GetStop().Z())/2<< " Time: "<<(i.GetStart().T()+i.GetStop().T())/2<<std::endl;
                  float z=((i.GetStart().Z()+i.GetStop().Z())/2);
                  float t=(i.GetStart().T()+i.GetStop().T())/2;
                  htz->Fill(z,t);
              }
          }
      }
    }
	htz->Draw("colz");
	
}



