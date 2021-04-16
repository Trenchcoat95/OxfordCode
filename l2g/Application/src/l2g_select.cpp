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

    TG4Event* ev = new TG4Event;
    t->SetBranchAddress("Event",&ev);
  
    
    TFile fout(foutname,"RECREATE");
    //TTree tout("tSelect","Selection");
    TTree *tout = t->CloneTree(0);
    //tout.Branch("cell","std::vector<cell>",&vec_cell);
    
    
    const int nev = t->GetEntries();
    
    std::cout << "Events: " << nev << " [";
    std::cout << std::setw(3) << int(0) << "%]" << std::flush;

    for(int i = 0; i < nev; i++)
    {

      t->GetEntry(i);

 
      bool copyevent= false;


      for (std::map<std::string,std::vector<TG4HitSegment> >::iterator it=ev->SegmentDetectors.begin();
            it!=ev->SegmentDetectors.end(); ++it)
      {
        //std::string s = it->first;
        //const char *cstr = it.c_str();
        //std::cout << "Segment detectors: " << it->first << std::endl;
        if (it->first!= "ArgonCube") copyevent=true;
      }
      
      //std::cout<<std::endl;
      if (copyevent) tout->Fill();    
      std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i)/nev*100) << "%]" << std::flush;
         
    }
    std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
    std::cout << std::endl;
    

    
    geo->Write();
    tout->Write();
    //t->CloneTree()->Write();
    TDirectory* cdir = fout.mkdir("DetSimPassThru");
    fout.cd("DetSimPassThru");
    gRooTracker->CloneTree()->Write();
    InputKinem->CloneTree()->Write();
    InputFiles->CloneTree()->Write();
    //fout.Close();
    
    //f.Close();
    //Tutto molto bello
	return 0;
}


int SelectGeo(const char* finname, const char* foutname)
{

  TVector3 LArCenter(0,0,6660);
  TVector3 ArgonCubeActiveDims(7146.999999999999,3010.2,5091);
  TVector3 LArUpLimit = LArCenter+0.5*ArgonCubeActiveDims;
  TVector3 LArLowLimit = LArCenter-0.5*ArgonCubeActiveDims;


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

    TG4Event* ev = new TG4Event;
    t->SetBranchAddress("Event",&ev);
  
    
    TFile fout(foutname,"RECREATE");
    //TTree tout("tSelect","Selection");
    TTree *tout = t->CloneTree(0);
    //tout.Branch("cell","std::vector<cell>",&vec_cell);
    
    
    const int nev = t->GetEntries();
    
    std::cout << "Events: " << nev << " [";
    std::cout << std::setw(3) << int(0) << "%]" << std::flush;

    for(int i = 0; i < nev; i++)
    {

      t->GetEntry(i);

 
      bool copyevent= false;


      size_t ntraj = ev->Trajectories.size();
      for(size_t j=0; j<ntraj; j++)
      {
        for(size_t k=0; k<ev->Trajectories.at(j).Points.size();k++)
        {
          TLorentzVector CurrentPositionL = ev->Trajectories.at(j).Points.at(k).GetPosition();
          TVector3 CurrentPosition(CurrentPositionL.X(),CurrentPositionL.Y(),CurrentPositionL.Z());
          TVector3 CurrentMomentum = ev->Trajectories.at(j).Points.at(k).GetMomentum();

          if((CurrentPosition.X()>LArUpLimit.X() || CurrentPosition.Y()>LArUpLimit.Y() || CurrentPosition.Z()>LArUpLimit.Z() 
             && CurrentPosition.X()<LArLowLimit.X() || CurrentPosition.Y()<LArLowLimit.Y() || CurrentPosition.Z()<LArLowLimit.Z()))
             {
              copyevent=true;
              break;
             }  
     
        }
        if(copyevent) break;
      }

      if (copyevent) tout->Fill();    
      std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i)/nev*100) << "%]" << std::flush;
         
    }
    std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
    std::cout << std::endl;
    

    
    geo->Write();
    tout->Write();
    //t->CloneTree()->Write();
    TDirectory* cdir = fout.mkdir("DetSimPassThru");
    fout.cd("DetSimPassThru");
    gRooTracker->CloneTree()->Write();
    InputKinem->CloneTree()->Write();
    InputFiles->CloneTree()->Write();
    //fout.Close();
    
    //f.Close();
    //Tutto molto bello
	return 0;
}

void help_select()
{
  std::cout << "Select <input file> <output file>" << std::endl;
  std::cout << "input file name could contain wild card" << std::endl;
} 

int main(int argc, char* argv[])
{

  if(argc != 3)
    help_select();
  else
    SelectGeo(argv[1], argv[2]);
}

