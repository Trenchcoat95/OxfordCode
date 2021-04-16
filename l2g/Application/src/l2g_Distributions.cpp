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



int Distribution(const char* finname, const char* foutname)
{

    //Units mm different from gdml which is in cm
	  //TChain* t = new TChain("EDepSimEvents","EDepSimEvents");
	  //t->Add(finname);
	  //TFile f(t->GetListOfFiles()->At(0)->GetTitle());

    TVector3 LArCenter(0,0,6660);
    TVector3 ArgonCubeActiveDims(7146.999999999999,3010.2,5091);
    TVector3 LArUpLimit = LArCenter+0.5*ArgonCubeActiveDims;
    TVector3 LArLowLimit = LArCenter-0.5*ArgonCubeActiveDims;

    //TVector3 ArgonCubeActiveDims(7146.999999999999,3010.2,5091);
    TH2* hxy = new TH2F("hxy", "hxy title", 1000, -10000, 10000, 1000, -10000, 10000);
    TH2* hyz = new TH2F("hyz", "hyz title", 1000, -10000, 10000, 1000, -10000, 10000);
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
    
    TG4Event* ev = new TG4Event;
    t->SetBranchAddress("Event",&ev);
  
    //TG4Event evout;
    std::vector<int> PDG;
    std::vector<double> X;
    std::vector<double> Y;
    std::vector<double> Z;
    std::vector<double> PX;
    std::vector<double> PY;
    std::vector<double> PZ;
    

    TFile fout(foutname,"RECREATE");
    TTree tout("tDistribution","Distribution");
    //TTree *tout = t->CloneTree(0);
    //tout.Branch("Event","TG4Event",&evout);
    tout.Branch("PDG","vector<int>",&PDG);
    tout.Branch("X","vector<double>",&X);
    tout.Branch("Y","vector<double>",&Y);
    tout.Branch("Z","vector<double>",&Z);
    tout.Branch("PX","vector<double>",&PX);
    tout.Branch("PY","vector<double>",&PY);
    tout.Branch("PZ","vector<double>",&PZ);
    
    const int nev = t->GetEntries();
    
    std::cout << "Events: " << nev << " [";
    std::cout << std::setw(3) << int(0) << "%]" << std::flush;

    for(int i = 0; i < nev; i++)
    {

      t->GetEntry(i);

      

      size_t ntraj = ev->Trajectories.size();
      for(size_t j=0; j<ntraj; j++)
      {
        TVector3 FirstPositionOutofLAr(0,0,0);
        TVector3 FirstMomentunOutofLAr(0,0,0);
        int PDGcode = ev->Trajectories.at(j).GetPDGCode();
        bool WasInLAr=false;
        bool GotOutofLAr=false;
        for(size_t k=0; k<ev->Trajectories.at(j).Points.size();k++)
        {
          TLorentzVector CurrentPositionL = ev->Trajectories.at(j).Points.at(k).GetPosition();
          TVector3 CurrentPosition(CurrentPositionL.X(),CurrentPositionL.Y(),CurrentPositionL.Z());
          TVector3 CurrentMomentum = ev->Trajectories.at(j).Points.at(k).GetMomentum();


          if(CurrentPosition.X()<LArUpLimit.X() && CurrentPosition.Y()<LArUpLimit.Y() && CurrentPosition.Z()<LArUpLimit.Z() 
             && CurrentPosition.X()>LArLowLimit.X() && CurrentPosition.Y()>LArLowLimit.Y() && CurrentPosition.Z()>LArLowLimit.Z())
             {WasInLAr=true;}

          if((CurrentPosition.X()>LArUpLimit.X() || CurrentPosition.Y()>LArUpLimit.Y() || CurrentPosition.Z()>LArUpLimit.Z() 
             && CurrentPosition.X()<LArLowLimit.X() || CurrentPosition.Y()<LArLowLimit.Y() || CurrentPosition.Z()<LArLowLimit.Z())&& WasInLAr && !GotOutofLAr)
             {
              FirstPositionOutofLAr=CurrentPosition;
              FirstMomentunOutofLAr=CurrentMomentum;
              GotOutofLAr=true;
             }  
     
        }
        if(WasInLAr && GotOutofLAr)
          {
          //std::cout<<"FirstPositionOutofLAr (t,x,y,z) "<<FirstPositionOutofLAr.X()<<" "<<FirstPositionOutofLAr.Y()<<" "<<FirstPositionOutofLAr.Z()<<std::endl;
          PDG.push_back(PDGcode);
          X.push_back(FirstPositionOutofLAr.X());
          Y.push_back(FirstPositionOutofLAr.Y());
          Z.push_back(FirstPositionOutofLAr.Z());
          PX.push_back(FirstMomentunOutofLAr.X());
          PY.push_back(FirstMomentunOutofLAr.Y());
          PZ.push_back(FirstMomentunOutofLAr.Z());
          hxy->Fill(FirstPositionOutofLAr.X(),FirstPositionOutofLAr.Y());
          hyz->Fill(FirstPositionOutofLAr.Y(),FirstPositionOutofLAr.Z());
          }
        
      }
          
      std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i)/nev*100) << "%]" << std::flush;
         
      tout.Fill();

      PDG.clear();
      X.clear();
      Y.clear();
      Z.clear();
      PX.clear();
      PY.clear();
      PZ.clear();
    }
    std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
    std::cout << std::endl;
    
    //fout.cd();
    
    
    //hxy->Write();
    //hyz->Write();
    tout.Write();
    //t->CloneTree()->Write();
	return 0;
}

void help_select()
{
  std::cout << "GetDistribution <input file> <output file>" << std::endl;
  std::cout << "input file name could contain wild card" << std::endl;
} 

int main(int argc, char* argv[])
{

  if(argc != 3)
    help_select();
  else
    Distribution(argv[1], argv[2]);
}

