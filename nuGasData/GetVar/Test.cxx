#include "stdio.h"
#include "TFile.h"
#include "TH1I.h"
#include "TSystem.h"
#include "TTree.h"
#include "TROOT.h"
#include "./include/GetVar.h"


void Test()
{
   float t= GetVar::GetX0();
   std::cout<<"t "<<t<<std::endl;
}