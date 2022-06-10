#pragma once
#ifndef MATERIAL_PART_UTILS_H_
#define MATERIAL_PART_UTILS_H_
#include <iostream>
#include "TVectorD.h"
#include "TMatrix.h"
#include "TMath.h"
#include "TVector3.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TF2.h"
#include "Math/Vector3D.h"
#include <algorithm>
#include <vector>

using namespace std;

using namespace ROOT::Math;



namespace utils
{
    XYZVector GArCenter(0,-150.473,1486); 
    double GAr_r = 349.9;
    double GAr_L = 669.6;
    double Plane_thick = 4;
    double  B=-0.5; //in Tesla
    XYZVector Z_axis(0,0,1);
    XYZVector Y_axis(0,1,0);
    double ZA =0.54141;
    double Ipar=64.7e-9;      ///GeV
    double rho= 1.032;   //g/cm^3
    double X1=2.49;
    double X0=0.1469;
    double muon_mass=0.1056583755; //GeV/c^2   ///material properties are for typical plastic polymere scintillator: polyvinyltoluene
    double a=0.1610;
    double m=3.24;
    double hw=21.54e-9;
    double Z=0.085+6*0.915;
    double xx0 = 42.54;  //radiation length in cm (different from Alice's code where it's xx0=1/X0 in cm^-1)
    const float mK  = 0.307075e-3; // [GeV*cm^2/g]
    const float me  = 0.511e-3;    // [GeV/c^2]
    const Double_t kAlmost1=1. - Double_t(FLT_EPSILON);
    const Double_t kAlmost0=Double_t(FLT_MIN);
}

#endif