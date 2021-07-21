#ifndef _GASELEMENT_
#define _GASELEMENT_

#include <iostream>
#include "TString.h"
#include "GasUtils.h"

using namespace std;
using namespace GasUtils;

class gasElement
{
 public:
  gasElement();
  gasElement(const double zz, const double aa);
  gasElement(const gasElement &obj);

  gasElement operator+(const gasElement & obj) const;
  gasElement operator*(const int mul) const;
  gasElement operator*(const double mul) const;

  double getZ() const{return fZ;}
  double getA() const{return fAA;}
  TString getName() const{return fName;}
  double getH() const {return fNElement[kHydrogen];}//free proton
  double getB() const {return fZ-getH();}//bound proton
  double getNElement(const int typ) const{return fNElement[typ];}

  void setName(const TString nn){fName=nn;}
  void Print(const TString info="") const;

 protected:
  void iniElements();
  void addElement(const int zz, const int aa);

  double fZ;
  double fAA;
  TString fName;
  double fNElement[100];//need to check default constructor
};

#endif

/*
//for testing
gSystem->Load("libgasElement.so");
gasElement gC(6, 12), gH(1, 1);
gH.Print();
gC.Print();
gasElement gCH=gC+gH;
gCH.Print();
gasElement gCH4=gC+gH*4;
gCH4.Print();
gasElement gAr(18, 40);
gAr.Print();
gasElement gP10 = gAr*0.9 + gCH4*0.1;
gP10.Print();

 */
