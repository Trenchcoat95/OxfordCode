#ifndef _GASMIXER_
#define _GASMIXER_

#include <iostream>
#include "TMath.h"
#include "TString.h"
#include "gasElement.h"
#include "GasUtils.h"

using namespace std;
using namespace GasUtils;

class gasMixer: public gasElement
{
public:
 gasMixer():gasElement(){}
 gasMixer(const double zz, const double aa): gasElement(zz, aa){}
 gasMixer(const int zz, const int pp, const int qq, const double xx);
 gasMixer(const gasElement ge):gasElement(ge){}//needed for gasMixer"="gasElement
  //gasMixer(const gasMixer g1, const double f1, const gasMixer g2, const double f2, const double pres, const bool kSimpleName=false);
  gasMixer(const gasMixer gg,  const double fracg, const double pres, const gasMixer *ptrPatcher=0x0);
  gasMixer(const gasMixer gbase, const gasMixer gtop,  const double fractop, const double pres, const gasMixer *basePatcherPtr=0x0, const gasMixer *topPatcherPtr=0x0);


  double getTPCHMole() const;
  double getProtonFree2Bound() const;
  double getHindex() const;
  double getTPCArMass() const;
  double getTPCArMassFraction() const;

  double getX0() const;
  //double getTheta0(const double pressure) const;
  double getMS() const;

  void Print() const;

  TString getBase(TString gasname="") const;
  TString getTop(TString gasname="") const;

  static gasMixer GetCH(const int nC, const int nH);//for CH
  static gasMixer GetCHF(const int nC, const int nH, const int nF);//for CHF

 private:
  /*
  //gas mixture (1-x)(Z, A) + x C_p H_q
  int fZ;
  int fP;
  int fQ;
  double fX;
  */

  //double getGasDensity(const double p) const;
  double getAOverX0() const;

  //double getHToH2() const;
  //double getNonHToH2() const;

  //double getPurityToCH() const;  
  //double get10BarTPCto3DSTH() const;
};

#endif
