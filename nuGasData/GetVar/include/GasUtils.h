#ifndef _GASUTILS_
#define _GASUTILS_

#include <iostream>
#include "TLine.h"
#include "TMath.h"
#include "TString.h"
#include "TGraphErrors.h"

using namespace std;

namespace GasUtils
{
  enum ElementType{
    kHydrogen=0,
    kHelium,
    kCarbon,
    kNitrogen,
    kOxygen,
    kFluorine,
    kArgon,
    kNElement,
    kUnknown
  };

  double molUnit();
  double nevtUnit();
  double NA();
  double RESXsec();
  double PressureLimit(const TString gasname);
  double GetElementX0(const int eleType);
  double GetElementA(const int eleType);
  TString GetElementName(const int type);

  double DUNEPOT();
  double DUNEFlux();

  double IdealGasMolarVolume1Bar();
  double GetTPCVolume();
  double GetTPCMole(const double pres);
  double GetDUNERESEventRate(const double nmol);
  double HmolToPSton(const double hmol);
  double Get3DSTHs();

  TString ALICEGas();
  TString ALICEGasEP();
  TString T2KGas();
  TString T2KGasEP();

  TLine * ZeroYLine(const double xmin, const double xmax, const int lc);
  TLine * PercentLine(const double xmin, const double xmax, const int lc);
  TLine * CHLine(const double xmin, const double xmax, const int lc);
  TLine * ALICEX0Line(const double xmin, const double xmax, const int lc);
  TLine * ALICEDvLine(const double xmin, const double xmax, const int lc);
  TLine * ALICEDiffusionLine(const double xmin, const double xmax, const int lc);
  TGraphErrors * ALICETownsendCurve();
  TLine * T2KX0Line(const double xmin, const double xmax, const int lc);
  TLine * T2KDvLine(const double xmin, const double xmax, const int lc);
  TLine * T2KDiffusionLine(const double xmin, const double xmax, const int lc);
//  TGraphErrors * T2KTownsendCurve();  // has not correct rpenning
  TLine * ThermalDiffusionLimit(const double xmin, const double xmax, const int lc, const double E);

}

#endif
