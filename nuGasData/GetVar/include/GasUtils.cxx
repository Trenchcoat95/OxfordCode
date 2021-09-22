#include "GasUtils.h"
#include "TFile.h"
#include "TTree.h"

TString GasUtils::GetElementName(const int type)
{
  if(type==kHydrogen){
    return "H";
  }
  else if(type==kHelium){
    return "He";
  }
  else if(type==kCarbon){
    return "C";
  }
  else if(type==kNitrogen){
    return "N";
  }
  else if(type==kOxygen){
    return "O";
  }
  else if(type==kFluorine){
    return "F";
  }
  else if(type==kArgon){
    return "Ar";
  }
  else{
    printf("gasElement::getElementName unknown type %d\n", type); exit(1);
  }

  return "";
}

double GasUtils::GetElementX0(const int eleType)
{
  //in g cm-2
  switch(eleType){
  case kHydrogen:
    return 63.05;
  case kHelium:
    return 94.32;
  case kCarbon:
    return 42.70;
  case kArgon:
    return 19.55;
  default:
    printf("GasUtils::GetElementX0 unknown element %d\n", eleType); 
    exit(1);
  }
}


double GasUtils::GetElementA(const int eleType)
{
  //in g cm-2
  switch(eleType){
  case kHydrogen:
    return 1;
  case kHelium:
    return 4;
  case kCarbon:
    return 12;
  case kArgon:
    return 40;
  default:
    printf("GasUtils::GetElementA unknown element %d\n", eleType); 
    exit(1);
  }
}

TString GasUtils::ALICEGas()
{
  return "90% Ne + 10% CO_{2}";
}

//static TString ALICEGasEP(){return ALICEGas()+" @ 400 V/cm, 1 atm";}
TString GasUtils::ALICEGasEP()
{
  return Form("#splitline{%s}{@ 400 V/cm, 1 atm}", ALICEGas().Data());
}

TString GasUtils::T2KGas()
{
  return "95% Ar + 3% CF_{4} + 2% iC_{4}H_{10}";
}

//static TString T2KGasEP(){return T2KGas()+" @ 275 V/cm, 1 atm";}

TString GasUtils::T2KGasEP()
{
  return Form("#splitline{%s}{@ 275 V/cm, 1 atm}", T2KGas().Data());
}

double GasUtils::molUnit()
{
  //return 1;//1E5;
  return 1E3;
}

double GasUtils::nevtUnit()
{
  return 1;//1E3;
}

double GasUtils::NA()
{
  //https://en.wikipedia.org/wiki/Mole_(unit)
  //6.02214076×1023
  return 6.022E23;// 1/mol
}

double GasUtils::RESXsec()
{
  //https://journals.aps.org/rmp/pdf/10.1103/RevModPhys.84.1307
  //neutrino cross section:
  //return 0.5E-42; //m2 at 1 GeV

  //antineutrino one is about 1/3 of it

  //take nearest 10s
  return 0.1E-42; //m2 at 1 GeV
}

double GasUtils::DUNEPOT()
{
  //Report from the DUNE Near Detector Concept Study Group
  //May 23, 2018
  //Assuming the ocial 3-horn conguration and 1.5E21 pot/year, with
  //return 1.5e21;//POT/year
  //take nearest 10s
  return 1e21;//POT/year
}


double GasUtils::DUNEFlux()
{
  /*
  https://home.fnal.gov/~ljf26/DUNEFluxes/
  https://home.fnal.gov/~ljf26/DUNEFluxes/OptimizedEngineeredNov2017/

  The histograms have units of nus / m^2 / POT. POT estimates assume 1.1e21 POT/year at 120 GeV.

  xlu@ubuntu:~/Desktop$ r histos_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_neutrino_LBNEND_fastmc.root 
  root [0] 
  Attaching file histos_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_neutrino_LBNEND_fastmc.root as _file0...
  root [1] numu_flux->Integral(9,100000)
  (const Double_t)9.28911881599179764e-04
  root [2] numu_flux->GetBinCenter(9)
  (const Double_t)1.06250000000000000e+00
  root [3] numu_flux->Scale(1,"width")
  root [4] numu_flux->Integral(9,100000,"width")
  (const Double_t)9.28911881599179764e-04
  */

  //ony above 1 GeV
  //return 9.2891e-04;//nus / m^2 / POT
  //take the nearest 10s
  return 10e-04;//nus / m^2 / POT
}


double GasUtils::GetTPCVolume()
{
  //At 10 bars, the active volume of the cylindrical argon-gas TPC in a new-build dipole magnet can be as big as 5.2 m in diameter and 5 m in length, giving an active mass of about 1.8 t. 
 //in m
  //const double rr = 2.5;//m
  const double rr = 2.6;//m
  const double ll = 5;//m
  const double vol = TMath::Pi()*rr*rr*ll;//m3

  return vol;
}

double GasUtils::IdealGasMolarVolume1Bar()
{
  //The molar volume of an ideal gas at 1 atmosphere of pressure is 0.024465 m3/mol at 25 °C. https://en.wikipedia.org/wiki/Molar_volume
  //const double molv = 0.024465;//m3/mol

  //https://en.wikipedia.org/wiki/Molar_volume
  //The molar volume of an ideal gas at 100 kPa (1 bar) is
  //0.022710980(38) m3/mol at 0 °C,
  //0.024789598(42) m3/mol at 25 °C.
  const double molv = 0.024789598;//m3/mol at 1 bar 25c

  return molv;
}

double GasUtils::GetTPCMole(const double pres)
{
  //pres in Bar
  return GetTPCVolume()/IdealGasMolarVolume1Bar()*pres;
}

double GasUtils::GetDUNERESEventRate(const double nmol)
{
  //k event per year
  return nmol*molUnit()*NA()*RESXsec()*DUNEFlux()*DUNEPOT()/nevtUnit();
}

double GasUtils::HmolToPSton(const double hmol)
{
  //http://www.polymerprocessing.com/polymers/PS.html
  //repeat unit C8H8
  //Molecular weight of repeat unit: 104.1 g/mol

  const double molm = 104.1;//g/mol

  const double chton = (hmol*molUnit()/8)*molm*1E-6;

  return chton;
}

double GasUtils::Get3DSTHs()
{
  const double vol = 200*200*200;//cm^3

  //Polystyrene ([C6H5CHCH2]n) 0.53768 57.5 81.7 43.79 1.936 1.06 PDG
  //density 1.06 {g cm−3}

  const double density = 1.06; //g/cm3

  const double mass = density * vol; //g

  //http://www.polymerprocessing.com/polymers/PS.html
  //repeat unit C8H8
  //Molecular weight of repeat unit: 104.1 g/mol

  const double molm = 104.1;//g/mol

  const double nmol = mass/molm;

  return 8*nmol/molUnit();
}


double GasUtils::PressureLimit(const TString gasname)
{
  //https://docs.google.com/spreadsheets/d/15sbDTghJfsXEjltO0P75jo5FVywGq2kE8NCiohpEOsk/edit?usp=sharing
  //for supercritical state, check density at 1 bar and 10 bar, compare density/A/bar. If similar around 4E-2 kg/m^3, then still gas
  //density/A/bar at 1bar  25C (unit 1E-2 kg/m^3) := D1
  //density/A/bar at 10bar 25C (unit 1E-2 kg/m^3) := D10
  //D1=D10 -> still gas at 10bar
  const double SCpressure = 10.001;

  const double verylarge = 1e10;

  if(gasname=="H2"){
    //https://www.wolframalpha.com/input/?i=vapour+pressure+of+hydrogen+at+25C
    //supercritical
    //D1=4.06, D2=4.04
    return SCpressure;
  }
  else if(gasname=="CH4"){
    //at 25C, temperature above cirtical point, pressure still below critical point
    //https://www.wolframalpha.com/input/?i=CH4+25C+10+atm
    //D1=4.05, D2=4.12
    return SCpressure;
  }
  else if(gasname=="C2H6"){
    //at 25C, pressure still below phase boundary
    //https://www.wolframalpha.com/input/?i=C2H6+25C+10+atm
    //https://www.wolframalpha.com/input/?i=vapour+pressure+of+C2H6+at+25C
    return 41.9;
  }
  else if(gasname=="C3H8"){
    //10 atm, 25C, f=93%, 9.3 atm gas, 9.4 atm liquid
    //9.5bar 
    //https://www.wolframalpha.com/input/?i=vapour+pressure+of+C3H8+at+25C
    //9.399 atm (atmospheres) = 9.5235367 bar
    return 9.5;//bar
  }
  else if(gasname=="C4H10"){
    //10atm, 25C, f = 34%, partial pressure 3.4 atm, gas; 3.5 atm liquid
    //https://www.wolframalpha.com/input/?i=vapour+pressure+of+isobutane+at+25C
    //3.461 atm (atmospheres) = 3.5068583 bar
    //isobutane
    return 3.5;//bar
  }
  else if(gasname=="He"){
    //D1 4.04, D2 4.02
    return SCpressure;
  }
  else if(gasname=="Ar"){
    //D1 4.03, D2 4.05
    return SCpressure;
  }
  else if(gasname=="CH2F2"){
    //https://www.wolframalpha.com/input/?i=vapour+pressure+of+CH2F2+at+25C
    return 17.31;//bar
  }
  else if(gasname=="CHF3"){
    //https://www.wolframalpha.com/input/?i=vapour+pressure+of+CHF3+at+25C
    return 45.98;//bar
  }
  else if(gasname=="C2H2F4"){
    //C2H2F4 1,1,1,2-Tetrafluoroethane (R134a), Norflurane
    //https://webbook.nist.gov/cgi/fluid.cgi?TLow=0&THigh=40&TInc=0.1&Applet=on&Digits=5&ID=C811972&Action=Load&Type=SatP&TUnit=C&PUnit=bar&DUnit=g%2Fml&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF
    //5.7 bar at 10C
    //6.65 bar at 25C
    return 6;//bar
  }
  else if(gasname=="CH3OCH3" || gasname=="DME"){//DME
    //https://www.wolframalpha.com/input/?i=vapour+pressure+of+Dimethyl+ether+at+25C
    return 5.90;//bar
  }
  else if(gasname=="LAr" || gasname=="GAr" || gasname == "CH"){//these are from coredraw for range
    return verylarge;
  }

  cout<<"Unknown gasname! "<<gasname<<endl; 
  exit(1);
}

TLine * GasUtils::ZeroYLine(const double xmin, const double xmax, const int lc)
{
  const double zeroy = 0;
  TLine *lnY = new TLine(xmin, zeroy, xmax, zeroy);
  lnY->SetLineWidth(1);
  lnY->SetLineStyle(kDashed);
  lnY->SetLineColor(lc);

  return lnY;
}

TLine * GasUtils::PercentLine(const double xmin, const double xmax, const int lc)
{
  const double limit_dy = 1;
  TLine *lnY = new TLine(xmin, limit_dy, xmax, limit_dy);
  lnY->SetLineWidth(1);
  lnY->SetLineStyle(kDashDotted);
  lnY->SetLineColor(lc);

  return lnY;
}

TLine * GasUtils::CHLine(const double xmin, const double xmax, const int lc)
{
  const double CHy = 1./6.;
  TLine *lnY = new TLine(xmin, CHy, xmax, CHy);
  lnY->SetLineWidth(2);
  lnY->SetLineStyle(kSolid);
  lnY->SetLineColor(lc);
  return lnY;
}

TLine * GasUtils::ALICEX0Line(const double xmin, const double xmax, const int lc)
{
  // is actually A/X0
  const double ALICEy = 0.7437524851; //90Ne+10CO2	20	10	0.1	1	12	6	2	16	8	28.93	42.7	34.24	18	1.2	3.2	0.8035714286	0.05357142857	0.1428571429	30.11754643	0.7437524851	
  TLine *lnY = new TLine(xmin, ALICEy, xmax, ALICEy);
  lnY->SetLineWidth(1);
  lnY->SetLineStyle(kSolid);
  lnY->SetLineColor(lc);
  return lnY;
}

TLine * GasUtils::ALICEDvLine(const double xmin, const double xmax, const int lc)      
{
  const double ALICEy = 2.83;//cm/us
  TLine *lnY = new TLine(xmin, ALICEy, xmax, ALICEy);
  lnY->SetLineWidth(1);
  lnY->SetLineStyle(kSolid);
  lnY->SetLineColor(lc);
  return lnY;
}

TLine * GasUtils::ALICEDiffusionLine(const double xmin, const double xmax, const int lc)
{
  // dt = dl
  const double ALICEy = 220;//um/sqrt(cm)
  TLine *lnY = new TLine(xmin, ALICEy, xmax, ALICEy);
  lnY->SetLineWidth(1);
  lnY->SetLineStyle(kSolid);
  lnY->SetLineColor(lc);

  return lnY;
}

TGraphErrors * GasUtils::ALICETownsendCurve()
{
  TFile* fin = new TFile("../gasSim/ALICE-gain.root", "READ");
  TTree* gastree = (TTree*) fin->Get("gastree");
  if(!gastree){
    printf("no tree!\n"); exit(1);
  }

  float alpharp, eta, efield, rpenning, pressure;
  vector<double> vx, vy;

  gastree->SetBranchAddress("alpharp", &alpharp);
  gastree->SetBranchAddress("eta", &eta);
  gastree->SetBranchAddress("E", &efield);
  gastree->SetBranchAddress("rpenning", &rpenning);
  gastree->SetBranchAddress("p", &pressure);

  for (int i=0; i<gastree->GetEntries(); i++) {
    gastree->GetEntry(i);

    if (fabs(pressure-1e3)>1 || rpenning<55 || rpenning>57.5)
      continue; 

    vx.push_back(efield*1e-3); // axis is in kV/cm
    // std::cerr<<efield<<"\t";
    vy.push_back(alpharp-eta);
  }
  // std::cerr << std::endl;
  
  TGraphErrors* gr = new TGraphErrors(vx.size(), &vx[0], &vy[0]);
  gr->Sort();
  gr->SetLineColor(kRed);
  gr->SetLineWidth(2);
  
  fin->Close();

  return gr;
}

TLine * GasUtils::T2KX0Line(const double xmin, const double xmax, const int lc)
{
  // ar (gas) http://pdg.lbl.gov/2009/AtomicNuclearProperties/HTML_PAGES/018.html -> 19.55 g cm-2
  // cf4 http://pdg.lbl.gov/2009/AtomicNuclearProperties/HTML_PAGES/326.html      -> 33.99 g cm-2
  // butane http://pdg.lbl.gov/2009/AtomicNuclearProperties/HTML_PAGES/124.html   -> 45.23 g cm-2

  const double T2Ky =  (0.95 * 40 +  0.03 * (12 + 4*19) + 0.02 * (4*12 + 10*1)) * ( 0.95 / 19.55 + 0.03 / 33.99 + 0.02 / 45.23); // A/X0 in cm2/g 
  TLine *lnY = new TLine(xmin, T2Ky, xmax, T2Ky);
  lnY->SetLineWidth(1);
  lnY->SetLineStyle(kSolid);
  lnY->SetLineColor(lc);
  return lnY;
}

TLine * GasUtils::T2KDvLine(const double xmin, const double xmax, const int lc)      
{
  const double T2Ky = 7.8;  // cm/us
  TLine *lnY = new TLine(xmin, T2Ky, xmax, T2Ky);
  lnY->SetLineWidth(1);
  lnY->SetLineStyle(kSolid);
  lnY->SetLineColor(lc);
  return lnY;
}

TLine * GasUtils::T2KDiffusionLine(const double xmin, const double xmax, const int lc)
{
  const double ALICEy = 265;// um/sqrt(cm) - transverse
  TLine *lnY = new TLine(xmin, ALICEy, xmax, ALICEy);
  lnY->SetLineWidth(1);
  lnY->SetLineStyle(kSolid);
  lnY->SetLineColor(lc);

  return lnY;
}

TLine * GasUtils::ThermalDiffusionLimit(const double xmin, const double xmax, const int lc, const double E)
{
  // see [Eq. 41] in arXiv:1710.01018
  // or Blum-Rolandi [Eq. 2.64]
  const double kB = TMath::K();   // J/K
  const double qe = TMath::Qe();  // C
  const double T = 298;           // K
  const double limit = std::sqrt( (2.*kB/qe) * T * (1./E) ) * 1e6 / std::sqrt(1e2);  // m/sqrt(m) -> um/sqrt(cm)

  TLine *lnY = new TLine(xmin, limit, xmax, limit);
  lnY->SetLineWidth(1);
  lnY->SetLineStyle(kDashDotted);
  lnY->SetLineColor(lc);

  return lnY;
}
