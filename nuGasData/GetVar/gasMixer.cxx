#include "gasMixer.h"

gasMixer gasMixer::GetCH(const int nC, const int nH)
{
  const gasMixer gC(6, 12), gH(1,1);
  return  gC*nC + gH*nH;
}

gasMixer gasMixer::GetCHF(const int nC, const int nH, const int nF)
{
  const gasMixer gF(9,19);

  return GetCH(nC, nH) + gF*nF;
}

gasMixer::gasMixer(const gasMixer gg,  const double fracg, const double pres, const gasMixer *ptrPatcher)
{
  gasMixer gCook;

  const double maxfg  = PressureLimit(gg.getName())/pres;

  if(fracg<maxfg){
    gCook = gg*fracg;
  }
  else if(ptrPatcher){//only correct over-pressure if there is a patcher
    gCook = gg*maxfg + (*ptrPatcher)*(fracg-maxfg);
    
    printf("\nbefore patching: "); 
    gCook.Print();
    
    gCook.setName(Form("[%.0f]%s", fracg*100, gg.getName().Data()));
    
    printf("\nafter patching: "); 
    gCook.Print();
  }

  (*this) = gCook;
}


gasMixer::gasMixer(const gasMixer gbase, const gasMixer gtop,  const double fractop, const double pres, const gasMixer *basePatcherPtr, const gasMixer *topPatcherPtr)
{
  const gasMixer gCookBase(gbase, 1-fractop, pres, basePatcherPtr);
  const gasMixer gCookTop(gtop, fractop, pres, topPatcherPtr);

  gasMixer gc;
  if(gCookBase.getName()!=""&&gCookTop.getName()!=""){
    gc = gCookBase + gCookTop;
  }
  (*this) = gc;
}

/*
gasMixer::gasMixer(const gasMixer g1, const double f1, const gasMixer g2, const double f2, const double pres, const bool kSimpleName)
{
  //if still gas
  //return g1*f1+g2*f2
  //otherwise return empty mix with name ""

  const double p1max = PressureLimit(g1.getName());
  const double p2max = PressureLimit(g2.getName());

  gasMixer tmp;

  if( (pres*f1)<p1max && (pres*f2)<p2max ){
    tmp = g1*f1+g2*f2;

    if(kSimpleName){
      tmp.setName(g1.getName()+g2.getName());
    }
  }

  (*this) = tmp;
}
*/

gasMixer::gasMixer(const int zz, const int pp, const int qq, const double xx)
{
  int tmpAA = zz*2;
  if(zz==18){
    tmpAA = 40;
  }
  //gas mixture (1-x)(Z, A) + x C_p H_q
  const gasElement gH(  1,  1);
  const gasElement gC(  6, 12);
  const gasElement gAZ( zz, tmpAA);

  const gasElement gCH=gC*pp + gH*qq;
  const gasElement gZCH=gAZ*(1-xx) + gCH*xx;

  (*this) = gZCH;
}

double gasMixer::getTPCHMole() const
{
  const double pres = 10;  //at 10 bar
  return getH()*GetTPCMole(pres)/molUnit();//mol
}

double gasMixer::getTPCArMass() const
{
  const double pres = 10;  //at 10 bar
  const double nmol = getNElement(kArgon)*GetTPCMole(pres);//mol

  //https://www.wolframalpha.com/input/?i=argon+molar+mass
  //39.948 g/mol (grams per mole)

  //const double molmass = 39.948;// g/mol;

  return nmol * 40/1E6;//in ton
}


double gasMixer::getTPCArMassFraction() const
{
  //in %
  const double pres = 10;  //at 10 bar
  const double ALLmass = getA()*GetTPCMole(pres)/1E6;//in t
  const double amass = getTPCArMass();//in t

  return amass/ALLmass*100;
}

/*

double gasMixer::get10BarTPCto3DSTH() const
{
  return getTPCHMole()/Get3DSTHs();
}

double gasMixer::getHToH2() const
{
  return getH()/2;
}

double gasMixer::getNonHToH2() const
{
  return getB() / 2;
}

double gasMixer::getPurityToCH() const
{
  const double CHpur = 0.166667;
  return getHToH2()/(getNonHToH2()+1E-10)/CHpur;
}
*/

double gasMixer::getProtonFree2Bound() const
{
  return getH()/(getB()+1E-10);
}

double gasMixer::getHindex() const
{
  return getProtonFree2Bound() * getH();
}

double gasMixer::getMS() const
{
  return getAOverX0();
}

double gasMixer::getAOverX0() const
{
  //in g cm-2
  double sumAoverX0 = 0;
  for(int ii=0; ii< kNElement; ii++){
    const double nele = getNElement(ii);
    if(nele==0){
      continue;
    }
    const double tmpx0 = GetElementX0(ii);
    const double tmpA = GetElementA(ii);
    const double toadd = nele*tmpA/tmpx0;
    sumAoverX0 += toadd;
    //printf("test %d %s %f %f %f\n", ii, GetElementName(ii).Data(), tmpx0, tmpA, toadd);
  }

  return sumAoverX0;
}

double gasMixer::getX0() const
{
  return getA()/getAOverX0();
}

/*
double gasMixer::getGasDensity(const double p) const
{
  // returns an estimated density in g/cm3 for gas mix at given pressure p in mbar
    
  // constants in g / cm3 - from wolframalpha at STP := 1000mbar, 298K
  double noble_stp;
  double quencher_stp;

  switch (this->fZ)
  {
  case 18:  // Ar
    noble_stp = 16.13e-4;
    break;
  
  case 2:  // He
    noble_stp = 1.615e-4;
  
  default:
    std::cout << "Unknown noble gas requested in gasMixer::getGasDensity - aborting" << std::endl;
    exit(1);
    break;
  }

  switch (this->fP)
  {
  case 1:  // CH4
    quencher_stp = 6.486e-4;
    break;
  case 2:  // C2H6
    quencher_stp = 12.23e-4;
    break;
  case 3:  // C3H8
    quencher_stp = 18.09e-4;  // could be replaced by a lambda that scales linearly between stp and 10bar
    break;

  default:
    std::cout << "Unknown quenching gas requested in gasMixer::getGasDensity - aborting" << std::endl;
    exit(1);
    break;
  }

  // scale using ideal gas law
  double density = (1-this->fX) * noble_stp + this->fX * quencher_stp;
  density *= p * 1e-3;

  return density;
}

double gasMixer::getTheta0(const double pressure) const
{
  // see PDG chapter 33.3 'Multiple scattering through small angles'
  // using eq 33.15 - p451 in PDG2018
  // need to agree on an energy an a particle (actually only if electron or not and energy)

  // constants
  const double x = 500;                             // traversed distance in material in cm 
  const double Joule_eV = 1.602176e-19;             // in J/eV
  double density = this->getGasDensity(pressure);   // in g/cm3

  const double x0 = this->getX0() / density;        // in cm
  const double c_vac = 1;                           // in natural units

  // the particle - a proton at 100 MeV
  const double z = 1;           // charge number of incident particle
  const double momentum = 100;  // in MeV
  const double mass = 938;      // in MeV
  const double energy = std::sqrt( mass * mass + momentum * momentum);  // in MeV
  const double beta = momentum / energy;

  // terms of eq 33.15
  const double term_momentum = 13.6e6 * Joule_eV / beta / c_vac / momentum;
  const double term_distance = std::sqrt( x / x0 );
  const double term_expansion = 1 + 0.038 * std::log( x * z * z / x0 / beta / beta);

  return term_momentum * z * term_distance * term_expansion;
}
*/

void gasMixer::Print() const
{
  gasElement::Print(Form("purity: %f H-index: %f", getProtonFree2Bound(), getHindex()));
}

TString gasMixer::getBase(TString gasname) const
{
  //internal name
  if(gasname==""){
    gasname = fName;
  }

  TString gout; 
  if(gasname.Contains("[")){
    gout = gasname(gasname.First("]")+1, gasname.Length());
    gout = gout(0, gout.First("["));
  }

  return gout;
}

TString gasMixer::getTop(TString gasname) const
{
  //internal name
  if(gasname==""){
    gasname = fName;
  }

  TString gout; 
  if(gasname.Contains("[")){
    gout=gasname(gasname.First("]")+1, gasname.Length())+"[";//need a padding [ for gasBase to work
    gout=getBase(gout);
  }
  return gout;
}
