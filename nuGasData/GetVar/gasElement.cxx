#include "gasElement.h"


gasElement::gasElement(): fZ(0), fAA(0), fName("")
{ 
  iniElements(); 
}

gasElement::gasElement(const double zz, const double aa): fZ(zz), fAA(aa), fName("") 
{ 
  if(zz>aa){
    printf("gasElement::gasElement zz>aa %f %f\n", zz, aa); exit(1);
  }

  iniElements();
  //allow zz=0 for null element in e.g. [100]Ar[0]CH4
  if(zz>1E-10){
    addElement(zz, aa); 
  }
}

gasElement::gasElement(const gasElement &obj)
{
  fAA = obj.fAA;
  fZ = obj.fZ;
  fName = obj.fName;
  for(unsigned int ii=0; ii<sizeof(fNElement)/sizeof(double); ii++){
    fNElement[ii] = obj.fNElement[ii];
  }
}

void gasElement::iniElements()
{
  const int nsize = sizeof(fNElement)/sizeof(double);
  for(unsigned int ii=0; ii<nsize; ii++){
    fNElement[ii]=0;
  }
}

void gasElement::addElement(const int zz, const int aa)
{
  ElementType kelmt = kUnknown;
  if(zz==1 && aa==1){
    kelmt = kHydrogen;
  }
  else if(zz==2 && aa==4){
    kelmt = kHelium;
  }
  else if(zz==6 && aa==12){
    kelmt = kCarbon;
  }
  else if(zz==7 && aa==14){
    kelmt = kNitrogen;
  }
  else if(zz==8 && aa==16){
    kelmt = kOxygen;
  }
  else  if(zz==9 && aa==19){
    kelmt = kFluorine;
  }
  else if(zz==18 && aa==40){
    kelmt = kArgon;
  }
  else{
    printf("gasElement::addElement element unknown %d %d\n", zz, aa); exit(1);
  }

  fName = GetElementName(kelmt);
  fNElement[kelmt]++;
}

gasElement gasElement::operator+(const gasElement & obj) const
{
  gasElement outp(*this);
  outp.fZ += obj.fZ;
  outp.fAA += obj.fAA;
  outp.fName += obj.fName;
  for(unsigned int ii=0; ii<kNElement; ii++){
    outp.fNElement[ii] += obj.fNElement[ii];
  }
  return outp;
}

gasElement gasElement::operator*(const double mul) const
{
  gasElement outp(*this);
  outp.fZ *= mul;
  outp.fAA *= mul;
  for(unsigned int ii=0; ii<kNElement; ii++){
    outp.fNElement[ii] *= mul;
  }
  outp.fName = Form("[%.0f]%s", mul*100, fName.Data());

  return outp;
}

gasElement gasElement::operator*(const int mul) const
{
  gasElement outp = (*this)*((double)mul);
  outp.fName = Form("%s%d", fName.Data(), mul);
  return outp;
}

void gasElement::Print(const TString info) const 
{
  printf("\ngasElement::Print Name %20s\tZ %10.2f\tA %10.2f\tH %10.2f\tB %10.2f", fName.Data(), fZ, fAA, getH(), getB());
  for(int ii=0; ii<kNElement; ii++){
    if(fNElement[ii]>1E-10){
      printf(",\tElement %2s%2d\t%.2f", GetElementName(ii).Data(), ii, fNElement[ii]);
    }
  }
  if(info.Length()){
    printf("\t%s", info.Data());
  }
  printf("\n\n");
}
