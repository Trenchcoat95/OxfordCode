#include "TROOT.h"
#include "TRandom.h"
#include "iostream"
#include "TMath.h"

void makerandommuons(float alow, float ahigh, int nmuons=1000)
{
  for (int i=0; i<nmuons; ++i)
    {
      float x = 0;
      float y = 0;
      float z = -500;
      float p = 3;
      float degangle = alow + gRandom->Uniform(ahigh-alow);
      float angle = degangle * 3.14159265 / 180 ;
      float px = p*sin(angle);
      float py = 0;
      float pz = p*cos(angle);
      float m = 0.10566;
      float e = TMath::Sqrt(p*p + m*m);
      std::cout << i << " " <<  1 << std::endl;
      std::cout << "1 13 0 0 0 0 " << px << " " << py << " " << pz << " " << e << " "  << m << " " << x << " " << y << " " << z << " " << 0 << std::endl; 
    }
}
