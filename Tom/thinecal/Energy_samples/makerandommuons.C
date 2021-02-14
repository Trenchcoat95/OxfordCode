#include "TROOT.h"
#include "TRandom.h"
#include "iostream"
#include "TMath.h"

void makerandommuons(float plow, float phigh, int nmuons=1000)
{
  for (int i=0; i<nmuons; ++i)
    {
      float p = plow + gRandom->Uniform(phigh-plow);
      float x = 0;
      float y = 0;
      float z = -500;
      float px = 0;
      float py = 0;
      float pz = p;
      float m = 0.10566;
      float e = TMath::Sqrt(p*p + m*m);
      std::cout << i << " " <<  1 << std::endl;
      std::cout << "1 13 0 0 0 0 " << px << " " << py << " " << pz << " " << e << " "  << m << " " << x << " " << y << " " << z << " " << 0 << std::endl; 
    }
}
