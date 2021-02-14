#include "TROOT.h"
#include "TRandom.h"
#include "iostream"
#include "TMath.h"

void makerandommuons(float plow, float phigh, int nmuons=1000)
{
  for (int i=0; i<nmuons; ++i)
    {
      float p = plow + gRandom->Uniform(phigh-plow);
      float x = (-2.0 + gRandom->Uniform(4.0))*100.0;
      float y = (-2.0 + gRandom->Uniform(2.0))*100.0;
      float z = -190;
      float px = 0;
      float py = 0;
      float pz = p;
      float m = 0.10566;
      float e = TMath::Sqrt(p*p + m*m);
      std::cout << i << " " <<  1 << std::endl;
      std::cout << "1 13 0 0 0 0 " << px << " " << py << " " << pz << " " << e << " "  << m << " " << x << " " << y << " " << z << " " << 0 << std::endl; 
    }
}
