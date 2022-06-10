#include "TROOT.h"
#include "TRandom.h"
#include "iostream"
#include "TMath.h"

void makerandommuons(float plow, float phigh, int nmuons=1000)
{
  //freopen("output.txt","w",stdout);
  gRandom->SetSeed(0);
  gSystem->RedirectOutput("muons.txt","w");

  double  Plane_XY_fid[4] = {-200, +200, -300, 0};
  double  Planes_Z[] =   {1179,1183,1199,1203,1219,1223,1239,1243,1339,1343,1539,1543};

  for (int i=0; i<nmuons; ++i)
    {
      float p = plow + gRandom->Rndm()*(phigh-plow);
      double px=-0.5+gRandom->Rndm();
      double py=-0.5+gRandom->Rndm();
      double pz=1.0+gRandom->Rndm();
      double norm=p/sqrt(px*px+py*py+pz*pz);
      px*=norm;
      py*=norm;
      pz*=norm;
      double pdg = (gRandom->Rndm()<0.5) ? -13 : 13;
      float x = gRandom->Rndm()*(Plane_XY_fid[1]-Plane_XY_fid[0])+Plane_XY_fid[0];
      float y = gRandom->Rndm()*(Plane_XY_fid[3]-Plane_XY_fid[2])+Plane_XY_fid[2];
      float z = gRandom->Rndm()*10+(Planes_Z[0]-10);
  

      
      float m = 0.10566;
      float e = TMath::Sqrt(p*p + m*m);
      std::cout << i << " " <<  1 << std::endl;
      std::cout << "1 "<<pdg<<" 0 0 0 0 " << px << " " << py << " " << pz << " " << e << " "  << m << " " << x << " " << y << " " << z << " " << 0 << std::endl; 
    }
  gApplication->Terminate(0);
  
}
