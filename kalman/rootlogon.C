{ // rootlogon.C
  //gInterpreter->AddIncludePath("/opt/root/6.22.02-install/include/"); // where "DDRec/Vector3D.h" resides
  //gInterpreter->GenerateDictionary("std::vector<ROOT::Math::Vector3D>", "vector;Math/Vector3D.h");
  gInterpreter->GenerateDictionary("std::vector<TVectorT<double>>", "vector;TVectorT.h");
  gInterpreter->GenerateDictionary("std::vector<TMatrixT<double>>", "vector;TMatrixT.h");
  // gSystem->Load("/path/to/libMyClasses");
  // gInterpreter->AddIncludePath("/path/to/MyClasses/include");
  //gROOT->ProcessLine(".L AutoDict_std__vector_ROOT__TVectorF_.cxx+");
}