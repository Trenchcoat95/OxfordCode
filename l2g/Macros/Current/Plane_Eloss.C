#include "TGraph2D.h"
#include "TGraph.h"
#include "l2g_trackmatch_t2.C"

//   In a ROOT session, you can do:
  //      root> .L Plane_Eloss.C
  //      root> garana t
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.Plane_Eloss();       // Use the trackmatching function


float BetheBlochGeant(float bg,
         float kp0,
         float kp1,
         float kp2,
         float kp3,
         float kp4) {
  //
  // This is the parameterization of the Bethe-Bloch formula inspired by Geant.
  //
  // bg  - beta*gamma
  // kp0 - density [g/cm^3]
  // kp1 - density effect first junction point
  // kp2 - density effect second junction point
  // kp3 - mean excitation energy [GeV]
  // kp4 - mean Z/A
  
  //
  // The default values for the kp* parameters are for silicon. 
  // The returned value is in [GeV/(g/cm^2)].
  // 

  const float mK  = 0.307075e-3; // [GeV*cm^2/g]
  const float me  = 0.511e-3;    // [GeV/c^2]
  const float rho = kp0;
  const float x0  = kp1*2.303;
  const float x1  = kp2*2.303;
  const float mI  = kp3;
  const float mZA = kp4;
  const float bg2 = bg*bg;
  const float maxT= 2*me*bg2;    // neglecting the electron mass
  
  //*** Density effect
  float d2=0.; 
  const float x=TMath::Log(bg);
  const float lhwI=TMath::Log(28.816*1e-9*TMath::Sqrt(rho*mZA)/mI);
  if (x > x1) {
    d2 = lhwI + x - 0.5;
  } else if (x > x0) {
    const float r=(x1-x)/(x1-x0);
    d2 = lhwI + x - 0.5 + (0.5 - lhwI - x0)*r*r*r;
  }

  return mK*mZA*(1+bg2)/bg2*
         (0.5*TMath::Log(2*me*bg2*maxT/(mI*mI)) - bg2/(1+bg2) - d2);
}

float BetheBlochLeo(float bg,
         float kp0,
         float kp1,
         float kp2,
         float kp3,
         float kp4,
         float kp5,
         float kp6,
         float kp7,
         float kp8,
         float kp9) {
  //
  // This is the parameterization of the Bethe-Bloch formula inspired by Geant.
  //
  // bg  - beta*gamma
  // kp0 - density [g/cm^3]
  // kp1 - density effect first junction point
  // kp2 - density effect second junction point
  // kp3 - mean excitation energy [GeV]
  // kp4 - mean Z/A
  // kp5 - C
  // kp6 - Z
  //
  // The default values for the kp* parameters are for silicon. 
  // The returned value is in [GeV/(g/cm^2)].
  // 

  //const float mK  = 0.1535e-3; // [GeV*cm^2/g]
  
  const float me  = 0.511e-3;    // [GeV/c^2]
  const float mK=0.1535e-3;
  const float rho = kp0;
  const float x0  = kp1;
  const float x1  = kp2;
  const float mI  = kp3;
  const float mZA = kp4;
  const float hw = kp5;
  const float a = kp6;
  const float m = kp7;
  const float Z = kp8;
  const float Mass = kp9;
  const float bg2 = bg*bg;
  const float maxT= 2*me*bg2/(1+2*(me/Mass)*sqrt(1+bg2)+pow((me/Mass),2));    // neglecting the electron mass
  
  //*** Density effect
  float d2=0.; 
  float x=TMath::Log10(bg);
  float C=2*log(mI/hw)+1;

  
  if (x > x1) {
    d2 = 4.6052*x + C;
  } else if (x > x0) {
    d2 = 4.6052*x + C + a*(pow(x1-x,m));
  }
  //std::cout<<"dEdx"<<d2<<std::endl;

  return mK*mZA*(1+bg2)/bg2*rho*
         (log(2*me*bg2*maxT/(mI*mI)) - 2*bg2/(1+bg2) - d2 + 2*(C/Z));
}


void garana::Plane_Eloss()
{
    float ZA =0.54141;
    float I=64.7e-9;      ///GeV
    float rho= 1.032;   //g/cm^3
    float X1=2.49;
    float X0=0.1469;
    float muon_mass=0.1056583755; //GeV/c^2
    float a=0.1610;
    float m=3.24;
    float hw=21.54e-9;
    float Z=0.085+6*0.915;


    float  First_Plane[6] = {-300, +300, -400, 100, 1244, 1248 }; //all in cm
    float  Plane_thick = First_Plane[5]-First_Plane[4];
 
    Int_t nentries = fChain->GetEntries();

    gROOT->SetBatch(); ///Run in batch mode


    #pragma region "Plot Declaration"

        /////////////////////////////////////1D plots of the three main variables and 2D correlations
        TH1F* h_deltap = new TH1F("h_deltap", "Momentum loss in plane", 100, 0.0, 0.04);
        TH1F* h_deltaE = new TH1F("h_deltaE", "Energy loss in plane", 100, 0.0, 0.1);
        TH1F* h_deltaEB = new TH1F("h_deltaEB", "Expected energy loss in plane", 100, 0.0, 0.1);
        TH1F* h_deltaEBGeant = new TH1F("h_deltaEBGeant", "Expected energy loss in plane", 100, 0.0, 0.1);
        TH1F* h_deltaERes = new TH1F("h_deltaERes", "Expected energy loss in plane", 100, -2, 2);
        TH1F* h_deltaEResGeant = new TH1F("h_deltaEResGeant", "Expected energy loss in plane", 100, -2, 2);
        TH1F* h_pcontact = new TH1F("h_pcontact", "Momentum at impact point", 50, 0, 8);
        TH1F* h_pzcontact = new TH1F("h_pzcontact", "Momentum z at impact point", 50, 0, 8);
        TH1F* h_costheta = new TH1F("h_costheta", "cos#theta at impact point with plane", 50, 0, 1);

        TH2F* h_xybef = new TH2F("h_xybef", "cos#theta at impact point with plane", 300, -310, 310, 300, -410,110);
        TH2F* h_zybef = new TH2F("h_zybef", "cos#theta at impact point with plane",300,1100,1400,300, -410,110);
        TH2F* h_xypost = new TH2F("h_xypost", "cos#theta at impact point with plane", 310, -310, 300, 300, -410,110);
        TH2F* h_zypost = new TH2F("h_zypost", "cos#theta at impact point with plane",300,1100,1400,300, -410,110);
        
        TH2F* h_pcontactVScostheta = new TH2F("h_pcontactVScostheta ","pcontactVScostheta",50, 0, 1, 50, 0, 8);
        TH2F* h_deltapVScostheta = new TH2F("h_deltapVScostheta ","deltapVScostheta",50, 0, 1, 50, 0, 0.04);
        TH2F* h_deltapVSpcontact = new TH2F("h_deltapVSpcontact ","deltapVSpcontact",50, 0, 4, 50, 0, 0.04);
        TH2F* h_deltaEVScostheta = new TH2F("h_deltaEVScostheta ","deltapVScostheta",50, 0, 1, 50, 0, 0.1);
        TH2F* h_deltaEexpVScostheta = new TH2F("h_deltaEexpVScostheta ","deltapVScostheta",50, 0, 1, 50, 0, 0.1);
       
        TH2F* P_deltapVScostheta = new TH2F("P_deltapVScostheta ","deltapVScostheta",50, 0, 1, 50, 0, 0.04);
        TH2F* P_deltapVSpcontact = new TH2F("P_deltapVSpcontact ","deltapVSpcontact",50, 0, 4, 50, 0, 0.04);
        TH2F* P_deltaEVScostheta = new TH2F("P_deltaEVScostheta ","deltapVScostheta",50, 0, 1, 50, 0, 0.1);
        TH2F* P_deltaEexpVScostheta = new TH2F("P_deltaEexpVScostheta ","deltapVScostheta",50, 0, 1, 50, 0, 0.1);
       

        TGraph* g_deltapVScostheta = new TGraph;
        TGraph* g_deltapVSpcontact = new TGraph;

        TGraph2D* g_deltapVScosthetaVSpcontact = new TGraph2D;
        
        
    #pragma endregion    

    bool showprog = true;  
    if(showprog==true) std::cout<<"Progress:  "<<std::endl;

    Int_t pp=0;

    for (Int_t i=0; i<nentries; i++) 
    {
        fChain->GetEntry(i);
        int prog = 100*i/nentries;
        std::string strprog = std::to_string(prog);
        if(showprog==true) std::cout<<strprog<<"%";

        //std::cout<<i<<std::endl;
        
        

        if(showprog==true) std::cout << std::string(strprog.length(),'\b')<<"\b";

        //fiducial cut in LAr
        int nTracks = TrackStartX->size();
        //std::cout<<"n of trajectories: "<<PDG->size()<<std::endl;
        //std::cout<<"n of tracks: "<<nTracks<<std::endl;
        
        int InGArcount = 0;
        //Cycle over MC Particle
        for (Int_t j=0; j<PDG->size(); j++) 
        {
           if((PDG->at(j)==13 || PDG->at(j)==-13) && PDGMother->at(j)==0) //consider only primary muons
           {
               TVector3 p_pre;
               TVector3 p_post;
               TVector3 xyz_pre;
               TVector3 xyz_post;
               TVector3 normal(0,0,1);
             
               float deltap = 0;
               float deltaE =0;
               float deltaEBethe =0;
               float deltaEBetheG=0;
               float E_pre=0;
               float E_post=0;
               float Efrac =0;
               float EfracG =0;
               float p_contact = 0;
               float thetaRad = 0;
               float costheta = 0;

               int ID = MCTrkID->at(j);
               //std::cout<<MCTrkID->at(j)<<std::endl;
               
               int BefPlane=0;
               int PostPlane=0;

               //Cycle over all trajectory points, Find the ones corresponding to the current MC Particle and Find first point in GAr
               for (Int_t k=0; k<TrajMCPX->size(); k++) 
                {
                    if (TrajMCPTrackID->at(k)==ID)
                    {
                    
                    if (TrajMCPX->at(k)>(First_Plane[0]+10) && TrajMCPX->at(k)<(First_Plane[1]-10) 
                        && TrajMCPY->at(k)>(First_Plane[2]+10) && TrajMCPY->at(k)<(First_Plane[3]-10)
                        && TrajMCPZ->at(k)<First_Plane[4] && TrajMCPZ->at(k)>(First_Plane[4]-10))
                    {
                        BefPlane++;
                        p_pre.SetXYZ(TrajMCPPX->at(k),TrajMCPPY->at(k),TrajMCPPZ->at(k));
                        xyz_pre.SetXYZ(TrajMCPX->at(k),TrajMCPY->at(k),TrajMCPZ->at(k));
                    }

                    if (TrajMCPX->at(k)>(First_Plane[0]+10) && TrajMCPX->at(k)<(First_Plane[1]-10) 
                        && TrajMCPY->at(k)>(First_Plane[2]+10) && TrajMCPY->at(k)<(First_Plane[3]-10)
                        && TrajMCPZ->at(k)>First_Plane[5]&& TrajMCPZ->at(k)<(First_Plane[5]+10)&&BefPlane>0)
                    {
                        PostPlane++;
                        p_post.SetXYZ(TrajMCPPX->at(k),TrajMCPPY->at(k),TrajMCPPZ->at(k));
                        xyz_post.SetXYZ(TrajMCPX->at(k),TrajMCPY->at(k),TrajMCPZ->at(k));
                        break;
                    }
                    }
                }

                
                
                //std::cout<<"InGAr"<<InGAr<<std::endl;
                //Fill the MC Particle in GAr distributions and Fill the nhits per track in GAr distribution
                if(BefPlane>0 && PostPlane>0 && (p_pre.Mag()-p_post.Mag())>0) 
                {
                    deltap = p_pre.Mag()-p_post.Mag();
                    E_pre =sqrt(p_pre.Mag2()+muon_mass*muon_mass);
                    E_post =sqrt(p_post.Mag2()+muon_mass*muon_mass);
                    deltaE =E_pre-E_post;
                
                    
                    p_contact = p_pre.Mag();
                    thetaRad = p_pre.Angle(normal);
                    costheta = cos(thetaRad);

                    deltaEBethe= BetheBlochLeo((p_pre.Mag()/muon_mass),rho,X0,X1,I,ZA,hw,a,m,Z,muon_mass) * rho * Plane_thick * (1/costheta);
                    deltaEBetheG= BetheBlochGeant((p_pre.Mag()/muon_mass),rho,X0,X1,I,ZA) * rho * Plane_thick * (1/costheta);
                    //std::cout<<deltaEBethe<<std::endl;

                    Efrac=(deltaE-deltaEBethe)/deltaE;
                    EfracG=(deltaE-deltaEBetheG)/deltaE;

                    h_deltaE->Fill(deltaE);
                    h_deltaEB->Fill(deltaEBethe);
                    h_deltaEBGeant->Fill(deltaEBetheG);
                    h_deltaERes->Fill(Efrac);
                    h_deltaEResGeant->Fill(EfracG);

                    h_xybef->Fill(xyz_pre.X(),xyz_pre.Y());
                    h_zybef->Fill(xyz_pre.Z(),xyz_pre.Y());
                    h_xypost->Fill(xyz_post.X(),xyz_post.Y());
                    h_zypost->Fill(xyz_post.Z(),xyz_post.Y());

                    h_deltap->Fill(deltap);
                    h_pcontact->Fill(p_contact);
                    h_costheta->Fill(costheta);
                    h_pzcontact->Fill(p_pre.Z());
                    h_pcontactVScostheta->Fill(costheta,p_contact);
                    h_deltapVScostheta->Fill(costheta,deltap);
                    h_deltapVSpcontact->Fill(p_contact,deltap);

                    h_deltaEVScostheta->Fill(costheta,deltaE*costheta);
                    h_deltaEexpVScostheta->Fill(costheta,deltaEBetheG*costheta);

                    if(costheta<=1&&costheta>=0&&deltap<=0.3&&deltap>=0.0&&p_contact>=0&&p_contact<=8)
                    {
                    g_deltapVScosthetaVSpcontact->SetPoint(pp,p_contact,costheta,deltap);
                    g_deltapVScostheta->SetPoint(pp,costheta,deltap);
                    g_deltapVSpcontact->SetPoint(pp,p_contact,deltap);
                    pp++;
                    }
                }
               
            
            
               
           }
        }
        
    //std::cout<<"InGArCount "<<InGArcount<<std::endl; 
    //std::cout<<"InGArCount "<<htheta_InGAr_Sample->GetEntries()<<std::endl;  

    }

    GetProbabilityPlot(h_deltapVSpcontact,P_deltapVSpcontact);
    GetProbabilityPlot(h_deltapVScostheta,P_deltapVScostheta);

    GetProbabilityPlot(h_deltaEVScostheta,P_deltaEVScostheta);
    GetProbabilityPlot(h_deltaEexpVScostheta,P_deltaEexpVScostheta);
  

    std::string Formula = "abs(0.39894228040143*"+std::to_string(h_deltap->GetBinWidth(0))+"*([0]/([2]*x))*(exp(-0.5*((log(x)-[1])/[2])^2)))";
    TF1 *LogNormal = new TF1("LogNormal",Formula.c_str(),0.0,0.3);
    LogNormal->SetParameters(h_deltap->GetEntries(),log(h_deltap->GetMean()),log(h_deltap->GetRMS()));
    
    #pragma region "Plotting"

        ///////////////////////////////////////////////1D distributions for three main variables
        gStyle->SetOptStat(1);
        gStyle->SetOptFit(1);
        TCanvas *mccanvasdeltap = new TCanvas("mccanvasdeltap","",1000,800);
        h_deltap->SetTitle("Momentum loss in plane;#Delta p [GeV/c];n");
        h_deltap->Fit("LogNormal");
        h_deltap->Draw();
        mccanvasdeltap->Print("E_loss/deltap.png");

        TCanvas *mccanvasdeltapL = new TCanvas("mccanvasdeltapL","",1000,800);
        h_deltap->SetTitle("Momentum loss in plane;#Delta p [GeV/c];n");
        h_deltap->Fit("landau","","",0.003,0.03);
        h_deltap->Draw();
        mccanvasdeltapL->Print("E_loss/deltapLandau.png");

        TCanvas *mccanvasdeltaE = new TCanvas("mccanvasdeltaE","",1000,800);
        h_deltaE->SetTitle("Energy loss in plane;#Delta E [GeV/c^{2}];n");
        h_deltaE->Draw();
        mccanvasdeltaE->Print("E_loss/deltaE.png");

        TCanvas *mccanvasdeltaEB = new TCanvas("mccanvasdeltaEB","",1000,800);
        h_deltaEB->SetTitle("Expected energy loss in plane;#Delta E [GeV/c^{2}];n");
        h_deltaEB->Draw();
        mccanvasdeltaEB->Print("E_loss/deltaEB.png");

        TCanvas *mccanvasdeltaEEB = new TCanvas("mccanvasdeltaEEB","",1000,800);
        h_deltaERes->SetTitle("Energy loss Residual;Res;n");
        h_deltaERes->Draw();
        mccanvasdeltaEEB->Print("E_loss/deltaERes.png");

        TCanvas *mccanvasdeltaEBG = new TCanvas("mccanvasdeltaEBG","",1000,800);
        h_deltaEBGeant->SetTitle("Expected energy loss in plane;#Delta E [GeV/c^{2}];n");
        h_deltaEBGeant->Draw();
        mccanvasdeltaEBG->Print("E_loss/deltaEBG.png");

        TCanvas *mccanvasdeltaEEBG = new TCanvas("mccanvasdeltaEEBG","",1000,800);
        h_deltaEResGeant->SetTitle("Energy loss Residual;Res;n");
        h_deltaEResGeant->Fit("gaus");
        h_deltaEResGeant->Draw();
        mccanvasdeltaEEBG->Print("E_loss/deltaEResG.png");





        gStyle->SetOptStat(0);

        TCanvas *mccanvasp = new TCanvas("mccanvasp","",1000,800);
        h_pcontact->SetTitle("Momentum at impact point;p_{cont} [GeV/c];n");
        h_pcontact->Draw();
        mccanvasp->Print("E_loss/pcontact.png");
       

        /*
        TCanvas *mccanvaspz = new TCanvas("mccanvaspz","",1000,800);
        h_pzcontact->SetTitle("Momentum at impact point z;pz_{cont} [GeV/c];n");
        h_pzcontact->Draw();
        mccanvaspz->Print("E_loss/pzcontact.png");
        */

        TCanvas *mccanvascos = new TCanvas("mccanvascos","",1000,800);
        h_costheta->SetTitle("cos#theta at impact point with plane;cos#theta;n");
        h_costheta->Draw();
        mccanvascos->Print("E_loss/costheta.png");

        ///////////////////////////////////////////////2D distributions 
       
    
        TCanvas *mccanvaspVScos = new TCanvas("mccanvaspVScos","",1000,800);
        h_pcontactVScostheta->SetTitle("Momentum as function of cos#theta;cos#theta;p [GeV/c]");
        h_pcontactVScostheta->Draw("COLZ");
        mccanvaspVScos->Print("E_loss/pVScos.png");

        TCanvas *mccanvasdeltapVScos = new TCanvas("mccanvasdeltapVScos","",1000,800);
        h_deltapVScostheta->SetTitle("Momentum loss as function of cos#theta;cos#theta;#Delta p [GeV/c]");
        h_deltapVScostheta->Draw("COLZ");
        mccanvasdeltapVScos->Print("E_loss/deltapVScos.png");

        TCanvas *mccanvasdeltapVSp = new TCanvas("mccanvasdeltapVSp","",1000,800);
        h_deltapVSpcontact->SetTitle("Momentum loss as function of p;p_{cont} [GeV/c];#Delta p [GeV/c]");
        h_deltapVSpcontact->Draw("COLZ");
        mccanvasdeltapVSp->Print("E_loss/deltapVSpcontact.png");

        TCanvas *mccanvasdeltaEVScos = new TCanvas("mccanvasdeltaEVScos","",1000,800);
        h_deltaEVScostheta->SetTitle("Energy loss as function times cos#theta of cos#theta;cos#theta;#Delta E #times cos#theta[GeV/c^{2}]");
        h_deltaEVScostheta->Draw("COLZ");
        mccanvasdeltaEVScos->Print("E_loss/deltaEVScos.png");

         TCanvas *mccanvasdeltaEexpVScos = new TCanvas("mccanvasdeltaEexpVScos","",1000,800);
        h_deltaEexpVScostheta->SetTitle("Expected energy loss times cos#theta as function of cos#theta;cos#theta;#Delta E #times cos#theta[GeV/c^{2}]");
        h_deltaEexpVScostheta->Draw("COLZ");
        mccanvasdeltaEexpVScos->Print("E_loss/deltaEexpVScos.png");



        ///////////////////////////////////////////////2D distributions XYZ

        TCanvas *mccanvasXY = new TCanvas("mccanvasXY","",1000,800);
        h_xybef->SetTitle("XY before and after plane;x[cm];y[cm]");
        h_xybef->SetMarkerStyle(3);
        h_xybef->SetMarkerColor(kRed);
        h_xybef->Draw();
        h_xypost->SetMarkerStyle(3);
        h_xypost->SetMarkerColor(kBlue);
        h_xypost->Draw("SAME");
        mccanvasXY->Print("E_loss/XY.png");


        TCanvas *mccanvasZY = new TCanvas("mccanvasZY","",1000,800);
        h_zybef->SetTitle("ZY before and after plane;z[cm];y[cm]");
        h_zybef->SetMarkerStyle(3);
        h_zybef->SetMarkerColor(kRed);
        h_zybef->Draw();
        h_zypost->SetMarkerStyle(3);
        h_zypost->SetMarkerColor(kBlue);
        h_zypost->Draw("SAME");
        mccanvasZY->Print("E_loss/ZY.png");


        //////////////////////////////////////////////////Probability plots

        TCanvas *mccanvasdeltapVScosP = new TCanvas("mccanvasdeltapVScosP","",1000,800);
        P_deltapVScostheta->SetTitle("P(#Deltap|cos#theta);cos#theta;#Delta p [GeV/c]");
        P_deltapVScostheta->Draw("COLZ");
        mccanvasdeltapVScosP->Print("E_loss/PdeltapVScos.png");

        TCanvas *mccanvasdeltapVSpP = new TCanvas("mccanvasdeltapVSpP","",1000,800);
        P_deltapVSpcontact->SetTitle("P(#Deltap|p);p_{cont} [GeV/c];#Delta p [GeV/c]");
        P_deltapVSpcontact->Draw("COLZ");
        mccanvasdeltapVSpP->Print("E_loss/PdeltapVSpcontact.png");

        TCanvas *PmccanvasdeltaEVScos = new TCanvas("PmccanvasdeltaEVScos","",1000,800);
        P_deltaEVScostheta->SetTitle("P(#DeltaE*cos#theta|cos#theta);cos#theta;#Delta E #times cos#theta[GeV/c^{2}]");
        P_deltaEVScostheta->Draw("COLZ");
        PmccanvasdeltaEVScos->Print("E_loss/PdeltaEVScos.png");

         TCanvas *PmccanvasdeltaEexpVScos = new TCanvas("PmccanvasdeltaEexpVScos","",1000,800);
        P_deltaEexpVScostheta->SetTitle("P(#DeltaE_{exp}*cos#theta|cos#theta);cos#theta;#Delta E_{exp} #times cos#theta[GeV/c^{2}]");
        P_deltaEexpVScostheta->Draw("COLZ");
        PmccanvasdeltaEVScos->Print("E_loss/PdeltaEexpVScos.png");


        ///////////////////////////////////////////////TGraphs 

        TCanvas *mccanvasdeltapVSpg = new TCanvas("mccanvasdeltapVSpg","",1000,800);
        g_deltapVSpcontact->SetTitle("Momentum loss as function of p ;p_{cont} [GeV/c]; #Delta p [GeV/c]");
        g_deltapVSpcontact->Draw("AP");
        mccanvasdeltapVSpg->Print("E_loss/G_deltapVSpcontact.png");

         TCanvas *mccanvasdeltapVScosg = new TCanvas("mccanvasdeltapVScosg","",1000,800);
        g_deltapVScostheta->SetTitle("Momentum loss as functioncos#theta;cos#theta;#Delta p [GeV/c]");
        g_deltapVScostheta->Draw("AP");
        mccanvasdeltapVScosg->Print("E_loss/G_deltapVScos.png");

        TCanvas *mccanvasdeltapVSpVScos = new TCanvas("mccanvasdeltapVSpVScos","",1000,800);
        g_deltapVScosthetaVSpcontact->SetTitle("Momentum loss as function of p and cos#theta;p_{cont} [GeV/c];cos#theta;#Delta p [GeV/c]");
        g_deltapVScosthetaVSpcontact->Draw("COLZ");
        mccanvasdeltapVSpVScos->Print("E_loss/deltapVSpcontactVScos.png");

    #pragma endregion

}