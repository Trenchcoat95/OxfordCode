#include "/nashome/b/bellanto/util/C--/globbing.h"
#include "/nashome/b/bellanto/util/C--/concat.h"
// #include "/nashome/b/bellanto/util/C--/unorderedUnique.h"

#include "/nashome/b/bellanto/util/rootcrap/includeROOT.h"
#include "/nashome/b/bellanto/util/rootcrap/treeReader.h"
#include "/nashome/b/bellanto/util/rootcrap/PlotMap.h"
#include "../AnaUtils/analysisFuncs.h"
#include "../AnaUtils/PDGutils.h"

#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/IDNumberGen.h"

void MakeHists(TTree* Oak, PlotMap plots);



//==============================================================================
//==============================================================================
string OutName = "spatialResiduals";

enum CounterTag{nMCpart,	nPartInFid,		nMatchTrack,	
				nTPCCluster,nMatchTraj,
nCounterTag}; // Must be last enum
int Counter[nCounterTag];



enum ROCregion{CROC,IROC,IOROC,OOROC};
// Hack-in the FCL resolution parameters from trackfit2 (default, actually)
double const fTPCClusterResolYZ = 1.0;
double const fTPCClusterResolX  = 0.5;





//==============================================================================
//==============================================================================
int main(int argc , const char* argv[]){
	for (int i=0; i<nCounterTag; ++i) Counter[i] = 0;



	int maxNfiles = -1;	
	bool haveInputs = false;
	string globListFile;

	int iArg = 1;
	while ((iArg < argc) && (argv[iArg][0] == '-')){
		switch (argv[iArg][1]) {
		case 'i':		case 'I':
			haveInputs = true;
			++iArg;
			globListFile = argv[iArg];
			break;
		case 'f':		case 'F':
			++iArg;
			maxNfiles = atoi(argv[iArg]);
			if (maxNfiles <= 1) {
				cout << "maxNfiles " << maxNfiles << "; that's not right!" << endl;
				exit(1);
			}
			break;
		// Special top secret case for top secret things.  Like debugging
		case 'L':
			cout << "No L option today!" << endl; exit(1);
			break;
		case 'h':
		default:
			cout << "-iI filename	File listing input globbing patterns.  1stline in" << endl;
			cout << "               file is globbing pattern for input files; 2nd line" << endl;
			cout << "               is ignored.  PNFS space is from xrootd, e.g." << endl;
			cout << "               \"dCache */OverlaidDST_*.root\" then will get all the";
			cout << "               files matching */OverlaidDST_*.root from" << endl;
			cout << "               /pnfs/dune/persistent/users/bellanto/" << endl;
			cout << "[-fF n]		Analyze n background files, max" << endl;
			cout << "[-h]			This message" << endl;
			exit(0);
		}
		++iArg;
    }
	if (!haveInputs) {
		cout << "You need an -i or -I argument in the command line" << endl;
		exit(1);
	}



	// Get list of input files
	string LookFor;
	std::ifstream globListInStream(globListFile);
	int			  	   nFiles;
	std::vector<string> FileList;

	globListInStream >> LookFor;	
	if (LookFor == "dCache") {
		globListInStream >> LookFor;
		LookFor = "/pnfs/dune/persistent/users/bellanto/" + LookFor;
		FileList = globbing(LookFor);
		string newAccess = "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr";
		std::vector<string>::iterator it = FileList.begin();
		for (; it<FileList.end(); ++it) {
			it->erase(0,5);
			*it = newAccess + *it;
		}
	} else {
		FileList = globbing(LookFor);
	}
	globListInStream.close();
	cout << ( nFiles = (int)FileList.size() ) << " file(s) found.\n";
	if (nFiles==0) exit(2);

	if ( maxNfiles>0 && maxNfiles<nFiles) {
		nFiles = maxNfiles;
		cout << "Processing only " << nFiles << " files." << endl;
	}
	cout << endl;





	// Output files and the plots therein
	PlotMap plots;
	Int_t nBinsX,nBinsY;		Double_t loX,loY, hiX, hiY;

	// Matching MC to reco tracks
	nBinsX = 100;		loX = 0.95;			hiX =  +1.0;
	nBinsY =  90;		loY = 0.00;			hiY = +20.0;
	plots.addTH2D("RECOmatch", "#delta(x) vs. cos(#theta) for signal MC-track match",
				nBinsX,loX,hiX, nBinsY,loY,hiY);

	// A reco'd momentum plot
	nBinsX = 25;		loX = 0.0;			hiX = +5.0;
	plots.addTH1D("RECO_P",		"Momentum of matched tracks",
				nBinsX,loX,hiX);

	// What's the drift distance to the closest traj point?
	nBinsX = 100;		loX = 0.0;			hiX = +10.0;
	plots.addTH1D("TrajToClus","Drift distance TPCCluster to closest MC traj. point",
				nBinsX,loX,hiX);

	// What's the spacing between consecutive trajectory points?
	nBinsX = 100;		loX = 0.0;			hiX =  +1.0;
	plots.addTH1D("TrajDeltaX",	"Drift distance spacing consecutive MC traj. points",
				nBinsX,loX,hiX);
	nBinsX = 100;		loX = 0.0;			hiX = +10.0;
	plots.addTH1D("TrajDeltaT",	"Perp. place distance spacing consecutive MC traj. points",
				nBinsX,loX,hiX);

	// Plots for track residual studies
	nBinsX = 45;		loX = 0.0;			hiX = M_PI/4.0;
	plots.addTH1D("impactAngle","MC track angle re long axis of pad",
				nBinsX,loX,hiX);
	
	nBinsX = nBinsY = 100;		loX = loY = -5.00;		hiX = hiY = +5.00;
	plots.addTH2D("residuaL_2D","(z,y) TPCCluster residuals large scale",
				nBinsX,loX,hiX, nBinsY,loY,hiY);

	nBinsX = nBinsY = 80;		loX = loY = -0.40;		hiX = hiY = +0.40;
	plots.addTH2D("resid_CROC_2D","(z,y) TPCCluster residuals,  CROC (cm)",
				nBinsX,loX,hiX, nBinsY,loY,hiY);
	plots.addTH2D("resid_IROC_2D","(z,y) TPCCluster residuals,  IROC (cm)",
				nBinsX,loX,hiX, nBinsY,loY,hiY);
	plots.addTH2D("residIOROC_2D","(z,y) TPCCluster residuals, IOROC (cm)",
				nBinsX,loX,hiX, nBinsY,loY,hiY);
	plots.addTH2D("residOOROC_2D","(z,y) TPCCluster residuals, OOROC (cm)",
				nBinsX,loX,hiX, nBinsY,loY,hiY);

	nBinsX = 60;				loX = 0.00;				hiX = +0.60;
	plots.addTH1D("resid_CROC_1D_0","TPCCluster residuals,  CROC (cm) 0.0 < #alpha < 0.1",
				nBinsX,loX,hiX);
	plots.addTH1D("resid_IROC_1D_0","TPCCluster residuals,  IROC (cm) 0.0 < #alpha < 0.1",
				nBinsX,loX,hiX);
	plots.addTH1D("residIOROC_1D_0","TPCCluster residuals, IOROC (cm) 0.0 < #alpha < 0.1",
				nBinsX,loX,hiX);
	plots.addTH1D("residOOROC_1D_0","TPCCluster residuals, OOROC (cm) 0.0 < #alpha < 0.1",
				nBinsX,loX,hiX);

	plots.addTH1D("resid_CROC_1D_1","TPCCluster residuals,  CROC (cm) 0.1 < #alpha < 0.2",
				nBinsX,loX,hiX);
	plots.addTH1D("resid_IROC_1D_1","TPCCluster residuals,  IROC (cm) 0.1 < #alpha < 0.2",
				nBinsX,loX,hiX);
	plots.addTH1D("residIOROC_1D_1","TPCCluster residuals, IOROC (cm) 0.1 < #alpha < 0.2",
				nBinsX,loX,hiX);
	plots.addTH1D("residOOROC_1D_1","TPCCluster residuals, OOROC (cm) 0.1 < #alpha < 0.2",
				nBinsX,loX,hiX);

	plots.addTH1D("resid_CROC_1D_2","TPCCluster residuals,  CROC (cm) 0.2 < #alpha < 0.3",
				nBinsX,loX,hiX);
	plots.addTH1D("resid_IROC_1D_2","TPCCluster residuals,  IROC (cm) 0.2 < #alpha < 0.3",
				nBinsX,loX,hiX);
	plots.addTH1D("residIOROC_1D_2","TPCCluster residuals, IOROC (cm) 0.2 < #alpha < 0.3",
				nBinsX,loX,hiX);
	plots.addTH1D("residOOROC_1D_2","TPCCluster residuals, OOROC (cm) 0.2 < #alpha < 0.3",
				nBinsX,loX,hiX);

	plots.addTH1D("resid_CROC_1D_3","TPCCluster residuals,  CROC (cm) 0.3 < #alpha < 0.4",
				nBinsX,loX,hiX);
	plots.addTH1D("resid_IROC_1D_3","TPCCluster residuals,  IROC (cm) 0.3 < #alpha < 0.4",
				nBinsX,loX,hiX);
	plots.addTH1D("residIOROC_1D_3","TPCCluster residuals, IOROC (cm) 0.3 < #alpha < 0.4",
				nBinsX,loX,hiX);
	plots.addTH1D("residOOROC_1D_3","TPCCluster residuals, OOROC (cm) 0.3 < #alpha < 0.4",
				nBinsX,loX,hiX);

	plots.addTH1D("resid_CROC_1D_4","TPCCluster residuals,  CROC (cm) 0.4 < #alpha < 0.5",
				nBinsX,loX,hiX);
	plots.addTH1D("resid_IROC_1D_4","TPCCluster residuals,  IROC (cm) 0.4 < #alpha < 0.5",
				nBinsX,loX,hiX);
	plots.addTH1D("residIOROC_1D_4","TPCCluster residuals, IOROC (cm) 0.4 < #alpha < 0.5",
				nBinsX,loX,hiX);
	plots.addTH1D("residOOROC_1D_4","TPCCluster residuals, OOROC (cm) 0.4 < #alpha < 0.5",
				nBinsX,loX,hiX);

	plots.addTH1D("resid_CROC_1D_5","TPCCluster residuals,  CROC (cm) 0.5 < #alpha < 0.6",
				nBinsX,loX,hiX);
	plots.addTH1D("resid_IROC_1D_5","TPCCluster residuals,  IROC (cm) 0.5 < #alpha < 0.6",
				nBinsX,loX,hiX);
	plots.addTH1D("residIOROC_1D_5","TPCCluster residuals, IOROC (cm) 0.5 < #alpha < 0.6",
				nBinsX,loX,hiX);
	plots.addTH1D("residOOROC_1D_5","TPCCluster residuals, OOROC (cm) 0.5 < #alpha < 0.6",
				nBinsX,loX,hiX);

	plots.addTH1D("resid_CROC_1D_6","TPCCluster residuals,  CROC (cm) 0.6 < #alpha < 0.7",
				nBinsX,loX,hiX);
	plots.addTH1D("resid_IROC_1D_6","TPCCluster residuals,  IROC (cm) 0.6 < #alpha < 0.7",
				nBinsX,loX,hiX);
	plots.addTH1D("residIOROC_1D_6","TPCCluster residuals, IOROC (cm) 0.6 < #alpha < 0.7",
				nBinsX,loX,hiX);
	plots.addTH1D("residOOROC_1D_6","TPCCluster residuals, OOROC (cm) 0.6 < #alpha < 0.7",
				nBinsX,loX,hiX);

	plots.addTH1D("resid_CROC_1D_7","TPCCluster residuals,  CROC (cm) 0.7 < #alpha < 0.8",
				nBinsX,loX,hiX);
	plots.addTH1D("resid_IROC_1D_7","TPCCluster residuals,  IROC (cm) 0.7 < #alpha < 0.8",
				nBinsX,loX,hiX);
	plots.addTH1D("residIOROC_1D_7","TPCCluster residuals, IOROC (cm) 0.7 < #alpha < 0.8",
				nBinsX,loX,hiX);
	plots.addTH1D("residOOROC_1D_7","TPCCluster residuals, OOROC (cm) 0.7 < #alpha < 0.8",
				nBinsX,loX,hiX);

	plots.addTH1D("resid_CROC_1D_8","TPCCluster residuals,  CROC (cm) 0.8 < #alpha < 0.9",
				nBinsX,loX,hiX);
	plots.addTH1D("resid_IROC_1D_8","TPCCluster residuals,  IROC (cm) 0.8 < #alpha < 0.9",
				nBinsX,loX,hiX);
	plots.addTH1D("residIOROC_1D_8","TPCCluster residuals, IOROC (cm) 0.8 < #alpha < 0.9",
				nBinsX,loX,hiX);
	plots.addTH1D("residOOROC_1D_8","TPCCluster residuals, OOROC (cm) 0.8 < #alpha < 0.9",
				nBinsX,loX,hiX);

	plots.addTH1D("resid_CROC_1D_9","TPCCluster residuals,  CROC (cm) 0.9 < #alpha < 1.0",
				nBinsX,loX,hiX);
	plots.addTH1D("resid_IROC_1D_9","TPCCluster residuals,  IROC (cm) 0.9 < #alpha < 1.0",
				nBinsX,loX,hiX);
	plots.addTH1D("residIOROC_1D_9","TPCCluster residuals, IOROC (cm) 0.9 < #alpha < 1.0",
				nBinsX,loX,hiX);
	plots.addTH1D("residOOROC_1D_9","TPCCluster residuals, OOROC (cm) 0.9 < #alpha < 1.0",
				nBinsX,loX,hiX);

	plots.addTH1D("resid_CROC_1D_A","TPCCluster residuals,  CROC (cm) 1.0 < #alpha < 1.1",
				nBinsX,loX,hiX);
	plots.addTH1D("resid_IROC_1D_A","TPCCluster residuals,  IROC (cm) 1.0 < #alpha < 1.1",
				nBinsX,loX,hiX);
	plots.addTH1D("residIOROC_1D_A","TPCCluster residuals, IOROC (cm) 1.0 < #alpha < 1.1",
				nBinsX,loX,hiX);
	plots.addTH1D("residOOROC_1D_A","TPCCluster residuals, OOROC (cm) 1.0 < #alpha < 1.1",
				nBinsX,loX,hiX);

	plots.addTH1D("resid_CROC_1D_B","TPCCluster residuals,  CROC (cm) 1.1 < #alpha < 1.2",
				nBinsX,loX,hiX);
	plots.addTH1D("resid_IROC_1D_B","TPCCluster residuals,  IROC (cm) 1.1 < #alpha < 1.2",
				nBinsX,loX,hiX);
	plots.addTH1D("residIOROC_1D_B","TPCCluster residuals, IOROC (cm) 1.1 < #alpha < 1.2",
				nBinsX,loX,hiX);
	plots.addTH1D("residOOROC_1D_B","TPCCluster residuals, OOROC (cm) 1.1 < #alpha < 1.2",
				nBinsX,loX,hiX);

	plots.addTH1D("resid_CROC_1D_C","TPCCluster residuals,  CROC (cm) 1.2 < #alpha < 1.3",
				nBinsX,loX,hiX);
	plots.addTH1D("resid_IROC_1D_C","TPCCluster residuals,  IROC (cm) 1.2 < #alpha < 1.3",
				nBinsX,loX,hiX);
	plots.addTH1D("residIOROC_1D_C","TPCCluster residuals, IOROC (cm) 1.2 < #alpha < 1.3",
				nBinsX,loX,hiX);
	plots.addTH1D("residOOROC_1D_C","TPCCluster residuals, OOROC (cm) 1.2 < #alpha < 1.3",
				nBinsX,loX,hiX);

	plots.addTH1D("resid_CROC_1D_D","TPCCluster residuals,  CROC (cm) 1.3 < #alpha < 1.4",
				nBinsX,loX,hiX);
	plots.addTH1D("resid_IROC_1D_D","TPCCluster residuals,  IROC (cm) 1.3 < #alpha < 1.4",
				nBinsX,loX,hiX);
	plots.addTH1D("residIOROC_1D_D","TPCCluster residuals, IOROC (cm) 1.3 < #alpha < 1.4",
				nBinsX,loX,hiX);
	plots.addTH1D("residOOROC_1D_D","TPCCluster residuals, OOROC (cm) 1.3 < #alpha < 1.4",
				nBinsX,loX,hiX);

	plots.addTH1D("resid_CROC_1D_E","TPCCluster residuals,  CROC (cm) 1.4 < #alpha < 1.5",
				nBinsX,loX,hiX);
	plots.addTH1D("resid_IROC_1D_E","TPCCluster residuals,  IROC (cm) 1.4 < #alpha < 1.5",
				nBinsX,loX,hiX);
	plots.addTH1D("residIOROC_1D_E","TPCCluster residuals, IOROC (cm) 1.4 < #alpha < 1.5",
				nBinsX,loX,hiX);
	plots.addTH1D("residOOROC_1D_E","TPCCluster residuals, OOROC (cm) 1.4 < #alpha < 1.5",
				nBinsX,loX,hiX);

	plots.addTH1D("resid_CROC_1D_F","TPCCluster residuals,  CROC (cm) 1.5 < #alpha < #pi",
				nBinsX,loX,hiX);
	plots.addTH1D("resid_IROC_1D_F","TPCCluster residuals,  IROC (cm) 1.5 < #alpha < #pi",
				nBinsX,loX,hiX);
	plots.addTH1D("residIOROC_1D_F","TPCCluster residuals, IOROC (cm) 1.5 < #alpha < #pi",
				nBinsX,loX,hiX);
	plots.addTH1D("residOOROC_1D_F","TPCCluster residuals, OOROC (cm) 1.5 < #alpha < #pi",
				nBinsX,loX,hiX);




	
	time_t unixT;	struct tm* localT;
	time(&unixT);	localT = localtime(&unixT);
	cout << "Starting at " << asctime(localT) << endl;
	
	for (int iFile=0; iFile<nFiles; ++iFile) {
		cout << "Reading file " << iFile+1 << ": " << FileList[iFile].c_str() << " ";
		TFile* f = TFile::Open(FileList[iFile].c_str(),"READ");
		if (!f || !f->IsOpen()) {
			cout << "...but that file can not be opened!" << endl;
			exit(3);
		}
		TTree* treeInFile;		
		TDirectory* dir = (TDirectory*)(f)->Get("anatree");
		dir->GetObject("GArAnaTree",treeInFile);
		MakeHists(treeInFile, plots);
		f->Close();
	}
	time(&unixT);	localT = localtime(&unixT);
	cout << "Stopping at " << asctime(localT) << endl;





	string TextOutName = OutName;
	TextOutName += ".txt";
	// could change to TextFileOut if you want!
	#define OUTDEV cout
	//std::ofstream TextFileOut(TextOutName.c_str());

	OUTDEV << Counter[nMCpart]      << "\tMC primaries analyzed, with\n";
	OUTDEV << Counter[nPartInFid]   << "\tof them in the fiducial, and\n";
	OUTDEV << Counter[nMatchTrack]	<< "\tare well matched to tracks.\n\n";
	OUTDEV << Counter[nTPCCluster]	<< "\tTPC Clusters, of which\n";
	OUTDEV << Counter[nMatchTraj]   << "\tmatch well to MC trajectories.\n";

	//TextFileOut.close();





	plots.writePlots(OutName);

    exit(0);
}



//==============================================================================
//==============================================================================
void MakeHists(TTree* Oak, PlotMap plots) {


	Long64_t iEntry;

	vectorFromTree<Int_t>		MCP_PDGMother		 (Oak,"PDGMother",&iEntry);
	vectorFromTree<Int_t>		MCP_PDG				 (Oak,"PDG",       &iEntry);
	vectorFromTree<Float_t>		MCP_X				 (Oak,"MCPStartX", &iEntry);
	vectorFromTree<Float_t>		MCP_Y				 (Oak,"MCPStartY", &iEntry);
	vectorFromTree<Float_t>		MCP_Z				 (Oak,"MCPStartZ", &iEntry);
	vectorFromTree<Float_t>		MCP_PX				 (Oak,"MCPStartPX",&iEntry);
	vectorFromTree<Float_t>		MCP_PY				 (Oak,"MCPStartPY",&iEntry);
	vectorFromTree<Float_t>		MCP_PZ				 (Oak,"MCPStartPZ",&iEntry);

	vectorFromTree<Float_t>		TrajMCPX			 (Oak,"TrajMCPX",&iEntry);
	vectorFromTree<Float_t>		TrajMCPY			 (Oak,"TrajMCPY",&iEntry);
	vectorFromTree<Float_t>		TrajMCPZ			 (Oak,"TrajMCPZ",&iEntry);
	vectorFromTree<Int_t>		TrajMCPIndex		 (Oak,"TrajMCPIndex",&iEntry);

	vectorFromTree<ULong64_t>	TrackIDNumber		 (Oak,"TrackIDNumber",&iEntry);
	vectorFromTree<Float_t>		TrackStartX			 (Oak,"TrackStartX",  &iEntry);
	vectorFromTree<Float_t>		TrackStartY			 (Oak,"TrackStartY",  &iEntry);
	vectorFromTree<Float_t>		TrackStartZ			 (Oak,"TrackStartZ",  &iEntry);
	vectorFromTree<Float_t>		TrackStartPX		 (Oak,"TrackStartPX", &iEntry);
	vectorFromTree<Float_t>		TrackStartPY		 (Oak,"TrackStartPY", &iEntry);
	vectorFromTree<Float_t>		TrackStartPZ		 (Oak,"TrackStartPZ", &iEntry);
	vectorFromTree<Float_t>		TrackEndX			 (Oak,"TrackEndX",  &iEntry);
	vectorFromTree<Float_t>		TrackEndY			 (Oak,"TrackEndY",  &iEntry);
	vectorFromTree<Float_t>		TrackEndZ			 (Oak,"TrackEndZ",  &iEntry);
	vectorFromTree<Float_t>		TrackEndPX			 (Oak,"TrackEndPX", &iEntry);
	vectorFromTree<Float_t>		TrackEndPY			 (Oak,"TrackEndPY", &iEntry);
	vectorFromTree<Float_t>		TrackEndPZ			 (Oak,"TrackEndPZ", &iEntry);

	vectorFromTree<ULong64_t>	NTPCClustersOnTrack  (Oak,"NTPCClustersOnTrack",&iEntry);
	vectorFromTree<Float_t>		TPCClusterX          (Oak,"TPCClusterX",&iEntry);
	vectorFromTree<Float_t>		TPCClusterY          (Oak,"TPCClusterY",&iEntry);
	vectorFromTree<Float_t>		TPCClusterZ          (Oak,"TPCClusterZ",&iEntry);
	vectorFromTree<ULong64_t>	TPCClusterTrkIDNumber(Oak,"TPCClusterTrkIDNumber",&iEntry);
	vectorFromTree<Float_t>		TPCClusterCovYY      (Oak,"TPCClusterCovYY",&iEntry);
	vectorFromTree<Float_t>		TPCClusterCovYZ      (Oak,"TPCClusterCovYZ",&iEntry);
	vectorFromTree<Float_t>		TPCClusterCovZZ      (Oak,"TPCClusterCovZZ",&iEntry);





	Long64_t nEntries = Oak->GetEntries();
	cout << "nEntries:\t" << nEntries << endl;
	int const hullo = 20000;
	for (iEntry=0; iEntry<nEntries; ++iEntry) {
		if (nEntries>hullo && iEntry%hullo==0 && iEntry!=0)
			cout << "\t\t...iEntry = " << iEntry << endl;


		vector<int> PDGsToUse = {muonPDG,pichPDG};
		for (int& iPDG : PDGsToUse) {

			int nMCPs = MCP_X.size();
			int iMCpart = -1;		float MC_X,MC_Y,MC_Z;
			for (int iMCP=0; iMCP<nMCPs; ++iMCP) {

				// Primaries only;
				if ( MCP_PDGMother(iMCP) != 0 ) continue;

				// Start with muons; then pions
				int PDG = MCP_PDG(iMCP);
				if ( abs(PDG) != iPDG  ) continue;
				++Counter[nMCpart];

				// inFid means track starts in fiducial
				MC_X = MCP_X(iMCP);
				MC_Y = MCP_Y(iMCP);
				MC_Z = MCP_Z(iMCP);
				if ( !inFiducial(MC_X, MC_Y, MC_Z) ) continue;
				++Counter[nPartInFid];

				iMCpart = iMCP;
				break;
			}
			if (iMCpart==-1) continue;
			TVector3 MCpart( MCP_PX(iMCpart),MCP_PY(iMCpart),MCP_PZ(iMCpart) );



			// Find a track that matches the MC muon.  Still done the hard way.
			int nTracks = TrackIDNumber.size();			int iRECOpart = -1;
			TVector3 MCpartHat = MCpart.Unit();			double cosTmax = -BIG;
			gar::rec::TrackEnd whichEnd;				TVector3 RECO_P;
			for (int iTrack=0; iTrack<nTracks; ++iTrack) {
				TVector3 RECO_P_beg(TrackStartPX(iTrack),
									TrackStartPY(iTrack),TrackStartPZ(iTrack));
				TVector3 RECO_P_end(TrackEndPX(iTrack),
									TrackEndPY(iTrack),  TrackEndPZ(iTrack));

				// Direction matching
				double cosTbeg = MCpartHat.Dot(RECO_P_beg)/RECO_P_beg.Mag();
				double cosTend = MCpartHat.Dot(RECO_P_end)/RECO_P_end.Mag();
				if (cosTbeg > cosTend) {
					if (cosTbeg > cosTmax) {
						cosTmax   = cosTbeg;
						iRECOpart = iTrack;
						whichEnd  = gar::rec::TrackEndBeg;
						RECO_P	  = RECO_P_beg;
					}
				} else {
					if (cosTend > cosTmax) {
						cosTmax   = cosTend;
						iRECOpart = iTrack;
						whichEnd  = gar::rec::TrackEndEnd;
						RECO_P	  = RECO_P_beg;
					}
				}
			} // End loop over all Tracks
			if (iRECOpart == -1) continue;

			// Plot offset of reco track
			TVector3 vecX;
			if ( whichEnd == gar::rec::TrackEndBeg ) {
				vecX.SetXYZ(TrackStartX(iRECOpart) -MC_X,
							TrackStartY(iRECOpart) -MC_Y,
							TrackStartZ(iRECOpart) -MC_Z);
			} else {
				vecX.SetXYZ(TrackEndX(iRECOpart)   -MC_X,
							TrackEndY(iRECOpart)   -MC_Y,
							TrackEndZ(iRECOpart)   -MC_Z);
			}
			Float_t delX = vecX.Cross(MCpartHat).Mag();
			plots["RECOmatch"]->Fill(cosTmax,delX);

			// And cut to match
			if ( cosTmax <= 0.997 ) continue;
			if (  delX   >= 3.0   ) continue;
			++Counter[nMatchTrack];
			ULong64_t iRECO_IDNumber = TrackIDNumber(iRECOpart);

			// Track momentum might be interesting here.
			plots["RECO_P"]->Fill( RECO_P.Mag() );

			// Get the relevant TPCclusters
			TVector3 aTPCCluster;		vector<TVector3> ourTPCClusters;
			int nTPCClusters = TPCClusterTrkIDNumber.size();
			for (int iTPCCluster=0; iTPCCluster<nTPCClusters; ++iTPCCluster) {
				if ( TPCClusterTrkIDNumber(iTPCCluster) == iRECO_IDNumber ) {
					aTPCCluster.SetXYZ( TPCClusterX(iTPCCluster),
										TPCClusterY(iTPCCluster),
										TPCClusterZ(iTPCCluster) );
					ourTPCClusters.push_back(aTPCCluster);
				}
			}

			if ( ourTPCClusters.size() != NTPCClustersOnTrack(iRECOpart) ) {
				cout << "Huh." << ourTPCClusters.size() << " vs. " <<
					NTPCClustersOnTrack(iRECOpart) << " at iEntry = " << iEntry
					<< endl;
				continue;
			}



			// Get the corresponding trajectory points
			TVector3 aTrajPoint;		vector<TVector3> ourTrajPoints;
			int nTrajPoints = TrajMCPIndex.size();
			for (int iTrajPoint=0; iTrajPoint<nTrajPoints; ++iTrajPoint) {
				if ( TrajMCPIndex(iTrajPoint) == iMCpart ) {
					aTrajPoint.SetXYZ( TrajMCPX(iTrajPoint),
									   TrajMCPY(iTrajPoint),
									   TrajMCPZ(iTrajPoint) );
					ourTrajPoints.push_back(aTrajPoint);
				}
			}



			// Match TPCclusters to points
			int nOurTPCClusters = ourTPCClusters.size();
			int nOurTrajPoints  = ourTrajPoints.size();
			if ( nOurTrajPoints <2 ) {
				cout << "Really? " << ourTrajPoints.size() << " trajectory points? "
					<< " at iEntry = " << iEntry<< endl;
				continue;
			}



			for (int iTPCCluster=0; iTPCCluster<nOurTPCClusters; ++iTPCCluster) {
				// Find trajectory point closest to this cluster in 3d and then which
				// of it's 2 neighbors is next closest.  Don't get the 2nd closest 
				// trajectory point as (in principle) that could be on another loop
				// of the helix.
				vector<double> distances(nOurTrajPoints);
				for (int iDist=0; iDist<nOurTrajPoints; ++iDist) {
					distances[iDist] = ( ourTrajPoints[iDist] - ourTPCClusters[iTPCCluster] ).Mag();
				}
				int iClosest  = std::min_element(distances.begin(),distances.end())
							   -distances.begin();
				plots["TrajToClus"]->Fill( distances[iClosest] );

				++Counter[nTPCCluster];
				// Following cut would bias the residuals!
				if ( distances[iClosest] >= BIG ) continue;
				++Counter[nMatchTraj];

				int iSecondest;
				if ( iClosest == 0 ) {
					iSecondest = 1;
				} else if ( iClosest == (nOurTrajPoints-1) ) {
					iSecondest = iClosest -1;
				} else {
					iSecondest = ( distances[iClosest-1] < distances[iClosest+1] ) ?
						iClosest-1 : iClosest+1;
				}

				// Plot some info about MC trajectory point spacing
				double dex = fabs( ourTrajPoints[iSecondest].X() -ourTrajPoints[iClosest].X() );
				plots["TrajDeltaX"]->Fill(dex);
				double dep = Qadd(ourTrajPoints[iSecondest].Z() -ourTrajPoints[iClosest].Z(),
								  ourTrajPoints[iSecondest].Y() -ourTrajPoints[iClosest].Y());
				plots["TrajDeltaT"]->Fill(dep);



				// Linear interp between these two trajectory points to get (y,z) of
				// trajectory at point closest to measured TPCCluster position, using the 
				// definition of closest that comes out of Reco/tpctrackfit2_module.cc.
				double smallXi = ourTrajPoints[iClosest].X();
				double Yslope = (ourTrajPoints[iSecondest].Y() -ourTrajPoints[iClosest].Y())
							   /(ourTrajPoints[iSecondest].X() -smallXi);
				double Yinter = ourTrajPoints[iClosest].Y();
				double Zslope = (ourTrajPoints[iSecondest].Z() -ourTrajPoints[iClosest].Z())
							   /(ourTrajPoints[iSecondest].X() -smallXi);
				double Zinter = ourTrajPoints[iClosest].Z();

				double sigmaX2  = fTPCClusterResolX *fTPCClusterResolX;
				double sigmaYZ2 = fTPCClusterResolYZ*fTPCClusterResolYZ;

				double xClus = ourTPCClusters[iTPCCluster].X();
				double yClus = ourTPCClusters[iTPCCluster].Y();
				double zClus = ourTPCClusters[iTPCCluster].Z();

				double XcloseNum, XcloseDen, Xclosest, Yclosest, Zclosest;
				XcloseNum  = Yslope*(Yinter -Yslope*smallXi) +Zslope*(Zinter -Zslope*smallXi)
							-(yClus*Yslope +zClus*Zslope);
				XcloseNum /= sigmaYZ2;
				XcloseNum -= xClus/sigmaX2;
				XcloseNum *= -1.0;
				XcloseDen  = 1.0/sigmaX2 + (Yslope*Yslope + Zslope*Zslope)/sigmaYZ2;
				Xclosest   = XcloseNum / XcloseDen;
				Yclosest   = Yslope*(Xclosest -smallXi) +Yinter;
				Zclosest   = Zslope*(Xclosest -smallXi) +Zinter;

				// And now you have the residuals
				float Yresid = ourTPCClusters[iTPCCluster].Y() -Yclosest;
				float Zresid = ourTPCClusters[iTPCCluster].Z() -Zclosest;
				float residual = Qadd(Zresid,Yresid);




				// Where to plot this residual?
				double impactAngle;
				TVector3 trajPerp     = ourTrajPoints[iClosest];
				trajPerp.SetX(0.0);			double magTrajPerp     = trajPerp.Mag();
				TVector3 trajStepPerp = ourTrajPoints[iSecondest] -ourTrajPoints[iClosest];
				trajStepPerp.SetX(0.0);		double magTrajStepPerp = trajStepPerp.Mag();
				impactAngle = trajPerp.Dot(trajStepPerp) / (magTrajPerp*magTrajStepPerp);
				impactAngle = acos(abs(impactAngle));
				plots["impactAngle"]->Fill(impactAngle);

				ROCregion InROC;
				if (magTrajPerp < 84.0) {
					InROC = CROC;
				} else if (magTrajPerp < 133.5) {
					InROC = IROC;
				} else if (magTrajPerp < 200.0) {
					InROC = IOROC;
				} else {
					InROC = OOROC;
				}

				plots["residuaL_2D"]->Fill(Zresid,Yresid);

				int alpha = impactAngle * 10;	// Truncate
				if (InROC == CROC) {
					plots["resid_CROC_2D"]->Fill(Zresid,Yresid);
					switch (alpha) {
						case 0:		plots["resid_CROC_1D_0"]->Fill(residual);	break;
						case 1:		plots["resid_CROC_1D_1"]->Fill(residual);	break;
						case 2:		plots["resid_CROC_1D_2"]->Fill(residual);	break;
						case 3:		plots["resid_CROC_1D_3"]->Fill(residual);	break;
						case 4:		plots["resid_CROC_1D_4"]->Fill(residual);	break;
						case 5:		plots["resid_CROC_1D_5"]->Fill(residual);	break;
						case 6:		plots["resid_CROC_1D_6"]->Fill(residual);	break;
						case 7:		plots["resid_CROC_1D_7"]->Fill(residual);	break;
						case 8:		plots["resid_CROC_1D_8"]->Fill(residual);	break;
						case 9:		plots["resid_CROC_1D_9"]->Fill(residual);	break;
						case 10:	plots["resid_CROC_1D_A"]->Fill(residual);	break;
						case 11:	plots["resid_CROC_1D_B"]->Fill(residual);	break;
						case 12:	plots["resid_CROC_1D_C"]->Fill(residual);	break;
						case 13:	plots["resid_CROC_1D_D"]->Fill(residual);	break;
						case 14:	plots["resid_CROC_1D_E"]->Fill(residual);	break;
						case 15:	plots["resid_CROC_1D_F"]->Fill(residual);	break;
					}

				} else if (InROC == IROC) {
					plots["resid_IROC_2D"]->Fill(Zresid,Yresid);
					switch (alpha) {
						case 0:		plots["resid_IROC_1D_0"]->Fill(residual);	break;
						case 1:		plots["resid_IROC_1D_1"]->Fill(residual);	break;
						case 2:		plots["resid_IROC_1D_2"]->Fill(residual);	break;
						case 3:		plots["resid_IROC_1D_3"]->Fill(residual);	break;
						case 4:		plots["resid_IROC_1D_4"]->Fill(residual);	break;
						case 5:		plots["resid_IROC_1D_5"]->Fill(residual);	break;
						case 6:		plots["resid_IROC_1D_6"]->Fill(residual);	break;
						case 7:		plots["resid_IROC_1D_7"]->Fill(residual);	break;
						case 8:		plots["resid_IROC_1D_8"]->Fill(residual);	break;
						case 9:		plots["resid_IROC_1D_9"]->Fill(residual);	break;
						case 10:	plots["resid_IROC_1D_A"]->Fill(residual);	break;
						case 11:	plots["resid_IROC_1D_B"]->Fill(residual);	break;
						case 12:	plots["resid_IROC_1D_C"]->Fill(residual);	break;
						case 13:	plots["resid_IROC_1D_D"]->Fill(residual);	break;
						case 14:	plots["resid_IROC_1D_E"]->Fill(residual);	break;
						case 15:	plots["resid_IROC_1D_F"]->Fill(residual);	break;
					}

				} else if (InROC == IOROC) {
					plots["residIOROC_2D"]->Fill(Zresid,Yresid);
					switch (alpha) {
						case 0:		plots["residIOROC_1D_0"]->Fill(residual);	break;
						case 1:		plots["residIOROC_1D_1"]->Fill(residual);	break;
						case 2:		plots["residIOROC_1D_2"]->Fill(residual);	break;
						case 3:		plots["residIOROC_1D_3"]->Fill(residual);	break;
						case 4:		plots["residIOROC_1D_4"]->Fill(residual);	break;
						case 5:		plots["residIOROC_1D_5"]->Fill(residual);	break;
						case 6:		plots["residIOROC_1D_6"]->Fill(residual);	break;
						case 7:		plots["residIOROC_1D_7"]->Fill(residual);	break;
						case 8:		plots["residIOROC_1D_8"]->Fill(residual);	break;
						case 9:		plots["residIOROC_1D_9"]->Fill(residual);	break;
						case 10:	plots["residIOROC_1D_A"]->Fill(residual);	break;
						case 11:	plots["residIOROC_1D_B"]->Fill(residual);	break;
						case 12:	plots["residIOROC_1D_C"]->Fill(residual);	break;
						case 13:	plots["residIOROC_1D_D"]->Fill(residual);	break;
						case 14:	plots["residIOROC_1D_E"]->Fill(residual);	break;
						case 15:	plots["residIOROC_1D_F"]->Fill(residual);	break;
					}

				} else {
					plots["residOOROC_2D"]->Fill(Zresid,Yresid);
					switch (alpha) {
						case 0:		plots["residOOROC_1D_0"]->Fill(residual);	break;
						case 1:		plots["residOOROC_1D_1"]->Fill(residual);	break;
						case 2:		plots["residOOROC_1D_2"]->Fill(residual);	break;
						case 3:		plots["residOOROC_1D_3"]->Fill(residual);	break;
						case 4:		plots["residOOROC_1D_4"]->Fill(residual);	break;
						case 5:		plots["residOOROC_1D_5"]->Fill(residual);	break;
						case 6:		plots["residOOROC_1D_6"]->Fill(residual);	break;
						case 7:		plots["residOOROC_1D_7"]->Fill(residual);	break;
						case 8:		plots["residOOROC_1D_8"]->Fill(residual);	break;
						case 9:		plots["residOOROC_1D_9"]->Fill(residual);	break;
						case 10:	plots["residOOROC_1D_A"]->Fill(residual);	break;
						case 11:	plots["residOOROC_1D_B"]->Fill(residual);	break;
						case 12:	plots["residOOROC_1D_C"]->Fill(residual);	break;
						case 13:	plots["residOOROC_1D_D"]->Fill(residual);	break;
						case 14:	plots["residOOROC_1D_E"]->Fill(residual);	break;
						case 15:	plots["residOOROC_1D_F"]->Fill(residual);	break;
					}
				}
			}	// end loop on TPCclusters for this track/MCpart.
		}	// end loop on mu, pi	


	}	// end loop on iEntry
	return;
}
