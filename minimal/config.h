#include <iostream>


using namespace std;




enum selectionType{UNCUT,DILEP,PHOTON,SEL};
enum Histograms1D{PT1,PT2,ETA1,ETA2,PHI1,PHI2,MLL,NJETS,NPHOTONS,ETMISS,HT,GENHT,NVTX,ETAG1,PHIG1,PTG1,SIGMAIETAIETAG1};
enum Histograms2D{ISRVFSR};
enum particleType{E,M};

map<Histograms1D,string> histoNames;
map<Histograms2D,string> histoNames2D;
void setHistoNames(){
histoNames[PT1]= "Pt1";
histoNames[PT2]= "Pt2";
histoNames[ETA1]= "Eta1";
histoNames[ETA2]= "Eta2";
histoNames[PHI1]= "Phi1";
histoNames[PHI2]= "Phi2";
histoNames[MLL]= "m_ll";
histoNames[NJETS]= "n_jets";
histoNames[NPHOTONS]= "n_photons";
histoNames[ETMISS]= "met";
histoNames[HT]= "ht";
histoNames[GENHT]= "gen_ht";
histoNames[NVTX]= "n_vtx";
histoNames[ETAG1]= "eta_g1";
histoNames[PHIG1]= "phi_g1";
histoNames[PTG1]= "pt_g1";
histoNames[SIGMAIETAIETAG1]= "sigmaIetaIeta_g1";
histoNames2D[ISRVFSR] = "ISRvFSR";
};


