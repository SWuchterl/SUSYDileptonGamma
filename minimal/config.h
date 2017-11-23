#include <iostream>
using namespace std;



enum selectionType{UNCUT,DILEP,PHOTON,SEL,ONZ,TRIGDILEP,TRIGSEL,TRIGONZ,TRIGDILEP_ptcuts,TRIGONZ_ptcuts,TRIGSEL_ptcuts};
enum Histograms1D{PT1,PT2,ETA1,ETA2,PHI1,PHI2,MLL,NJETS,NPHOTONS,ETMISS,HT,GENHT,NVTX,ETAG1,PHIG1,PTG1,SIGMAIETAIETAG1,DeltaEtaLL,DeltaPhiLL,DeltaEtaLLG,DeltaPhiLLG,DeltaRLL,DeltaRLLG,ST,STG,STMET,ZPT,MTLL,MTLLG};
enum Histograms2D{ISRVFSR};
enum particleType{E,M};




map<Histograms1D,string> histoNames;
map<Histograms2D,string> histoNames2D;
void setHistoNames(){
histoNames[PT1]= "pt1";
histoNames[PT2]= "pt2";
histoNames[ETA1]= "eta1";
histoNames[ETA2]= "eta2";
histoNames[PHI1]= "phi1";
histoNames[PHI2]= "phi2";
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
histoNames[DeltaEtaLL] = "deltaEtaLL";
histoNames[DeltaPhiLL] = "deltaPhiLL";
histoNames[DeltaEtaLLG] = "deltaEtaLLG";
histoNames[DeltaPhiLLG] = "deltaPhiLLG";
histoNames[DeltaRLL] = "deltaRLL";
histoNames[DeltaRLLG] = "deltaRLLG";
histoNames[ST] = "st";
histoNames[STG] = "stg";
histoNames[STMET] = "stmet";
histoNames[ZPT] = "zpt";
histoNames[MTLL] = "mtll";
histoNames[MTLLG] = "mtllg";
};


