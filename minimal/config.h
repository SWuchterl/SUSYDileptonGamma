#include <iostream>
using namespace std;



enum selectionType{UNCUT,DILEP,PHOTON,SEL,ONZ,TRIGDILEP,TRIGSEL,TRIGONZ,
    TRIGDILEP_ptcuts,TRIGONZ_ptcuts,TRIGSEL_ptcuts,TRIGDILEP_pt1cut,TRIGDILEP_pt2cut,ONZMET,ONZG,ABOVEZG,EXO,EGRegression};
enum Histograms1D{PT1,PT2,ETA1,ETA2,PHI1,PHI2,MLL,NJETS,NPHOTONS,ETMISS,
    HT,GENHT,NVTX,ETAG1,PHIG1,PTG1,SIGMAIETAIETAG1,DeltaEtaLL,DeltaPhiLL,
    DeltaEtaLLG,DeltaPhiLLG,DeltaRLL,DeltaRLLG,ST,STG,STMET,ZPT,MTLL,MTLLG,
    CUTFLOW,CUTFLOW_fine,HOVERE,R9,SIGMAIPHIIPHIG1,DELTARGL1,DELTARGL2,MTL1MET,MTL2MET,MTGMET,MTLLMET,MTLLGMET,MT2,MOTHERID,MLLG,PT_llg,MZG_exo,gammaMotherID,genPhotonPT};
enum Histograms2D{ISRVFSR,PTGvsMLLG};
enum particleType{E,M};

enum cutFlowFlags{TRIGGERED,LEPTONID_leading,LEPTONPT_leading,LEPTONID_trailing,LEPTONPT_trailing,PHOTON1,PHOTON1ID,PHOTON1PT,PHOTON1DR,M50,ZMASS,DIELECTRON,DIMUON,GENVETO};


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
histoNames[SIGMAIPHIIPHIG1]= "sigmaIphiIphi_g1";
histoNames[DELTARGL1]= "deltaR1_g1";
histoNames[DELTARGL2]= "deltaR2_g1";
histoNames[HOVERE]= "hOverE_g1";
histoNames[R9]= "r9_g1";
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
histoNames[MTL1MET] = "mtl1met";
histoNames[MTL2MET] = "mtl2met";
histoNames[MTGMET] = "mtgmet";
histoNames[MTLLMET] = "mtllmet";
histoNames[MTLLGMET] = "mtllgmet";
histoNames[CUTFLOW] = "cutflow";
histoNames[CUTFLOW_fine] = "cutflow_fine";
histoNames[MT2] = "mt2";
histoNames[MOTHERID] = "motherID";
histoNames[MLLG] = "m_llg";
histoNames[PT_llg] = "pt_llg";
//histoNames[PTGvsMLLG] = "ptg_mllg";
histoNames[MZG_exo] = "mzg_exo";
histoNames[gammaMotherID] = "gammaMotherID";
histoNames[genPhotonPT] = "genPhotonPT";

histoNames2D[PTGvsMLLG] = "ptg_mllg";
histoNames2D[ISRVFSR] = "ISRvFSR";
};


