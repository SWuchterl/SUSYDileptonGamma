#include <iostream>
using namespace std;



enum selectionType{UNCUT,DILEP,PHOTON,SEL,ONZ,TRIGDILEP,TRIGSEL,TRIGONZ,
    TRIGDILEP_ptcuts,TRIGONZ_ptcuts,TRIGSEL_ptcuts,TRIGDILEP_pt1cut,TRIGDILEP_pt2cut,
    ONZMET,ONZG,ABOVEZG,EXO,EGRegression,ONZMET100,ONZMET200,ONZMET100200,ONZMET200300,ONZMET100300,ONZMET0100,ONZMET100150,ONZMET150,
    ControlRegionDY,ControlRegionTT,ControlRegionTT080,ControlRegionTT80,ValidationRegion,ValidationRegion080,ValidationRegion80,ControlRegionZZ,LooseLeptons,ControlRegionWZ,ControlRegionWW};
    //CRDY,CRTT,CRTT080,CRTT80,VR,VR080,VR80,CRZZ,LooseLeptons,CRWZ,CRWW};
enum Histograms1D{PT1,PT2,ETA1,ETA2,PHI1,PHI2,MLL,NJETS,NPHOTONS,ETMISS,
    HT,GENHT,NVTX,ETAG1,PHIG1,PTG1,SIGMAIETAIETAG1,DeltaEtaLL,DeltaPhiLL,
    DeltaEtaLLG,DeltaPhiLLG,DeltaRLL,DeltaRLLG,ST,STG,STMET,ZPT,MTLL,MTLLG,
    CUTFLOW,CUTFLOW_fine,HOVERE,R9,SIGMAIPHIIPHIG1,DELTARGL1,DELTARGL2,MTL1MET,MTL2MET,MTGMET,MTLLMET,MTLLGMET,MT2,MOTHERID,MLLG,PT_llg,MZG_exo,gammaMotherID,genPhotonPT,
    genPhotonPT_Veto,PTG1_Veto,genPhotonPT_NoVeto,PTG1_NoVeto,VetoCompare,DeltaPhiLLMet,DeltaEtaLLMet,DeltaRLLMet,
    WEIGHT_TOPPT,WEIGHT_NISR,WEIGHT_EWKINOPAIRPT,WEIGHT_LEPTONPAIRPT,
    NElectrons,NMuons,
    PT3,PT4,ETA3,ETA4,PHI3,PHI4,MLL2,ZPT2,MTL3MET,MLLLL,
    JetPt1,JetPt2,JetPt3,JetPt4,JetPhi1,JetPhi2,JetPhi3,JetPhi4,JetEta1,JetEta2,JetEta3,JetEta4,
    DeltaEtaLL_neg,DeltaPhiLL_neg,DeltaRLL_neg};
enum Histograms2D{ISRVFSR,PTGvsMLLG};
enum particleType{E,M,DUMMYPARTICLE,EM,ME};

enum cutFlowFlags{TRIGGERED,LEPTONID_leading,LEPTONPT_leading,LEPTONID_trailing,LEPTONPT_trailing,PHOTON1,PHOTON1ID,PHOTON1PT,PHOTON1DR,M50,ZMASS,DIELECTRON,DIMUON,GENVETO,
    LEPTONIDPure_leading,LEPTONIDIso_leading,LEPTONIDImpact_leading,LEPTONIDEta_leading,LEPTONIDPure_trailing,LEPTONIDIso_trailing,LEPTONIDImpact_trailing,LEPTONIDEta_trailing,
    LEPTONIDDeltaR_leading,LEPTONIDDeltaR_trailing,genZLL,PHOTON1SEED,PHOTON1ETA,EMUON};


map<Histograms1D,string> histoNames;
map<Histograms2D,string> histoNames2D;
void setHistoNames(){
histoNames[PT1]= "pt1";
histoNames[PT2]= "pt2";
histoNames[PT3]= "pt3";
histoNames[PT4]= "pt4";
histoNames[ETA1]= "eta1";
histoNames[ETA2]= "eta2";
histoNames[ETA3]= "eta3";
histoNames[ETA4]= "eta4";
histoNames[PHI1]= "phi1";
histoNames[PHI2]= "phi2";
histoNames[PHI3]= "phi3";
histoNames[PHI4]= "phi4";
histoNames[MLL]= "m_ll";
histoNames[MLL2]= "m_ll2";
histoNames[MLLLL]= "m_llll";
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
histoNames[DeltaEtaLL_neg] = "deltaEtaLL_neg";
histoNames[DeltaPhiLL_neg] = "deltaPhiLL_neg";
histoNames[DeltaEtaLLG] = "deltaEtaLLG";
histoNames[DeltaPhiLLG] = "deltaPhiLLG";
histoNames[DeltaRLL] = "deltaRLL";
histoNames[DeltaRLL_neg] = "deltaRLL_neg";
histoNames[DeltaRLLG] = "deltaRLLG";
histoNames[ST] = "st";
histoNames[STG] = "stg";
histoNames[STMET] = "stmet";
histoNames[ZPT] = "zpt";
histoNames[ZPT2] = "zpt2";
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
histoNames[MZG_exo] = "mzg_exo";
histoNames[gammaMotherID] = "gammaMotherID";
histoNames[genPhotonPT] = "genPhotonPT";
histoNames[genPhotonPT_Veto] = "genPhotonPT_Veto";
histoNames[PTG1_Veto] = "PhotonPT_Veto";
histoNames[genPhotonPT_NoVeto] = "genPhotonPT_NoVeto";
histoNames[PTG1_NoVeto] = "PhotonPT_NoVeto";
histoNames[VetoCompare] = "VetoCompare";
histoNames[DeltaPhiLLMet] = "DeltaPhiLLMet";
histoNames[DeltaEtaLLMet] = "DeltaEtaLLMet";
histoNames[DeltaRLLMet] = "DeltaRLLMet";
histoNames[JetPt1] = "jetPt1";
histoNames[JetPt2] = "jetPt2";
histoNames[JetPt3] = "jetPt3";
histoNames[JetPt4] = "jetPt4";
histoNames[JetPhi1] = "jetPhi1";
histoNames[JetPhi2] = "jetPhi2";
histoNames[JetPhi3] = "jetPhi3";
histoNames[JetPhi4] = "jetPhi4";
histoNames[JetEta1] = "jetEta1";
histoNames[JetEta2] = "jetEta2";
histoNames[JetEta3] = "jetEta3";
histoNames[JetEta4] = "jetEta4";



histoNames[WEIGHT_NISR] = "weight_nISR";
histoNames[WEIGHT_EWKINOPAIRPT] = "weight_EWKinoPairPt";
histoNames[WEIGHT_LEPTONPAIRPT] = "weight_leptonPairPt";
histoNames[WEIGHT_TOPPT] = "weight_topPt";

histoNames[NElectrons] = "nElectrons";
histoNames[NMuons] = "nMuons";
histoNames[MTL3MET] = "mTL3Met";

histoNames2D[PTGvsMLLG] = "ptg_mllg";
histoNames2D[ISRVFSR] = "ISRvFSR";
};


