#include <iostream>
using namespace std;



enum selectionType{UNCUT=0,DILEP,PHOTON,SEL,ONZ,TRIGDILEP,TRIGSEL,TRIGONZ,
    TRIGDILEP_ptcuts,TRIGONZ_ptcuts,TRIGSEL_ptcuts,TRIGDILEP_pt1cut,TRIGDILEP_pt2cut,
    ONZMET,ONZG,ABOVEZG,EXO,EGRegression,ONZMET100,ONZMET200,ONZMET100200,ONZMET200300,ONZMET100300,ONZMET0100,ONZMET100150,ONZMET150,
    ControlRegionDY,ControlRegionTT,ControlRegionTT080,ControlRegionTT80,ValidationRegion,ValidationRegion080,ValidationRegion80,ControlRegionZZ,LooseLeptons,ControlRegionWZ,ControlRegionWW};
    //CRDY,CRTT,CRTT080,CRTT80,VR,VR080,VR80,CRZZ,LooseLeptons,CRWZ,CRWW};
enum Histograms1D{PT1=0,PT2,ETA1,ETA2,PHI1,PHI2,MLL,NJETS,NPHOTONS,ETMISS,
    HT,GENHT,NVTX,ETAG1,PHIG1,PTG1,SIGMAIETAIETAG1,DeltaEtaLL,DeltaPhiLL,
    DeltaEtaLLG,DeltaPhiLLG,DeltaRLL,DeltaRLLG,ST,STG,STMET,ZPT,MTLL,MTLLG,
    CUTFLOW,CUTFLOW_fine,HOVERE,R9,SIGMAIPHIIPHIG1,DELTARGL1,DELTARGL2,MTL1MET,MTL2MET,MTGMET,MTLLMET,MTLLGMET,MT2,MOTHERID,MLLG,PT_llg,MZG_exo,gammaMotherID,genPhotonPT,
    genPhotonPT_Veto,PTG1_Veto,genPhotonPT_NoVeto,PTG1_NoVeto,VetoCompare,DeltaPhiLLMet,DeltaEtaLLMet,DeltaRLLMet,
    WEIGHT_TOPPT,WEIGHT_NISR,WEIGHT_EWKINOPAIRPT,WEIGHT_LEPTONPAIRPT,WEIGHT_PDF,
    NElectrons,NMuons,
    PT3,PT4,ETA3,ETA4,PHI3,PHI4,MLL2,ZPT2,MTL3MET,MLLLL,
    JetPt1,JetPt2,JetPt3,JetPt4,JetPhi1,JetPhi2,JetPhi3,JetPhi4,JetEta1,JetEta2,JetEta3,JetEta4,
    DeltaEtaLL_neg,DeltaPhiLL_neg,DeltaRLL_neg,NBJETS,DeltaPhiGMet,DeltaRGMet,
    FakeElectron,FakeJet,FakePhoton,Fakes,
    DeltaPhiL1G,DeltaPhiL2G,ML1G,ML2G,justMCPUWeight,withNISRWeight,withTopPtWeight,withEWKWeight,withPDFWeight,ETMISSRAW};
enum Histograms2D{ISRVFSR=0,PTGvsMLLG,
    MetZpt,MetDeltaPhiLL,MetMllg,MetMt2};
enum particleType{E=0,M,DUMMYPARTICLE,ElMu,MuEl};

enum cutFlowFlags{TRIGGERED=0,TRIGGEREDMATCHED,LEPTONID_leading,LEPTONPT_leading,LEPTONID_trailing,LEPTONPT_trailing,PHOTON1,PHOTON1ID,PHOTON1PT,PHOTON1DR,M50,ZMASS,DIELECTRON,DIMUON,GENVETO,
    LEPTONIDPure_leading,LEPTONIDIso_leading,LEPTONIDImpact_leading,LEPTONIDEta_leading,LEPTONIDPure_trailing,LEPTONIDIso_trailing,LEPTONIDImpact_trailing,LEPTONIDEta_trailing,
    LEPTONIDDeltaR_leading,LEPTONIDDeltaR_trailing,genZLL,PHOTON1SEED,PHOTON1ETA,EMUON,ex2LEP,MT2Cut,MetCut};


enum cutFlowMapName{
    cutFlow_Fine_onZEE=0,
    cutFlow_Fine_onZMM,
    cutFlow_Fine_onZEM,
    cutFlow_onZEE,
    cutFlow_onZMM,
    cutFlow_onZEM
    };

enum selectionFolderName{ 
    dilep=0,DILEPEE,DILEPMM,DILEPEM,DILEPLL,
    photon,PHOTONEE,PHOTONMM,PHOTONEM,PHOTONLL,
    sel,SELEE,SELMM,SELEM,SELLL,
    onz,ONZEE,ONZMM,ONZEM,ONZLL,
    trigdilep,TRIGDILEPEE,TRIGDILEPMM,TRIGDILEPEM,TRIGDILEPLL,
    trigsel,TRIGSELEE,TRIGSELMM,TRIGSELEM,TRIGSELLL,
    trigonz,TRIGONZEE,TRIGONZMM,TRIGONZEM,TRIGONZLL,
    //TRIGDILEP,TRIGDILEPEE,TRIGDILEPMM,TRIGDILEPEM,TRIGDILEPLL,
    trigdilep_ptcuts,TRIGDILEP_ptcutsEE,TRIGDILEP_ptcutsMM,TRIGDILEP_ptcutsEM,TRIGDILEP_ptcutsLL,
    trigonz_ptcuts,TRIGONZ_ptcutsEE,TRIGONZ_ptcutsMM,TRIGONZ_ptcutsEM,TRIGONZ_ptcutsLL,
    trigsel_ptcuts,TRIGSEL_ptcutsEE,TRIGSEL_ptcutsMM,TRIGSEL_ptcutsEM,TRIGSEL_ptcutsLL,
    onzmet,ONZMETEE,ONZMETMM,ONZMETEM,ONZMETLL,
    onzg,ONZGEE,ONZGMM,ONZGEM,ONZGLL,
    abovezg,ABOVEZGEE,ABOVEZGMM,ABOVEZGEM,ABOVEZGLL,
    //EXO,EXOEE,EXOMM,EXOEM,EXOLL,
    //EGRegression,EGRegressionEE,EGRegressionMM,EGRegressionEM,EGRegressionLL,
    onzmet100,ONZMET100EE,ONZMET100MM,ONZMET100EM,ONZMET100LL,
    onzmet200,ONZMET200EE,ONZMET200MM,ONZMET200EM,ONZMET200LL,
    onzmet100200,ONZMET100200EE,ONZMET100200MM,ONZMET100200EM,ONZMET100200LL,
    onzmet200300,ONZMET200300EE,ONZMET200300MM,ONZMET200300EM,ONZMET200300LL,
    onzmet100300,ONZMET100300EE,ONZMET100300MM,ONZMET100300EM,ONZMET100300LL,
    onzmet0100,ONZMET0100EE,ONZMET0100MM,ONZMET0100EM,ONZMET0100LL,
    onzmet100150,ONZMET100150EE,ONZMET100150MM,ONZMET100150EM,ONZMET100150LL,
    onzmet150,ONZMET150EE,ONZMET150MM,ONZMET150EM,ONZMET150LL,
    controlregionDY,ControlRegionDYEE,ControlRegionDYMM,ControlRegionDYEM,ControlRegionDYLL,
    controlregionTT,ControlRegionTTEE,ControlRegionTTMM,ControlRegionTTEM,ControlRegionTTLL,
    controlregionTT080,ControlRegionTT080EE,ControlRegionTT080MM,ControlRegionTT080EM,ControlRegionTT080LL,
    controlregionTT80,ControlRegionTT80EE,ControlRegionTT80MM,ControlRegionTT80EM,ControlRegionTT80LL,
    validationregion,ValidationRegionEE,ValidationRegionMM,ValidationRegionEM,ValidationRegionLL,
    validationregion080,ValidationRegion080EE,ValidationRegion080MM,ValidationRegion080EM,ValidationRegion080LL,
    validationregion80,ValidationRegion80EE,ValidationRegion80MM,ValidationRegion80EM,ValidationRegion80LL,
    controlregionZZ,ControlRegionZZEE,ControlRegionZZMM,ControlRegionZZEM,ControlRegionZZLL,
    controlregionWZ,ControlRegionWZEE,ControlRegionWZMM,ControlRegionWZEM,ControlRegionWZLL,
    JESu,JESd,JERu,JERd,
    LEPSFUP,LEPSFDOWN,PHOTONSFUP,PHOTONSFDOWN,PUUP,PUDOWN,ISRUP,ISRDOWN,EWKUP,EWKDOWN,NOPUUP,NOPUDOWN,NOPU,
    EE,LL,MM,EM,nom,
    PDF0,PDF1,PDF2,PDF3,PDF4,PDF5,PDF6,PDF7,PDF8,PDF9,
    PDF10,PDF11,PDF12,PDF13,PDF14,PDF15,PDF16,PDF17,PDF18,PDF19,
    PDF20,PDF21,PDF22,PDF23,PDF24,PDF25,PDF26,PDF27,PDF28,PDF29,
    PDF30,PDF31,PDF32,PDF33,PDF34,PDF35,PDF36,PDF37,PDF38,PDF39,
    PDF40,PDF41,PDF42,PDF43,PDF44,PDF45,PDF46,PDF47,PDF48,PDF49,
    PDF50,PDF51,PDF52,PDF53,PDF54,PDF55,PDF56,PDF57,PDF58,PDF59,
    PDF60,PDF61,PDF62,PDF63,PDF64,PDF65,PDF66,PDF67,PDF68,PDF69,
    PDF70,PDF71,PDF72,PDF73,PDF74,PDF75,PDF76,PDF77,PDF78,PDF79,
    PDF80,PDF81,PDF82,PDF83,PDF84,PDF85,PDF86,PDF87,PDF88,PDF89,
    PDF90,PDF91,PDF92,PDF93,PDF94,PDF95,PDF96,PDF97,PDF98,PDF99,
    PDF100,PDF101,PDF102,PDF103,PDF104,PDF105,PDF106,PDF107,PDF108,PDF109,
    sig,sig_zz,sig_gz,sig_gg,sig80,sig080,GENMET,sigMt2,sig80Mt2,sig080Mt2
    //PDF110,PDF111,PDF112,PDF113,PDF114,PDF115,PDF116,PDF117,PDF118,PDF119,
    };

selectionFolderName PDFNAMES[110] = {
    PDF0,PDF1,PDF2,PDF3,PDF4,PDF5,PDF6,PDF7,PDF8,PDF9,
    PDF10,PDF11,PDF12,PDF13,PDF14,PDF15,PDF16,PDF17,PDF18,PDF19,
    PDF20,PDF21,PDF22,PDF23,PDF24,PDF25,PDF26,PDF27,PDF28,PDF29,
    PDF30,PDF31,PDF32,PDF33,PDF34,PDF35,PDF36,PDF37,PDF38,PDF39,
    PDF40,PDF41,PDF42,PDF43,PDF44,PDF45,PDF46,PDF47,PDF48,PDF49,
    PDF50,PDF51,PDF52,PDF53,PDF54,PDF55,PDF56,PDF57,PDF58,PDF59,
    PDF60,PDF61,PDF62,PDF63,PDF64,PDF65,PDF66,PDF67,PDF68,PDF69,
    PDF70,PDF71,PDF72,PDF73,PDF74,PDF75,PDF76,PDF77,PDF78,PDF79,
    PDF80,PDF81,PDF82,PDF83,PDF84,PDF85,PDF86,PDF87,PDF88,PDF89,
    PDF90,PDF91,PDF92,PDF93,PDF94,PDF95,PDF96,PDF97,PDF98,PDF99,
    PDF100,PDF101,PDF102,PDF103,PDF104,PDF105,PDF106,PDF107,PDF108,PDF109,
    };

enum changemet{normal=0,JESUP,JESDOWN,JERUP,JERDOWN,METGEN};
enum changepu{normalPU=0,upPU,downPU,noPU,noPUDown,noPUUp};
enum changeLEPSF{normalLEPSF=0,upLEPSF,downLEPSF};
enum changePHOTONSF{normalPHOTONSF=0,upPHOTONSF,downPHOTONSF};
enum changeISR{normalISR=0,upISR,downISR};
enum changeEWK{normalEWK=0,upEWK,downEWK};

map<cutFlowMapName,string> cutFlowMapNameString;

map<selectionFolderName,string> selectionFolderNameString;
void setFolderNames(){
    cutFlowMapNameString[cutFlow_Fine_onZEE]="cutFlow_Fine_onZEE";
    cutFlowMapNameString[cutFlow_Fine_onZMM]="cutFlow_Fine_onZMM";
    cutFlowMapNameString[cutFlow_Fine_onZEM]="cutFlow_Fine_onZEM";
    cutFlowMapNameString[cutFlow_onZEE]="cutFlow_onZEE";
    cutFlowMapNameString[cutFlow_onZMM]="cutFlow_onZMM";
    cutFlowMapNameString[cutFlow_onZEM]="cutFlow_onZEM";
    
selectionFolderNameString[dilep]="dilep"; selectionFolderNameString[DILEPEE]="dilepEE"; selectionFolderNameString[DILEPMM]="dilepMM"; selectionFolderNameString[DILEPEM]="dilepEM"; selectionFolderNameString[DILEPLL]="dilepLL";    
//selectionFolderNameString[PHOTON]="Dilep"; selectionFolderNameString[PHOTONEE]="DilepEE"; selectionFolderNameString[PHOTONMM]="DilepMM"; selectionFolderNameString[PHOTONEM]="EM"; selectionFolderNameString[PHOTONLL]="DilepLL";    
selectionFolderNameString[sel]="sel"; selectionFolderNameString[SELEE]="selEE"; selectionFolderNameString[SELMM]="selMM"; selectionFolderNameString[SELEM]="selEM"; selectionFolderNameString[SELLL]="selLL";    
selectionFolderNameString[onz]="onZ"; selectionFolderNameString[ONZEE]="onZEE"; selectionFolderNameString[ONZMM]="onZMM"; selectionFolderNameString[ONZEM]="onZEM"; selectionFolderNameString[ONZLL]="onZLL";    
selectionFolderNameString[trigdilep]="trigDilep"; selectionFolderNameString[TRIGDILEPEE]="trigDilepEE"; selectionFolderNameString[TRIGDILEPMM]="trigDilepMM"; selectionFolderNameString[TRIGDILEPEM]="trigDilepEM"; selectionFolderNameString[TRIGDILEPLL]="trigDilepLL";    
selectionFolderNameString[trigdilep_ptcuts]="trigDilep_ptcuts"; selectionFolderNameString[TRIGDILEP_ptcutsEE]="trigDilep_ptcutsEE"; selectionFolderNameString[TRIGDILEP_ptcutsMM]="trigDilep_ptcutsMM"; selectionFolderNameString[TRIGDILEP_ptcutsEM]="trigDilep_ptcutsEM"; selectionFolderNameString[TRIGDILEP_ptcutsLL]="trigDilep_ptcutsLL";    
selectionFolderNameString[onzmet0100]="onZMet0100"; selectionFolderNameString[ONZMET0100EE]="onZMet0100EE"; selectionFolderNameString[ONZMET0100MM]="onZMet0100MM"; selectionFolderNameString[ONZMET0100EM]="onZMet0100EM"; selectionFolderNameString[ONZMET0100LL]="onZMet0100LL";    
selectionFolderNameString[onzmet100150]="onZMet100150"; selectionFolderNameString[ONZMET100150EE]="onZMet100150EE"; selectionFolderNameString[ONZMET100150MM]="oonZMet100150MM"; selectionFolderNameString[ONZMET100150EM]="onZMet100150EM"; selectionFolderNameString[ONZMET100150LL]="onZMet100150LL";    
selectionFolderNameString[onzmet150]="onZMet150"; selectionFolderNameString[ONZMET150EE]="onZMet150EE"; selectionFolderNameString[ONZMET150MM]="onZMet150MM"; selectionFolderNameString[ONZMET150EM]="onZMet150EM"; selectionFolderNameString[ONZMET150LL]="onZMet150LL";    
selectionFolderNameString[controlregionDY]="CRDY"; selectionFolderNameString[ControlRegionDYEE]="CRDYEE"; selectionFolderNameString[ControlRegionDYMM]="CRDYMM"; selectionFolderNameString[ControlRegionDYEM]="CRDYEM"; selectionFolderNameString[ControlRegionDYLL]="CRDYLL";    
selectionFolderNameString[controlregionTT]="CRTT"; selectionFolderNameString[ControlRegionTTEE]="CRTTEE"; selectionFolderNameString[ControlRegionTTMM]="CRTTMM"; selectionFolderNameString[ControlRegionTTEM]="CRTTEM"; selectionFolderNameString[ControlRegionTTLL]="CRTTLL";    
selectionFolderNameString[controlregionTT080]="CRTT080"; selectionFolderNameString[ControlRegionTT080EE]="CRTT080EE"; selectionFolderNameString[ControlRegionTT080MM]="CRTT080MM"; selectionFolderNameString[ControlRegionTT080EM]="CRTT080EM"; selectionFolderNameString[ControlRegionTT080LL]="CRTT080LL";    
selectionFolderNameString[controlregionTT80]="CRTT80"; selectionFolderNameString[ControlRegionTT80EE]="CRTT80EE"; selectionFolderNameString[ControlRegionTT80MM]="CRTT80MM"; selectionFolderNameString[ControlRegionTT080EM]="CRTT80EM"; selectionFolderNameString[ControlRegionTT080LL]="CRTT80LL";    
selectionFolderNameString[validationregion]="VR"; selectionFolderNameString[ValidationRegionEE]="VREE"; selectionFolderNameString[ValidationRegionMM]="VRMM"; selectionFolderNameString[ValidationRegionEM]="VREM"; selectionFolderNameString[ValidationRegionLL]="VRLL";    
selectionFolderNameString[validationregion080]="VR080"; selectionFolderNameString[ValidationRegion080EE]="VR080EE"; selectionFolderNameString[ValidationRegion080MM]="VR080MM"; selectionFolderNameString[ValidationRegion080EM]="VR080EM"; selectionFolderNameString[ValidationRegion080LL]="VR080LL";    
selectionFolderNameString[validationregion80]="VR80"; selectionFolderNameString[ValidationRegion80EE]="VR80EE"; selectionFolderNameString[ValidationRegion80MM]="VR80MM"; selectionFolderNameString[ValidationRegion80EM]="VR80EM"; selectionFolderNameString[ValidationRegion80LL]="VR80LL";    
selectionFolderNameString[controlregionZZ]="CRZZ"; selectionFolderNameString[ControlRegionZZEE]="CRZZEE"; selectionFolderNameString[ControlRegionZZMM]="CRZZMM"; selectionFolderNameString[ControlRegionZZEM]="CRZZEM"; selectionFolderNameString[ControlRegionZZLL]="CRZZLL";    
selectionFolderNameString[controlregionWZ]="CRWZ"; selectionFolderNameString[ControlRegionWZEE]="CRWZEE"; selectionFolderNameString[ControlRegionWZMM]="CRWZMM"; selectionFolderNameString[ControlRegionWZEM]="CRWZEM"; selectionFolderNameString[ControlRegionWZLL]="CRWZLL";    
selectionFolderNameString[JESu]="JESu"; 
selectionFolderNameString[JESd]="JESd"; 
selectionFolderNameString[JERu]="JERu"; 
selectionFolderNameString[JERd]="JERd";
 
selectionFolderNameString[LEPSFUP]="lepSFu"; 
selectionFolderNameString[LEPSFDOWN]="lepSFd"; 
selectionFolderNameString[PHOTONSFUP]="photonSFu"; 
selectionFolderNameString[PHOTONSFDOWN]="photonSFd"; 
selectionFolderNameString[PUUP]="PUu"; 
selectionFolderNameString[NOPUUP]="NoPUu"; 
selectionFolderNameString[NOPU]="NoPU"; 
selectionFolderNameString[PUDOWN]="PUd"; 
selectionFolderNameString[NOPUDOWN]="NoPUd"; 
selectionFolderNameString[ISRUP]="ISRu"; 
selectionFolderNameString[ISRDOWN]="ISRd"; 
selectionFolderNameString[EWKUP]="EWKu"; 
selectionFolderNameString[EWKDOWN]="EWKd"; 
selectionFolderNameString[GENMET]="genmet"; 

selectionFolderNameString[EE]="EE"; selectionFolderNameString[MM]="MM"; selectionFolderNameString[EM]="EM"; selectionFolderNameString[LL]="LL";
selectionFolderNameString[nom]="nom"; 
selectionFolderNameString[PDF0]="0"; 
selectionFolderNameString[PDF1]="1"; 
selectionFolderNameString[PDF2]="2"; 
selectionFolderNameString[PDF3]="3"; 
selectionFolderNameString[PDF4]="4"; 
selectionFolderNameString[PDF5]="5"; 
selectionFolderNameString[PDF6]="6"; 
selectionFolderNameString[PDF7]="7"; 
selectionFolderNameString[PDF8]="8"; 
selectionFolderNameString[PDF9]="9"; 
selectionFolderNameString[PDF10]="10"; 
selectionFolderNameString[PDF11]="11"; 
selectionFolderNameString[PDF12]="12"; 
selectionFolderNameString[PDF13]="13"; 
selectionFolderNameString[PDF14]="14"; 
selectionFolderNameString[PDF15]="15"; 
selectionFolderNameString[PDF16]="16"; 
selectionFolderNameString[PDF17]="17"; 
selectionFolderNameString[PDF18]="18"; 
selectionFolderNameString[PDF19]="19"; 
selectionFolderNameString[PDF20]="20"; 
selectionFolderNameString[PDF21]="21"; 
selectionFolderNameString[PDF22]="22"; 
selectionFolderNameString[PDF23]="23"; 
selectionFolderNameString[PDF24]="24"; 
selectionFolderNameString[PDF25]="25"; 
selectionFolderNameString[PDF26]="26"; 
selectionFolderNameString[PDF27]="27"; 
selectionFolderNameString[PDF28]="28"; 
selectionFolderNameString[PDF29]="29"; 
selectionFolderNameString[PDF30]="30"; 
selectionFolderNameString[PDF31]="31"; 
selectionFolderNameString[PDF32]="32"; 
selectionFolderNameString[PDF33]="33"; 
selectionFolderNameString[PDF34]="34"; 
selectionFolderNameString[PDF35]="35"; 
selectionFolderNameString[PDF36]="36"; 
selectionFolderNameString[PDF37]="37"; 
selectionFolderNameString[PDF38]="38"; 
selectionFolderNameString[PDF39]="39"; 
selectionFolderNameString[PDF40]="40"; 
selectionFolderNameString[PDF41]="41"; 
selectionFolderNameString[PDF42]="42"; 
selectionFolderNameString[PDF43]="43"; 
selectionFolderNameString[PDF44]="44"; 
selectionFolderNameString[PDF45]="45"; 
selectionFolderNameString[PDF46]="46"; 
selectionFolderNameString[PDF47]="47"; 
selectionFolderNameString[PDF48]="48"; 
selectionFolderNameString[PDF49]="49"; 
selectionFolderNameString[PDF50]="50"; 
selectionFolderNameString[PDF51]="51"; 
selectionFolderNameString[PDF52]="52"; 
selectionFolderNameString[PDF53]="53"; 
selectionFolderNameString[PDF54]="54"; 
selectionFolderNameString[PDF55]="55"; 
selectionFolderNameString[PDF56]="56"; 
selectionFolderNameString[PDF57]="57"; 
selectionFolderNameString[PDF58]="58"; 
selectionFolderNameString[PDF59]="59"; 
selectionFolderNameString[PDF60]="60"; 
selectionFolderNameString[PDF61]="61"; 
selectionFolderNameString[PDF62]="62"; 
selectionFolderNameString[PDF63]="63"; 
selectionFolderNameString[PDF64]="64"; 
selectionFolderNameString[PDF65]="65"; 
selectionFolderNameString[PDF66]="66"; 
selectionFolderNameString[PDF67]="67"; 
selectionFolderNameString[PDF68]="68"; 
selectionFolderNameString[PDF69]="69"; 
selectionFolderNameString[PDF70]="70"; 
selectionFolderNameString[PDF71]="71"; 
selectionFolderNameString[PDF72]="72"; 
selectionFolderNameString[PDF73]="73"; 
selectionFolderNameString[PDF74]="74"; 
selectionFolderNameString[PDF75]="75"; 
selectionFolderNameString[PDF76]="76"; 
selectionFolderNameString[PDF77]="77"; 
selectionFolderNameString[PDF78]="78"; 
selectionFolderNameString[PDF79]="79"; 
selectionFolderNameString[PDF80]="80"; 
selectionFolderNameString[PDF81]="81"; 
selectionFolderNameString[PDF82]="82"; 
selectionFolderNameString[PDF83]="83"; 
selectionFolderNameString[PDF84]="84"; 
selectionFolderNameString[PDF85]="85"; 
selectionFolderNameString[PDF86]="86"; 
selectionFolderNameString[PDF87]="87"; 
selectionFolderNameString[PDF88]="88"; 
selectionFolderNameString[PDF89]="89"; 
selectionFolderNameString[PDF90]="90"; 
selectionFolderNameString[PDF91]="91"; 
selectionFolderNameString[PDF92]="92"; 
selectionFolderNameString[PDF93]="93"; 
selectionFolderNameString[PDF94]="94"; 
selectionFolderNameString[PDF95]="95"; 
selectionFolderNameString[PDF96]="96"; 
selectionFolderNameString[PDF97]="97"; 
selectionFolderNameString[PDF98]="98"; 
selectionFolderNameString[PDF99]="99"; 
selectionFolderNameString[PDF100]="100"; 
selectionFolderNameString[PDF101]="101"; 
selectionFolderNameString[PDF102]="102"; 
selectionFolderNameString[PDF103]="103"; 
selectionFolderNameString[PDF104]="104"; 
selectionFolderNameString[PDF105]="105"; 
selectionFolderNameString[PDF106]="106"; 
selectionFolderNameString[PDF107]="107"; 
selectionFolderNameString[PDF108]="108"; 
selectionFolderNameString[PDF109]="109"; 
selectionFolderNameString[sig]="sig"; 
selectionFolderNameString[sig80]="sig80"; 
selectionFolderNameString[sig080]="sig080"; 
selectionFolderNameString[sigMt2]="sigMt2"; 
selectionFolderNameString[sig80Mt2]="sig80Mt2"; 
selectionFolderNameString[sig080Mt2]="sig080Mt2"; 
selectionFolderNameString[sig_gg]="sig_gg"; 
selectionFolderNameString[sig_zz]="sig_zz"; 
selectionFolderNameString[sig_gz]="sig_gz"; 
}



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
histoNames[NBJETS]= "n_bjets";
histoNames[NPHOTONS]= "n_photons";
histoNames[ETMISS]= "met";
histoNames[ETMISSRAW]= "met_RAW";
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
histoNames[DeltaPhiGMet] = "DeltaPhiGMet";
histoNames[DeltaRGMet] = "DeltaRGMet";
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
histoNames[DeltaPhiL1G] = "DeltaPhiL1G";
histoNames[DeltaPhiL2G] = "DeltaPhiL2G";
histoNames[ML1G] = "m_l1g";
histoNames[ML2G] = "m_l2g";



histoNames[WEIGHT_NISR] = "weight_nISR";
histoNames[WEIGHT_EWKINOPAIRPT] = "weight_EWKinoPairPt";
histoNames[WEIGHT_LEPTONPAIRPT] = "weight_leptonPairPt";
histoNames[WEIGHT_TOPPT] = "weight_topPt";
histoNames[WEIGHT_PDF] = "weight_PDF";

histoNames[NElectrons] = "nElectrons";
histoNames[NMuons] = "nMuons";
histoNames[MTL3MET] = "mTL3Met";

histoNames[FakeElectron] = "FakeElectron";
histoNames[FakeJet] = "FakeJet";
histoNames[FakePhoton] = "FakePhoton";
histoNames[Fakes] = "Fakes";

histoNames2D[PTGvsMLLG] = "ptg_mllg";
histoNames2D[ISRVFSR] = "ISRvFSR";
histoNames2D[MetZpt] = "MetZpt";
histoNames2D[MetDeltaPhiLL] = "MetDeltaPhiLL";
histoNames2D[MetMllg] = "MetMllg";
histoNames2D[MetMt2] = "MetMt2";
};


