#include <math.h>
#include <regex>
#include <time.h>

#include "TROOT.h"
#include "TFile.h"
#include "TSelector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TEfficiency.h"
#include "TRandom2.h"
#include "TTree.h"

#include "TKey.h"
#include "TList.h"
//#include "TIter.h"
#include "TObject.h"

#include "TreeParticles.hpp"
#include "UserFunctions.h"
#include "Weighter.h"
//#include "CutFlow.h"
//#include "Resolution.h"

#include "MT2Functor.h"

#include <iostream>

#include "config.h"

#include <chrono> //for sleeping
#include <thread> // --do--
#include <cstdlib>//for random increments
#include <ctime>// --do--

#include <boost/algorithm/string.hpp>


#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "particleClasses.hpp"

using namespace std;




class treeSplitter : public TSelector {
public:

treeSplitter();
virtual ~treeSplitter() {
}

//Functions
virtual void Init(TTree *tree);
virtual void SlaveBegin(TTree *tree);
virtual Bool_t Process(Long64_t entry);
virtual void Terminate();
virtual Int_t Version() const {
        return 2;
}

bool Cleaning();
bool CleaningTriggerStudies();
//void CalculateVariables(const tree::Lepton& l1, const tree::Lepton& l2, const tree::Photon& g,const tree::Particle met ,const particleType particle);
//void CalculateVariables(const tree::Lepton& l1, const tree::Lepton& l2,const tree::Particle met ,const particleType particle);
//void CalculateVariables(const tree::Lepton& l1, const tree::Lepton& l2, const tree::Photon& g,const tree::Particle met ,const particleType particle);
void CalculateVariables(const tree::Lepton& l1, const tree::Lepton& l2,const particleType particle);
void ClearVariables();
bool CheckParticles();
bool Check2Ele();
bool Check2Mu();
bool CheckEMu();

bool matchGenParticle(const tree::Particle& pa);

bool FindGenPhotonMatch(const selPhoton& pa);
tree::GenParticle& GetGenPhotonMatch(const selPhoton& pa);
bool FindGenPhotonMatch(const tree::Photon& pa);
tree::GenParticle& GetGenPhotonMatch(const tree::Photon& pa);
bool FindRecoPhotonMatch(const tree::GenParticle pa);

bool matchSelMuon(const selMuon& pa);
bool matchSelElectron(const selElectron& pa);
//bool matchLepton(const tree::Lepton& pa, selEvent& ev, bool isEle);
bool matchLepton(const selElectron& pa);
bool matchLepton(const selMuon& pa);

bool matchRecoPhoton(const selPhoton& pho,const string collection);

bool SelectEvent(selectionType selection);
bool SelectEventTriggerStudies(selectionType selection);

bool testSelection(const tree::Electron & pa, selectionType,bool leading);
bool testSelection(const tree::Muon & pa, selectionType,bool leading);
bool testSelection(const selPhoton & pa, selectionType);
bool testSelection(const selMuon & pa, selectionType);
bool testSelection(const selElectron & pa, selectionType);
bool testSelection(const selJet& pa);







void InitScaleFactors();
void InitScaleFactorsAlternative();
void InitScaleFactorsFinal();
void InitScaleFactorsFinalMedium();

float GetScaleFactorAndError(float pt, float eta,bool isFastSim, bool isEle);
float GetScaleFactorAndErrorAlternative(float pt, float eta,bool isFastSim, bool isEle, int runNr);
float GetScaleFactorAndErrorFinal(float pt, float eta,bool isFastSim, bool isEle,bool error=false);
//float GetScaleFactorAndErrorFinalMedium(float pt, float eta,bool isFastSim, bool isEle,bool error=false);
float GetScaleFactorAndErrorPhotons(vector<selPhoton>& vecGamma,bool error=false);




//void InitCompressedTree();
void InitTree();
void InitHTTree();
void SaveTree();
void SaveHTTree();


bool trigMET;

bool GenPhotonVeto(const int a);

// config.ini
boost::property_tree::ptree propertyTree;
int config_veto;
float config_eventpercentage;
bool config_doHt;
bool config_doMET;
bool config_doHtPure;
string config_outputfolder;


//Tree Variables
TTreeReader fReader;
TTreeReaderValue<std::vector<tree::Photon> > photons;
TTreeReaderValue<std::vector<tree::Jet> > jets;
TTreeReaderValue<std::vector<tree::Electron> > electrons;
TTreeReaderValue<std::vector<tree::Muon> > muons;
TTreeReaderValue<std::vector<tree::Particle> > genJets;
TTreeReaderValue<std::vector<tree::GenParticle> > genParticles;
TTreeReaderValue<std::vector<tree::IntermediateGenParticle> > intermediateGenParticles;
TTreeReaderValue<tree::MET> met;
TTreeReaderValue<tree::MET> met_gen;
TTreeReaderValue<tree::MET> metRaw;
TTreeReaderValue<tree::MET> met_JESu;
TTreeReaderValue<tree::MET> met_JESd;
TTreeReaderValue<tree::MET> met_JERu;
TTreeReaderValue<tree::MET> met_JERd;
TTreeReaderValue<Float_t> pu_weight;
TTreeReaderValue<Char_t> mc_weight;
TTreeReaderValue<std::vector<Float_t> > pdf_weights;
TTreeReaderValue<Int_t> nGoodVertices;
TTreeReaderValue<Int_t> nTracksPV;
TTreeReaderValue<Float_t> genHt;
TTreeReaderValue<Float_t> ht;
TTreeReaderValue<Float_t> rho;
TTreeReaderValue<Int_t> nTruePV;
TTreeReaderValue<Int_t> nISR;

TTreeReaderValue<Float_t> EWKinoPairPt;
TTreeReaderValue<Float_t> leptonPairPt;
TTreeReaderValue<Float_t> topPt1;
TTreeReaderValue<Float_t> topPt2;

TTreeReaderValue<ULong64_t> evtNo;
TTreeReaderValue<UInt_t> runNo;
TTreeReaderValue<UInt_t> lumNo;


TTreeReaderValue<Bool_t> hlt_met110;
TTreeReaderValue<Bool_t> hlt_met120;
TTreeReaderValue<Bool_t> hlt_met170_Noise;
TTreeReaderValue<Bool_t> hlt_met170_HBHE;
TTreeReaderValue<Bool_t> hlt_met170_Jet;
TTreeReaderValue<Bool_t> hlt_met170_Not;
TTreeReaderValue<Bool_t> hlt_met300;
TTreeReaderValue<Bool_t> hlt_met400;
TTreeReaderValue<Bool_t> hlt_met500;
TTreeReaderValue<Bool_t> hlt_met600;


TTreeReaderValue<Bool_t> hlt_ht200;
TTreeReaderValue<Bool_t> hlt_ht250;
TTreeReaderValue<Bool_t> hlt_ht300;
TTreeReaderValue<Bool_t> hlt_ht350;
TTreeReaderValue<Bool_t> hlt_ht400;
TTreeReaderValue<Bool_t> hlt_ht475;
TTreeReaderValue<Bool_t> hlt_ht600;
TTreeReaderValue<Bool_t> hlt_ht650;
TTreeReaderValue<Bool_t> hlt_ht800;
TTreeReaderValue<Bool_t> hlt_ele17_ele12_iso;
TTreeReaderValue<Bool_t> hlt_ele23_ele12_iso;
TTreeReaderValue<Bool_t> hlt_mu17_mu8_iso;
TTreeReaderValue<Bool_t> hlt_mu17_tkMu8_iso;
TTreeReaderValue<Bool_t> hlt_mu17_mu8_iso_dz;
TTreeReaderValue<Bool_t> hlt_mu17_tkMu8_iso_dz;
TTreeReaderValue<Bool_t> hlt_tkMu17_tkMu8_iso_dz;

TTreeReaderValue<Bool_t> hlt_mu17_ele12_iso;
TTreeReaderValue<Bool_t> hlt_mu23_ele8_iso;
TTreeReaderValue<Bool_t> hlt_mu23_ele8_iso_dz;
TTreeReaderValue<Bool_t> hlt_mu23_ele12_iso;
TTreeReaderValue<Bool_t> hlt_mu23_ele12_iso_dz;
TTreeReaderValue<Bool_t> hlt_mu8_ele17_iso;
TTreeReaderValue<Bool_t> hlt_mu8_ele23_iso;
TTreeReaderValue<Bool_t> hlt_mu8_ele23_iso_dz;
TTreeReaderValue<Bool_t> hlt_mu12_ele23_iso;
TTreeReaderValue<Bool_t> hlt_mu12_ele23_iso_dz;

TTreeReaderValue<Bool_t> hlt_doubleEle33;
TTreeReaderValue<Bool_t> hlt_doubleEle33_mw;
TTreeReaderValue<Bool_t> hlt_mu27_tkMu8;
TTreeReaderValue<Bool_t> hlt_mu30_tkMu11;

TTreeReaderValue<Bool_t> hlt_mu30_ele30;
TTreeReaderValue<Bool_t> hlt_mu33_ele33;


//Trigger Objects
//EE
TTreeReaderValue<std::vector<tree::Particle> > trigObj_ele17_ele12;
TTreeReaderValue<std::vector<tree::Particle> > trigObj_ele23_ele12;
TTreeReaderValue<std::vector<tree::Particle> > trigObj_ele33_ele33;
TTreeReaderValue<std::vector<tree::Particle> > trigObj_ele33_ele33_mw;

TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu17_mu8;
TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu17_mu8tk;
TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu17_mu8_dz;
TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu17_mu8tk_dz;
TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu17tk_mu8tk_dz;
TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu27_mu8;
TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu30_mu11;

TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu17_ele12_eleLeg;
TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu17_ele12_muLeg;

TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu23_ele8_eleLeg;
TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu23_ele8_muLeg;

TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu23_ele8_dz;

TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu23_ele12_eleLeg;
TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu23_ele12_muLeg;

TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu23_ele12_dz;

TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu8_ele17_eleLeg;
TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu8_ele17_muLeg;

TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu8_ele23_eleLeg;
TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu8_ele23_muLeg;

TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu8_ele23_dz;

TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu12_ele23_muLeg;
TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu12_ele23_eleLeg;

TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu12_ele23_dz;

TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu30_ele30_eleLeg;
TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu30_ele30_muLeg;

TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu33_ele33_eleLeg;
TTreeReaderValue<std::vector<tree::Particle> > trigObj_mu33_ele33_muLeg;


// signal scan
TTreeReaderValue<UShort_t> nBinos;
TTreeReaderValue<UShort_t> signal_m1;
TTreeReaderValue<UShort_t> signal_m2;

//Selections
//vector<tree::Photon*> selPhotons;
//vector<tree::Jet*> selJets;
//vector<tree::Jet*> selBJets;
//vector<tree::Jet*> selHEJets;
//vector<tree::Electron*> selElectrons;
//vector<tree::Muon*> selMuons;




vector<tree::Photon> temp_photons;
vector<tree::Jet> temp_jets;
vector<tree::Electron> temp_electrons;
vector<tree::Muon> temp_muons;
vector<tree::Particle> temp_genJets;
vector<tree::GenParticle> temp_genParticles;
vector<tree::IntermediateGenParticle> temp_intermediateGenParticles;   //eeTree->Branch("intermediateGenParticles",&intermediateGenParticles);
tree::MET temp_met;   //eeTree->Branch("met",&met);
tree::MET temp_met_gen;   //eeTree->Branch("met",&met);
tree::MET temp_metRaw;   //eeTree->Branch("metRaw",&metRaw);
tree::MET temp_met_JESu;   //eeTree->Branch("metJESu",&met_JESu);
tree::MET temp_met_JESd;  //eeTree->Branch("metJESd",&met_JESd);
tree::MET temp_met_JERu;   //eeTree->Branch("metJERu",&met_JERu);
tree::MET temp_met_JERd;   //eeTree->Branch("metJERd",&met_JERd);
Int_t temp_nGoodVertices;  //eeTree->Branch("nGoodVertices",&nGoodVertices);
Int_t temp_nTracksPV;   //eeTree->Branch("nTracksPV",&nTracksPV);

Float_t temp_pu_weight;   //eeTree->Branch("pu_weight",&pu_weight);
Float_t temp_mc_weight;   //eeTree->Branch("mc_weight",&mc_weight);
vector<Float_t> temp_pdf_weights;   //eeTree->Branch("pdf_weights",&pdf_weights);

Float_t temp_genHt;   //eeTree->Branch("genHt",&genHt);
Float_t temp_ht;   //eeTree->Branch("ht",&ht);
Int_t temp_nTruePV;   //eeTree->Branch("nTruePV",&nTruePV);

Int_t temp_nISR;   //eeTree->Branch("nISR",&nISR);
Float_t temp_EWKinoPairPt;   //eeTree->Branch("EWKinoPairPt",&EWKinoPairPt);
Float_t temp_leptonPairPt;   //eeTree->Branch("leptonPairPt",&leptonPairPt);
Float_t temp_topPt1;   //eeTree->Branch("topPt1",&topPt1);
Float_t temp_topPt2;   //eeTree->Branch("topPt2",&topPt2);

UInt_t temp_runNo;   //eeTree->Branch("runNo",&runNo);
UInt_t temp_lumNo;   //eeTree->Branch("lumNo",&lumNo);
ULong64_t temp_evtNo;   //eeTree->Branch("evtNo",&evtNo);

UShort_t temp_nBinos;
UShort_t temp_signal_m1;
UShort_t temp_signal_m2;

float temp_pu_weightUp;
float temp_pu_weightDown;

vector<selMuon> temp_negMuons;
vector<selMuon> temp_posMuons;
vector<selElectron> temp_negElectrons;
vector<selElectron> temp_posElectrons;

int temp_countNegCharge;
int temp_countPosCharge;


Weighter puWeighterUp;
Weighter puWeighterDown;


int genZToLLCounter;

TTree* eeTree;
TTree* htTree;
//string outputPathCompressed;
string outputFilenameCompressed;

TFile* treeFile;
TFile* treeHTFile;

int nEntries;


//map<string,map<Histograms1D,TH1F>> h1Maps;
//map<string,map<Histograms2D,TH2F>> h2Maps;

//map<string,map<Histograms1D,TEfficiency>> eff1Maps;
//map<string,map<Histograms2D,TEfficiency>> eff2Maps;

//map<string,map<Histograms1D,TH1F>> c1Maps;


//map<string,map<Histograms1D,TH1F>> s1Maps;
//map<string,map<string,map<Histograms1D,TH1F>>> s1Maps;

//CR
//map<string,map<string,map<Histograms1D,TH1F>>> cr1Maps;


//map <cutFlowFlags, bool> decisionMapCutFlowFine;
//map <cutFlowFlags, float> decisionMapCutFlowFine_weight;
//void clearCutFlowMap();

TH1F cutFlow;
map<string,TH1F> cutFlowMap;
TH1F weightHisto;
map<string,TH1F> weightHistoMap;
float nGen;
string inputName;


selEvent selectedEvent;


//float countGen=0;
//float countReco=0;
//TH1F genHist=TH1F("gen", ";p_{T}^{#gamma} (GeV)", 5000, 0, 5000);
//TH1F recoHist=TH1F("reco", ";p_{T}^{#gamma} (GeV)", 5000, 0, 5000);

int nWeights;



float totalWeight;



//bool evtHasGenPhotonVeto;

//cutflow booleans
//bool cutflowIsTriggered;
//bool cutflow2Leptons;
//bool cutflowMll50;
//bool cutflow1Photon;
//bool cutflowOnZ;
//bool cutflowDiEle;
//bool cutflowDiMu;


//scale factors
Weighter DiEleWeighterID;
Weighter DiEleWeighterIso;
Weighter DiEleWeighterConv;
Weighter DiEleWeighterTrack;

Weighter DiMuWeighterID;
Weighter DiMuWeighterID_BCDEF;
Weighter DiMuWeighterID_GH;
Weighter DiMuWeighterIso;
Weighter DiMuWeighterIso_BCDEF;
Weighter DiMuWeighterIso_GH;
Weighter DiMuWeighterIP2D;
Weighter DiMuWeighterSIP3D;
Weighter DiMuWeighterTrack;
Weighter DiMuWeighterTrack_BCDEF;
Weighter DiMuWeighterTrack_GH;

//Weighter FastSimDiEleWeighterID;
//Weighter FastSimDiEleWeighterIso;
//Weighter FastSimDiEleWeighterConv;
//Weighter FastSimDiMuWeighterID;
//Weighter FastSimDiMuWeighterIso;
//Weighter FastSimDiMuWeighterIP2D;
//Weighter FastSimDiMuWeighterSIP3D;

Weighter PhotonIDWeighter;
Weighter PhotonConversionWeighter;

Weighter electronMllWeighter;

bool isData;
bool isSignal;
bool isTotalSignal;

MT2Functor fctMT2_;

bool doWeights_TopPt;
bool doWeights_nISR;
bool doWeights_EWKinoPairPt;
bool doWeights_LeptonPairPt;

bool noPromptPhotons;
bool isZGammaInclusive;

double startTime;
ClassDef(treeSplitter, 1)
};
