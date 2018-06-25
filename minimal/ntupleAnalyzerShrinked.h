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

#include "TreeParticles.hpp"
#include "UserFunctions.h"
#include "Weighter.h"
//#include "CutFlow.h"
//#include "Resolution.h"

//#include "MT2Functor.h"

#include <iostream>

#include "config.h"

#include <chrono> //for sleeping
#include <thread> // --do--
#include <cstdlib>//for random increments 
#include <ctime>// --do--

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>


using namespace std;


struct selPhoton : public tree::Photon{
   public:
   void setAll(const tree::Photon& g){
      p=g.p;
      sigmaIetaIeta=g.sigmaIetaIeta; // full 5x5
      sigmaIphiIphi=g.sigmaIphiIphi;
      hOverE=g.hOverE;
      hasPixelSeed=g.hasPixelSeed;
      passElectronVeto=g.passElectronVeto;
      r9=g.r9;
      sigmaPt=g.sigmaPt;
      hasGainSwitch=g.hasGainSwitch;
//
      cIso=g.cIso;
      nIso=g.nIso;
      pIso=g.pIso;
      cIsoWorst=g.cIsoWorst;
//
      isTrue=g.isTrue;
      isTrueAlternative=g.isTrueAlternative;
      pMultifit=g.pMultifit;
      pUncorrected=g.pUncorrected;

      isLoose=g.isLoose;
      isMedium=g.isMedium;
      isTight=g.isTight;

      isMediumMVA=g.isMediumMVA;
      mvaValue=g.mvaValue;
      mvaCategory=g.mvaCategory;
   }
   TLorentzVector vec;
   float deltaR1;
   float deltaR2;
   bool matchedToPhoton;
   bool matchedToJet;
   bool matchedToElectron;
};

struct selElectron : public tree::Electron{
   public:
   void setAll(const tree::Electron& e){
     
      p=e.p;
      
      charge=e.charge; // +/- 1
      rIso=e.rIso;
      passImpactParameter=e.passImpactParameter;
      d0=e.d0;
      dZ=e.dZ;
      SIP3D=e.SIP3D;
      miniIso=e.miniIso;
      
      isVetoID=e.isVetoID;
      isLoose=e.isLoose;
      isMedium=e.isMedium;
      isTight=e.isTight;
      isMediumMVA=e.isMediumMVA;
      isTightMVA=e.isTightMVA;
      isTightMVASlope=e.isTightMVASlope;
      mvaValue=e.mvaValue;
      mvaCategory=e.mvaCategory;
      seedCrystalE=e.seedCrystalE;
      isPassConvVeto=e.isPassConvVeto;
      pUncorrected=e.pUncorrected;
   }
   TLorentzVector vec;
   bool matched=false;
   float deltaR1;
   float deltaR2;
   
   int chargeInt;
   
};
struct selMuon : public tree::Muon{
   public:
   void setAll(const tree::Muon& m){
     
      p=m.p;
      
      charge=m.charge; // +/- 1
      rIso=m.rIso;
      passImpactParameter=m.passImpactParameter;
      d0=m.d0;
      dZ=m.dZ;
      SIP3D=m.SIP3D;
      miniIso=m.miniIso;
      
      isTight=m.isTight;
      isMedium=m.isMedium;
      isMediumRun=m.isMediumRun;
      nTrkLayers=m.nTrkLayers;
   }
   TLorentzVector vec;
   bool matched=false;
   float deltaR1;
   float deltaR2;
   
   int chargeInt;
   
};
struct selLepton : public tree::Lepton{
   public:
   void setAll(const selMuon& m){
     
      p=m.p;
      
      charge=m.charge; // +/- 1
      rIso=m.rIso;
      passImpactParameter=m.passImpactParameter;
      d0=m.d0;
      dZ=m.dZ;
      SIP3D=m.SIP3D;
      miniIso=m.miniIso;
      
      vec=m.vec;
      matched=m.matched;
      deltaR1=m.deltaR1;
      deltaR2=m.deltaR2;
     
      chargeInt=m.chargeInt;
      
   }
   void setAll(const selElectron& e){
     
      p=e.p;
      
      charge=e.charge; // +/- 1
      rIso=e.rIso;
      passImpactParameter=e.passImpactParameter;
      d0=e.d0;
      dZ=e.dZ;
      SIP3D=e.SIP3D;
      miniIso=e.miniIso;
      
      vec=e.vec;
      matched=e.matched;
      deltaR1=e.deltaR1;
      deltaR2=e.deltaR2;
     
      chargeInt=e.chargeInt;
   }
   TLorentzVector vec;
   bool matched=false;
   float deltaR1;
   float deltaR2;
   
   int chargeInt;
   
};

struct selJet : public tree::Jet{
   public:
   void setAll(const tree::Jet& j){
      p=j.p;
      isLoose = j.isLoose;
      hasPhotonMatch= j.hasPhotonMatch;
      hasElectronMatch=j.hasElectronMatch;
      hasMuonMatch=j.hasMuonMatch;
      bDiscriminator=j.bDiscriminator;
      uncert=j.uncert;
      chf=j.chf;
      nhf=j.nhf;
      cef=j.cef;
      nef=j.nef;
      nch=j.nch;
      nconstituents=j.nconstituents;
      ptRes=j.ptRes;
      phiRes=j.phiRes;
      sfRes=j.sfRes;
      sfResUp=j.sfResUp;
      sfResDn=j.sfResDn;
      uncorJecFactor=j.uncorJecFactor; // uncorrected jet momentum over corrected jet momentum
   }
   TLorentzVector vec;
   float deltaR1;
   float deltaR2;
   bool bTag;
};


class selEvent {
    public:
    
    float totalWeight=1.;
    
    bool isDiElectron=false;
    bool isDiMuon=false;
    bool isMuonElectron=false;
    bool isElectronMuon=false;
    //TriggerDecisions(sum)
    bool trigDiEle=false;
    bool trigDiMu=false;
    bool trigMuEle=false;
    bool trigHt=false;
    //TriggerObjectMatchingResults
    bool trigDiEleMatch=false;
    bool trigDiMuMatch=false;
    bool trigMuEleMatch=false;

    //additional variables
    //leptons
    TLorentzVector l1;
    TLorentzVector l2;
    float pt1=-10000.;  
    float pt2=-10000.;  
    float phi1=-1000.;  
    float phi2=-1000;  
    float eta1=-1000.;  
    float eta2=-1000.;
    float miniIso1=1000.;  
    float miniIso2=1000.;
    //float charge;
    float chargeProduct;
    float deltaRll=0.;  
    float deltaRl1e=0.;  
    float deltaRl2e=0.;  
    float mll=-1000.0;
    
    //photon
    vector<selPhoton> selPhotons;
    vector<selJet> selJets;
    vector<selElectron> selElectrons;
    vector<selMuon> selMuons;
    //MET
    float ETmiss=-10000.;
    TLorentzVector ETmiss_vec;
    float MT2_val=-10000.;
    
    float calcHt=0.;
    
    bool evtHasGenPhotonVeto=false;
    
    int matchedEleSize=0;
    int matchedMuSize=0;
    int matchedLeptonSize=0;
    int selLeptonSize=0;
    int selMuonSize=0;
    int selElectronSize=0;
    
    int nselJets;
    int nselBJets;
    
};









class ntupleAnalyzerShrinked : public TSelector {
 public:

  ntupleAnalyzerShrinked();
  virtual ~ntupleAnalyzerShrinked() { }

    //Functions
  virtual void Init(TTree *tree);
  virtual void SlaveBegin(TTree *tree);
  virtual Bool_t Process(Long64_t entry);
  virtual void Terminate();
  virtual Int_t Version() const { return 2; }
  
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
  bool SelectEventZZ(selectionType selection);

  bool testSelection(const tree::Electron& pa, selectionType,bool leading);
  bool testSelection(const tree::Muon& pa, selectionType,bool leading);
  bool testSelection(const selPhoton& pa, selectionType);
  bool testSelection(const selMuon& pa, selectionType);
  bool testSelection(const selElectron& pa, selectionType);
  bool testSelection(const selJet& pa);
  
  
  
  //void Filler(selEvent& ev, map<Histograms1D,TH1F>& m,bool withPhoton);
  //tuple<int,string>
  void Filler(selEvent& ev, map<Histograms1D,TH1F>& m,bool withPhoton,bool slimmed=false,int changePDF=9999,string changeMET="N");
  //void FillerZZ(selEvent& ev, map<Histograms1D,TH1F>& m,bool withPhoton,selLepton& l1, selLepton& l2, selLepton& l3, selLepton& l4);
  //void FillerWZ(selEvent& ev, map<Histograms1D,TH1F>& m,bool withPhoton,selLepton& l1, selLepton& l2, selLepton& l3);
  void FillerZZ(selEvent& ev, map<Histograms1D,TH1F>& m,bool withPhoton,selLepton& l1, selLepton& l2, selLepton& l3, selLepton& l4,bool slimmed=false,int changePDF=9999,string changeMET="N");
  void FillerWZ(selEvent& ev, map<Histograms1D,TH1F>& m,bool withPhoton,selLepton& l1, selLepton& l2, selLepton& l3,bool slimmed=false,int changePDF=9999,string changeMET="N");
  void Filler2D(selEvent& ev, map<Histograms2D,TH2F>& m,bool withPhoton);
  void FillerSignal(selEvent& ev, map<Histograms1D,TH1F>& m,float divideFactor=1.);
  //void FillerSignal(selEvent& ev, map<Histograms1D,TH1F>& m,float divideFactor);
  void FillerTrigger(selEvent& ev, map<Histograms1D,TEfficiency>& m,bool withPhoton,bool TriggerBool);
  
  void InitAllHistos();
  void InitTriggerStudiesHistos();
  void InitCutFlowHistos();
  void InitCutFlowHistos_Fine();
  void InitSignalScanHistos(string masspoint);
  
  void InitScaleFactors();
  void InitScaleFactorsAlternative();
  void InitScaleFactorsFinal();
  
  float GetScaleFactorAndError(float pt, float eta,bool isFastSim, bool isEle);
  float GetScaleFactorAndErrorAlternative(float pt, float eta,bool isFastSim, bool isEle, int runNr);
  float GetScaleFactorAndErrorFinal(float pt, float eta,bool isFastSim, bool isEle);
  float GetScaleFactorAndErrorPhotons(vector<selPhoton>& vecGamma);
  
  map<Histograms1D,TH1F> InitHistograms(const selectionType selection);
  map<Histograms1D,TH1F> InitHistogramsSlimmed(const selectionType selection);
  map<Histograms2D,TH2F> Init2DHistograms(const selectionType selection);
  map<Histograms1D,TEfficiency> InitTriggerStudies(const selectionType selection);
  map<Histograms1D,TH1F> InitCutFlowHistograms(const selectionType selection);
  map<Histograms1D,TH1F> InitCutFlowHistograms_Fine(const selectionType selection);
  map<Histograms1D,TH1F> InitSignalScanHistograms(const selectionType selection);
  void InitWeightHistos(map<string,map<string,map<Histograms1D,TH1F>>>& map_,const selectionType selection_, string name_);
  
  
  void InitCompressedTree();
  
  void FillHistograms();
  void FillHistograms2D();
  void FillTriggerStudies();
  void FillCutFlowHistograms();
  void SetCutFlowHistogramsStatus();
  
  void FillSignalHistograms();

  void FillCutFlowHistograms_Fine();
  

  bool GenPhotonVeto(const int a);

  // config.ini
  boost::property_tree::ptree propertyTree;
  map<selectionType,bool> config_selectionsToProcessMap;
  int config_veto;
  bool config_docutflow;
  bool config_docutflowfine;
  bool config_dosignalscan;
  bool config_dosignalscanSplit;
  float config_eventpercentage;
  string config_outputfolder;


    //Tree Variables
  TTreeReader fReader;
  TTreeReaderValue<std::vector<tree::Photon>> photons;
  TTreeReaderValue<std::vector<tree::Jet>> jets;
  TTreeReaderValue<std::vector<tree::Electron>> electrons;
  TTreeReaderValue<std::vector<tree::Muon>> muons;
  TTreeReaderValue<std::vector<tree::Particle>> genJets;
  TTreeReaderValue<std::vector<tree::GenParticle>> genParticles;
  TTreeReaderValue<std::vector<tree::IntermediateGenParticle>> intermediateGenParticles;
  TTreeReaderValue<tree::MET> met;
  TTreeReaderValue<tree::MET> metRaw;
  TTreeReaderValue<tree::MET> met_JESu;
  TTreeReaderValue<tree::MET> met_JESd;
  TTreeReaderValue<tree::MET> met_JERu;
  TTreeReaderValue<tree::MET> met_JERd;
  TTreeReaderValue<Float_t> pu_weight;
  TTreeReaderValue<Char_t> mc_weight;
  TTreeReaderValue<std::vector<Float_t>> pdf_weights;
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
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_ele17_ele12;
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_ele23_ele12;
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_ele33_ele33;
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_ele33_ele33_mw;
  
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu17_mu8;
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu17_mu8tk;
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu17_mu8_dz;
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu17_mu8tk_dz;
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu17tk_mu8tk_dz;
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu27_mu8;
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu30_mu11;

  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu17_ele12_eleLeg;
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu17_ele12_muLeg;
  
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu23_ele8_eleLeg;
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu23_ele8_muLeg;
  
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu23_ele8_dz;
  
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu23_ele12_eleLeg;
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu23_ele12_muLeg;
  
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu23_ele12_dz;
  
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu8_ele17_eleLeg;
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu8_ele17_muLeg;
  
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu8_ele23_eleLeg;
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu8_ele23_muLeg;
  
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu8_ele23_dz;
  
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu12_ele23_muLeg;
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu12_ele23_eleLeg;
    
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu12_ele23_dz;
  
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu30_ele30_eleLeg;
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu30_ele30_muLeg;
  
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu33_ele33_eleLeg;
  TTreeReaderValue<std::vector<tree::Particle>> trigObj_mu33_ele33_muLeg;


  // signal scan
  TTreeReaderValue<UShort_t> nBinos;
  TTreeReaderValue<UShort_t> signal_m1;
  TTreeReaderValue<UShort_t> signal_m2;

    //Selections
  vector<tree::Photon*> selPhotons;
  vector<tree::Jet*> selJets;
  vector<tree::Jet*> selBJets;
  vector<tree::Jet*> selHEJets;
  vector<tree::Electron*> selElectrons;
  vector<tree::Muon*> selMuons;
  
  
  
  //vector<tree::Muon> myMuons;

  int nEntries;


  map<string,map<Histograms1D,TH1F>> h1Maps;
  map<string,map<Histograms2D,TH2F>> h2Maps;
  
  map<string,map<Histograms1D,TEfficiency>> eff1Maps;
  map<string,map<Histograms2D,TEfficiency>> eff2Maps;
  
  map<string,map<Histograms1D,TH1F>> c1Maps;


  //map<string,map<Histograms1D,TH1F>> s1Maps;
  map<string,map<string,map<Histograms1D,TH1F>>> s1Maps;

  //CR
  map<string,map<string,map<Histograms1D,TH1F>>> cr1Maps;


  map <cutFlowFlags, bool> decisionMapCutFlowFine;
  map <cutFlowFlags, float> decisionMapCutFlowFine_weight;
  void clearCutFlowMap();
  
  TH1F cutFlow;
  float nGen;
  string inputName;


  selEvent selectedEvent;

  bool cleaned=0;


  int nWeights;

  //TriggerDecisions(sum)
  //bool trigDiEle;
  //bool trigDiMu;
  //bool trigMuEle;
  //bool trigHt;
  
  float totalWeight;
  

  //additional variables
  //bool isDiElectron;
  //bool isDiMuon;
  //bool isMuonElectron;
  //bool isElectronMuon;
  //leptons
  //TLorentzVector lep1;
  //TLorentzVector lep2;
  //float pt1;  
  //float pt2;  
  //float pt1UnCor;  
  //float pt2UnCor;  
  //float eta1UnCor;  
  //float eta2UnCor;  
  //float phi1UnCor;  
  //float phi2UnCor;  
  //float phi1;  
  //float phi2;  
  //float eta1;  
  //float eta2;
  //float chargeProduct;
  //float deltaRll;  
  //float deltaRl1g;  
  //float deltaRl2g;  
  //float mll;
  //float miniIso1;
  //float miniIso2;
  //photon
  //float pt;
  //float eta;
  //float phi;
  //float miniIso;
  //MET
  //float ETmiss;
  //TLorentzVector ETmiss_vec;
  //float runNo;
  //float MT2_val;
  //MT2Functor fctMT2_;

  //bool evtHasGenPhotonVeto;

  //cutflow booleans
  bool cutflowIsTriggered;
  bool cutflow2Leptons;
  bool cutflowMll50;
  bool cutflow1Photon;
  bool cutflowOnZ;
  bool cutflowDiEle;
  bool cutflowDiMu;


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
  Weighter FastSimDiEleWeighterID;
  Weighter FastSimDiEleWeighterIso;
  Weighter FastSimDiEleWeighterConv;
  Weighter FastSimDiMuWeighterID;
  Weighter FastSimDiMuWeighterIso;
  Weighter FastSimDiMuWeighterIP2D;
  Weighter FastSimDiMuWeighterSIP3D;

  Weighter PhotonIDWeighter;
  Weighter PhotonConversionWeighter;

  Weighter electronMllWeighter;

  bool isData;
  bool isSignal;
  bool isTotalSignal;
  
  bool doWeights_TopPt;
  bool doWeights_nISR;
  bool doWeights_EWKinoPairPt;
  bool doWeights_LeptonPairPt;

  bool noPromptPhotons;
  bool isZGammaInclusive;

  double startTime;
  ClassDef(ntupleAnalyzerShrinked, 1)
};

