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

#include <iostream>

#include "config.h"


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
      //hasGainSwitch=g.hasGainSwitch;
//
      //cIso=g.cIso;
      //nIso=g.nIso;
      //pIso=g.pIso;
      //cIsoWorst=g.cIsoWorst;
//
      //isTrue=g.isTrue;
      //isTrueAlternative=g.isTrueAlternative;
      //pMultifit=g.pMultifit;
      //pUncorrected=g.pUncorrected;

      isLoose=g.isLoose;
      isMedium=g.isMedium;
      isTight=g.isTight;

      isMediumMVA=g.isMediumMVA;
      mvaValue=g.mvaValue;
      //mvaCategory=g.mvaCategory;
   }
   TLorentzVector vec;
   float deltaR1;
   float deltaR2;
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
      //nhf=j.nhf;
      //cef=j.cef;
      //nef=j.nef;
      //nch=j.nch;
      //nconstituents=j.nconstituents;
      //ptRes=j.ptRes;
      //phiRes=j.phiRes;
      //sfRes=j.sfRes;
      //sfResUp=j.sfResUp;
      //sfResDn=j.sfResDn;
      //uncorJecFactor=j.uncorJecFactor; // uncorrected jet momentum over corrected jet momentum
   }
   TLorentzVector vec;
   float deltaR1;
   float deltaR2;
};


class selEvent {
    public:
    
    float totalWeight;
    
    bool isDiElectron;
    //TriggerDecisions(sum)
    bool trigDiEle;
    bool trigDiMu;
    bool trigMuEle;
    bool trigHt;

    //additional variables
    //leptons
    TLorentzVector l1;
    TLorentzVector l2;
    float pt1;  
    float pt2;  
    float phi1;  
    float phi2;  
    float eta1;  
    float eta2;
    //float charge;
    float chargeProduct;
    float deltaRll;  
    float deltaRl1e;  
    float deltaRl2e;  
    float mll;
    float miniIso1;
    float miniIso2;
    //photon
    //vector<tree::Photon> selPhotons;
    vector<selPhoton> selPhotons;
    vector<selJet> selJets;
    //MET
    float ETmiss;
};









class HistogramProducer : public TSelector {
 public:

  HistogramProducer();
  virtual ~HistogramProducer() { }

    //Functions
  virtual void Init(TTree *tree);
  virtual void SlaveBegin(TTree *tree);
  virtual Bool_t Process(Long64_t entry);
  virtual void Terminate();
  virtual Int_t Version() const { return 2; }
  
  //void Cleaning();
  bool Cleaning();
  bool CleaningTriggerStudies();
  void CalculateVariables(const tree::Lepton& l1, const tree::Lepton& l2, const tree::Photon& g,const tree::Particle met ,const particleType particle);
  void CalculateVariables(const tree::Lepton& l1, const tree::Lepton& l2,const tree::Particle met ,const particleType particle);
  bool CheckParticles();
  bool Check2Ele();
  bool Check2Mu();
  
  bool matchGenParticle(const tree::Particle& pa);
  
  bool SelectEvent(selectionType selection);
  bool SelectEventTriggerStudies(selectionType selection);

  bool testSelection(const tree::Electron& pa, selectionType,bool leading);
  bool testSelection(const tree::Muon& pa, selectionType,bool leading);
  bool testSelection(const selPhoton& pa, selectionType);
  bool testSelection(const selJet& pa, selectionType);
  
  
  void Filler(selEvent& ev, map<Histograms1D,TH1F>& m,bool withPhoton);
  void FillerTrigger(selEvent& ev, map<Histograms1D,TEfficiency>& m,bool withPhoton,bool TriggerBool);
  
  void InitAllHistos();
  void InitTriggerStudiesHistos();
  void InitScaleFactors();
  void InitScaleFactorsAlternative();
  float GetScaleFactorAndError(float pt, float eta,bool isFastSim, bool isEle);
  float GetScaleFactorAndErrorAlternative(float pt, float eta,bool isFastSim, bool isEle, int runNr);
  map<Histograms1D,TH1F> InitHistograms(const selectionType selection);
  map<Histograms2D,TH2F> Init2DHistograms(const selectionType selection);
  //void InitTriggerStudies();
  map<Histograms1D,TEfficiency> InitTriggerStudies(const selectionType selection);
  void FillHistograms();
  void FillHistograms2D();
  void FillTriggerStudies();

  bool GenPhotonVeto();


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
  //TTreeReaderValue<Bool_t> hlt_mu17_ele12_iso;
  //TTreeReaderValue<Bool_t> hlt_mu23_ele8_iso;
  //TTreeReaderValue<Bool_t> hlt_mu23_ele8_iso_dz;
  //TTreeReaderValue<Bool_t> hlt_mu23_ele12_iso;
  //TTreeReaderValue<Bool_t> hlt_mu23_ele12_iso_dz;
  //TTreeReaderValue<Bool_t> hlt_mu8_ele17_iso;
  //TTreeReaderValue<Bool_t> hlt_mu8_ele23_iso;
  //TTreeReaderValue<Bool_t> hlt_mu8_ele23_iso_dz;
  //TTreeReaderValue<Bool_t> hlt_mu12_ele23_iso;
  //TTreeReaderValue<Bool_t> hlt_mu12_ele23_iso_dz;
  TTreeReaderValue<Bool_t> hlt_doubleEle33;
  TTreeReaderValue<Bool_t> hlt_doubleEle33_mw;
  TTreeReaderValue<Bool_t> hlt_mu27_tkMu8;
  TTreeReaderValue<Bool_t> hlt_mu30_tkMu11;
  //TTreeReaderValue<Bool_t> hlt_mu30_ele30;
  //TTreeReaderValue<Bool_t> hlt_mu33_ele33;

  // signal scan
  TTreeReaderValue<UShort_t> signal_nBinos;
  TTreeReaderValue<UShort_t> signal_m1;
  TTreeReaderValue<UShort_t> signal_m2;

    //Selections
  vector<tree::Photon*> selPhotons;
  vector<tree::Jet*> selJets;
  vector<tree::Jet*> selBJets;
  vector<tree::Jet*> selHEJets;
  vector<tree::Electron*> selElectrons;
  vector<tree::Muon*> selMuons;

  vector<tree::Photon> artificialPhotons;

  int nEntries;


  map<string,map<Histograms1D,TH1F>> h1Maps;
  map<string,map<Histograms2D,TH2F>> h2Maps;
  
  map<string,map<Histograms1D,TEfficiency>> eff1Maps;
  map<string,map<Histograms2D,TEfficiency>> eff2Maps;
  
  TH1F cutFlow;
  string inputName;


  selEvent selectedEvent;



  //TriggerDecisions(sum)
  bool trigDiEle;
  bool trigDiMu;
  bool trigMuEle;
  bool trigHt;
  
  float totalWeight;
  
  //Decisions for Selection
  //bool l1_

  //additional variables
  bool isDiElectron;
  //leptons
  TLorentzVector lep1;
  TLorentzVector lep2;
  float pt1;  
  float pt2;  
  float phi1;  
  float phi2;  
  float eta1;  
  float eta2;
  float chargeProduct;
  float deltaRll;  
  float deltaRl1g;  
  float deltaRl2g;  
  float mll;
  float miniIso1;
  float miniIso2;
  //photon
  float pt;
  float eta;
  float phi;
  float miniIso;
  //MET
  float ETmiss;
  //float runNo;
    

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

  bool isData;
  bool isSignal;

  bool noPromptPhotons;

  double startTime;
  ClassDef(HistogramProducer, 1)
};

