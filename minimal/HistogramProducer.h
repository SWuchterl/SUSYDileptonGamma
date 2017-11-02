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
//#include "Weighter.h"
//#include "CutFlow.h"
//#include "Resolution.h"

#include <iostream>

using namespace std;





class selEvent {
    public:
    
    float totalWeight;
    
    bool isDiElectron;
    //TriggerDecisions(sum)
    bool trigDiEle;
    bool trigDiMu;
    bool trigMuEle;
    bool trigHt;
  
    //Decisions for Selection
    //bool l1_

    //additional variables
    //leptons
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
    vector<tree::Photon> selPhotons;
    //float pt;
    //float eta;
    //float phi;
    //float miniIso;
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
  bool Cleaning(string selection);
  void CalculateVariables(const tree::Lepton& l1, const tree::Lepton& l2, const tree::Photon& g,const tree::Particle met ,const string particle);
  void CalculateVariables(const tree::Lepton& l1, const tree::Lepton& l2,const tree::Particle met ,const string particle);
  bool CheckParticles(string selection);
  bool Check2Ele();
  bool Check2Mu();
  
  bool matchGenParticle(const tree::Particle& pa);
  
  bool SelectEvent(string selection);

  bool testSelection(const tree::Electron& pa);
  bool testSelection(const tree::Muon& pa);
  bool testSelection(const tree::Photon& pa);
  
  void InitHistograms();
  void FillHistograms();

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

  map<string,TH1F> h1Map;
  TH1F cutFlow;
  string inputName;

  //vector<selEvent> selectedEvents;
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
    

  bool isData;
  bool isSignal;

  double startTime;
  ClassDef(HistogramProducer, 1)
};

