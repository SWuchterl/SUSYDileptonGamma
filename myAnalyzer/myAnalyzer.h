#include <math.h>
#include <regex>
#include <time.h>

#include "TROOT.h"
#include "TFile.h"
#include "TSelector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TEfficiency.h"
#include "TRandom2.h"
#include "TTree.h"

#include "TLorentzVector.h"

#include "TreeParticles.hpp"
#include "particleClasses.hpp"
#include "UserFunctions.h"
#include "Weighter.h"
//#include "CutFlow.h"
//#include "Resolution.h"


#include <iostream>

#include "config.h"

#include <chrono> //for sleeping
#include <thread> // --do--
#include <cstdlib>//for random increments 
#include <ctime>// --do--

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>


using namespace std;




typedef pair<UShort_t,UShort_t> SignalPoint;






class myAnalyzer : public TSelector {
 public:

  myAnalyzer();
  virtual ~myAnalyzer() { }

    //Functions
  virtual void Init(TTree *tree);
  virtual void SlaveBegin(TTree *tree);
  virtual Bool_t Process(Long64_t entry);
  virtual void Terminate();
  virtual Int_t Version() const { return 2; }
  
  //bool Cleaning();
  //void CalculateVariables(const tree::Lepton& l1, const tree::Lepton& l2,const particleType particle);
  //void ClearVariables();
  bool CheckParticles();
  bool Check2Ele();
  bool Check2Mu();
  bool CheckEMu();
  
  //bool matchGenParticle(const tree::Particle& pa);
  
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

  //bool testSelection(const tree::Electron& pa, selectionType,bool leading);
  //bool testSelection(const tree::Muon& pa, selectionType,bool leading);
  bool testSelection(const selElectron& pa, selectionType,bool leading);
  bool testSelection(const selMuon& pa, selectionType,bool leading);
  bool testSelection(const selPhoton& pa, selectionType);
  bool testSelection(const selMuon& pa, selectionType);
  bool testSelection(const selElectron& pa, selectionType);
  bool testSelection(const selJet& pa);
  
  
  
  //void Filler(selEvent& ev, map<Histograms1D,TH1F>& m,bool withPhoton);
  //tuple<int,string>
  //void Filler(selEvent& ev, map<Histograms1D,TH1F>& m,bool withPhoton,bool slimmed=false,int changePDF=9999,string changeMET="N");
  void Filler(map<Histograms1D,TH1F>& m,bool withPhoton,bool slimmed=false,int changePDF=9999,changemet changeMET=normal,
    changepu changePU=normalPU,
    changeLEPSF changeLepSF=normalLEPSF,
    changePHOTONSF changePhotonSF=normalPHOTONSF,
    changeISR changeisr=normalISR,
    changeEWK changeewk=normalEWK);
  //void FillerZZ(selEvent& ev, map<Histograms1D,TH1F>& m,bool withPhoton,selLepton& l1, selLepton& l2, selLepton& l3, selLepton& l4);
  //void FillerWZ(selEvent& ev, map<Histograms1D,TH1F>& m,bool withPhoton,selLepton& l1, selLepton& l2, selLepton& l3);
  //void FillerZZ(selEvent& ev, map<Histograms1D,TH1F>& m,bool withPhoton,selLepton& l1, selLepton& l2, selLepton& l3, selLepton& l4,bool slimmed=false,int changePDF=9999,string changeMET="N");
  //void FillerWZ(selEvent& ev, map<Histograms1D,TH1F>& m,bool withPhoton,selLepton& l1, selLepton& l2, selLepton& l3,bool slimmed=false,int changePDF=9999,string changeMET="N");
  //void Filler2D(selEvent& ev, map<Histograms2D,TH2F>& m,bool withPhoton);
  //void FillerSignal(selEvent& ev, map<Histograms1D,TH1F>& m,float divideFactor=1.);
  void FillerZZ(map<Histograms1D,TH1F>& m,bool withPhoton,selLepton& le1, selLepton& le2, selLepton& le3, selLepton& le4,bool slimmed=false,int changePDF=9999,changemet changeMET=normal,
    changepu changePU=normalPU,
    changeLEPSF changeLepSF=normalLEPSF,
    changePHOTONSF changePhotonSF=normalPHOTONSF,
    changeISR changeisr=normalISR,
    changeEWK changeewk=normalEWK);
  void FillerWZ(map<Histograms1D,TH1F>& m,bool withPhoton,selLepton& le1, selLepton& le2, selLepton& le3,bool slimmed=false,int changePDF=9999,changemet changeMET=normal,
    changepu changePU=normalPU,
    changeLEPSF changeLepSF=normalLEPSF,
    changePHOTONSF changePhotonSF=normalPHOTONSF,
    changeISR changeisr=normalISR,
    changeEWK changeewk=normalEWK);
  void Filler2D(map<Histograms2D,TH2F>& m,bool withPhoton);
  //void FillerSignal(map<Histograms1D,TH1F>& m,float divideFactor=1.);
  void FillerSignal(map<Histograms1D,TH1F>& m,float divideFactor=1.,int changePDF=9999,
    changemet changeMET=normal,
    changepu changePU=normalPU,
    changeLEPSF changeLepSF=normalLEPSF,
    changePHOTONSF changePhotonSF=normalPHOTONSF,
    changeISR changeisr=normalISR,
    changeEWK changeewk=normalEWK);
    
    
    //Weighter puWeighterUp;
    //Weighter puWeighterDown;
    
    
    
  //void FillerSignal(selEvent& ev, map<Histograms1D,TH1F>& m,float divideFactor);
  //void FillerTrigger(selEvent& ev, map<Histograms1D,TEfficiency>& m,bool withPhoton,bool TriggerBool);
  void FillerTrigger(map<Histograms1D,TEfficiency>& m,bool withPhoton,bool TriggerBool);
  
  void InitAllHistos();
  void InitTriggerStudiesHistos();
  void InitCutFlowHistos();
  void InitCutFlowHistos_Fine();
  //void InitSignalScanHistos(string masspoint);
  void InitSignalScanHistos(SignalPoint masspoint);
  

  

  
  map<Histograms1D,TH1F> InitHistograms(const selectionType selection);
  map<Histograms1D,TH1F> InitHistogramsSlimmed(const selectionType selection);
  map<Histograms2D,TH2F> Init2DHistograms(const selectionType selection);
  map<Histograms1D,TEfficiency> InitTriggerStudies(const selectionType selection);
  map<Histograms1D,TH1F> InitCutFlowHistograms(const selectionType selection);
  map<Histograms1D,TH1F> InitCutFlowHistograms_Fine(const selectionType selection);
  map<Histograms1D,TH1F> InitSignalScanHistograms(const selectionType selection);
  //void InitWeightHistos(map<string,map<string,map<Histograms1D,TH1F>>>& map_,const selectionType selection_, string name_);
  void InitWeightHistos(map<selectionFolderName,map<selectionFolderName,map<selectionFolderName,map<Histograms1D,TH1F>>>>& map_,const selectionType selection_, selectionFolderName name_);
  
  
  
  void FillHistograms();
  void FillHistograms2D();
  void FillTriggerStudies();
  void FillCutFlowHistograms();
  void SetCutFlowHistogramsStatus();
  
  void FillSignalHistograms();

  void FillCutFlowHistograms_Fine();
  


  // config.ini
  boost::property_tree::ptree propertyTree;
  map<selectionType,bool> config_selectionsToProcessMap;
  int config_veto;
  bool config_docutflow;
  bool config_docutflowfine;
  bool config_dosignalscan;
  bool config_dosignalscanSplit;
  float config_eventpercentage;
  bool config_doTrigger;
  string config_outputfolder;


    //Tree Variables
  TTreeReader fReader;
  
  TTreeReaderValue<std::vector<selPhoton>> selPhotons;
  TTreeReaderValue<std::vector<selJet>> selJets;
  TTreeReaderValue<std::vector<selMuon>> selMuons;
  TTreeReaderValue<std::vector<selElectron>> selElectrons;
  
  TTreeReaderValue<float> chargeProduct;
  TTreeReaderValue<float> deltaRll;
  TTreeReaderValue<float> mll;
  TTreeReaderValue<float> miniIso1;
  TTreeReaderValue<float> miniIso2;
  TTreeReaderValue<float> ETmiss;
  //TTreeReaderValue<TLorentzVector> ETmiss_vec;
  TTreeReaderValue<tree::Particle4Vector> ETmiss_vec;
  
  TTreeReaderValue<float> pt1;
  TTreeReaderValue<float> pt2;
  TTreeReaderValue<float> phi1;
  TTreeReaderValue<float> phi2;
  TTreeReaderValue<float> eta1;
  TTreeReaderValue<float> eta2;
  //TTreeReaderValue<TLorentzVector> l1;
  //TTreeReaderValue<TLorentzVector> l2;
  TTreeReaderValue<tree::Particle4Vector> l1;
  TTreeReaderValue<tree::Particle4Vector> l2;
  
  TTreeReaderValue<int> selMuonSize;
  TTreeReaderValue<int> selElectronSize;
  TTreeReaderValue<int> selLeptonSize;
  TTreeReaderValue<int> selPhotonSize;
  TTreeReaderValue<int> selJetSize;
  TTreeReaderValue<int> selBJetSize;
  
  TTreeReaderValue<int> matchedEleSize;
  TTreeReaderValue<int> matchedMuSize;
  TTreeReaderValue<int> matchedLeptonSize;
  
  TTreeReaderValue<bool> isDiElectron;
  TTreeReaderValue<bool> isDiMuon;
  TTreeReaderValue<bool> isMuonElectron;
  TTreeReaderValue<bool> isElectronMuon;
  
  TTreeReaderValue<bool> evtHasGenPhotonVeto;
  
  TTreeReaderValue<bool> trigHt;
  TTreeReaderValue<bool> trigDiEle;
  TTreeReaderValue<bool> trigDiMu;
  TTreeReaderValue<bool> trigMuEle;
  TTreeReaderValue<bool> trigDiEleMatch;
  TTreeReaderValue<bool> trigDiMuMatch;
  TTreeReaderValue<bool> trigMuEleMatch;
  
  TTreeReaderValue<float> calcHt;
  
  TTreeReaderValue<float> lepSF_weight;
  TTreeReaderValue<float> lepSF_weightUp;
  TTreeReaderValue<float> lepSF_weightDown;
  TTreeReaderValue<float> photonSF_weight;
  TTreeReaderValue<float> photonSF_weightUp;
  TTreeReaderValue<float> photonSF_weightDown;
  TTreeReaderValue<float> topPt_weight;
  TTreeReaderValue<float> isr_weight;
  TTreeReaderValue<float> isr_weightUp;
  TTreeReaderValue<float> isr_weightDown;
  TTreeReaderValue<float> ewk_weight;
  TTreeReaderValue<float> ewk_weightUp;
  TTreeReaderValue<float> ewk_weightDown;
  
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
  TTreeReaderValue<float> pu_weightUp;
  TTreeReaderValue<float> pu_weightDown;
  TTreeReaderValue<float> mc_weight;
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

  // signal scan
  TTreeReaderValue<UShort_t> nBinos;
  TTreeReaderValue<UShort_t> signal_m1;
  TTreeReaderValue<UShort_t> signal_m2;


  TTreeReaderValue<int> countNegCharge;
  TTreeReaderValue<int> countPosCharge;

  TTreeReaderValue<vector<selElectron>> negElectrons;
  TTreeReaderValue<vector<selElectron>> posElectrons;
  TTreeReaderValue<vector<selMuon>> negMuons;
  TTreeReaderValue<vector<selMuon>> posMuons;


  //int nEntries;
  long nEntries;
  SignalPoint sp_;

  //map<string,map<Histograms1D,TH1F>> h1Maps;
  //map<string,map<Histograms2D,TH2F>> h2Maps;
  
  //map<string,map<Histograms1D,TEfficiency>> eff1Maps;
  //map<string,map<Histograms2D,TEfficiency>> eff2Maps;
  
  //map<string,map<Histograms1D,TH1F>> c1Maps;


  //map<string,map<string,map<Histograms1D,TH1F>>> s1Maps;

  //CR
  //map<string,map<string,map<Histograms1D,TH1F>>> cr1Maps;
  map<selectionFolderName,map<selectionFolderName,map<Histograms1D,TH1F>>> h1Maps;
  map<selectionFolderName,map<selectionFolderName,map<Histograms2D,TH2F>>> h2Maps;
  
  map<selectionFolderName,map<selectionFolderName,map<Histograms1D,TEfficiency>>> eff1Maps;
  map<selectionFolderName,map<selectionFolderName,map<Histograms2D,TEfficiency>>> eff2Maps;
  
  //map<selectionFolderName,map<selectionFolderName,map<Histograms1D,TH1F>>> c1Maps;
  map<string,map<Histograms1D,TH1F>> c1Maps;



  //map<selectionFolderName,map<selectionFolderName,map<string,map<Histograms1D,TH1F>>>> s1Maps;
  //map<string,map<selectionFolderName,map<selectionFolderName,map<Histograms1D,TH1F>>>> s1Maps;
  //map<string,map<selectionFolderName,map<selectionFolderName,map<selectionFolderName,map<Histograms1D,TH1F>>>>> s1Maps;
  map<SignalPoint,map<selectionFolderName,map<selectionFolderName,map<selectionFolderName,map<Histograms1D,TH1F>>>>> s1Maps;

  //CR
  map<selectionFolderName,map<selectionFolderName,map<selectionFolderName,map<Histograms1D,TH1F>>>> cr1Maps;


  map <cutFlowFlags, bool> decisionMapCutFlowFine;
  map <cutFlowFlags, float> decisionMapCutFlowFine_weight;
  void clearCutFlowMap();
  
  TH1F cutFlow;
  float nGen;
  string inputName;


  int nWeights;

  
  float totalWeight;
  float totalWeightCalc;
  


  //cutflow booleans
  bool cutflowIsTriggered;
  bool cutflow2Leptons;
  bool cutflowMll50;
  bool cutflow1Photon;
  bool cutflowOnZ;
  bool cutflowDiEle;
  bool cutflowDiMu;

  bool doWeights_TopPt;
  bool doWeights_nISR;
  bool doWeights_EWKinoPairPt;
  bool doWeights_LeptonPairPt;

  bool isData;
  bool isSignal;
  bool isTotalSignal;
  

  double startTime;
  ClassDef(myAnalyzer, 1)
};

