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

#include "TKey.h"
#include "TList.h"
//#include "TIter.h"
#include "TObject.h"

#include "TLorentzVector.h"

#include "TreeParticles.hpp"
#include "particleClasses.hpp"
#include "UserFunctions.h"
#include "Weighter.h"
//#include "CutFlow.h"
//#include "Resolution.h"

	
#include <boost/algorithm/string.hpp>
//#include <boost/algorithm/string/split.hpp>
#include <iostream>
//using namespace boost;
#include "config.h"

#include <chrono> //for sleeping
#include <thread> // --do--
#include <cstdlib>//for random increments 
#include <ctime>// --do--

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

//#include <unordered_map>

using namespace std;




typedef pair<UShort_t,UShort_t> SignalPoint;



//selectionType,Histograms1D,Histograms2D,cutFlowFlags,selectionFolderName,changemet,changepu,changeLEPSF,changePHOTONSF,changeISR,changeEWK,
//namespace std{
	//template <>
	//struct hash<selectionType>{
		//size_t operator()(const selectionType& v) const{
			//return hash<int>()(v);
		//}
	//};
	//template <>
	//struct hash<Histograms1D>{
		//size_t operator()(const Histograms1D& v) const{
			//return hash<int>()(v);
		//}
	//};
	//template <>
	//struct hash<Histograms2D>{
		//size_t operator()(const Histograms2D& v) const{
			//return hash<int>()(v);
		//}
	//};
	//template <>
	//struct hash<cutFlowFlags>{
		//size_t operator()(const cutFlowFlags& v) const{
			//return hash<int>()(v);
		//}
	//};
	//template <>
	//struct hash<selectionFolderName>{
		//size_t operator()(const selectionFolderName& v) const{
			//return hash<int>()(v);
		//}
	//};
	//template <>
	//struct hash<changemet>{
		//size_t operator()(const changemet& v) const{
			//return hash<int>()(v);
		//}
	//};
	//template <>
	//struct hash<changepu>{
		//size_t operator()(const changepu& v) const{
			//return hash<int>()(v);
		//}
	//};
	//template <>
	//struct hash<changeLEPSF>{
		//size_t operator()(const changeLEPSF& v) const{
			//return hash<int>()(v);
		//}
	//};
	//template <>
	//struct hash<changePHOTONSF>{
		//size_t operator()(const changePHOTONSF& v) const{
			//return hash<int>()(v);
		//}
	//};
	//template <>
	//struct hash<changeISR>{
		//size_t operator()(const changeISR& v) const{
			//return hash<int>()(v);
		//}
	//};
	//template <>
	//struct hash<changeEWK>{
		//size_t operator()(const changeEWK& v) const{
			//return hash<int>()(v);
		//}
	//};
	//template <>
	//struct hash<SignalPoint>{
		//size_t operator()(const SignalPoint& v) const{
        //size_t h1 = std::hash<unsigned short>()(v.first);
        //size_t h2 = std::hash<unsigned short>()(v.second);
        //return h1 ^ ( h2 << 1 );
		//}
	//};
//}




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
  

  bool CheckParticles();
  bool Check2Ele();
  bool Check2Mu();
  bool CheckEMu();
  
  
  bool FindGenPhotonMatch(const selPhoton& pa);
  tree::GenParticle& GetGenPhotonMatch(const selPhoton& pa);
  bool FindGenPhotonMatch(const tree::Photon& pa);
  tree::GenParticle& GetGenPhotonMatch(const tree::Photon& pa);
  bool FindRecoPhotonMatch(const tree::GenParticle pa);

  bool matchSelMuon(const selMuon& pa);
  bool matchSelElectron(const selElectron& pa);
  bool matchLepton(const selElectron& pa);
  bool matchLepton(const selMuon& pa);
  
  bool matchRecoPhoton(const selPhoton& pho,const string collection);
  
  bool SelectEvent(selectionType selection);
  bool SelectEventTriggerStudies(selectionType selection);
  bool SelectEventZZ(selectionType selection);

  bool testSelection(const selElectron& pa, selectionType,bool leading);
  bool testSelection(const selMuon& pa, selectionType,bool leading);
  bool testSelection(const selPhoton& pa, selectionType);
  bool testSelection(const selMuon& pa, selectionType);
  bool testSelection(const selElectron& pa, selectionType);
  bool testSelection(const selJet& pa);
  
  
  void Filler(map<Histograms1D,TH1F>& m,bool withPhoton,bool slimmed=false,int changePDF=9999,changemet changeMET=normal,
    changepu changePU=normalPU,
    changeLEPSF changeLepSF=normalLEPSF,
    changePHOTONSF changePhotonSF=normalPHOTONSF,
    changeISR changeisr=normalISR,
    changeEWK changeewk=normalEWK);
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
  //void Filler2D(map<Histograms2D,TH2F>& m,bool withPhoton);
  void Filler2D(map<Histograms2D,TH2F>& m,bool withPhoton,bool slimmed=false,int changePDF=9999,changemet changeMET=normal,
    changepu changePU=normalPU,
    changeLEPSF changeLepSF=normalLEPSF,
    changePHOTONSF changePhotonSF=normalPHOTONSF,
    changeISR changeisr=normalISR,
    changeEWK changeewk=normalEWK);
  //void FillerSignal(map<Histograms1D,TH1F>& m,float divideFactor=1.);
  void FillerSignal(map<Histograms1D,TH1F>& m,float divideFactor=1.,int changePDF=9999,
    changemet changeMET=normal,
    changepu changePU=normalPU,
    changeLEPSF changeLepSF=normalLEPSF,
    changePHOTONSF changePhotonSF=normalPHOTONSF,
    changeISR changeisr=normalISR,
    changeEWK changeewk=normalEWK);
  void FillerSignal(map<Histograms2D,TH2F>& m,float divideFactor=1.,int changePDF=9999,
    changemet changeMET=normal,
    changepu changePU=normalPU,
    changeLEPSF changeLepSF=normalLEPSF,
    changePHOTONSF changePhotonSF=normalPHOTONSF,
    changeISR changeisr=normalISR,
    changeEWK changeewk=normalEWK);
    
  void FillerTrigger(map<Histograms1D,TEfficiency>& m,bool withPhoton,bool TriggerBool);
  
  void InitAllHistos();
  void InitTriggerStudiesHistos();
  void InitCutFlowHistos();
  void InitCutFlowHistos_Fine();
  void InitSignalScanHistos(SignalPoint masspoint);
  

  

  
  map<Histograms1D,TH1F> InitHistograms(const selectionType selection);
  map<Histograms1D,TH1F> InitHistogramsSlimmed(const selectionType selection);
  map<Histograms2D,TH2F> Init2DHistograms(const selectionType selection);
  map<Histograms1D,TEfficiency> InitTriggerStudies(const selectionType selection);
  map<Histograms1D,TH1F> InitCutFlowHistograms(const selectionType selection);
  map<Histograms1D,TH1F> InitCutFlowHistograms_Fine(const selectionType selection);
  map<Histograms1D,TH1F> InitSignalScanHistograms(const selectionType selection);
  map<Histograms2D,TH2F> InitSignalScanHistograms2D(const selectionType selection);
  void InitWeightHistos(map<selectionFolderName,map<selectionFolderName,map<selectionFolderName,map<Histograms1D,TH1F>>>>& map_,const selectionType selection_, selectionFolderName name_);
  
  
  
  void FillHistograms();
  void FillHistograms2D();
  void FillTriggerStudies();
  void FillCutFlowHistograms();
  void SetCutFlowHistogramsStatus();
  
  void FillSignalHistograms();

  void FillCutFlowHistograms_Fine();
  
  bool GenPhotonTTVeto();
  float GetZZKFactor();

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
  bool config_doTriggerPure;
  bool config_doTriggerMET;
  
  
  bool config_doPdfSplit;
  int config_pdfStart;
  int config_pdfStep;
  
  float config_eeWeight;
  float config_emWeight;
  float config_mmWeight;
  float config_eeEff;
  float config_emEff;
  float config_mmEff;
  
  
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
  TTreeReaderValue<tree::Particle4Vector> ETmiss_vec;
  
  TTreeReaderValue<float> pt1;
  TTreeReaderValue<float> pt2;
  TTreeReaderValue<float> phi1;
  TTreeReaderValue<float> phi2;
  TTreeReaderValue<float> eta1;
  TTreeReaderValue<float> eta2;
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
  TTreeReaderValue<bool> trigMET;
  TTreeReaderValue<bool> trigDiEle;
  TTreeReaderValue<bool> trigDiMu;
  TTreeReaderValue<bool> trigMuEle;
  TTreeReaderValue<bool> trigDiEleMatch;
  TTreeReaderValue<bool> trigDiMuMatch;
  TTreeReaderValue<bool> trigMuEleMatch;
  
  TTreeReaderValue<float> calcHt;
  TTreeReaderValue<float> Mt2;
  
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
  TTreeReaderValue<tree::MET> met_gen;
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


  long nEntries;
  SignalPoint sp_;

  //float puWeights[3];
  float puWeights[6];
  float lepSfWeights[3];
  float photonSfWeights[3];
  float ewkWeights[3];
  float isrWeights[3];

  const char *emptyLabelPtr="";
  const char *weightLabelPtr="weight";
  const char *metLabelPtr=";#it{p}_{T}^{miss} (GeV)";



  map<selectionFolderName,map<selectionFolderName,map<Histograms1D,TH1F>>> h1Maps;
  map<selectionFolderName,map<selectionFolderName,map<Histograms2D,TH2F>>> h2Maps;
  
  map<selectionFolderName,map<selectionFolderName,map<Histograms1D,TEfficiency>>> eff1Maps;
  map<selectionFolderName,map<selectionFolderName,map<Histograms2D,TEfficiency>>> eff2Maps;
  

  map<cutFlowMapName,map<Histograms1D,TH1F>> c1Maps;

  map<SignalPoint,map<selectionFolderName,map<selectionFolderName,map<selectionFolderName,map<Histograms1D,TH1F>>>>> s1Maps;
  map<SignalPoint,map<selectionFolderName,map<selectionFolderName,map<selectionFolderName,map<Histograms2D,TH2F>>>>> s2Maps;

  map<selectionFolderName,map<selectionFolderName,map<selectionFolderName,map<Histograms1D,TH1F>>>> cr1Maps;


  map<cutFlowFlags, bool> decisionMapCutFlowFine;
  map<cutFlowFlags, float> decisionMapCutFlowFine_weight;
  void clearCutFlowMap();
  
  TH1F cutFlow;
  TH1F weightHisto;
  map<string,TH1F> weightHistoMap;
  float nGen;
  string inputName;


  int nWeights;

  
  float totalWeight;
  float totalWeightCalc;
  
  float zzKFactor;
  
  
  //const char *emptyLabelPtr="";
  //const char *weightLabelPtr="weight";
  //const char *metLabelPtr=";#it{p}_{T}^{miss} (GeV)";


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
  
  bool isTTInclusive;
  bool isZZ4L;
  bool isZZ2L;
  

  double startTime;
  ClassDef(myAnalyzer, 1)
};

