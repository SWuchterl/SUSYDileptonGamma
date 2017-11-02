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
#include "CutFlow.h"
#include "Resolution.h"


class SignalScan : public TSelector {
 public:

  SignalScan();
  virtual ~SignalScan() { }

  virtual void Init(TTree *tree);
  virtual void SlaveBegin(TTree *tree){};
  virtual Bool_t Process(Long64_t entry);
  virtual void Terminate();
  virtual Int_t Version() const { return 2; }

  float getPhotonWeight(const tree::Photon& p);
  void fillSignalSelection(const string p, string const& s, float weight, bool isSel, bool, float, float);

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
  TTreeReaderValue<tree::MET> metGen;
  TTreeReaderValue<Float_t> pu_weight;
  TTreeReaderValue<Char_t> mc_weight;
  TTreeReaderValue<std::vector<Float_t>> pdf_weights;

  TTreeReaderValue<Int_t> nGoodVertices;
  TTreeReaderValue<Float_t> rho;
  TTreeReaderValue<ULong64_t> eventNo;
  TTreeReaderValue<UInt_t> runNo;
  TTreeReaderValue<UInt_t> lumNo;
  TTreeReaderValue<Int_t> nTruePV;
  TTreeReaderValue<Int_t> nISR;

  // signal scan
  TTreeReaderValue<UShort_t> signal_nBinos;
  TTreeReaderValue<UShort_t> signal_m1;
  TTreeReaderValue<UShort_t> signal_m2;

  map<string,map<string,map<string,TH1F>>> h1MapsMaps;
  map<string,unsigned int> nEventMap;

  string inputName;

  vector<tree::Photon*> selPhotons;
  vector<tree::Jet*> selJets;
  vector<tree::Electron*> selElectrons;
  vector<tree::Muon*> selMuons;

  map<string,Weighter> weighters;

  bool puSummer16;
  bool electroweak;
  double startTime;
  ClassDef(SignalScan, 1)
};

SignalScan::SignalScan():
  photons(fReader, "photons"),
  jets(fReader, "jets"),
//  electrons(fReader, "electrons"),
//  muons(fReader, "muons"),
  genJets(fReader, "genJets"),
  genParticles(fReader, "genParticles"),
//  intermediateGenParticles(fReader, "intermediateGenParticles"),
 met(fReader, "met"),
//  metRaw(fReader, "met_raw"),
  met_JESu(fReader, "met_JESu"),
  met_JESd(fReader, "met_JESd"),
  met_JERu(fReader, "met_JERu"),
  met_JERd(fReader, "met_JERd"),
  metGen(fReader, "met_gen"),
//  nGoodVertices(fReader, "nGoodVertices"),
  pu_weight(fReader, "pu_weight"),
  mc_weight(fReader, "mc_weight"),
  pdf_weights(fReader, "pdf_weights"),
//  rho(fReader, "rho"),
//  runNo(fReader, "runNo"),
  nTruePV(fReader, "true_nPV"),
  nISR(fReader, "nISR"),
  signal_nBinos(fReader, "signal_nBinos"),
  signal_m1(fReader, "signal_m1"),
  signal_m2(fReader, "signal_m2"),
  puSummer16(false),
  electroweak(false),
  startTime(time(NULL))
{
  weighters["sf_photon_id_loose"] = Weighter("../plotter/data/dataMcScaleFactors_80X.root", "EGamma_SF2D");
  weighters["sf_photon_pixel"] = Weighter("../plotter/data/ScalingFactors_80X_Summer16.root", "Scaling_Factors_HasPix_R9 Inclusive");
  weighters.at("sf_photon_id_loose").fillOverflow2d();
  weighters.at("sf_photon_pixel").fillOverflow2d();
}

void SignalScan::Init(TTree *tree)
{
  fReader.SetTree(tree);
  inputName = fReader.GetTree()->GetCurrentFile()->GetName();
  electroweak = inputName.find("TChi") != string::npos;
  puSummer16 = inputName.find("T5bbbbZg") != string::npos ||
               inputName.find("T5ttttZg") != string::npos ||
               inputName.find("T5qqqqHg") != string::npos ||
               inputName.find("GMSB") != string::npos;

  string puUp = "pileupWeightUp_mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU";
  string puDn = "pileupWeightDown_mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU";
  if (puSummer16) {
    puUp = "pileupWeightUp_mix_2016_25ns_Moriond17MC_PoissonOOTPU";
    puDn = "pileupWeightDown_mix_2016_25ns_Moriond17MC_PoissonOOTPU";
    weighters["puWeight"] = Weighter("../../CMSSW/treewriter/CMSSW_8_0_25/src/TreeWriter/PUreweighting/data/puWeights.root", "pileupWeight_mix_2016_25ns_Moriond17MC_PoissonOOTPU");
  }
  weighters["puWeightUp"] = Weighter("../../CMSSW/treewriter/CMSSW_8_0_25/src/TreeWriter/PUreweighting/data/puWeights.root", puUp);
  weighters["puWeightDn"] = Weighter("../../CMSSW/treewriter/CMSSW_8_0_25/src/TreeWriter/PUreweighting/data/puWeights.root", puDn);

}

map<string,TH1F> initSignalHistograms() {
  map<string,TH1F> hMap;
  hMap["met"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_isr"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_isrUp"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_isrDn"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["metGen"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_nopu"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_nopuUp"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_nopuDn"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_puUp"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_puDn"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_jesUp"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_jesDn"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_jerUp"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_jerDn"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_weight1"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_weight2"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_weight3"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_weight4"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_weight5"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_weight6"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_weight7"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_weight8"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_contamination"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  return hMap;
}

void SignalScan::fillSignalSelection(const string p, string const& s, float weight, bool isSel=false, bool controlRegion=false, float shift=1., float scale=0.)
{
  if (!h1MapsMaps.count(p)) {
    h1MapsMaps[p] = {{"signal_lowEMHT", initSignalHistograms()}, {"signal_highEMHT", initSignalHistograms()},
                     {"signal_allAngles_lowEMHT", initSignalHistograms()}, {"signal_allAngles_highEMHT", initSignalHistograms()}};
  }
  auto m1 = &h1MapsMaps[p][s];
  auto _met = met->p.Pt();
  auto isrPt = 0;
  if (electroweak) {
    isrPt = getPropagatorPt(*genParticles);
  }
  auto isrWeight = electroweak ? isrReweightingEWK(isrPt) : isrReweighting(*nISR);
  auto isrWeightE = electroweak ? isrReweightingEWK(isrPt, true) : isrReweighting(*nISR, true);

  m1->at("met").Fill(isSel?_met:-1, weight);
  m1->at("met_isr").Fill(isSel?_met:-1, weight*isrWeight);
  m1->at("met_isrUp").Fill(isSel?_met:-1, weight*(isrWeight+isrWeightE));
  m1->at("met_isrDn").Fill(isSel?_met:-1, weight*(isrWeight-isrWeightE));
  m1->at("metGen").Fill(isSel?metGen->p.Pt():-1, weight);
  m1->at("met_nopu").Fill(isSel?_met:-1, weight/ *pu_weight);
  if (*nTruePV>=20) m1->at("met_nopuUp").Fill(isSel?_met:-1, weight/ *pu_weight);
  else             m1->at("met_nopuDn").Fill(isSel?_met:-1, weight/ *pu_weight);
  m1->at("met_puUp").Fill(isSel?_met:-1, weight*weighters.at("puWeightUp").getWeight(*nTruePV)/ *pu_weight);
  m1->at("met_puDn").Fill(isSel?_met:-1, weight*weighters.at("puWeightDn").getWeight(*nTruePV)/ *pu_weight);
  m1->at("met_jesUp").Fill(isSel?met_JESu->p.Pt():-1, weight);
  m1->at("met_jesDn").Fill(isSel?met_JESd->p.Pt():-1, weight);
  m1->at("met_jerUp").Fill(isSel?met_JERu->p.Pt():-1, weight);
  m1->at("met_jerDn").Fill(isSel?met_JERd->p.Pt():-1, weight);
  if (pdf_weights->size()>8) {
  m1->at("met_weight1").Fill(isSel?_met:-1, weight*pdf_weights->at(1));
  m1->at("met_weight2").Fill(isSel?_met:-1, weight*pdf_weights->at(2));
  m1->at("met_weight3").Fill(isSel?_met:-1, weight*pdf_weights->at(3));
  m1->at("met_weight4").Fill(isSel?_met:-1, weight*pdf_weights->at(4));
  m1->at("met_weight5").Fill(isSel?_met:-1, weight*pdf_weights->at(5));
  m1->at("met_weight6").Fill(isSel?_met:-1, weight*pdf_weights->at(6));
  m1->at("met_weight7").Fill(isSel?_met:-1, weight*pdf_weights->at(7));
  m1->at("met_weight8").Fill(isSel?_met:-1, weight*pdf_weights->at(8));
  }

  m1->at("met_contamination").Fill(controlRegion?_met*shift:-1, weight*scale);
}


float SignalScan::getPhotonWeight(const tree::Photon& p) {
  float weight = 1.;
  auto pt = p.p.Pt();
  auto eta = p.p.Eta();
  // sf for id and electron veto
  weight = weighters.at("sf_photon_id_loose").getWeight(eta, pt) * weighters.at("sf_photon_pixel").getWeight(fabs(eta), pt);
  // trigger efficiency
  weight *= fabs(eta)<photonsEtaMaxBarrel ? 0.964 : 0.94;

  return weight;
}

Bool_t SignalScan::Process(Long64_t entry)
{
  if (!(entry%int(2e6))) cout << 1.*entry / fReader.GetEntries(false) << endl;
  fReader.SetLocalEntry(entry);
  auto pn = getSignalPointName(1, *signal_m1, *signal_m2);
  auto pnGG = getSignalPointName(2, *signal_m1, *signal_m2);
  //auto pn = getSignalPointName(*signal_nBinos, *signal_m1, *signal_m2);

  if (!nEventMap.count(pn)) {
    unsigned int nGen = 0;
    string cutFlowName = "TreeWriter/hCutFlow";
    inputName = fReader.GetTree()->GetCurrentFile()->GetName();
    if (inputName.find("1600_100") != string::npos) { // testing one point for synchronization
      TH1F* cutFlow = (TH1F*)fReader.GetTree()->GetCurrentFile()->Get(cutFlowName.c_str());
      if (cutFlow) {
        nGen = cutFlow->GetBinContent(2);
      } else {
        cout << "Could not read cut flow histogram" << endl;
      }
    } else {
      if (inputName.find("T5ttttZg") != string::npos) {
        cutFlowName += "T5ttttZg";
      } else if (inputName.find("T5bbbbZg") != string::npos) {
        cutFlowName += "T5bbbbZg";
      } else if (inputName.find("T5qqqqHg") != string::npos) {
        cutFlowName += "T5qqqqHg";
      } else if (inputName.find("T5") != string::npos) {
        cutFlowName += "T5Wg";
      } else if (inputName.find("T6") != string::npos) {
        cutFlowName += "T6Wg";
      } else if (inputName.find("TChiWG") != string::npos) {
        cutFlowName += "TChiWG";
      } else if (inputName.find("TChiNG") != string::npos) {
        cutFlowName += "TChiNG";
      } else if (inputName.find("GMSB") != string::npos) {
        cutFlowName += "GMSB";
      } else {
        cout << "Could not find correct cutflowname for " << inputName << endl;
      }
      cutFlowName += "_"+to_string(*signal_m1);
      if (*signal_m2) cutFlowName += "_"+to_string(*signal_m2);
      TH1F* cutFlow = (TH1F*)fReader.GetTree()->GetCurrentFile()->Get(cutFlowName.c_str());
      if (cutFlow) {
        nGen = cutFlow->GetBinContent(2);
      } else {
        cout << "Could not read cutFlow histogram " << cutFlowName << endl;
      }
    }
    nEventMap[pn] = nGen;
    nEventMap[pnGG] = nGen/4;
  }

  if (unMatchedSuspiciousJet(*jets, *genJets)) {
    nEventMap[pn] -= 1;
    return kTRUE;
  }

  selPhotons.clear();
  selJets.clear();
  for (auto& photon : *photons) {
    if (photon.isLoose && !photon.hasPixelSeed && photon.p.Pt() > 100 && fabs(photon.p.Eta()) < photonsEtaMaxBarrel) {
      selPhotons.push_back(&photon);
    }
  }

 bool orthogonal = false;
  if (selPhotons.size()) {
    auto dPhi = fabs(met->p.DeltaPhi(selPhotons.at(0)->p));
    orthogonal = .3<dPhi && dPhi<2.84;
  }

  for (auto& jet : *jets) {
    if (!jet.isLoose
//      || jet.hasPhotonMatch || jet.hasElectronMatch || jet.hasMuonMatch
      || indexOfMatchedParticle<tree::Photon*>(jet, selPhotons, .4) >= 0
      || jet.p.Pt() < 30 || fabs(jet.p.Eta()) > 3) continue;
    selJets.push_back(&jet);
  }
  float emht=0;
  for (auto& p : selPhotons) emht += p->p.Pt();
  for (auto& p : selJets) emht += p->p.Pt();

  if (puSummer16) {
    *pu_weight = weighters.at("puWeight").getWeight(*nTruePV);
  }
  auto selW = *mc_weight * *pu_weight;
  if (selPhotons.size()) selW *= getPhotonWeight(*selPhotons.at(0));
  bool passSelection = orthogonal && selPhotons.size() && emht>700;
  fillSignalSelection(pn, "signal_lowEMHT", selW, passSelection && emht<2000, !selPhotons.size() && emht>700 && emht<2000, 0.905, 0.001271);
  fillSignalSelection(pn, "signal_highEMHT", selW, passSelection && emht>=2000, !selPhotons.size() && emht>=2000, 0.85, 0.003678);
  if (*signal_nBinos == 2) {
    fillSignalSelection(pnGG, "signal_lowEMHT", selW, passSelection && emht<2000, !selPhotons.size() && emht>700 && emht<2000, 0.905, 0.001271);
    fillSignalSelection(pnGG, "signal_highEMHT", selW, passSelection && emht>=2000, !selPhotons.size() && emht>=2000, 0.85, 0.003678);
  }

  // do not make orthogonal cut
  bool passSelectionAllAngles = selPhotons.size() && emht>700;
  fillSignalSelection(pn, "signal_allAngles_lowEMHT", selW, passSelectionAllAngles && emht<2000, !selPhotons.size() && emht>700 && emht<2000, 0.905, 0.001271);
  fillSignalSelection(pn, "signal_allAngles_highEMHT", selW, passSelectionAllAngles && emht>=2000, !selPhotons.size() && emht>=2000, 0.85, 0.003678);
  if (*signal_nBinos == 2) {
    fillSignalSelection(pnGG, "signal_allAngles_lowEMHT", selW, passSelectionAllAngles && emht<2000, !selPhotons.size() && emht>700 && emht<2000, 0.905, 0.001271);
    fillSignalSelection(pnGG, "signal_allAngles_highEMHT", selW, passSelectionAllAngles && emht>=2000, !selPhotons.size() && emht>=2000, 0.85, 0.003678);
  }

  return kTRUE;
}

void cdNewDir(TFile& file, const string& name) {
  if (!file.Get(name.c_str())) {
    file.mkdir(name.c_str());
  }
  file.cd(name.c_str());
}

template<typename T>
void save2File(map<string,map<string,map<string,T>>>& hMapsMaps, TFile& file, map<string,unsigned>& nEventMap )
{
  for (auto& hMapMapIt: hMapsMaps) {
    file.mkdir(hMapMapIt.first.c_str());
    for (auto& hMapIt: hMapMapIt.second) {
      file.mkdir((hMapMapIt.first+"/"+hMapIt.first).c_str());
      file.cd((hMapMapIt.first+"/"+hMapIt.first).c_str());
      for (auto& h: hMapIt.second) {
        h.second.Scale(1./nEventMap[hMapMapIt.first]);
        h.second.Write(h.first.c_str(), TObject::kWriteDelete);
      }
    }
  }
}

void SignalScan::Terminate()
{
  auto outputName = getOutputFilename(inputName, "signalScan");
  TFile file(outputName.c_str(), "RECREATE");
  save2File(h1MapsMaps, file, nEventMap);
  file.Close();
  cout << "Created " << outputName << " in " << (time(NULL) - startTime)/60 << " min" << endl;
}

