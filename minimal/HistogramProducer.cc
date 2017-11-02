#include "HistogramProducer.h"

HistogramProducer::HistogramProducer():
   photons(fReader, "photons"),
   jets(fReader, "jets"),
   electrons(fReader, "electrons"),
   muons(fReader, "muons"),
   genJets(fReader, "genJets"),
   genParticles(fReader, "genParticles"),
   intermediateGenParticles(fReader, "intermediateGenParticles"),
   met(fReader, "met"),
   metRaw(fReader, "met_raw"),
   met_JESu(fReader, "met_JESu"),
   met_JESd(fReader, "met_JESd"),
   met_JERu(fReader, "met_JERu"),
   met_JERd(fReader, "met_JERd"),
   nGoodVertices(fReader, "nGoodVertices"),
   nTracksPV(fReader, "nTracksPV"),
   pu_weight(fReader, "pu_weight"),
   mc_weight(fReader, "mc_weight"),
   //pdf_weights(fReader, "pdf_weights"),
   genHt(fReader, "genHt"),
   ht(fReader, "ht"),
   nTruePV(fReader, "true_nPV"),
   nISR(fReader, "nISR"),
   runNo(fReader, "runNo"),
   lumNo(fReader, "lumNo"),
   evtNo(fReader, "evtNo"),
   hlt_ht200(fReader, "HLT_PFHT200_v"),
   hlt_ht250(fReader, "HLT_PFHT250_v"),
   hlt_ht300(fReader, "HLT_PFHT300_v"),
   hlt_ht350(fReader, "HLT_PFHT350_v"),
   hlt_ht400(fReader, "HLT_PFHT400_v"),
   hlt_ht475(fReader, "HLT_PFHT475_v"),
   hlt_ht600(fReader, "HLT_PFHT600_v"),
   hlt_ht650(fReader, "HLT_PFHT650_v"),
   hlt_ht800(fReader, "HLT_PFHT800"),
   hlt_ele17_ele12_iso(fReader, "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"),
   hlt_ele23_ele12_iso(fReader, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"),
   hlt_mu17_mu8_iso(fReader, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v"),
   hlt_mu17_tkMu8_iso(fReader, "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v"),
   hlt_mu17_mu8_iso_dz(fReader, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"),
   hlt_mu17_tkMu8_iso_dz(fReader, "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"),
   hlt_tkMu17_tkMu8_iso_dz(fReader, "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"),
   hlt_mu17_ele12_iso(fReader, "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"),
   hlt_mu23_ele8_iso(fReader, "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v"),
   hlt_mu23_ele8_iso_dz(fReader, "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v"),
   hlt_mu23_ele12_iso(fReader, "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"),
   hlt_mu23_ele12_iso_dz(fReader, "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"),
   hlt_mu8_ele17_iso(fReader, "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v"),
   hlt_mu8_ele23_iso(fReader, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"),
   hlt_mu8_ele23_iso_dz(fReader, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"),
   hlt_mu12_ele23_iso(fReader, "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"),
   hlt_mu12_ele23_iso_dz(fReader, "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"),
   hlt_doubleEle33(fReader, "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v"),
   hlt_doubleEle33_mw(fReader, "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v"),
   hlt_mu27_tkMu8(fReader, "HLT_Mu27_TkMu8_v"),
   hlt_mu30_tkMu11(fReader, "HLT_Mu30_TkMu11_v"),
   hlt_mu30_ele30(fReader, "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v"),
   hlt_mu33_ele33(fReader, "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v"),
   startTime(time(NULL))//,
   //rand()
{
}

/*
void HistogramProducer::CalculateVariables(const tree::Lepton& l1, const tree::Lepton& l2, const tree::Photon& g,const tree::Particle mett ,const string particle)
{
   //additional variables
   //leptons
   bool isEle;
   if (particle=="e"){
      isEle=true;
   }else{
      isEle=false;
   }
   isDiElectron=isEle;
   float leptonMass = (isEle)?  0.51099895e-3 : 0.1056584;
   pt1 = l1.p.Pt();  
   pt2 = l2.p.Pt();  
   phi1 = l1.p.Phi();  
   phi2 = l2.p.Phi();  
   eta1 = l1.p.Eta();  
   eta2 = l2.p.Eta();
   int charge1 = (int) l1.charge;
   int charge2 = (int) l2.charge;
   chargeProduct = charge1*charge2;
   TLorentzVector lepton1 (0.,0.,0.,0.);
   TLorentzVector lepton2 (0.,0.,0.,0.);
   lepton1.SetPtEtaPhiM(pt1,eta1,phi1,leptonMass);
   lepton2.SetPtEtaPhiM(pt2,eta2,phi2,leptonMass);
   mll = (lepton1+lepton2).M();
   miniIso1 = l1.miniIso;
   miniIso2 = l2.miniIso;
   //photon
   //TLorentzVector photon (0.,0.,0.,0.);
   //photon.SetPtEtaPhiM(pt,eta,phi,0.);
   //pt = g.p.Pt();
   //eta = g.p.Eta();
   //phi = g.p.Phi();
   //Delta R between all three partciles
   deltaRll=lepton1.DeltaR(lepton2);
   //deltaRl1g=photon.DeltaR(lepton1);
   //deltaRl2g=photon.DeltaR(lepton2);
   //MET
   ETmiss = mett.p.Pt();
}
*/
void HistogramProducer::CalculateVariables(const tree::Lepton& l1, const tree::Lepton& l2, const tree::Particle mett ,const string particle)
{
   //additional variables
   //leptons
   bool isEle;
   if (particle=="e"){
      isEle=true;
   }else{
      isEle=false;
   }
   isDiElectron=isEle;
   //float leptonMass = (isEle)?  0.51099895e-3 : 0.1056584;
   float leptonMass =0.51099895e-3;
   pt1 = l1.p.Pt();  
   pt2 = l2.p.Pt();  
   phi1 = l1.p.Phi();  
   phi2 = l2.p.Phi();  
   eta1 = l1.p.Eta();  
   eta2 = l2.p.Eta();
   int charge1 = (int) l1.charge;
   int charge2 = (int) l2.charge;
   chargeProduct = charge1*charge2;
   TLorentzVector lepton1 (0.,0.,0.,0.);
   TLorentzVector lepton2 (0.,0.,0.,0.);
   lepton1.SetPtEtaPhiM(pt1,eta1,phi1,leptonMass);
   lepton2.SetPtEtaPhiM(pt2,eta2,phi2,leptonMass);
   mll = (lepton1+lepton2).M();
   miniIso1 = l1.miniIso;
   miniIso2 = l2.miniIso;
   //photon
   //TLorentzVector photon (0.,0.,0.,0.);
   //photon.SetPtEtaPhiM(pt,eta,phi,0.);
   //pt = -1.;
   //eta = -1.;
   //phi = -1.;
   //Delta R between all three partciles
   deltaRll=lepton1.DeltaR(lepton2);
   //deltaRl1g=-1.;
   //deltaRl2g=-1.;
   //MET
   ETmiss = mett.p.Pt();
}



bool HistogramProducer::CheckParticles(string selection="DY"){
   if (selection==""){
      return (((muons->size()>=2.) || (electrons->size()>=2.)) && (photons->size()>=1.));
   }
   if (selection=="DY"){
      return ((muons->size()>=2.) || (electrons->size()>=2.));      
   }
   return false;
}
bool HistogramProducer::Check2Ele(){
   return ((electrons->size()>=2.) && (photons->size()>=1.));
}
bool HistogramProducer::Check2Mu(){
   return ((muons->size()>=2.) && (photons->size()>=1.));
}


bool HistogramProducer::Cleaning(string selection="") //return if event is valid (2 Leptons exist)
{
   trigDiEle = *hlt_ele17_ele12_iso || *hlt_ele23_ele12_iso || *hlt_doubleEle33 || *hlt_doubleEle33_mw;
   trigDiMu  = *hlt_mu17_mu8_iso || *hlt_mu17_tkMu8_iso || *hlt_mu17_mu8_iso_dz || *hlt_mu17_tkMu8_iso_dz
               || *hlt_tkMu17_tkMu8_iso_dz || *hlt_mu27_tkMu8 || *hlt_mu30_tkMu11;
   trigMuEle = *hlt_mu17_ele12_iso || *hlt_mu23_ele8_iso || *hlt_mu23_ele8_iso_dz || *hlt_mu23_ele12_iso
               || *hlt_mu23_ele12_iso_dz || *hlt_mu8_ele17_iso || *hlt_mu8_ele23_iso
               || *hlt_mu8_ele23_iso_dz || *hlt_mu12_ele23_iso || *hlt_mu12_ele23_iso_dz
               || *hlt_mu30_ele30 || *hlt_mu33_ele33;
   trigHt    = *hlt_ht200 || *hlt_ht250 || *hlt_ht300 || *hlt_ht350 || *hlt_ht400
               || *hlt_ht475 || *hlt_ht600 || *hlt_ht650 || *hlt_ht800;
   
   if(selection==""){
      if (isData){           
         //find Selection (Dataset)
         bool isEESelection = inputName.find("DoubleEG") != string::npos;
         bool isMuMuSelection = inputName.find("DoubleMuon") != string::npos;
         bool isEMuSelection = inputName.find("MuonEG") != string::npos;
        
         //auto g  = photons->at(0);
         auto missingET = *met;
        
        
         if (isEESelection && trigDiEle && Check2Ele()) {
            auto e1 = electrons->at(0);
            auto e2 = electrons->at(1);
            if (trigDiMu && Check2Mu()) {
               auto m1 = muons->at(0);
               auto m2 = muons->at(1);
               if (min(e1.p.Pt(), e2.p.Pt()) > min(m1.p.Pt(), m2.p.Pt())){
                  //CalculateVariables(e1, e2, g, missingET,"e");
                  CalculateVariables(e1, e2, missingET,"e");
                  return true;
               }
            }else{
               //CalculateVariables(e1, e2, g, missingET,"e");
               CalculateVariables(e1, e2, missingET,"e");
               return true;
            }
         }
         if (isMuMuSelection && trigDiMu && Check2Mu()) {
            auto m1 = muons->at(0);
            auto m2 = muons->at(1);
            if (trigDiEle && Check2Ele()) {
               auto e1 = electrons->at(0);
               auto e2 = electrons->at(1);
               if (min(m1.p.Pt(), m2.p.Pt()) > min(e1.p.Pt(), e2.p.Pt())){
                  //CalculateVariables(m1, m2, g, missingET,"m");
                  CalculateVariables(m1, m2, missingET,"m");
                  return true;
               }
            }else{
               //CalculateVariables(m1, m2, g, missingET,"m");
               CalculateVariables(m1, m2, missingET,"m");
               return true;
            }
         }
      }else{
      //auto g  = photons->at(0);
      auto missingET = *met;
      
      if (trigDiEle && Check2Ele()) {
         auto e1 = electrons->at(0);
         auto e2 = electrons->at(1);
         if (trigDiMu && Check2Mu()) {
            auto m1 = muons->at(0);
            auto m2 = muons->at(1);
            if (min(e1.p.Pt(), e2.p.Pt()) > min(m1.p.Pt(), m2.p.Pt())){
               //CalculateVariables(e1, e2, g, missingET,"e");
               CalculateVariables(e1, e2, missingET,"e");
               return true;
            }
         }else{
            //CalculateVariables(e1, e2, g, missingET,"e");
            CalculateVariables(e1, e2, missingET,"e");
            return true;
         }
      }
      if (trigDiMu && Check2Mu()) {
         auto m1 = muons->at(0);
         auto m2 = muons->at(1);
         if (trigDiEle && Check2Ele()) {
            auto e1 = electrons->at(0);
            auto e2 = electrons->at(1);
            if (min(m1.p.Pt(), m2.p.Pt()) > min(e1.p.Pt(), e2.p.Pt())){
               //CalculateVariables(m1, m2, g, missingET,"m");
               CalculateVariables(m1, m2, missingET,"m");
               return true;
            }
         }
         else{
            //CalculateVariables(m1, m2, g, missingET,"m");
            CalculateVariables(m1, m2, missingET,"m");
            return true;
         }
      }
   }
   }
   
   if(selection=="DY"){
      if (isData){           
         //find Selection (Dataset)
         bool isEESelection = inputName.find("DoubleEG") != string::npos;
         bool isMuMuSelection = inputName.find("DoubleMuon") != string::npos;
         bool isEMuSelection = inputName.find("MuonEG") != string::npos;
        
         //auto g  = photons->at(0);
         auto missingET = *met;
        
        
         if (isEESelection && trigDiEle && Check2Ele()) {
            auto e1 = electrons->at(0);
            auto e2 = electrons->at(1);
            if (trigDiMu && Check2Mu()) {
               auto m1 = muons->at(0);
               auto m2 = muons->at(1);
               if (min(e1.p.Pt(), e2.p.Pt()) > min(m1.p.Pt(), m2.p.Pt())){
                  //CalculateVariables(e1, e2, g, missingET,"e");
                  CalculateVariables(e1, e2, missingET,"e");
                  return true;
               }
            }else{
               CalculateVariables(e1, e2, missingET,"e");
               return true;
            }
         }
         if (isMuMuSelection && trigDiMu && Check2Mu()) {
            auto m1 = muons->at(0);
            auto m2 = muons->at(1);
            if (trigDiEle && Check2Ele()) {
               auto e1 = electrons->at(0);
               auto e2 = electrons->at(1);
               if (min(m1.p.Pt(), m2.p.Pt()) > min(e1.p.Pt(), e2.p.Pt())){
                  CalculateVariables(m1, m2, missingET,"m");
                  return true;
               }
            }else{
               CalculateVariables(m1, m2, missingET,"m");
               return true;
            }
         }
      }else
      {
      //auto g  = photons->at(0);
      auto missingET = *met;
      
      if (trigDiEle && Check2Ele()) {
         auto e1 = electrons->at(0);
         auto e2 = electrons->at(1);
         if (trigDiMu && Check2Mu()) {
            auto m1 = muons->at(0);
            auto m2 = muons->at(1);
            if (min(e1.p.Pt(), e2.p.Pt()) > min(m1.p.Pt(), m2.p.Pt())){
               //CalculateVariables(e1, e2, g, missingET,"e");
               CalculateVariables(e1, e2, missingET,"e");
               return true;
            }
         }else{
            CalculateVariables(e1, e2, missingET,"e");
            return true;
         }
      }
      if (trigDiMu && Check2Mu()) {
         auto m1 = muons->at(0);
         auto m2 = muons->at(1);
         if (trigDiEle && Check2Ele()) {
            auto e1 = electrons->at(0);
            auto e2 = electrons->at(1);
            if (min(m1.p.Pt(), m2.p.Pt()) > min(e1.p.Pt(), e2.p.Pt())){
               CalculateVariables(m1, m2,  missingET,"m");
               return true;
            }
         }
         else{
            CalculateVariables(m1, m2, missingET,"m");
            return true;
         }
      }
   }
   }
   
   
   
   return false;
}


bool HistogramProducer::testSelection(const tree::Electron& pa){
   bool decision;
   decision = pa.isPassConvVeto && pa.passImpactParameter && (pa.p.Pt()>20.) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (pa.isMediumMVA);
   return decision;
}
bool HistogramProducer::testSelection(const tree::Muon& pa){
   bool decision;
   decision = pa.passImpactParameter && (pa.p.Pt()>20.) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (pa.isTight);
   return decision;
}
bool HistogramProducer::testSelection(const tree::Photon& pa){
   bool decision;
   decision = pa.passElectronVeto && (pa.p.Pt()>5.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isLoose);
   return decision;
}


//bool matchGenParticle(const tree::Particle& pa){
//
//for{genParticles::iterator = genParticles.begin(), iterator}

//}


void HistogramProducer::Init(TTree *tree)
{
   fReader.SetTree(tree);
   inputName = fReader.GetTree()->GetCurrentFile()->GetName();
   cutFlow = *((TH1F*)fReader.GetTree()->GetCurrentFile()->Get("TreeWriter/hCutFlow"));
   fReader.GetEntries(true); // jumps to last file
   string lastInputName = fReader.GetTree()->GetCurrentFile()->GetName();
  
   nEntries=fReader.GetTree()->GetEntries();
   cout<<"Processing "<<nEntries<<" Events!"<<endl;
   isData = inputName.find("Run201") != string::npos;
   isSignal = inputName.find("SMS") != string::npos;
  
   if (inputName!=lastInputName) {
      // adds cut flow of last file. This makes only sense if there is for 1 or two files
      cutFlow.Add((TH1F*)fReader.GetTree()->GetCurrentFile()->Get("TreeWriter/hCutFlow"));
   }
   float nGen = cutFlow.GetBinContent(2);

   InitHistograms();

}



void HistogramProducer::InitHistograms(){
   h1Map["met"] = TH1F("", ";#it{p}_{T}^{miss} (GeV)", 200, 0, 1000);
   h1Map["pt1"] = TH1F("", ";#it{p}_{T}^{leading} (GeV)", 200, 0, 1000);
   h1Map["pt2"] = TH1F("", ";#it{p}_{T}^{trailing} (GeV)", 200, 0, 1000);
   h1Map["pt1_ele"] = TH1F("", ";#it{p}_{T}^{leading e} (GeV)", 200, 0, 1000);
   h1Map["pt1_mu"] = TH1F("", ";#it{p}_{T}^{leading #mu} (GeV)", 200, 0, 1000);
   h1Map["pt2_ele"] = TH1F("", ";#it{p}_{T}^{trailing e} (GeV)", 200, 0, 1000);
   h1Map["pt2_mu"] = TH1F("", ";#it{p}_{T}^{trailing #mu} (GeV)", 200, 0, 1000);
   h1Map["pt"] = TH1F("", ";#it{p}_{T}^{#gamma} (GeV)", 200, 0, 1000);
   h1Map["m_ll"] = TH1F("", ";#it{m}_{ll} (GeV)", 200, 0, 1000);
   h1Map["m_ll_e"] = TH1F("", ";#it{m}_{ee} (GeV)", 200, 0, 1000);
   h1Map["m_ll_m"] = TH1F("", ";#it{m}_{#mu#mu} (GeV)", 200, 0, 1000);
   
   //h1Map["m_ll_m_gen"]=TH1F("", ";#it{m}_{#mu#mu} (GeV)", 200, 0, 1000);
}


void HistogramProducer::SlaveBegin(TTree *tree)
{
}

bool HistogramProducer::SelectEvent(string selection=""){
    
   if (CheckParticles()){//only Fill Histograms if nMu>=2, nEle>=2, nGamma>=1
      if (Cleaning()){
         selEvent L;
         L.isDiElectron=isDiElectron;
         L.trigDiEle=trigDiEle;
         L.trigDiMu=trigDiMu;
         L.trigMuEle=trigMuEle;
         L.trigHt=trigHt;

         L.pt1=pt1;
         L.pt2=pt2;
         L.phi1=phi1;
         L.phi2=phi2;
         L.eta1=eta1;
         L.eta2=eta2;
         L.chargeProduct=chargeProduct;
         L.deltaRll=deltaRll;
         //L.deltaRl1e=deltaRl1e;
         //L.deltaRl2e=deltaRl2e;
         L.mll=mll;
         L.miniIso1=miniIso1;
         L.miniIso2=miniIso2;
         //L.pt=pt;
         //L.eta=eta;
         //L.phi=phi;
         //L.miniIso=miniIso;
         L.ETmiss=ETmiss;
         L.totalWeight=totalWeight;
         
         if (selection==""){
            if (isDiElectron){
               auto e1 = electrons->at(0);
               auto e2 = electrons->at(1);
               auto g = photons->at(0);
               if(testSelection(e1) && testSelection(e2)){
                  for (vector<tree::Photon>::iterator it = photons->begin(); it != photons->end(); ++it){
                     if (testSelection(*it)){
                        L.selPhotons.push_back(*it);
                     }
                  }
                  selectedEvent=L;
                  return true;
               }
            }
            else{
               auto m1 = muons->at(0); 
               auto m2 = muons->at(1); 
               auto g = photons->at(0);
               if(testSelection(m1) && testSelection(m2)){
                  selectedEvent=L;
                  return true;
               }
            }
         }
         if (selection=="DY"){
            if (isDiElectron){
               auto e1 = electrons->at(0);
               auto e2 = electrons->at(1);
               //auto g = photons->at(0);
               if(testSelection(e1) && testSelection(e2)){
                  selectedEvent=L;
                  return true;
               }
            }
            else{
               auto m1 = muons->at(0); 
               auto m2 = muons->at(1); 
               //auto g = photons->at(0);
               if(testSelection(m1) && testSelection(m2)){
                  selectedEvent=L;
                  return true;
               }
            }
         }
      }
   }
  
  
  return false;
}


void HistogramProducer::FillHistograms(){
   if(SelectEvent("")){
      //cout<<"Accepted"<<endl;
      h1Map["met"].Fill(selectedEvent.ETmiss, selectedEvent.totalWeight);
      //h1Map["pt"].Fill(selectedEvent.pt, selectedEvent.totalWeight);
      h1Map["pt1"].Fill(selectedEvent.pt1, selectedEvent.totalWeight);
      h1Map["pt2"].Fill(selectedEvent.pt2, selectedEvent.totalWeight);
      h1Map["m_ll"].Fill(selectedEvent.mll,selectedEvent.totalWeight);   
      (selectedEvent.isDiElectron)? h1Map["m_ll_e"].Fill(selectedEvent.mll,selectedEvent.totalWeight) : h1Map["m_ll_m"].Fill(selectedEvent.mll,selectedEvent.totalWeight);
      (selectedEvent.isDiElectron)? h1Map["pt1_ele"].Fill(selectedEvent.pt1,selectedEvent.totalWeight) : h1Map["pt1_mu"].Fill(selectedEvent.pt1,selectedEvent.totalWeight);
      (selectedEvent.isDiElectron)? h1Map["pt2_ele"].Fill(selectedEvent.pt2,selectedEvent.totalWeight) : h1Map["pt2_mu"].Fill(selectedEvent.pt2,selectedEvent.totalWeight);
   }
   
   
   //auto 
   
}


Bool_t HistogramProducer::Process(Long64_t entry){
   if (entry>1000000){
      return kTRUE;
   }
   fReader.SetLocalEntry(entry);
   //auto totalWeight = *mc_weight * *pu_weight;
   totalWeight = *mc_weight * *pu_weight;



   //print out progression
   double currentEventNo = (double) entry;
   double totalEventNo = (double) nEntries;
   double tempVerbose = currentEventNo/totalEventNo*100.;
   if(entry%500000==0){
      int tempVerbose2 = (int)tempVerbose;
      cout<<entry<<" / "<<nEntries<<" | "<<tempVerbose2<<"%"<<endl;
   }

   FillHistograms();
   
   return kTRUE;
}

template<typename T>
void save2File(const map<string,T>& hMap, TFile& file)
{
   for (auto& h : hMap) {
      h.second.Write(h.first.c_str(), TObject::kWriteDelete);
   }
}

void HistogramProducer::Terminate()
{
   auto outputName = getOutputFilename(inputName);
   TFile file(outputName.c_str(), "RECREATE");
   save2File(h1Map, file);
   cutFlow.Write("hCutFlow");
   file.Close();
   cout << "Created " << outputName << " in " << (time(NULL) - startTime)/60 << " min" << endl;
}
//
//void HistogramProducer::resetSelection() {
   //selPhotons.clear();
   //selJets.clear();
   //selBJets.clear();
   //selHEJets.clear();
   //selElectrons.clear();
   //selMuons.clear();
   //artificialPhotons.clear();
//}
