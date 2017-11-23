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
   hlt_ht800(fReader, "HLT_PFHT800_v"),
   hlt_ele17_ele12_iso(fReader, "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"),
   hlt_ele23_ele12_iso(fReader, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"),
   hlt_mu17_mu8_iso(fReader, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v"),
   hlt_mu17_tkMu8_iso(fReader, "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v"),
   hlt_mu17_mu8_iso_dz(fReader, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"),
   hlt_mu17_tkMu8_iso_dz(fReader, "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"),
   hlt_tkMu17_tkMu8_iso_dz(fReader, "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"),
   //hlt_mu17_ele12_iso(fReader, "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"),
   //hlt_mu23_ele8_iso(fReader, "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v"),
   //hlt_mu23_ele8_iso_dz(fReader, "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v"),
   //hlt_mu23_ele12_iso(fReader, "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"),
   //hlt_mu23_ele12_iso_dz(fReader, "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"),
   //hlt_mu8_ele17_iso(fReader, "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v"),
   //hlt_mu8_ele23_iso(fReader, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"),
   //hlt_mu8_ele23_iso_dz(fReader, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"),
   //hlt_mu12_ele23_iso(fReader, "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"),
   //hlt_mu12_ele23_iso_dz(fReader, "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"),
   hlt_doubleEle33(fReader, "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v"),
   hlt_doubleEle33_mw(fReader, "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v"),
   hlt_mu27_tkMu8(fReader, "HLT_Mu27_TkMu8_v"),
   hlt_mu30_tkMu11(fReader, "HLT_Mu30_TkMu11_v"),
   //hlt_mu30_ele30(fReader, "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v"),
   //hlt_mu33_ele33(fReader, "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v"),
   startTime(time(NULL))//,
   //rand()
{
}



bool HistogramProducer::GenPhotonVeto(){
   bool isGenPhoton=false;
   if(!isData){
      for (vector<tree::GenParticle>::iterator it = genParticles->begin(); it != genParticles->end(); ++it){
         if(it->pdgId==22){
            //isGenPhoton=isGenPhoton || it->fromHardProcess;
            isGenPhoton=isGenPhoton || true;
         }
      }
   }
   //bool cutPrompt = noPromptPhotons && count_if(genParticles->begin(), genParticles->end(), [] (const tree::GenParticle& p) { return p.pdgId==22 && p.promptStatus == DIRECTPROMPT;});
   //return cutPrompt;
   //return false;
   return isGenPhoton && noPromptPhotons;
}


void HistogramProducer::CalculateVariables(const tree::Lepton& l1, const tree::Lepton& l2, const tree::Particle mett ,const particleType particle)
{
   //leptons
   bool isEle;
   if (particle==E){
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
   lep1=lepton1;
   lep2=lepton2;
   mll = (lepton1+lepton2).M();
   miniIso1 = l1.miniIso;
   miniIso2 = l2.miniIso;
   //Delta R between all three partciles
   deltaRll=lepton1.DeltaR(lepton2);
   //MET
   ETmiss = mett.p.Pt();
}



bool HistogramProducer::CheckParticles(){
   return ((muons->size()>=2.) || (electrons->size()>=2.));      
}
bool HistogramProducer::Check2Ele(){
   return (electrons->size()>=2.);
}
bool HistogramProducer::Check2Mu(){
   return (muons->size()>=2.);
}


bool HistogramProducer::Cleaning() //return if event is valid (2 Leptons exist)
{
   //cout<<"signal?"<<isSignal<<endl;
   if(!isSignal){
      trigDiEle = *hlt_ele17_ele12_iso || *hlt_ele23_ele12_iso || *hlt_doubleEle33 || *hlt_doubleEle33_mw;
      trigDiMu  = *hlt_mu17_mu8_iso || *hlt_mu17_tkMu8_iso || *hlt_mu17_mu8_iso_dz || *hlt_mu17_tkMu8_iso_dz
                  || *hlt_tkMu17_tkMu8_iso_dz || *hlt_mu27_tkMu8 || *hlt_mu30_tkMu11;
      //trigMuEle = *hlt_mu17_ele12_iso || *hlt_mu23_ele8_iso || *hlt_mu23_ele8_iso_dz || *hlt_mu23_ele12_iso
                  //|| *hlt_mu23_ele12_iso_dz || *hlt_mu8_ele17_iso || *hlt_mu8_ele23_iso
                  //|| *hlt_mu8_ele23_iso_dz || *hlt_mu12_ele23_iso || *hlt_mu12_ele23_iso_dz
                  //|| *hlt_mu30_ele30 || *hlt_mu33_ele33;
      trigMuEle = false;
      trigHt    = *hlt_ht200 || *hlt_ht250 || *hlt_ht300 || *hlt_ht350 || *hlt_ht400
                  || *hlt_ht475 || *hlt_ht600 || *hlt_ht650 || *hlt_ht800;
      //trigHt    = true;
   }else{
      trigDiEle = true;
      trigDiMu = true;
      trigMuEle = false;
      trigHt = true;
   }
   

   if (isData){           
      //find Selection (Dataset)
      bool isEESelection = inputName.find("DoubleEG") != string::npos;
      bool isMuMuSelection = inputName.find("DoubleMuon") != string::npos;
      bool isEMuSelection = inputName.find("MuonEG") != string::npos;
      bool isHTSelection = inputName.find("JetHT")!= string::npos;
      auto missingET = *met;
        
      if (isEESelection && trigDiEle && Check2Ele()) {
         auto e1 = electrons->at(0);
         auto e2 = electrons->at(1);
         if (trigDiMu && Check2Mu()) {
            auto m1 = muons->at(0);
            auto m2 = muons->at(1);
            if (min(e1.p.Pt(), e2.p.Pt()) > min(m1.p.Pt(), m2.p.Pt())){
               CalculateVariables(e1, e2, missingET,E);
               return true;
            }else{return false;}
         }else{
            CalculateVariables(e1, e2, missingET,E);
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
               CalculateVariables(m1, m2, missingET,M);
               return true;
            }else{return false;}
         }else{
            CalculateVariables(m1, m2, missingET,M);
            return true;
         }
      }
      return false;
      }else{
         auto missingET = *met;
      
         if (trigDiEle && Check2Ele()) {
            auto e1 = electrons->at(0);
            auto e2 = electrons->at(1);
            if (trigDiMu && Check2Mu()) {
               auto m1 = muons->at(0);
               auto m2 = muons->at(1);
               if (min(e1.p.Pt(), e2.p.Pt()) > min(m1.p.Pt(), m2.p.Pt())){
                  CalculateVariables(e1, e2, missingET,E);
                  return true;
               }else{return false;}
            }else{
            CalculateVariables(e1, e2, missingET,E);
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
                  CalculateVariables(m1, m2, missingET,M);
                  return true;
               }else{return false;}
            }else{
            CalculateVariables(m1, m2, missingET,M);
            return true;
            }
         }
      return false;
   }
   return false;
}

bool HistogramProducer::CleaningTriggerStudies() //return if event is valid (2 Leptons exist)
{
   //cout<<"signal?"<<isSignal<<endl;
   if(!isSignal){
      trigDiEle = *hlt_ele17_ele12_iso || *hlt_ele23_ele12_iso || *hlt_doubleEle33 || *hlt_doubleEle33_mw;
      trigDiMu  = *hlt_mu17_mu8_iso || *hlt_mu17_tkMu8_iso || *hlt_mu17_mu8_iso_dz || *hlt_mu17_tkMu8_iso_dz
                  || *hlt_tkMu17_tkMu8_iso_dz || *hlt_mu27_tkMu8 || *hlt_mu30_tkMu11;
      //trigMuEle = *hlt_mu17_ele12_iso || *hlt_mu23_ele8_iso || *hlt_mu23_ele8_iso_dz || *hlt_mu23_ele12_iso
                  //|| *hlt_mu23_ele12_iso_dz || *hlt_mu8_ele17_iso || *hlt_mu8_ele23_iso
                  //|| *hlt_mu8_ele23_iso_dz || *hlt_mu12_ele23_iso || *hlt_mu12_ele23_iso_dz
                  //|| *hlt_mu30_ele30 || *hlt_mu33_ele33;
      trigMuEle = false;
      trigHt    = *hlt_ht200 || *hlt_ht250 || *hlt_ht300 || *hlt_ht350 || *hlt_ht400
                  || *hlt_ht475 || *hlt_ht600 || *hlt_ht650 || *hlt_ht800;
      //trigHt    = true;
   }else{
      trigDiEle = true;
      trigDiMu = true;
      trigMuEle = false;
      trigHt = true;
   }
   

   if (isData){           
      //find Selection (Dataset)
      bool isEESelection = inputName.find("DoubleEG") != string::npos;
      bool isMuMuSelection = inputName.find("DoubleMuon") != string::npos;
      bool isEMuSelection = inputName.find("MuonEG") != string::npos;
      bool isHTSelection = inputName.find("JetHT")!= string::npos;
      auto missingET = *met;
        
      if (isHTSelection && Check2Ele()) {
         auto e1 = electrons->at(0);
         auto e2 = electrons->at(1);
         if (Check2Mu()) {
            auto m1 = muons->at(0);
            auto m2 = muons->at(1);
            if (min(e1.p.Pt(), e2.p.Pt()) > min(m1.p.Pt(), m2.p.Pt())){
               CalculateVariables(e1, e2, missingET,E);
               return true;
            }else{return false;}
         }else{
            CalculateVariables(e1, e2, missingET,E);
            return true;
         }
      }
      if (isMuMuSelection && Check2Mu()) {
         auto m1 = muons->at(0);
         auto m2 = muons->at(1);
         if (Check2Ele()) {
            auto e1 = electrons->at(0);
            auto e2 = electrons->at(1);
            if (min(m1.p.Pt(), m2.p.Pt()) > min(e1.p.Pt(), e2.p.Pt())){
               CalculateVariables(m1, m2, missingET,M);
               return true;
            }else{return false;}
         }else{
            CalculateVariables(m1, m2, missingET,M);
            return true;
         }
      }
      return false;
   }else{
      auto missingET = *met;
   
      if (trigHt && Check2Ele()) {
         auto e1 = electrons->at(0);
         auto e2 = electrons->at(1);
         if (Check2Mu()) {
            auto m1 = muons->at(0);
            auto m2 = muons->at(1);
            if (min(e1.p.Pt(), e2.p.Pt()) > min(m1.p.Pt(), m2.p.Pt())){
               CalculateVariables(e1, e2, missingET,E);
               return true;
            }else{return false;}
         }else{
         CalculateVariables(e1, e2, missingET,E);
         return true;
         }
      }
      if (trigHt && Check2Mu()) {
         auto m1 = muons->at(0);
         auto m2 = muons->at(1);
         if (Check2Ele()) {
            auto e1 = electrons->at(0);
            auto e2 = electrons->at(1);
            if (min(m1.p.Pt(), m2.p.Pt()) > min(e1.p.Pt(), e2.p.Pt())){
               CalculateVariables(m1, m2, missingET,M);
               return true;
            }else{return false;}
         }else{
         CalculateVariables(m1, m2, missingET,M);
         return true;
         }
      }
   return false;
   }
   return false;
}




bool HistogramProducer::testSelection(const tree::Electron& pa,selectionType selection ,bool leading){
   bool decision=false;
   if ((selection==UNCUT)||(selection==PHOTON)){
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && ElectronTightMVA(pa.p.Eta(),pa.p.Pt(),pa.mvaValue);
      decision = pa.isPassConvVeto && pa.passImpactParameter && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && pa.isTightMVA;
   }
   if(selection==SEL || selection==ONZ){
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && ElectronTightMVA(pa.p.Eta(),pa.p.Pt(),pa.mvaValue);
      decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVA;
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && ElectronTightMVA(pa.p.Eta(),pa.p.Pt(),pa.mvaValue);
   }
   if(selection==DILEP){
      decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVA ;
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && ElectronTightMVA(pa.p.Eta(),pa.p.Pt(),pa.mvaValue) ;
   }
   if(selection==TRIGSEL || selection==TRIGDILEP || selection==TRIGONZ){
      decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>0.) : (pa.p.Pt()>0.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVA;
   }
   return decision;
}
bool HistogramProducer::testSelection(const tree::Muon& pa, selectionType selection, bool leading){
   bool decision=false;
   if ((selection==UNCUT)||(selection==PHOTON)){
      decision = pa.passImpactParameter  && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (pa.isTight);
   }
   if (selection==SEL || selection == ONZ){
      decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (deltaRll>0.1);
   }
   if(selection==DILEP){
      decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (deltaRll>0.1);
   }
   if(selection==TRIGSEL || selection==TRIGDILEP || selection==TRIGONZ){
      decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>0.) : (pa.p.Pt()>0.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (deltaRll>0.1);
   }
   return decision;
}
bool HistogramProducer::testSelection(const selPhoton& pa, selectionType selection){
   bool decision=false;
   if ((selection==UNCUT)||(selection==PHOTON)){
      decision = (pa.passElectronVeto) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isLoose) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3);
   }
   if(selection==SEL || selection==DILEP ||selection==ONZ || selection==TRIGDILEP || selection==TRIGSEL || selection==TRIGONZ){
      decision = (pa.passElectronVeto) && (pa.p.Pt()>20.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isLoose) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3); //study Delta R cut ;
   }
   return decision;
}
bool HistogramProducer::testSelection(const selJet& pa, selectionType selection){
   bool decision=false;
   if ((selection==UNCUT)||(selection==PHOTON)){
      decision = (fabs(pa.p.Eta())<2.4) && (pa.isLoose) && !(pa.hasElectronMatch) && !(pa.hasPhotonMatch) && !(pa.hasMuonMatch);
   }
   if(selection==SEL||selection==DILEP || selection==ONZ || selection==TRIGDILEP || selection==TRIGSEL || selection==TRIGONZ){
      decision = (pa.p.Pt()>30.) && (fabs(pa.p.Eta())<2.4) && (pa.isLoose) && !(pa.hasElectronMatch) && !(pa.hasPhotonMatch) && !(pa.hasMuonMatch); //&& (pa.deltaR1>0.3) && (pa.deltaR2>0.3) //study Delta R cut ;
   }
   return decision;
}

//bool matchGenParticle(const tree::Particle& pa){
   //float deltaR_temp,pt_temp;
   //for{genParticles::iterator = genParticles.begin(), iterator}
   //
//}


void HistogramProducer::Init(TTree *tree)
{
   fReader.SetTree(tree);
   inputName = fReader.GetTree()->GetCurrentFile()->GetName();
   cutFlow = *((TH1F*)fReader.GetTree()->GetCurrentFile()->Get("TreeWriter/hCutFlow"));
   fReader.GetEntries(true); // jumps to last file
   string lastInputName = fReader.GetTree()->GetCurrentFile()->GetName();
  
   nEntries=fReader.GetTree()->GetEntries();
   //cout<<"Processing "<<nEntries<<" Events!"<<endl;
   isData = inputName.find("Run201") != string::npos;
   isSignal = inputName.find("SMS") != string::npos;
  
  setHistoNames();
  
   if (inputName!=lastInputName) {
      // adds cut flow of last file. This makes only sense if there is for 1 or two files
      cutFlow.Add((TH1F*)fReader.GetTree()->GetCurrentFile()->Get("TreeWriter/hCutFlow"));
   }
   float nGen = cutFlow.GetBinContent(2);


   noPromptPhotons = inputName.find("DYJets") != string::npos || inputName.find("TTTo2L2Nu") != string::npos;
   //noPromptPhotons = inputName.find("DYJets") != string::npos ;


   InitAllHistos();
   InitTriggerStudiesHistos();
   
   //InitScaleFactors();
   InitScaleFactorsAlternative();

}

void HistogramProducer::InitScaleFactors(){
   DiEleWeighterID=Weighter("scaleFactors/ScaleFactorsElectrons.root","GsfElectronToMVATightTightIP2DSIP3D4");// |eta| vs pt
   DiEleWeighterID.fillOverflow2d();
   DiEleWeighterIso=Weighter("scaleFactors/ScaleFactorsElectrons.root","MVAVLooseElectronToMini");// |eta| vs pt
   DiEleWeighterIso.fillOverflow2d();
   DiEleWeighterConv=Weighter("scaleFactors/ScaleFactorsElectrons.root","MVATightElectronToConvVetoIHit0"); // |eta| vs pt
   DiEleWeighterConv.fillOverflow2d();
   DiEleWeighterTrack=Weighter("scaleFactors/TrackScaleFactorsElectrons.root","EGamma_SF2D"); //pt vs |eta|
   DiEleWeighterTrack.fillOverflow2d();
   
   DiMuWeighterID = Weighter("scaleFactors/ScaleFactorMuonID.root","SF");// |eta| vs pt
   DiMuWeighterID.fillOverflow2d();
   DiMuWeighterIso = Weighter("scaleFactors/ScaleFactorMuonMiniIso.root","SF");// |eta| vs pt
   DiMuWeighterIso.fillOverflow2d();
   DiMuWeighterIP2D = Weighter("scaleFactors/ScaleFactorMuonIP2D.root","SF");// |eta| vs pt
   DiMuWeighterIP2D.fillOverflow2d();
   DiMuWeighterSIP3D = Weighter("scaleFactors/ScaleFactorMuonSIP3D.root","SF");// |eta| vs pt
   DiMuWeighterSIP3D.fillOverflow2d();
   DiMuWeighterTrack = Weighter("scaleFactors/TrackScaleFactorsMuons.root","muonTrackScaleFactorEtaHisto");// |eta|
   DiMuWeighterTrack.fillOverflow2d();
   
   FastSimDiEleWeighterID=Weighter("scaleFactors/FastSimScaleFactorElectronID.root","histo2D");// |eta| vs pt
   FastSimDiEleWeighterID.fillOverflow2d();
   FastSimDiEleWeighterIso=Weighter("scaleFactors/FastSimScaleFactorElectronIso.root","histo2D");// |eta| vs pt
   FastSimDiEleWeighterIso.fillOverflow2d();
   FastSimDiEleWeighterConv=Weighter("scaleFactors/FastSimScaleFactorElectronConvVeto.root","histo2D");// |eta| vs pt
   FastSimDiEleWeighterConv.fillOverflow2d();
   
   FastSimDiMuWeighterID=Weighter("scaleFactors/FastSimScaleFactorMuonID.root","histo2D");// |eta| vs pt
   FastSimDiMuWeighterID.fillOverflow2d();
   FastSimDiMuWeighterIso=Weighter("scaleFactors/FastSimScaleFactorMuonIso.root","histo2D");// |eta| vs pt
   FastSimDiMuWeighterIso.fillOverflow2d();
   FastSimDiMuWeighterIP2D=Weighter("scaleFactors/FastSimScaleFactorMuonIP2D.root","histo2D");// |eta| vs pt
   FastSimDiMuWeighterIP2D.fillOverflow2d();
   FastSimDiMuWeighterSIP3D=Weighter("scaleFactors/FastSimScaleFactorMuonSIP3D.root","histo2D");// |eta| vs pt
   FastSimDiMuWeighterSIP3D.fillOverflow2d();
   }
   
void HistogramProducer::InitScaleFactorsAlternative(){
   DiEleWeighterID=Weighter("scaleFactors/Alternative/ScaleFactorElectronID.root","EGamma_SF2D");// pt vs |eta|
   DiEleWeighterID.fillOverflow2d();
   DiEleWeighterTrack=Weighter("scaleFactors/Alternative/egammaEffi.txt_EGM2D.root","EGamma_SF2D"); //pt vs |eta|
   DiEleWeighterTrack.fillOverflow2d();
   
   DiMuWeighterID_BCDEF = Weighter("scaleFactors/Alternative/EfficienciesAndSF_BCDEF.root","MC_NUM_MediumID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio");// |eta| vs pt
   DiMuWeighterID_BCDEF.fillOverflow2d();
   DiMuWeighterID_GH = Weighter("scaleFactors/Alternative/EfficienciesAndSF_GH.root","MC_NUM_MediumID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio");// |eta| vs pt
   DiMuWeighterID_GH.fillOverflow2d();
   DiMuWeighterIso_BCDEF = Weighter("scaleFactors/Alternative/IsoEfficienciesAndSF_BCDEF.root","TightISO_MediumID_pt_eta/pt_abseta_ratio");// |eta| vs pt
   DiMuWeighterIso_BCDEF.fillOverflow2d();
   DiMuWeighterIso_GH = Weighter("scaleFactors/Alternative/IsoEfficienciesAndSF_GH.root","TightISO_MediumID_pt_eta/pt_abseta_ratio");// |eta| vs pt
   DiMuWeighterIso_GH.fillOverflow2d();
   DiMuWeighterTrack = Weighter("scaleFactors/TrackScaleFactorsMuons.root","muonTrackScaleFactorEtaHisto");// |eta|
   DiMuWeighterTrack.fillOverflow2d();
   }

float HistogramProducer::GetScaleFactorAndError(float pt, float eta,bool isFastSim,bool isEle){
   //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSLeptonSF#Data_leading_order_FullSim_M_AN1
   float DiEleSF =1.;
   float DiMuSF = 1.;
   float DiEleSFErr = 0.;
   float DiMuSFErr = 0.;
   float tempPt,tempPtTrack,tempEta;
   float absEta=fabs(eta);

   float FastSimDiEleSF = 1.;
   float FastSimDiEleSFErr = 0.;
   float FastSimDiMuSF = 1.;
   float FastSimDiMuSFErr = 0.;

   if(isEle){
      DiEleSF = DiEleSF * DiEleWeighterID.getWeight(pt,absEta);
      DiEleSF = DiEleSF * DiEleWeighterIso.getWeight(pt,absEta);
      DiEleSF = DiEleSF * DiEleWeighterConv.getWeight(pt,absEta);
      DiEleSF = DiEleSF * DiEleWeighterTrack.getWeight(absEta,pt); 
      
      DiEleSFErr = DiEleWeighterID.getError(pt,absEta);
      DiEleSFErr = pow((pow(DiEleSFErr,2.) + pow(DiEleWeighterIso.getError(pt,absEta),2.)),0.5);
      DiEleSFErr = pow((pow(DiEleSFErr,2.) + pow(DiEleWeighterConv.getError(pt,absEta),2.)),0.5);
      DiEleSFErr = pow((pow(DiEleSFErr,2.) + pow(DiEleWeighterTrack.getError(absEta,pt),2.)),0.5);
      DiEleSFErr = pow((pow(DiEleSFErr,2.) + pow((0.01*DiEleSF*(pt>20. && pt<80.)),2.)),0.5); //additional 1% uncertainty for all (normal just pt<20 and >80)
      
      if(isFastSim){
         FastSimDiEleSF = DiEleSF * FastSimDiEleWeighterID.getWeight(pt,absEta);
         FastSimDiEleSF = FastSimDiEleSF * FastSimDiEleWeighterIso.getWeight(pt,absEta);
         FastSimDiEleSF = FastSimDiEleSF * FastSimDiEleWeighterConv.getWeight(pt,absEta);

         FastSimDiEleSFErr = pow((pow(DiEleSFErr,2.) + pow(FastSimDiEleWeighterID.getError(pt,absEta),2.)),0.5);
         FastSimDiEleSFErr = pow((pow(FastSimDiEleSFErr,2.) + pow(FastSimDiEleWeighterIso.getError(pt,absEta),2.)),0.5);
         FastSimDiEleSFErr = pow((pow(FastSimDiEleSFErr,2.) + pow(FastSimDiEleWeighterConv.getError(absEta,pt),2.)),0.5);
         FastSimDiEleSFErr = pow((pow(FastSimDiEleSFErr,2.) + pow((0.02*FastSimDiEleSF),2.)),0.5); //additional 2% uncertainty
   
      }
   }else{
      DiMuSF = DiMuSF * DiMuWeighterID.getWeight(pt,absEta);
      DiMuSF = DiMuSF * DiMuWeighterIso.getWeight(pt,absEta);
      DiMuSF = DiMuSF * DiMuWeighterIP2D.getWeight(pt,absEta);
      DiMuSF = DiMuSF * DiMuWeighterSIP3D.getWeight(absEta,pt); 
      DiMuSF = DiMuSF * DiMuWeighterTrack.getWeight(absEta); 
      
      DiMuSFErr = DiMuWeighterID.getError(pt,absEta);
      DiMuSFErr = pow((pow(DiMuSFErr,2.) + pow(DiMuWeighterIso.getError(pt,absEta),2.)),0.5);
      DiMuSFErr = pow((pow(DiMuSFErr,2.) + pow(DiMuWeighterIP2D.getError(pt,absEta),2.)),0.5);
      DiMuSFErr = pow((pow(DiMuSFErr,2.) + pow(DiMuWeighterSIP3D.getError(absEta,pt),2.)),0.5);
      DiMuSFErr = pow((pow(DiMuSFErr,2.) + pow(DiMuWeighterTrack.getError(absEta),2.)),0.5);
      DiMuSFErr = pow((pow(DiMuSFErr,2.) + pow((0.03*DiMuSF),2.)),0.5); //additional 3% 
      
      if(isFastSim){
         FastSimDiMuSF = DiMuSF * FastSimDiMuWeighterID.getWeight(pt,absEta);
         FastSimDiMuSF = FastSimDiMuSF * FastSimDiMuWeighterIso.getWeight(pt,absEta);
         FastSimDiMuSF = FastSimDiMuSF * FastSimDiMuWeighterIP2D.getWeight(pt,absEta);
         FastSimDiMuSF = FastSimDiMuSF * FastSimDiMuWeighterSIP3D.getWeight(pt,absEta);

         FastSimDiMuSFErr = pow((pow(DiMuSFErr,2.) + pow(FastSimDiMuWeighterID.getError(pt,absEta),2.)),0.5);
         FastSimDiMuSFErr = pow((pow(FastSimDiMuSFErr,2.) + pow(FastSimDiMuWeighterIso.getError(pt,absEta),2.)),0.5);
         FastSimDiMuSFErr = pow((pow(FastSimDiMuSFErr,2.) + pow(FastSimDiMuWeighterIP2D.getError(absEta,pt),2.)),0.5);   
         FastSimDiMuSFErr = pow((pow(FastSimDiMuSFErr,2.) + pow(FastSimDiMuWeighterSIP3D.getError(absEta,pt),2.)),0.5);
         FastSimDiMuSFErr = pow((pow(FastSimDiMuSFErr,2.) + pow((0.02*FastSimDiMuSF),2.)),0.5); //additional 2% uncertainty
      }
   }
   float SF=1.;
   if(isFastSim){
      SF = isEle ?  FastSimDiEleSF : FastSimDiMuSF;
   }else{
      SF = isEle ?  DiEleSF : DiMuSF;
   }
   return SF;
}
float HistogramProducer::GetScaleFactorAndErrorAlternative(float pt, float eta,bool isFastSim,bool isEle,int runNr){
   //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSLeptonSF#Data_leading_order_FullSim_M_AN1
   float DiEleSF =1.;
   float DiMuSF_BCDEF = 1.;
   float DiMuSF_GH = 1.;
   float DiEleSFErr = 0.;
   float DiMuSFErr = 0.;
   float tempPt,tempPtTrack,tempEta;
   float absEta=fabs(eta);



   if(isEle){
      DiEleSF = DiEleSF * DiEleWeighterID.getWeight(absEta,pt);
      DiEleSF = DiEleSF * DiEleWeighterTrack.getWeight(absEta,pt); 
      
      //DiEleSFErr = DiEleWeighterID.getError(pt,absEta);
      //DiEleSFErr = pow((pow(DiEleSFErr,2.) + pow(DiEleWeighterIso.getError(pt,absEta),2.)),0.5);
      //DiEleSFErr = pow((pow(DiEleSFErr,2.) + pow(DiEleWeighterConv.getError(pt,absEta),2.)),0.5);
      //DiEleSFErr = pow((pow(DiEleSFErr,2.) + pow(DiEleWeighterTrack.getError(absEta,pt),2.)),0.5);
      //DiEleSFErr = pow((pow(DiEleSFErr,2.) + pow((0.01*DiEleSF*(pt>20. && pt<80.)),2.)),0.5); //additional 1% uncertainty for all (normal just pt<20 and >80)
      
   }else{
      DiMuSF_BCDEF = DiMuSF_BCDEF * DiMuWeighterID_BCDEF.getWeight(pt,absEta);
      DiMuSF_BCDEF = DiMuSF_BCDEF * DiMuWeighterIso_BCDEF.getWeight(pt,absEta);
      DiMuSF_BCDEF = DiMuSF_BCDEF * DiMuWeighterTrack.getWeight(absEta); 
      DiMuSF_GH = DiMuSF_GH * DiMuWeighterID_GH.getWeight(pt,absEta);
      DiMuSF_GH = DiMuSF_GH * DiMuWeighterIso_GH.getWeight(pt,absEta);
      DiMuSF_GH = DiMuSF_GH * DiMuWeighterTrack.getWeight(absEta); 
      
      //DiMuSFErr = DiMuWeighterID.getError(pt,absEta);
      //DiMuSFErr = pow((pow(DiMuSFErr,2.) + pow(DiMuWeighterIso.getError(pt,absEta),2.)),0.5);
      //DiMuSFErr = pow((pow(DiMuSFErr,2.) + pow(DiMuWeighterIP2D.getError(pt,absEta),2.)),0.5);
      //DiMuSFErr = pow((pow(DiMuSFErr,2.) + pow(DiMuWeighterSIP3D.getError(absEta,pt),2.)),0.5);
      //DiMuSFErr = pow((pow(DiMuSFErr,2.) + pow(DiMuWeighterTrack.getError(absEta),2.)),0.5);
      //DiMuSFErr = pow((pow(DiMuSFErr,2.) + pow((0.03*DiMuSF),2.)),0.5); //additional 3% 
      //
   }
   float SF=1.;
   if(isRunABCDEF(runNr)){
      SF = isEle? DiEleSF : DiMuSF_BCDEF;
   }else{
      SF = isEle ?  DiEleSF : DiMuSF_GH;
   }
   return SF;
}


map<Histograms1D,TH1F> HistogramProducer::InitHistograms(const selectionType selection_){
    
   map<Histograms1D,TH1F> hMap;
//enum Histograms1D{PT1,PT2,ETA1,ETA2,PHI1,PHI2,MLL,NJETS,NPHOTONS,ETMISS,HT,GENHT,NVTX,ETAG1,PHIG1,PTG1,SIGMAIETAIETAG1,DeltaEtaLL,DeltaPhiLL,DeltaEtaLLG,DeltaPhiLLG,DeltaRLL,DeltaRLLG,};

   hMap[ETMISS] = TH1F("", ";#it{p}_{T}^{miss} (GeV)", 200, 0, 1000);
   hMap[PT1] = TH1F("", ";#it{p}_{T}^{leading} (GeV)", 200, 0, 1000);
   hMap[PT2] = TH1F("", ";#it{p}_{T}^{trailing} (GeV)", 200, 0, 1000);
   hMap[MLL] = TH1F("", ";#it{m}_{ll} (GeV)", 200, 0, 1000);
   hMap[NPHOTONS] = TH1F("","n_#gamma",10,0,10);
   hMap[NVTX] = TH1F("","n_{Vtx}",60,0,60);
   hMap[HT] = TH1F("","#it{H}_{T} (GeV)",200,0,1000);
   hMap[GENHT] = TH1F("","#it{H}_{T}^{gen} (GeV)",200,0,1000);
   hMap[NJETS] = TH1F("","n_{Jets}",20,0,20);
   hMap[ETA1] = TH1F("", ";|#eta_{trailing}|", 260, 0, 2.6);
   hMap[ETA2] = TH1F("", ";|#eta_{leading}|", 260, 0, 2.6);
   hMap[PHI1] = TH1F("", ";|#phi_{trailing}|", 350, 0, 3.5);
   hMap[PHI2] = TH1F("", ";|#phi_{leading}|", 350, 0, 3.5);
   hMap[DeltaEtaLL] = TH1F("", ";#Delta#Eta_{ll}", 600, 0, 6.);
   hMap[DeltaPhiLL] = TH1F("", ";#Delta#Phi_{ll}", 600, 0, 6.);
   hMap[DeltaRLL] = TH1F("", ";#DeltaR_{ll}", 5000, 0, 6.);

   
   if ((selection_==PHOTON)||(selection_==SEL)||(selection_==ONZ)){
      hMap[PTG1] = TH1F("", ";#it{p}_{T}^{#gamma 1} (GeV)", 200, 0, 1000);
      hMap[ETAG1] = TH1F("", ";|#eta_{#gamma 1}|", 260, 0, 2.6);
      hMap[PHIG1] = TH1F("", ";|#phi_{#gamma 1}|", 350, 0, 3.5);
      hMap[SIGMAIETAIETAG1] = TH1F("", ";#sigma_{i#etai#eta}^{#gamma 1}", 400, 0, 0.04);
      hMap[DeltaRLLG] = TH1F("", ";#DeltaR_{ll,#gamma}", 6000, 0, 6.);
      hMap[DeltaEtaLLG] = TH1F("", ";#Delta#Eta_{ll,#gamma}", 600, 0, 6.);
      hMap[DeltaPhiLLG] = TH1F("", ";#Delta#Phi_{ll,#gamma}", 600, 0, 6.);
   }

   return hMap;

}

map<Histograms2D,TH2F> HistogramProducer::Init2DHistograms(const selectionType selection_){
    
   map<Histograms2D,TH2F> h2Map;
   
   if ((selection_==PHOTON)||(selection_==SEL)){
      h2Map[ISRVFSR] = TH2F("",";#it{m}_{ll};#it{m}{ll#gamma}",200,0,200,500,0,500);   
   }

   return h2Map;

}

map<Histograms1D,TEfficiency> HistogramProducer::InitTriggerStudies(const selectionType selection_){
    
   map<Histograms1D,TEfficiency> hMap;
    
   hMap[ETMISS] = TEfficiency("", ";#it{p}_{T}^{miss} (GeV)", 200, 0, 1000);
   hMap[PT1] = TEfficiency("", ";#it{p}_{T}^{leading} (GeV)", 1000, 0, 1000);
   hMap[PT2] = TEfficiency("", ";#it{p}_{T}^{trailing} (GeV)", 1000, 0, 1000);
   hMap[MLL] = TEfficiency("", ";#it{m}_{ll} (GeV)", 200, 0, 1000);
   hMap[NPHOTONS] = TEfficiency("","n_#gamma",10,0,10);
   hMap[NVTX] = TEfficiency("","n_{Vtx}",60,0,60);
   hMap[HT] = TEfficiency("","#it{H}_{T} (GeV)",200,0,1000);
   hMap[GENHT] = TEfficiency("","#it{H}_{T}^{gen} (GeV)",200,0,1000);
   hMap[NJETS] = TEfficiency("","n_{Jets}",20,0,20);
   hMap[ETA1] = TEfficiency("", ";|#eta_{trailing}|", 260, 0, 2.6);
   hMap[ETA2] = TEfficiency("", ";|#eta_{leading}|", 260, 0, 2.6);
   hMap[PHI1] = TEfficiency("", ";|#phi_{trailing}|", 350, 0, 3.5);
   hMap[PHI2] = TEfficiency("", ";|#phi_{leading}|", 350, 0, 3.5);
   
   if ((selection_==TRIGONZ)||(selection_==TRIGSEL)){
      hMap[PTG1] = TEfficiency("", ";#it{p}_{T}^{#gamma 1} (GeV)", 200, 0, 1000);
      hMap[ETAG1] = TEfficiency("", ";|#eta_{#gamma 1}|", 260, 0, 2.6);
      hMap[PHIG1] = TEfficiency("", ";|#phi_{#gamma 1}|", 350, 0, 3.5);
      hMap[SIGMAIETAIETAG1] = TEfficiency("", ";#sigma_{i#etai#eta}^{#gamma 1}", 400, 0, 0.04);
   }

   return hMap;

}




void HistogramProducer::SlaveBegin(TTree *tree)
{
}

bool HistogramProducer::SelectEvent(selectionType selection){
   if (CheckParticles()){//only Fill Histograms if nMu>=2, nEle>=2, nGamma>=0
      if (Cleaning()){//decide in events with DiEle and DiMu where to contribute
         selEvent L;
         L.isDiElectron=isDiElectron;
         L.trigDiEle=trigDiEle;
         L.trigDiMu=trigDiMu;
         L.trigMuEle=trigMuEle;
         L.trigHt=trigHt;
         //if(trigHt&&isDiElectron){
            //cout<<"  "<<trigDiEle<<endl;
         //}
         L.l1=lep1;
         L.l2=lep2;
         L.pt1=pt1;
         L.pt2=pt2;
         L.phi1=phi1;
         L.phi2=phi2;
         L.eta1=eta1;
         L.eta2=eta2;
         L.chargeProduct=chargeProduct;
         L.deltaRll=deltaRll;
         L.mll=mll;
         L.miniIso1=miniIso1;
         L.miniIso2=miniIso2;
         L.ETmiss=ETmiss;
         if(!isData){
            //L.totalWeight = totalWeight * GetScaleFactorAndError(pt1,eta1,isSignal,isDiElectron)*GetScaleFactorAndError(pt2,eta2,isSignal,isDiElectron);
            L.totalWeight = totalWeight * GetScaleFactorAndErrorAlternative(pt1,eta1,isSignal,isDiElectron,*runNo)*GetScaleFactorAndErrorAlternative(pt2,eta2,isSignal,isDiElectron,*runNo);
         }else{
            L.totalWeight=totalWeight;
         }
         if (L.isDiElectron){
            auto e1 = electrons->at(0);
            auto e2 = electrons->at(1);
            if(testSelection(e1,selection,true) && testSelection(e2,selection,false) && ((selection==UNCUT)? true : (mll>50.)) ){

               for (vector<tree::Photon>::iterator it = photons->begin(); it != photons->end(); ++it){
                  auto g = *it;
                  selPhoton gamma;
                  gamma.setAll(g);
                  gamma.vec.SetPtEtaPhiM(gamma.p.Pt(),gamma.p.Eta(),gamma.p.Phi(),0.);
                  gamma.deltaR1=gamma.vec.DeltaR(L.l1);
                  gamma.deltaR2=gamma.vec.DeltaR(L.l2);
                  if (testSelection(gamma,selection)){
                     L.selPhotons.push_back(gamma);
                  }
               }
               for (vector<tree::Jet>::iterator it = jets->begin(); it != jets->end(); ++it){
                  auto j = *it;
                  selJet jet;
                  jet.setAll(j);
                  jet.vec.SetPtEtaPhiM(jet.p.Pt(),jet.p.Eta(),jet.p.Phi(),0.);
                  jet.deltaR1=jet.vec.DeltaR(L.l1);
                  jet.deltaR2=jet.vec.DeltaR(L.l2);
                  if (testSelection(jet,selection)){
                     L.selJets.push_back(jet);
                  }
               }
               selectedEvent=L;
               if(!GenPhotonVeto()){
                  return true;
               }else{
                  return false;
               }
            }
         }
         else{
            auto m1 = muons->at(0); 
            auto m2 = muons->at(1);
            if(testSelection(m1,selection,true) && testSelection(m2,selection,false) && ((selection==UNCUT)? true : (mll>50.))){
               for (vector<tree::Photon>::iterator it = photons->begin(); it != photons->end(); ++it){
                  auto g = *it;
                  selPhoton gamma;
                  gamma.setAll(g);
                  gamma.vec.SetPtEtaPhiM(gamma.p.Pt(),gamma.p.Eta(),gamma.p.Phi(),0.);
                  gamma.deltaR1=gamma.vec.DeltaR(L.l1);
                  gamma.deltaR2=gamma.vec.DeltaR(L.l2);
                  if (testSelection(gamma,selection)){
                     L.selPhotons.push_back(gamma);
                  }
               }
               for (vector<tree::Jet>::iterator it = jets->begin(); it != jets->end(); ++it){
                  auto j = *it;
                  selJet jet;
                  jet.setAll(j);
                  jet.vec.SetPtEtaPhiM(jet.p.Pt(),jet.p.Eta(),jet.p.Phi(),0.);
                  jet.deltaR1=jet.vec.DeltaR(L.l1);
                  jet.deltaR2=jet.vec.DeltaR(L.l2);
                  if (testSelection(jet,selection)){
                     L.selJets.push_back(jet);
                  }
               }
               selectedEvent=L;
               if(!GenPhotonVeto()){
                  return true;
               }else{
                  return false;
               }
            }
         }
      }else{return false;}
   }else{return false;}
  
  
  return false;
}

bool HistogramProducer::SelectEventTriggerStudies(selectionType selection){
   if (CheckParticles()){//only Fill Histograms if nMu>=2, nEle>=2, nGamma>=0
      if (CleaningTriggerStudies()){//decide in events with DiEle and DiMu where to contribute
         selEvent L;
         L.isDiElectron=isDiElectron;
         L.trigDiEle=trigDiEle;
         L.trigDiMu=trigDiMu;
         L.trigMuEle=trigMuEle;
         L.trigHt=trigHt;
         //if(trigHt&&isDiElectron){
            //cout<<"  "<<trigDiEle<<endl;
         //}
         L.l1=lep1;
         L.l2=lep2;
         L.pt1=pt1;
         L.pt2=pt2;
         L.phi1=phi1;
         L.phi2=phi2;
         L.eta1=eta1;
         L.eta2=eta2;
         L.chargeProduct=chargeProduct;
         L.deltaRll=deltaRll;
         L.mll=mll;
         L.miniIso1=miniIso1;
         L.miniIso2=miniIso2;
         L.ETmiss=ETmiss;
         if(!isData){
            //L.totalWeight = totalWeight * GetScaleFactorAndError(pt1,eta1,isSignal,isDiElectron)*GetScaleFactorAndError(pt2,eta2,isSignal,isDiElectron);
            L.totalWeight = totalWeight * GetScaleFactorAndErrorAlternative(pt1,eta1,isSignal,isDiElectron,*runNo)*GetScaleFactorAndErrorAlternative(pt2,eta2,isSignal,isDiElectron,*runNo);
         }else{
            L.totalWeight=totalWeight;
         }
         if (L.isDiElectron){
            auto e1 = electrons->at(0);
            auto e2 = electrons->at(1);
            if(testSelection(e1,selection,true) && testSelection(e2,selection,false) && ((selection==UNCUT)? true : (mll>50.)) ){

               for (vector<tree::Photon>::iterator it = photons->begin(); it != photons->end(); ++it){
                  auto g = *it;
                  selPhoton gamma;
                  gamma.setAll(g);
                  gamma.vec.SetPtEtaPhiM(gamma.p.Pt(),gamma.p.Eta(),gamma.p.Phi(),0.);
                  gamma.deltaR1=gamma.vec.DeltaR(L.l1);
                  gamma.deltaR2=gamma.vec.DeltaR(L.l2);
                  if (testSelection(gamma,selection)){
                     L.selPhotons.push_back(gamma);
                  }
               }
               for (vector<tree::Jet>::iterator it = jets->begin(); it != jets->end(); ++it){
                  auto j = *it;
                  selJet jet;
                  jet.setAll(j);
                  jet.vec.SetPtEtaPhiM(jet.p.Pt(),jet.p.Eta(),jet.p.Phi(),0.);
                  jet.deltaR1=jet.vec.DeltaR(L.l1);
                  jet.deltaR2=jet.vec.DeltaR(L.l2);
                  if (testSelection(jet,selection)){
                     L.selJets.push_back(jet);
                  }
               }
               selectedEvent=L;
               if(!GenPhotonVeto()){
                  return true;
               }else{
                  return false;
               }
            }
         }
         else{
            auto m1 = muons->at(0); 
            auto m2 = muons->at(1);
            if(testSelection(m1,selection,true) && testSelection(m2,selection,false) && ((selection==UNCUT)? true : (mll>50.))){
               for (vector<tree::Photon>::iterator it = photons->begin(); it != photons->end(); ++it){
                  auto g = *it;
                  selPhoton gamma;
                  gamma.setAll(g);
                  gamma.vec.SetPtEtaPhiM(gamma.p.Pt(),gamma.p.Eta(),gamma.p.Phi(),0.);
                  gamma.deltaR1=gamma.vec.DeltaR(L.l1);
                  gamma.deltaR2=gamma.vec.DeltaR(L.l2);
                  if (testSelection(gamma,selection)){
                     L.selPhotons.push_back(gamma);
                  }
               }
               for (vector<tree::Jet>::iterator it = jets->begin(); it != jets->end(); ++it){
                  auto j = *it;
                  selJet jet;
                  jet.setAll(j);
                  jet.vec.SetPtEtaPhiM(jet.p.Pt(),jet.p.Eta(),jet.p.Phi(),0.);
                  jet.deltaR1=jet.vec.DeltaR(L.l1);
                  jet.deltaR2=jet.vec.DeltaR(L.l2);
                  if (testSelection(jet,selection)){
                     L.selJets.push_back(jet);
                  }
               }
               selectedEvent=L;
               if(!GenPhotonVeto()){
                  return true;
               }else{
                  return false;
               }
            }
         }
      }else{return false;}
   }else{return false;}
  
  
  return false;
}


void HistogramProducer::InitAllHistos(){
   //h1Maps["uncutEE"]=InitHistograms("uncut"); 
   //h1Maps["uncutMM"]=InitHistograms("uncut");
   //h1Maps["uncut"]=InitHistograms("uncut");
   h1Maps["dilepEE"]=InitHistograms(DILEP); 
   h1Maps["dilepMM"]=InitHistograms(DILEP);
   h1Maps["dilep"]=InitHistograms(DILEP);
   //h1Maps["1photonEE"]=InitHistograms("1photon"); 
   //h1Maps["1photonMM"]=InitHistograms("1photon"); 
   //h1Maps["1photon"]=InitHistograms("1photon"); 
   h1Maps["selEE"]=InitHistograms(SEL);
   h1Maps["selMM"]=InitHistograms(SEL); 
   h1Maps["sel"]=InitHistograms(SEL); 
   h1Maps["onZEE"]=InitHistograms(ONZ);
   h1Maps["onZMM"]=InitHistograms(ONZ); 
   h1Maps["onZ"]=InitHistograms(ONZ); 
   //h2Maps["selEE"]=Init2DHistograms("sel");
   //h2Maps["selMM"]=Init2DHistograms("sel"); 
   //h2Maps["sel"]=Init2DHistograms("sel"); 
}
   
void HistogramProducer::InitTriggerStudiesHistos(){
   eff1Maps["trigDilep"]=InitTriggerStudies(TRIGDILEP);
   eff1Maps["trigDilepEE"]=InitTriggerStudies(TRIGDILEP);
   eff1Maps["trigDilepMM"]=InitTriggerStudies(TRIGDILEP);
   eff1Maps["trigSelEE"]=InitTriggerStudies(TRIGSEL);
   eff1Maps["trigSelMM"]=InitTriggerStudies(TRIGSEL);
   eff1Maps["trigOnZEE"]=InitTriggerStudies(TRIGONZ);
   eff1Maps["trigOnZMM"]=InitTriggerStudies(TRIGONZ);

}

void HistogramProducer::FillHistograms2D(){
   //if(SelectEvent(SEL)){
   if(SelectEvent(SEL)||SelectEvent(ONZ)){
      if (selectedEvent.selPhotons.size()!=0){
         auto m1 = &h2Maps["selEE"];
         auto m2 = &h2Maps["selMM"];
         auto m3 = &h2Maps["sel"];
         m3->at(ISRVFSR).Fill(selectedEvent.mll,(selectedEvent.l1+selectedEvent.l2+selectedEvent.selPhotons.at(0).vec).M());         
         if (selectedEvent.isDiElectron){ 
            m1->at(ISRVFSR).Fill(selectedEvent.mll,(selectedEvent.l1+selectedEvent.l2+selectedEvent.selPhotons.at(0).vec).M());    
         }else{
            m2->at(ISRVFSR).Fill(selectedEvent.mll,(selectedEvent.l1+selectedEvent.l2+selectedEvent.selPhotons.at(0).vec).M());          
         }
      }   
   
   }
}
void HistogramProducer::FillHistograms(){
//enum Histograms1D{PT1,PT2,ETA1,ETA2,PHI1,PHI2,MLL,NJETS,NPHOTONS,ETMISS,HT,GENHT,NVTX,ETAG1,PHIG1,PTG1,SIGMAIETAIETAG1,DeltaEtaLL,DeltaPhiLL,DeltaEtaLLG,DeltaPhiLLG,DeltaRLL,DeltaRLLG,};

   /*if(SelectEvent(UNCUT)){
      auto m1 = &h1Maps["uncutEE"];
      auto m2 = &h1Maps["uncutMM"];
      auto m3 = &h1Maps["uncut"];
      m3->at("met").Fill(selectedEvent.ETmiss, selectedEvent.totalWeight);
      m3->at("pt1").Fill(selectedEvent.pt1, selectedEvent.totalWeight);
      m3->at("pt2").Fill(selectedEvent.pt2, selectedEvent.totalWeight);
      m3->at("m_ll").Fill(selectedEvent.mll,selectedEvent.totalWeight);   
      m3->at("n_photons").Fill(selectedEvent.selPhotons.size(),selectedEvent.totalWeight);
      m3->at("n_vtx").Fill(*nGoodVertices,selectedEvent.totalWeight);
      m3->at("ht").Fill(*ht,selectedEvent.totalWeight);
      m3->at("gen_ht").Fill(*genHt,selectedEvent.totalWeight);
      m3->at("n_jets").Fill(selectedEvent.selJets.size(),selectedEvent.totalWeight);
      m3->at("eta1").Fill(fabs(selectedEvent.eta1),selectedEvent.totalWeight);
      m3->at("eta2").Fill(fabs(selectedEvent.eta2),selectedEvent.totalWeight);
      m3->at("phi1").Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
      m3->at("phi2").Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight); 
      if (selectedEvent.isDiElectron){ 
         m1->at("met").Fill(selectedEvent.ETmiss, selectedEvent.totalWeight);
         m1->at("pt1").Fill(selectedEvent.pt1, selectedEvent.totalWeight);
         m1->at("pt2").Fill(selectedEvent.pt2, selectedEvent.totalWeight);
         m1->at("m_ll").Fill(selectedEvent.mll,selectedEvent.totalWeight);   
         m1->at("n_photons").Fill(selectedEvent.selPhotons.size(),selectedEvent.totalWeight);
         m1->at("n_vtx").Fill(*nGoodVertices,selectedEvent.totalWeight);
         m1->at("ht").Fill(*ht,selectedEvent.totalWeight);
         m1->at("gen_ht").Fill(*genHt,selectedEvent.totalWeight);
         m1->at("n_jets").Fill(selectedEvent.selJets.size(),selectedEvent.totalWeight);
         m1->at("eta1").Fill(fabs(selectedEvent.eta1),selectedEvent.totalWeight);
         m1->at("eta2").Fill(fabs(selectedEvent.eta2),selectedEvent.totalWeight);
         m1->at("phi1").Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
         m1->at("phi2").Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
      }else{
         m2->at("met").Fill(selectedEvent.ETmiss, selectedEvent.totalWeight);
         m2->at("pt1").Fill(selectedEvent.pt1, selectedEvent.totalWeight);
         m2->at("pt2").Fill(selectedEvent.pt2, selectedEvent.totalWeight);
         m2->at("m_ll").Fill(selectedEvent.mll,selectedEvent.totalWeight);   
         m2->at("n_photons").Fill(selectedEvent.selPhotons.size(),selectedEvent.totalWeight);
         m2->at("n_vtx").Fill(*nGoodVertices,selectedEvent.totalWeight);
         m2->at("ht").Fill(*ht,selectedEvent.totalWeight);
         m2->at("gen_ht").Fill(*genHt,selectedEvent.totalWeight);
         m2->at("n_jets").Fill(selectedEvent.selJets.size(),selectedEvent.totalWeight);
         m2->at("eta1").Fill(fabs(selectedEvent.eta1),selectedEvent.totalWeight);
         m2->at("eta2").Fill(fabs(selectedEvent.eta2),selectedEvent.totalWeight);
         m2->at("phi1").Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
         m2->at("phi2").Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);      
      }
   }*/
   if(SelectEvent(DILEP)){
      auto m1 = &h1Maps["dilepEE"];
      auto m2 = &h1Maps["dilepMM"];
      auto m3 = &h1Maps["dilep"];
      m3->at(ETMISS).Fill(selectedEvent.ETmiss, selectedEvent.totalWeight);
      m3->at(PT1).Fill(selectedEvent.pt1, selectedEvent.totalWeight);
      m3->at(PT2).Fill(selectedEvent.pt2, selectedEvent.totalWeight);
      m3->at(MLL).Fill(selectedEvent.mll,selectedEvent.totalWeight);   
      m3->at(NPHOTONS).Fill(selectedEvent.selPhotons.size(),selectedEvent.totalWeight);
      m3->at(NVTX).Fill(*nGoodVertices,selectedEvent.totalWeight);
      m3->at(HT).Fill(*ht,selectedEvent.totalWeight);
      m3->at(GENHT).Fill(*genHt,selectedEvent.totalWeight);
      m3->at(NJETS).Fill(selectedEvent.selJets.size(),selectedEvent.totalWeight);
      m3->at(ETA1).Fill(fabs(selectedEvent.eta1),selectedEvent.totalWeight);
      m3->at(ETA2).Fill(fabs(selectedEvent.eta2),selectedEvent.totalWeight);
      m3->at(PHI1).Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
      m3->at(PHI2).Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
      m3->at(DeltaEtaLL).Fill(fabs(selectedEvent.eta1-selectedEvent.eta2),selectedEvent.totalWeight);
      m3->at(DeltaPhiLL).Fill(fabs(selectedEvent.phi1-selectedEvent.phi2),selectedEvent.totalWeight);
      m3->at(DeltaRLL).Fill(fabs(selectedEvent.deltaRll),selectedEvent.totalWeight);
      if (selectedEvent.isDiElectron){ 
        m1->at(ETMISS).Fill(selectedEvent.ETmiss, selectedEvent.totalWeight);
        m1->at(PT1).Fill(selectedEvent.pt1, selectedEvent.totalWeight);
        m1->at(PT2).Fill(selectedEvent.pt2, selectedEvent.totalWeight);
        m1->at(MLL).Fill(selectedEvent.mll,selectedEvent.totalWeight);   
        m1->at(NPHOTONS).Fill(selectedEvent.selPhotons.size(),selectedEvent.totalWeight);
        m1->at(NVTX).Fill(*nGoodVertices,selectedEvent.totalWeight);
        m1->at(HT).Fill(*ht,selectedEvent.totalWeight);
        m1->at(GENHT).Fill(*genHt,selectedEvent.totalWeight);
        m1->at(NJETS).Fill(selectedEvent.selJets.size(),selectedEvent.totalWeight);
        m1->at(ETA1).Fill(fabs(selectedEvent.eta1),selectedEvent.totalWeight);
        m1->at(ETA2).Fill(fabs(selectedEvent.eta2),selectedEvent.totalWeight);
        m1->at(PHI1).Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
        m1->at(PHI2).Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
        m1->at(DeltaEtaLL).Fill(fabs(selectedEvent.eta1-selectedEvent.eta2),selectedEvent.totalWeight);
        m1->at(DeltaPhiLL).Fill(fabs(selectedEvent.phi1-selectedEvent.phi2),selectedEvent.totalWeight);
        m1->at(DeltaRLL).Fill(fabs(selectedEvent.deltaRll),selectedEvent.totalWeight);
      }else{
        m2->at(ETMISS).Fill(selectedEvent.ETmiss, selectedEvent.totalWeight);
        m2->at(PT1).Fill(selectedEvent.pt1, selectedEvent.totalWeight);
        m2->at(PT2).Fill(selectedEvent.pt2, selectedEvent.totalWeight);
        m2->at(MLL).Fill(selectedEvent.mll,selectedEvent.totalWeight);   
        m2->at(NPHOTONS).Fill(selectedEvent.selPhotons.size(),selectedEvent.totalWeight);
        m2->at(NVTX).Fill(*nGoodVertices,selectedEvent.totalWeight);
        m2->at(HT).Fill(*ht,selectedEvent.totalWeight);
        m2->at(GENHT).Fill(*genHt,selectedEvent.totalWeight);
        m2->at(NJETS).Fill(selectedEvent.selJets.size(),selectedEvent.totalWeight);
        m2->at(ETA1).Fill(fabs(selectedEvent.eta1),selectedEvent.totalWeight);
        m2->at(ETA2).Fill(fabs(selectedEvent.eta2),selectedEvent.totalWeight);
        m2->at(PHI1).Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
        m2->at(PHI2).Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
        m2->at(DeltaEtaLL).Fill(fabs(selectedEvent.eta1-selectedEvent.eta2),selectedEvent.totalWeight);
        m2->at(DeltaPhiLL).Fill(fabs(selectedEvent.phi1-selectedEvent.phi2),selectedEvent.totalWeight);
        m2->at(DeltaRLL).Fill(fabs(selectedEvent.deltaRll),selectedEvent.totalWeight);
      }
   }
   /*if(SelectEvent(PHOTON)){
      if (selectedEvent.selPhotons.size()!=0){
         auto m1 = &h1Maps["1photonEE"];
         auto m2 = &h1Maps["1photonMM"];
         auto m3 = &h1Maps["1photon"];
         m3->at("met").Fill(selectedEvent.ETmiss, selectedEvent.totalWeight);
         m3->at("pt1").Fill(selectedEvent.pt1, selectedEvent.totalWeight);
         m3->at("pt2").Fill(selectedEvent.pt2, selectedEvent.totalWeight);
         m3->at("m_ll").Fill(selectedEvent.mll,selectedEvent.totalWeight);   
         m3->at("n_photons").Fill(selectedEvent.selPhotons.size(),selectedEvent.totalWeight);
         m3->at("n_vtx").Fill(*nGoodVertices,selectedEvent.totalWeight);
         m3->at("ht").Fill(*ht,selectedEvent.totalWeight);
         m3->at("gen_ht").Fill(*genHt,selectedEvent.totalWeight);
         m3->at("n_jets").Fill(selectedEvent.selJets.size(),selectedEvent.totalWeight);
         m3->at("eta1").Fill(fabs(selectedEvent.eta1),selectedEvent.totalWeight);
         m3->at("eta2").Fill(fabs(selectedEvent.eta2),selectedEvent.totalWeight);
         m3->at("phi1").Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
         m3->at("phi2").Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
         m3->at("pt_g1").Fill(selectedEvent.selPhotons.at(0).p.Pt(),selectedEvent.totalWeight);
         m3->at("phi_g1").Fill(fabs(selectedEvent.selPhotons.at(0).p.Phi()),selectedEvent.totalWeight);
         m3->at("eta_g1").Fill(fabs(selectedEvent.selPhotons.at(0).p.Eta()),selectedEvent.totalWeight);
         m3->at("sigmaIetaIeta_g1").Fill(selectedEvent.selPhotons.at(0).sigmaIetaIeta,selectedEvent.totalWeight);         
         if (selectedEvent.isDiElectron){ 
            m1->at("met").Fill(selectedEvent.ETmiss, selectedEvent.totalWeight);
            m1->at("pt1").Fill(selectedEvent.pt1, selectedEvent.totalWeight);
            m1->at("pt2").Fill(selectedEvent.pt2, selectedEvent.totalWeight);
            m1->at("m_ll").Fill(selectedEvent.mll,selectedEvent.totalWeight);   
            m1->at("n_photons").Fill(selectedEvent.selPhotons.size(),selectedEvent.totalWeight);
            m1->at("n_vtx").Fill(*nGoodVertices,selectedEvent.totalWeight);
            m1->at("ht").Fill(*ht,selectedEvent.totalWeight);
            m1->at("gen_ht").Fill(*genHt,selectedEvent.totalWeight);
            m1->at("n_jets").Fill(selectedEvent.selJets.size(),selectedEvent.totalWeight);
            m1->at("eta1").Fill(fabs(selectedEvent.eta1),selectedEvent.totalWeight);
            m1->at("eta2").Fill(fabs(selectedEvent.eta2),selectedEvent.totalWeight);
            m1->at("phi1").Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
            m1->at("phi2").Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
            m1->at("pt_g1").Fill(selectedEvent.selPhotons.at(0).p.Pt(),selectedEvent.totalWeight);
            m1->at("phi_g1").Fill(fabs(selectedEvent.selPhotons.at(0).p.Phi()),selectedEvent.totalWeight);
            m1->at("eta_g1").Fill(fabs(selectedEvent.selPhotons.at(0).p.Eta()),selectedEvent.totalWeight);
            m1->at("sigmaIetaIeta_g1").Fill(selectedEvent.selPhotons.at(0).sigmaIetaIeta,selectedEvent.totalWeight);
         }else{
            m2->at("met").Fill(selectedEvent.ETmiss, selectedEvent.totalWeight);
            m2->at("pt1").Fill(selectedEvent.pt1, selectedEvent.totalWeight);
            m2->at("pt2").Fill(selectedEvent.pt2, selectedEvent.totalWeight);
            m2->at("m_ll").Fill(selectedEvent.mll,selectedEvent.totalWeight);   
            m2->at("n_photons").Fill(selectedEvent.selPhotons.size(),selectedEvent.totalWeight);
            m2->at("n_vtx").Fill(*nGoodVertices,selectedEvent.totalWeight);
            m2->at("ht").Fill(*ht,selectedEvent.totalWeight);
            m2->at("gen_ht").Fill(*genHt,selectedEvent.totalWeight);
            m2->at("n_jets").Fill(selectedEvent.selJets.size(),selectedEvent.totalWeight);
            m2->at("eta1").Fill(fabs(selectedEvent.eta1),selectedEvent.totalWeight);
            m2->at("eta2").Fill(fabs(selectedEvent.eta2),selectedEvent.totalWeight);
            m2->at("phi1").Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
            m2->at("phi2").Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
            m2->at("pt_g1").Fill(selectedEvent.selPhotons.at(0).p.Pt(),selectedEvent.totalWeight);
            m2->at("phi_g1").Fill(fabs(selectedEvent.selPhotons.at(0).p.Phi()),selectedEvent.totalWeight);
            m2->at("eta_g1").Fill(fabs(selectedEvent.selPhotons.at(0).p.Eta()),selectedEvent.totalWeight);
            m2->at("sigmaIetaIeta_g1").Fill(selectedEvent.selPhotons.at(0).sigmaIetaIeta,selectedEvent.totalWeight);      
         }
      }
   }*/
   if(SelectEvent(SEL)){
      if (selectedEvent.selPhotons.size()!=0){
         auto m1 = &h1Maps["selEE"];
         auto m2 = &h1Maps["selMM"];
         auto m3 = &h1Maps["sel"];
         m3->at(ETMISS).Fill(selectedEvent.ETmiss, selectedEvent.totalWeight);
         m3->at(PT1).Fill(selectedEvent.pt1, selectedEvent.totalWeight);
         m3->at(PT2).Fill(selectedEvent.pt2, selectedEvent.totalWeight);
         m3->at(MLL).Fill(selectedEvent.mll,selectedEvent.totalWeight);   
         m3->at(NPHOTONS).Fill(selectedEvent.selPhotons.size(),selectedEvent.totalWeight);
         m3->at(NVTX).Fill(*nGoodVertices,selectedEvent.totalWeight);
         m3->at(HT).Fill(*ht,selectedEvent.totalWeight);
         m3->at(GENHT).Fill(*genHt,selectedEvent.totalWeight);
         m3->at(NJETS).Fill(selectedEvent.selJets.size(),selectedEvent.totalWeight);
         m3->at(ETA1).Fill(fabs(selectedEvent.eta1),selectedEvent.totalWeight);
         m3->at(ETA2).Fill(fabs(selectedEvent.eta2),selectedEvent.totalWeight);
         m3->at(PHI1).Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
         m3->at(PHI1).Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
         m3->at(PTG1).Fill(selectedEvent.selPhotons.at(0).p.Pt(),selectedEvent.totalWeight);
         m3->at(PHIG1).Fill(selectedEvent.selPhotons.at(0).p.Phi(),selectedEvent.totalWeight);
         m3->at(ETAG1).Fill(selectedEvent.selPhotons.at(0).p.Eta(),selectedEvent.totalWeight);
         m3->at(SIGMAIETAIETAG1).Fill(selectedEvent.selPhotons.at(0).sigmaIetaIeta,selectedEvent.totalWeight);
         m3->at(DeltaEtaLL).Fill(fabs(selectedEvent.eta1-selectedEvent.eta2),selectedEvent.totalWeight);
         m3->at(DeltaPhiLL).Fill(fabs(selectedEvent.phi1-selectedEvent.phi2),selectedEvent.totalWeight);
         m3->at(DeltaRLL).Fill(fabs(selectedEvent.deltaRll),selectedEvent.totalWeight);
         m3->at(DeltaEtaLLG).Fill(fabs((selectedEvent.l1+selectedEvent.l2).Eta()-selectedEvent.selPhotons.at(0).vec.Eta()), selectedEvent.totalWeight);
         m3->at(DeltaPhiLLG).Fill(fabs((selectedEvent.l1+selectedEvent.l2).Phi()-selectedEvent.selPhotons.at(0).vec.Phi()),selectedEvent.totalWeight);
         m3->at(DeltaRLLG).Fill(fabs((selectedEvent.l1+selectedEvent.l2).DeltaR(selectedEvent.selPhotons.at(0).vec)),selectedEvent.totalWeight);
         if (selectedEvent.isDiElectron){ 
            m1->at(ETMISS).Fill(selectedEvent.ETmiss, selectedEvent.totalWeight);
            m1->at(PT1).Fill(selectedEvent.pt1, selectedEvent.totalWeight);
            m1->at(PT2).Fill(selectedEvent.pt2, selectedEvent.totalWeight);
            m1->at(MLL).Fill(selectedEvent.mll,selectedEvent.totalWeight);   
            m1->at(NPHOTONS).Fill(selectedEvent.selPhotons.size(),selectedEvent.totalWeight);
            m1->at(NVTX).Fill(*nGoodVertices,selectedEvent.totalWeight);
            m1->at(HT).Fill(*ht,selectedEvent.totalWeight);
            m1->at(GENHT).Fill(*genHt,selectedEvent.totalWeight);
            m1->at(NJETS).Fill(selectedEvent.selJets.size(),selectedEvent.totalWeight);
            m1->at(ETA1).Fill(fabs(selectedEvent.eta1),selectedEvent.totalWeight);
            m1->at(ETA2).Fill(fabs(selectedEvent.eta2),selectedEvent.totalWeight);
            m1->at(PHI1).Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
            m1->at(PHI2).Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
            m1->at(PTG1).Fill(selectedEvent.selPhotons.at(0).p.Pt(),selectedEvent.totalWeight);
            m1->at(PHIG1).Fill(fabs(selectedEvent.selPhotons.at(0).p.Phi()),selectedEvent.totalWeight);
            m1->at(ETAG1).Fill(fabs(selectedEvent.selPhotons.at(0).p.Eta()),selectedEvent.totalWeight);
            m1->at(SIGMAIETAIETAG1).Fill(selectedEvent.selPhotons.at(0).sigmaIetaIeta,selectedEvent.totalWeight);
            m1->at(DeltaEtaLL).Fill(fabs(selectedEvent.eta1-selectedEvent.eta2),selectedEvent.totalWeight);
            m1->at(DeltaPhiLL).Fill(fabs(selectedEvent.phi1-selectedEvent.phi2),selectedEvent.totalWeight);
            m1->at(DeltaRLL).Fill(fabs(selectedEvent.deltaRll),selectedEvent.totalWeight);
            m1->at(DeltaEtaLLG).Fill(fabs((selectedEvent.l1+selectedEvent.l2).Eta()-selectedEvent.selPhotons.at(0).vec.Eta()), selectedEvent.totalWeight);
            m1->at(DeltaPhiLLG).Fill(fabs((selectedEvent.l1+selectedEvent.l2).Phi()-selectedEvent.selPhotons.at(0).vec.Phi()),selectedEvent.totalWeight);
            m1->at(DeltaRLLG).Fill(fabs((selectedEvent.l1+selectedEvent.l2).DeltaR(selectedEvent.selPhotons.at(0).vec)),selectedEvent.totalWeight);
         }else{
            m2->at(ETMISS).Fill(selectedEvent.ETmiss, selectedEvent.totalWeight);
            m2->at(PT1).Fill(selectedEvent.pt1, selectedEvent.totalWeight);
            m2->at(PT2).Fill(selectedEvent.pt2, selectedEvent.totalWeight);
            m2->at(MLL).Fill(selectedEvent.mll,selectedEvent.totalWeight);   
            m2->at(NPHOTONS).Fill(selectedEvent.selPhotons.size(),selectedEvent.totalWeight);
            m2->at(NVTX).Fill(*nGoodVertices,selectedEvent.totalWeight);
            m2->at(HT).Fill(*ht,selectedEvent.totalWeight);
            m2->at(GENHT).Fill(*genHt,selectedEvent.totalWeight);
            m2->at(NJETS).Fill(selectedEvent.selJets.size(),selectedEvent.totalWeight);
            m2->at(ETA1).Fill(fabs(selectedEvent.eta1),selectedEvent.totalWeight);
            m2->at(ETA2).Fill(fabs(selectedEvent.eta2),selectedEvent.totalWeight);
            m2->at(PHI1).Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
            m2->at(PHI2).Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
            m2->at(PTG1).Fill(selectedEvent.selPhotons.at(0).p.Pt(),selectedEvent.totalWeight);
            m2->at(PHIG1).Fill(fabs(selectedEvent.selPhotons.at(0).p.Phi()),selectedEvent.totalWeight);
            m2->at(ETAG1).Fill(fabs(selectedEvent.selPhotons.at(0).p.Eta()),selectedEvent.totalWeight);
            m2->at(SIGMAIETAIETAG1).Fill(selectedEvent.selPhotons.at(0).sigmaIetaIeta,selectedEvent.totalWeight);
            m2->at(DeltaEtaLL).Fill(fabs(selectedEvent.eta1-selectedEvent.eta2),selectedEvent.totalWeight);
            m2->at(DeltaPhiLL).Fill(fabs(selectedEvent.phi1-selectedEvent.phi2),selectedEvent.totalWeight);
            m2->at(DeltaRLL).Fill(fabs(selectedEvent.deltaRll),selectedEvent.totalWeight);
            m2->at(DeltaEtaLLG).Fill(fabs((selectedEvent.l1+selectedEvent.l2).Eta()-selectedEvent.selPhotons.at(0).vec.Eta()), selectedEvent.totalWeight);
            m2->at(DeltaPhiLLG).Fill(fabs((selectedEvent.l1+selectedEvent.l2).Phi()-selectedEvent.selPhotons.at(0).vec.Phi()),selectedEvent.totalWeight);
            m2->at(DeltaRLLG).Fill(fabs((selectedEvent.l1+selectedEvent.l2).DeltaR(selectedEvent.selPhotons.at(0).vec)),selectedEvent.totalWeight);
         }
      }
   }
  
   if(SelectEvent(ONZ)){
      if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81 && selectedEvent.mll<101)){
         auto m1 = &h1Maps["onZEE"];
         auto m2 = &h1Maps["onZMM"];
         auto m3 = &h1Maps["onZ"];
         m3->at(ETMISS).Fill(selectedEvent.ETmiss, selectedEvent.totalWeight);
         m3->at(PT1).Fill(selectedEvent.pt1, selectedEvent.totalWeight);
         m3->at(PT2).Fill(selectedEvent.pt2, selectedEvent.totalWeight);
         m3->at(MLL).Fill(selectedEvent.mll,selectedEvent.totalWeight);   
         m3->at(NPHOTONS).Fill(selectedEvent.selPhotons.size(),selectedEvent.totalWeight);
         m3->at(NVTX).Fill(*nGoodVertices,selectedEvent.totalWeight);
         m3->at(HT).Fill(*ht,selectedEvent.totalWeight);
         m3->at(GENHT).Fill(*genHt,selectedEvent.totalWeight);
         m3->at(NJETS).Fill(selectedEvent.selJets.size(),selectedEvent.totalWeight);
         m3->at(ETA1).Fill(fabs(selectedEvent.eta1),selectedEvent.totalWeight);
         m3->at(ETA2).Fill(fabs(selectedEvent.eta2),selectedEvent.totalWeight);
         m3->at(PHI1).Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
         m3->at(PHI1).Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
         m3->at(PTG1).Fill(selectedEvent.selPhotons.at(0).p.Pt(),selectedEvent.totalWeight);
         m3->at(PHIG1).Fill(selectedEvent.selPhotons.at(0).p.Phi(),selectedEvent.totalWeight);
         m3->at(ETAG1).Fill(selectedEvent.selPhotons.at(0).p.Eta(),selectedEvent.totalWeight);
         m3->at(SIGMAIETAIETAG1).Fill(selectedEvent.selPhotons.at(0).sigmaIetaIeta,selectedEvent.totalWeight);
         m3->at(DeltaEtaLL).Fill(fabs(selectedEvent.eta1-selectedEvent.eta2),selectedEvent.totalWeight);
         m3->at(DeltaPhiLL).Fill(fabs(selectedEvent.phi1-selectedEvent.phi2),selectedEvent.totalWeight);
         m3->at(DeltaRLL).Fill(fabs(selectedEvent.deltaRll),selectedEvent.totalWeight);
         m3->at(DeltaEtaLLG).Fill(fabs((selectedEvent.l1+selectedEvent.l2).Eta()-selectedEvent.selPhotons.at(0).vec.Eta()), selectedEvent.totalWeight);
         m3->at(DeltaPhiLLG).Fill(fabs((selectedEvent.l1+selectedEvent.l2).Phi()-selectedEvent.selPhotons.at(0).vec.Phi()),selectedEvent.totalWeight);
         m3->at(DeltaRLLG).Fill(fabs((selectedEvent.l1+selectedEvent.l2).DeltaR(selectedEvent.selPhotons.at(0).vec)),selectedEvent.totalWeight);
         if (selectedEvent.isDiElectron){ 
            m1->at(ETMISS).Fill(selectedEvent.ETmiss, selectedEvent.totalWeight);
            m1->at(PT1).Fill(selectedEvent.pt1, selectedEvent.totalWeight);
            m1->at(PT2).Fill(selectedEvent.pt2, selectedEvent.totalWeight);
            m1->at(MLL).Fill(selectedEvent.mll,selectedEvent.totalWeight);   
            m1->at(NPHOTONS).Fill(selectedEvent.selPhotons.size(),selectedEvent.totalWeight);
            m1->at(NVTX).Fill(*nGoodVertices,selectedEvent.totalWeight);
            m1->at(HT).Fill(*ht,selectedEvent.totalWeight);
            m1->at(GENHT).Fill(*genHt,selectedEvent.totalWeight);
            m1->at(NJETS).Fill(selectedEvent.selJets.size(),selectedEvent.totalWeight);
            m1->at(ETA1).Fill(fabs(selectedEvent.eta1),selectedEvent.totalWeight);
            m1->at(ETA2).Fill(fabs(selectedEvent.eta2),selectedEvent.totalWeight);
            m1->at(PHI1).Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
            m1->at(PHI2).Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
            m1->at(PTG1).Fill(selectedEvent.selPhotons.at(0).p.Pt(),selectedEvent.totalWeight);
            m1->at(PHIG1).Fill(fabs(selectedEvent.selPhotons.at(0).p.Phi()),selectedEvent.totalWeight);
            m1->at(ETAG1).Fill(fabs(selectedEvent.selPhotons.at(0).p.Eta()),selectedEvent.totalWeight);
            m1->at(SIGMAIETAIETAG1).Fill(selectedEvent.selPhotons.at(0).sigmaIetaIeta,selectedEvent.totalWeight);
            m1->at(DeltaEtaLL).Fill(fabs(selectedEvent.eta1-selectedEvent.eta2),selectedEvent.totalWeight);
            m1->at(DeltaPhiLL).Fill(fabs(selectedEvent.phi1-selectedEvent.phi2),selectedEvent.totalWeight);
            m1->at(DeltaRLL).Fill(fabs(selectedEvent.deltaRll),selectedEvent.totalWeight);
            m1->at(DeltaEtaLLG).Fill(fabs((selectedEvent.l1+selectedEvent.l2).Eta()-selectedEvent.selPhotons.at(0).vec.Eta()), selectedEvent.totalWeight);
            m1->at(DeltaPhiLLG).Fill(fabs((selectedEvent.l1+selectedEvent.l2).Phi()-selectedEvent.selPhotons.at(0).vec.Phi()),selectedEvent.totalWeight);
            m1->at(DeltaRLLG).Fill(fabs((selectedEvent.l1+selectedEvent.l2).DeltaR(selectedEvent.selPhotons.at(0).vec)),selectedEvent.totalWeight);
         }else{
            m2->at(ETMISS).Fill(selectedEvent.ETmiss, selectedEvent.totalWeight);
            m2->at(PT1).Fill(selectedEvent.pt1, selectedEvent.totalWeight);
            m2->at(PT2).Fill(selectedEvent.pt2, selectedEvent.totalWeight);
            m2->at(MLL).Fill(selectedEvent.mll,selectedEvent.totalWeight);   
            m2->at(NPHOTONS).Fill(selectedEvent.selPhotons.size(),selectedEvent.totalWeight);
            m2->at(NVTX).Fill(*nGoodVertices,selectedEvent.totalWeight);
            m2->at(HT).Fill(*ht,selectedEvent.totalWeight);
            m2->at(GENHT).Fill(*genHt,selectedEvent.totalWeight);
            m2->at(NJETS).Fill(selectedEvent.selJets.size(),selectedEvent.totalWeight);
            m2->at(ETA1).Fill(fabs(selectedEvent.eta1),selectedEvent.totalWeight);
            m2->at(ETA2).Fill(fabs(selectedEvent.eta2),selectedEvent.totalWeight);
            m2->at(PHI1).Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
            m2->at(PHI2).Fill(fabs(selectedEvent.phi2),selectedEvent.totalWeight);
            m2->at(PTG1).Fill(selectedEvent.selPhotons.at(0).p.Pt(),selectedEvent.totalWeight);
            m2->at(PHIG1).Fill(fabs(selectedEvent.selPhotons.at(0).p.Phi()),selectedEvent.totalWeight);
            m2->at(ETAG1).Fill(fabs(selectedEvent.selPhotons.at(0).p.Eta()),selectedEvent.totalWeight);
            m2->at(SIGMAIETAIETAG1).Fill(selectedEvent.selPhotons.at(0).sigmaIetaIeta,selectedEvent.totalWeight);
            m2->at(DeltaEtaLL).Fill(fabs(selectedEvent.eta1-selectedEvent.eta2),selectedEvent.totalWeight);
            m2->at(DeltaPhiLL).Fill(fabs(selectedEvent.phi1-selectedEvent.phi2),selectedEvent.totalWeight);
            m2->at(DeltaRLL).Fill(fabs(selectedEvent.deltaRll),selectedEvent.totalWeight);
            m2->at(DeltaEtaLLG).Fill(fabs((selectedEvent.l1+selectedEvent.l2).Eta()-selectedEvent.selPhotons.at(0).vec.Eta()), selectedEvent.totalWeight);
            m2->at(DeltaPhiLLG).Fill(fabs((selectedEvent.l1+selectedEvent.l2).Phi()-selectedEvent.selPhotons.at(0).vec.Phi()),selectedEvent.totalWeight);
            m2->at(DeltaRLLG).Fill(fabs((selectedEvent.l1+selectedEvent.l2).DeltaR(selectedEvent.selPhotons.at(0).vec)),selectedEvent.totalWeight);
         }
      }
   }
}


void HistogramProducer::FillTriggerStudies(){
   if(SelectEventTriggerStudies(TRIGDILEP)){
      auto m1 = &eff1Maps["trigDilepEE"];
      auto m2 = &eff1Maps["trigDilepMM"];
      if (selectedEvent.trigHt){//baselineTrigger
         if (selectedEvent.isDiElectron){ 
           m1->at(ETMISS).Fill(selectedEvent.trigDiEle,selectedEvent.ETmiss);
           m1->at(PT1).Fill(selectedEvent.trigDiEle,selectedEvent.pt1);
           m1->at(PT2).Fill(selectedEvent.trigDiEle,selectedEvent.pt2);
           m1->at(MLL).Fill(selectedEvent.trigDiEle,selectedEvent.mll);   
           m1->at(NPHOTONS).Fill(selectedEvent.trigDiEle,selectedEvent.selPhotons.size());
           m1->at(NVTX).Fill(selectedEvent.trigDiEle,*nGoodVertices);
           m1->at(HT).Fill(selectedEvent.trigDiEle,*ht);
           m1->at(GENHT).Fill(selectedEvent.trigDiEle,*genHt);
           m1->at(NJETS).Fill(selectedEvent.trigDiEle,selectedEvent.selJets.size());
           m1->at(ETA1).Fill(selectedEvent.trigDiEle,fabs(selectedEvent.eta1));
           m1->at(ETA2).Fill(selectedEvent.trigDiEle,fabs(selectedEvent.eta2));
           m1->at(PHI1).Fill(selectedEvent.trigDiEle,fabs(selectedEvent.phi2));
           m1->at(PHI2).Fill(selectedEvent.trigDiEle,fabs(selectedEvent.phi2));
         }else{
           m2->at(ETMISS).Fill(selectedEvent.trigDiMu,selectedEvent.ETmiss);
           m2->at(PT1).Fill(selectedEvent.trigDiMu,selectedEvent.pt1);
           m2->at(PT2).Fill(selectedEvent.trigDiMu,selectedEvent.pt2);
           m2->at(MLL).Fill(selectedEvent.trigDiMu,selectedEvent.mll);   
           m2->at(NPHOTONS).Fill(selectedEvent.trigDiMu,selectedEvent.selPhotons.size());
           m2->at(NVTX).Fill(selectedEvent.trigDiMu,*nGoodVertices);
           m2->at(HT).Fill(selectedEvent.trigDiMu,*ht);
           m2->at(GENHT).Fill(selectedEvent.trigDiMu,*genHt);
           m2->at(NJETS).Fill(selectedEvent.trigDiMu,selectedEvent.selJets.size());
           m2->at(ETA1).Fill(selectedEvent.trigDiMu,fabs(selectedEvent.eta1));
           m2->at(ETA2).Fill(selectedEvent.trigDiMu,fabs(selectedEvent.eta2));
           m2->at(PHI1).Fill(selectedEvent.trigDiMu,fabs(selectedEvent.phi2));
           m2->at(PHI2).Fill(selectedEvent.trigDiMu,fabs(selectedEvent.phi2));      
         }
      }
   }
   if(SelectEvent(TRIGSEL)){
      if (selectedEvent.selPhotons.size()!=0){
         auto m1 = &eff1Maps["trigSelEE"];
         auto m2 = &eff1Maps["trigSelMM"];
         if (selectedEvent.trigHt){//baselineTrigger
            if (selectedEvent.isDiElectron){ 
               m1->at(ETMISS).Fill(selectedEvent.trigDiEle,selectedEvent.ETmiss);
               m1->at(PT1).Fill(selectedEvent.trigDiEle,selectedEvent.pt1);
               m1->at(PT2).Fill(selectedEvent.trigDiEle,selectedEvent.pt2);
               m1->at(MLL).Fill(selectedEvent.trigDiEle,selectedEvent.mll);   
               m1->at(NPHOTONS).Fill(selectedEvent.trigDiEle,selectedEvent.selPhotons.size());
               m1->at(NVTX).Fill(selectedEvent.trigDiEle,*nGoodVertices);
               m1->at(HT).Fill(selectedEvent.trigDiEle,*ht);
               m1->at(GENHT).Fill(selectedEvent.trigDiEle,*genHt);
               m1->at(NJETS).Fill(selectedEvent.trigDiEle,selectedEvent.selJets.size());
               m1->at(ETA1).Fill(selectedEvent.trigDiEle,fabs(selectedEvent.eta1));
               m1->at(ETA2).Fill(selectedEvent.trigDiEle,fabs(selectedEvent.eta2));
               m1->at(PHI1).Fill(selectedEvent.trigDiEle,fabs(selectedEvent.phi2));
               m1->at(PHI2).Fill(selectedEvent.trigDiEle,fabs(selectedEvent.phi2));
               m1->at(PTG1).Fill(selectedEvent.trigDiEle,selectedEvent.selPhotons.at(0).p.Pt());
               m1->at(PHIG1).Fill(selectedEvent.trigDiEle,fabs(selectedEvent.selPhotons.at(0).p.Phi()));
               m1->at(ETAG1).Fill(selectedEvent.trigDiEle,fabs(selectedEvent.selPhotons.at(0).p.Eta()));
               m1->at(SIGMAIETAIETAG1).Fill(selectedEvent.trigDiEle,selectedEvent.selPhotons.at(0).sigmaIetaIeta);
            }else{
               m2->at(ETMISS).Fill(selectedEvent.trigDiMu,selectedEvent.ETmiss);
               m2->at(PT1).Fill(selectedEvent.trigDiMu,selectedEvent.pt1);
               m2->at(PT2).Fill(selectedEvent.trigDiMu,selectedEvent.pt2);
               m2->at(MLL).Fill(selectedEvent.trigDiMu,selectedEvent.mll);   
               m2->at(NPHOTONS).Fill(selectedEvent.trigDiMu,selectedEvent.selPhotons.size());
               m2->at(NVTX).Fill(selectedEvent.trigDiMu,*nGoodVertices);
               m2->at(HT).Fill(selectedEvent.trigDiMu,*ht);
               m2->at(GENHT).Fill(selectedEvent.trigDiMu,*genHt);
               m2->at(NJETS).Fill(selectedEvent.trigDiMu,selectedEvent.selJets.size());
               m2->at(ETA1).Fill(selectedEvent.trigDiMu,fabs(selectedEvent.eta1));
               m2->at(ETA2).Fill(selectedEvent.trigDiMu,fabs(selectedEvent.eta2));
               m2->at(PHI1).Fill(selectedEvent.trigDiMu,fabs(selectedEvent.phi2));
               m2->at(PHI2).Fill(selectedEvent.trigDiMu,fabs(selectedEvent.phi2));
               m2->at(PTG1).Fill(selectedEvent.trigDiMu,selectedEvent.selPhotons.at(0).p.Pt());
               m2->at(PHIG1).Fill(selectedEvent.trigDiMu,fabs(selectedEvent.selPhotons.at(0).p.Phi()));
               m2->at(ETAG1).Fill(selectedEvent.trigDiMu,fabs(selectedEvent.selPhotons.at(0).p.Eta()));
               m2->at(SIGMAIETAIETAG1).Fill(selectedEvent.trigDiMu,selectedEvent.selPhotons.at(0).sigmaIetaIeta);
            }
         }
      }
   }
   if(SelectEvent(TRIGONZ)){
      if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81 && selectedEvent.mll<101)){
         auto m1 = &eff1Maps["trigOnZEE"];
         auto m2 = &eff1Maps["trigOnZMM"];
         if (selectedEvent.trigHt){//baselineTrigger
            if (selectedEvent.isDiElectron){ 
               m1->at(ETMISS).Fill(selectedEvent.trigDiEle,selectedEvent.ETmiss);
               m1->at(PT1).Fill(selectedEvent.trigDiEle,selectedEvent.pt1);
               m1->at(PT2).Fill(selectedEvent.trigDiEle,selectedEvent.pt2);
               m1->at(MLL).Fill(selectedEvent.trigDiEle,selectedEvent.mll);   
               m1->at(NPHOTONS).Fill(selectedEvent.trigDiEle,selectedEvent.selPhotons.size());
               m1->at(NVTX).Fill(selectedEvent.trigDiEle,*nGoodVertices);
               m1->at(HT).Fill(selectedEvent.trigDiEle,*ht);
               m1->at(GENHT).Fill(selectedEvent.trigDiEle,*genHt);
               m1->at(NJETS).Fill(selectedEvent.trigDiEle,selectedEvent.selJets.size());
               m1->at(ETA1).Fill(selectedEvent.trigDiEle,fabs(selectedEvent.eta1));
               m1->at(ETA2).Fill(selectedEvent.trigDiEle,fabs(selectedEvent.eta2));
               m1->at(PHI1).Fill(selectedEvent.trigDiEle,fabs(selectedEvent.phi2));
               m1->at(PHI2).Fill(selectedEvent.trigDiEle,fabs(selectedEvent.phi2));
               m1->at(PTG1).Fill(selectedEvent.trigDiEle,selectedEvent.selPhotons.at(0).p.Pt());
               m1->at(PHIG1).Fill(selectedEvent.trigDiEle,fabs(selectedEvent.selPhotons.at(0).p.Phi()));
               m1->at(ETAG1).Fill(selectedEvent.trigDiEle,fabs(selectedEvent.selPhotons.at(0).p.Eta()));
               m1->at(SIGMAIETAIETAG1).Fill(selectedEvent.trigDiEle,selectedEvent.selPhotons.at(0).sigmaIetaIeta);
            }else{
               m2->at(ETMISS).Fill(selectedEvent.trigDiMu,selectedEvent.ETmiss);
               m2->at(PT1).Fill(selectedEvent.trigDiMu,selectedEvent.pt1);
               m2->at(PT2).Fill(selectedEvent.trigDiMu,selectedEvent.pt2);
               m2->at(MLL).Fill(selectedEvent.trigDiMu,selectedEvent.mll);   
               m2->at(NPHOTONS).Fill(selectedEvent.trigDiMu,selectedEvent.selPhotons.size());
               m2->at(NVTX).Fill(selectedEvent.trigDiMu,*nGoodVertices);
               m2->at(HT).Fill(selectedEvent.trigDiMu,*ht);
               m2->at(GENHT).Fill(selectedEvent.trigDiMu,*genHt);
               m2->at(NJETS).Fill(selectedEvent.trigDiMu,selectedEvent.selJets.size());
               m2->at(ETA1).Fill(selectedEvent.trigDiMu,fabs(selectedEvent.eta1));
               m2->at(ETA2).Fill(selectedEvent.trigDiMu,fabs(selectedEvent.eta2));
               m2->at(PHI1).Fill(selectedEvent.trigDiMu,fabs(selectedEvent.phi2));
               m2->at(PHI2).Fill(selectedEvent.trigDiMu,fabs(selectedEvent.phi2));
               m2->at(PTG1).Fill(selectedEvent.trigDiMu,selectedEvent.selPhotons.at(0).p.Pt());
               m2->at(PHIG1).Fill(selectedEvent.trigDiMu,fabs(selectedEvent.selPhotons.at(0).p.Phi()));
               m2->at(ETAG1).Fill(selectedEvent.trigDiMu,fabs(selectedEvent.selPhotons.at(0).p.Eta()));
               m2->at(SIGMAIETAIETAG1).Fill(selectedEvent.trigDiMu,selectedEvent.selPhotons.at(0).sigmaIetaIeta);
            }
         }
      }
   }
}


Bool_t HistogramProducer::Process(Long64_t entry){
   //float tempPercentage = (float) entry/ (float)nEntries;
   //if (tempPercentage>0.05){
      //return kTRUE;
   //}
   fReader.SetLocalEntry(entry);
   totalWeight = *mc_weight * *pu_weight;
   
   //print out progression
   //double currentEventNo = (double) entry;
   //double totalEventNo = (double) nEntries;
   //double tempVerbose = currentEventNo/totalEventNo*100.;
   //if(entry%1000000==0){
      //int tempVerbose2 = (int)tempVerbose;
      //cout<<"output/"+getOutputFilename(inputName)<<" | "<<tempVerbose2<<"%"<<endl;
   //}

   FillHistograms();
   //FillHistograms2D();
   
   FillTriggerStudies();
   
   return kTRUE;
}



template<typename T>
void save2File(const map<string,map<Histograms1D,T>>& hMaps, TFile& file)
{
  for (auto& hMapIt : hMaps) {
    if (!file.Get(hMapIt.first.c_str())) {
      file.mkdir(hMapIt.first.c_str());
    }
    file.cd(hMapIt.first.c_str());
    for (auto& h : hMapIt.second) {
      h.second.Write(histoNames[h.first].c_str(), TObject::kWriteDelete);
    }
    file.cd();
  }
}
template<typename T>
void save2File(const map<string,map<Histograms2D,T>>& hMaps, TFile& file)
{
  for (auto& hMapIt : hMaps) {
    if (!file.Get(hMapIt.first.c_str())) {
      file.mkdir(hMapIt.first.c_str());
    }
    file.cd(hMapIt.first.c_str());
    for (auto& h : hMapIt.second) {
      h.second.Write(histoNames2D[h.first].c_str(), TObject::kWriteDelete);
    }
    file.cd();
  }
}




void HistogramProducer::Terminate()
{
   //auto outputName = "output_01/"+getOutputFilename(inputName);
   auto outputName = "output/"+getOutputFilename(inputName);
   TFile file(outputName.c_str(), "RECREATE");
   //save2File(h1Map, file);
   save2File(h1Maps, file);
   save2File(eff1Maps, file);
   //save2File(h2Maps, file);
   cutFlow.Write("hCutFlow");
   file.Close();
   cout << "Created " << outputName << " in " << (time(NULL) - startTime)/60 << " min" << endl;
   //cout << "Created " << getOutputFilename(inputName) << " in " << (time(NULL) - startTime)/60 << " min" << endl;
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
