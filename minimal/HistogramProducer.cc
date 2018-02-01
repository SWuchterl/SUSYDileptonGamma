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
   
   signal_m1(fReader,"signal_m1"),
   signal_m2(fReader,"signal_m2"),
   nBinos(fReader,"signal_nBinos"),
   startTime(time(NULL))//,
   //rand()
{
}



bool HistogramProducer::GenPhotonVeto(const int a){
   bool isGenPhoton=false;
   bool vetoPt130ZG=false;
   if(!isData){
      //for (vector<tree::GenParticle>::iterator it = genParticles->begin(); it != genParticles->end(); ++it){
      for (vector<tree::GenParticle>::iterator it = genParticles->begin(); it != genParticles->end(); it++){
         if(it->pdgId==22){
            
            if(isZGammaInclusive){
               if(it->p.Pt()>=130.){
                  vetoPt130ZG = vetoPt130ZG || true;
               }
            }
            
            switch (a)
            {
            case 0: 
               isGenPhoton=isGenPhoton || true; //veto all genPhotons
               break;
            
            case 1:
               isGenPhoton=isGenPhoton || (abs(it->motherId)<16 && abs(it->motherId)>10); //veto genPhotons from leptons (FSR)
               break;
            case 2:
               isGenPhoton=isGenPhoton || it->fromHardProcess || (abs(it->motherId)<16 && abs(it->motherId)>10);
               break;
            case 3:
               isGenPhoton=isGenPhoton || it->isPrompt ;
               break;
            case 4:
               isGenPhoton=isGenPhoton || !(it->promptStatus==DIRECTPROMPT) ;
               break;
            case 5:
               isGenPhoton=isGenPhoton || ((it->statusID==1) && (it->isPrompt)) ;
               break;
            case 6:
               isGenPhoton=isGenPhoton || ((it->isPromptFinalState)) ;
               break;
            default:
               isGenPhoton=false;
               break;
            }
            //isGenPhoton=isGenPhoton || true;
            //isGenPhoton=isGenPhoton || (abs(it->motherId)<7);
            //isGenPhoton = isGenPhoton || (abs(it->motherId)<7) || (abs(it->motherId)>16);
         }
      }
   }
   return (isGenPhoton && noPromptPhotons)||(vetoPt130ZG);
   //return false;
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
   bool isMuon;
   if (particle==M){
      isMuon=true;
   }else{
      isMuon=false;
   }
   isDiMuon=isMuon;
   float leptonMass = (isEle)?  0.51099895e-3 : 0.1056584;
   float qter = 1.0;
   
   int charge1 = (int) l1.charge;
   int charge2 = (int) l2.charge;
   //
   //pt1UnCor = (isEle)? electrons->at(0).pUncorrected.Pt() : l1.p.Pt();
   //pt2UnCor = (isEle)? electrons->at(1).pUncorrected.Pt() : l1.p.Pt();
   //eta1UnCor = (isEle)? electrons->at(0).pUncorrected.Eta() : l1.p.Eta();
   //eta2UnCor = (isEle)? electrons->at(1).pUncorrected.Eta() : l1.p.Eta();
   //phi1UnCor = (isEle)? electrons->at(0).pUncorrected.Phi() : l1.p.Phi();
   //phi2UnCor = (isEle)? electrons->at(1).pUncorrected.Phi() : l1.p.Phi();
   
   TLorentzVector lepton1 (0.,0.,0.,0.);
   TLorentzVector lepton2 (0.,0.,0.,0.);
   lepton1.SetPtEtaPhiM(l1.p.Pt(),l1.p.Eta(),l1.p.Phi(),leptonMass);
   lepton2.SetPtEtaPhiM(l2.p.Pt(),l2.p.Eta(),l2.p.Phi(),leptonMass);
   
   float pt1_old= l1.p.Pt();
   float pt2_old= l2.p.Pt();
   
   
   //if(isData && !isEle){
   if(isData && isMuon){
      rmcor.momcor_data(lepton1, charge1, 0, qter);
      rmcor.momcor_data(lepton2, charge2, 0, qter);
      miniIso1 = l1.miniIso*pt1_old/lepton1.Pt();
      miniIso2 = l2.miniIso*pt2_old/lepton2.Pt();
   }else{
      //if(!isEle && !isData){
      if(!isEle && isMuon){
         //rmcor.momcor_mc(lepton1, charge1, l1.nTrkLayers, qter);
         //rmcor.momcor_mc(lepton2, charge2, l2.nTrkLayers, qter);
         rmcor.momcor_mc(lepton1, charge1, muons->at(0).nTrkLayers, qter);
         rmcor.momcor_mc(lepton2, charge2, muons->at(0).nTrkLayers, qter);
      }
   }
   
   //pt1 = l1.p.Pt();  
   //pt2 = l2.p.Pt();  
   //phi1 = l1.p.Phi();  
   //phi2 = l2.p.Phi();  
   //eta1 = l1.p.Eta();  
   //eta2 = l2.p.Eta();
   pt1 = lepton1.Pt();  
   pt2 = lepton2.Pt();  
   phi1 = lepton1.Phi();  
   phi2 = lepton2.Phi();  
   eta1 = lepton1.Eta();  
   eta2 = lepton2.Eta();
   
   //for uncorrected plots only !!!!!!!!!!   //for uncorrected plots only !!!!!!!!!!   //for uncorrected plots only !!!!!!!!!!
   //for uncorrected plots only !!!!!!!!!!   //for uncorrected plots only !!!!!!!!!!   //for uncorrected plots only !!!!!!!!!!
   //for uncorrected plots only !!!!!!!!!!   //for uncorrected plots only !!!!!!!!!!   //for uncorrected plots only !!!!!!!!!!
   //pt1 = pt1UnCor;
   //pt2 = pt2UnCor;  
   //phi1 = phi1UnCor;  
   //phi2 = phi2UnCor;  
   //eta1 = eta1UnCor;  
   //eta2 = eta2UnCor;  
   //lepton1.SetPtEtaPhiM(pt1,eta1,phi1,leptonMass);
   //lepton2.SetPtEtaPhiM(pt2,eta2,phi2,leptonMass);

   
   //cout<<pt1_old<<" | "<<pt1<<endl;
   
   //int charge1 = (int) l1.charge;
   //int charge2 = (int) l2.charge;
   chargeProduct = charge1*charge2;
   //TLorentzVector lepton1 (0.,0.,0.,0.);
   //TLorentzVector lepton2 (0.,0.,0.,0.);
   //lepton1.SetPtEtaPhiM(pt1,eta1,phi1,leptonMass);
   //lepton2.SetPtEtaPhiM(pt2,eta2,phi2,leptonMass);
   lep1=lepton1;
   lep2=lepton2;
   mll = (lepton1+lepton2).M();
   //miniIso1 = l1.miniIso;
   //miniIso2 = l2.miniIso;
   //Delta R between all three partciles
   deltaRll=lepton1.DeltaR(lepton2);
   //MET
   ETmiss = mett.p.Pt();
   
   TLorentzVector missing (0.,0.,0.,0.);
   missing.SetPtEtaPhiM(mett.p.Pt(),mett.p.Eta(),mett.p.Phi(),0.);
   ETmiss_vec = missing;
   
   //double pa[3] = {aVec.M(),aVec.Px(),aVec.Py()};
   //double pb[3] = {bVec.M(),bVec.Px(),bVec.Py()};
   double pa[3] = {lep1.M(),lep1.Px(),lep1.Py()};
   double pb[3] = {lep2.M(),lep2.Px(),lep2.Py()};
  
   //double pmiss[3] = {0.,met.Px(),met.Py()};
   double pmiss[3] = {0.,missing.Px(),missing.Py()};
   fctMT2_.set_mn(0.);
   //fctMT2_.set_momenta(pa,pb,pmiss,*(longIntBranches_[treeName]["eventNr"]));
   fctMT2_.set_momenta(pa,pb,pmiss,*evtNo);
   //*(floatBranches_[treeName]["MT2"]) = static_cast<float>(fctMT2_.get_mt2()); 
   MT2_val = static_cast<float>(fctMT2_.get_mt2()); 
   
   evtHasGenPhotonVeto = GenPhotonVeto(6);
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
   }else{
      trigDiEle = true;
      trigDiMu = true;
      trigMuEle = false;
      trigHt = true;
   }
   

   if (isData){           
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
        
      if (isHTSelection && trigHt && Check2Ele()) {
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
      if (isHTSelection && trigHt && Check2Mu()) {
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


void HistogramProducer::clearCutFlowMap(){
   decisionMapCutFlowFine[TRIGGERED]=false;
   decisionMapCutFlowFine[LEPTONID_leading]=false;
   decisionMapCutFlowFine[LEPTONID_trailing]=false;
   decisionMapCutFlowFine[LEPTONPT_leading]=false;
   decisionMapCutFlowFine[LEPTONPT_trailing]=false;
   decisionMapCutFlowFine[PHOTON1]=false;
   decisionMapCutFlowFine[PHOTON1ID]=false;
   decisionMapCutFlowFine[PHOTON1PT]=false;
   decisionMapCutFlowFine[PHOTON1DR]=false;
   decisionMapCutFlowFine[M50]=false;
   decisionMapCutFlowFine[DIELECTRON]=false;
   decisionMapCutFlowFine[DIMUON]=false;
   decisionMapCutFlowFine[ZMASS]=false;
   decisionMapCutFlowFine[GENVETO]=false;
   
   decisionMapCutFlowFine_weight[TRIGGERED]=totalWeight;
   decisionMapCutFlowFine_weight[LEPTONID_leading]=totalWeight;
   decisionMapCutFlowFine_weight[LEPTONID_trailing]=totalWeight;
   decisionMapCutFlowFine_weight[LEPTONPT_leading]=totalWeight;
   decisionMapCutFlowFine_weight[LEPTONPT_trailing]=totalWeight;
   decisionMapCutFlowFine_weight[PHOTON1]=totalWeight;
   decisionMapCutFlowFine_weight[PHOTON1ID]=totalWeight;
   decisionMapCutFlowFine_weight[PHOTON1PT]=totalWeight;
   decisionMapCutFlowFine_weight[PHOTON1DR]=totalWeight;
   decisionMapCutFlowFine_weight[M50]=totalWeight;
   decisionMapCutFlowFine_weight[DIELECTRON]=totalWeight;
   decisionMapCutFlowFine_weight[DIMUON]=totalWeight;
   decisionMapCutFlowFine_weight[ZMASS]=totalWeight;
   decisionMapCutFlowFine_weight[GENVETO]=totalWeight;
}



//ELECTRON
bool HistogramProducer::testSelection(const tree::Electron& pa,selectionType selection ,bool leading){
   bool decision=false;
   if ((selection==UNCUT)||(selection==PHOTON)){
      decision = pa.isPassConvVeto && pa.passImpactParameter && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && pa.isTightMVA;
      //if (pa.isPassConvVeto && pa.passImpactParameter && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && pa.isTightMVA){
         //decisionMapCutFlowFine[LEPTON2ID]=true;
      //}
   }
   if(selection==SEL){
      decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVA;
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVASlope;
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isLoose;
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTight;
      //if (pa.isPassConvVeto && pa.passImpactParameter && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVA){
         //if(selection==ONZ) leading? decisionMapCutFlowFine[LEPTONID_leading]=true : decisionMapCutFlowFine[LEPTONID_trailing]=true;
         //if((leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.))){
            //if(selection==ONZ) leading? decisionMapCutFlowFine[LEPTONPT_leading]=true : decisionMapCutFlowFine[LEPTONPT_trailing]=true;
         //}
      //}
   
   }
   if(selection==ONZ){
      decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVA;
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVASlope;
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isLoose;
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTight;
      if (pa.isPassConvVeto && pa.passImpactParameter && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVA){
         leading? decisionMapCutFlowFine[LEPTONID_leading]=true : decisionMapCutFlowFine[LEPTONID_trailing]=true;
         if((leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.))){
            leading? decisionMapCutFlowFine[LEPTONPT_leading]=true : decisionMapCutFlowFine[LEPTONPT_trailing]=true;
         }
      }
   
   }
   if(selection==DILEP){
      decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVA ;
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVASlope ;
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isLoose;
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTight;
      //if(pa.isPassConvVeto && pa.passImpactParameter && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVA){
         //leading? decisionMapCutFlowFine[LEPTONID_leading]=true : decisionMapCutFlowFine[LEPTONID_trailing]=true;
         //if((leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.))){
            //leading? decisionMapCutFlowFine[LEPTONPT_leading]=true : decisionMapCutFlowFine[LEPTONPT_trailing]=true;
         //}
      //}
   }
   if(selection==TRIGSEL || selection==TRIGDILEP || selection==TRIGONZ){
      decision =(*ht>200.) && pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>0.) : (pa.p.Pt()>0.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVA;
   }
   if(selection==TRIGSEL_ptcuts || selection==TRIGDILEP_ptcuts || selection==TRIGONZ_ptcuts){
      decision =(*ht>200.) && pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVA;
   }
   if(selection==EXO){
      decision = pa.isPassConvVeto && pa.passImpactParameter && (leading? (fabs(pa.p.Eta())<1.4442 ? (pa.p.Pt()>65.) : (pa.p.Pt()>70.)) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.5) && (pa.miniIso<0.1) && pa.isTightMVA;
   }
   if(selection==EGRegression){
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading? (pa.p.Pt()>32.) : (pa.p.Pt()>25.)) && (fabs(pa.p.Eta())<2.5) && (pa.miniIso<0.1) && pa.isLoose;
      decision = (leading? (pa.p.Pt()>32.) : (pa.p.Pt()>25.)) && (fabs(pa.p.Eta())<2.5) && pa.isLoose;
   }
   return decision;
}
//MUON
bool HistogramProducer::testSelection(const tree::Muon& pa, selectionType selection, bool leading){
   bool decision=false;
   if ((selection==UNCUT)||(selection==PHOTON)){
      decision = pa.passImpactParameter  && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (pa.isTight);
   }
   if (selection==SEL){
      decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (deltaRll>0.1);
      //if(pa.passImpactParameter && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (deltaRll>0.1)){
         //if(selection==ONZ) leading? decisionMapCutFlowFine[LEPTONID_leading]=true : decisionMapCutFlowFine[LEPTONID_trailing]=true;
         //if((leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.))){
            //if(selection==ONZ) leading? decisionMapCutFlowFine[LEPTONPT_leading]=true : decisionMapCutFlowFine[LEPTONPT_trailing]=true;
         //}
      //}
   }
   if (selection == ONZ){
      decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (deltaRll>0.1);
      if(pa.passImpactParameter && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (deltaRll>0.1)){
         leading? decisionMapCutFlowFine[LEPTONID_leading]=true : decisionMapCutFlowFine[LEPTONID_trailing]=true;
         if((leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.))){
            leading? decisionMapCutFlowFine[LEPTONPT_leading]=true : decisionMapCutFlowFine[LEPTONPT_trailing]=true;
         }
      }
   }
   if(selection==DILEP || selection==EGRegression){
      decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (deltaRll>0.1);
      //if (pa.passImpactParameter && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (deltaRll>0.1)){
         //leading? decisionMapCutFlowFine[LEPTONID_leading]=true : decisionMapCutFlowFine[LEPTONID_trailing]=true;
         //if((leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.))){
            //leading? decisionMapCutFlowFine[LEPTONPT_leading]=true : decisionMapCutFlowFine[LEPTONPT_trailing]=true;
         //}
      //}
   }
   if(selection==TRIGSEL || selection==TRIGDILEP || selection==TRIGONZ){
      decision = (*ht>200.) && pa.passImpactParameter && (leading ? (pa.p.Pt()>0.) : (pa.p.Pt()>0.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (deltaRll>0.1);
   }
   if(selection==TRIGSEL_ptcuts || selection==TRIGDILEP_ptcuts || selection==TRIGONZ_ptcuts){
      decision = (*ht>200.)&& pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (deltaRll>0.1);
   }
   if(selection==TRIGDILEP_pt1cut){
      decision = (*ht>200.)&& pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>0.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (deltaRll>0.1);
   }
   if(selection==TRIGDILEP_pt2cut){
      decision = (*ht>200.)&& pa.passImpactParameter && (leading ? (pa.p.Pt()>0.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (deltaRll>0.1);
   }
   if(selection==EXO){
      decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>52.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (pa.isMedium);
   }
   return decision;
}
//PHOTON
bool HistogramProducer::testSelection(const selPhoton& pa, selectionType selection){
   bool decision=false;
   if ((selection==UNCUT)||(selection==PHOTON)){
      //decision = (pa.passElectronVeto) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isMediumMVA) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3);
      //decision = (pa.passElectronVeto) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isMediumMVA) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3);
      //decision = !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isMediumMVA) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3);
      decision = !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isLoose) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3);
   }
   if(selection==SEL || selection==DILEP ||selection==ONZ || selection==EGRegression){
      //decision = (pa.passElectronVeto) && (pa.p.Pt()>20.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isMediumMVA) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3); //study Delta R cut ;
      //decision = (pa.passElectronVeto) && (pa.p.Pt()>25.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isMediumMVA) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3); //study Delta R cut ;
      //decision = (pa.p.Pt()>25.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isMediumMVA) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3); //study Delta R cut ;
      decision = (pa.p.Pt()>25.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isLoose) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3); //study Delta R cut ;
      if(!(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isLoose)){
         if(selection==ONZ) decisionMapCutFlowFine[PHOTON1ID]=true;
         if((pa.p.Pt()>25.)){
            if(selection==ONZ) decisionMapCutFlowFine[PHOTON1PT]=true;
            if((pa.deltaR1>0.3) && (pa.deltaR2>0.3)){
               if(selection==ONZ) decisionMapCutFlowFine[PHOTON1DR]=true;
            }
         }
      }
   }
   if(selection==EXO){
      //decision = (isDiElectron? (pa.p.Pt()>(65. + 5.* ((1.566 < fabs(pa.p.Eta())) && (fabs(pa.p.Eta()) < 2.5))  )) : (pa.p.Pt()>40.)) && (pa.passElectronVeto) && (fabs(pa.p.Eta())<2.5) && (pa.isMediumMVA) && (pa.deltaR1>0.4) && (pa.deltaR2>0.4); //study Delta R cut ;
      decision = (isDiElectron? (pa.p.Pt()>(65. + 5.* ((1.566 < fabs(pa.p.Eta())) && (fabs(pa.p.Eta()) < 2.5))  )) : (pa.p.Pt()>40.)) && (pa.passElectronVeto) && (fabs(pa.p.Eta())<2.5) && (pa.isLoose) && (pa.deltaR1>0.4) && (pa.deltaR2>0.4); //study Delta R cut ;
   }
   if(selection==TRIGDILEP || selection==TRIGSEL || selection==TRIGONZ){
      //decision = (pa.passElectronVeto) && (pa.p.Pt()>20.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isMediumMVA) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3); //study Delta R cut ;
      //decision = (pa.passElectronVeto) && (pa.p.Pt()>25.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isMediumMVA) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3); //study Delta R cut ;
      //decision = (pa.p.Pt()>25.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isMediumMVA) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3); //study Delta R cut ;
      decision = (pa.p.Pt()>25.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isLoose) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3); //study Delta R cut ;
   }
   return decision;
}
//JET
bool HistogramProducer::testSelection(const selJet& pa, selectionType selection){
   bool decision=false;
   if ((selection==UNCUT)||(selection==PHOTON)){
      decision = (fabs(pa.p.Eta())<2.4) && (pa.isLoose) && !(pa.hasElectronMatch) && !(pa.hasPhotonMatch) && !(pa.hasMuonMatch);
   }
   if(selection==SEL||selection==DILEP || selection==ONZ || selection==TRIGDILEP || selection==TRIGSEL || selection==TRIGONZ || selection==EXO || selection==EGRegression){
      decision = (pa.p.Pt()>30.) && (fabs(pa.p.Eta())<2.4) && (pa.isLoose) && !(pa.hasElectronMatch) && !(pa.hasPhotonMatch) && !(pa.hasMuonMatch); //&& (pa.deltaR1>0.3) && (pa.deltaR2>0.3) //study Delta R cut ;
   }
   return decision;
}



void HistogramProducer::Init(TTree *tree)
{
   fReader.SetTree(tree);
   inputName = fReader.GetTree()->GetCurrentFile()->GetName();
   
   isTotalSignal = (inputName.find("SMS-TChiNG_BF") != string::npos) || (inputName.find("SMS-T5bbbbZg_nTuple")!=string::npos);
   isData = inputName.find("Run201") != string::npos;
   isSignal = inputName.find("SMS") != string::npos;
   
   boost::property_tree::ini_parser::read_ini("example.ini", propertyTree);

   config_selectionsToProcessMap[UNCUT] = propertyTree.get<bool>("selectionsToProcess.uncut");
   config_selectionsToProcessMap[PHOTON] = propertyTree.get<bool>("selectionsToProcess.photon");
   config_selectionsToProcessMap[ONZ] = propertyTree.get<bool>("selectionsToProcess.onz");
   config_selectionsToProcessMap[DILEP] = propertyTree.get<bool>("selectionsToProcess.dilep");
   config_selectionsToProcessMap[ONZMET] = propertyTree.get<bool>("selectionsToProcess.onzmet");
   config_selectionsToProcessMap[EGRegression] = propertyTree.get<bool>("selectionsToProcess.egregression");
   config_selectionsToProcessMap[EXO] = propertyTree.get<bool>("selectionsToProcess.exo");
   config_selectionsToProcessMap[SEL] = propertyTree.get<bool>("selectionsToProcess.sel");
   config_selectionsToProcessMap[ABOVEZG] = propertyTree.get<bool>("selectionsToProcess.abovezg");
   config_selectionsToProcessMap[ONZG] = propertyTree.get<bool>("selectionsToProcess.onzg");
   
   config_veto = propertyTree.get<int>("generalSettings.veto");
   config_docutflow = propertyTree.get<bool>("generalSettings.docutflow");
   config_docutflowfine = propertyTree.get<bool>("generalSettings.docutflowfine");
   config_dosignalscan = propertyTree.get<bool>("generalSettings.dosignalscan");
   config_eventpercentage = propertyTree.get<float>("generalSettings.eventpercentage");
   config_outputfolder = propertyTree.get<string>("generalSettings.outputfolder");
   
   //cout<<"Loaded example.ini config file succesfully."<<endl;
   
   if(! isTotalSignal){
      cutFlow = *((TH1F*)fReader.GetTree()->GetCurrentFile()->Get("TreeWriter/hCutFlow"));
   }
   fReader.GetEntries(true); // jumps to last file
   //string lastInputName = fReader.GetTree()->GetCurrentFile()->GetName();
  
   nEntries=fReader.GetTree()->GetEntries();
   //isData = inputName.find("Run201") != string::npos;
   //isSignal = inputName.find("SMS") != string::npos;
  
  
  //rochester muon pt corrections
   rochcor2016 *rmcor = new rochcor2016();
  
  
  setHistoNames();
  
   //if (inputName!=lastInputName) {
       ////////adds cut flow of last file. This makes only sense if there is for 1 or two files
      //cutFlow.Add((TH1F*)fReader.GetTree()->GetCurrentFile()->Get("TreeWriter/hCutFlow"));
   //}
   //float nGen = cutFlow.GetBinContent(2);


   noPromptPhotons = inputName.find("DYJets") != string::npos || inputName.find("TTTo2L2Nu") != string::npos;
   isZGammaInclusive = inputName.find("ZGTo2LG_ext") != string::npos;

   if (config_docutflow) InitCutFlowHistos();
   if (config_docutflowfine) InitCutFlowHistos_Fine();
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
   
   //PhotonIDWeighter = Weighter("scaleFactors/ScaleFactorsPhotonsID.root","EGamma_SF2D"); // pt vs eta
   PhotonIDWeighter = Weighter("scaleFactors/ScaleFactorsPhotonsID_LooseCut.root","EGamma_SF2D"); // pt vs eta
   PhotonIDWeighter.fillOverflow2d();
   
   //electronMllWeighter = Weighter("electronMllWeights.root","769ec98b61da6715"); // vs mll
   //electronMllWeighter.fillOverflow2d();
   
   }
   
void HistogramProducer::InitScaleFactorsAlternative(){
   //DiEleWeighterID=Weighter("scaleFactors/Alternative/ScaleFactorElectronID.root","EGamma_SF2D");// pt vs |eta|
   DiEleWeighterID=Weighter("scaleFactors/Alternative/ScaleFactorElectronID_MVATight.root","EGamma_SF2D");// pt vs |eta|
   //DiEleWeighterID=Weighter("scaleFactors/Alternative/ScaleFactorElectronID_LooseCut.root","EGamma_SF2D");// pt vs |eta|
   //DiEleWeighterID=Weighter("scaleFactors/Alternative/ScaleFactorElectronID_MediumCut.root","EGamma_SF2D");// pt vs |eta|
   //DiEleWeighterID=Weighter("scaleFactors/Alternative/ScaleFactorElectronID_TightCut.root","EGamma_SF2D");// pt vs |eta|
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
   
   //PhotonIDWeighter = Weighter("scaleFactors/ScaleFactorsPhotonsID.root","EGamma_SF2D"); // pt vs eta
   PhotonIDWeighter = Weighter("scaleFactors/ScaleFactorsPhotonsID_LooseCut.root","EGamma_SF2D"); // pt vs eta
   PhotonIDWeighter.fillOverflow2d();
   
   //electronMllWeighter = Weighter("electronMllWeights.root","769ec98b61da6715"); // vs mll
   //electronMllWeighter.fillOverflow2d();
   
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


float HistogramProducer::GetScaleFactorAndErrorPhotons(vector<selPhoton>& vecGamma){
   float totalSF = 1.;
   //for (vector<selPhoton>::iterator it = vecGamma.begin(); it != vecGamma.end(); ++it){
   for (vector<selPhoton>::iterator it = vecGamma.begin(); it != vecGamma.end(); it++){
      float tempPt=it->p.Pt();
      float tempEta=it->p.Eta();
      //0.9938 +- 0.0119 https://twiki.cern.ch/twiki/pub/CMS/EgammaIDRecipesRun2/EleVetoScalingFactors_Moriond17.pdf
      totalSF = totalSF * PhotonIDWeighter.getWeight(fabs(tempEta),tempPt)*0.9938;
   } 
   return totalSF;

}


map<Histograms1D,TH1F> HistogramProducer::InitHistograms(const selectionType selection_){
   map<Histograms1D,TH1F> hMap;

   hMap[ETMISS] = TH1F("", ";#it{p}_{T}^{miss} (GeV)", 5000, 0, 5000);
   hMap[PT1] = TH1F("", ";#it{p}_{T}^{leading} (GeV)", 5000, 0, 5000);
   hMap[PT2] = TH1F("", ";#it{p}_{T}^{trailing} (GeV)", 5000, 0, 5000);
   hMap[MLL] = TH1F("", ";#it{m}_{ll} (GeV)", 5000, 0, 5000);
   hMap[NPHOTONS] = TH1F("","n_#gamma",10,0,10);
   hMap[NVTX] = TH1F("","n_{Vtx}",60,0,60);
   hMap[HT] = TH1F("","#it{H}_{T} (GeV)",5000,0,5000);
   hMap[GENHT] = TH1F("","#it{H}_{T}^{gen} (GeV)",5000,0,5000);
   hMap[NJETS] = TH1F("","n_{Jets}",20,0,20);
   hMap[ETA1] = TH1F("", ";|#eta_{trailing}|", 260, 0, 2.6);
   hMap[ETA2] = TH1F("", ";|#eta_{leading}|", 260, 0, 2.6);
   hMap[PHI1] = TH1F("", ";|#phi_{trailing}|", 350, 0, 3.5);
   hMap[PHI2] = TH1F("", ";|#phi_{leading}|", 350, 0, 3.5);
   hMap[DeltaEtaLL] = TH1F("", ";#Delta#Eta_{ll}", 6000, 0, 6.);
   hMap[DeltaPhiLL] = TH1F("", ";#Delta#Phi_{ll}", 6000, 0, 6.);
   hMap[DeltaRLL] = TH1F("", ";#DeltaR_{ll}", 5000, 0, 6.);
   hMap[ZPT] = TH1F("", ";Z_{p_T}", 5000, 0, 5000);
   hMap[MTLL] = TH1F("", ";m_{T}^{ll}", 5000, 0, 5000);
   hMap[ST] = TH1F("", ";S_T", 5000, 0, 5000.);
   hMap[MT2] = TH1F("", ";M_{T2}", 50000, 0, 5000.);
   hMap[DeltaPhiLLMet] = TH1F("", "#Delta#Phi_{ll,MET}", 6000, 0, 6.);
   hMap[DeltaEtaLLMet] = TH1F("", ";#Delta#Eta_{ll,MET}", 6000, 0, 6.);
   hMap[DeltaRLLMet] = TH1F("", "#DeltaR_{ll,MET}", 6000, 0, 6.);
   hMap[VetoCompare] = TH1F("", "", 2, 0, 2);

   
   if ((selection_==PHOTON)||(selection_==SEL)||(selection_==ONZ)||(selection_==EXO)){
      hMap[PTG1] = TH1F("", ";#it{p}_{T}^{#gamma 1} (GeV)", 5000, 0, 5000);
      hMap[ETAG1] = TH1F("", ";|#eta_{#gamma 1}|", 260, 0, 2.6);
      hMap[PHIG1] = TH1F("", ";|#phi_{#gamma 1}|", 350, 0, 3.5);
      hMap[SIGMAIETAIETAG1] = TH1F("", ";#sigma_{i#etai#eta}^{#gamma 1}", 400, 0, 0.04);
      hMap[SIGMAIPHIIPHIG1] = TH1F("", ";#sigma_{i#phii#phi}^{#gamma 1}", 2000, 0, 0.2);
      hMap[R9] = TH1F("", ";r9", 1500, 0, 1.5);
      hMap[HOVERE] = TH1F("", ";H/E", 1000, 0, 0.1);
      hMap[DELTARGL1] = TH1F("", ";#DeltaR_{l1,#gamma}", 6000, 0, 6.);
      hMap[DELTARGL2] = TH1F("", ";#DeltaR_{l2,#gamma}", 6000, 0, 6.);
      hMap[DeltaRLLG] = TH1F("", ";#DeltaR_{ll,#gamma}", 6000, 0, 6.);
      hMap[DeltaEtaLLG] = TH1F("", ";#Delta#Eta_{ll,#gamma}", 600, 0, 6.);
      hMap[DeltaPhiLLG] = TH1F("", ";#Delta#Phi_{ll,#gamma}", 600, 0, 6.);
      hMap[STG] = TH1F("", ";S_T", 5000, 0, 5000.);
      hMap[STMET] = TH1F("", ";S_T + #it{p}_{T}^{miss} (GeV)", 5000, 0, 5000.);   
      hMap[MTLLG] = TH1F("", ";m_{T}^{ll#gamma}", 5000, 0, 5000);
      hMap[MTL1MET] = TH1F("", ";m_{T}^{l1,met}", 5000, 0, 5000);
      hMap[MTL2MET] = TH1F("", ";m_{T}^{l2,met}", 5000, 0, 5000);
      hMap[MTGMET] = TH1F("", ";m_{T}^{#gamma,met}", 5000, 0, 5000);
      hMap[MTLLMET] = TH1F("", ";m_{T}^{ll,met}", 5000, 0, 5000);
      hMap[MTLLGMET] = TH1F("", ";m_{T}^{ll#gamma,met}", 5000, 0, 5000);
      hMap[MOTHERID] = TH1F("", ";ID_{mother}", 200, 0, 200);
      hMap[MLLG] = TH1F("", ";ID_{ll#gamma}", 5000, 0, 5000);
      hMap[PT_llg] = TH1F("", ";p_{T}^{ll#gamma}", 5000, 0, 5000);
      hMap[MZG_exo] = TH1F("", ";m_{Z#gamma}", 5000, 0, 5000);
      hMap[gammaMotherID] = TH1F("", ";motherID_{#gamma}", 5000, 0, 5000);
      hMap[genPhotonPT] = TH1F("", ";gen p_T^{#gamma,matched}", 5000, 0, 5000);
      hMap[genPhotonPT_Veto] = TH1F("", ";gen p_T^{#gamma,matched,veto}", 5000, 0, 5000);
      hMap[PTG1_Veto] = TH1F("", ";gen p_T^{#gamma,veto}", 5000, 0, 5000);
      hMap[genPhotonPT_NoVeto] = TH1F("", ";gen p_T^{#gamma,matched,Noveto}", 5000, 0, 5000);
      hMap[PTG1_NoVeto] = TH1F("", ";gen p_T^{#gamma,Noveto}", 5000, 0, 5000);
      //hMap[VetoCompare] = TH1F("", "", 2, 0, 2);
   }

   return hMap;

}

map<Histograms1D,TH1F> HistogramProducer::InitSignalScanHistograms(const selectionType selection_){
   map<Histograms1D,TH1F> sMap;
   
   sMap[ETMISS] = TH1F("", ";#it{p}_{T}^{miss} (GeV)", 5000, 0, 5000);
   
   return sMap;
}


map<Histograms2D,TH2F> HistogramProducer::Init2DHistograms(const selectionType selection_){
    
   map<Histograms2D,TH2F> h2Map;
   
   if ((selection_==PHOTON)||(selection_==SEL)||(selection_==ONZ)||(selection_==EXO) ){
      h2Map[PTGvsMLLG] = TH2F("",";#it{m}_{ll#gamma};#it{p}_{T}^{gamma}",1000,0,1000,1000,0,1000);   
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
   hMap[ETA1] = TEfficiency("", ";|#eta_{leading}|", 260, 0, 2.6);
   hMap[ETA2] = TEfficiency("", ";|#eta_{trailing}|", 260, 0, 2.6);
   hMap[PHI1] = TEfficiency("", ";|#phi_{leading}|", 350, 0, 3.5);
   hMap[PHI2] = TEfficiency("", ";|#phi_{trailing}|", 350, 0, 3.5);

   if ((selection_==TRIGONZ)||(selection_==TRIGSEL)||(selection_==TRIGSEL_ptcuts)||(selection_==TRIGSEL_ptcuts)){
      hMap[PTG1] = TEfficiency("", ";#it{p}_{T}^{#gamma 1} (GeV)", 200, 0, 1000);
      hMap[ETAG1] = TEfficiency("", ";|#eta_{#gamma 1}|", 260, 0, 2.6);
      hMap[PHIG1] = TEfficiency("", ";|#phi_{#gamma 1}|", 350, 0, 3.5);
      hMap[SIGMAIETAIETAG1] = TEfficiency("", ";#sigma_{i#etai#eta}^{#gamma 1}", 400, 0, 0.04);
   }

   return hMap;

}


map<Histograms1D,TH1F> HistogramProducer::InitCutFlowHistograms(const selectionType selection){
   map<Histograms1D,TH1F> cMap;
   cMap[CUTFLOW] = TH1F("","",5,0,5);
   cMap[CUTFLOW].Fill("triggered", 0);
   cMap[CUTFLOW].Fill("2leptons", 0);
   cMap[CUTFLOW].Fill("m50", 0);
   cMap[CUTFLOW].Fill("1photon", 0);
   cMap[CUTFLOW].Fill("Z", 0);
   return cMap;
}


map<Histograms1D,TH1F> HistogramProducer::InitCutFlowHistograms_Fine(const selectionType selection){
   map<Histograms1D,TH1F> cMap;
   cMap[CUTFLOW_fine] = TH1F("","",10,0,10);
   cMap[CUTFLOW_fine].Fill("triggered", 0);
   cMap[CUTFLOW_fine].Fill("GenPhotonVeto", 0);
   cMap[CUTFLOW_fine].Fill("2LeptonID", 0);
   cMap[CUTFLOW_fine].Fill("2LeptonPT", 0);
   cMap[CUTFLOW_fine].Fill("m50", 0);
   cMap[CUTFLOW_fine].Fill("1Photon", 0);
   cMap[CUTFLOW_fine].Fill("1PhotonID", 0);
   cMap[CUTFLOW_fine].Fill("1PhotonPT", 0);
   cMap[CUTFLOW_fine].Fill("1PhotonDeltaR", 0);
   cMap[CUTFLOW_fine].Fill("Z", 0);
   return cMap;
}

void HistogramProducer::SlaveBegin(TTree *tree)
{
}

bool HistogramProducer::SelectEvent(selectionType selection){
   clearCutFlowMap();
   if (CheckParticles()){//only Fill Histograms if nMu>=2, nEle>=2, nGamma>=0
      if (Cleaning()){//decide in events with DiEle and DiMu where to contribute
         selEvent L;
         L.isDiElectron=isDiElectron;
         L.isDiMuon=isDiMuon;
         L.trigDiEle=trigDiEle;
         L.trigDiMu=trigDiMu;
         L.trigMuEle=trigMuEle;
         L.trigHt=trigHt;
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
         L.ETmiss_vec=ETmiss_vec;
         L.MT2_val = MT2_val;
         L.evtHasGenPhotonVeto = evtHasGenPhotonVeto;
         if(!isData){
            //L.totalWeight = totalWeight * GetScaleFactorAndError(pt1,eta1,isSignal,isDiElectron)*GetScaleFactorAndError(pt2,eta2,isSignal,isDiElectron);
            L.totalWeight = totalWeight * GetScaleFactorAndErrorAlternative(pt1,eta1,isSignal,isDiElectron,*runNo)*GetScaleFactorAndErrorAlternative(pt2,eta2,isSignal,isDiElectron,*runNo);
            //L.totalWeight = totalWeight * GetScaleFactorAndErrorAlternative(pt1UnCor,eta1UnCor,isSignal,isDiElectron,*runNo)*GetScaleFactorAndErrorAlternative(pt2UnCor,eta2UnCor,isSignal,isDiElectron,*runNo);
            //if ((selection==SEL||selection==ONZ||selection==DILEP)&&(isDiElectron)) totalWeight = totalWeight * electronMllWeighter.getWeight(mll);
            //if ((selection==SEL||ONZ||DILEP)&&(isDiElectron)) L.totalWeight = L.totalWeight * electronMllWeighter.getWeight(mll);
         }else{
            L.totalWeight=totalWeight;
         }
         //fill cutflowInfo DiEle or DiMu?
         cutflowDiEle=isDiElectron;
         cutflowDiMu=isDiMuon;
         decisionMapCutFlowFine[DIMUON]=isDiMuon;
         decisionMapCutFlowFine[DIELECTRON]=isDiElectron;
         
         if (L.isDiElectron){
            
            //Fill cutflow triggered bool
            if(selection==ONZ) cutflowIsTriggered=true;
            //if(selection==ONZ) decisionMapCutFlowFine[TRIGGERED]=true;
            if(selection==ONZ) decisionMapCutFlowFine[TRIGGERED]=trigDiEle;
            if(selection==ONZ) decisionMapCutFlowFine_weight[TRIGGERED]=L.totalWeight;
            
            auto e1 = electrons->at(0);
            auto e2 = electrons->at(1);
            
            //if(!GenPhotonVeto(5)){
            if(!GenPhotonVeto(config_veto)){
            //if(!GenPhotonVeto(config_veto+10)){
               if(selection==ONZ) decisionMapCutFlowFine[GENVETO]=true;
               if(selection==ONZ) decisionMapCutFlowFine_weight[GENVETO]=L.totalWeight;
            
               if(testSelection(e1,selection,true) && testSelection(e2,selection,false)){
                  //Fill cutflow 2 leptons bool
                  if(selection==ONZ) cutflow2Leptons=true;
                  if((selection==UNCUT)? true : (mll>50.)){
                     //Fill cutflow mll>50 bool
                     if(selection==ONZ) cutflowMll50=true;
                     if(selection==ONZ) decisionMapCutFlowFine[M50]=true;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONID_leading]=L.totalWeight;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONID_trailing]=L.totalWeight;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONPT_leading]=L.totalWeight;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONPT_trailing]=L.totalWeight;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[M50]=L.totalWeight;

                     //for (vector<tree::Photon>::iterator it = photons->begin(); it != photons->end(); ++it){
                     for (vector<tree::Photon>::iterator it = photons->begin(); it != photons->end(); it++){
                        auto g = *it;
                        selPhoton gamma;
                        gamma.setAll(g);
                        gamma.vec.SetPtEtaPhiM(gamma.p.Pt(),gamma.p.Eta(),gamma.p.Phi(),0.);
                        gamma.deltaR1=gamma.vec.DeltaR(L.l1);
                        gamma.deltaR2=gamma.vec.DeltaR(L.l2);
                        if (testSelection(gamma,selection)){
                           L.selPhotons.push_back(gamma);
                        }
                        //if(!isData){
                           //L.totalWeight=L.totalWeight*GetScaleFactorAndErrorPhotons(L.selPhotons);
                        //}
                     }
                     if(selection==ONZ) decisionMapCutFlowFine[PHOTON1]=true;
                     if(selection==ONZ) decisionMapCutFlowFine[PHOTON1]=L.totalWeight;
                     if(!isData){
                           L.totalWeight=L.totalWeight*GetScaleFactorAndErrorPhotons(L.selPhotons);
                     }
                     if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1ID]=L.totalWeight;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1PT]=L.totalWeight;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1DR]=L.totalWeight;
                     
                     //for (vector<tree::Jet>::iterator it = jets->begin(); it != jets->end(); ++it){
                     for (vector<tree::Jet>::iterator it = jets->begin(); it != jets->end(); it++){
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
                        if(selectedEvent.selPhotons.size()!=0){
                           if(selection==ONZ) cutflow1Photon=true;
                           //if(selection==ONZ) decisionMapCutFlowFine[PHOTON1]=true;
                           if(L.mll<101. && L.mll>81.){
                              if(selection==ONZ) cutflowOnZ=true;
                              if(selection==ONZ) decisionMapCutFlowFine[ZMASS]=true;
                              if(selection==ONZ) decisionMapCutFlowFine_weight[ZMASS]=L.totalWeight;
                           }
                        }
                     return true;
                     //else{
                        //return false;
                     //}
                  }
               }
            }
         }
         else{
            
            if(selection==ONZ) cutflowIsTriggered=true;
            //if(selection==ONZ) decisionMapCutFlowFine[TRIGGERED]=true;
            if(selection==ONZ) decisionMapCutFlowFine[TRIGGERED]=trigDiMu;

            auto m1 = muons->at(0); 
            auto m2 = muons->at(1);
            m1.p.SetPtEtaPhi(pt1,eta1,phi1);
            m2.p.SetPtEtaPhi(pt2,eta2,phi2);
            m1.miniIso=miniIso1;
            m2.miniIso=miniIso2;
            
            
            //if(!GenPhotonVeto(5)){
            if(!GenPhotonVeto(config_veto)){
            //if(!GenPhotonVeto(config_veto+10)){
               
               if(selection==ONZ) decisionMapCutFlowFine[GENVETO]=true;

            
               if(testSelection(m1,selection,true) && testSelection(m2,selection,false)){
                  if(selection==ONZ) cutflow2Leptons=true;
                  if((selection==UNCUT)? true : (mll>50.)){
                     if(selection==ONZ) cutflowMll50=true;
                     decisionMapCutFlowFine[M50]=true;
                     //for (vector<tree::Photon>::iterator it = photons->begin(); it != photons->end(); ++it){
                     for (vector<tree::Photon>::iterator it = photons->begin(); it != photons->end(); it++){
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
                     if(selection==ONZ) decisionMapCutFlowFine[PHOTON1]=true;
                     if(!isData){
                           L.totalWeight=L.totalWeight*GetScaleFactorAndErrorPhotons(L.selPhotons);
                     }
                     //for (vector<tree::Jet>::iterator it = jets->begin(); it != jets->end(); ++it){
                     for (vector<tree::Jet>::iterator it = jets->begin(); it != jets->end(); it++){
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
                        
                     if(L.selPhotons.size()>0){
                        if(selection==ONZ) cutflow1Photon=true;
                        //decisionMapCutFlowFine[PHOTON1]=true;
                        if(L.mll<101. && L.mll>81.){
                           if(selection==ONZ) cutflowOnZ=true;
                           decisionMapCutFlowFine[ZMASS]=true;
                        }
                     }
                     return true;
                     //else{
                        //return false;
                     //}
                  }
               }
            }
         }
      }else{return false;}
   }else{return false;}
  
  
  return false;
}

bool HistogramProducer::SelectEventTriggerStudies(selectionType selection){
   if (CheckParticles()){//only Fill Histograms if nMu>=2 or nEle>=2, nGamma>=0
      if (CleaningTriggerStudies()){//decide in events with DiEle and DiMu where to contribute
         selEvent L;
         L.isDiElectron=isDiElectron;
         L.trigDiEle=trigDiEle;
         L.trigDiMu=trigDiMu;
         L.trigMuEle=trigMuEle;
         L.trigHt=trigHt;
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
            L.totalWeight = totalWeight * GetScaleFactorAndError(pt1,eta1,isSignal,isDiElectron)*GetScaleFactorAndError(pt2,eta2,isSignal,isDiElectron);
            //L.totalWeight = totalWeight * GetScaleFactorAndErrorAlternative(pt1,eta1,isSignal,isDiElectron,*runNo)*GetScaleFactorAndErrorAlternative(pt2,eta2,isSignal,isDiElectron,*runNo);
         }else{
            L.totalWeight=totalWeight;
         }
         if (L.isDiElectron){
            auto e1 = electrons->at(0);
            auto e2 = electrons->at(1);
            if(testSelection(e1,selection,true) && testSelection(e2,selection,false) && (chargeProduct < 0.) && ((selection==UNCUT)? true : (mll>50.)) ){

               //for (vector<tree::Photon>::iterator it = photons->begin(); it != photons->end(); ++it){
               for (vector<tree::Photon>::iterator it = photons->begin(); it != photons->end(); it++){
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
               //for (vector<tree::Jet>::iterator it = jets->begin(); it != jets->end(); ++it){
               for (vector<tree::Jet>::iterator it = jets->begin(); it != jets->end(); it++){
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
               if(!GenPhotonVeto(config_veto)){
                  return true;
               }else{
                  return false;
               }
            }
         }
         else{
            auto m1 = muons->at(0); 
            auto m2 = muons->at(1);
            if(testSelection(m1,selection,true) && testSelection(m2,selection,false) && (chargeProduct < 0.) && ((selection==UNCUT)? true : (mll>50.))){
               //for (vector<tree::Photon>::iterator it = photons->begin(); it != photons->end(); ++it){
               for (vector<tree::Photon>::iterator it = photons->begin(); it != photons->end(); it++){
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
               //for (vector<tree::Jet>::iterator it = jets->begin(); it != jets->end(); ++it){
               for (vector<tree::Jet>::iterator it = jets->begin(); it != jets->end(); it++){
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
               if(!GenPhotonVeto(config_veto)){
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
   if (config_selectionsToProcessMap[UNCUT]==true){
      h1Maps["uncutEE"]=InitHistograms(UNCUT); 
      h1Maps["uncutMM"]=InitHistograms(UNCUT);
      h1Maps["uncut"]=InitHistograms(UNCUT);
   }
   if (config_selectionsToProcessMap[DILEP]==true){
      h1Maps["dilepEE"]=InitHistograms(DILEP); 
      h1Maps["dilepMM"]=InitHistograms(DILEP);
      h1Maps["dilep"]=InitHistograms(DILEP);
   }
   if (config_selectionsToProcessMap[PHOTON]==true){
      h1Maps["1photonEE"]=InitHistograms(PHOTON); 
      h1Maps["1photonMM"]=InitHistograms(PHOTON); 
      h1Maps["1photon"]=InitHistograms(PHOTON); 
   }
   if (config_selectionsToProcessMap[SEL]==true){
      h1Maps["selEE"]=InitHistograms(SEL);
      h1Maps["selMM"]=InitHistograms(SEL); 
      h1Maps["sel"]=InitHistograms(SEL); 
   }
   if (config_selectionsToProcessMap[ONZ]==true){
      h1Maps["onZEE"]=InitHistograms(ONZ);
      h1Maps["onZMM"]=InitHistograms(ONZ); 
      h1Maps["onZ"]=InitHistograms(ONZ);
   }
   if (config_selectionsToProcessMap[ONZMET]==true){
      h1Maps["onZMetEE"]=InitHistograms(ONZ);
      h1Maps["onZMetMM"]=InitHistograms(ONZ); //>300
      h1Maps["onZMet"]=InitHistograms(ONZ); 
      
      
      h1Maps["onZMet0100EE"]=InitHistograms(ONZ); //<100
      h1Maps["onZMet0100MM"]=InitHistograms(ONZ); 
      h1Maps["onZMet0100"]=InitHistograms(ONZ); 
      
      h1Maps["onZMet0200EE"]=InitHistograms(ONZ); //<200
      h1Maps["onZMet0200MM"]=InitHistograms(ONZ); 
      h1Maps["onZMet0200"]=InitHistograms(ONZ); 
      
      h1Maps["onZMet0300EE"]=InitHistograms(ONZ); //<200
      h1Maps["onZMet0300MM"]=InitHistograms(ONZ); 
      h1Maps["onZMet0300"]=InitHistograms(ONZ); 
      
      h1Maps["onZMet100EE"]=InitHistograms(ONZ); //>100
      h1Maps["onZMet100MM"]=InitHistograms(ONZ); 
      h1Maps["onZMet100"]=InitHistograms(ONZ); 
      
      h1Maps["onZMet200EE"]=InitHistograms(ONZ);//>200
      h1Maps["onZMet200MM"]=InitHistograms(ONZ); 
      h1Maps["onZMet200"]=InitHistograms(ONZ); 
      
      h1Maps["onZMet100200EE"]=InitHistograms(ONZ);//>100 <200
      h1Maps["onZMet100200MM"]=InitHistograms(ONZ); 
      h1Maps["onZMet100200"]=InitHistograms(ONZ); 
      
      h1Maps["onZMet100300EE"]=InitHistograms(ONZ);//>100 <300
      h1Maps["onZMet100300MM"]=InitHistograms(ONZ); 
      h1Maps["onZMet100300"]=InitHistograms(ONZ); 
      
      h1Maps["onZMet200300EE"]=InitHistograms(ONZ);//>200 <300
      h1Maps["onZMet200300MM"]=InitHistograms(ONZ); 
      h1Maps["onZMet200300"]=InitHistograms(ONZ); 
   }
   if (config_selectionsToProcessMap[ONZG]==true){  
      h1Maps["onZGEE"]=InitHistograms(ONZ);
      h1Maps["onZGMM"]=InitHistograms(ONZ); 
      h1Maps["onZG"]=InitHistograms(ONZ);
   }
   if (config_selectionsToProcessMap[ABOVEZG]==true){ 
      h1Maps["mllG110EE"]=InitHistograms(ONZ);
      h1Maps["mllG110MM"]=InitHistograms(ONZ); 
      h1Maps["mllG110"]=InitHistograms(ONZ); 
   }
   if (config_selectionsToProcessMap[EXO]==true){
      h1Maps["exoEE"]=InitHistograms(EXO);
      h1Maps["exoMM"]=InitHistograms(EXO); 
      h1Maps["exo"]=InitHistograms(EXO); 
   }
   if (config_selectionsToProcessMap[EGRegression]==true){
      h1Maps["EGRegressionEE"]=InitHistograms(EGRegression);
      h1Maps["EGRegressionMM"]=InitHistograms(EGRegression); 
      h1Maps["EGRegression"]=InitHistograms(EGRegression); 
   }
   if (config_selectionsToProcessMap[SEL]==true){
      h2Maps["selEE"]=Init2DHistograms(SEL);
      h2Maps["selMM"]=Init2DHistograms(SEL); 
      h2Maps["sel"]=Init2DHistograms(SEL); 
   }
   if (config_selectionsToProcessMap[ONZ]==true){
      h2Maps["onZEE"]=Init2DHistograms(ONZ);
      h2Maps["onZMM"]=Init2DHistograms(ONZ); 
      h2Maps["onZ"]=Init2DHistograms(ONZ); 
   }
   if (config_selectionsToProcessMap[DILEP]==true){
      h2Maps["dilepEE"]=Init2DHistograms(DILEP);
      h2Maps["dilepMM"]=Init2DHistograms(DILEP); 
      h2Maps["dilep"]=Init2DHistograms(DILEP); 
   }
   if (config_selectionsToProcessMap[EXO]==true){
      h2Maps["exoEE"]=Init2DHistograms(EXO);
      h2Maps["exoMM"]=Init2DHistograms(EXO); 
      h2Maps["exo"]=Init2DHistograms(EXO); 
   }
}

void HistogramProducer::InitSignalScanHistos(string masspoint){
   s1Maps[masspoint]["signal"]=InitSignalScanHistograms(ONZ);
   s1Maps[masspoint]["signalEE"]=InitSignalScanHistograms(ONZ);
   s1Maps[masspoint]["signalMM"]=InitSignalScanHistograms(ONZ);
}


void HistogramProducer::InitCutFlowHistos(){
   //c1Maps["dilepEE"]=InitCutFlowHistograms(DILEP);
   //c1Maps["dilepMM"]=InitCutFlowHistograms(DILEP);
   //c1Maps["selEE"]=InitCutFlowHistograms(SEL);
   //c1Maps["selMM"]=InitCutFlowHistograms(SEL);
   c1Maps["cutFlow_onZEE"]=InitCutFlowHistograms(ONZ);
   c1Maps["cutFlow_onZMM"]=InitCutFlowHistograms(ONZ);
}


void HistogramProducer::InitCutFlowHistos_Fine(){
   //c1Maps["dilepEE"]=InitCutFlowHistograms(DILEP);
   //c1Maps["dilepMM"]=InitCutFlowHistograms(DILEP);
   //c1Maps["cutFlow_Fine_selEE"]=InitCutFlowHistograms(SEL);
   //c1Maps["cutFlow_Fine_selMM"]=InitCutFlowHistograms(SEL);
   c1Maps["cutFlow_Fine_onZEE"]=InitCutFlowHistograms_Fine(ONZ);
   c1Maps["cutFlow_Fine_onZMM"]=InitCutFlowHistograms_Fine(ONZ);
}

   
void HistogramProducer::InitTriggerStudiesHistos(){
   eff1Maps["trigDilepEE"]=InitTriggerStudies(TRIGDILEP);
   eff1Maps["trigDilepMM"]=InitTriggerStudies(TRIGDILEP);
   eff1Maps["trigSelEE"]=InitTriggerStudies(TRIGSEL);
   eff1Maps["trigSelMM"]=InitTriggerStudies(TRIGSEL);
   eff1Maps["trigOnZEE"]=InitTriggerStudies(TRIGONZ);
   eff1Maps["trigOnZMM"]=InitTriggerStudies(TRIGONZ);
   eff1Maps["trigDilepEE_ptcuts"]=InitTriggerStudies(TRIGDILEP_ptcuts);
   eff1Maps["trigDilepMM_ptcuts"]=InitTriggerStudies(TRIGDILEP_ptcuts);
   eff1Maps["trigDilepEE_pt1cut"]=InitTriggerStudies(TRIGDILEP_ptcuts);
   eff1Maps["trigDilepMM_pt1cut"]=InitTriggerStudies(TRIGDILEP_ptcuts);
   eff1Maps["trigDilepEE_pt2cut"]=InitTriggerStudies(TRIGDILEP_ptcuts);
   eff1Maps["trigDilepMM_pt2cut"]=InitTriggerStudies(TRIGDILEP_ptcuts);
   eff1Maps["trigSelEE_ptcuts"]=InitTriggerStudies(TRIGSEL_ptcuts);
   eff1Maps["trigSelMM_ptcuts"]=InitTriggerStudies(TRIGSEL_ptcuts);
   eff1Maps["trigOnZEE_ptcuts"]=InitTriggerStudies(TRIGONZ_ptcuts);
   eff1Maps["trigOnZMM_ptcuts"]=InitTriggerStudies(TRIGONZ_ptcuts);

}

void HistogramProducer::FillHistograms2D(){
   //if(SelectEvent(SEL)){
   //if(SelectEvent(SEL)||SelectEvent(ONZ)){
      //if (selectedEvent.selPhotons.size()!=0){
         //auto m1 = &h2Maps["selEE"];
         //auto m2 = &h2Maps["selMM"];
         //auto m3 = &h2Maps["sel"];
         //m3->at(ISRVFSR).Fill(selectedEvent.mll,(selectedEvent.l1+selectedEvent.l2+selectedEvent.selPhotons.at(0).vec).M());         
         //if (selectedEvent.isDiElectron){ 
            //m1->at(ISRVFSR).Fill(selectedEvent.mll,(selectedEvent.l1+selectedEvent.l2+selectedEvent.selPhotons.at(0).vec).M());    
         //}else{
            //m2->at(ISRVFSR).Fill(selectedEvent.mll,(selectedEvent.l1+selectedEvent.l2+selectedEvent.selPhotons.at(0).vec).M());          
         //}
      //}   
   //
   //}

   if(SelectEvent(SEL)){
      if (selectedEvent.selPhotons.size()!=0){
         Filler2D(selectedEvent,h2Maps["sel"],true);
         if (selectedEvent.isDiElectron){
            Filler2D(selectedEvent,h2Maps["selEE"],true); 
         }else{
            if (selectedEvent.isDiMuon) Filler2D(selectedEvent,h2Maps["selMM"],true); 
         }
      }
   }
   
   if(SelectEvent(ONZ)){
      if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)){
         Filler2D(selectedEvent,h2Maps["onZ"],true);
         if (selectedEvent.isDiElectron){ 
            Filler2D(selectedEvent,h2Maps["onZEE"],true);
         }else{
            if (selectedEvent.isDiMuon) Filler2D(selectedEvent,h2Maps["onZMM"],true);
         }
      }
   }
   if(SelectEvent(EXO)){
      if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>50. && selectedEvent.mll<130.)){
         Filler2D(selectedEvent,h2Maps["exo"],true);
         if (selectedEvent.isDiElectron){ 
            Filler2D(selectedEvent,h2Maps["exoEE"],true);
         }else{
            if (selectedEvent.isDiMuon) Filler2D(selectedEvent,h2Maps["exoMM"],true);
         }
      }
   }
   
}


void HistogramProducer::Filler(selEvent& ev, map<Histograms1D,TH1F>& m,bool withPhoton){
   m.at(ETMISS).Fill(ev.ETmiss, ev.totalWeight);
   m.at(PT1).Fill(ev.pt1, ev.totalWeight);
   m.at(PT2).Fill(ev.pt2, ev.totalWeight);
   m.at(MLL).Fill(ev.mll,ev.totalWeight);   
   m.at(NPHOTONS).Fill(ev.selPhotons.size(),ev.totalWeight);
   m.at(NVTX).Fill(*nGoodVertices,ev.totalWeight);
   m.at(HT).Fill(*ht,ev.totalWeight);
   m.at(GENHT).Fill(*genHt,ev.totalWeight);
   m.at(NJETS).Fill(ev.selJets.size(),ev.totalWeight);
   m.at(ETA1).Fill(fabs(ev.eta1),ev.totalWeight);
   m.at(ETA2).Fill(fabs(ev.eta2),ev.totalWeight);
   m.at(PHI1).Fill(fabs(ev.phi2),ev.totalWeight);
   m.at(PHI2).Fill(fabs(ev.phi2),ev.totalWeight);
   m.at(DeltaEtaLL).Fill(fabs(ev.eta1-ev.eta2),ev.totalWeight);
   m.at(DeltaPhiLL).Fill(fabs(ev.phi1-ev.phi2),ev.totalWeight);
   m.at(DeltaRLL).Fill(fabs(ev.deltaRll),ev.totalWeight);
   m.at(ZPT).Fill((ev.l1+ev.l2).Pt(),ev.totalWeight);
   m.at(MTLL).Fill((ev.l1+ev.l2).Mt(),ev.totalWeight);
   m.at(ST).Fill(ev.pt1+ev.pt2,ev.totalWeight);
   m.at(MT2).Fill(ev.MT2_val,ev.totalWeight);
   m.at(DeltaPhiLLMet).Fill(fabs((ev.l1+ev.l2).Phi() - ev.ETmiss_vec.Phi()),ev.totalWeight);
   m.at(DeltaEtaLLMet).Fill(fabs((ev.l1+ev.l2).Eta() - ev.ETmiss_vec.Eta()),ev.totalWeight);
   m.at(DeltaRLLMet).Fill(fabs((ev.l1+ev.l2).DeltaR(ev.ETmiss_vec)),ev.totalWeight);
   
   if(ev.evtHasGenPhotonVeto) m.at(VetoCompare).Fill(1.,ev.totalWeight);
   if(!ev.evtHasGenPhotonVeto) m.at(VetoCompare).Fill(0.,ev.totalWeight);
   
   if(withPhoton){
      m.at(PTG1).Fill(ev.selPhotons.at(0).p.Pt(),ev.totalWeight);
      m.at(PHIG1).Fill(ev.selPhotons.at(0).p.Phi(),ev.totalWeight);
      m.at(ETAG1).Fill(ev.selPhotons.at(0).p.Eta(),ev.totalWeight);
      m.at(SIGMAIETAIETAG1).Fill(ev.selPhotons.at(0).sigmaIetaIeta,ev.totalWeight);
      m.at(SIGMAIPHIIPHIG1).Fill(ev.selPhotons.at(0).sigmaIphiIphi,ev.totalWeight);
      m.at(R9).Fill(ev.selPhotons.at(0).r9,ev.totalWeight);
      m.at(HOVERE).Fill(ev.selPhotons.at(0).hOverE,ev.totalWeight);
      m.at(DELTARGL1).Fill(ev.selPhotons.at(0).deltaR1,ev.totalWeight);
      m.at(DELTARGL2).Fill(ev.selPhotons.at(0).deltaR2,ev.totalWeight);
      m.at(DeltaEtaLLG).Fill(fabs((ev.l1+ev.l2).Eta()-ev.selPhotons.at(0).vec.Eta()), ev.totalWeight);
      m.at(DeltaPhiLLG).Fill(fabs((ev.l1+ev.l2).Phi()-ev.selPhotons.at(0).vec.Phi()),ev.totalWeight);
      m.at(DeltaRLLG).Fill(fabs((ev.l1+ev.l2).DeltaR(ev.selPhotons.at(0).vec)),ev.totalWeight);
      m.at(MTLLG).Fill((ev.l1+ev.l2+ev.selPhotons.at(0).vec).Mt(),ev.totalWeight);
      m.at(MTL1MET).Fill((ev.l1+ev.ETmiss_vec).Mt(),ev.totalWeight);
      m.at(MTL2MET).Fill((ev.l2+ev.ETmiss_vec).Mt(),ev.totalWeight);
      m.at(MTGMET).Fill((ev.selPhotons.at(0).vec+ev.ETmiss_vec).Mt(),ev.totalWeight);
      m.at(MTLLMET).Fill((ev.l1+ev.l2+ev.ETmiss_vec).Mt(),ev.totalWeight);
      m.at(MTLLGMET).Fill((ev.l1+ev.l2+ev.selPhotons.at(0).vec+ev.ETmiss_vec).Mt(),ev.totalWeight);
      m.at(STG).Fill(ev.pt1+ev.pt2+ev.selPhotons.at(0).p.Pt(),ev.totalWeight);
      m.at(STMET).Fill(ev.pt1+ev.pt2+ev.selPhotons.at(0).p.Pt()+ev.ETmiss,ev.totalWeight);
      m.at(MLLG).Fill((ev.l1+ev.l2+ev.selPhotons.at(0).vec).M(),ev.totalWeight);
      //m.at(PT_llg).Fill((ev.l1+ev.l2+ev.selPhotons.at(0).vec).Pt(),ev.totalWeight);
      //m.at(MZG_exo).Fill(40./150. * (ev.l1+ev.l2+ev.selPhotons.at(0).vec).M(),ev.totalWeight);
      
      if(FindGenPhotonMatch(ev.selPhotons.at(0))) m.at(genPhotonPT).Fill(GetGenPhotonMatch(ev.selPhotons.at(0)).p.Pt(),ev.totalWeight);
      
      if(ev.evtHasGenPhotonVeto) if(FindGenPhotonMatch(ev.selPhotons.at(0))) m.at(genPhotonPT_Veto).Fill(GetGenPhotonMatch(ev.selPhotons.at(0)).p.Pt(),ev.totalWeight);
      if(ev.evtHasGenPhotonVeto) m.at(PTG1_Veto).Fill(ev.selPhotons.at(0).p.Pt(),ev.totalWeight);
      if(!ev.evtHasGenPhotonVeto) if(FindGenPhotonMatch(ev.selPhotons.at(0))) m.at(genPhotonPT_NoVeto).Fill(GetGenPhotonMatch(ev.selPhotons.at(0)).p.Pt(),ev.totalWeight);
      if(!ev.evtHasGenPhotonVeto) m.at(PTG1_NoVeto).Fill(ev.selPhotons.at(0).p.Pt(),ev.totalWeight);
      //if(ev.evtHasGenPhotonVeto) m.at(VetoCompare).Fill(1.,ev.totalWeight);
      //if(!ev.evtHasGenPhotonVeto) m.at(VetoCompare).Fill(0.,ev.totalWeight);
      
      //for (vector<tree::GenParticle>::iterator it = genParticles->begin(); it != genParticles->end(); it++){
         //int Id_mother=-1;
         //if (abs(it->pdgId)==22){
            //Id_mother = abs(it->motherId);
            //m.at(gammaMotherID).Fill(Id_mother,ev.totalWeight);
         //}
      //}
      
   }
}
void HistogramProducer::Filler2D(selEvent& ev, map<Histograms2D,TH2F>& m,bool withPhoton){
   if(withPhoton){
      //m.at(PTG1).Fill(ev.selPhotons.at(0).p.Pt(),ev.totalWeight);
      m.at(ISRVFSR).Fill(ev.mll,(ev.l1+selectedEvent.l2+ev.selPhotons.at(0).vec).M(),ev.totalWeight);  
      m.at(PTGvsMLLG).Fill((ev.l1+selectedEvent.l2+ev.selPhotons.at(0).vec).M(),ev.selPhotons.at(0).p.Pt(),ev.totalWeight);  
   }
}

void HistogramProducer::FillerTrigger(selEvent& ev, map<Histograms1D,TEfficiency>& m,bool withPhoton,bool TriggerBool){
   m.at(ETMISS).Fill(TriggerBool,ev.ETmiss);
   m.at(PT1).Fill(TriggerBool,ev.pt1);
   m.at(PT2).Fill(TriggerBool,ev.pt2);
   m.at(MLL).Fill(TriggerBool,ev.mll);   
   m.at(NPHOTONS).Fill(TriggerBool,ev.selPhotons.size());
   m.at(NVTX).Fill(TriggerBool,*nGoodVertices);
   m.at(HT).Fill(TriggerBool,*ht);
   m.at(GENHT).Fill(TriggerBool,*genHt);
   m.at(NJETS).Fill(TriggerBool,ev.selJets.size());
   m.at(ETA1).Fill(TriggerBool,fabs(ev.eta1));
   m.at(ETA2).Fill(TriggerBool,fabs(ev.eta2));
   m.at(PHI1).Fill(TriggerBool,fabs(ev.phi2));
   m.at(PHI2).Fill(TriggerBool,fabs(ev.phi2));
   if(withPhoton){
      m.at(PTG1).Fill(TriggerBool,ev.selPhotons.at(0).p.Pt());
      m.at(PHIG1).Fill(TriggerBool,fabs(ev.selPhotons.at(0).p.Phi()));
      m.at(ETAG1).Fill(TriggerBool,fabs(ev.selPhotons.at(0).p.Eta()));
      m.at(SIGMAIETAIETAG1).Fill(TriggerBool,ev.selPhotons.at(0).sigmaIetaIeta);
   }
}

void HistogramProducer::FillerSignal(selEvent& ev, map<Histograms1D,TH1F>& m){
   m.at(ETMISS).Fill(ev.ETmiss, ev.totalWeight*1./nGen);
}

void HistogramProducer::FillHistograms(){

   if(SelectEvent(UNCUT)&&config_selectionsToProcessMap[UNCUT]){
      Filler(selectedEvent,h1Maps["uncut"],false); 
      if (selectedEvent.isDiElectron){ 
         Filler(selectedEvent,h1Maps["uncutEE"],false); 
      }else{
         Filler(selectedEvent,h1Maps["uncutMM"],false);      
      }
   }
   if(SelectEvent(DILEP)&&config_selectionsToProcessMap[DILEP]){
      Filler(selectedEvent,h1Maps["dilep"],false);
      if (selectedEvent.isDiElectron){ 
         Filler(selectedEvent,h1Maps["dilepEE"],false);
      }else{
        if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["dilepMM"],false);
      }
   }
   if(SelectEvent(PHOTON)&&config_selectionsToProcessMap[PHOTON]){
      if (selectedEvent.selPhotons.size()!=0){
         Filler(selectedEvent,h1Maps["1photon"],true);         
         if (selectedEvent.isDiElectron){ 
            Filler(selectedEvent,h1Maps["1photonEE"],true);
         }else{
            Filler(selectedEvent,h1Maps["1photonMM"],true); 
         }
      }
   }
   if(SelectEvent(SEL)&&config_selectionsToProcessMap[SEL]){
      if (selectedEvent.selPhotons.size()!=0){
         Filler(selectedEvent,h1Maps["sel"],true);
         if (selectedEvent.isDiElectron){
            Filler(selectedEvent,h1Maps["selEE"],true); 
         }else{
            if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["selMM"],true); 
         }
      }
   }
   clearCutFlowMap();
   if(SelectEvent(ONZ)&&config_selectionsToProcessMap[ONZ]){
        
      if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)){
         Filler(selectedEvent,h1Maps["onZ"],true);
         if (selectedEvent.isDiElectron){ 
            Filler(selectedEvent,h1Maps["onZEE"],true);
         }else{
            if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["onZMM"],true);
         }
      }
   }
   clearCutFlowMap();
   if(SelectEvent(EXO)&&config_selectionsToProcessMap[EXO]){
        
      if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>50. && selectedEvent.mll<130.)){
         Filler(selectedEvent,h1Maps["exo"],true);
         if (selectedEvent.isDiElectron){ 
            Filler(selectedEvent,h1Maps["exoEE"],true);
         }else{
            if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["exoMM"],true);
         }
      }
   }
   if(SelectEvent(EGRegression)&&config_selectionsToProcessMap[EGRegression]){
      Filler(selectedEvent,h1Maps["EGRegression"],false);
      if (selectedEvent.isDiElectron){ 
         Filler(selectedEvent,h1Maps["EGRegressionEE"],false);
      }else{
         if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["EGRegressionMM"],false);
      }
   }
   if(SelectEvent(ONZ)&&config_selectionsToProcessMap[ONZMET]){ //ONZ+MET>300
        
      if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)&&(selectedEvent.ETmiss>=300.)){
         Filler(selectedEvent,h1Maps["onZMet"],true);
         if (selectedEvent.isDiElectron){ 
            Filler(selectedEvent,h1Maps["onZMetEE"],true);
         }else{
            if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["onZMetMM"],true);
         }
      }

      if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)&&(selectedEvent.ETmiss<100.)){ //ONZ+MET<100
         Filler(selectedEvent,h1Maps["onZMet0100"],true);
         if (selectedEvent.isDiElectron){ 
            Filler(selectedEvent,h1Maps["onZMet0100EE"],true);
         }else{
            if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["onZMet0100MM"],true);
         }
      }
      if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)&&(selectedEvent.ETmiss<200.)){ //ONZ+MET>200
         Filler(selectedEvent,h1Maps["onZMet0200"],true);
         if (selectedEvent.isDiElectron){ 
            Filler(selectedEvent,h1Maps["onZMet0200EE"],true);
         }else{
            if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["onZMet0200MM"],true);
         }
      }
      if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)&&(selectedEvent.ETmiss<300.)){ //ONZ+MET>300
         Filler(selectedEvent,h1Maps["onZMet0300"],true);
         if (selectedEvent.isDiElectron){ 
            Filler(selectedEvent,h1Maps["onZMet0200EE"],true);
         }else{
            if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["onZMet0300MM"],true);
         }
      }
      
      if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)&&(selectedEvent.ETmiss>=100.)){ //ONZ+MET>100
         Filler(selectedEvent,h1Maps["onZMet100"],true);
         if (selectedEvent.isDiElectron){ 
            Filler(selectedEvent,h1Maps["onZMet100EE"],true);
         }else{
            if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["onZMet100MM"],true);
         }
      }      
      
       if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)&&(selectedEvent.ETmiss>=200.)){ //ONZ+MET>200
         Filler(selectedEvent,h1Maps["onZMet200"],true);
         if (selectedEvent.isDiElectron){ 
            Filler(selectedEvent,h1Maps["onZMet200EE"],true);
         }else{
            if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["onZMet200MM"],true);
         }
      }    
      
       if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)&&(selectedEvent.ETmiss>=100.)&&(selectedEvent.ETmiss<200.)){ //ONZ+MET>100<200
         Filler(selectedEvent,h1Maps["onZMet100200"],true);
         if (selectedEvent.isDiElectron){ 
            Filler(selectedEvent,h1Maps["onZMet100200EE"],true);
         }else{
            if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["onZMet100200MM"],true);
         }
      }
      
      if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)&&(selectedEvent.ETmiss>=100.)&&(selectedEvent.ETmiss<300.)){//ONZ+MET>100<300
         Filler(selectedEvent,h1Maps["onZMet100300"],true);
         if (selectedEvent.isDiElectron){ 
            Filler(selectedEvent,h1Maps["onZMet100300EE"],true);
         }else{
            if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["onZMet100300MM"],true);
         }
      }      
      
       if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)&&(selectedEvent.ETmiss>=200.)&&(selectedEvent.ETmiss<300.)){ //ONZ+MET>200<300
         Filler(selectedEvent,h1Maps["onZMet200300"],true);
         if (selectedEvent.isDiElectron){ 
            Filler(selectedEvent,h1Maps["onZMet200300EE"],true);
         }else{
            if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["onZMet200300MM"],true);
         }
      }     
   }
   if(SelectEvent(ONZG)&&config_selectionsToProcessMap[ONZG]){ //ONZG
        
      if ((selectedEvent.selPhotons.size()!=0)&&((selectedEvent.l1+selectedEvent.l2+selectedEvent.selPhotons.at(0).vec).M()>81.)&&((selectedEvent.l1+selectedEvent.l2+selectedEvent.selPhotons.at(0).vec).M()<101.)){
         Filler(selectedEvent,h1Maps["onZG"],true);
         if (selectedEvent.isDiElectron){ 
            Filler(selectedEvent,h1Maps["onZGEE"],true);
         }else{
            if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["onZGMM"],true);
         }
      }
   }
   if(SelectEvent(ONZ)&&config_selectionsToProcessMap[ABOVEZG]){//MLLG>110
        
      if ((selectedEvent.selPhotons.size()!=0)&&((selectedEvent.l1+selectedEvent.l2+selectedEvent.selPhotons.at(0).vec).M()>110.)){
         Filler(selectedEvent,h1Maps["mllG110"],true);
         if (selectedEvent.isDiElectron){ 
            Filler(selectedEvent,h1Maps["mllG110EE"],true);
         }else{
            if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["mllG110MM"],true);
         }
      }
   }
}

void HistogramProducer::FillSignalHistograms(){
   if(SelectEvent(ONZ)){   
      if((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)){
         string signalPoint = getSignalPointName(*nBinos,*signal_m1,*signal_m2);
         if(!(s1Maps.count(signalPoint)>0)){
            InitSignalScanHistos(signalPoint);
         }
         FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("signal"));
         if (selectedEvent.isDiElectron){
            FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("signalEE"));
         }else{
            if (selectedEvent.isDiMuon) FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("signalMM"));
         }
      }
   }
}


void HistogramProducer::FillCutFlowHistograms(){
   if(cutflowDiEle && !(cutflowDiMu)){
      if (cutflowIsTriggered) c1Maps["cutFlow_onZEE"].at(CUTFLOW).Fill("triggered",selectedEvent.totalWeight);
      if (cutflow2Leptons) c1Maps["cutFlow_onZEE"].at(CUTFLOW).Fill("2leptons",selectedEvent.totalWeight);
      if (cutflowMll50) c1Maps["cutFlow_onZEE"].at(CUTFLOW).Fill("m50",selectedEvent.totalWeight);
      if (cutflow1Photon) c1Maps["cutFlow_onZEE"].at(CUTFLOW).Fill("1photon",selectedEvent.totalWeight);
      if (cutflowOnZ) c1Maps["cutFlow_onZEE"].at(CUTFLOW).Fill("Z",selectedEvent.totalWeight);
   }
   if(cutflowDiMu && !(cutflowDiEle)){
      if (cutflowIsTriggered) c1Maps["cutFlow_onZMM"].at(CUTFLOW).Fill("triggered",selectedEvent.totalWeight);
      if (cutflow2Leptons) c1Maps["cutFlow_onZMM"].at(CUTFLOW).Fill("2leptons",selectedEvent.totalWeight);
      if (cutflowMll50) c1Maps["cutFlow_onZMM"].at(CUTFLOW).Fill("m50",selectedEvent.totalWeight);
      if (cutflow1Photon) c1Maps["cutFlow_onZMM"].at(CUTFLOW).Fill("1photon",selectedEvent.totalWeight);
      if (cutflowOnZ) c1Maps["cutFlow_onZMM"].at(CUTFLOW).Fill("Z",selectedEvent.totalWeight);
   }

}
void HistogramProducer::SetCutFlowHistogramsStatus(){
   cutflowIsTriggered=false;
   cutflow2Leptons=false;
   cutflowMll50=false;
   cutflow1Photon=false;
   cutflowOnZ=false;
}




void HistogramProducer::FillCutFlowHistograms_Fine(){
   auto m=decisionMapCutFlowFine;
   auto mW=decisionMapCutFlowFine_weight;
   if(m[DIELECTRON] && !(m[DIMUON])){
      //if (m[TRIGGERED]) c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("triggered",selectedEvent.totalWeight);
      //if (m[TRIGGERED]&&m[GENVETO]) c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("GenPhotonVeto",selectedEvent.totalWeight);
      //if (m[TRIGGERED]&&m[GENVETO]&&m[LEPTONID_leading]&&m[LEPTONID_trailing]) c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("2LeptonID",selectedEvent.totalWeight);
      //if (m[TRIGGERED]&&m[GENVETO]&&m[LEPTONID_leading]&&m[LEPTONID_trailing]&&m[LEPTONPT_leading]&&m[LEPTONPT_trailing]) c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("2LeptonPT",selectedEvent.totalWeight);
      //if (m[M50]) c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("m50",selectedEvent.totalWeight);
      //if (m[PHOTON1]) c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("1Photon",selectedEvent.totalWeight);
      //if (m[PHOTON1ID]) c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("1PhotonID",selectedEvent.totalWeight);
      //if (m[PHOTON1ID]&&m[PHOTON1PT]) c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("1PhotonPT",selectedEvent.totalWeight);
      //if (m[PHOTON1ID]&&m[PHOTON1PT]&&m[PHOTON1DR]) c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("1PhotonDeltaR",selectedEvent.totalWeight);
      //if (m[ZMASS]) c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("Z",selectedEvent.totalWeight);
   
      if (m[TRIGGERED]){
         c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("triggered",m[TRIGGERED]);
         if(m[GENVETO]){
            c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("GenPhotonVeto",mW[GENVETO]);
            if(m[LEPTONID_leading] && m[LEPTONID_trailing]){
               c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("2LeptonID",mW[LEPTONID_leading]);
               if(m[LEPTONPT_leading] && m[LEPTONPT_trailing]){
                  c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("2LeptonPT",mW[LEPTONPT_leading]);
                  if(m[M50]){
                     c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("m50",mW[M50]);
                     if(m[PHOTON1]){
                        c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("1Photon",mW[PHOTON1]);
                        if(m[PHOTON1ID]){
                           c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("1PhotonID",mW[PHOTON1ID]);
                           if(m[PHOTON1PT]){
                              c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("1PhotonPT",mW[PHOTON1PT]);
                              if(m[PHOTON1DR]){
                                 c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("1PhotonDeltaR",mW[PHOTON1DR]);
                                 if(m[ZMASS]){
                                    c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("Z",mW[ZMASS]);
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   if(m[DIMUON] && !(m[DIELECTRON])){
      //if (m[TRIGGERED]) c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("triggered",selectedEvent.totalWeight);
      //if (m[GENVETO]) c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("GenPhotonVeto",selectedEvent.totalWeight);
      //if (m[TRIGGERED]&&m[GENVETO]) c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("GenPhotonVeto",selectedEvent.totalWeight);
      //if (m[LEPTONID_leading]&&m[LEPTONID_trailing]) c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("2LeptonID",selectedEvent.totalWeight);
      //if (m[TRIGGERED]&&m[GENVETO]&&m[LEPTONID_leading]&&m[LEPTONID_trailing]) c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("2LeptonID",selectedEvent.totalWeight);
      //if (m[LEPTONPT_leading]&&m[LEPTONPT_trailing]) c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("2LeptonPT",selectedEvent.totalWeight);
      //if (m[TRIGGERED]&&m[GENVETO]&&m[LEPTONID_leading]&&m[LEPTONID_trailing]&&m[LEPTONPT_leading]&&m[LEPTONPT_trailing]) c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("2LeptonPT",selectedEvent.totalWeight);
      //if (m[TRIGGERED]&&m[GENVETO]&&m[LEPTONID_leading]&&m[LEPTONID_trailing]&&m[LEPTONPT_leading]&&m[LEPTONPT_trailing]&&m[M50]) c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("m50",selectedEvent.totalWeight);
      //if (m[TRIGGERED]&&m[GENVETO]&&m[LEPTONID_leading]&&m[LEPTONID_trailing]&&m[LEPTONPT_leading]&&m[LEPTONPT_trailing]&&m[M50]&&m[PHOTON1]) c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("1Photon",selectedEvent.totalWeight);
      //if (m[TRIGGERED]&&m[GENVETO]&&m[LEPTONID_leading]&&m[LEPTONID_trailing]&&m[LEPTONPT_leading]&&m[LEPTONPT_trailing]&&m[M50]&&m[PHOTON1]&&m[PHOTON1ID]) c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("1PhotonID",selectedEvent.totalWeight);
      //if (m[TRIGGERED]&&m[GENVETO]&&m[LEPTONID_leading]&&m[LEPTONID_trailing]&&m[LEPTONPT_leading]&&m[LEPTONPT_trailing]&&m[M50]&&m[PHOTON1]&&m[PHOTON1ID]&&m[PHOTON1PT]) c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("1PhotonPT",selectedEvent.totalWeight);
      //if (m[TRIGGERED]&&m[GENVETO]&&m[LEPTONID_leading]&&m[LEPTONID_trailing]&&m[LEPTONPT_leading]&&m[LEPTONPT_trailing]&&m[M50]&&m[PHOTON1]&&m[PHOTON1ID]&&m[PHOTON1PT]&&m[PHOTON1DR]) c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("1PhotonDeltaR",selectedEvent.totalWeight);
      //if (m[TRIGGERED]&&m[GENVETO]&&m[LEPTONID_leading]&&m[LEPTONID_trailing]&&m[LEPTONPT_leading]&&m[LEPTONPT_trailing]&&m[M50]&&m[PHOTON1]&&m[PHOTON1ID]&&m[PHOTON1PT]&&m[PHOTON1DR]&&m[ZMASS]) c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("Z",selectedEvent.totalWeight);
   
      if (m[TRIGGERED]){
         c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("triggered",selectedEvent.totalWeight);
         if(m[GENVETO]){
            c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("GenPhotonVeto",selectedEvent.totalWeight);
            if(m[LEPTONID_leading] && m[LEPTONID_trailing]){
               c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("2LeptonID",selectedEvent.totalWeight);
               if(m[LEPTONPT_leading] && m[LEPTONPT_trailing]){
                  c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("2LeptonPT",selectedEvent.totalWeight);
                  if(m[M50]){
                     c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("m50",selectedEvent.totalWeight);
                     if(m[PHOTON1]){
                        c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("1Photon",selectedEvent.totalWeight);
                        if(m[PHOTON1ID]){
                           c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("1PhotonID",selectedEvent.totalWeight);
                           if(m[PHOTON1PT]){
                              c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("1PhotonPT",selectedEvent.totalWeight);
                              if(m[PHOTON1DR]){
                                 c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("1PhotonDeltaR",selectedEvent.totalWeight);
                                 if(m[ZMASS]){
                                    c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("Z",selectedEvent.totalWeight);
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   
   
   
   }
}



void HistogramProducer::FillTriggerStudies(){
   if(SelectEventTriggerStudies(TRIGDILEP)){
      //auto m1 = &eff1Maps["trigDilepEE"];
      //auto m2 = &eff1Maps["trigDilepMM"];
      if (selectedEvent.trigHt){//baselineTrigger
         if (selectedEvent.isDiElectron){ 
           FillerTrigger(selectedEvent,eff1Maps["trigDilepEE"],false,selectedEvent.trigDiEle);
         }else{
            FillerTrigger(selectedEvent,eff1Maps["trigDilepMM"],false,selectedEvent.trigDiMu);  
         }
      }
   }
   if(SelectEventTriggerStudies(TRIGDILEP_ptcuts)){
      if (selectedEvent.trigHt){//baselineTrigger
         if (selectedEvent.isDiElectron){ 
           FillerTrigger(selectedEvent,eff1Maps["trigDilepEE_ptcuts"],false,selectedEvent.trigDiEle);
         }else{
            FillerTrigger(selectedEvent,eff1Maps["trigDilepMM_ptcuts"],false,selectedEvent.trigDiMu);  
         }
      }
   }
   if(SelectEventTriggerStudies(TRIGDILEP_pt1cut)){
      if (selectedEvent.trigHt){//baselineTrigger
         if (selectedEvent.isDiElectron){ 
           FillerTrigger(selectedEvent,eff1Maps["trigDilepEE_pt1cut"],false,selectedEvent.trigDiEle);
         }else{
            FillerTrigger(selectedEvent,eff1Maps["trigDilepMM_pt1cut"],false,selectedEvent.trigDiMu);  
         }
      }
   }
   if(SelectEventTriggerStudies(TRIGDILEP_pt2cut)){
      if (selectedEvent.trigHt){//baselineTrigger
         if (selectedEvent.isDiElectron){ 
           FillerTrigger(selectedEvent,eff1Maps["trigDilepEE_pt2cut"],false,selectedEvent.trigDiEle);
         }else{
            FillerTrigger(selectedEvent,eff1Maps["trigDilepMM_pt2cut"],false,selectedEvent.trigDiMu);  
         }
      }
   }
   if(SelectEventTriggerStudies(TRIGSEL)){
      if (selectedEvent.selPhotons.size()!=0){
         if (selectedEvent.trigHt){//baselineTrigger
            if (selectedEvent.isDiElectron){
               FillerTrigger(selectedEvent,eff1Maps["trigSelEE"],true,selectedEvent.trigDiEle); 
            }else{
               FillerTrigger(selectedEvent,eff1Maps["trigSelMM"],true,selectedEvent.trigDiMu); 
            }
         }
      }
   }
   if(SelectEventTriggerStudies(TRIGSEL_ptcuts)){
      if (selectedEvent.selPhotons.size()!=0){
         if (selectedEvent.trigHt){//baselineTrigger
            if (selectedEvent.isDiElectron){
               FillerTrigger(selectedEvent,eff1Maps["trigSelEE_ptcuts"],true,selectedEvent.trigDiEle); 
            }else{
               FillerTrigger(selectedEvent,eff1Maps["trigSelMM_ptcuts"],true,selectedEvent.trigDiMu); 
            }
         }
      }
   }
   if(SelectEventTriggerStudies(TRIGONZ)){
      if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81 && selectedEvent.mll<101)){
         if (selectedEvent.trigHt){//baselineTrigger
            if (selectedEvent.isDiElectron){
               FillerTrigger(selectedEvent,eff1Maps["trigOnZEE"],true,selectedEvent.trigDiEle);  
            }else{
               FillerTrigger(selectedEvent,eff1Maps["trigOnZMM"],true,selectedEvent.trigDiMu); 
            }
         }
      }
   }
   if(SelectEventTriggerStudies(TRIGONZ_ptcuts)){
      if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81 && selectedEvent.mll<101)){
         if (selectedEvent.trigHt){//baselineTrigger
            if (selectedEvent.isDiElectron){
               FillerTrigger(selectedEvent,eff1Maps["trigOnZEE_ptcuts"],true,selectedEvent.trigDiEle);  
            }else{
               FillerTrigger(selectedEvent,eff1Maps["trigOnZMM_ptcuts"],true,selectedEvent.trigDiMu); 
            }
         }
      }
   }
}

bool HistogramProducer::FindGenPhotonMatch(const selPhoton& pa){
   bool genMatch = false;
   for (auto& gen : *genParticles) {
      //if (gen.pdgId == 22 &&  gen.p.DeltaR(pa.p) < 0.3 && fabs(gen.p.Pt() - pa.p.Pt()) < 10. ) {
      if ((abs(gen.pdgId) == 22)  && (gen.statusID==1) &&  (gen.p.DeltaR(pa.p) < 0.2)) {
         genMatch = true;
         break;
      }
   }
   return genMatch;
}

tree::GenParticle& HistogramProducer::GetGenPhotonMatch(const selPhoton& pa){
   bool genMatch = false;
   auto& foundGen = genParticles->at(0);
   for (auto& gen : *genParticles) {
      if ((abs(gen.pdgId) == 22)  && (gen.statusID==1) &&  (gen.p.DeltaR(pa.p) < 0.2)) {
         genMatch = true;
         foundGen=gen;
         break;
      }
   }
   return foundGen;
}


Bool_t HistogramProducer::Process(Long64_t entry){
   float tempPercentage = (float) entry/ (float)nEntries;
   if(!(abs(config_eventpercentage-100.)<0.1)){
      if (tempPercentage>config_eventpercentage/100.){
         return kTRUE;
      }
   }
   fReader.SetLocalEntry(entry);
   totalWeight = *mc_weight * *pu_weight;
   
   //cout<<entry<<" | "<<tempPercentage<<" | "<<*evtNo<<endl;


   //int progress = tempPercentage*100.;
   //if(entry%100000==0){
   //
		//std::cout<<"[";
		//for(int i=0;i<100;i++)
			//if(i<progress)
				//std::cout<<'=';
			//else if(i==progress)
				//std::cout<<'>';
			//else
				//std::cout<<' ';
		//std::cout<<"] "<<progress<<" %"<<" "<<getOutputFilename(inputName)<<'\r';
		//std::cout.flush();
	//
//}

   string cutFlowName = "TreeWriter/hCutFlow";
   
   if(isTotalSignal){
      if(inputName.find("TChiNG")!=string::npos){
         cutFlowName+="TChiNG";
      }else{
         if(inputName.find("T5bbbbZg")!=string::npos){
            cutFlowName+="T5bbbbZg";
         }
      }
      cutFlowName += "_"+to_string(*signal_m1);
      if (*signal_m2) cutFlowName += "_"+to_string(*signal_m2);
      cutFlow = *((TH1F*)fReader.GetTree()->GetCurrentFile()->Get(cutFlowName.c_str()));
      if (cutFlow.GetEntries()) {
         nGen = cutFlow.GetBinContent(2);
      }else{
         cout << "Could not read cutFlow histogram " << cutFlowName << endl;
      }
   }

   SetCutFlowHistogramsStatus();
   clearCutFlowMap();
   //cout<<"init "<<decisionMapCutFlowFine[TRIGGERED]<<endl;
   FillHistograms();
   //FillHistograms2D();
   //cout<<cutflowIsTriggered<<" | "<<decisionMapCutFlowFine[TRIGGERED]<<endl;
   if(inputName.find("JetHT")!= string::npos){
      FillTriggerStudies();
   }
   
   if (config_docutflow) FillCutFlowHistograms();
   if (config_docutflowfine) FillCutFlowHistograms_Fine();
   clearCutFlowMap();
   if (config_dosignalscan) FillSignalHistograms();
   
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
template<typename T>
void save2File(const map<string,map<string,map<Histograms1D,T>>>& hMaps, TFile& file)
{
   for (auto& hMapIt : hMaps) {
      if (!file.Get(hMapIt.first.c_str())) {
         file.mkdir(hMapIt.first.c_str());
      }
      file.cd(hMapIt.first.c_str());
      for(auto& hMapIt2 : hMapIt.second){
         if (!file.Get(hMapIt2.first.c_str())) {
            file.mkdir((hMapIt.first+"/"+hMapIt2.first).c_str());
         }
         file.cd((hMapIt.first+"/"+hMapIt2.first).c_str());
         for (auto h: hMapIt2.second){
            //h.second.Scale(1./nEventMap[hMapMapIt.first]);
            h.second.Write(histoNames[h.first].c_str(), TObject::kWriteDelete);
         }
         file.cd();
      }
   }
}



void HistogramProducer::Terminate()
{
   //auto outputName = "output/"+getOutputFilename(inputName);
   //auto outputName = "output_EGRegression/"+getOutputFilename(inputName);
   //auto outputName = "output_noVeto/"+getOutputFilename(inputName);
   //auto outputName = "output_mllWeight/"+getOutputFilename(inputName);
   //auto outputName = "output_FSRVeto/"+getOutputFilename(inputName);
   //auto outputName = "output_FSR+HardVeto/"+getOutputFilename(inputName);
   
   auto outputName = "output"+config_outputfolder+"/"+getOutputFilename(inputName);
   
   TFile file(outputName.c_str(), "RECREATE");
   save2File(h1Maps, file);
   save2File(eff1Maps, file);
   save2File(h2Maps, file);
   save2File(c1Maps, file);
   save2File(s1Maps, file);
   if(! isTotalSignal){
      cutFlow.Write("hCutFlow");
   }
   file.Close();
   cout << "Created " << outputName << " in " << (time(NULL) - startTime)/60 << " min" << endl;
   //cout << "Created " << getOutputFilename(inputName) << " in " << (time(NULL) - startTime)/60 << " min" << endl;
}

