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
   EWKinoPairPt(fReader, "EWKinoPairPt"),
   leptonPairPt(fReader, "leptonPairPt"),
   topPt1(fReader, "topPt1"),
   topPt2(fReader, "topPt2"),
   
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





void HistogramProducer::CalculateVariables(const tree::Lepton& l1, const tree::Lepton& l2, const tree::Particle mett ,const particleType particle=DUMMYPARTICLE)
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
   bool isEMuon;
   bool isMuonE;
   if (particle==EM){
      isEMuon=true;
   }else{
      isEMuon=false;
   }
   if (particle==ME){
      isMuonE=true;
   }else{
      isMuonE=false;
   }
   //isMuonElectron=isEMuon;
   //isElectronMuon=isMuonE;
   isMuonElectron=isMuonE;
   isElectronMuon=isEMuon;
   float leptonMass1 = (isEle||isEMuon)?  0.51099895e-3 : 0.1056584;
   float leptonMass2 = (isMuon||isMuonE)?  0.1056584 : 0.51099895e-3;

   
   
   
   
   int charge1 = (int) l1.charge;
   int charge2 = (int) l2.charge;
   
   
   //TLorentzVector lepton1 (0.,0.,0.,0.);
   //TLorentzVector lepton2 (0.,0.,0.,0.);
   //TLorentzVector lepton1;
   //TLorentzVector lepton2;
   //lepton1.SetPtEtaPhiM(l1.p.Pt(),l1.p.Eta(),l1.p.Phi(),leptonMass1);
   //lepton2.SetPtEtaPhiM(l2.p.Pt(),l2.p.Eta(),l2.p.Phi(),leptonMass2);
   lep1.SetPtEtaPhiM(l1.p.Pt(),l1.p.Eta(),l1.p.Phi(),leptonMass1);
   lep2.SetPtEtaPhiM(l2.p.Pt(),l2.p.Eta(),l2.p.Phi(),leptonMass2);
   

   //pt1 = lepton1.Pt();  
   //pt2 = lepton2.Pt();  
   //phi1 = lepton1.Phi();  
   //phi2 = lepton2.Phi();  
   //eta1 = lepton1.Eta();  
   //eta2 = lepton2.Eta();
   pt1 = lep1.Pt();  
   pt2 = lep2.Pt();  
   phi1 = lep1.Phi();  
   phi2 = lep2.Phi();  
   eta1 = lep1.Eta();  
   eta2 = lep2.Eta();
   
   
   

   chargeProduct = charge1*charge2;
   //lep1=lepton1;
   //lep2=lepton2;
   //mll = (lepton1+lepton2).M();
   mll = (lep1+lep2).M();
   miniIso1 = l1.miniIso;
   miniIso2 = l2.miniIso;
   //Delta R between all three partciles
   //deltaRll=lepton1.DeltaR(lepton2);
   deltaRll=lep1.DeltaR(lep2);
   //MET
   ETmiss = mett.p.Pt();
   
   //TLorentzVector missing (0.,0.,0.,0.);
   //TLorentzVector missing;
   //missing.SetPtEtaPhiM(mett.p.Pt(),mett.p.Eta(),mett.p.Phi(),0.);
   //ETmiss_vec = missing;
   //TLorentzVector missing;
   //missing.
   ETmiss_vec.SetPtEtaPhiM(mett.p.Pt(),mett.p.Eta(),mett.p.Phi(),0.);;
   

   //double pa[3] = {lep1.M(),lep1.Px(),lep1.Py()};
   //double pb[3] = {lep2.M(),lep2.Px(),lep2.Py()};
   //double pmiss[3] = {0.,missing.Px(),missing.Py()};
   //fctMT2_.set_mn(0.);
   //fctMT2_.set_momenta(pa,pb,pmiss,*evtNo);
   //MT2_val = static_cast<float>(fctMT2_.get_mt2()); 
   
   evtHasGenPhotonVeto = GenPhotonVeto(6);
   

}

void HistogramProducer::ClearVariables(){

   isDiElectron=false;
   isDiMuon=false;
   isMuonElectron=false;
   isElectronMuon=false;
   float leptonMass1 = 0.;
   float leptonMass2 = 0.;
   float qter = 1.0;
   
   trigDiEle=false;
   trigDiMu=false;
   trigMuEle=false;
   trigHt=false;
   
   int charge1 = 0;
   int charge2 = 0;
   
   //TLorentzVector lepton1 (0.,0.,0.,0.);
   //TLorentzVector lepton2 (0.,0.,0.,0.);
   lep1.SetPtEtaPhiM(0.,0.,0.,0.);
   lep2.SetPtEtaPhiM(0.,0.,0.,0.);

   pt1 = 0.;  
   pt2 = 0.;  
   phi1 = 0.;  
   phi2 = 0.;  
   eta1 = 0.;  
   eta2 = 0.;

   miniIso1=1000.;
   miniIso2=1000.;

   chargeProduct = charge1*charge2;

   //lep1=lepton1;
   //lep2=lepton2;
   //mll = (lepton1+lepton2).M();
   mll = 0.;
   deltaRll=0.;
   //MET
   ETmiss = 0.;
   
   //TLorentzVector missing (0.,0.,0.,0.);
   //ETmiss_vec = missing;
   ETmiss_vec.SetPtEtaPhiM(0.,0.,0.,0.);
   
//;
   //double pa[3] = {lep1.M(),lep1.Px(),lep1.Py()};
   //double pb[3] = {lep2.M(),lep2.Px(),lep2.Py()};
  //
   //double pmiss[3] = {0.,missing.Px(),missing.Py()};
   //fctMT2_.set_mn(0.);
   //fctMT2_.set_momenta(pa,pb,pmiss,*evtNo);
   //MT2_val = static_cast<float>(fctMT2_.get_mt2()); 
   
   evtHasGenPhotonVeto = false;   
   
   selEvent L;
   L.isDiElectron=isDiElectron;
   L.isDiMuon=isDiMuon;
   L.isMuonElectron=isMuonElectron;
   L.isElectronMuon=isElectronMuon;
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
   //L.MT2_val = MT2_val;
   L.evtHasGenPhotonVeto = evtHasGenPhotonVeto;
   

   selectedEvent=L;
   
   }

bool HistogramProducer::CheckParticles(){
   return ((muons->size()>=2.) || (electrons->size()>=2.) || ((muons->size()>=1)&&(electrons->size()>=1)) );      
   //return ((myMuons.size()>=2.) || (electrons->size()>=2.) || ((myMuons.size()>=1)&&(electrons->size()>=1)) );      
}
bool HistogramProducer::Check2Ele(){
   return (electrons->size()>=2.);
}
bool HistogramProducer::Check2Mu(){
   return (muons->size()>=2.);
   //return (myMuons.size()>=2.);
}
bool HistogramProducer::CheckEMu(){
   return ((muons->size()>=1.) && (electrons->size()>=1.));
   //return ((myMuons.size()>=1.) && (electrons->size()>=1.));
}


void HistogramProducer::CorrectAllMuonPt(){
   float qter=1.;
   if(isData){
      for(auto& it : *muons){
         //cout<<it.p.Pt()<<endl;
         tree::Muon tempMuon;

         tempMuon.p.SetPtEtaPhi(it.p.Pt(),it.p.Eta(),it.p.Phi());
         tempMuon.charge=it.charge; // +/- 1
         tempMuon.rIso=it.rIso;
         tempMuon.passImpactParameter=it.passImpactParameter;
         tempMuon.d0=it.d0;
         tempMuon.dZ=it.dZ;
         tempMuon.SIP3D=it.SIP3D;
         tempMuon.miniIso=it.miniIso;         
         tempMuon.isTight=it.isTight;
         tempMuon.isMedium=it.isMedium;
         tempMuon.isMediumRun=it.isMediumRun;
         tempMuon.nTrkLayers=it.nTrkLayers;
      
         
         
         //int chargeTemp = (int) tempMuon.charge;
         
         //TLorentzVector leptonDummy (0.,0.,0.,0.);
         TLorentzVector leptonDummy;
         float isotemp=tempMuon.miniIso;
         float ptOld=tempMuon.p.Pt();
         leptonDummy.SetPtEtaPhiM(tempMuon.p.Pt(),tempMuon.p.Eta(),tempMuon.p.Phi(),0.1056584);
         //rmcor.momcor_data(leptonDummy, (int) tempMuon.charge, 0, qter);
         tempMuon.p.SetPtEtaPhi(leptonDummy.Pt(),leptonDummy.Eta(),leptonDummy.Phi());
         tempMuon.miniIso = isotemp*ptOld/leptonDummy.Pt();
         //myMuons.push_back(tempMuon);
         //cout<<tempMuon.p.Pt()<<endl;
      }
   }
   if(!isData){
      for(auto& it : *muons){
         //cout<<it.p.Pt()<<endl;
         //TLorentzVector leptonDummy (0.,0.,0.,0.);
         TLorentzVector leptonDummy;
         float isotemp=it.miniIso;
         float ptOld=it.p.Pt();
         tree::Muon tempMuon;

         tempMuon.p.SetPtEtaPhi(it.p.Pt(),it.p.Eta(),it.p.Phi());
         tempMuon.charge=it.charge; // +/- 1
         tempMuon.rIso=it.rIso;
         tempMuon.passImpactParameter=it.passImpactParameter;
         tempMuon.d0=it.d0;
         tempMuon.dZ=it.dZ;
         tempMuon.SIP3D=it.SIP3D;
         tempMuon.miniIso=it.miniIso;         
         tempMuon.isTight=it.isTight;
         tempMuon.isMedium=it.isMedium;
         tempMuon.isMediumRun=it.isMediumRun;
         tempMuon.nTrkLayers=it.nTrkLayers;
      
         
         
         //int chargeTemp = (int) tempMuon.charge;
         
         leptonDummy.SetPtEtaPhiM(tempMuon.p.Pt(),tempMuon.p.Eta(),tempMuon.p.Phi(),0.1056584);
         //rmcor.momcor_mc(leptonDummy, (int) tempMuon.charge, tempMuon.nTrkLayers, qter);
         tempMuon.p.SetPtEtaPhi(leptonDummy.Pt(),leptonDummy.Eta(),leptonDummy.Phi());
         tempMuon.miniIso = isotemp*ptOld/leptonDummy.Pt();
         //myMuons.push_back(tempMuon);
         //cout<<tempMuon.p.Pt()<<endl;
      }
   }
}


bool HistogramProducer::Cleaning() //return if event is valid (2 Leptons exist)
{
   if(!isSignal){
      trigDiEle = *hlt_ele17_ele12_iso || *hlt_ele23_ele12_iso || *hlt_doubleEle33 || *hlt_doubleEle33_mw;
      trigDiMu  = *hlt_mu17_mu8_iso || *hlt_mu17_tkMu8_iso || *hlt_mu17_mu8_iso_dz || *hlt_mu17_tkMu8_iso_dz
                  || *hlt_tkMu17_tkMu8_iso_dz || *hlt_mu27_tkMu8 || *hlt_mu30_tkMu11;
      trigMuEle = *hlt_mu17_ele12_iso || *hlt_mu23_ele8_iso || *hlt_mu23_ele8_iso_dz || *hlt_mu23_ele12_iso
                  || *hlt_mu23_ele12_iso_dz || *hlt_mu8_ele17_iso || *hlt_mu8_ele23_iso
                  || *hlt_mu8_ele23_iso_dz || *hlt_mu12_ele23_iso || *hlt_mu12_ele23_iso_dz
                  || *hlt_mu30_ele30 || *hlt_mu33_ele33;
      trigHt    = *hlt_ht200 || *hlt_ht250 || *hlt_ht300 || *hlt_ht350 || *hlt_ht400
                  || *hlt_ht475 || *hlt_ht600 || *hlt_ht650 || *hlt_ht800;
   }else{
      trigDiEle = true;
      trigDiMu = true;
      //trigMuEle = false;
      trigMuEle = true;
      trigHt = true;
   }
   

   if (isData){
      bool isEESelection = inputName.find("DoubleEG") != string::npos;
      bool isMuMuSelection = inputName.find("DoubleMuon") != string::npos;
      bool isEMuSelection = inputName.find("MuonEG") != string::npos;
      bool isHTSelection = inputName.find("JetHT")!= string::npos;
      auto missingET = *met;
        
      //if (isEESelection && trigDiEle && Check2Ele()) {
         //auto e1 = electrons->at(0);
         //auto e2 = electrons->at(1);
         //if (trigDiMu && Check2Mu()) {
            //auto m1 = muons->at(0);
            //auto m2 = muons->at(1);
            //if (min(e1.p.Pt(), e2.p.Pt()) > min(m1.p.Pt(), m2.p.Pt())){
               //CalculateVariables(e1, e2, missingET,E);
               //return true;
            //}else{return false;}
         //}else{
            //CalculateVariables(e1, e2, missingET,E);
            //return true;
         //}
      //}
      //if (isMuMuSelection && trigDiMu && Check2Mu()) {
         //auto m1 = muons->at(0);
         //auto m2 = muons->at(1);
         //if (trigDiEle && Check2Ele()) {
            //auto e1 = electrons->at(0);
            //auto e2 = electrons->at(1);
            //if (min(m1.p.Pt(), m2.p.Pt()) > min(e1.p.Pt(), e2.p.Pt())){
               //CalculateVariables(m1, m2, missingET,M);
               //return true;
            //}else{return false;}
         //}else{
            //CalculateVariables(m1, m2, missingET,M);
            //return true;
         //}
      //}
      
      //if(isMuMuSelection && Check2Mu() && trigDiMu){
         //auto m1 = muons->at(0);
         //auto m2 = muons->at(1);
         //CalculateVariables(m1,m2,missingET,M);
         //return true;
      //}
      //if(isEESelection && Check2Ele() && !trigDiMu && trigDiEle){
         //auto e1 = electrons->at(0);
         //auto e2 = electrons->at(1);
         //CalculateVariables(e1,e2,missingET,E);
         //return true;
      //}
      //if(isEMuSelection && CheckEMu() && !trigDiMu && !trigDiEle && trigMuEle){
         //auto e0 = electrons->at(0);
         //auto m0 = muons->at(0);
         //if (e0.p.Pt()>m0.p.Pt()){
            //CalculateVariables(e0,m0,missingET,EM);
         //}else{
            //CalculateVariables(m0,e0,missingET,ME);
         //}
         //return true;
      //}
      
      //third approach
      //if(isMuMuSelection && Check2Mu() && trigDiMu){
         //if (!(trigDiEle || trigMuEle)){
            //auto m1 = muons->at(0);
            //auto m2 = muons->at(1);
            //CalculateVariables(m1,m2,missingET,M);
            //return true;
         //}
      //}
      //if(isEESelection && Check2Ele() && trigDiEle){
         //if (!(trigDiMu || trigMuEle)){
            //auto e1 = electrons->at(0);
            //auto e2 = electrons->at(1);
            //CalculateVariables(e1,e2,missingET,E);
            //return true;
         //}
      //}
      //if(isEMuSelection && CheckEMu() && trigMuEle){
         //if (!(trigDiMu || trigDiEle)){
            //auto e0 = electrons->at(0);
            //auto m0 = muons->at(0);
            //if(e0.p.Pt() > m0.p.Pt()){
               //CalculateVariables(e0,m0,missingET,EM);
            //}else{
               //CalculateVariables(m0,e0,missingET,ME);
            //}
            //return true;
         //}
      //}
      //cout<<"here"<<endl;
      //->only combinations leftover
      //-->look for leading and trailing lepton, indepentdently of flavor
      //--->put in this category
      //---->only CalculateVariables() if one is in the right Dataset
      
      float pt1_ = 0.;
      float pt2_ = 0.;
      float tol = 0.0001;
  
      std::string leptonFlavor1 = "N";
      std::string leptonFlavor2 = "N";
      
     for(auto& it : *electrons){
        if ((it.p.Pt()-tol) > pt1_){
           pt1_ = it.p.Pt();
           leptonFlavor1 = "E";
        }
      }
      for(auto& it : *muons){
      //for(auto& it : myMuons){
        if ((it.p.Pt()-tol) > pt1_){
           pt1_ = it.p.Pt();
           leptonFlavor1 = "M";
        }
      }
      
      for(auto& it : *electrons){
        if ((it.p.Pt() < (pt1_-tol)) && ((it.p.Pt()-tol) > pt2_)){
           pt2_ = it.p.Pt();
           leptonFlavor2 = "E";
        }
      }
      for(auto& it : *muons){
      //for(auto& it : myMuons){
        if ((it.p.Pt() < (pt1_-tol)) && ((it.p.Pt()-tol) > pt2_)){
           pt2_ = it.p.Pt();
           leptonFlavor2 = "M";
        }
      }
      if (leptonFlavor1 == "E" && leptonFlavor2 == "E") {
         if(isEESelection && Check2Ele() && trigDiEle){ CalculateVariables(electrons->at(0),electrons->at(1),missingET,E); return true;}
      }
      else if (leptonFlavor1 == "M" && leptonFlavor2 == "M") {
         if(isMuMuSelection && Check2Mu() && trigDiMu){ CalculateVariables(muons->at(0),muons->at(1),missingET,M); return true;}      
      }
      else if (leptonFlavor1 == "E" && leptonFlavor2 == "M") {
         if(isEMuSelection && CheckEMu() && trigMuEle){ CalculateVariables(electrons->at(0),muons->at(0),missingET,EM); return true;}      
      }
      else if (leptonFlavor1 == "M" && leptonFlavor2 == "E") {
         if(isEMuSelection && CheckEMu() && trigMuEle){ CalculateVariables(muons->at(0),electrons->at(0),missingET,ME); return true;}      
      }
      
      
      //fourth approach
      //if (isEESelection){
         //if(trigDiEle){
            //if(Check2Ele()){
               //CalculateVariables(electrons->at(0),electrons->at(1),missingET,E);
               //return true;
            //}
         //}
      //}
      //if (isMuMuSelection){
         //if(trigDiMu && !trigDiEle){
            //if(Check2Mu()){
               //CalculateVariables(muons->at(0),muons->at(1),missingET,M);
               //return true;
            //}
         //}
      //}
      //if (isEMuSelection){
         //if(trigMuEle && !trigDiEle && !trigDiMu){
         //if(trigMuEle){
            //if(CheckEMu()){
               //if(muons->at(0).p.Pt() > electrons->at(0).p.Pt()){
                  //CalculateVariables(muons->at(0),electrons->at(0),missingET,ME);
                  //return true;                  
               //}else{
                  //CalculateVariables(electrons->at(0),muons->at(0),missingET,EM);
                  //return true;
               //}
            //}
         //}
      //}
      
      
      
      return false;
      }else{
         auto missingET = *met;
      //
         //if (trigDiEle && Check2Ele()) {
            //auto e1 = electrons->at(0);
            //auto e2 = electrons->at(1);
            //if (trigDiMu && Check2Mu()) {
               //auto m1 = muons->at(0);
               //auto m2 = muons->at(1);
               //if (min(e1.p.Pt(), e2.p.Pt()) > min(m1.p.Pt(), m2.p.Pt())){
                  //CalculateVariables(e1, e2, missingET,E);
                  //return true;
               //}else{return false;}
            //}else{
            //CalculateVariables(e1, e2, missingET,E);
            //return true;
            //}
         //}
         //if (trigDiMu && Check2Mu()) {
            //auto m1 = muons->at(0);
            //auto m2 = muons->at(1);
            //if (trigDiEle && Check2Ele()) {
               //auto e1 = electrons->at(0);
               //auto e2 = electrons->at(1);
               //if (min(m1.p.Pt(), m2.p.Pt()) > min(e1.p.Pt(), e2.p.Pt())){
                  //CalculateVariables(m1, m2, missingET,M);
                  //return true;
               //}else{return false;}
            //}else{
            //CalculateVariables(m1, m2, missingET,M);
            //return true;
            //}
         //}
         //if(Check2Mu() && (isSignal||trigDiMu)){
            //auto m1 = muons->at(0);
            //auto m2 = muons->at(1);
            //CalculateVariables(m1,m2,missingET,M);
            //return true;
         //}
         //if(Check2Ele() && (isSignal || !trigDiMu) && (isSignal||trigDiEle)){
            //auto e1 = electrons->at(0);
            //auto e2 = electrons->at(1);
            //CalculateVariables(e1,e2,missingET,E);
            //return true;
         //}
         //if(CheckEMu() && (isSignal || !trigDiMu) && (isSignal || !trigDiEle) && (isSignal || trigMuEle)){
            //auto e0 = electrons->at(0);
            //auto m0 = muons->at(0);
            //CalculateVariables(e0,m0,missingET,EM);
            //return true;
         //}
         
         //third approach
         //cout<<"1"<<endl;
         //if(Check2Mu() && trigDiMu){
            //if (!(trigDiEle || trigMuEle)){
               //auto m1 = muons->at(0);
               //auto m2 = muons->at(1);
               //CalculateVariables(m1,m2,missingET,M);
               //return true;
            //}
         //}
         //cout<<"2"<<endl;
         //if(Check2Ele() && trigDiEle){
            //if (!(trigDiMu || trigMuEle)){
               //auto e1 = electrons->at(0);
               //auto e2 = electrons->at(1);
               //CalculateVariables(e1,e2,missingET,E);
               //return true;
            //}
         //}
         //cout<<"3"<<endl;
         //if(CheckEMu() && trigMuEle){
            //if (!(trigDiMu || trigDiEle)){
               //auto e0 = electrons->at(0);
               //auto m0 = muons->at(0);
               //if(e0.p.Pt() > m0.p.Pt()){
                  //CalculateVariables(e0,m0,missingET,EM);
               //}else{
                  //CalculateVariables(m0,e0,missingET,ME);
               //}
               //return true;
            //}
         //}
         //cout<<"4"<<endl;
         //->only combinations leftover or signal....
         //-->look for leading and trailing lepton, indepentdently of flavor
         //--->put in this category
         //---->only CalculateVariables() if one is in the right Dataset
         float pt1_ = 0.;
         float pt2_ = 0.;
         float tol = 0.0001;
     
         std::string leptonFlavor1 = "N";
         std::string leptonFlavor2 = "N";
         
        for(auto& it : *electrons){
           if ((it.p.Pt()-tol) > pt1_){
              pt1_ = it.p.Pt();
              leptonFlavor1 = "E";
           }
         }
         for(auto& it : *muons){
         //for(auto& it : myMuons){
           if ((it.p.Pt()-tol) > pt1_){
              pt1_ = it.p.Pt();
              leptonFlavor1 = "M";
           }
         }
         for(auto& it : *electrons){
           if ((it.p.Pt() < (pt1_-tol)) && ((it.p.Pt()-tol) > pt2_)){
              pt2_ = it.p.Pt();
              leptonFlavor2 = "E";
           }
         }
         for(auto& it : *muons){
         //for(auto& it : myMuons){
           if ((it.p.Pt() < (pt1_-tol)) && ((it.p.Pt()-tol) > pt2_)){
              pt2_ = it.p.Pt();
              leptonFlavor2 = "M";
           }
         }
         if (leptonFlavor1 == "E" && leptonFlavor2 == "E") {
            if(Check2Ele() && trigDiEle){CalculateVariables(electrons->at(0),electrons->at(1),missingET,E); return true;}
         }
         else if (leptonFlavor1 == "M" && leptonFlavor2 == "M") {
            if(Check2Mu() && trigDiMu){ CalculateVariables(muons->at(0),muons->at(1),missingET,M); return true;}
         }
         else if (leptonFlavor1 == "E" && leptonFlavor2 == "M") {
            if(CheckEMu() && trigMuEle){ CalculateVariables(electrons->at(0),muons->at(0),missingET,EM); return true;}     
         }
         else if (leptonFlavor1 == "M" && leptonFlavor2 == "E") {
            if(CheckEMu() && trigMuEle){ CalculateVariables(muons->at(0),electrons->at(0),missingET,ME); return true;}     
         }
         
               //fourth approach
            //if(trigMuEle){
               //if(CheckEMu()){
                  //if(muons->at(0).p.Pt() > electrons->at(0).p.Pt()){
                     //CalculateVariables(muons->at(0),electrons->at(0),missingET,ME);
                     //return true;                  
                  //}else{
                     //CalculateVariables(electrons->at(0),muons->at(0),missingET,EM);
                     //return true;
                  //}
               //}
            //}
            //
            //if(trigDiEle){
               //if(Check2Ele()){
                  //CalculateVariables(electrons->at(0),electrons->at(1),missingET,E);
                  //return true;
               //}
            //}
         

            //if(trigDiMu && !trigDiEle){
            //if(trigDiMu){
               //if(Check2Mu()){
                  //CalculateVariables(muons->at(0),muons->at(1),missingET,M);
                  //return true;
               //}
            //}
         

            //if(trigMuEle && !trigDiEle && !trigDiMu){
            //if(trigMuEle){
               //if(CheckEMu()){
                  //if(muons->at(0).p.Pt() > electrons->at(0).p.Pt()){
                     //CalculateVariables(muons->at(0),electrons->at(0),missingET,ME);
                     //return true;                  
                  //}else{
                     //CalculateVariables(electrons->at(0),muons->at(0),missingET,EM);
                     //return true;
                  //}
               //}
            //}
         
         
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
      trigMuEle = *hlt_mu17_ele12_iso || *hlt_mu23_ele8_iso || *hlt_mu23_ele8_iso_dz || *hlt_mu23_ele12_iso
                  || *hlt_mu23_ele12_iso_dz || *hlt_mu8_ele17_iso || *hlt_mu8_ele23_iso
                  || *hlt_mu8_ele23_iso_dz || *hlt_mu12_ele23_iso || *hlt_mu12_ele23_iso_dz
                  || *hlt_mu30_ele30 || *hlt_mu33_ele33;
      //trigMuEle = false;
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
        
      //if (isHTSelection && trigHt && Check2Ele()) {
         //auto e1 = electrons->at(0);
         //auto e2 = electrons->at(1);
         //if (Check2Mu()) {
            //auto m1 = muons->at(0);
            //auto m2 = muons->at(1);
            //if (min(e1.p.Pt(), e2.p.Pt()) > min(m1.p.Pt(), m2.p.Pt())){
               //CalculateVariables(e1, e2, missingET,E);
               //return true;
            //}else{return false;}
         //}else{
            //CalculateVariables(e1, e2, missingET,E);
            //return true;
         //}
      //}
      //if (isHTSelection && trigHt && Check2Mu()) {
         //auto m1 = muons->at(0);
         //auto m2 = muons->at(1);
         //if (Check2Ele()) {
            //auto e1 = electrons->at(0);
            //auto e2 = electrons->at(1);
            //if (min(m1.p.Pt(), m2.p.Pt()) > min(e1.p.Pt(), e2.p.Pt())){
               //CalculateVariables(m1, m2, missingET,M);
               //return true;
            //}else{return false;}
         //}else{
            //CalculateVariables(m1, m2, missingET,M);
            //return true;
         //}
      //}
      //return false;

      float pt1_ = 0.;
      float pt2_ = 0.;
      float tol = 0.0001;
  
      std::string leptonFlavor1 = "N";
      std::string leptonFlavor2 = "N";
      
     for(auto& it : *electrons){
        if ((it.p.Pt()-tol) > pt1_){
           pt1_ = it.p.Pt();
           leptonFlavor1 = "E";
        }
      }
      for(auto& it : *muons){
      //for(auto& it : myMuons){
        if ((it.p.Pt()-tol) > pt1_){
           pt1_ = it.p.Pt();
           leptonFlavor1 = "M";
        }
      }
      for(auto& it : *electrons){
        if ((it.p.Pt() < (pt1_-tol)) && ((it.p.Pt()-tol) > pt2_)){
           pt2_ = it.p.Pt();
           leptonFlavor2 = "E";
        }
      }
      for(auto& it : *muons){
      //for(auto& it : myMuons){
        if ((it.p.Pt() < (pt1_-tol)) && ((it.p.Pt()-tol) > pt2_)){
           pt2_ = it.p.Pt();
           leptonFlavor2 = "M";
        }
      }
      

      
      if (leptonFlavor1 == "E" && leptonFlavor2 == "E") {
         if(isHTSelection && Check2Ele() && trigHt){ CalculateVariables(electrons->at(0),electrons->at(1),missingET,E); return true;}
      }
      else if (leptonFlavor1 == "M" && leptonFlavor2 == "M") {
         if(isHTSelection && Check2Mu() && trigHt){ CalculateVariables(muons->at(0),muons->at(1),missingET,M); return true;}      
      }
      else if (leptonFlavor1 == "E" && leptonFlavor2 == "M") {
         if(isHTSelection && CheckEMu() && trigHt){ CalculateVariables(electrons->at(0),muons->at(0),missingET,EM); return true;}      
      }
      else if (leptonFlavor1 == "M" && leptonFlavor2 == "E") {
         if(isHTSelection && CheckEMu() && trigHt){ CalculateVariables(muons->at(0),electrons->at(0),missingET,ME); return true;}      
      }
      return false;      

   }else{
      auto missingET = *met;
   
      //if (trigHt && Check2Ele()) {
         //auto e1 = electrons->at(0);
         //auto e2 = electrons->at(1);
         //if (Check2Mu()) {
            //auto m1 = muons->at(0);
            //auto m2 = muons->at(1);
            //if (min(e1.p.Pt(), e2.p.Pt()) > min(m1.p.Pt(), m2.p.Pt())){
               //CalculateVariables(e1, e2, missingET,E);
               //return true;
            //}else{return false;}
         //}else{
         //CalculateVariables(e1, e2, missingET,E);
         //return true;
         //}
      //}
      //if (trigHt && Check2Mu()) {
         //auto m1 = muons->at(0);
         //auto m2 = muons->at(1);
         //if (Check2Ele()) {
            //auto e1 = electrons->at(0);
            //auto e2 = electrons->at(1);
            //if (min(m1.p.Pt(), m2.p.Pt()) > min(e1.p.Pt(), e2.p.Pt())){
               //CalculateVariables(m1, m2, missingET,M);
               //return true;
            //}else{return false;}
         //}else{
         //CalculateVariables(m1, m2, missingET,M);
         //return true;
         //}
      //}
      
         float pt1_ = 0.;
         float pt2_ = 0.;
         float tol = 0.0001;
     
         std::string leptonFlavor1 = "N";
         std::string leptonFlavor2 = "N";
         
        for(auto& it : *electrons){
           if ((it.p.Pt()-tol) > pt1_){
              pt1_ = it.p.Pt();
              leptonFlavor1 = "E";
           }
         }
         for(auto& it : *muons){
         //for(auto& it : myMuons){
           if ((it.p.Pt()-tol) > pt1_){
              pt1_ = it.p.Pt();
              leptonFlavor1 = "M";
           }
         }
         for(auto& it : *electrons){
           if ((it.p.Pt() < (pt1_-tol)) && ((it.p.Pt()-tol) > pt2_)){
              pt2_ = it.p.Pt();
              leptonFlavor2 = "E";
           }
         }
         for(auto& it : *muons){
         //for(auto& it : myMuons){
           if ((it.p.Pt() < (pt1_-tol)) && ((it.p.Pt()-tol) > pt2_)){
              pt2_ = it.p.Pt();
              leptonFlavor2 = "M";
           }
         }
         if (leptonFlavor1 == "E" && leptonFlavor2 == "E") {
            if(Check2Ele() && trigHt){CalculateVariables(electrons->at(0),electrons->at(1),missingET,E); return true;}
         }
         else if (leptonFlavor1 == "M" && leptonFlavor2 == "M") {
            if(Check2Mu() && trigHt){ CalculateVariables(muons->at(0),muons->at(1),missingET,M); return true;}
            //if(Check2Mu() && trigHt){ CalculateVariables(myMuons.at(0),myMuons.at(1),missingET,M); return true;}
         }
         else if (leptonFlavor1 == "E" && leptonFlavor2 == "M") {
            if(CheckEMu() && trigHt){ CalculateVariables(electrons->at(0),muons->at(0),missingET,EM); return true;}     
            //if(CheckEMu() && trigHt){ CalculateVariables(electrons->at(0),myMuons.at(0),missingET,EM); return true;}     
         }
         else if (leptonFlavor1 == "M" && leptonFlavor2 == "E") {
            if(CheckEMu() && trigHt){ CalculateVariables(muons->at(0),electrons->at(0),missingET,ME); return true;}     
            //if(CheckEMu() && trigHt){ CalculateVariables(myMuons.at(0),electrons->at(0),missingET,ME); return true;}     
         }
      
   return false;
   }
   return false;
}


void HistogramProducer::clearCutFlowMap(){
   decisionMapCutFlowFine[TRIGGERED]=false;
   decisionMapCutFlowFine[genZLL]=false;
   //decisionMapCutFlowFine[LEPTONID_leading]=false;
   //decisionMapCutFlowFine[LEPTONID_trailing]=false;
   decisionMapCutFlowFine[LEPTONIDPure_leading]=false;
   decisionMapCutFlowFine[LEPTONIDPure_trailing]=false;
   decisionMapCutFlowFine[LEPTONIDImpact_leading]=false;
   decisionMapCutFlowFine[LEPTONIDImpact_trailing]=false;
   decisionMapCutFlowFine[LEPTONIDIso_leading]=false;
   decisionMapCutFlowFine[LEPTONIDIso_trailing]=false;
   decisionMapCutFlowFine[LEPTONIDDeltaR_leading]=false;
   decisionMapCutFlowFine[LEPTONIDDeltaR_trailing]=false;
   decisionMapCutFlowFine[LEPTONIDEta_leading]=false;
   decisionMapCutFlowFine[LEPTONIDEta_trailing]=false;
   decisionMapCutFlowFine[LEPTONPT_leading]=false;
   decisionMapCutFlowFine[LEPTONPT_trailing]=false;
   decisionMapCutFlowFine[PHOTON1]=false;
   decisionMapCutFlowFine[PHOTON1ID]=false;
   decisionMapCutFlowFine[PHOTON1SEED]=false;
   decisionMapCutFlowFine[PHOTON1ETA]=false;
   decisionMapCutFlowFine[PHOTON1PT]=false;
   decisionMapCutFlowFine[PHOTON1DR]=false;
   decisionMapCutFlowFine[M50]=false;
   decisionMapCutFlowFine[DIELECTRON]=false;
   decisionMapCutFlowFine[DIMUON]=false;
   decisionMapCutFlowFine[EMUON]=false;
   decisionMapCutFlowFine[ZMASS]=false;
   decisionMapCutFlowFine[GENVETO]=false;
   
   decisionMapCutFlowFine_weight[TRIGGERED]=totalWeight;
   decisionMapCutFlowFine_weight[genZLL]=totalWeight;
   //decisionMapCutFlowFine_weight[LEPTONID_leading]=totalWeight;
   //decisionMapCutFlowFine_weight[LEPTONID_trailing]=totalWeight;
   decisionMapCutFlowFine_weight[LEPTONIDPure_leading]=totalWeight;
   decisionMapCutFlowFine_weight[LEPTONIDPure_trailing]=totalWeight;
   decisionMapCutFlowFine_weight[LEPTONIDImpact_leading]=totalWeight;
   decisionMapCutFlowFine_weight[LEPTONIDImpact_trailing]=totalWeight;
   decisionMapCutFlowFine_weight[LEPTONIDEta_leading]=totalWeight;
   decisionMapCutFlowFine_weight[LEPTONIDEta_trailing]=totalWeight;
   decisionMapCutFlowFine_weight[LEPTONIDIso_leading]=totalWeight;
   decisionMapCutFlowFine_weight[LEPTONIDIso_trailing]=totalWeight;
   decisionMapCutFlowFine_weight[LEPTONIDDeltaR_leading]=totalWeight;
   decisionMapCutFlowFine_weight[LEPTONIDDeltaR_trailing]=totalWeight;
   decisionMapCutFlowFine_weight[LEPTONPT_leading]=totalWeight;
   decisionMapCutFlowFine_weight[LEPTONPT_trailing]=totalWeight;
   decisionMapCutFlowFine_weight[PHOTON1]=totalWeight;
   decisionMapCutFlowFine_weight[PHOTON1ID]=totalWeight;
   decisionMapCutFlowFine_weight[PHOTON1SEED]=totalWeight;
   decisionMapCutFlowFine_weight[PHOTON1ETA]=totalWeight;
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
   }
   if(selection==SEL){
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVA;
      decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVA;
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && pa.isTightMVA;
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVASlope;
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (deltaRll>0.1) && ElectronTightMVA(pa.p.Eta(),pa.p.Pt(),pa.mvaValue);
   
   }
   if(selection==LooseLeptons){
      decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && pa.isLoose;
   }
   if(selection==ONZ || selection==ControlRegionZZ || selection==ControlRegionWZ){
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVA;
      decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVA;
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && pa.isTightMVA;
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (deltaRll>0.1) && ElectronTightMVA(pa.p.Eta(),pa.p.Pt(),pa.mvaValue);
      //if (pa.isPassConvVeto && pa.isTightMVA){
      if (pa.isPassConvVeto && pa.isTightMVASlope){
         leading? decisionMapCutFlowFine[LEPTONIDPure_leading]=true : decisionMapCutFlowFine[LEPTONIDPure_trailing]=true;
         if (pa.passImpactParameter){
            leading? decisionMapCutFlowFine[LEPTONIDImpact_leading]=true : decisionMapCutFlowFine[LEPTONIDImpact_trailing]=true;
            //if((fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6)){
            if((fabs(pa.p.Eta())<2.4)){
               leading? decisionMapCutFlowFine[LEPTONIDEta_leading]=true : decisionMapCutFlowFine[LEPTONIDEta_trailing]=true;
               if(pa.miniIso<0.1){
                  leading? decisionMapCutFlowFine[LEPTONIDIso_leading]=true : decisionMapCutFlowFine[LEPTONIDIso_trailing]=true;
                  if(deltaRll>0.1){
                  //if(deltaRll>0.){
                     leading? decisionMapCutFlowFine[LEPTONIDDeltaR_leading]=true : decisionMapCutFlowFine[LEPTONIDDeltaR_trailing]=true;
                     if((leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.))){
                        leading? decisionMapCutFlowFine[LEPTONPT_leading]=true : decisionMapCutFlowFine[LEPTONPT_trailing]=true;
                     }
                  }
               }
            }
         }
      }
   
   }
   if(selection==DILEP){
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVA ;
      decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4)  && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVA ;
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4)  && (pa.miniIso<0.1) && pa.isTightMVA ;
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4)  && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVASlope ;
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4)  && (pa.miniIso<0.1) && (deltaRll>0.1) && ElectronTightMVA(pa.p.Eta(),pa.p.Pt(),pa.mvaValue) ;
   }
   if(selection==TRIGSEL || selection==TRIGDILEP || selection==TRIGONZ){
      decision =(*ht>200.) && pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>0.) : (pa.p.Pt()>0.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVA;
      //decision =(*ht>200.) && pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>0.) : (pa.p.Pt()>0.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVASlope;
      //decision =(*ht>200.) && pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>0.) : (pa.p.Pt()>0.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && ElectronTightMVA(pa.p.Eta(),pa.p.Pt(),pa.mvaValue);
   }
   if(selection==TRIGSEL_ptcuts || selection==TRIGDILEP_ptcuts || selection==TRIGONZ_ptcuts){
      decision =(*ht>200.) && pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVA;
      //decision =(*ht>200.) && pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && pa.isTightMVASlope;
      //decision =(*ht>200.) && pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (deltaRll>0.1) && ElectronTightMVA(pa.p.Eta(),pa.p.Pt(),pa.mvaValue);
   }
   if(selection==EXO){
      decision = pa.isPassConvVeto && pa.passImpactParameter && (leading? (fabs(pa.p.Eta())<1.4442 ? (pa.p.Pt()>65.) : (pa.p.Pt()>70.)) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.5) && (pa.miniIso<0.1) && pa.isTightMVA;
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading? (fabs(pa.p.Eta())<1.4442 ? (pa.p.Pt()>65.) : (pa.p.Pt()>70.)) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.5) && (pa.miniIso<0.1) && ElectronTightMVA(pa.p.Eta(),pa.p.Pt(),pa.mvaValue);
   }
   if(selection==EGRegression){
      decision = (leading? (pa.p.Pt()>32.) : (pa.p.Pt()>25.)) && (fabs(pa.p.Eta())<2.5) && pa.isLoose;
   }
   return decision;
}
//MUON
bool HistogramProducer::testSelection(const tree::Muon& pa, selectionType selection_, bool leading){
   bool decision=false;
   if ((selection_==UNCUT)||(selection_==PHOTON)){
      decision = pa.passImpactParameter  && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (pa.isTight);
   }
   if (selection_==SEL){
      //decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (deltaRll>0.1);
      decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.2) && (pa.isMedium) && (deltaRll>0.1);
      //decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.2) && (pa.isMedium) ;
   }
   if (selection_==LooseLeptons){
      decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.2) ;
   }
   if (selection_ == ONZ || selection_==DILEP || selection_==ControlRegionZZ || selection_==ControlRegionWZ){
      //decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (deltaRll>0.1);
      decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.2) && (pa.isMedium) && (deltaRll>0.1);
      //decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.2) && (pa.isMedium) ;
      if((pa.isMedium)){
         leading? decisionMapCutFlowFine[LEPTONIDPure_leading]=true : decisionMapCutFlowFine[LEPTONIDPure_trailing]=true;
         if(pa.passImpactParameter){
            leading? decisionMapCutFlowFine[LEPTONIDImpact_leading]=true : decisionMapCutFlowFine[LEPTONIDImpact_trailing]=true;
            //if((fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6)){
            if((fabs(pa.p.Eta())<2.4) ){
               leading? decisionMapCutFlowFine[LEPTONIDEta_leading]=true : decisionMapCutFlowFine[LEPTONIDEta_trailing]=true;
               if(pa.miniIso<0.2){
                  leading? decisionMapCutFlowFine[LEPTONIDIso_leading]=true : decisionMapCutFlowFine[LEPTONIDIso_trailing]=true;
                  if((deltaRll>0.1)){
                  //if((deltaRll>0.)){
                     leading? decisionMapCutFlowFine[LEPTONIDDeltaR_leading]=true : decisionMapCutFlowFine[LEPTONIDDeltaR_trailing]=true;
                     if((leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.))){
                        leading? decisionMapCutFlowFine[LEPTONPT_leading]=true : decisionMapCutFlowFine[LEPTONPT_trailing]=true;
                     }
                  }
               }
            }
         }
      }
   }
   //if(selection==DILEP || selection==EGRegression){
   //if(selection_==DILEP){
      //decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (deltaRll>0.1);
      //if (pa.passImpactParameter && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (deltaRll>0.1)){
         //leading? decisionMapCutFlowFine[LEPTONID_leading]=true : decisionMapCutFlowFine[LEPTONID_trailing]=true;
         //if((leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.))){
            //leading? decisionMapCutFlowFine[LEPTONPT_leading]=true : decisionMapCutFlowFine[LEPTONPT_trailing]=true;
         //}
      //}
   //}
   if(selection_==TRIGSEL || selection_==TRIGDILEP || selection_==TRIGONZ){
      decision = (*ht>200.) && pa.passImpactParameter && (leading ? (pa.p.Pt()>0.) : (pa.p.Pt()>0.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (deltaRll>0.1);
   }
   if(selection_==TRIGSEL_ptcuts || selection_==TRIGDILEP_ptcuts || selection_==TRIGONZ_ptcuts){
      decision = (*ht>200.)&& pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (deltaRll>0.1);
   }
   if(selection_==TRIGDILEP_pt1cut){
      decision = (*ht>200.)&& pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>0.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (deltaRll>0.1);
   }
   if(selection_==TRIGDILEP_pt2cut){
      decision = (*ht>200.)&& pa.passImpactParameter && (leading ? (pa.p.Pt()>0.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (deltaRll>0.1);
   }
   if(selection_==EXO){
      decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>52.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (pa.isMedium);
   }
   return decision;
}
//PHOTON
bool HistogramProducer::testSelection(const selPhoton& pa, selectionType selection){
   bool decision=false;
   if ((selection==UNCUT)||(selection==PHOTON)){
      decision = !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isLoose) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3);
      //decision = !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isLoose) ;
   }
   if(selection==SEL || selection==DILEP ||selection==ONZ || selection==EGRegression || selection==ControlRegionZZ){
      decision = (pa.p.Pt()>25.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isLoose) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3); //study Delta R cut ;
      //decision = (pa.p.Pt()>25.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isLoose) ; //study Delta R cut ;
      if((pa.isLoose)){
         if(selection==ONZ) decisionMapCutFlowFine[PHOTON1ID]=true;
         if(!(pa.hasPixelSeed)){
            if(selection==ONZ) decisionMapCutFlowFine[PHOTON1SEED]=true;
            if(fabs(pa.p.Eta()<1.4442)){
               if(selection==ONZ) decisionMapCutFlowFine[PHOTON1ETA]=true;
               if((pa.p.Pt()>25.)){
                  if(selection==ONZ) decisionMapCutFlowFine[PHOTON1PT]=true;
                  if((pa.deltaR1>0.3) && (pa.deltaR2>0.3)){
                  //if((pa.deltaR1>0.) && (pa.deltaR2>0.)){
                     if(selection==ONZ) decisionMapCutFlowFine[PHOTON1DR]=true;
                  }
               }
            }
         }
      }
   }
   if(selection==EXO){
      decision = (isDiElectron? (pa.p.Pt()>(65. + 5.* ((1.566 < fabs(pa.p.Eta())) && (fabs(pa.p.Eta()) < 2.5))  )) : (pa.p.Pt()>40.)) && (pa.passElectronVeto) && (fabs(pa.p.Eta())<2.5) && (pa.isLoose) && (pa.deltaR1>0.4) && (pa.deltaR2>0.4); //study Delta R cut ;
   }
   if(selection==TRIGDILEP || selection==TRIGSEL || selection==TRIGONZ){
      decision = (pa.p.Pt()>25.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isLoose) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3); //study Delta R cut ;
   }
   return decision;
}
//SELMUON
bool HistogramProducer::testSelection(const selMuon& pa, selectionType selection_){
   bool decision=false;
   if (selection_==LooseLeptons){
      //decision = pa.passImpactParameter && (pa.p.Pt()>20.) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.2) ;
      //decision = pa.passImpactParameter && (pa.p.Pt()>20.) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.4) ;
      decision = pa.passImpactParameter && (pa.p.Pt()>20.) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.2) && (pa.isMedium);
      //decision = pa.passImpactParameter && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.2) ;
   }
   return decision;
}
//SELELECTRON
bool HistogramProducer::testSelection(const selElectron& pa, selectionType selection_){
   bool decision=false;
   if(selection_==LooseLeptons){
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (pa.p.Pt()>20.) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && pa.isLoose;
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (pa.p.Pt()>20.) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.4) && pa.isLoose;
      decision = pa.isPassConvVeto && pa.passImpactParameter && (pa.p.Pt()>20.) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && pa.isTightMVA;
      //decision = pa.isPassConvVeto && pa.passImpactParameter && (fabs(pa.p.Eta())<2.4) && pa.isLoose;
   }
   return decision;
}
//JET
bool HistogramProducer::testSelection(const selJet& pa, selectionType selection){
   bool decision=false;
   if ((selection==UNCUT)||(selection==PHOTON)){
      decision = (fabs(pa.p.Eta())<2.4) && (pa.isLoose) && !(pa.hasElectronMatch) && !(pa.hasPhotonMatch) && !(pa.hasMuonMatch);
   }
   if(selection==SEL||selection==DILEP || selection==ONZ || selection==TRIGDILEP || selection==TRIGSEL || selection==TRIGONZ || selection==EXO || selection==EGRegression|| selection==ControlRegionZZ){
      decision = (pa.p.Pt()>30.) && (fabs(pa.p.Eta())<2.4) && (pa.isLoose) && !(pa.hasElectronMatch) && !(pa.hasPhotonMatch) && !(pa.hasMuonMatch); //&& (pa.deltaR1>0.3) && (pa.deltaR2>0.3) //study Delta R cut ;
   }
   return decision;
}



void HistogramProducer::Init(TTree *tree)
{
   fReader.SetTree(tree);
   inputName = fReader.GetTree()->GetCurrentFile()->GetName();
   
   isTotalSignal = (inputName.find("SMS-TChiNG_BF") != string::npos) || (inputName.find("SMS-T6ttZg_nTuple")!=string::npos) || (inputName.find("SMS-T5bbbbZg_nTuple")!=string::npos) || (inputName.find("GGM")!=string::npos);
   isData = inputName.find("Run201") != string::npos;
   isSignal = (inputName.find("SMS") != string::npos) || (inputName.find("GGM") != string::npos);
   
   boost::property_tree::ini_parser::read_ini("example.ini", propertyTree);


   doWeights_TopPt = (inputName.find("TTTo") != string::npos) || (inputName.find("TTGamma") != string::npos) || (inputName.find("ST_") != string::npos);
   //doWeights_TopPt = false;
   doWeights_EWKinoPairPt = (inputName.find("TChi") != string::npos);
   //doWeights_LeptonPairPt = (inputName.find("ZG") != string::npos) || (inputName.find("DY") != string::npos);
   doWeights_LeptonPairPt = false;
   doWeights_nISR = (inputName.find("T5") != string::npos);
   //doWeights_nISR = (inputName.find("TTTo") != string::npos) || (inputName.find("TTGamma") != string::npos);
   //doWeights_nISR = false;


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
   config_selectionsToProcessMap[ControlRegionDY] = propertyTree.get<bool>("selectionsToProcess.control");
   config_selectionsToProcessMap[ControlRegionTT] = propertyTree.get<bool>("selectionsToProcess.control");
   config_selectionsToProcessMap[ControlRegionZZ] = propertyTree.get<bool>("selectionsToProcess.control");
   config_selectionsToProcessMap[ControlRegionWZ] = propertyTree.get<bool>("selectionsToProcess.control");
   config_selectionsToProcessMap[ControlRegionWW] = propertyTree.get<bool>("selectionsToProcess.control");
   config_selectionsToProcessMap[ValidationRegion] = propertyTree.get<bool>("selectionsToProcess.validation");
   
   config_veto = propertyTree.get<int>("generalSettings.veto");
   config_docutflow = propertyTree.get<bool>("generalSettings.docutflow");
   config_docutflowfine = propertyTree.get<bool>("generalSettings.docutflowfine");
   config_dosignalscan = propertyTree.get<bool>("generalSettings.dosignalscan");
   config_dosignalscanSplit = propertyTree.get<bool>("generalSettings.dosignalscantchingsplit");
   config_dosignalscanSplit = config_dosignalscanSplit && config_dosignalscan && inputName.find("TChiNG") != string::npos;
   //config_dosignalscanSplit = false;
   config_eventpercentage = propertyTree.get<float>("generalSettings.eventpercentage");
   config_outputfolder = propertyTree.get<string>("generalSettings.outputfolder");
   
   
   if(isTotalSignal){
      config_selectionsToProcessMap[UNCUT] = false;
      config_selectionsToProcessMap[PHOTON] = false;
      config_selectionsToProcessMap[ONZ] = false;
      config_selectionsToProcessMap[DILEP] = false;
      config_selectionsToProcessMap[ONZMET] = false;
      config_selectionsToProcessMap[EGRegression] = false;
      config_selectionsToProcessMap[EXO] = false;
      config_selectionsToProcessMap[SEL] = false;
      config_selectionsToProcessMap[ABOVEZG] = false;
      config_selectionsToProcessMap[ONZG] = false;
      config_selectionsToProcessMap[ControlRegionDY]=false;
      config_selectionsToProcessMap[ControlRegionTT]=false;
      config_selectionsToProcessMap[ControlRegionZZ]=false;
      config_selectionsToProcessMap[ControlRegionWZ]=false;
      config_selectionsToProcessMap[ControlRegionWW]=false;
      config_selectionsToProcessMap[ValidationRegion]=false;
      config_docutflow = false;
      config_docutflowfine = false;
   }
   
   if(! isTotalSignal){
      cutFlow = *((TH1F*)fReader.GetTree()->GetCurrentFile()->Get("TreeWriter/hCutFlow"));
   }
   fReader.GetEntries(true); // jumps to last file
   //string lastInputName = fReader.GetTree()->GetCurrentFile()->GetName();
  
   nEntries=fReader.GetTree()->GetEntries();

  
  
  //rochester muon pt corrections
   //rochcor2016 *rmcor = new rochcor2016();
  
  
  setHistoNames();
  
   //if (inputName!=lastInputName) {
       ////////adds cut flow of last file. This makes only sense if there is for 1 or two files
      //cutFlow.Add((TH1F*)fReader.GetTree()->GetCurrentFile()->Get("TreeWriter/hCutFlow"));
   //}
   //float nGen = cutFlow.GetBinContent(2);


   //noPromptPhotons = inputName.find("DYJets") != string::npos || inputName.find("TTTo2L2Nu") != string::npos;
   noPromptPhotons = inputName.find("DYJets") != string::npos;
   isZGammaInclusive = (inputName.find("ZGTo2LG_ext") != string::npos) || (inputName.find("ZGTo2LG_Total") != string::npos);

   if (config_docutflow) InitCutFlowHistos();
   if (config_docutflowfine) InitCutFlowHistos_Fine();
   if(!isTotalSignal){
      InitAllHistos();
      InitTriggerStudiesHistos();
   }
   //InitScaleFactors();
   //InitScaleFactorsAlternative();
   InitScaleFactorsFinal();
   
   //InitCompressedTree();
   
   dummy1=0;
   dummy2=0;
   dummy3=0;
   dummy4=0;
   dummy5=0;
   dummy6=0;
   dummy7=0;
   dummy8=0;
   dummy9=0;
   dummy10=0;
   dummy11=0;
   dummy12=0;
   dummy13=0;
   dummy14=0;
   dummy15=0;
   dummy16=0;

}


void HistogramProducer::InitCompressedTree(){

   newtreeCompressed = new TTree("compressedTree","NNData");

   
   newtreeCompressed->Branch("weight",&comp_weight);
   
   newtreeCompressed->Branch("Lepton1Pt",&comp_Lepton1Pt);
   newtreeCompressed->Branch("Lepton1Px",&comp_Lepton1Px);
   newtreeCompressed->Branch("Lepton1Py",&comp_Lepton1Py);
   newtreeCompressed->Branch("Lepton1Pz",&comp_Lepton1Py);
   newtreeCompressed->Branch("Lepton1Eta",&comp_Lepton1Eta);
   newtreeCompressed->Branch("Lepton1Phi",&comp_Lepton1Phi);
   newtreeCompressed->Branch("Lepton1MiniIso",&comp_Lepton1MiniIso);
   newtreeCompressed->Branch("Lepton1IsElectron",&comp_Lepton1IsElectron);
   newtreeCompressed->Branch("Lepton1IsMuon",&comp_Lepton1IsMuon);
   
   newtreeCompressed->Branch("Lepton2Pt",&comp_Lepton2Pt);
   newtreeCompressed->Branch("Lepton2Px",&comp_Lepton2Px);
   newtreeCompressed->Branch("Lepton2Py",&comp_Lepton2Py);
   newtreeCompressed->Branch("Lepton2Pz",&comp_Lepton2Py);
   newtreeCompressed->Branch("Lepton2Eta",&comp_Lepton2Eta);
   newtreeCompressed->Branch("Lepton2Phi",&comp_Lepton2Phi);
   newtreeCompressed->Branch("Lepton2MiniIso",&comp_Lepton2MiniIso);
   newtreeCompressed->Branch("Lepton2IsElectron",&comp_Lepton2IsElectron);
   newtreeCompressed->Branch("Lepton2IsMuon",&comp_Lepton2IsMuon);
   
   newtreeCompressed->Branch("mll",&comp_mll);
   newtreeCompressed->Branch("deltaRll",&comp_deltaRll);

   newtreeCompressed->Branch("NPhoton",&comp_NPhoton);   
   newtreeCompressed->Branch("PhotonPt",&comp_PhotonPt);
   newtreeCompressed->Branch("PhotonPx",&comp_PhotonPx);
   newtreeCompressed->Branch("PhotonPy",&comp_PhotonPx);
   newtreeCompressed->Branch("PhotonPz",&comp_PhotonPz);
   newtreeCompressed->Branch("PhotonEta",&comp_PhotonEta);
   newtreeCompressed->Branch("PhotonPhi",&comp_PhotonPhi);
   
   newtreeCompressed->Branch("Photon1Pt",&comp_Photon1Pt);
   newtreeCompressed->Branch("Photon1Px",&comp_Photon1Px);
   newtreeCompressed->Branch("Photon1Py",&comp_Photon1Py);
   newtreeCompressed->Branch("Photon1Pz",&comp_Photon1Pz);
   newtreeCompressed->Branch("Photon1Eta",&comp_Photon1Eta);
   newtreeCompressed->Branch("Photon1Phi",&comp_Photon1Phi);
   
   
   newtreeCompressed->Branch("NJet",&comp_NJet);
   newtreeCompressed->Branch("JetPt",&comp_JetPt);
   newtreeCompressed->Branch("JetPx",&comp_JetPx);
   newtreeCompressed->Branch("JetPy",&comp_JetPx);
   newtreeCompressed->Branch("JetPz",&comp_JetPz);
   newtreeCompressed->Branch("JetEta",&comp_JetEta);
   newtreeCompressed->Branch("JetPhi",&comp_JetPhi);
   newtreeCompressed->Branch("JetBDiscriminator",&comp_JetBDiscriminator);
   newtreeCompressed->Branch("JetChf",&comp_JetChf);
   newtreeCompressed->Branch("JetNhf",&comp_JetNhf);
   newtreeCompressed->Branch("JetNConstituents",&comp_JetNConstituents);
   
   newtreeCompressed->Branch("MetPt",&comp_MetPt);
   newtreeCompressed->Branch("MetPx",&comp_MetPx);
   newtreeCompressed->Branch("MetPy",&comp_MetPy);
   newtreeCompressed->Branch("MetPz",&comp_MetPz);
   newtreeCompressed->Branch("MetEta",&comp_MetEta);
   newtreeCompressed->Branch("MetPhi",&comp_MetPhi);
   
   newtreeCompressed->Branch("ht",&comp_ht);
   
   newtreeCompressed->Branch("m1",&comp_m1);
   newtreeCompressed->Branch("m2",&comp_m2);
   
   newtreeCompressed->Branch("eventNo",&comp_eventNo);
   
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
   
   PhotonIDWeighter = Weighter("scaleFactors/ScaleFactorsPhotonsID_LooseCut.root","EGamma_SF2D"); // pt vs eta
   PhotonIDWeighter.fillOverflow2d();
   
   
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
   
   PhotonConversionWeighter = Weighter("scaleFactors/ScaleFactorPhotonConversion.root","Scaling_Factors_HasPix_R9 Inclusive"); // pt vs eta
   PhotonConversionWeighter.fillOverflow2d();
   
   
   }
   
void HistogramProducer::InitScaleFactorsFinal(){
   //DiEleWeighterID=Weighter("newScaleFactors/electron/scaleFactor_ele_ID.root","GsfElectronToMVATightTightIP2DSIP3D4");// |eta| vs pt
   //DiEleWeighterID.fillOverflow2d();
   //DiEleWeighterIso=Weighter("newScaleFactors/electron/scaleFactor_ele_ID.root","MVAVLooseElectronToMini");// |eta| vs pt
   //DiEleWeighterIso.fillOverflow2d();
   //DiEleWeighterConv=Weighter("newScaleFactors/electron/scaleFactor_ele_ID.root","MVATightElectronToConvVetoIHit0"); // |eta| vs pt
   //DiEleWeighterConv.fillOverflow2d();
   //DiEleWeighterTrack=Weighter("newScaleFactors/electron/scaleFactor_ele_reco.root","EGamma_SF2D"); //pt vs |eta|
   //DiEleWeighterTrack.fillOverflow2d();
   
   DiEleWeighterID=Weighter("scaleFactors/Alternative/ScaleFactorElectronID_MVATight.root","EGamma_SF2D");// pt vs |eta|
   DiEleWeighterID.fillOverflow2d();
   DiEleWeighterTrack=Weighter("scaleFactors/Alternative/egammaEffi.txt_EGM2D.root","EGamma_SF2D"); //pt vs |eta|
   DiEleWeighterTrack.fillOverflow2d();  
   
   DiMuWeighterID = Weighter("newScaleFactors/muon/scaleFactor_muon_ID.root","SF");// |eta| vs pt
   DiMuWeighterID.fillOverflow2d();
   DiMuWeighterIso = Weighter("newScaleFactors/muon/scaleFactor_muon_iso.root","SF");// |eta| vs pt
   DiMuWeighterIso.fillOverflow2d();
   DiMuWeighterIP2D = Weighter("newScaleFactors/muon/scaleFactor_muon_IP2D.root","SF");// |eta| vs pt
   DiMuWeighterIP2D.fillOverflow2d();
   DiMuWeighterSIP3D = Weighter("newScaleFactors/muon/scaleFactor_muon_SIP3D.root","SF");// |eta| vs pt
   DiMuWeighterSIP3D.fillOverflow2d();
   DiMuWeighterTrack = Weighter("newScaleFactors/muon/TrackScaleFactorsMuons.root","muonTrackScaleFactorEtaHisto");// |eta|
   DiMuWeighterTrack.fillOverflow2d();
   
   //FastSimDiEleWeighterID=Weighter("newScaleFactors/electron/scaleFactor_ele_FastSim_ID.root","histo2D");// |eta| vs pt
   //FastSimDiEleWeighterID.fillOverflow2d();
   //FastSimDiEleWeighterIso=Weighter("newScaleFactors/electron/scaleFactor_ele_FastSim_iso.root","histo2D");// |eta| vs pt
   //FastSimDiEleWeighterIso.fillOverflow2d();
   //FastSimDiEleWeighterConv=Weighter("newScaleFactors/electron/scaleFactor_ele_FastSim_ConvVeto.root","histo2D");// |eta| vs pt
   //FastSimDiEleWeighterConv.fillOverflow2d();
   
   //FastSimDiMuWeighterID=Weighter("newScaleFactors/muon/scaleFactor_muon_FastSim_ID.root","histo2D");// |eta| vs pt
   //FastSimDiMuWeighterID.fillOverflow2d();
   //FastSimDiMuWeighterIso=Weighter("newScaleFactors/muon/scaleFactor_muon_FastSim_iso.root","histo2D");// |eta| vs pt
   //FastSimDiMuWeighterIso.fillOverflow2d();
   //FastSimDiMuWeighterIP2D=Weighter("newScaleFactors/muon/scaleFactor_muon_FastSim_IP2D.root","histo2D");// |eta| vs pt
   //FastSimDiMuWeighterIP2D.fillOverflow2d();
   //FastSimDiMuWeighterSIP3D=Weighter("newScaleFactors/muon/scaleFactor_muon_FastSim_SIP3D.root","histo2D");// |eta| vs pt
   //FastSimDiMuWeighterSIP3D.fillOverflow2d();
   
   //PhotonIDWeighter = Weighter("scaleFactors/ScaleFactorsPhotonsID.root","EGamma_SF2D"); // pt vs eta
   PhotonIDWeighter = Weighter("newScaleFactors/photon/scaleFactor_photon_ID.root","EGamma_SF2D"); // pt vs eta
   PhotonIDWeighter.fillOverflow2d();
   
   PhotonConversionWeighter = Weighter("newScaleFactors/photon/scaleFactor_photon_conv.root","Scaling_Factors_HasPix_R9 Inclusive"); // pt vs eta
   PhotonConversionWeighter.fillOverflow2d();
   
   
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


float HistogramProducer::GetScaleFactorAndErrorFinal(float pt, float eta,bool isFastSim,bool isEle){
   //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSLeptonSF#Data_leading_order_FullSim_M_AN1
   float DiEleSF =1.;
   float DiMuSF = 1.;
   float DiEleSFErr = 0.;
   float DiMuSFErr = 0.;
   //float tempPt,tempPtTrack,tempEta;
   float absEta=fabs(eta);

   //float FastSimDiEleSF = 1.;
   //float FastSimDiEleSFErr = 0.;
   //float FastSimDiMuSF = 1.;
   //float FastSimDiMuSFErr = 0.;

   if(isEle){
      //DiEleSF = DiEleSF * DiEleWeighterID.getWeight(pt,absEta);
      //DiEleSF = DiEleSF * DiEleWeighterIso.getWeight(pt,absEta);
      //DiEleSF = DiEleSF * DiEleWeighterConv.getWeight(pt,absEta);
      //DiEleSF = DiEleSF * DiEleWeighterTrack.getWeight(absEta,pt); 
      //
      //DiEleSFErr = DiEleWeighterID.getError(pt,absEta);
      //DiEleSFErr = pow((pow(DiEleSFErr,2.) + pow(DiEleWeighterIso.getError(pt,absEta),2.)),0.5);
      //DiEleSFErr = pow((pow(DiEleSFErr,2.) + pow(DiEleWeighterConv.getError(pt,absEta),2.)),0.5);
      //DiEleSFErr = pow((pow(DiEleSFErr,2.) + pow(DiEleWeighterTrack.getError(absEta,pt),2.)),0.5);
      //DiEleSFErr = pow((pow(DiEleSFErr,2.) + pow((0.01*DiEleSF*(pt<20. && pt>80.)),2.)),0.5); //additional 1% uncertainty for all (normal just pt<20 and >80)
      
      
      DiEleSF = DiEleSF * DiEleWeighterID.getWeight(absEta,pt);
      DiEleSF = DiEleSF * DiEleWeighterTrack.getWeight(absEta,pt); 
      
      //if(isFastSim){
         //FastSimDiEleSF = DiEleSF * FastSimDiEleWeighterID.getWeight(pt,absEta);
         //FastSimDiEleSF = FastSimDiEleSF * FastSimDiEleWeighterIso.getWeight(pt,absEta);
         //FastSimDiEleSF = FastSimDiEleSF * FastSimDiEleWeighterConv.getWeight(pt,absEta);
//
         //FastSimDiEleSFErr = pow((pow(DiEleSFErr,2.) + pow(FastSimDiEleWeighterID.getError(pt,absEta),2.)),0.5);
         //FastSimDiEleSFErr = pow((pow(FastSimDiEleSFErr,2.) + pow(FastSimDiEleWeighterIso.getError(pt,absEta),2.)),0.5);
         //FastSimDiEleSFErr = pow((pow(FastSimDiEleSFErr,2.) + pow(FastSimDiEleWeighterConv.getError(absEta,pt),2.)),0.5);
         //FastSimDiEleSFErr = pow((pow(FastSimDiEleSFErr,2.) + pow((0.02*FastSimDiEleSF),2.)),0.5); //additional 2% uncertainty
   
      //}
   }else{
      DiMuSF = DiMuSF * DiMuWeighterID.getWeight(pt,absEta);
      DiMuSF = DiMuSF * DiMuWeighterIso.getWeight(pt,absEta);
      DiMuSF = DiMuSF * DiMuWeighterIP2D.getWeight(pt,absEta);
      //DiMuSF = DiMuSF * DiMuWeighterSIP3D.getWeight(absEta,pt); 
      DiMuSF = DiMuSF * DiMuWeighterSIP3D.getWeight(pt,absEta); 
      DiMuSF = DiMuSF * DiMuWeighterTrack.getWeight(absEta); 
      
      //DiMuSFErr = DiMuWeighterID.getError(pt,absEta);
      //DiMuSFErr = pow((pow(DiMuSFErr,2.) + pow(DiMuWeighterIso.getError(pt,absEta),2.)),0.5);
      //DiMuSFErr = pow((pow(DiMuSFErr,2.) + pow(DiMuWeighterIP2D.getError(pt,absEta),2.)),0.5);
      //DiMuSFErr = pow((pow(DiMuSFErr,2.) + pow(DiMuWeighterSIP3D.getError(pt,absEta),2.)),0.5);
      //DiMuSFErr = pow((pow(DiMuSFErr,2.) + pow(DiMuWeighterTrack.getError(absEta),2.)),0.5);
      //DiMuSFErr = pow((pow(DiMuSFErr,2.) + pow((0.03*DiMuSF),2.)),0.5); //additional 3% 
      
      //if(isFastSim){
         //FastSimDiMuSF = DiMuSF * FastSimDiMuWeighterID.getWeight(pt,absEta);
         //FastSimDiMuSF = FastSimDiMuSF * FastSimDiMuWeighterIso.getWeight(pt,absEta);
         //FastSimDiMuSF = FastSimDiMuSF * FastSimDiMuWeighterIP2D.getWeight(pt,absEta);
         //FastSimDiMuSF = FastSimDiMuSF * FastSimDiMuWeighterSIP3D.getWeight(pt,absEta);
//
         //FastSimDiMuSFErr = pow((pow(DiMuSFErr,2.) + pow(FastSimDiMuWeighterID.getError(pt,absEta),2.)),0.5);
         //FastSimDiMuSFErr = pow((pow(FastSimDiMuSFErr,2.) + pow(FastSimDiMuWeighterIso.getError(pt,absEta),2.)),0.5);
         //FastSimDiMuSFErr = pow((pow(FastSimDiMuSFErr,2.) + pow(FastSimDiMuWeighterIP2D.getError(absEta,pt),2.)),0.5);   
         //FastSimDiMuSFErr = pow((pow(FastSimDiMuSFErr,2.) + pow(FastSimDiMuWeighterSIP3D.getError(absEta,pt),2.)),0.5);
         //FastSimDiMuSFErr = pow((pow(FastSimDiMuSFErr,2.) + pow((0.02*FastSimDiMuSF),2.)),0.5); //additional 2% uncertainty
      //}
   }
   float SF=1.;
   //if(isFastSim){
      //SF = isEle ?  FastSimDiEleSF : FastSimDiMuSF;
   //}else{
      SF = isEle ?  DiEleSF : DiMuSF;
   //}
   return SF;
}




float HistogramProducer::GetScaleFactorAndErrorAlternative(float pt, float eta,bool isFastSim,bool isEle,int runNr){
   float DiEleSF =1.;
   float DiMuSF_BCDEF = 1.;
   float DiMuSF_GH = 1.;
   float DiEleSFErr = 0.;
   float DiMuSFErr = 0.;
   float tempPt,tempPtTrack,tempEta;
   float absEta=fabs(eta);

   //float DiMuAverage=1.;
   //float DiMuAverageErr=0.;

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
      
      //DiMuAverage = 
      
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
      totalSF = totalSF * PhotonIDWeighter.getWeight(fabs(tempEta),tempPt)*PhotonConversionWeighter.getWeight(fabs(tempEta),tempPt);
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
   hMap[ZPT] = TH1F("", ";Z_{p_T}", 5000, 0, 5000);
   hMap[ST] = TH1F("", ";S_T", 5000, 0, 5000.);
   hMap[MT2] = TH1F("", ";M_{T2}", 50000, 0, 5000.);
   hMap[VetoCompare] = TH1F("", "", 2, 0, 2);
   hMap[NElectrons] = TH1F("","nElectrons",20,0,20);
   hMap[NMuons] = TH1F("","nMuons",20,0,20);

   hMap[JetPt1] = TH1F("", ";#it{p}_{T}^{Jet 1} (GeV)", 5000, 0, 5000);
   hMap[JetPt2] = TH1F("", ";#it{p}_{T}^{Jet 2} (GeV)", 5000, 0, 5000);
   hMap[JetPt3] = TH1F("", ";#it{p}_{T}^{Jet 3} (GeV)", 5000, 0, 5000);
   hMap[JetPt4] = TH1F("", ";#it{p}_{T}^{Jet 4} (GeV)", 5000, 0, 5000);
   hMap[JetPhi1] = TH1F("", ";|#phi_{Jet 1}|", 350, 0, 3.5);
   hMap[JetPhi2] = TH1F("", ";|#phi_{Jet 2}|", 350, 0, 3.5);
   hMap[JetPhi3] = TH1F("", ";|#phi_{Jet 3}|", 350, 0, 3.5);
   hMap[JetPhi4] = TH1F("", ";|#phi_{Jet 4}|", 350, 0, 3.5);
   hMap[JetEta1] = TH1F("", ";|#eta_{Jet 1}|", 260, 0, 2.6);
   hMap[JetEta2] = TH1F("", ";|#eta_{Jet 2}|", 260, 0, 2.6);
   hMap[JetEta3] = TH1F("", ";|#eta_{Jet 3}|", 260, 0, 2.6);
   hMap[JetEta4] = TH1F("", ";|#eta_{Jet 4}|", 260, 0, 2.6);

   hMap[WEIGHT_NISR] = TH1F("","weight",10000,0,10);
   hMap[WEIGHT_TOPPT] = TH1F("","weight",10000,0,10);
   hMap[WEIGHT_EWKINOPAIRPT] = TH1F("","weight",10000,0,10);
   hMap[WEIGHT_LEPTONPAIRPT] = TH1F("","weight",10000,0,10);




      hMap[DeltaEtaLL] = TH1F("", ";#Delta#Eta_{ll}", 6000, 0, 6.);
      hMap[DeltaPhiLL] = TH1F("", ";#Delta#Phi_{ll}", 6000, 0, 6.);
      hMap[DeltaRLL] = TH1F("", ";#DeltaR_{ll}", 6000, 0, 6.);

   hMap[DeltaEtaLL_neg] = TH1F("", ";#Delta#Eta_{ll}", 12000, -6., 6.);
   hMap[DeltaPhiLL_neg] = TH1F("", ";#Delta#Phi_{ll}", 12000, -6., 6.);
   hMap[DeltaRLL_neg] = TH1F("", ";#DeltaR_{ll}", 12000, -6., 6.);
         //
   if(selection_!=ControlRegionZZ){
   //hMap[DeltaEtaLL] = TH1F("", ";#Delta#Eta_{ll}", 6000, 0, 6.);
   //hMap[DeltaPhiLL] = TH1F("", ";#Delta#Phi_{ll}", 6000, 0, 6.);
   //hMap[DeltaRLL] = TH1F("", ";#DeltaR_{ll}", 6000, 0, 6.);
   //hMap[DeltaEtaLL_neg] = TH1F("", ";#Delta#Eta_{ll}", 12000, -6., 6.);
   //hMap[DeltaPhiLL_neg] = TH1F("", ";#Delta#Phi_{ll}", 12000, -6., 6.);
   //hMap[DeltaRLL_neg] = TH1F("", ";#DeltaR_{ll}", 12000, -6., 6.);
   hMap[DeltaPhiLLMet] = TH1F("", "#Delta#Phi_{ll,MET}", 6000, 0, 6.);
   hMap[DeltaEtaLLMet] = TH1F("", ";#Delta#Eta_{ll,MET}", 6000, 0, 6.);
   hMap[DeltaRLLMet] = TH1F("", "#DeltaR_{ll,MET}", 6000, 0, 6.);
   hMap[MTLL] = TH1F("", ";m_{T}^{ll}", 5000, 0, 5000); 
   }

   if(selection_==ControlRegionZZ){
      //hMap[DeltaEtaLL] = TH1F("", ";#Delta#Eta_{ll}", 6000, 0, 6.);
      //hMap[DeltaPhiLL] = TH1F("", ";#Delta#Phi_{ll}", 6000, 0, 6.);
      //hMap[DeltaRLL] = TH1F("", ";#DeltaR_{ll}", 6000, 0, 6.);
      
      hMap[PT3] = TH1F("", ";#it{p}_{T}^{trailing} (GeV)", 5000, 0, 5000);
      hMap[PT4] = TH1F("", ";#it{p}_{T}^{trailing} (GeV)", 5000, 0, 5000);
      hMap[ETA3] = TH1F("", ";|#eta_{leading}|", 260, 0, 2.6);
      hMap[ETA4] = TH1F("", ";|#eta_{leading}|", 260, 0, 2.6);
      hMap[PHI3] = TH1F("", ";|#phi_{leading}|", 350, 0, 3.5);
      hMap[PHI4] = TH1F("", ";|#phi_{leading}|", 350, 0, 3.5);
      hMap[ZPT2] = TH1F("", ";Z_{p_T}", 5000, 0, 5000);
      hMap[MLL2] = TH1F("", ";#it{m}_{ll} (GeV)", 5000, 0, 5000);
      hMap[MTL3MET] = TH1F("", ";#it{m}_{T}^{l2,Met} (GeV)", 5000, 0, 5000);
      hMap[MLLLL] = TH1F("", ";#it{m}_{llll} (GeV)", 5000, 0, 5000);
      
   }
   
   
      
   //if ((selection_==PHOTON)||(selection_==SEL)||(selection_==ONZ)||(selection_==EXO)||selection_==ControlRegionZZ){
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
   
   sMap[WEIGHT_NISR] = TH1F("","weight",10000,0,10);
   sMap[WEIGHT_TOPPT] = TH1F("","weight",10000,0,10);
   sMap[WEIGHT_EWKINOPAIRPT] = TH1F("","weight",10000,0,10);
   sMap[WEIGHT_LEPTONPAIRPT] = TH1F("","weight",10000,0,10);
   
   
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
   cMap[CUTFLOW_fine] = TH1F("","",17,0,17);
   cMap[CUTFLOW_fine].Fill("genZToLL",0);
   cMap[CUTFLOW_fine].Fill("triggered", 0);
   cMap[CUTFLOW_fine].Fill("GenPhotonVeto", 0);
   //cMap[CUTFLOW_fine].Fill("2LeptonID", 0);
   cMap[CUTFLOW_fine].Fill("2LeptonIDPure", 0);
   cMap[CUTFLOW_fine].Fill("2LeptonIDImpact", 0);
   cMap[CUTFLOW_fine].Fill("2LeptonIDEta", 0);
   cMap[CUTFLOW_fine].Fill("2LeptonIDIso", 0);
   cMap[CUTFLOW_fine].Fill("2LeptonIDDeltaR", 0);
   cMap[CUTFLOW_fine].Fill("2LeptonPT", 0);
   cMap[CUTFLOW_fine].Fill("m50", 0);
   cMap[CUTFLOW_fine].Fill("1Photon", 0);
   cMap[CUTFLOW_fine].Fill("1PhotonID", 0);
   cMap[CUTFLOW_fine].Fill("1PhotonSeed", 0);
   cMap[CUTFLOW_fine].Fill("1PhotonEta", 0);
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
   
   //int matchedEleCounter=0;
   //int matchedMuCounter=0;
   
   
   if (CheckParticles()){
      if (Cleaning()){//decide category in dilep events
         selEvent L;
         L.isDiElectron=isDiElectron;
         L.isDiMuon=isDiMuon;
         L.isMuonElectron=isMuonElectron;
         L.isElectronMuon=isElectronMuon;
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
         //L.MT2_val = MT2_val;
         L.evtHasGenPhotonVeto = evtHasGenPhotonVeto;
         
         //if(selection==ControlRegionZZ) cout<<"CR ZZ start"<<endl;
         
         if(!isData){
            if(L.isDiElectron){
               //L.totalWeight = totalWeight * GetScaleFactorAndErrorAlternative(pt1,eta1,isSignal,isDiElectron,*runNo)*GetScaleFactorAndErrorAlternative(pt2,eta2,isSignal,isDiElectron,*runNo);
               //L.totalWeight = totalWeight * GetScaleFactorAndErrorAlternative(L.pt1,L.eta1,isSignal,true,*runNo)*GetScaleFactorAndErrorAlternative(L.pt2,L.eta2,isSignal,true,*runNo);
               L.totalWeight = totalWeight * GetScaleFactorAndErrorFinal(L.pt1,L.eta1,isSignal,true)*GetScaleFactorAndErrorFinal(L.pt2,L.eta2,isSignal,true);
            }else{
               if(L.isDiMuon){
                  //L.totalWeight = totalWeight * GetScaleFactorAndErrorAlternative(pt1,eta1,isSignal,false,*runNo)*GetScaleFactorAndErrorAlternative(L.pt2,L.eta2,isSignal,false,*runNo);
                  L.totalWeight = totalWeight * GetScaleFactorAndErrorFinal(pt1,eta1,isSignal,false)*GetScaleFactorAndErrorFinal(L.pt2,L.eta2,isSignal,false);
               }
               else{
                  if(L.isMuonElectron||L.isElectronMuon){
                     //L.totalWeight = totalWeight * GetScaleFactorAndErrorAlternative(muons->at(0).p.Pt(),muons->at(0).p.Eta(),isSignal,false,*runNo)*GetScaleFactorAndErrorAlternative(electrons->at(0).p.Pt(),electrons->at(0).p.Eta(),isSignal,true,*runNo);
                     L.totalWeight = totalWeight * GetScaleFactorAndErrorFinal(muons->at(0).p.Pt(),muons->at(0).p.Eta(),isSignal,false)*GetScaleFactorAndErrorFinal(electrons->at(0).p.Pt(),electrons->at(0).p.Eta(),isSignal,true);
                     //L.totalWeight = totalWeight * GetScaleFactorAndErrorFinal(myMuons.at(0).p.Pt(),myMuons.at(0).p.Eta(),isSignal,false)*GetScaleFactorAndErrorFinal(electrons->at(0).p.Pt(),electrons->at(0).p.Eta(),isSignal,true);
                  }
               }
            }
         }else{
            L.totalWeight=totalWeight;
         }
         //fill cutflowInfo DiEle or DiMu?
         cutflowDiEle=isDiElectron;
         cutflowDiMu=isDiMuon;
         decisionMapCutFlowFine[DIMUON]=isDiMuon;
         decisionMapCutFlowFine[DIELECTRON]=isDiElectron;
         decisionMapCutFlowFine[EMUON]=isMuonElectron||isElectronMuon;
         
         int genNNToZ=0;
         int genNNToG=0;
         //int genNNToH=0;
         int genZToLL=0;
         
         if((selection==ONZ) && (!isTotalSignal)){
            //get "triggered" signal events for cutflow
            for(auto& intermediate: *intermediateGenParticles){
               for(auto& daughter: intermediate.daughters){
                  if((abs(daughter.pdgId)==23)&&(intermediate.pdgId>1000021)){
                     genNNToZ=genNNToZ+1;
                  }
                  if((abs(daughter.pdgId)==22)&&(intermediate.pdgId>1000021)){
                     genNNToG=genNNToG+1;
                  }
                  //if((abs(daughter.pdgId)==25)&&(intermediate.pdgId>1000021)){
                     //genNNToH=genNNToH+1;
                  //}
                  if(((abs(daughter.pdgId)==11)||(abs(daughter.pdgId)==13))&&(intermediate.pdgId==23)){
                     genZToLL=genZToLL+1;
                  }
               }
            }
         }
         

         
         if (L.isDiElectron){
            //if(selection==ONZ) decisionMapCutFlowFine[genZLL]=(counterGenEle>=2);
            if(selection==ONZ) decisionMapCutFlowFine_weight[genZLL]=totalWeight;
            
            //Fill cutflow triggered bool
            if(selection==ONZ) cutflowIsTriggered=true;
            if(selection==ONZ){
               if(isSignal){
                  decisionMapCutFlowFine[TRIGGERED]=((genNNToG==1)&&(genNNToZ==1)&&(genZToLL==2));
                  //cout<<"top!!"<<endl;
               }else{
                  decisionMapCutFlowFine[TRIGGERED]=trigDiEle;
               }
            }
            if(selection==ONZ) decisionMapCutFlowFine_weight[TRIGGERED]=totalWeight;
            
            auto e1 = electrons->at(0);
            auto e2 = electrons->at(1);
            
            
            
            if(!GenPhotonVeto(config_veto)){
               if(selection==ONZ) decisionMapCutFlowFine[GENVETO]=true;
               if(selection==ONZ) decisionMapCutFlowFine_weight[GENVETO]=totalWeight;
            
            
               if(testSelection(e1,selection,true) && testSelection(e2,selection,false) && (chargeProduct < 0.) ){

                  //Fill cutflow 2 leptons bool
                  if(selection==ONZ) cutflow2Leptons=true;
                  if((selection==UNCUT)? true : (L.mll>50.)){
                     //Fill cutflow mll>50 bool
                     if(selection==ONZ) cutflowMll50=true;
                     if(selection==ONZ) decisionMapCutFlowFine[M50]=true;
                     //if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONID_leading]=L.totalWeight;
                     //if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONID_trailing]=L.totalWeight;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDPure_leading]=L.totalWeight;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDPure_trailing]=L.totalWeight;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDImpact_leading]=L.totalWeight;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDImpact_trailing]=L.totalWeight;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDEta_leading]=L.totalWeight;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDEta_trailing]=L.totalWeight;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDIso_leading]=L.totalWeight;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDIso_trailing]=L.totalWeight;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDDeltaR_leading]=L.totalWeight;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDDeltaR_trailing]=L.totalWeight;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONPT_leading]=L.totalWeight;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONPT_trailing]=L.totalWeight;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[M50]=L.totalWeight;


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
                     for (vector<tree::Muon>::iterator it = muons->begin(); it != muons->end(); it++){
                     //for (vector<tree::Muon>::iterator it = myMuons.begin(); it != myMuons.end(); it++){
                        auto g = *it;
                        selMuon muon;
                        muon.setAll(g);
                        muon.vec.SetPtEtaPhiM(muon.p.Pt(),muon.p.Eta(),muon.p.Phi(),0.1056584);
                        if (testSelection(muon,LooseLeptons)){
                           muon.deltaR1 = muon.vec.DeltaR(L.l1);
                           muon.deltaR2 = muon.vec.DeltaR(L.l2);
                           muon.chargeInt = (int) g.charge;
                           if(matchLepton(muon,L)){
                              muon.matched=true;
                              L.matchedMuSize+=1;
                              L.matchedLeptonSize+=1;
                           }
                              L.selMuons.push_back(muon);
                        }
                     }
                     for (vector<tree::Electron>::iterator it = electrons->begin(); it != electrons->end(); it++){
                        auto g = *it;
                        selElectron ele;
                        ele.setAll(g);
                        ele.vec.SetPtEtaPhiM(ele.p.Pt(),ele.p.Eta(),ele.p.Phi(),0.51099895e-3);
                        if (testSelection(ele,LooseLeptons)){
                           ele.deltaR1 = ele.vec.DeltaR(L.l1);
                           ele.deltaR2 = ele.vec.DeltaR(L.l2);
                           ele.chargeInt = (int) g.charge;
                           if(matchLepton(ele,L)){
                              ele.matched=true;
                              L.matchedEleSize+=1;
                              L.matchedLeptonSize+=1;
                           }
                              L.selElectrons.push_back(ele);
                        }
                     }
                     
                     L.selLeptonSize=(L.selMuons.size()+L.selElectrons.size());
                     L.selMuonSize=(L.selMuons.size());
                     L.selElectronSize=(L.selElectrons.size());
                     
                     
                     //TLorentzVector temp4Vec(0.,0.,0.,0.);
                     //for (vector<selElectron>::iterator it = L.selElectrons.begin(); it != L.selElectrons.end(); it++){
                        //if(it->matched==false){
                           //temp4Vec=temp4Vec+it->vec;
                        //}
                     //}
                     //for (vector<selMuon>::iterator it = L.selMuons.begin(); it != L.selMuons.end(); it++){
                        //if(it->matched==false){
                           //temp4Vec=temp4Vec+it->vec;
                        //}
                     //}
                     //
                     //L.invAddLeptMass=temp4Vec.M();
                     
                     
                     if((selection==ONZ)&&(photons->size()!=0)) decisionMapCutFlowFine[PHOTON1]=true;
                     if((selection==ONZ)&&(photons->size()!=0)) decisionMapCutFlowFine_weight[PHOTON1]=L.totalWeight;
                     
                     
                     if(!isData){
                           L.totalWeight=L.totalWeight*GetScaleFactorAndErrorPhotons(L.selPhotons);
                     }
                     
                     if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1ID]=L.totalWeight;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1SEED]=L.totalWeight;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1ETA]=L.totalWeight;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1PT]=L.totalWeight;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1DR]=L.totalWeight;
                     
                     for (vector<tree::Jet>::iterator it = jets->begin(); it != jets->end(); it++){
                        auto j = *it;
                        selJet jet;
                        jet.setAll(j);
                        jet.vec.SetPtEtaPhiM(jet.p.Pt(),jet.p.Eta(),jet.p.Phi(),0.);
                        jet.deltaR1=jet.vec.DeltaR(L.l1);
                        jet.deltaR2=jet.vec.DeltaR(L.l2);
                        if (testSelection(jet,selection)){
                           L.selJets.push_back(jet);
                           L.calcHt += jet.vec.Pt();
                        }
                     }
                     
                     
                     selectedEvent=L;
                        if(selectedEvent.selPhotons.size()!=0){
                           if(selection==ONZ) cutflow1Photon=true;
                           if((L.mll>81.) && (L.mll<101.)){
                              if(selection==ONZ) cutflowOnZ=true;
                              if(selection==ONZ) decisionMapCutFlowFine[ZMASS]=true;
                              if(selection==ONZ) decisionMapCutFlowFine_weight[ZMASS]=L.totalWeight;
                           }
                        }
                     return true;
                  }
               }
            }
         }
            if(L.isDiMuon){

               //if(selection==ONZ) decisionMapCutFlowFine[genZLL]=(counterGenMu>=2);
               if(selection==ONZ) decisionMapCutFlowFine_weight[genZLL]=totalWeight;
               
               
               if(selection==ONZ) cutflowIsTriggered=true;
               if(selection==ONZ){
                  if(isSignal){
                     decisionMapCutFlowFine[TRIGGERED]=((genNNToG==1)&&(genNNToZ==1)&&(genZToLL==2));
                  }else{
                     decisionMapCutFlowFine[TRIGGERED]=trigDiMu;
                  }
               }
               //if(selection==ONZ) decisionMapCutFlowFine[TRIGGERED]=trigDiMu;
               if(selection==ONZ) decisionMapCutFlowFine_weight[TRIGGERED]=totalWeight;

               auto m1 = muons->at(0); 
               auto m2 = muons->at(1);
               //auto m1 = myMuons.at(0); 
               //auto m2 = myMuons.at(1);
               
               if(!GenPhotonVeto(config_veto)){
                  
                  if(selection==ONZ) decisionMapCutFlowFine[GENVETO]=true;
                  if(selection==ONZ) decisionMapCutFlowFine_weight[GENVETO]=totalWeight;


                  if(testSelection(m1,selection,true) && testSelection(m2,selection,false) && (chargeProduct < 0.)){

                     if(selection==ONZ) cutflow2Leptons=true;
                     if((selection==UNCUT)? true : (L.mll>50.)){
                        if(selection==ONZ) {cutflowMll50=true;}
                        if(selection==ONZ) decisionMapCutFlowFine[M50]=true;
                        //if(selection==ONZ) decisionMapCutFlowFine_weight[M50]=L.totalWeight;
                        //if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONID_leading]=L.totalWeight;
                        //if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONID_trailing]=L.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDPure_leading]=L.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDPure_trailing]=L.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDImpact_leading]=L.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDImpact_trailing]=L.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDEta_leading]=L.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDEta_trailing]=L.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDIso_leading]=L.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDIso_trailing]=L.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDDeltaR_leading]=L.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDDeltaR_trailing]=L.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONPT_leading]=L.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONPT_trailing]=L.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[M50]=L.totalWeight;
                        //if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONPT_leading]=L.totalWeight;
                        //if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONPT_trailing]=L.totalWeight;
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
                        for (vector<tree::Muon>::iterator it = muons->begin(); it != muons->end(); it++){
                        //for (vector<tree::Muon>::iterator it = myMuons.begin(); it != myMuons.end(); it++){
                           auto g = *it;
                           selMuon muon;
                           muon.setAll(g);
                           muon.vec.SetPtEtaPhiM(muon.p.Pt(),muon.p.Eta(),muon.p.Phi(),0.1056584);
                           if (testSelection(muon,LooseLeptons)){
                              muon.deltaR1 = muon.vec.DeltaR(L.l1);
                              muon.deltaR2 = muon.vec.DeltaR(L.l2);
                              muon.chargeInt = (int) g.charge;
                              if(matchLepton(muon,L)){
                                 muon.matched=true;
                                 L.matchedMuSize+=1;
                                 L.matchedLeptonSize+=1;
                              }
                                 L.selMuons.push_back(muon);
                           }
                        }
                        for (vector<tree::Electron>::iterator it = electrons->begin(); it != electrons->end(); it++){
                           auto g = *it;
                           selElectron ele;
                           ele.setAll(g);
                           ele.vec.SetPtEtaPhiM(ele.p.Pt(),ele.p.Eta(),ele.p.Phi(),0.51099895e-3);
                           if (testSelection(ele,LooseLeptons)){
                              ele.deltaR1 = ele.vec.DeltaR(L.l1);
                              ele.deltaR2 = ele.vec.DeltaR(L.l2);
                              ele.chargeInt = (int) g.charge;
                              if(matchLepton(ele,L)){
                                 ele.matched=true;
                                 L.matchedEleSize+=1;
                                 L.matchedLeptonSize+=1;
                              }
                                 L.selElectrons.push_back(ele);
                           }
                        }
                                            
                     //TLorentzVector temp4Vec(0.,0.,0.,0.);
                     //for (vector<selElectron>::iterator it = L.selElectrons.begin(); it != L.selElectrons.end(); it++){
                        //if(it->matched==false){
                           //temp4Vec=temp4Vec+it->vec;
                        //}
                     //}
                     //for (vector<selMuon>::iterator it = L.selMuons.begin(); it != L.selMuons.end(); it++){
                        //if(it->matched==false){
                           //temp4Vec=temp4Vec+it->vec;
                        //}
                     //}
                     //L.invAddLeptMass=temp4Vec.M();

                     
                     L.selLeptonSize=(L.selMuons.size()+L.selElectrons.size());
                     L.selMuonSize=(L.selMuons.size());
                     L.selElectronSize=(L.selElectrons.size());

                        
                        if(selection==ONZ) if(photons->size()!=0) decisionMapCutFlowFine[PHOTON1]=true;
                        if(selection==ONZ) if(photons->size()!=0) decisionMapCutFlowFine_weight[PHOTON1]=L.totalWeight;
                        if(!isData){
                              L.totalWeight=L.totalWeight*GetScaleFactorAndErrorPhotons(L.selPhotons);
                        }
                        
                        if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1ID]=L.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1SEED]=L.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1ETA]=L.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1PT]=L.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1DR]=L.totalWeight;
                        
                        for (vector<tree::Jet>::iterator it = jets->begin(); it != jets->end(); it++){
                           auto j = *it;
                           selJet jet;
                           jet.setAll(j);
                           jet.vec.SetPtEtaPhiM(jet.p.Pt(),jet.p.Eta(),jet.p.Phi(),0.);
                           jet.deltaR1=jet.vec.DeltaR(L.l1);
                           jet.deltaR2=jet.vec.DeltaR(L.l2);
                           if (testSelection(jet,selection)){
                              L.selJets.push_back(jet);
                              L.calcHt += jet.vec.Pt();
                           }
                        }
                     //for (auto& j : selJets) {
                        //L.calcHt += j->p.Pt();
                     //}                        
                        
                        selectedEvent=L;
                           
                        if(L.selPhotons.size()>0){
                           if(selection==ONZ) cutflow1Photon=true;
                           if(L.mll<101. && L.mll>81.){
                              if(selection==ONZ) cutflowOnZ=true;
                              if(selection==ONZ) decisionMapCutFlowFine[ZMASS]=true;
                              if(selection==ONZ) decisionMapCutFlowFine_weight[ZMASS]=L.totalWeight;
                           }
                        }
                        return true;
                     }
                  }
               }
            }
            //else{
               
               if(L.isElectronMuon || L.isMuonElectron){
               
                  if(selection==ONZ) cutflowIsTriggered=true;
                  if(selection==ONZ){
                     if(isSignal){
                        decisionMapCutFlowFine[TRIGGERED]=((genNNToG==1)&&(genNNToZ==1)&&(genZToLL==2));
                     }else{
                        decisionMapCutFlowFine[TRIGGERED]=trigMuEle;
                     }
                  }
                  //if(selection==ONZ) decisionMapCutFlowFine[TRIGGERED]=trigMuEle;
                  if(selection==ONZ) decisionMapCutFlowFine_weight[TRIGGERED]=totalWeight;

                  auto e0 = electrons->at(0); 
                  auto m0 = muons->at(0);
                  //auto m0 = myMuons.at(0);

                  
                  if(!GenPhotonVeto(config_veto)){
                     
                     if(selection==ONZ) decisionMapCutFlowFine[GENVETO]=true;
                     if(selection==ONZ) decisionMapCutFlowFine_weight[GENVETO]=totalWeight;


                     if(testSelection(m0,selection,L.isMuonElectron) && testSelection(e0,selection,L.isElectronMuon) && (chargeProduct < 0.)){
                        if(selection==ONZ) cutflow2Leptons=true;
                        if((selection==UNCUT)? true : (mll>50.)){
                           if(selection==ONZ) cutflowMll50=true;
                           if(selection==ONZ) decisionMapCutFlowFine[M50]=true;
                           if(selection==ONZ) decisionMapCutFlowFine_weight[M50]=L.totalWeight;
                           if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONID_leading]=L.totalWeight;
                           if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONID_trailing]=L.totalWeight;
                           if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONPT_leading]=L.totalWeight;
                           if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONPT_trailing]=L.totalWeight;

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
                           for (vector<tree::Muon>::iterator it = muons->begin(); it != muons->end(); it++){
                           //for (vector<tree::Muon>::iterator it = myMuons.begin(); it != myMuons.end(); it++){
                              auto g = *it;
                              selMuon muon;
                              muon.setAll(g);
                              muon.vec.SetPtEtaPhiM(muon.p.Pt(),muon.p.Eta(),muon.p.Phi(),0.1056584);
                              if (testSelection(muon,LooseLeptons)){
                                 muon.deltaR1 = muon.vec.DeltaR(L.l1);
                                 muon.deltaR2 = muon.vec.DeltaR(L.l2);
                                 muon.chargeInt = (int) g.charge;
                                 if(matchLepton(muon,L)){
                                    muon.matched=true;
                                    L.matchedMuSize+=1;
                                    L.matchedLeptonSize+=1;
                                 }
                                 L.selMuons.push_back(muon);

                              }
                           }
                           for (vector<tree::Electron>::iterator it = electrons->begin(); it != electrons->end(); it++){
                              auto g = *it;
                              selElectron ele;
                              ele.setAll(g);
                              ele.vec.SetPtEtaPhiM(ele.p.Pt(),ele.p.Eta(),ele.p.Phi(),0.51099895e-3);
                              if (testSelection(ele,LooseLeptons)){
                                 ele.deltaR1 = ele.vec.DeltaR(L.l1);
                                 ele.deltaR2 = ele.vec.DeltaR(L.l2);
                                 ele.chargeInt = (int) g.charge;
                                 if(matchLepton(ele,L)){
                                    ele.matched=true;
                                    L.matchedEleSize+=1;
                                    L.matchedLeptonSize+=1;
                                 }
                                 L.selElectrons.push_back(ele);
                              }
                           }

                     L.selLeptonSize=(L.selMuons.size()+L.selElectrons.size());
                     L.selMuonSize=(L.selMuons.size());
                     L.selElectronSize=(L.selElectrons.size());

                     //TLorentzVector temp4Vec(0.,0.,0.,0.);
                     //for (vector<selElectron>::iterator it = L.selElectrons.begin(); it != L.selElectrons.end(); it++){
                        //if(it->matched==false){
                           //temp4Vec=temp4Vec+it->vec;
                        //}
                     //}
                     //for (vector<selMuon>::iterator it = L.selMuons.begin(); it != L.selMuons.end(); it++){
                        //if(it->matched==false){
                           //temp4Vec=temp4Vec+it->vec;
                        //}
                     //}
                     //L.invAddLeptMass=temp4Vec.M();

                           
                           if(selection==ONZ) if(photons->size()!=0) decisionMapCutFlowFine[PHOTON1]=true;
                           if(selection==ONZ) if(photons->size()!=0) decisionMapCutFlowFine_weight[PHOTON1]=L.totalWeight;
                           if(!isData){
                                 L.totalWeight=L.totalWeight*GetScaleFactorAndErrorPhotons(L.selPhotons);
                           }
                           
                           if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1ID]=L.totalWeight;
                           if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1PT]=L.totalWeight;
                           if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1DR]=L.totalWeight;
                           
                           for (vector<tree::Jet>::iterator it = jets->begin(); it != jets->end(); it++){
                              auto j = *it;
                              selJet jet;
                              jet.setAll(j);
                              jet.vec.SetPtEtaPhiM(jet.p.Pt(),jet.p.Eta(),jet.p.Phi(),0.);
                              jet.deltaR1=jet.vec.DeltaR(L.l1);
                              jet.deltaR2=jet.vec.DeltaR(L.l2);
                              if (testSelection(jet,selection)){
                                 L.selJets.push_back(jet);
                                 L.calcHt += jet.vec.Pt();
                              }
                           }
                           
                           //for (auto& j : selJets) {
                              //L.calcHt += j->p.Pt();
                           //}
                           
                           selectedEvent=L;
                              
                           if(L.selPhotons.size()>0){
                              if(selection==ONZ) cutflow1Photon=true;
                              if(L.mll<101. && L.mll>81.){
                                 if(selection==ONZ) cutflowOnZ=true;
                                 if(selection==ONZ) decisionMapCutFlowFine[ZMASS]=true;
                                 if(selection==ONZ) decisionMapCutFlowFine_weight[ZMASS]=L.totalWeight;
                              }
                           }
                           return true;
                        }
                     }
                  }
               }
               return false;       
            //}
            return false;
         //}
         return false;
      }else{return false;}
   }else{return false;}
  
  
  return false;
}

bool HistogramProducer::SelectEventTriggerStudies(selectionType selection){
   if (CheckParticles()){//only Fill Histograms if nMu>=2 or nEle>=2, nGamma>=0
      if (CleaningTriggerStudies()){//decide in events with DiEle and DiMu where to contribute
         selEvent L;
         L.isDiElectron=isDiElectron;
         L.isDiMuon=isDiMuon;
         L.isMuonElectron=isMuonElectron;
         L.isElectronMuon=isElectronMuon;
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
         //L.MT2_val = MT2_val;
         L.evtHasGenPhotonVeto = evtHasGenPhotonVeto;
         

         
         if(!isData){
            if(L.isDiElectron){
               //L.totalWeight = totalWeight * GetScaleFactorAndErrorAlternative(pt1,eta1,isSignal,isDiElectron,*runNo)*GetScaleFactorAndErrorAlternative(pt2,eta2,isSignal,isDiElectron,*runNo);
               //L.totalWeight = totalWeight * GetScaleFactorAndErrorAlternative(L.pt1,L.eta1,isSignal,true,*runNo)*GetScaleFactorAndErrorAlternative(L.pt2,L.eta2,isSignal,true,*runNo);
               L.totalWeight = totalWeight * GetScaleFactorAndErrorFinal(L.pt1,L.eta1,isSignal,true)*GetScaleFactorAndErrorFinal(L.pt2,L.eta2,isSignal,true);
            }else{
               if(L.isDiMuon){
                  //L.totalWeight = totalWeight * GetScaleFactorAndErrorAlternative(pt1,eta1,isSignal,false,*runNo)*GetScaleFactorAndErrorAlternative(L.pt2,L.eta2,isSignal,false,*runNo);
                  L.totalWeight = totalWeight * GetScaleFactorAndErrorFinal(pt1,eta1,isSignal,false)*GetScaleFactorAndErrorFinal(L.pt2,L.eta2,isSignal,false);
               }
               else{
                  if(L.isMuonElectron||L.isElectronMuon){
                     L.totalWeight = totalWeight * GetScaleFactorAndErrorFinal(muons->at(0).p.Pt(),muons->at(0).p.Eta(),isSignal,false)*GetScaleFactorAndErrorFinal(electrons->at(0).p.Pt(),electrons->at(0).p.Eta(),isSignal,true);
                     //L.totalWeight = totalWeight * GetScaleFactorAndErrorFinal(myMuons.at(0).p.Pt(),myMuons.at(0).p.Eta(),isSignal,false)*GetScaleFactorAndErrorFinal(electrons->at(0).p.Pt(),electrons->at(0).p.Eta(),isSignal,true);
                  }
               }
            }
         }else{
            L.totalWeight=totalWeight;
         }
         if (L.isDiElectron){
            auto e1 = electrons->at(0);
            auto e2 = electrons->at(1);
            
            if(!GenPhotonVeto(config_veto)){
            
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
                  
                  
                   for (vector<tree::Muon>::iterator it = muons->begin(); it != muons->end(); it++){
                   //for (vector<tree::Muon>::iterator it = myMuons.begin(); it != myMuons.end(); it++){
                           auto g = *it;
                           selMuon muon;
                           muon.setAll(g);
                           muon.vec.SetPtEtaPhiM(muon.p.Pt(),muon.p.Eta(),muon.p.Phi(),0.1056584);
                           if (testSelection(muon,LooseLeptons)){
                              muon.deltaR1 = muon.vec.DeltaR(L.l1);
                              muon.deltaR2 = muon.vec.DeltaR(L.l2);
                              if(matchLepton(muon,L)){
                                 muon.matched=true;
                                 L.matchedMuSize+=1;
                                 L.matchedLeptonSize+=1;
                              }
                                 L.selMuons.push_back(muon);
                           }
                        }
                        for (vector<tree::Electron>::iterator it = electrons->begin(); it != electrons->end(); it++){
                           auto g = *it;
                           selElectron ele;
                           ele.setAll(g);
                           ele.vec.SetPtEtaPhiM(ele.p.Pt(),ele.p.Eta(),ele.p.Phi(),0.51099895e-3);
                           if (testSelection(ele,LooseLeptons)){
                              ele.deltaR1 = ele.vec.DeltaR(L.l1);
                              ele.deltaR2 = ele.vec.DeltaR(L.l2);
                              if(matchLepton(ele,L)){
                                 ele.matched=true;
                                 L.matchedEleSize+=1;
                                 L.matchedLeptonSize+=1;
                              }
                                 L.selElectrons.push_back(ele);
                           }
                        }
                  
                  L.selLeptonSize=(L.selMuons.size()+L.selElectrons.size());
                  
                  
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
      }
      else{
         if (L.isDiMuon){
            

            
         auto m1 = muons->at(0); 
         auto m2 = muons->at(1);
         //auto m1 = myMuons.at(0); 
         //auto m2 = myMuons.at(1);
            
         if(!GenPhotonVeto(config_veto)){
            
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
               
                                       for (vector<tree::Muon>::iterator it = muons->begin(); it != muons->end(); it++){
                                       //for (vector<tree::Muon>::iterator it = myMuons.begin(); it != myMuons.end(); it++){
                           auto g = *it;
                           selMuon muon;
                           muon.setAll(g);
                           muon.vec.SetPtEtaPhiM(muon.p.Pt(),muon.p.Eta(),muon.p.Phi(),0.1056584);
                           if (testSelection(muon,LooseLeptons)){
                              muon.deltaR1 = muon.vec.DeltaR(L.l1);
                              muon.deltaR2 = muon.vec.DeltaR(L.l2);
                              if(matchLepton(muon,L)){
                                 muon.matched=true;
                                 L.matchedMuSize+=1;
                                 L.matchedLeptonSize+=1;
                              }
                                 L.selMuons.push_back(muon);
                           }
                        }
                        for (vector<tree::Electron>::iterator it = electrons->begin(); it != electrons->end(); it++){
                           auto g = *it;
                           selElectron ele;
                           ele.setAll(g);
                           ele.vec.SetPtEtaPhiM(ele.p.Pt(),ele.p.Eta(),ele.p.Phi(),0.51099895e-3);
                           if (testSelection(ele,LooseLeptons)){
                              ele.deltaR1 = ele.vec.DeltaR(L.l1);
                              ele.deltaR2 = ele.vec.DeltaR(L.l2);
                              if(matchLepton(ele,L)){
                                 ele.matched=true;
                                 L.matchedEleSize+=1;
                                 L.matchedLeptonSize+=1;
                              }
                                 L.selElectrons.push_back(ele);
                           }
                        }
               
               L.selLeptonSize=(L.selMuons.size()+L.selElectrons.size());
               
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
      }else{
            if (L.isElectronMuon || L.isMuonElectron){
            
                  
               auto e0 = electrons->at(0); 
               auto m0 = muons->at(0);
               //auto m0 = myMuons.at(0);
               
               if(!GenPhotonVeto(config_veto)){
               
                  if(testSelection(e0,selection,isElectronMuon) && testSelection(m0,selection,isMuonElectron) && (chargeProduct < 0.) && ((selection==UNCUT)? true : (mll>50.))){

                     
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
                     
                     
                                             for (vector<tree::Muon>::iterator it = muons->begin(); it != muons->end(); it++){
                                             //for (vector<tree::Muon>::iterator it = myMuons.begin(); it != myMuons.end(); it++){
                           auto g = *it;
                           selMuon muon;
                           muon.setAll(g);
                           muon.vec.SetPtEtaPhiM(muon.p.Pt(),muon.p.Eta(),muon.p.Phi(),0.1056584);
                           if (testSelection(muon,LooseLeptons)){
                              muon.deltaR1 = muon.vec.DeltaR(L.l1);
                              muon.deltaR2 = muon.vec.DeltaR(L.l2);
                              if(matchLepton(muon,L)){
                                 muon.matched=true;
                                 L.matchedMuSize+=1;
                                 L.matchedLeptonSize+=1;
                              }
                                 L.selMuons.push_back(muon);
                           }
                        }
                        for (vector<tree::Electron>::iterator it = electrons->begin(); it != electrons->end(); it++){
                           auto g = *it;
                           selElectron ele;
                           ele.setAll(g);
                           ele.vec.SetPtEtaPhiM(ele.p.Pt(),ele.p.Eta(),ele.p.Phi(),0.51099895e-3);
                           if (testSelection(ele,LooseLeptons)){
                              ele.deltaR1 = ele.vec.DeltaR(L.l1);
                              ele.deltaR2 = ele.vec.DeltaR(L.l2);
                              if(matchLepton(ele,L)){
                                 ele.matched=true;
                                 L.matchedEleSize+=1;
                                 L.matchedLeptonSize+=1;
                              }
                                 L.selElectrons.push_back(ele);
                           }
                        }
                     
                     L.selLeptonSize=(L.selMuons.size()+L.selElectrons.size());
                     
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
      h1Maps["dilepEM"]=InitHistograms(DILEP);
      h1Maps["dilep"]=InitHistograms(DILEP);
   }
   if (config_selectionsToProcessMap[PHOTON]==true){
      h1Maps["1photonEE"]=InitHistograms(PHOTON); 
      h1Maps["1photonMM"]=InitHistograms(PHOTON); 
      h1Maps["1photonEM"]=InitHistograms(PHOTON); 
      h1Maps["1photon"]=InitHistograms(PHOTON); 
   }
   if (config_selectionsToProcessMap[SEL]==true){
      h1Maps["selEE"]=InitHistograms(SEL);
      h1Maps["selMM"]=InitHistograms(SEL); 
      h1Maps["selEM"]=InitHistograms(SEL); 
      h1Maps["sel"]=InitHistograms(SEL); 
   }
   if (config_selectionsToProcessMap[ONZ]==true){
      h1Maps["onZEE"]=InitHistograms(ONZ);
      h1Maps["onZMM"]=InitHistograms(ONZ); 
      h1Maps["onZEM"]=InitHistograms(ONZ); 
      h1Maps["onZ"]=InitHistograms(ONZ);
   }
   if (config_selectionsToProcessMap[ONZMET]==true){

      h1Maps["onZMet0100EE"]=InitHistograms(ONZ); //<100
      h1Maps["onZMet0100MM"]=InitHistograms(ONZ); 
      h1Maps["onZMet0100EM"]=InitHistograms(ONZ); 
      h1Maps["onZMet0100"]=InitHistograms(ONZ); 
      
      h1Maps["onZMet100EE"]=InitHistograms(ONZ); //>100
      h1Maps["onZMet100MM"]=InitHistograms(ONZ); 
      h1Maps["onZMet100EM"]=InitHistograms(ONZ); 
      h1Maps["onZMet100"]=InitHistograms(ONZ); 
      
      h1Maps["onZMet150EE"]=InitHistograms(ONZ); //>150
      h1Maps["onZMet150MM"]=InitHistograms(ONZ); 
      h1Maps["onZMet150EM"]=InitHistograms(ONZ); 
      h1Maps["onZMet150"]=InitHistograms(ONZ); 
      
      h1Maps["onZMet100150EE"]=InitHistograms(ONZ); //>100<150
      h1Maps["onZMet100150MM"]=InitHistograms(ONZ); 
      h1Maps["onZMet100150EM"]=InitHistograms(ONZ); 
      h1Maps["onZMet100150"]=InitHistograms(ONZ); 
   }
   if (config_selectionsToProcessMap[ONZG]==true){  
      h1Maps["onZGEE"]=InitHistograms(ONZ);
      h1Maps["onZGMM"]=InitHistograms(ONZ); 
      h1Maps["onZGEM"]=InitHistograms(ONZ); 
      h1Maps["onZG"]=InitHistograms(ONZ);
   }
   if (config_selectionsToProcessMap[ABOVEZG]==true){ 
      h1Maps["mllG110EE"]=InitHistograms(ONZ);
      h1Maps["mllG110MM"]=InitHistograms(ONZ); 
      h1Maps["mllG110EM"]=InitHistograms(ONZ); 
      h1Maps["mllG110"]=InitHistograms(ONZ); 
   }
   if (config_selectionsToProcessMap[EXO]==true){
      h1Maps["exoEE"]=InitHistograms(EXO);
      h1Maps["exoMM"]=InitHistograms(EXO); 
      h1Maps["exoEM"]=InitHistograms(EXO); 
      h1Maps["exo"]=InitHistograms(EXO); 
   }
   if (config_selectionsToProcessMap[EGRegression]==true){
      h1Maps["EGRegressionEE"]=InitHistograms(EGRegression);
      h1Maps["EGRegressionMM"]=InitHistograms(EGRegression); 
      h1Maps["EGRegressionEM"]=InitHistograms(EGRegression); 
      h1Maps["EGRegression"]=InitHistograms(EGRegression); 
   }
   if (config_selectionsToProcessMap[ControlRegionDY]==true){
      h1Maps["CRDYEE"]=InitHistograms(ONZ);
      h1Maps["CRDYMM"]=InitHistograms(ONZ); 
      h1Maps["CRDYEM"]=InitHistograms(ONZ); 
      h1Maps["CRDY"]=InitHistograms(ONZ); 
   }
   if (config_selectionsToProcessMap[ControlRegionTT]==true){
      h1Maps["CRTTEE"]=InitHistograms(ONZ);
      h1Maps["CRTTMM"]=InitHistograms(ONZ); 
      h1Maps["CRTTEM"]=InitHistograms(ONZ); 
      h1Maps["CRTT"]=InitHistograms(ONZ); 
      h1Maps["CRTT080EE"]=InitHistograms(ONZ);
      h1Maps["CRTT080MM"]=InitHistograms(ONZ); 
      h1Maps["CRTT080EM"]=InitHistograms(ONZ); 
      h1Maps["CRTT080"]=InitHistograms(ONZ); 
      h1Maps["CRTT80EE"]=InitHistograms(ONZ);
      h1Maps["CRTT80MM"]=InitHistograms(ONZ); 
      h1Maps["CRTT80EM"]=InitHistograms(ONZ); 
      h1Maps["CRTT80"]=InitHistograms(ONZ); 
   }
   if (config_selectionsToProcessMap[ControlRegionZZ]==true){
      h1Maps["CRZZEE"]=InitHistograms(ControlRegionZZ);
      h1Maps["CRZZMM"]=InitHistograms(ControlRegionZZ); 
      h1Maps["CRZZEM"]=InitHistograms(ControlRegionZZ); 
      h1Maps["CRZZ"]=InitHistograms(ControlRegionZZ); 
   }
   if (config_selectionsToProcessMap[ControlRegionWZ]==true){
      h1Maps["CRWZEE"]=InitHistograms(ControlRegionZZ);
      h1Maps["CRWZMM"]=InitHistograms(ControlRegionZZ); 
      h1Maps["CRWZEM"]=InitHistograms(ControlRegionZZ); 
      h1Maps["CRWZ"]=InitHistograms(ControlRegionZZ); 
   }
   if (config_selectionsToProcessMap[ControlRegionWW]==true){
      h1Maps["CRWWEE"]=InitHistograms(DILEP);
      h1Maps["CRWWLL"]=InitHistograms(DILEP);
      h1Maps["CRWWMM"]=InitHistograms(DILEP); 
      h1Maps["CRWWEM"]=InitHistograms(DILEP); 
      h1Maps["CRWW"]=InitHistograms(DILEP); 
   }
   if (config_selectionsToProcessMap[ValidationRegion]==true){
      h1Maps["VREE"]=InitHistograms(ONZ);
      h1Maps["VRMM"]=InitHistograms(ONZ); 
      h1Maps["VREM"]=InitHistograms(ONZ); 
      h1Maps["VR"]=InitHistograms(ONZ); 
      h1Maps["VR080EE"]=InitHistograms(ONZ);
      h1Maps["VR080MM"]=InitHistograms(ONZ); 
      h1Maps["VR080EM"]=InitHistograms(ONZ); 
      h1Maps["VR080"]=InitHistograms(ONZ); 
      h1Maps["VR80EE"]=InitHistograms(ONZ);
      h1Maps["VR80MM"]=InitHistograms(ONZ); 
      h1Maps["VR80EM"]=InitHistograms(ONZ); 
      h1Maps["VR80"]=InitHistograms(ONZ); 
   }
   if (config_selectionsToProcessMap[SEL]==true){
      h2Maps["selEE"]=Init2DHistograms(SEL);
      h2Maps["selMM"]=Init2DHistograms(SEL); 
      h2Maps["selEM"]=Init2DHistograms(SEL); 
      h2Maps["sel"]=Init2DHistograms(SEL); 
   }
   if (config_selectionsToProcessMap[ONZ]==true){
      h2Maps["onZEE"]=Init2DHistograms(ONZ);
      h2Maps["onZMM"]=Init2DHistograms(ONZ); 
      h2Maps["onZEM"]=Init2DHistograms(ONZ); 
      h2Maps["onZ"]=Init2DHistograms(ONZ); 
   }
   if (config_selectionsToProcessMap[DILEP]==true){
      h2Maps["dilepEE"]=Init2DHistograms(DILEP);
      h2Maps["dilepMM"]=Init2DHistograms(DILEP); 
      h2Maps["dilepEM"]=Init2DHistograms(DILEP); 
      h2Maps["dilep"]=Init2DHistograms(DILEP); 
   }
   if (config_selectionsToProcessMap[EXO]==true){
      h2Maps["exoEE"]=Init2DHistograms(EXO);
      h2Maps["exoMM"]=Init2DHistograms(EXO); 
      h2Maps["exoEM"]=Init2DHistograms(EXO); 
      h2Maps["exo"]=Init2DHistograms(EXO); 
   }
}

void HistogramProducer::InitSignalScanHistos(string masspoint){
   s1Maps[masspoint]["sig"]=InitSignalScanHistograms(ONZ);
   s1Maps[masspoint]["sigEE"]=InitSignalScanHistograms(ONZ);
   s1Maps[masspoint]["sigMM"]=InitSignalScanHistograms(ONZ);
   s1Maps[masspoint]["sigEM"]=InitSignalScanHistograms(ONZ);
   s1Maps[masspoint]["sig080"]=InitSignalScanHistograms(ONZ);
   s1Maps[masspoint]["sig080EE"]=InitSignalScanHistograms(ONZ);
   s1Maps[masspoint]["sig080MM"]=InitSignalScanHistograms(ONZ);
   s1Maps[masspoint]["sig080EM"]=InitSignalScanHistograms(ONZ);
   s1Maps[masspoint]["sig80"]=InitSignalScanHistograms(ONZ);
   s1Maps[masspoint]["sig80EE"]=InitSignalScanHistograms(ONZ);
   s1Maps[masspoint]["sig80MM"]=InitSignalScanHistograms(ONZ);
   s1Maps[masspoint]["sig80EM"]=InitSignalScanHistograms(ONZ);
   if(config_dosignalscanSplit){
      s1Maps[masspoint]["sig_gg"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_ggEE"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_ggMM"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_ggEM"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_zz"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_zzEE"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_zzMM"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_zzEM"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_gz"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_gzEE"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_gzMM"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_gzEM"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_gg080"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_gg080EE"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_gg080MM"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_gg080EM"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_zz080"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_zz080EE"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_zz080MM"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_zz080EM"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_gz080"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_gz080EE"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_gz080MM"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_gz080EM"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_gg80"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_gg80EE"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_gg80MM"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_gg80EM"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_zz80"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_zz80EE"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_zz80MM"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_zz80EM"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_gz80"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_gz80EE"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_gz80MM"]=InitSignalScanHistograms(ONZ);
      s1Maps[masspoint]["sig_gz80EM"]=InitSignalScanHistograms(ONZ);
   }
   //s1Maps[masspoint]["HIGHsignal"]=InitSignalScanHistograms(ONZ);
   //s1Maps[masspoint]["HIGHsignalEE"]=InitSignalScanHistograms(ONZ);
   //s1Maps[masspoint]["HIGHsignalMM"]=InitSignalScanHistograms(ONZ);
   //s1Maps[masspoint]["HIGHsignalEM"]=InitSignalScanHistograms(ONZ);
}


void HistogramProducer::InitCutFlowHistos(){
   c1Maps["cutFlow_onZEE"]=InitCutFlowHistograms(ONZ);
   c1Maps["cutFlow_onZMM"]=InitCutFlowHistograms(ONZ);
   c1Maps["cutFlow_onZEM"]=InitCutFlowHistograms(ONZ);
}


void HistogramProducer::InitCutFlowHistos_Fine(){
   c1Maps["cutFlow_Fine_onZEE"]=InitCutFlowHistograms_Fine(ONZ);
   c1Maps["cutFlow_Fine_onZMM"]=InitCutFlowHistograms_Fine(ONZ);
   c1Maps["cutFlow_Fine_onZEM"]=InitCutFlowHistograms_Fine(ONZ);
}

   
void HistogramProducer::InitTriggerStudiesHistos(){
   eff1Maps["trigDilepEE"]=InitTriggerStudies(TRIGDILEP);
   eff1Maps["trigDilepMM"]=InitTriggerStudies(TRIGDILEP);
   eff1Maps["trigDilepEM"]=InitTriggerStudies(TRIGDILEP);
   eff1Maps["trigSelEE"]=InitTriggerStudies(TRIGSEL);
   eff1Maps["trigSelMM"]=InitTriggerStudies(TRIGSEL);
   eff1Maps["trigSelEM"]=InitTriggerStudies(TRIGSEL);
   eff1Maps["trigOnZEE"]=InitTriggerStudies(TRIGONZ);
   eff1Maps["trigOnZMM"]=InitTriggerStudies(TRIGONZ);
   eff1Maps["trigOnZEM"]=InitTriggerStudies(TRIGONZ);
   eff1Maps["trigDilepEE_ptcuts"]=InitTriggerStudies(TRIGDILEP_ptcuts);
   eff1Maps["trigDilepMM_ptcuts"]=InitTriggerStudies(TRIGDILEP_ptcuts);
   eff1Maps["trigDilepEM_ptcuts"]=InitTriggerStudies(TRIGDILEP_ptcuts);
   //eff1Maps["trigDilepEE_pt1cut"]=InitTriggerStudies(TRIGDILEP_ptcuts);
   //eff1Maps["trigDilepMM_pt1cut"]=InitTriggerStudies(TRIGDILEP_ptcuts);
   //eff1Maps["trigDilepEM_pt1cut"]=InitTriggerStudies(TRIGDILEP_ptcuts);
   //eff1Maps["trigDilepEE_pt2cut"]=InitTriggerStudies(TRIGDILEP_ptcuts);
   //eff1Maps["trigDilepMM_pt2cut"]=InitTriggerStudies(TRIGDILEP_ptcuts);
   //eff1Maps["trigDilepEM_pt2cut"]=InitTriggerStudies(TRIGDILEP_ptcuts);
   eff1Maps["trigSelEE_ptcuts"]=InitTriggerStudies(TRIGSEL_ptcuts);
   eff1Maps["trigSelMM_ptcuts"]=InitTriggerStudies(TRIGSEL_ptcuts);
   eff1Maps["trigSelEM_ptcuts"]=InitTriggerStudies(TRIGSEL_ptcuts);
   eff1Maps["trigOnZEE_ptcuts"]=InitTriggerStudies(TRIGONZ_ptcuts);
   eff1Maps["trigOnZMM_ptcuts"]=InitTriggerStudies(TRIGONZ_ptcuts);
   eff1Maps["trigOnZEM_ptcuts"]=InitTriggerStudies(TRIGONZ_ptcuts);

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
   
   float weight_TopPt=1.;
   float weight_LeptonPt=1.;
   float weight_nIsr=1.;
   float weight_EWKinoPt=1.;
   
   float tempWeight=ev.totalWeight;
   
   if(doWeights_TopPt){
      weight_TopPt = topPtReweighting(*topPt1,*topPt2);
   }
   if(doWeights_nISR){
      weight_nIsr = isrReweighting(*nISR,false);
   }
   if(doWeights_EWKinoPairPt){
      weight_EWKinoPt = isrReweightingEWK(*EWKinoPairPt,false);
   }
   if(doWeights_LeptonPairPt){
      weight_LeptonPt = isrReweightingEWK(*leptonPairPt,false);
   }
   
   
   
   tempWeight = tempWeight*weight_TopPt*weight_LeptonPt*weight_nIsr*weight_EWKinoPt;
   
   m.at(WEIGHT_EWKINOPAIRPT).Fill(weight_EWKinoPt);
   m.at(WEIGHT_LEPTONPAIRPT).Fill(weight_LeptonPt);
   m.at(WEIGHT_NISR).Fill(weight_nIsr);
   m.at(WEIGHT_TOPPT).Fill(weight_TopPt);
   
   
   m.at(ETMISS).Fill(ev.ETmiss, tempWeight);
   m.at(PT1).Fill(ev.pt1, tempWeight);
   m.at(PT2).Fill(ev.pt2, tempWeight);
   m.at(MLL).Fill(ev.mll,tempWeight);   
   m.at(NPHOTONS).Fill(ev.selPhotons.size(),tempWeight);
   m.at(NVTX).Fill(*nGoodVertices,tempWeight);
   //m.at(HT).Fill(*ht,tempWeight);
   m.at(HT).Fill(ev.calcHt,tempWeight);
   m.at(GENHT).Fill(*genHt,tempWeight);
   m.at(NJETS).Fill(ev.selJets.size(),tempWeight);
   m.at(ETA1).Fill(fabs(ev.eta1),tempWeight);
   m.at(ETA2).Fill(fabs(ev.eta2),tempWeight);
   m.at(PHI1).Fill(fabs(ev.phi2),tempWeight);
   m.at(PHI2).Fill(fabs(ev.phi2),tempWeight);
   m.at(DeltaEtaLL).Fill(fabs(ev.eta1-ev.eta2),tempWeight);
   m.at(DeltaPhiLL).Fill(fabs(ev.phi1-ev.phi2),tempWeight);
   m.at(DeltaRLL).Fill(fabs(ev.deltaRll),tempWeight);
   m.at(DeltaEtaLL_neg).Fill(ev.eta1-ev.eta2,tempWeight);
   m.at(DeltaPhiLL_neg).Fill(ev.phi1-ev.phi2,tempWeight);
   m.at(DeltaRLL_neg).Fill(ev.deltaRll,tempWeight);
   m.at(ZPT).Fill((ev.l1+ev.l2).Pt(),tempWeight);
   m.at(MTLL).Fill((ev.l1+ev.l2).Mt(),tempWeight);
   m.at(ST).Fill(ev.pt1+ev.pt2,tempWeight);
   //m.at(MT2).Fill(ev.MT2_val,tempWeight);
   m.at(DeltaPhiLLMet).Fill(fabs((ev.l1+ev.l2).Phi() - ev.ETmiss_vec.Phi()),tempWeight);
   m.at(DeltaEtaLLMet).Fill(fabs((ev.l1+ev.l2).Eta() - ev.ETmiss_vec.Eta()),tempWeight);
   m.at(DeltaRLLMet).Fill(fabs((ev.l1+ev.l2).DeltaR(ev.ETmiss_vec)),tempWeight);
   
   int jetSize=ev.selJets.size();

   //cout<<"here"<<endl;

   if(jetSize>0){
      m.at(JetPt1).Fill(ev.selJets.at(0).p.Pt(),tempWeight);
      m.at(JetPhi1).Fill(ev.selJets.at(0).p.Phi(),tempWeight);
      m.at(JetEta1).Fill(ev.selJets.at(0).p.Eta(),tempWeight);
      if(jetSize>1){
         m.at(JetPt2).Fill(ev.selJets.at(1).p.Pt(),tempWeight);
         m.at(JetPhi2).Fill(ev.selJets.at(1).p.Phi(),tempWeight);
         m.at(JetEta2).Fill(ev.selJets.at(1).p.Eta(),tempWeight);
         if(jetSize>2){
            m.at(JetPt3).Fill(ev.selJets.at(2).p.Pt(),tempWeight);
            m.at(JetPhi3).Fill(ev.selJets.at(2).p.Phi(),tempWeight);
            m.at(JetEta3).Fill(ev.selJets.at(2).p.Eta(),tempWeight);
            if(jetSize>3){
               m.at(JetPt4).Fill(ev.selJets.at(3).p.Pt(),tempWeight);
               m.at(JetPhi4).Fill(ev.selJets.at(3).p.Phi(),tempWeight);
               m.at(JetEta4).Fill(ev.selJets.at(3).p.Eta(),tempWeight);
            }
         }
      }
   }
   //cout<<"here2"<<endl;
   m.at(NElectrons).Fill(ev.selElectrons.size(),tempWeight);
   m.at(NMuons).Fill(ev.selMuons.size(),tempWeight);
   
   
   if(ev.evtHasGenPhotonVeto) m.at(VetoCompare).Fill(1.,tempWeight);
   if(!ev.evtHasGenPhotonVeto) m.at(VetoCompare).Fill(0.,tempWeight);
   
   if(withPhoton){
      m.at(PTG1).Fill(ev.selPhotons.at(0).p.Pt(),tempWeight);
      m.at(PHIG1).Fill(ev.selPhotons.at(0).p.Phi(),tempWeight);
      m.at(ETAG1).Fill(ev.selPhotons.at(0).p.Eta(),tempWeight);
      m.at(SIGMAIETAIETAG1).Fill(ev.selPhotons.at(0).sigmaIetaIeta,tempWeight);
      m.at(SIGMAIPHIIPHIG1).Fill(ev.selPhotons.at(0).sigmaIphiIphi,tempWeight);
      m.at(R9).Fill(ev.selPhotons.at(0).r9,tempWeight);
      m.at(HOVERE).Fill(ev.selPhotons.at(0).hOverE,tempWeight);
      m.at(DELTARGL1).Fill(ev.selPhotons.at(0).deltaR1,tempWeight);
      m.at(DELTARGL2).Fill(ev.selPhotons.at(0).deltaR2,tempWeight);
      m.at(DeltaEtaLLG).Fill(fabs((ev.l1+ev.l2).Eta()-ev.selPhotons.at(0).vec.Eta()), tempWeight);
      m.at(DeltaPhiLLG).Fill(fabs((ev.l1+ev.l2).Phi()-ev.selPhotons.at(0).vec.Phi()),tempWeight);
      m.at(DeltaRLLG).Fill(fabs((ev.l1+ev.l2).DeltaR(ev.selPhotons.at(0).vec)),tempWeight);
      m.at(MTLLG).Fill((ev.l1+ev.l2+ev.selPhotons.at(0).vec).Mt(),tempWeight);
      m.at(MTL1MET).Fill((ev.l1+ev.ETmiss_vec).Mt(),tempWeight);
      m.at(MTL2MET).Fill((ev.l2+ev.ETmiss_vec).Mt(),tempWeight);
      m.at(MTGMET).Fill((ev.selPhotons.at(0).vec+ev.ETmiss_vec).Mt(),tempWeight);
      m.at(MTLLMET).Fill((ev.l1+ev.l2+ev.ETmiss_vec).Mt(),tempWeight);
      m.at(MTLLGMET).Fill((ev.l1+ev.l2+ev.selPhotons.at(0).vec+ev.ETmiss_vec).Mt(),tempWeight);
      m.at(STG).Fill(ev.pt1+ev.pt2+ev.selPhotons.at(0).p.Pt(),tempWeight);
      m.at(STMET).Fill(ev.pt1+ev.pt2+ev.selPhotons.at(0).p.Pt()+ev.ETmiss,tempWeight);
      m.at(MLLG).Fill((ev.l1+ev.l2+ev.selPhotons.at(0).vec).M(),tempWeight);
      //m.at(PT_llg).Fill((ev.l1+ev.l2+ev.selPhotons.at(0).vec).Pt(),tempWeight);
      //m.at(MZG_exo).Fill(40./150. * (ev.l1+ev.l2+ev.selPhotons.at(0).vec).M(),tempWeight);
      
      if(FindGenPhotonMatch(ev.selPhotons.at(0))) m.at(genPhotonPT).Fill(GetGenPhotonMatch(ev.selPhotons.at(0)).p.Pt(),tempWeight);
      
      if(ev.evtHasGenPhotonVeto) if(FindGenPhotonMatch(ev.selPhotons.at(0))) m.at(genPhotonPT_Veto).Fill(GetGenPhotonMatch(ev.selPhotons.at(0)).p.Pt(),tempWeight);
      if(ev.evtHasGenPhotonVeto) m.at(PTG1_Veto).Fill(ev.selPhotons.at(0).p.Pt(),tempWeight);
      if(!ev.evtHasGenPhotonVeto) if(FindGenPhotonMatch(ev.selPhotons.at(0))) m.at(genPhotonPT_NoVeto).Fill(GetGenPhotonMatch(ev.selPhotons.at(0)).p.Pt(),tempWeight);
      if(!ev.evtHasGenPhotonVeto) m.at(PTG1_NoVeto).Fill(ev.selPhotons.at(0).p.Pt(),tempWeight);
      //if(ev.evtHasGenPhotonVeto) m.at(VetoCompare).Fill(1.,tempWeight);
      //if(!ev.evtHasGenPhotonVeto) m.at(VetoCompare).Fill(0.,tempWeight);
      
      //for (vector<tree::GenParticle>::iterator it = genParticles->begin(); it != genParticles->end(); it++){
         //int Id_mother=-1;
         //if (abs(it->pdgId)==22){
            //Id_mother = abs(it->motherId);
            //m.at(gammaMotherID).Fill(Id_mother,tempWeight);
         //}
      //}
      
   }
}
void HistogramProducer::FillerZZ(selEvent& ev, map<Histograms1D,TH1F>& m,bool withPhoton,selLepton& l1, selLepton& l2, selLepton& l3, selLepton& l4){
   
   float weight_TopPt=1.;
   float weight_LeptonPt=1.;
   float weight_nIsr=1.;
   float weight_EWKinoPt=1.;
   
   float tempWeight=ev.totalWeight;
   
   if(doWeights_TopPt){
      weight_TopPt = topPtReweighting(*topPt1,*topPt2);
   }
   if(doWeights_nISR){
      weight_nIsr = isrReweighting(*nISR,false);
   }
   if(doWeights_EWKinoPairPt){
      weight_EWKinoPt = isrReweightingEWK(*EWKinoPairPt,false);
   }
   if(doWeights_LeptonPairPt){
      weight_LeptonPt = isrReweightingEWK(*leptonPairPt,false);
   }
   
   tempWeight = tempWeight*weight_TopPt*weight_LeptonPt*weight_nIsr*weight_EWKinoPt;
   
   m.at(WEIGHT_EWKINOPAIRPT).Fill(weight_EWKinoPt);
   m.at(WEIGHT_LEPTONPAIRPT).Fill(weight_LeptonPt);
   m.at(WEIGHT_NISR).Fill(weight_nIsr);
   m.at(WEIGHT_TOPPT).Fill(weight_TopPt);
   
   
   m.at(ETMISS).Fill(ev.ETmiss, tempWeight);
   m.at(PT1).Fill(l1.p.Pt(), tempWeight);
   m.at(PT2).Fill(l2.p.Pt(), tempWeight);
   m.at(PT3).Fill(l3.p.Pt(), tempWeight);
   m.at(PT4).Fill(l4.p.Pt(), tempWeight);
   m.at(MLL).Fill((l1.vec+l2.vec).M(),tempWeight);   
   m.at(MLL2).Fill((l3.vec+l4.vec).M(),tempWeight);   
   m.at(MLLLL).Fill((l1.vec+l2.vec+l3.vec+l4.vec).M(),tempWeight);   
   m.at(NPHOTONS).Fill(ev.selPhotons.size(),tempWeight);
   m.at(NVTX).Fill(*nGoodVertices,tempWeight);
   //m.at(HT).Fill(*ht,tempWeight);
   m.at(HT).Fill(ev.calcHt,tempWeight);
   m.at(GENHT).Fill(*genHt,tempWeight);
   m.at(NJETS).Fill(ev.selJets.size(),tempWeight);
   m.at(ETA1).Fill(fabs(l1.p.Eta()),tempWeight);
   m.at(ETA2).Fill(fabs(l2.p.Eta()),tempWeight);
   m.at(ETA3).Fill(fabs(l3.p.Eta()),tempWeight);
   m.at(ETA4).Fill(fabs(l4.p.Eta()),tempWeight);
   m.at(PHI1).Fill(fabs(l1.p.Phi()),tempWeight);
   m.at(PHI2).Fill(fabs(l2.p.Phi()),tempWeight);
   m.at(PHI3).Fill(fabs(l3.p.Phi()),tempWeight);
   m.at(PHI4).Fill(fabs(l4.p.Phi()),tempWeight);
   m.at(ZPT).Fill((l1.vec+l2.vec).Pt(),tempWeight);
   m.at(ZPT2).Fill((l3.vec+l4.vec).Pt(),tempWeight);
   m.at(ST).Fill(ev.pt1+ev.pt2,tempWeight);
   m.at(MT2).Fill(ev.MT2_val,tempWeight);



   int jetSize=ev.selJets.size();

   if(jetSize>0){
      m.at(JetPt1).Fill(ev.selJets.at(0).p.Pt(),tempWeight);
      m.at(JetPhi1).Fill(ev.selJets.at(0).p.Phi(),tempWeight);
      m.at(JetEta1).Fill(ev.selJets.at(0).p.Eta(),tempWeight);
      if(jetSize>1){
         m.at(JetPt2).Fill(ev.selJets.at(1).p.Pt(),tempWeight);
         m.at(JetPhi2).Fill(ev.selJets.at(1).p.Phi(),tempWeight);
         m.at(JetEta2).Fill(ev.selJets.at(1).p.Eta(),tempWeight);
         if(jetSize>2){
            m.at(JetPt3).Fill(ev.selJets.at(2).p.Pt(),tempWeight);
            m.at(JetPhi3).Fill(ev.selJets.at(2).p.Phi(),tempWeight);
            m.at(JetEta3).Fill(ev.selJets.at(2).p.Eta(),tempWeight);
            if(jetSize>3){
               m.at(JetPt4).Fill(ev.selJets.at(3).p.Pt(),tempWeight);
               m.at(JetPhi4).Fill(ev.selJets.at(3).p.Phi(),tempWeight);
               m.at(JetEta4).Fill(ev.selJets.at(3).p.Eta(),tempWeight);
            }
         }
      }
   }




   m.at(NElectrons).Fill(ev.selElectrons.size(),tempWeight);
   m.at(NMuons).Fill(ev.selMuons.size(),tempWeight);

   
   if(withPhoton){
      m.at(PTG1).Fill(ev.selPhotons.at(0).p.Pt(),tempWeight);
      m.at(PHIG1).Fill(ev.selPhotons.at(0).p.Phi(),tempWeight);
      m.at(ETAG1).Fill(ev.selPhotons.at(0).p.Eta(),tempWeight);
      m.at(SIGMAIETAIETAG1).Fill(ev.selPhotons.at(0).sigmaIetaIeta,tempWeight);
      m.at(SIGMAIPHIIPHIG1).Fill(ev.selPhotons.at(0).sigmaIphiIphi,tempWeight);
      m.at(R9).Fill(ev.selPhotons.at(0).r9,tempWeight);
      m.at(HOVERE).Fill(ev.selPhotons.at(0).hOverE,tempWeight);
      m.at(DELTARGL1).Fill(ev.selPhotons.at(0).deltaR1,tempWeight);
      m.at(DELTARGL2).Fill(ev.selPhotons.at(0).deltaR2,tempWeight);
      m.at(DeltaEtaLLG).Fill(fabs((ev.l1+ev.l2).Eta()-ev.selPhotons.at(0).vec.Eta()), tempWeight);
      m.at(DeltaPhiLLG).Fill(fabs((ev.l1+ev.l2).Phi()-ev.selPhotons.at(0).vec.Phi()),tempWeight);
      m.at(DeltaRLLG).Fill(fabs((ev.l1+ev.l2).DeltaR(ev.selPhotons.at(0).vec)),tempWeight);
      m.at(MTLLG).Fill((ev.l1+ev.l2+ev.selPhotons.at(0).vec).Mt(),tempWeight);
      m.at(MTL1MET).Fill((ev.l1+ev.ETmiss_vec).Mt(),tempWeight);
      m.at(MTL2MET).Fill((ev.l2+ev.ETmiss_vec).Mt(),tempWeight);
      m.at(MTGMET).Fill((ev.selPhotons.at(0).vec+ev.ETmiss_vec).Mt(),tempWeight);
      m.at(MTLLMET).Fill((ev.l1+ev.l2+ev.ETmiss_vec).Mt(),tempWeight);
      m.at(MTLLGMET).Fill((ev.l1+ev.l2+ev.selPhotons.at(0).vec+ev.ETmiss_vec).Mt(),tempWeight);
      m.at(STG).Fill(ev.pt1+ev.pt2+ev.selPhotons.at(0).p.Pt(),tempWeight);
      m.at(STMET).Fill(ev.pt1+ev.pt2+ev.selPhotons.at(0).p.Pt()+ev.ETmiss,tempWeight);
      m.at(MLLG).Fill((ev.l1+ev.l2+ev.selPhotons.at(0).vec).M(),tempWeight);
      
      if(FindGenPhotonMatch(ev.selPhotons.at(0))) m.at(genPhotonPT).Fill(GetGenPhotonMatch(ev.selPhotons.at(0)).p.Pt(),tempWeight);
      
      if(ev.evtHasGenPhotonVeto) if(FindGenPhotonMatch(ev.selPhotons.at(0))) m.at(genPhotonPT_Veto).Fill(GetGenPhotonMatch(ev.selPhotons.at(0)).p.Pt(),tempWeight);
      if(ev.evtHasGenPhotonVeto) m.at(PTG1_Veto).Fill(ev.selPhotons.at(0).p.Pt(),tempWeight);
      if(!ev.evtHasGenPhotonVeto) if(FindGenPhotonMatch(ev.selPhotons.at(0))) m.at(genPhotonPT_NoVeto).Fill(GetGenPhotonMatch(ev.selPhotons.at(0)).p.Pt(),tempWeight);
      if(!ev.evtHasGenPhotonVeto) m.at(PTG1_NoVeto).Fill(ev.selPhotons.at(0).p.Pt(),tempWeight);
      
   }
}

void HistogramProducer::FillerWZ(selEvent& ev, map<Histograms1D,TH1F>& m,bool withPhoton,selLepton& l1, selLepton& l2, selLepton& l3){
   
   float weight_TopPt=1.;
   float weight_LeptonPt=1.;
   float weight_nIsr=1.;
   float weight_EWKinoPt=1.;
   
   float tempWeight=ev.totalWeight;
   
   if(doWeights_TopPt){
      weight_TopPt = topPtReweighting(*topPt1,*topPt2);
   }
   if(doWeights_nISR){
      weight_nIsr = isrReweighting(*nISR,false);
   }
   if(doWeights_EWKinoPairPt){
      weight_EWKinoPt = isrReweightingEWK(*EWKinoPairPt,false);
   }
   if(doWeights_LeptonPairPt){
      weight_LeptonPt = isrReweightingEWK(*leptonPairPt,false);
   }
   
   tempWeight = tempWeight*weight_TopPt*weight_LeptonPt*weight_nIsr*weight_EWKinoPt;
   
   m.at(WEIGHT_EWKINOPAIRPT).Fill(weight_EWKinoPt);
   m.at(WEIGHT_LEPTONPAIRPT).Fill(weight_LeptonPt);
   m.at(WEIGHT_NISR).Fill(weight_nIsr);
   m.at(WEIGHT_TOPPT).Fill(weight_TopPt);
   
   
   m.at(ETMISS).Fill(ev.ETmiss, tempWeight);
   m.at(PT1).Fill(l1.p.Pt(), tempWeight);
   m.at(PT2).Fill(l2.p.Pt(), tempWeight);
   m.at(PT3).Fill(l3.p.Pt(), tempWeight);
   m.at(MLL).Fill((l1.vec+l2.vec).M(),tempWeight);   
   m.at(NPHOTONS).Fill(ev.selPhotons.size(),tempWeight);
   m.at(NVTX).Fill(*nGoodVertices,tempWeight);
   //m.at(HT).Fill(*ht,tempWeight);
   m.at(HT).Fill(ev.calcHt,tempWeight);
   m.at(GENHT).Fill(*genHt,tempWeight);
   m.at(NJETS).Fill(ev.selJets.size(),tempWeight);
   m.at(ETA1).Fill(fabs(l1.p.Eta()),tempWeight);
   m.at(ETA2).Fill(fabs(l2.p.Eta()),tempWeight);
   m.at(ETA3).Fill(fabs(l3.p.Eta()),tempWeight);
   m.at(PHI1).Fill(fabs(l1.p.Phi()),tempWeight);
   m.at(PHI2).Fill(fabs(l2.p.Phi()),tempWeight);
   m.at(PHI3).Fill(fabs(l3.p.Phi()),tempWeight);
   m.at(ZPT).Fill((l1.vec+l2.vec).Pt(),tempWeight);
   m.at(ST).Fill(ev.pt1+ev.pt2,tempWeight);
   m.at(MT2).Fill(ev.MT2_val,tempWeight);

   m.at(DeltaEtaLL).Fill(fabs(l1.p.Eta()-l2.p.Eta()),tempWeight);
   m.at(DeltaPhiLL).Fill(fabs(l1.p.Phi()-l2.p.Phi()),tempWeight);
   m.at(DeltaRLL).Fill(fabs(l1.vec.DeltaR(l2.vec)),tempWeight);

   m.at(DeltaEtaLL_neg).Fill(l1.p.Eta()-l2.p.Eta(),tempWeight);
   m.at(DeltaPhiLL_neg).Fill(l1.p.Phi()-l2.p.Phi(),tempWeight);
   m.at(DeltaRLL_neg).Fill(l1.vec.DeltaR(l2.vec),tempWeight);

   int jetSize=ev.selJets.size();

   if(jetSize>0){
      m.at(JetPt1).Fill(ev.selJets.at(0).p.Pt(),tempWeight);
      m.at(JetPhi1).Fill(ev.selJets.at(0).p.Phi(),tempWeight);
      m.at(JetEta1).Fill(ev.selJets.at(0).p.Eta(),tempWeight);
      if(jetSize>1){
         m.at(JetPt2).Fill(ev.selJets.at(1).p.Pt(),tempWeight);
         m.at(JetPhi2).Fill(ev.selJets.at(1).p.Phi(),tempWeight);
         m.at(JetEta2).Fill(ev.selJets.at(1).p.Eta(),tempWeight);
         if(jetSize>2){
            m.at(JetPt3).Fill(ev.selJets.at(2).p.Pt(),tempWeight);
            m.at(JetPhi3).Fill(ev.selJets.at(2).p.Phi(),tempWeight);
            m.at(JetEta3).Fill(ev.selJets.at(2).p.Eta(),tempWeight);
            if(jetSize>3){
               m.at(JetPt4).Fill(ev.selJets.at(3).p.Pt(),tempWeight);
               m.at(JetPhi4).Fill(ev.selJets.at(3).p.Phi(),tempWeight);
               m.at(JetEta4).Fill(ev.selJets.at(3).p.Eta(),tempWeight);
            }
         }
      }
   }



   m.at(MTL3MET).Fill((l1.vec+ev.ETmiss_vec).Mt(),tempWeight);

   m.at(NElectrons).Fill(ev.selElectrons.size(),tempWeight);
   m.at(NMuons).Fill(ev.selMuons.size(),tempWeight);

   
   if(withPhoton){
      m.at(PTG1).Fill(ev.selPhotons.at(0).p.Pt(),tempWeight);
      m.at(PHIG1).Fill(ev.selPhotons.at(0).p.Phi(),tempWeight);
      m.at(ETAG1).Fill(ev.selPhotons.at(0).p.Eta(),tempWeight);
      m.at(SIGMAIETAIETAG1).Fill(ev.selPhotons.at(0).sigmaIetaIeta,tempWeight);
      m.at(SIGMAIPHIIPHIG1).Fill(ev.selPhotons.at(0).sigmaIphiIphi,tempWeight);
      m.at(R9).Fill(ev.selPhotons.at(0).r9,tempWeight);
      m.at(HOVERE).Fill(ev.selPhotons.at(0).hOverE,tempWeight);
      m.at(DELTARGL1).Fill(ev.selPhotons.at(0).deltaR1,tempWeight);
      m.at(DELTARGL2).Fill(ev.selPhotons.at(0).deltaR2,tempWeight);
      m.at(DeltaEtaLLG).Fill(fabs((ev.l1+ev.l2).Eta()-ev.selPhotons.at(0).vec.Eta()), tempWeight);
      m.at(DeltaPhiLLG).Fill(fabs((ev.l1+ev.l2).Phi()-ev.selPhotons.at(0).vec.Phi()),tempWeight);
      m.at(DeltaRLLG).Fill(fabs((ev.l1+ev.l2).DeltaR(ev.selPhotons.at(0).vec)),tempWeight);
      m.at(MTLLG).Fill((ev.l1+ev.l2+ev.selPhotons.at(0).vec).Mt(),tempWeight);
      m.at(MTL1MET).Fill((ev.l1+ev.ETmiss_vec).Mt(),tempWeight);
      m.at(MTL2MET).Fill((ev.l2+ev.ETmiss_vec).Mt(),tempWeight);
      m.at(MTGMET).Fill((ev.selPhotons.at(0).vec+ev.ETmiss_vec).Mt(),tempWeight);
      m.at(MTLLMET).Fill((ev.l1+ev.l2+ev.ETmiss_vec).Mt(),tempWeight);
      m.at(MTLLGMET).Fill((ev.l1+ev.l2+ev.selPhotons.at(0).vec+ev.ETmiss_vec).Mt(),tempWeight);
      m.at(STG).Fill(ev.pt1+ev.pt2+ev.selPhotons.at(0).p.Pt(),tempWeight);
      m.at(STMET).Fill(ev.pt1+ev.pt2+ev.selPhotons.at(0).p.Pt()+ev.ETmiss,tempWeight);
      m.at(MLLG).Fill((ev.l1+ev.l2+ev.selPhotons.at(0).vec).M(),tempWeight);
      
      if(FindGenPhotonMatch(ev.selPhotons.at(0))) m.at(genPhotonPT).Fill(GetGenPhotonMatch(ev.selPhotons.at(0)).p.Pt(),tempWeight);
      
      if(ev.evtHasGenPhotonVeto) if(FindGenPhotonMatch(ev.selPhotons.at(0))) m.at(genPhotonPT_Veto).Fill(GetGenPhotonMatch(ev.selPhotons.at(0)).p.Pt(),tempWeight);
      if(ev.evtHasGenPhotonVeto) m.at(PTG1_Veto).Fill(ev.selPhotons.at(0).p.Pt(),tempWeight);
      if(!ev.evtHasGenPhotonVeto) if(FindGenPhotonMatch(ev.selPhotons.at(0))) m.at(genPhotonPT_NoVeto).Fill(GetGenPhotonMatch(ev.selPhotons.at(0)).p.Pt(),tempWeight);
      if(!ev.evtHasGenPhotonVeto) m.at(PTG1_NoVeto).Fill(ev.selPhotons.at(0).p.Pt(),tempWeight);
      
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
   
   //float weight_TopPt=1.;
   //float weight_LeptonPt=1.;
   //float weight_nIsr=1.;
   //float weight_EWKinoPt=1.;
   
    //if(doWeights_TopPt){
      //weight_TopPt = topPtReweighting(*topPt1,*topPt2);
   //}
   //if(doWeights_nISR){
      //weight_nIsr = isrReweighting(*nISR,false);
   //}
   //if(doWeights_EWKinoPairPt){
      //weight_EWKinoPt = isrReweightingEWK(*EWKinoPairPt,false);
   //}
   //if(doWeights_LeptonPairPt){
      //weight_LeptonPt = isrReweightingEWK(*leptonPairPt,false);
   //}
     //
   //ev.totalWeight = ev.totalWeight*weight_TopPt*weight_LeptonPt*weight_nIsr*weight_EWKinoPt;
   //
   //m.at(WEIGHT_EWKINOPAIRPT).Fill(weight_EWKinoPt);
   //m.at(WEIGHT_LEPTONPAIRPT).Fill(weight_LeptonPt);
   //m.at(WEIGHT_NISR).Fill(weight_nIsr);
   //m.at(WEIGHT_TOPPT).Fill(weight_TopPt);
     
   m.at(ETMISS).Fill(TriggerBool,ev.ETmiss);
   m.at(PT1).Fill(TriggerBool,ev.pt1);
   m.at(PT2).Fill(TriggerBool,ev.pt2);
   m.at(MLL).Fill(TriggerBool,ev.mll);   
   m.at(NPHOTONS).Fill(TriggerBool,ev.selPhotons.size());
   m.at(NVTX).Fill(TriggerBool,*nGoodVertices);
   //m.at(HT).Fill(TriggerBool,*ht);
   m.at(HT).Fill(TriggerBool,ev.calcHt);
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

void HistogramProducer::FillerSignal(selEvent& ev, map<Histograms1D,TH1F>& m, float divideFactor=1.){
   
   float weight_TopPt=1.;
   float weight_LeptonPt=1.;
   float weight_nIsr=1.;
   float weight_EWKinoPt=1.;
   
   float tempWeight=ev.totalWeight;
   
   if(doWeights_TopPt){
      weight_TopPt = topPtReweighting(*topPt1,*topPt2);
   }
   if(doWeights_nISR){
      weight_nIsr = isrReweighting(*nISR,false);
   }
   if(doWeights_EWKinoPairPt){
      weight_EWKinoPt = isrReweightingEWK(*EWKinoPairPt,false);
   }
   if(doWeights_LeptonPairPt){
      weight_LeptonPt = isrReweightingEWK(*leptonPairPt,false);
   }
   
   tempWeight = tempWeight*weight_TopPt*weight_LeptonPt*weight_nIsr*weight_EWKinoPt;
   
   m.at(WEIGHT_EWKINOPAIRPT).Fill(weight_EWKinoPt);
   m.at(WEIGHT_LEPTONPAIRPT).Fill(weight_LeptonPt);
   m.at(WEIGHT_NISR).Fill(weight_nIsr);
   m.at(WEIGHT_TOPPT).Fill(weight_TopPt);
   
   
   
   //m.at(ETMISS).Fill(ev.ETmiss, ev.totalWeight*1./(nGen/divideFactor));
   m.at(ETMISS).Fill(ev.ETmiss, tempWeight*1./(nGen/divideFactor));
}

void HistogramProducer::FillHistograms(){

   //if(SelectEvent(UNCUT)&&config_selectionsToProcessMap[UNCUT]){
      //Filler(selectedEvent,h1Maps["uncut"],false); 
      //if (selectedEvent.isDiElectron){ 
         //Filler(selectedEvent,h1Maps["uncutEE"],false); 
      //}else{
         //Filler(selectedEvent,h1Maps["uncutMM"],false);      
      //}
   //}

   ClearVariables();

   if(config_selectionsToProcessMap[DILEP]){
   //cout<<"DILEP"<<endl;
      if(SelectEvent(DILEP)){
         if(selectedEvent.selLeptonSize==selectedEvent.matchedLeptonSize){
            if((selectedEvent.isDiMuon || selectedEvent.isDiElectron) && (!(selectedEvent.isElectronMuon || selectedEvent.isMuonElectron))) Filler(selectedEvent,h1Maps["dilep"],false);
            if (selectedEvent.isDiElectron && (!selectedEvent.isDiMuon) && (!(selectedEvent.isElectronMuon || selectedEvent.isMuonElectron))){
               Filler(selectedEvent,h1Maps["dilepEE"],false);
            }else{
              if (selectedEvent.isDiMuon && !(selectedEvent.isMuonElectron||selectedEvent.isElectronMuon)) {Filler(selectedEvent,h1Maps["dilepMM"],false);}
              if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron) Filler(selectedEvent,h1Maps["dilepEM"],false);
            }
         }
      }
   }

   ClearVariables();
   if(config_selectionsToProcessMap[PHOTON]){
      if(SelectEvent(PHOTON)){
         if(selectedEvent.selLeptonSize==selectedEvent.matchedLeptonSize){
            if (selectedEvent.selPhotons.size()!=0){
               if(selectedEvent.isDiMuon || selectedEvent.isDiElectron) Filler(selectedEvent,h1Maps["1photon"],true);         
               if (selectedEvent.isDiElectron){ 
                  Filler(selectedEvent,h1Maps["1photonEE"],true);
               }else{
                  if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["1photonMM"],true);
                  if(selectedEvent.isElectronMuon || selectedEvent.isMuonElectron)  Filler(selectedEvent,h1Maps["1photonEM"],true);
               }
            }
         }
      }
   }
   ClearVariables();
   if(config_selectionsToProcessMap[SEL]){
         //cout<<"SEL"<<endl;

      if(SelectEvent(SEL)){
         if(selectedEvent.selLeptonSize==selectedEvent.matchedLeptonSize){
            if (selectedEvent.selPhotons.size()!=0){
               if(selectedEvent.isDiMuon || selectedEvent.isDiElectron) Filler(selectedEvent,h1Maps["sel"],true);
               if (selectedEvent.isDiElectron){
                  Filler(selectedEvent,h1Maps["selEE"],true); 
               }else{
                  if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["selMM"],true); 
                  if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron) Filler(selectedEvent,h1Maps["selEM"],true); 
               }
            }
         }
      }
   }
   //ClearVariables();
   //clearCutFlowMap();
   //if(config_selectionsToProcessMap[ONZ]){
      //if(SelectEvent(ONZ)){
         //
         //int sum = selectedEvent.isDiElectron+selectedEvent.isDiMuon+selectedEvent.isMuonElectron+selectedEvent.isElectronMuon;
         //if(!(sum==1)){
         //cout<<sum<<selectedEvent.isDiElectron<<selectedEvent.isDiMuon<<selectedEvent.isMuonElectron<<selectedEvent.isElectronMuon<<endl;
         //}
           //
         //if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)){
            //if(selectedEvent.isDiMuon || selectedEvent.isDiElectron){
               //Filler(selectedEvent,h1Maps["onZ"],true);
            //}
            //if (selectedEvent.isDiElectron){
               //Filler(selectedEvent,h1Maps["onZEE"],true);
            //}else{
               //if (selectedEvent.isDiMuon){
                  //Filler(selectedEvent,h1Maps["onZMM"],true);
                  //}
               //if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron) Filler(selectedEvent,h1Maps["onZEM"],true);
               //
            //}
         //}
      //}
   //}

   //clearCutFlowMap();
      ClearVariables();

   if(config_selectionsToProcessMap[EXO]){
      if(SelectEvent(EXO)){
         if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>50. && selectedEvent.mll<130.)){
            if(selectedEvent.isDiMuon || selectedEvent.isDiElectron) Filler(selectedEvent,h1Maps["exo"],true);
            if (selectedEvent.isDiElectron){ 
               Filler(selectedEvent,h1Maps["exoEE"],true);
            }else{
               if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["exoMM"],true);
               if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron) Filler(selectedEvent,h1Maps["exoEM"],true);
            }
         }
      }
   }
      ClearVariables();

   if(config_selectionsToProcessMap[EGRegression]){
      if(SelectEvent(EGRegression)){
         if(selectedEvent.isDiMuon || selectedEvent.isDiElectron) Filler(selectedEvent,h1Maps["EGRegression"],false);
         if (selectedEvent.isDiElectron){ 
            Filler(selectedEvent,h1Maps["EGRegressionEE"],false);
         }else{
            if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["EGRegressionMM"],false);
            if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron) Filler(selectedEvent,h1Maps["EGRegressionEM"],false);
         }
      }
   }
      ClearVariables();

   if(config_selectionsToProcessMap[ONZMET]){
      if(SelectEvent(ONZ)){ 
            //cout<<"ONZMET"<<endl;

         if(selectedEvent.selLeptonSize==selectedEvent.matchedLeptonSize){

            if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)&&(selectedEvent.ETmiss<100.)){ //ONZ+MET<100
               //cout<<"ONZMET0100 "<<selectedEvent.ETmiss<<" "<<selectedEvent.totalWeight<<endl;
               if(selectedEvent.isDiMuon || selectedEvent.isDiElectron) Filler(selectedEvent,h1Maps["onZMet0100"],true);
               if (selectedEvent.isDiElectron){ 
                  Filler(selectedEvent,h1Maps["onZMet0100EE"],true);
               }else{
                  if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["onZMet0100MM"],true);
                  if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron) Filler(selectedEvent,h1Maps["onZMet0100EM"],true);
               }
            }
            
            if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)&&(selectedEvent.ETmiss>=100.)){
               //cout<<"ONZMET100 "<<selectedEvent.ETmiss<<" "<<selectedEvent.<<" "<<selectedEvent.totalWeight<<endl; //ONZ+MET>100
               if(selectedEvent.isDiMuon || selectedEvent.isDiElectron){
                  Filler(selectedEvent,h1Maps["onZMet100"],true);
                  }
               if (selectedEvent.isDiElectron){
                  Filler(selectedEvent,h1Maps["onZMet100EE"],true);
               }//else{
               if (selectedEvent.isDiMuon){
                     //cout<<"ONZMET100 "<<selectedEvent.ETmiss<<" "<<selectedEvent.pt1<<" "<<selectedEvent.pt2<<" "<<selectedEvent.totalWeight<<endl; //ONZ+MET>100
                     Filler(selectedEvent,h1Maps["onZMet100MM"],true);
               }
               if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron){
                   Filler(selectedEvent,h1Maps["onZMet100EM"],true);
               }
            }
            
            if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)&&(selectedEvent.ETmiss>=150.)){ //ONZ+MET>150
               //cout<<"ONZMET150 "<<selectedEvent.ETmiss<<" "<<selectedEvent.totalWeight<<endl;
               if(selectedEvent.isDiMuon || selectedEvent.isDiElectron){
                  Filler(selectedEvent,h1Maps["onZMet150"],true);
                  //if(*signal_m1==1500 && *signal_m2==400) dummy6+=1;
                  if(*signal_m1==600) dummy6+=1;
                  }
               if (selectedEvent.isDiElectron){
                  Filler(selectedEvent,h1Maps["onZMet150EE"],true);
                  //if(*signal_m1==1500 && *signal_m2==400) dummy7+=1;
                  if(*signal_m1==600) dummy7+=1;
               }//else{
               if (selectedEvent.isDiMuon){
                  //cout<<"ONZMET150 "<<selectedEvent.ETmiss<<" "<<selectedEvent.pt1<<" "<<selectedEvent.pt2<<" "<<selectedEvent.totalWeight<<endl; //ONZ+MET>100
                     Filler(selectedEvent,h1Maps["onZMet150MM"],true);
                     //if(*signal_m1==1500 && *signal_m2==400) dummy8+=1;
                     if(*signal_m1==600) dummy8+=1;
               }
               if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron){
                   Filler(selectedEvent,h1Maps["onZMet150EM"],true);
                   //if(*signal_m1==1500 && *signal_m2==400) dummy9+=1;
                   if(*signal_m1==600) dummy9+=1;
               }
               //if(*signal_m1==1500 && *signal_m2==400) dummy10+=1;
               if(*signal_m1==600) dummy10+=1;
            }
            
            if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)&&(selectedEvent.ETmiss<150.)&&(selectedEvent.ETmiss>=100.)){ //ONZ+MET>100<150
               if(selectedEvent.isDiMuon || selectedEvent.isDiElectron){
                  Filler(selectedEvent,h1Maps["onZMet100150"],true);
                  }
               if (selectedEvent.isDiElectron){
                  Filler(selectedEvent,h1Maps["onZMet100150EE"],true);
               }//else{
               if (selectedEvent.isDiMuon){
                     Filler(selectedEvent,h1Maps["onZMet100150MM"],true);
               }
               if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron){
                   Filler(selectedEvent,h1Maps["onZMet100150EM"],true);
               
               }
            }
         }
      }
   }
      ClearVariables();

   if(config_selectionsToProcessMap[ONZG]){
      if(SelectEvent(ONZG)){ //ONZG
           
         if ((selectedEvent.selPhotons.size()!=0)&&((selectedEvent.l1+selectedEvent.l2+selectedEvent.selPhotons.at(0).vec).M()>81.)&&((selectedEvent.l1+selectedEvent.l2+selectedEvent.selPhotons.at(0).vec).M()<101.)){
            if(selectedEvent.isDiMuon || selectedEvent.isDiElectron) Filler(selectedEvent,h1Maps["onZG"],true);
            if (selectedEvent.isDiElectron){ 
               Filler(selectedEvent,h1Maps["onZGEE"],true);
            }else{
               if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["onZGMM"],true);
               if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron) Filler(selectedEvent,h1Maps["onZGEM"],true);
            }
         }
      }
   }
      ClearVariables();

   if(config_selectionsToProcessMap[ABOVEZG]){
      if(SelectEvent(ONZ)){//MLLG>110
           
         if ((selectedEvent.selPhotons.size()!=0)&&((selectedEvent.l1+selectedEvent.l2+selectedEvent.selPhotons.at(0).vec).M()>110.)){
            if(selectedEvent.isDiMuon || selectedEvent.isDiElectron) Filler(selectedEvent,h1Maps["mllG110"],true);
            if (selectedEvent.isDiElectron){ 
               Filler(selectedEvent,h1Maps["mllG110EE"],true);
            }else{
               if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["mllG110MM"],true);
               if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron) Filler(selectedEvent,h1Maps["mllG110EM"],true);
            }
         }
      }
   }
   
   ClearVariables();
   clearCutFlowMap();
   if(config_selectionsToProcessMap[ControlRegionDY]){
         //cout<<"CRDY"<<endl;

      if(SelectEvent(ONZ)){
         if(selectedEvent.selLeptonSize==selectedEvent.matchedLeptonSize){
            if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)&&(selectedEvent.ETmiss<100.)){ //ONZ+MET<100 ~ CR DY/Z(+gamma)
               if(selectedEvent.isDiMuon || selectedEvent.isDiElectron) Filler(selectedEvent,h1Maps["CRDY"],true);
               if (selectedEvent.isDiElectron){
                  Filler(selectedEvent,h1Maps["CRDYEE"],true);
               }else{                          
                  if (selectedEvent.isDiMuon){Filler(selectedEvent,h1Maps["CRDYMM"],true);}
                  if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron) Filler(selectedEvent,h1Maps["CRDYEM"],true);
               }
            }
         }
      }
   }
   if(config_selectionsToProcessMap[ControlRegionWW]){
         //cout<<"CRWW"<<endl;

      if(SelectEvent(DILEP)){
         if(selectedEvent.selLeptonSize==selectedEvent.matchedLeptonSize){
            if ((selectedEvent.mll<81. || selectedEvent.mll>101.)&&(selectedEvent.ETmiss<100. && selectedEvent.ETmiss>20.)){ //ONZ+MET<100 ~ CR DY/Z(+gamma)
               if((selectedEvent.selJets.size()==0)&&((selectedEvent.l1+selectedEvent.l2).Pt()>45.)){
               
                  if(selectedEvent.isDiMuon || selectedEvent.isDiElectron || selectedEvent.isElectronMuon || selectedEvent.isMuonElectron){
                     Filler(selectedEvent,h1Maps["CRWWLL"],false);
                   }
                  if(selectedEvent.isDiMuon || selectedEvent.isDiElectron){
                     Filler(selectedEvent,h1Maps["CRWW"],false);
                   }
                  if (selectedEvent.isDiElectron){
                     Filler(selectedEvent,h1Maps["CRWWEE"],false);
                  }                      
                  if (selectedEvent.isDiMuon){
                     Filler(selectedEvent,h1Maps["CRWWMM"],false);
                  }
                  if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron){
                     Filler(selectedEvent,h1Maps["CRWWEM"],false);
                  }
               }
            }
         }
      }
   }
   ClearVariables();
   clearCutFlowMap();
   if(config_selectionsToProcessMap[ValidationRegion]){
         //cout<<"VR"<<endl;

      if(SelectEvent(ONZ)){
         if(selectedEvent.selLeptonSize==selectedEvent.matchedLeptonSize){
            if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)&&(selectedEvent.ETmiss>=100.)&&(selectedEvent.ETmiss<150.)){ //ONZ+MET<150>100 ~ VR
               if(selectedEvent.isDiMuon || selectedEvent.isDiElectron){
                  Filler(selectedEvent,h1Maps["VR"],true);
                  if(selectedEvent.selPhotons.at(0).p.Pt()>=80.){
                     Filler(selectedEvent,h1Maps["VR80"],true);
                  }else{
                     Filler(selectedEvent,h1Maps["VR080"],true);
                  }
               }
               if (selectedEvent.isDiElectron){
                  Filler(selectedEvent,h1Maps["VREE"],true);
                  if(selectedEvent.selPhotons.at(0).p.Pt()>=80.){
                     Filler(selectedEvent,h1Maps["VR80EE"],true);
                  }else{
                     Filler(selectedEvent,h1Maps["VR080EE"],true);
                  }
               }
               if (selectedEvent.isDiMuon){
                  Filler(selectedEvent,h1Maps["VRMM"],true);
                  if(selectedEvent.selPhotons.at(0).p.Pt()>=80.){
                     Filler(selectedEvent,h1Maps["VR80MM"],true);
                  }else{
                     Filler(selectedEvent,h1Maps["VR080MM"],true);
                  }
               }
               if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron){
                  Filler(selectedEvent,h1Maps["VREM"],true);
                  if(selectedEvent.selPhotons.at(0).p.Pt()>=80.){
                     Filler(selectedEvent,h1Maps["VR80EM"],true);
                  }else{
                     Filler(selectedEvent,h1Maps["VR080EM"],true);
                  }
               }
            }
         }
      }
   }
   ClearVariables();
   clearCutFlowMap();
   if(config_selectionsToProcessMap[ControlRegionTT]){
         //cout<<"CRTT"<<endl;

      if(SelectEvent(SEL)){
         if(selectedEvent.selLeptonSize==selectedEvent.matchedLeptonSize){
            if ((selectedEvent.selPhotons.size()!=0)){ //Different Flavor + 1Photon ~ CR TT(+gamma)
               if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron){
                  Filler(selectedEvent,h1Maps["CRTTEM"],true);
                  if(selectedEvent.selPhotons.at(0).p.Pt()>=80.){
                     Filler(selectedEvent,h1Maps["CRTT80EM"],true);
                  }else{
                     Filler(selectedEvent,h1Maps["CRTT080EM"],true);
                  }
               }
            }
         }
      }
   }
   ClearVariables();
   clearCutFlowMap();
   if(config_selectionsToProcessMap[ControlRegionZZ]){
         //cout<<"CRZZ"<<endl;

      if(SelectEvent(ControlRegionZZ)){
         if((selectedEvent.selLeptonSize==4)&&(selectedEvent.matchedLeptonSize==2)){
            if((selectedEvent.selMuonSize==4)||(selectedEvent.selElectronSize==4)||(selectedEvent.selMuonSize==2 && selectedEvent.selElectronSize==2)){

               int countNegCharge=0;
               int countPosCharge=0;
               
               vector<selMuon> negMuons;
               vector<selMuon> posMuons;
               vector<selElectron> negElectrons;
               vector<selElectron> posElectrons;
               
               float ZMass = 91.1876;
               
               float tempTotalWeight=selectedEvent.totalWeight;
               
               float addWeight=1.;
               
               for (vector<selElectron>::iterator it = selectedEvent.selElectrons.begin(); it != selectedEvent.selElectrons.end(); it++){
                  if (it->chargeInt <0){
                     countNegCharge+=1;
                     negElectrons.push_back(*it);
                  }else{
                     countPosCharge+=1;
                     posElectrons.push_back(*it);
                  }
                  if(it->matched==false){
                     tempTotalWeight=tempTotalWeight*GetScaleFactorAndErrorFinal(it->p.Pt(),it->p.Eta(),isSignal,true);
                  }
               }
               for (vector<selMuon>::iterator it = selectedEvent.selMuons.begin(); it != selectedEvent.selMuons.end(); it++){
                  if (it->chargeInt <0){
                     countNegCharge+=1;
                     negMuons.push_back(*it);
                  }else{
                     countPosCharge+=1;
                     posMuons.push_back(*it);
                  }
                  if(it->matched==false){
                     tempTotalWeight=tempTotalWeight*GetScaleFactorAndErrorFinal(it->p.Pt(),it->p.Eta(),isSignal,false);
                  }
               }
               
               selectedEvent.totalWeight=tempTotalWeight;
               
               if(countNegCharge==2 && countPosCharge==2){
               
                  selLepton lep1,lep2,lep3,lep4;
               
                  if(selectedEvent.selMuonSize==4){
                     
                     if((negMuons.size()==2)&&(posMuons.size()==2)){
                        selMuon neg1 = negMuons.at(0);
                        selMuon neg2 = negMuons.at(1);
                        selMuon pos1 = posMuons.at(0);
                        selMuon pos2 = posMuons.at(1);                  
                        
                        //there are 2 combinations n1p1/n2p2 - n1p2/n2p1
                        float n1p1 = abs((neg1.vec+pos1.vec).M()-ZMass);
                        float n1p2 = abs((neg1.vec+pos2.vec).M()-ZMass);
                        float n2p1 = abs((neg2.vec+pos1.vec).M()-ZMass);
                        float n2p2 = abs((neg2.vec+pos2.vec).M()-ZMass);
                        
                        if((n1p1<n1p2)&&(n1p1<n2p1)&&(n1p1<n2p2)){
                           if(neg1.p.Pt()>pos1.p.Pt()){
                              lep1.setAll(neg1);
                              lep2.setAll(pos1);
                           }else{
                              lep1.setAll(pos1);
                              lep2.setAll(neg1);
                           }
                           if(neg2.p.Pt()>pos2.p.Pt()){
                              lep3.setAll(neg2);
                              lep4.setAll(pos2);
                           }else{
                              lep3.setAll(pos2);
                              lep4.setAll(neg2);
                           }
                        }else{
                           if((n1p2<n1p1)&&(n1p2<n2p1)&&(n1p2<n2p2)){
                              if(neg1.p.Pt()>pos2.p.Pt()){
                                 lep1.setAll(neg1);
                                 lep2.setAll(pos2);
                              }else{
                                 lep1.setAll(pos2);
                                 lep2.setAll(neg1);
                              }
                              if(neg2.p.Pt()>pos1.p.Pt()){
                                 lep3.setAll(neg2);
                                 lep4.setAll(pos1);
                              }else{
                                 lep3.setAll(pos1);
                                 lep4.setAll(neg2);
                              }
                           }else{
                              if((n2p1<n1p1)&&(n2p1<n1p2)&&(n2p1<n2p2)){
                                 if(neg2.p.Pt()>pos1.p.Pt()){
                                    lep1.setAll(neg2);
                                    lep2.setAll(pos1);
                                 }else{
                                    lep1.setAll(pos1);
                                    lep2.setAll(neg2);
                                 }
                                 if(neg1.p.Pt()>pos2.p.Pt()){
                                    lep3.setAll(neg1);
                                    lep4.setAll(pos2);
                                 }else{
                                    lep3.setAll(pos2);
                                    lep4.setAll(neg1);
                                 }
                              }else{
                                 if(neg2.p.Pt()>pos2.p.Pt()){
                                    lep1.setAll(neg2);
                                    lep2.setAll(pos2);
                                 }else{
                                    lep1.setAll(pos2);
                                    lep2.setAll(neg2);
                                 }
                                 if(neg1.p.Pt()>pos1.p.Pt()){
                                    lep3.setAll(neg1);
                                    lep4.setAll(pos1);
                                 }else{
                                    lep3.setAll(pos1);
                                    lep4.setAll(neg1);
                                 }
                              }
                           }
                        }
                    }
                    
                     float mll1=(lep1.vec+lep2.vec).M();
                     float mll2=(lep3.vec+lep4.vec).M();
                     if(mll1<106. && mll1>76. && mll2<130. && mll2>50.){
                        FillerZZ(selectedEvent,h1Maps["CRZZMM"],false,lep1,lep2,lep3,lep4);
                     }
                  }
                  
                  if(selectedEvent.selElectronSize==4){
                     
                     if((negElectrons.size()==2)&&(posElectrons.size()==2)){
                        selElectron neg1 = negElectrons.at(0);
                        selElectron neg2 = negElectrons.at(1);
                        selElectron pos1 = posElectrons.at(0);
                        selElectron pos2 = posElectrons.at(1);                  
                        
                        //there are 2 combinations n1p1/n2p2 - n1p2/n2p1
                        float n1p1 = abs((neg1.vec+pos1.vec).M()-ZMass);
                        float n1p2 = abs((neg1.vec+pos2.vec).M()-ZMass);
                        float n2p1 = abs((neg2.vec+pos1.vec).M()-ZMass);
                        float n2p2 = abs((neg2.vec+pos2.vec).M()-ZMass);
                        
                        
                        if((n1p1<n1p2)&&(n1p1<n2p1)&&(n1p1<n2p2)){
                           if(neg1.p.Pt()>pos1.p.Pt()){
                              lep1.setAll(neg1);
                              lep2.setAll(pos1);
                           }else{
                              lep1.setAll(pos1);
                              lep2.setAll(neg1);
                           }
                           if(neg2.p.Pt()>pos2.p.Pt()){
                              lep3.setAll(neg2);
                              lep4.setAll(pos2);
                           }else{
                              lep3.setAll(pos2);
                              lep4.setAll(neg2);
                           }
                        }else{
                           if((n1p2<n1p1)&&(n1p2<n2p1)&&(n1p2<n2p2)){
                              if(neg1.p.Pt()>pos2.p.Pt()){
                                 lep1.setAll(neg1);
                                 lep2.setAll(pos2);
                              }else{
                                 lep1.setAll(pos2);
                                 lep2.setAll(neg1);
                              }
                              if(neg2.p.Pt()>pos1.p.Pt()){
                                 lep3.setAll(neg2);
                                 lep4.setAll(pos1);
                              }else{
                                 lep3.setAll(pos1);
                                 lep4.setAll(neg2);
                              }
                           }else{
                              if((n2p1<n1p1)&&(n2p1<n1p2)&&(n2p1<n2p2)){
                                 if(neg2.p.Pt()>pos1.p.Pt()){
                                    lep1.setAll(neg2);
                                    lep2.setAll(pos1);
                                 }else{
                                    lep1.setAll(pos1);
                                    lep2.setAll(neg2);
                                 }
                                 if(neg1.p.Pt()>pos2.p.Pt()){
                                    lep3.setAll(neg1);
                                    lep4.setAll(pos2);
                                 }else{
                                    lep3.setAll(pos2);
                                    lep4.setAll(neg1);
                                 }
                              }else{
                                 if(neg2.p.Pt()>pos2.p.Pt()){
                                    lep1.setAll(neg2);
                                    lep2.setAll(pos2);
                                 }else{
                                    lep1.setAll(pos2);
                                    lep2.setAll(neg2);
                                 }
                                 if(neg1.p.Pt()>pos1.p.Pt()){
                                    lep3.setAll(neg1);
                                    lep4.setAll(pos1);
                                 }else{
                                    lep3.setAll(pos1);
                                    lep4.setAll(neg1);
                                 }
                              }
                           }
                        }
                     }
                     float mll1=(lep1.vec+lep2.vec).M();
                     float mll2=(lep3.vec+lep4.vec).M();
                     if(mll1<106. && mll1>76. && mll2<130. && mll2>50.){
                        FillerZZ(selectedEvent,h1Maps["CRZZEE"],false,lep1,lep2,lep3,lep4);
                     }
                     
                  }
                  
                  if(selectedEvent.selMuonSize==2 && selectedEvent.selElectronSize==2){
                     
                     if((negMuons.size()==1)&&(posMuons.size()==1)&&(negElectrons.size()==1)&&(posElectrons.size()==1)){
                     
                        selMuon m1 = negMuons.at(0);
                        selMuon m2 = posMuons.at(0);
                        selElectron e1 = negElectrons.at(0);
                        selElectron e2 = posElectrons.at(0);                  
                        
                        float m1m2 = abs((m1.vec+m2.vec).M()-ZMass);
                        float e1e2 = abs((e1.vec+e2.vec).M()-ZMass);
                        
                        if(m1m2<e1e2){
                           if(m1.p.Pt()>m2.p.Pt()){
                              lep1.setAll(m1);
                              lep2.setAll(m2);
                           }else{
                              lep1.setAll(m2);
                              lep2.setAll(m1);
                           }
                           if(e1.p.Pt()>e2.p.Pt()){
                              lep3.setAll(e1);
                              lep4.setAll(e2);
                           }else{
                              lep3.setAll(e2);
                              lep4.setAll(e1);
                           }
                           
                           float mll1=(lep1.vec+lep2.vec).M();
                           float mll2=(lep3.vec+lep4.vec).M();
                           if(mll1<106. && mll1>76. && mll2<130. && mll2>50.){
                              FillerZZ(selectedEvent,h1Maps["CRZZMM"],false,lep1,lep2,lep3,lep4);
                           }
                           
                        }else{
                           if(m1.p.Pt()>m2.p.Pt()){
                              lep3.setAll(m1);
                              lep4.setAll(m2);
                           }else{
                              lep3.setAll(m2);
                              lep4.setAll(m1);
                           }
                           if(e1.p.Pt()>e2.p.Pt()){
                              lep1.setAll(e1);
                              lep2.setAll(e2);
                           }else{
                              lep1.setAll(e2);
                              lep2.setAll(e1);
                           }
                           
                           float mll1=(lep1.vec+lep2.vec).M();
                           float mll2=(lep3.vec+lep4.vec).M();
                           if(mll1<106. && mll1>76. && mll2<130. && mll2>50.){
                              FillerZZ(selectedEvent,h1Maps["CRZZEE"],false,lep1,lep2,lep3,lep4);
                           }
                        }
                     }
                  }
                  float mll1=(lep1.vec+lep2.vec).M();
                  float mll2=(lep3.vec+lep4.vec).M();
                  if(mll1<106. && mll1>76. && mll2<130. && mll2>50.){
                     FillerZZ(selectedEvent,h1Maps["CRZZ"],false,lep1,lep2,lep3,lep4);
                  }
               }
            }
            
            //if ((selectedEvent.mll>61. && selectedEvent.mll<121.)){
               //if (selectedEvent.isDiMuon || selectedEvent.isDiElectron) Filler(selectedEvent,h1Maps["ControlRegionZZ"],false);
               //if (selectedEvent.isDiElectron) Filler(selectedEvent,h1Maps["ControlRegionZZEE"],false);                         
               //if (selectedEvent.isDiMuon)Filler(selectedEvent,h1Maps["ControlRegionZZMM"],false);
               //if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron) Filler(selectedEvent,h1Maps["ControlRegionZZEM"],false);
            //}
         }
      }
   }
   
ClearVariables();
   clearCutFlowMap();
   if(config_selectionsToProcessMap[ControlRegionWZ]){
         //cout<<"CRWZ"<<endl;

      //if(SelectEvent(ControlRegionWZ)){
      if(SelectEvent(ControlRegionZZ)){
         if((selectedEvent.selLeptonSize==3)&&(selectedEvent.matchedLeptonSize==2)){
            if((selectedEvent.selMuonSize==3)||(selectedEvent.selElectronSize==3)||(selectedEvent.selMuonSize==1 && selectedEvent.selElectronSize==2)||(selectedEvent.selMuonSize==2 && selectedEvent.selElectronSize==1)){

               int countNegCharge=0;
               int countPosCharge=0;
               
               vector<selMuon> negMuons;
               vector<selMuon> posMuons;
               vector<selElectron> negElectrons;
               vector<selElectron> posElectrons;
               
               float ZMass = 91.1876;
               
               float tempTotalWeight=selectedEvent.totalWeight;
               for (vector<selElectron>::iterator it = selectedEvent.selElectrons.begin(); it != selectedEvent.selElectrons.end(); it++){
                  if (it->chargeInt <0){
                     countNegCharge+=1;
                     negElectrons.push_back(*it);
                  }else{
                     countPosCharge+=1;
                     posElectrons.push_back(*it);
                  }
                  if(it->matched==false){
                     tempTotalWeight=tempTotalWeight*GetScaleFactorAndErrorFinal(it->p.Pt(),it->p.Eta(),isSignal,true);
                  }
               }
               for (vector<selMuon>::iterator it = selectedEvent.selMuons.begin(); it != selectedEvent.selMuons.end(); it++){
                  if (it->chargeInt <0){
                     countNegCharge+=1;
                     negMuons.push_back(*it);
                  }else{
                     countPosCharge+=1;
                     posMuons.push_back(*it);
                  }
                  if(it->matched==false){
                     tempTotalWeight=tempTotalWeight*GetScaleFactorAndErrorFinal(it->p.Pt(),it->p.Eta(),isSignal,false);
                  }
               }
               
               selectedEvent.totalWeight=tempTotalWeight;
               
               if((countNegCharge==1 && countPosCharge==2)||(countNegCharge==2 && countPosCharge==1)){
                  //find lepton combinations
               
               selLepton lep1,lep2,lep3;
               
                  if(selectedEvent.selMuonSize==3){
                     
                     if((negMuons.size()==2)&&(posMuons.size()==1)){
                     
                        selMuon neg1 = negMuons.at(0);
                        selMuon neg2 = negMuons.at(1);
                        selMuon pos1 = posMuons.at(0);
                        
                        //there are 2 combinations n1p1/n2 - n2p1/n1
                        float n1p1 = abs((neg1.vec+pos1.vec).M()-ZMass);
                        float n2p1 = abs((neg2.vec+pos1.vec).M()-ZMass);

                        if(n1p1<n2p1){
                           if(neg1.p.Pt()>pos1.p.Pt()){
                              lep1.setAll(neg1);
                              lep2.setAll(pos1);
                           }else{
                              lep1.setAll(pos1);
                              lep2.setAll(neg1);
                           }
                           lep3.setAll(neg2);
                        }else{
                           if(n2p1<n1p1){
                              if(neg2.p.Pt()>pos1.p.Pt()){
                                 lep1.setAll(neg2);
                                 lep2.setAll(pos1);
                              }else{
                                 lep1.setAll(pos1);
                                 lep2.setAll(neg2);
                              }
                              lep3.setAll(neg1);
                           }
                        }
                     }
                     if((negMuons.size()==1)&&(posMuons.size()==2)){
                     
                        selMuon neg1 = negMuons.at(0);
                        selMuon pos1 = posMuons.at(0);
                        selMuon pos2 = posMuons.at(1);
                        
                        //there are 2 combinations n1p1/p2 - n1p2/p1
                        float n1p1 = abs((neg1.vec+pos1.vec).M()-ZMass);
                        float n1p2 = abs((neg1.vec+pos2.vec).M()-ZMass);

                        if(n1p1<n1p2){
                           if(neg1.p.Pt()>pos1.p.Pt()){
                              lep1.setAll(neg1);
                              lep2.setAll(pos1);
                           }else{
                              lep1.setAll(pos1);
                              lep2.setAll(neg1);
                           }
                           lep3.setAll(pos2);
                        }else{
                           if(n1p2<n1p1){
                              if(neg1.p.Pt()>pos2.p.Pt()){
                                 lep1.setAll(neg1);
                                 lep2.setAll(pos2);
                              }else{
                                 lep1.setAll(pos2);
                                 lep2.setAll(neg1);
                              }
                              lep3.setAll(pos1);
                           }
                        }
                     }
                     
                     float mll1=(lep1.vec+lep2.vec).M();
                     if(mll1<106. && mll1>76.){
                        if( (lep1.vec+selectedEvent.ETmiss_vec).Mt() > 50. ){
                           if(selectedEvent.ETmiss>70.){
                              FillerWZ(selectedEvent,h1Maps["CRWZMM"],false,lep1,lep2,lep3);
                           }
                        }
                     }
                     
                  }
                  
                  
                  if(selectedEvent.selElectronSize==3){
                     
                     if((negElectrons.size()==2)&&(posElectrons.size()==1)){
                     
                        selElectron neg1 = negElectrons.at(0);
                        selElectron neg2 = negElectrons.at(1);
                        selElectron pos1 = posElectrons.at(0);
                        
                        //there are 2 combinations n1p1/n2 - n2p1/n1
                        float n1p1 = abs((neg1.vec+pos1.vec).M()-ZMass);
                        float n2p1 = abs((neg2.vec+pos1.vec).M()-ZMass);

                        if(n1p1<n2p1){
                           if(neg1.p.Pt()>pos1.p.Pt()){
                              lep1.setAll(neg1);
                              lep2.setAll(pos1);
                           }else{
                              lep1.setAll(pos1);
                              lep2.setAll(neg1);
                           }
                           lep3.setAll(neg2);
                        }else{
                           if(n2p1<n1p1){
                              if(neg2.p.Pt()>pos1.p.Pt()){
                                 lep1.setAll(neg2);
                                 lep2.setAll(pos1);
                              }else{
                                 lep1.setAll(pos1);
                                 lep2.setAll(neg2);
                              }
                              lep3.setAll(neg1);
                           }
                        }
                     }
                     if((negElectrons.size()==1)&&(posElectrons.size()==2)){
                     
                        selElectron neg1 = negElectrons.at(0);
                        selElectron pos1 = posElectrons.at(0);
                        selElectron pos2 = posElectrons.at(1);
                        
                        //there are 2 combinations n1p1/p2 - n1p2/p1
                        float n1p1 = abs((neg1.vec+pos1.vec).M()-ZMass);
                        float n1p2 = abs((neg1.vec+pos2.vec).M()-ZMass);

                        if(n1p1<n1p2){
                           if(neg1.p.Pt()>pos1.p.Pt()){
                              lep1.setAll(neg1);
                              lep2.setAll(pos1);
                           }else{
                              lep1.setAll(pos1);
                              lep2.setAll(neg1);
                           }
                           lep3.setAll(pos2);
                        }else{
                           if(n1p2<n1p1){
                              if(neg1.p.Pt()>pos2.p.Pt()){
                                 lep1.setAll(neg1);
                                 lep2.setAll(pos2);
                              }else{
                                 lep1.setAll(pos2);
                                 lep2.setAll(neg1);
                              }
                              lep3.setAll(pos1);
                           }
                        }
                     }
                     
                     
                     float mll1=(lep1.vec+lep2.vec).M();
                     if(mll1<106. && mll1>76.){
                        if( (lep1.vec+selectedEvent.ETmiss_vec).Mt() > 50. ){
                           if(selectedEvent.ETmiss>70.){
                              FillerWZ(selectedEvent,h1Maps["CRWZEE"],false,lep1,lep2,lep3);
                           }
                        }
                     }
                     
                  }
                  
                  if((selectedEvent.selElectronSize==2)&&(selectedEvent.selMuonSize==1)){
                  
                     if(negElectrons.size()==1 && posElectrons.size()==1 && posMuons.size()==1){
                        selElectron neg1 = negElectrons.at(0);
                        selElectron pos1 = posElectrons.at(0);
                        selMuon pos2 = posMuons.at(0);
                        
                        if(neg1.p.Pt()>pos1.p.Pt()){
                           lep1.setAll(neg1);
                           lep2.setAll(pos1);
                        }else{
                           lep1.setAll(pos1);
                           lep2.setAll(neg1);
                        }
                        lep3.setAll(pos2);
                     }
                     if(negElectrons.size()==1 && posElectrons.size()==1 && negMuons.size()==1){
                        selElectron neg1 = negElectrons.at(0);
                        selElectron pos1 = posElectrons.at(0);
                        selMuon neg2 = negMuons.at(0);
                        
                        if(neg1.p.Pt()>pos1.p.Pt()){
                           lep1.setAll(neg1);
                           lep2.setAll(pos1);
                        }else{
                           lep1.setAll(pos1);
                           lep2.setAll(neg1);
                        }
                        lep3.setAll(neg2);
                     }    
                     
                     float mll1=(lep1.vec+lep2.vec).M();
                     if(mll1<106. && mll1>76.){
                        if( (lep1.vec+selectedEvent.ETmiss_vec).Mt() > 50. ){
                           if(selectedEvent.ETmiss>70.){
                              FillerWZ(selectedEvent,h1Maps["CRWZEE"],false,lep1,lep2,lep3);
                           }
                        }
                     }
                     
                  }
                  
                  if((selectedEvent.selElectronSize==1)&&(selectedEvent.selMuonSize==2)){
                  
                     if(negMuons.size()==1 && posMuons.size()==1 && posElectrons.size()==1){
                        selMuon neg1 = negMuons.at(0);
                        selMuon pos1 = posMuons.at(0);
                        selElectron pos2 = posElectrons.at(0);
                        
                        if(neg1.p.Pt()>pos1.p.Pt()){
                           lep1.setAll(neg1);
                           lep2.setAll(pos1);
                        }else{
                           lep1.setAll(pos1);
                           lep2.setAll(neg1);
                        }
                        lep3.setAll(pos2);
                     }
                     if(negMuons.size()==1 && posMuons.size()==1 && negElectrons.size()==1){
                        selMuon neg1 = negMuons.at(0);
                        selMuon pos1 = posMuons.at(0);
                        selElectron neg2 = negElectrons.at(0);
                        
                        if(neg1.p.Pt()>pos1.p.Pt()){
                           lep1.setAll(neg1);
                           lep2.setAll(pos1);
                        }else{
                           lep1.setAll(pos1);
                           lep2.setAll(neg1);
                        }
                        lep3.setAll(neg2);
                     }   
                     
                     float mll1=(lep1.vec+lep2.vec).M();
                     if(mll1<106. && mll1>76.){
                        if( (lep1.vec+selectedEvent.ETmiss_vec).Mt() > 50. ){
                           if(selectedEvent.ETmiss>70.){
                              FillerWZ(selectedEvent,h1Maps["CRWZMM"],false,lep1,lep2,lep3);
                           }
                        }
                     }                  
                  }
                  
                               
                  
               float mll1=(lep1.vec+lep2.vec).M();
               if(mll1<106. && mll1>76.){
                  if( (lep1.vec+selectedEvent.ETmiss_vec).Mt() > 50. ){
                     if(selectedEvent.ETmiss>70.){
                        FillerWZ(selectedEvent,h1Maps["CRWZ"],false,lep1,lep2,lep3);
                     }
                  }
               }
               }
            }
            
         }
      }
   }
   
   ClearVariables();
   clearCutFlowMap();
   if(config_selectionsToProcessMap[ONZ]){
         //cout<<"ONZ"<<endl;

      if(SelectEvent(ONZ)){
            if(selectedEvent.selLeptonSize==selectedEvent.matchedLeptonSize){
            if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)){
               //cout<<"ONZ "<<selectedEvent.ETmiss<<" "<<selectedEvent.totalWeight<<endl;
               if(selectedEvent.isDiMuon || selectedEvent.isDiElectron){
                  Filler(selectedEvent,h1Maps["onZ"],true);
               }
               if (selectedEvent.isDiElectron){
                  Filler(selectedEvent,h1Maps["onZEE"],true);
               }else{
                  if (selectedEvent.isDiMuon){
                     //cout<<"ONZ "<<selectedEvent.ETmiss<<" "<<selectedEvent.pt1<<" "<<selectedEvent.pt2<<" "<<selectedEvent.totalWeight<<endl; //ONZ+MET>100
                     Filler(selectedEvent,h1Maps["onZMM"],true);
                     }
                  if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron) Filler(selectedEvent,h1Maps["onZEM"],true);
                  
               }
            }
         }
      }
   }
   ClearVariables();
}

void HistogramProducer::FillSignalHistograms(){
   ClearVariables();
   if(SelectEvent(ONZ)){
      if(selectedEvent.selLeptonSize==selectedEvent.matchedLeptonSize){
         if((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)&&(selectedEvent.ETmiss>=150.)){
            float tempPhotonLeadingPt = selectedEvent.selPhotons.at(0).p.Pt();
            string signalPoint;
            if(inputName.find("GGM")!=string::npos){
               signalPoint = getSignalPointName(10,*signal_m1,*signal_m2);
            }else{
               signalPoint = getSignalPointName(*nBinos,*signal_m1,*signal_m2);
            }
            if(!(s1Maps.count(signalPoint)>0)){
               InitSignalScanHistos(signalPoint);
            }
            if(selectedEvent.isDiMuon || selectedEvent.isDiElectron){
               FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig"));
               if(tempPhotonLeadingPt>=80.){
                  FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig80"));
               }else{
                  FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig080"));
               }
            }
            if (selectedEvent.isDiElectron){
               FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sigEE"));
               if(tempPhotonLeadingPt>=80.){
                  FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig80EE"));
               }else{
                  FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig080EE"));
               }
            }
            if (selectedEvent.isDiMuon){
               FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sigMM"));
               if(tempPhotonLeadingPt>=80.){
                  FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig80MM"));
               }else{
                  FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig080MM"));
               }
            }
            if(selectedEvent.isElectronMuon || selectedEvent.isMuonElectron){
               FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sigEM"));
               if(tempPhotonLeadingPt>=80.){
                  FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig80EM"));
               }else{
                  FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig080EM"));
               }
            }
            
            if(config_dosignalscanSplit){
               
               if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2){
                  if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 22) {
                     if (fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 22) {
                        if(selectedEvent.isDiMuon || selectedEvent.isDiElectron){
                           //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_TChiNG_gg"),4.0);
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gg"),4.0);
                           if(tempPhotonLeadingPt>=80.){
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gg80"),4.0);
                           }else{
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gg080"),4.0);
                           }
                        }
                        if (selectedEvent.isDiElectron){
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_ggEE"),4.0);
                           if(tempPhotonLeadingPt>=80.){
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gg80EE"),4.0);
                           }else{
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gg080EE"),4.0);
                           }
                        }
                        if (selectedEvent.isDiMuon){
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_ggMM"),4.0);
                           if(tempPhotonLeadingPt>=80.){
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gg80MM"),4.0);
                           }else{
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gg080MM"),4.0);
                           }
                        }
                        if(selectedEvent.isElectronMuon || selectedEvent.isMuonElectron){
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_ggEM"),4.0);
                           if(tempPhotonLeadingPt>=80.){
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gg80EM"),4.0);
                           }else{
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gg080EM"),4.0);
                           }
                        }
                     }
                  }
               }
               if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2){
                  if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 23) {
                     if (fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 23) {
                        if(selectedEvent.isDiMuon || selectedEvent.isDiElectron){
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zz"),16.0);
                           if(tempPhotonLeadingPt>=80.){
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zz80"),16.0);
                           }else{
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zz080"),16.0);
                           }
                        }
                        if (selectedEvent.isDiElectron){
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zzEE"),16.0);
                           if(tempPhotonLeadingPt>=80.){
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zz80EE"),16.0);
                           }else{
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zz080EE"),16.0);
                           }
                        }
                        if (selectedEvent.isDiMuon){
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zzMM"),16.0);
                           if(tempPhotonLeadingPt>=80.){
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zz80MM"),16.0);
                           }else{
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zz080MM"),16.0);
                           }
                        }
                        if(selectedEvent.isElectronMuon || selectedEvent.isMuonElectron){
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zzEM"),16.0);
                           if(tempPhotonLeadingPt>=80.){
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zz80EM"),16.0);
                           }else{
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zz080EM"),16.0);
                           }
                        }
                     }
                  }
               }
               if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2){
                  if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 23) {
                     if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 22) {
                        if(selectedEvent.isDiMuon || selectedEvent.isDiElectron){
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gz"),4.0);
                           if(tempPhotonLeadingPt>=80.){
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gz80"),4.0);
                           }else{
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gz080"),4.0);
                           }
                        }
                        if (selectedEvent.isDiElectron){
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gzEE"),4.0);
                           if(tempPhotonLeadingPt>=80.){
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gz80EE"),4.0);
                           }else{
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gz080EE"),4.0);
                           }
                        }
                        if (selectedEvent.isDiMuon){
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gzMM"),4.0);
                           if(tempPhotonLeadingPt>=80.){
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gz80MM"),4.0);
                           }else{
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gz080MM"),4.0);
                           }
                        }
                        if(selectedEvent.isElectronMuon || selectedEvent.isMuonElectron){
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gzEM"),4.0);
                           if(tempPhotonLeadingPt>=80.){
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gz80EM"),4.0);
                           }else{
                              FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gz080EM"),4.0);
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



void HistogramProducer::FillCompressedTree(){
   
   comp_weight = selectedEvent.totalWeight*(1./nGen);
   
   comp_Lepton1Pt = selectedEvent.pt1;
   comp_Lepton1Px = selectedEvent.l1.Px();
   comp_Lepton1Py = selectedEvent.l1.Py();
   comp_Lepton1Pz = selectedEvent.l1.Pz();
   comp_Lepton1Eta = selectedEvent.eta1;
   comp_Lepton1Phi = selectedEvent.phi1;
   comp_Lepton1MiniIso = selectedEvent.miniIso1;
   comp_Lepton1IsElectron = (selectedEvent.isDiElectron||selectedEvent.isElectronMuon);
   comp_Lepton1IsMuon = (selectedEvent.isDiMuon||selectedEvent.isMuonElectron);

   comp_Lepton2Pt = selectedEvent.pt1;
   comp_Lepton2Px = selectedEvent.l2.Px();
   comp_Lepton2Py = selectedEvent.l2.Py();
   comp_Lepton2Pz = selectedEvent.l2.Pz();
   comp_Lepton2Eta = selectedEvent.eta2;
   comp_Lepton2Phi = selectedEvent.phi2;
   comp_Lepton2MiniIso = selectedEvent.miniIso2;
   comp_Lepton2IsElectron = (selectedEvent.isDiElectron||selectedEvent.isMuonElectron);
   comp_Lepton2IsMuon = (selectedEvent.isDiMuon||selectedEvent.isElectronMuon);

   comp_deltaRll = selectedEvent.deltaRll;  
   comp_mll = selectedEvent.mll;

   comp_PhotonPt.clear();
   comp_PhotonPx.clear();
   comp_PhotonPy.clear();
   comp_PhotonPz.clear();
   comp_PhotonEta.clear();
   comp_PhotonPhi.clear();

   for (vector<selPhoton>::iterator it = selectedEvent.selPhotons.begin(); it != selectedEvent.selPhotons.end(); it++){
      comp_PhotonPt.push_back(it->vec.Pt());
      comp_PhotonPx.push_back(it->vec.Px());
      comp_PhotonPy.push_back(it->vec.Py());
      comp_PhotonPz.push_back(it->vec.Pz());
      comp_PhotonEta.push_back(it->vec.Eta());
      comp_PhotonPhi.push_back(it->vec.Phi());
   }
   comp_NPhoton = selectedEvent.selPhotons.size();

   comp_Photon1Pt = comp_PhotonPt.at(0);
   comp_Photon1Px = comp_PhotonPx.at(0);
   comp_Photon1Py = comp_PhotonPy.at(0);
   comp_Photon1Pz = comp_PhotonPz.at(0);
   comp_Photon1Eta = comp_PhotonEta.at(0);
   comp_Photon1Phi = comp_PhotonPhi.at(0);


   comp_JetPt.clear();
   comp_JetPx.clear();
   comp_JetPy.clear();
   comp_JetPz.clear();
   comp_JetEta.clear();
   comp_JetPhi.clear();
   comp_JetBDiscriminator.clear();
   comp_JetChf.clear();
   comp_JetNhf.clear();
   comp_JetNConstituents.clear();

   for (vector<selJet>::iterator it = selectedEvent.selJets.begin(); it != selectedEvent.selJets.end(); it++){
      comp_JetPt.push_back(it->vec.Pt());
      comp_JetPx.push_back(it->vec.Px());
      comp_JetPy.push_back(it->vec.Py());
      comp_JetPz.push_back(it->vec.Pz());
      comp_JetEta.push_back(it->vec.Eta());
      comp_JetPhi.push_back(it->vec.Phi());
      comp_JetBDiscriminator.push_back(it->bDiscriminator);
      comp_JetChf.push_back(it->chf);
      //comp_JetNhf.push_back(it->nhf);
      //comp_JetNConstituents.push_back(it->nconstituents);
   }
   comp_NJet = selectedEvent.selJets.size();

   comp_MetPt = selectedEvent.ETmiss_vec.Pt();
   comp_MetPx = selectedEvent.ETmiss_vec.Px();
   comp_MetPy = selectedEvent.ETmiss_vec.Py();
   comp_MetPz = selectedEvent.ETmiss_vec.Pz();
   comp_MetEta = selectedEvent.ETmiss_vec.Eta();
   comp_MetPhi = selectedEvent.ETmiss_vec.Phi();
   
   comp_ht = *ht;
   
   comp_m1=*signal_m1;
   comp_m2=*signal_m2;
   
   comp_eventNo=*evtNo;
   
   newtreeCompressed->Fill();
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
   if(m[DIELECTRON] && !(m[DIMUON]||m[EMUON])){
      
      //if(m[genZLL]) c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("genZToLL",mW[genZLL]);
      
      if (m[TRIGGERED]){
      //if (m[genZLL]){
         c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("triggered",mW[TRIGGERED]);
         //c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("triggered",mW[genZLL]);
         if(m[GENVETO]){
            c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("GenPhotonVeto",mW[GENVETO]);
            if(m[LEPTONIDPure_leading] && m[LEPTONIDPure_trailing]){
               c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("2LeptonIDPure",mW[LEPTONIDPure_leading]);
               if(m[LEPTONIDImpact_leading] && m[LEPTONIDImpact_trailing]){
                  c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("2LeptonIDImpact",mW[LEPTONIDImpact_leading]);
                  if(m[LEPTONIDEta_leading] && m[LEPTONIDEta_trailing]){
                     c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("2LeptonIDEta",mW[LEPTONIDEta_leading]);
                     if(m[LEPTONIDIso_leading] && m[LEPTONIDIso_trailing]){
                        c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("2LeptonIDIso",mW[LEPTONIDIso_leading]);
                        if(m[LEPTONIDDeltaR_leading] && m[LEPTONIDDeltaR_trailing]){
                           c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("2LeptonIDDeltaR",mW[LEPTONIDDeltaR_leading]);
                           if(m[LEPTONPT_leading] && m[LEPTONPT_trailing]){
                              c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("2LeptonPT",mW[LEPTONPT_leading]);
                              if(m[M50]){
                                 c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("m50",mW[M50]);
                                 if(m[PHOTON1]){
                                    c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("1Photon",mW[PHOTON1]);
                                    if(m[PHOTON1ID]){
                                       c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("1PhotonID",mW[PHOTON1ID]);
                                       if(m[PHOTON1SEED]){
                                          c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("1PhotonSeed",mW[PHOTON1SEED]);
                                          if(m[PHOTON1ETA]){
                                             c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("1PhotonEta",mW[PHOTON1ETA]);
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
                  }
               }
            }
         }
      }
   }else{
   if(m[DIMUON] && !(m[DIELECTRON]||m[EMUON])){
      //if(m[genZLL])  c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("genZToLL",mW[genZLL]);
      if (m[TRIGGERED]){
         c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("triggered",mW[TRIGGERED]);
      //if (m[genZLL]){
         //c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("triggered",mW[genZLL]);
         if(m[GENVETO]){
            c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("GenPhotonVeto",mW[GENVETO]);
            if(m[LEPTONIDPure_leading] && m[LEPTONIDPure_trailing]){
               c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("2LeptonIDPure",mW[LEPTONIDPure_leading]);
               if(m[LEPTONIDImpact_leading] && m[LEPTONIDImpact_trailing]){
                  c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("2LeptonIDImpact",mW[LEPTONIDImpact_leading]);
                  if(m[LEPTONIDEta_leading] && m[LEPTONIDEta_trailing]){
                     c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("2LeptonIDEta",mW[LEPTONIDEta_leading]);
                     if(m[LEPTONIDIso_leading] && m[LEPTONIDIso_trailing]){
                        c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("2LeptonIDIso",mW[LEPTONIDIso_leading]);
                        if(m[LEPTONIDDeltaR_leading] && m[LEPTONIDDeltaR_trailing]){
                           c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("2LeptonIDDeltaR",mW[LEPTONIDDeltaR_leading]);
                           if(m[LEPTONPT_leading] && m[LEPTONPT_trailing]){
                              c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("2LeptonPT",mW[LEPTONPT_leading]);
                              if(m[M50]){
                                 c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("m50",mW[M50]);
                                 if(m[PHOTON1]){
                                    c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("1Photon",mW[PHOTON1]);
                                    if(m[PHOTON1ID]){
                                       c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("1PhotonID",mW[PHOTON1ID]);
                                       if(m[PHOTON1SEED]){
                                          c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("1PhotonSeed",mW[PHOTON1SEED]);
                                          if(m[PHOTON1ETA]){
                                             c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("1PhotonEta",mW[PHOTON1ETA]);
                                             if(m[PHOTON1PT]){
                                                c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("1PhotonPT",mW[PHOTON1PT]);
                                                if(m[PHOTON1DR]){
                                                   c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("1PhotonDeltaR",mW[PHOTON1DR]);
                                                   if(m[ZMASS]){
                                                      c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("Z",mW[ZMASS]);
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
               }
            }
         }
      }
   }
   if(m[EMUON] && !(m[DIELECTRON]||m[DIMUON])){
      if(m[genZLL])  c1Maps["cutFlow_Fine_onZEM"].at(CUTFLOW_fine).Fill("genZToLL",mW[genZLL]);
      if (m[TRIGGERED]){
         c1Maps["cutFlow_Fine_onZEM"].at(CUTFLOW_fine).Fill("triggered",mW[TRIGGERED]);
      //if (m[genZLL]){
         //c1Maps["cutFlow_Fine_onZEM"].at(CUTFLOW_fine).Fill("triggered",mW[genZLL]);
         if(m[GENVETO]){
            c1Maps["cutFlow_Fine_onZEM"].at(CUTFLOW_fine).Fill("GenPhotonVeto",mW[GENVETO]);
            if(m[LEPTONIDPure_leading] && m[LEPTONIDPure_trailing]){
               c1Maps["cutFlow_Fine_onZEM"].at(CUTFLOW_fine).Fill("2LeptonIDPure",mW[LEPTONIDPure_leading]);
               if(m[LEPTONIDImpact_leading] && m[LEPTONIDImpact_trailing]){
                  c1Maps["cutFlow_Fine_onZEM"].at(CUTFLOW_fine).Fill("2LeptonIDImpact",mW[LEPTONIDImpact_leading]);
                  if(m[LEPTONIDEta_leading] && m[LEPTONIDEta_trailing]){
                     c1Maps["cutFlow_Fine_onZEM"].at(CUTFLOW_fine).Fill("2LeptonIDEta",mW[LEPTONIDEta_leading]);
                     if(m[LEPTONIDIso_leading] && m[LEPTONIDIso_trailing]){
                        c1Maps["cutFlow_Fine_onZEM"].at(CUTFLOW_fine).Fill("2LeptonIDIso",mW[LEPTONIDIso_leading]);
                        if(m[LEPTONIDDeltaR_leading] && m[LEPTONIDDeltaR_trailing]){
                           c1Maps["cutFlow_Fine_onZEM"].at(CUTFLOW_fine).Fill("2LeptonIDDeltaR",mW[LEPTONIDDeltaR_leading]);
                           if(m[LEPTONPT_leading] && m[LEPTONPT_trailing]){
                              c1Maps["cutFlow_Fine_onZEM"].at(CUTFLOW_fine).Fill("2LeptonPT",mW[LEPTONPT_leading]);
                              if(m[M50]){
                                 c1Maps["cutFlow_Fine_onZEM"].at(CUTFLOW_fine).Fill("m50",mW[M50]);
                                 if(m[PHOTON1]){
                                    c1Maps["cutFlow_Fine_onZEM"].at(CUTFLOW_fine).Fill("1Photon",mW[PHOTON1]);
                                    if(m[PHOTON1ID]){
                                       c1Maps["cutFlow_Fine_onZEM"].at(CUTFLOW_fine).Fill("1PhotonID",mW[PHOTON1ID]);
                                       if(m[PHOTON1SEED]){
                                          c1Maps["cutFlow_Fine_onZEM"].at(CUTFLOW_fine).Fill("1PhotonSeed",mW[PHOTON1SEED]);
                                          if(m[PHOTON1ETA]){
                                             c1Maps["cutFlow_Fine_onZEM"].at(CUTFLOW_fine).Fill("1PhotonEta",mW[PHOTON1ETA]);
                                             if(m[PHOTON1PT]){
                                                c1Maps["cutFlow_Fine_onZEM"].at(CUTFLOW_fine).Fill("1PhotonPT",mW[PHOTON1PT]);
                                                if(m[PHOTON1DR]){
                                                   c1Maps["cutFlow_Fine_onZEM"].at(CUTFLOW_fine).Fill("1PhotonDeltaR",mW[PHOTON1DR]);
                                                   if(m[ZMASS]){
                                                      c1Maps["cutFlow_Fine_onZEM"].at(CUTFLOW_fine).Fill("Z",mW[ZMASS]);
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
               }
            }
         }
      }
   }
}

}



void HistogramProducer::FillTriggerStudies(){
   ClearVariables();
   if(SelectEventTriggerStudies(TRIGDILEP)){
      if (selectedEvent.trigHt){//baselineTrigger
         if(selectedEvent.isDiElectron && !(selectedEvent.isDiMuon || selectedEvent.isMuonElectron || selectedEvent.isMuonElectron)){
            FillerTrigger(selectedEvent,eff1Maps["trigDilepEE"],false,selectedEvent.trigDiEle);
            }
         if(selectedEvent.isDiMuon && !(selectedEvent.isDiElectron || selectedEvent.isMuonElectron || selectedEvent.isMuonElectron)){
            FillerTrigger(selectedEvent,eff1Maps["trigDilepMM"],false,selectedEvent.trigDiMu);
            }
         if((selectedEvent.isMuonElectron||selectedEvent.isElectronMuon) && !(selectedEvent.isDiElectron || selectedEvent.isDiMuon)){
            FillerTrigger(selectedEvent,eff1Maps["trigDilepEM"],false,selectedEvent.trigMuEle);
            }
      }
   }
   ClearVariables();
   if(SelectEventTriggerStudies(TRIGDILEP_ptcuts)){
      if (selectedEvent.trigHt){//baselineTrigger
         if(selectedEvent.isDiElectron && !(selectedEvent.isDiMuon || selectedEvent.isMuonElectron || selectedEvent.isMuonElectron)){
            FillerTrigger(selectedEvent,eff1Maps["trigDilepEE_ptcuts"],false,selectedEvent.trigDiEle);
         }
         if(selectedEvent.isDiMuon && !(selectedEvent.isDiElectron || selectedEvent.isMuonElectron || selectedEvent.isMuonElectron)){
            FillerTrigger(selectedEvent,eff1Maps["trigDilepMM_ptcuts"],false,selectedEvent.trigDiMu);
         }
         if((selectedEvent.isMuonElectron||selectedEvent.isElectronMuon) && !(selectedEvent.isDiElectron || selectedEvent.isDiMuon)){
            FillerTrigger(selectedEvent,eff1Maps["trigDilepEM_ptcuts"],false,selectedEvent.trigMuEle);
         }
      }
   }
   //ClearVariables();
   //if(SelectEventTriggerStudies(TRIGDILEP_pt1cut)){
      //if (selectedEvent.trigHt){//baselineTrigger
         //if (selectedEvent.isDiElectron){ 
           //FillerTrigger(selectedEvent,eff1Maps["trigDilepEE_pt1cut"],false,selectedEvent.trigDiEle);
         //}else{
            //FillerTrigger(selectedEvent,eff1Maps["trigDilepMM_pt1cut"],false,selectedEvent.trigDiMu);  
         //}
      //}
   //}
   //ClearVariables();
   //if(SelectEventTriggerStudies(TRIGDILEP_pt2cut)){
      //if (selectedEvent.trigHt){//baselineTrigger
         //if (selectedEvent.isDiElectron){ 
           //FillerTrigger(selectedEvent,eff1Maps["trigDilepEE_pt2cut"],false,selectedEvent.trigDiEle);
         //}else{
            //FillerTrigger(selectedEvent,eff1Maps["trigDilepMM_pt2cut"],false,selectedEvent.trigDiMu);  
         //}
      //}
   //}
   //ClearVariables();
   //if(SelectEventTriggerStudies(TRIGSEL)){
      //if (selectedEvent.selPhotons.size()!=0){
         //if (selectedEvent.trigHt){//baselineTrigger
            //if (selectedEvent.isDiElectron){
               //FillerTrigger(selectedEvent,eff1Maps["trigSelEE"],true,selectedEvent.trigDiEle); 
            //}else{
               //FillerTrigger(selectedEvent,eff1Maps["trigSelMM"],true,selectedEvent.trigDiMu); 
            //}
         //}
      //}
   //}
   //ClearVariables();
   //if(SelectEventTriggerStudies(TRIGSEL_ptcuts)){
      //if (selectedEvent.selPhotons.size()!=0){
         //if (selectedEvent.trigHt){//baselineTrigger
            //if (selectedEvent.isDiElectron){
               //FillerTrigger(selectedEvent,eff1Maps["trigSelEE_ptcuts"],true,selectedEvent.trigDiEle); 
            //}else{
               //FillerTrigger(selectedEvent,eff1Maps["trigSelMM_ptcuts"],true,selectedEvent.trigDiMu); 
            //}
         //}
      //}
   //}
   //ClearVariables();
   //if(SelectEventTriggerStudies(TRIGONZ)){
      //if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81 && selectedEvent.mll<101)){
         //if (selectedEvent.trigHt){//baselineTrigger
            //if (selectedEvent.isDiElectron){
               //FillerTrigger(selectedEvent,eff1Maps["trigOnZEE"],true,selectedEvent.trigDiEle);  
            //}else{
               //FillerTrigger(selectedEvent,eff1Maps["trigOnZMM"],true,selectedEvent.trigDiMu); 
            //}
         //}
      //}
   //}   
   //ClearVariables();
   //if(SelectEventTriggerStudies(TRIGONZ_ptcuts)){
      //if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81 && selectedEvent.mll<101)){
         //if (selectedEvent.trigHt){//baselineTrigger
            //if (selectedEvent.isDiElectron){
               //FillerTrigger(selectedEvent,eff1Maps["trigOnZEE_ptcuts"],true,selectedEvent.trigDiEle);  
            //}else{
               //FillerTrigger(selectedEvent,eff1Maps["trigOnZMM_ptcuts"],true,selectedEvent.trigDiMu); 
            //}
         //}
      //}
   //}
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

bool HistogramProducer::FindRecoPhotonMatch(const tree::GenParticle pa){
   bool recoMatch = false;
   for(auto& reco : *photons) {
      if(reco.p.DeltaR(pa.p)<0.2){
         recoMatch=true;
         break;
      }
   }
   return recoMatch;
}


bool HistogramProducer::FindGenPhotonMatch(const tree::Photon& pa){
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

tree::GenParticle& HistogramProducer::GetGenPhotonMatch(const tree::Photon& pa){
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


bool HistogramProducer::matchLepton(const selElectron& pa, selEvent& ev){
   bool foundMatch=false;
   if((pa.deltaR1 < 0.1 ) || (pa.deltaR2 < 0.1 )){
      foundMatch=true;
   }else{
      foundMatch=false;
   }
   return foundMatch;
}
bool HistogramProducer::matchLepton(const selMuon& pa, selEvent& ev){
   bool foundMatch=false;
   TLorentzVector dummyLepton(0.,0.,0.,0.);
   dummyLepton.SetPtEtaPhiM(pa.p.Pt(),pa.p.Eta(),pa.p.Phi(),0.51099895e-3);
   if((pa.deltaR1 < 0.1 ) || (pa.deltaR2 < 0.1 )){
      foundMatch=true;
   }else{
      foundMatch=false;
   }
   return foundMatch;
}


Bool_t HistogramProducer::Process(Long64_t entry){
   
   float tempPercentage = (float) entry/ (float)nEntries;
   if(!(abs(config_eventpercentage-100.)<0.1)){
      if (tempPercentage>config_eventpercentage/100.){
         return kTRUE;
      }
   }

   

   fReader.SetLocalEntry(entry);
   //CorrectAllMuonPt();

   totalWeight = *mc_weight * *pu_weight;
   

   //cout<<pdf_weights->size()<<" "<<inputName<<endl;
   //cout<<"--------"<<endl;
   //cout<<entry<<" | "<<tempPercentage<<" | "<<*evtNo<<endl;


   int progress = tempPercentage*100.;
   if(entry%100000==0){
   
		std::cout<<"[";
		for(int i=0;i<100;i++)
			if(i<progress)
				std::cout<<'=';
			else if(i==progress)
				std::cout<<'>';
			else
				std::cout<<' ';
		std::cout<<"] "<<progress<<" %"<<" "<<getOutputFilename(inputName)<<'\r';
		std::cout.flush();
	
   }
   //
   //if (!(*signal_m1==1500 && *signal_m2==400)){
      //return kTRUE;
   //}
   //if (!(*signal_m1==600)){
      //return kTRUE;
   //}


   //myMuons.clear();
//
   //CorrectAllMuonPt();



   string cutFlowName = "TreeWriter/hCutFlow";
   
   if(isTotalSignal){
      //if(inputName.find("TChiNG")!=string::npos){
      if(inputName.find("TC")!=string::npos){
         cutFlowName+="TChiNG";
      }else{
         //if(inputName.find("T5bbbbZg")!=string::npos){
         if(inputName.find("T5")!=string::npos){
            cutFlowName+="T5bbbbZg";
         }else{
            if(inputName.find("GGM")!=string::npos){
               cutFlowName+="GGM";
            }else{
               //if(inputName.find("T6ttZg")!=string::npos){
               if(inputName.find("T6")!=string::npos){
                  cutFlowName+="T6ttZg";
               }
            }
         }
      }
      if(inputName.find("GGM")!=string::npos){
         cutFlowName += "_M1"+to_string(*signal_m1);
      }else{
         cutFlowName += "_"+to_string(*signal_m1);
      }
      if (*signal_m2){
         if(inputName.find("GGM")!=string::npos){
            //if(inputName.find("M1-200to1500_M2")!=string::npos) cutFlowName += "_M2"+to_string(*signal_m2);
            //if(inputName.find("M1-50to1500_M3")!=string::npos) cutFlowName += "_M3"+to_string(*signal_m2);
            if(inputName.find("M1-2")!=string::npos) cutFlowName += "_M2"+to_string(*signal_m2);
            if(inputName.find("M1-5")!=string::npos) cutFlowName += "_M3"+to_string(*signal_m2);
         }else{ 
            cutFlowName += "_"+to_string(*signal_m2);
         }
      }
      cutFlow = *((TH1F*)fReader.GetTree()->GetCurrentFile()->Get(cutFlowName.c_str()));
      if (cutFlow.GetEntries()) {
         nGen = cutFlow.GetBinContent(2);
      }else{
            cout << "Could not read cutFlow histogram " << cutFlowName << endl;
      }
   }else{
      cutFlow = *((TH1F*)fReader.GetTree()->GetCurrentFile()->Get(cutFlowName.c_str()));
      if(cutFlow.GetEntries()){
         nGen = cutFlow.GetBinContent(2);
      }
   }

   if (!isTotalSignal){
      SetCutFlowHistogramsStatus();
      clearCutFlowMap();
      ClearVariables();
   }
   if (!isTotalSignal){
      FillHistograms();
   }
   if (config_docutflow) FillCutFlowHistograms();
   if (config_docutflowfine) FillCutFlowHistograms_Fine();
   if (!isTotalSignal){
      clearCutFlowMap();
   }
   //FillHistograms2D();
   
   //if(SelectEvent(SEL)){
      //if(selectedEvent.selPhotons.size()!=0){
         //FillCompressedTree();
      //}
   //}
   
   if((inputName.find("JetHT")!= string::npos)||(!isData && !isTotalSignal)){
      FillTriggerStudies();
   }
   
   
   //if (config_docutflow) FillCutFlowHistograms();
   //if (config_docutflowfine) FillCutFlowHistograms_Fine();
   //clearCutFlowMap();
   if (config_dosignalscan){
       FillSignalHistograms();
   }
   if (!isTotalSignal){
      clearCutFlowMap();
   }

   //if(SelectEvent(SEL)){
      //if(selectedEvent.selPhotons.size()!=0){
         //FillCompressedTree();
      //}
   //}

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


void HistogramProducer::SaveCompressedTree(){
   outputPathCompressed = "/user/swuchterl/nTuple_compressed";
   
   outputFilenameCompressed = getOutputFilename(inputName);
   
   auto saveStr= outputPathCompressed+"/"+outputFilenameCompressed;
   TFile newfileCompressed(saveStr.c_str(),"RECREATE");
   newtreeCompressed->Write();
   
   newfileCompressed.Close();
}



void HistogramProducer::Terminate()
{
   auto outputName = "output"+config_outputfolder+"/"+getOutputFilename(inputName);
   
   //SaveCompressedTree();
   //cout<<endl;
   //cout<<"signal LL "<<dummy1<<endl;
   //cout<<"signal EE "<<dummy2<<endl;
   //cout<<"signal MM "<<dummy3<<endl;
   //cout<<"signal EM "<<dummy4<<endl;
   //cout<<"signal xx "<<dummy5<<endl;
   //cout<<"normal LL "<<dummy6<<endl;
   //cout<<"normal EE "<<dummy7<<endl;
   //cout<<"normal MM "<<dummy8<<endl;
   //cout<<"normal EM "<<dummy9<<endl;
   //cout<<"normal xx "<<dummy10<<endl;
   
   
   //cout<<dummy13<<endl;
   //cout<<dummy14<<endl;
   //cout<<dummy15<<endl;
   //cout<<dummy16<<endl;
   
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

