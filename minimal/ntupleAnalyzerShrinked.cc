#include "ntupleAnalyzerShrinked.h"

ntupleAnalyzerShrinked::ntupleAnalyzerShrinked():
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
   pdf_weights(fReader, "pdf_weights"),
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
   
   //Trigger Objects
   //EE
   trigObj_ele17_ele12(fReader,"hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter"),
   trigObj_ele23_ele12(fReader,"hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter"),
   trigObj_ele33_ele33(fReader,"hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter"),
   trigObj_ele33_ele33_mw(fReader,"hltDiEle33CaloIdLGsfTrkIdVLMWPMS2UnseededFilter"),
   
   trigObj_mu17_mu8(fReader,"hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4"),
   trigObj_mu17_mu8tk(fReader,"hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4"),
   trigObj_mu17_mu8_dz(fReader,"hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2"),
   trigObj_mu17_mu8tk_dz(fReader,"hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2"),
   trigObj_mu17tk_mu8tk_dz(fReader,"hltDiMuonTrk17Trk8RelTrkIsoFiltered0p4DzFiltered0p2"),
   trigObj_mu27_mu8(fReader,"hltDiMuonGlb27Trk8DzFiltered0p2"),
   trigObj_mu30_mu11(fReader,"hltDiMuonGlb30Trk11DzFiltered0p2"),
   
   trigObj_mu17_ele12_eleLeg(fReader,"hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"),
   trigObj_mu17_ele12_muLeg(fReader,"hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered17"),
   
   trigObj_mu23_ele8_eleLeg(fReader,"hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"),
   trigObj_mu23_ele8_muLeg(fReader,"hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23"), 
     
   trigObj_mu23_ele8_dz(fReader,"hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLDZFilter"),  
    
   trigObj_mu23_ele12_eleLeg(fReader,"hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"),   
   trigObj_mu23_ele12_muLeg(fReader,"hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23"),   
   
   trigObj_mu23_ele12_dz(fReader,"hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLDZFilter"),    
   
   trigObj_mu8_ele17_eleLeg(fReader,"hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"),   
   trigObj_mu8_ele17_muLeg(fReader,"hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8"),   
   
   trigObj_mu8_ele23_eleLeg(fReader,"hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"),   
   trigObj_mu8_ele23_muLeg(fReader,"hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8"),  
    
   trigObj_mu8_ele23_dz(fReader,"hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLDZFilter"),   
   
   trigObj_mu12_ele23_muLeg(fReader,"hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"),   
   trigObj_mu12_ele23_eleLeg(fReader,"hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered12"),   
   
   trigObj_mu12_ele23_dz(fReader,"hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLDZFilter"),   
   
   trigObj_mu30_ele30_eleLeg(fReader,"hltEle30CaloIdLGsfTrkIdVLDPhiUnseededFilter"),   
   trigObj_mu30_ele30_muLeg(fReader,"hltL3fL1sMu22orMu25orMu20EG15orMu5EG20L1f0L2f10QL3Filtered30Q"),   
   
   trigObj_mu33_ele33_eleLeg(fReader,"hltEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter"),   
   trigObj_mu33_ele33_muLeg(fReader,"hltL3fL1sMu22orMu25orMu20EG15orMu5EG20L1f0L2f10QL3Filtered33Q"),   
   
   
   
   signal_m1(fReader,"signal_m1"),
   signal_m2(fReader,"signal_m2"),
   nBinos(fReader,"signal_nBinos"),
   startTime(time(NULL))//,
   //rand()
{
}


void ntupleAnalyzerShrinked::Init(TTree *tree)
{
   fReader.SetTree(tree);
   inputName = fReader.GetTree()->GetCurrentFile()->GetName();
   
   isTotalSignal = (inputName.find("SMS-TChiNG_BF") != string::npos) || (inputName.find("GMSB") != string::npos) || (inputName.find("SMS-T6ttZg_nTuple")!=string::npos) || (inputName.find("SMS-T5bbbbZg_nTuple")!=string::npos) || (inputName.find("GGM")!=string::npos);
   isData = inputName.find("Run201") != string::npos;
   isSignal = (inputName.find("SMS") != string::npos) || (inputName.find("GGM") != string::npos) || (inputName.find("GMSB") != string::npos);
   
   boost::property_tree::ini_parser::read_ini("example.ini", propertyTree);

   doWeights_TopPt = (inputName.find("TTTo") != string::npos) || (inputName.find("TTGamma") != string::npos) || (inputName.find("ST_") != string::npos);
   doWeights_EWKinoPairPt = (inputName.find("TChi") != string::npos);
   doWeights_LeptonPairPt = false;
   doWeights_nISR = (inputName.find("T5") != string::npos);

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
   //config_selectionsToProcessMap[ControlRegionDY] = false;
   config_selectionsToProcessMap[ControlRegionTT] = propertyTree.get<bool>("selectionsToProcess.control");
   //config_selectionsToProcessMap[ControlRegionTT] = false;
   config_selectionsToProcessMap[ControlRegionZZ] = propertyTree.get<bool>("selectionsToProcess.control");
   //config_selectionsToProcessMap[ControlRegionZZ] = false;
   config_selectionsToProcessMap[ControlRegionWZ] = propertyTree.get<bool>("selectionsToProcess.control");
   config_selectionsToProcessMap[ControlRegionWZ] = false;
   //config_selectionsToProcessMap[ControlRegionWW] = propertyTree.get<bool>("selectionsToProcess.control");
   config_selectionsToProcessMap[ControlRegionWW] = false;
   config_selectionsToProcessMap[ValidationRegion] = propertyTree.get<bool>("selectionsToProcess.validation");
   
   config_veto = propertyTree.get<int>("generalSettings.veto");
   config_docutflow = propertyTree.get<bool>("generalSettings.docutflow");
   config_docutflowfine = propertyTree.get<bool>("generalSettings.docutflowfine");
   config_dosignalscan = propertyTree.get<bool>("generalSettings.dosignalscan");
   config_dosignalscanSplit = propertyTree.get<bool>("generalSettings.dosignalscantchingsplit");
   config_dosignalscanSplit = config_dosignalscanSplit && config_dosignalscan && inputName.find("TChiNG") != string::npos;
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
   fReader.GetEntries(true);
  
   nEntries=fReader.GetTree()->GetEntries();
  
  setHistoNames();
  
   if(isSignal){
      //if(inputName.find("GGM")!=string::npos){
      //nWeights=0;
      //}else{
      nWeights=9;
      //}
   }else{
      //if(inputName.find("ST_")!=string::npos){
         //nWeights=0;
      //}else{
         nWeights=isData? 0 : 110;
         //nWeights=0;
         //nWeights=110;

      //}
   }
  
   noPromptPhotons = inputName.find("DYJets") != string::npos;
   isZGammaInclusive = (inputName.find("ZGTo2LG_ext") != string::npos) || (inputName.find("ZGTo2LG_Total") != string::npos);

   if (config_docutflow) InitCutFlowHistos();
   if (config_docutflowfine) InitCutFlowHistos_Fine();
   if(!isTotalSignal){
      InitAllHistos();
      InitTriggerStudiesHistos();
   }
   InitScaleFactorsFinal();
}


void ntupleAnalyzerShrinked::SlaveBegin(TTree *tree)
{
}



Bool_t ntupleAnalyzerShrinked::Process(Long64_t entry){
   
   float tempPercentage = (float) entry/ (float)nEntries;
   if(!(abs(config_eventpercentage-100.)<0.1)){
      if (tempPercentage>config_eventpercentage/100.){
         return kTRUE;
      }
   }

   

   fReader.SetLocalEntry(entry);


   //for(int i=0; i<pdf_weights->size(); i++){
      //cout<<i<<" "<<pdf_weights->at(i)<<endl;
   //}


  //if(isSignal){
   //if(inputName.find("GGM")!=string::npos){
      //nWeights=0;
   //}else{
      //nWeights=9;
   //}
  //}else{
   //nWeights=isData? 0 : pdf_weights->size();
  //}

   //cout<<nWeights<<endl;

   totalWeight = *mc_weight * *pu_weight;
   
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


   string cutFlowName = "TreeWriter/hCutFlow";
   
   if(isTotalSignal){
      if(inputName.find("TChi")!=string::npos){
         cutFlowName+="TChiNG";
      }else{
         if(inputName.find("T5")!=string::npos){
            cutFlowName+="T5bbbbZg";
         }else{
            if(inputName.find("GGM")!=string::npos){
               cutFlowName+="GGM";
            }else{
               if(inputName.find("GMSB")!=string::npos){
                  cutFlowName+="GMSB";
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
   //ClearVariables();
   //cleaned=Cleaning();

   if (!isTotalSignal){
      SetCutFlowHistogramsStatus();
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
   
   
   if((inputName.find("JetHT")!= string::npos)||(!isData && !isTotalSignal)){
      FillTriggerStudies();
   }

   if (config_dosignalscan){
       FillSignalHistograms();
   }
   if (!isTotalSignal){
      clearCutFlowMap();
   }

   return kTRUE;
}

void ntupleAnalyzerShrinked::FillHistograms2D(){
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

void ntupleAnalyzerShrinked::FillSignalHistograms(){
   ClearVariables();
   if(SelectEvent(ONZ)){
      if(selectedEvent.selLeptonSize==selectedEvent.matchedLeptonSize){
         if((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)&&(selectedEvent.ETmiss>=150.)){
            //float tempPhotonLeadingPt = selectedEvent.selPhotons.at(0).p.Pt();
            string signalPoint;
            if(inputName.find("GGM")!=string::npos){
               signalPoint = getSignalPointName(10,*signal_m1,*signal_m2);
            }else{
               if(inputName.find("GMSB")!=string::npos){
                  signalPoint = getSignalPointName(11,*signal_m1,*signal_m2);
               }else{
                  signalPoint = getSignalPointName(*nBinos,*signal_m1,*signal_m2);
               }
            }
            //cout<<signalPoint<<endl;
            if(!(s1Maps.count(signalPoint)>0)){
               InitSignalScanHistos(signalPoint);
            }
            if(selectedEvent.isDiMuon || selectedEvent.isDiElectron){
               FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig"));
               //if(tempPhotonLeadingPt>=80.){
                  //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig80"));
               //}else{
                  //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig080"));
               //}
            }
            if (selectedEvent.isDiElectron){
               FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sigEE"));
               //if(tempPhotonLeadingPt>=80.){
                  //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig80EE"));
               //}else{
                  //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig080EE"));
               //}
            }
            if (selectedEvent.isDiMuon){
               FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sigMM"));
               //if(tempPhotonLeadingPt>=80.){
                  //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig80MM"));
               //}else{
                  //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig080MM"));
               //}
            }
            if(selectedEvent.isElectronMuon || selectedEvent.isMuonElectron){
               FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sigEM"));
               //if(tempPhotonLeadingPt>=80.){
                  //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig80EM"));
               //}else{
                  //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig080EM"));
               //}
            }
            
            if(config_dosignalscanSplit){
               
               if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2){
                  if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 22) {
                     if (fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 22) {
                        if(selectedEvent.isDiMuon || selectedEvent.isDiElectron){
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gg"),4.0);
                           //if(tempPhotonLeadingPt>=80.){
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gg80"),4.0);
                           //}else{
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gg080"),4.0);
                           //}
                        }
                        if (selectedEvent.isDiElectron){
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_ggEE"),4.0);
                           //if(tempPhotonLeadingPt>=80.){
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gg80EE"),4.0);
                           //}else{
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gg080EE"),4.0);
                           //}
                        }
                        if (selectedEvent.isDiMuon){
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_ggMM"),4.0);
                           //if(tempPhotonLeadingPt>=80.){
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gg80MM"),4.0);
                           //}else{
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gg080MM"),4.0);
                           //}
                        }
                        if(selectedEvent.isElectronMuon || selectedEvent.isMuonElectron){
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_ggEM"),4.0);
                           //if(tempPhotonLeadingPt>=80.){
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gg80EM"),4.0);
                           //}else{
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gg080EM"),4.0);
                           //}
                        }
                     }
                  }
               }
               if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2){
                  if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 23) {
                     if (fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 23) {
                        if(selectedEvent.isDiMuon || selectedEvent.isDiElectron){
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zz"),16.0);
                           //if(tempPhotonLeadingPt>=80.){
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zz80"),16.0);
                           //}else{
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zz080"),16.0);
                           //}
                        }
                        if (selectedEvent.isDiElectron){
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zzEE"),16.0);
                           //if(tempPhotonLeadingPt>=80.){
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zz80EE"),16.0);
                           //}else{
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zz080EE"),16.0);
                           //}
                        }
                        if (selectedEvent.isDiMuon){
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zzMM"),16.0);
                           //if(tempPhotonLeadingPt>=80.){
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zz80MM"),16.0);
                           //}else{
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zz080MM"),16.0);
                           //}
                        }
                        if(selectedEvent.isElectronMuon || selectedEvent.isMuonElectron){
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zzEM"),16.0);
                           //if(tempPhotonLeadingPt>=80.){
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zz80EM"),16.0);
                           //}else{
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_zz080EM"),16.0);
                           //}
                        }
                     }
                  }
               }
               if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2){
                  if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 23) {
                     if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 22) {
                        if(selectedEvent.isDiMuon || selectedEvent.isDiElectron){
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gz"),4.0);
                           //if(tempPhotonLeadingPt>=80.){
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gz80"),4.0);
                           //}else{
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gz080"),4.0);
                           //}
                        }
                        if (selectedEvent.isDiElectron){
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gzEE"),4.0);
                           //if(tempPhotonLeadingPt>=80.){
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gz80EE"),4.0);
                           //}else{
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gz080EE"),4.0);
                           //}
                        }
                        if (selectedEvent.isDiMuon){
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gzMM"),4.0);
                           //if(tempPhotonLeadingPt>=80.){
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gz80MM"),4.0);
                           //}else{
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gz080MM"),4.0);
                           //}
                        }
                        if(selectedEvent.isElectronMuon || selectedEvent.isMuonElectron){
                           FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gzEM"),4.0);
                           //if(tempPhotonLeadingPt>=80.){
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gz80EM"),4.0);
                           //}else{
                              //FillerSignal(selectedEvent,s1Maps.at(signalPoint).at("sig_gz080EM"),4.0);
                           //}
                        }
                     }
                  }
               }
            }
         }
      }
   }
}


void ntupleAnalyzerShrinked::FillTriggerStudies(){
   ClearVariables();
   if(SelectEventTriggerStudies(TRIGDILEP)){
      if (selectedEvent.trigHt){//baselineTrigger
         if(selectedEvent.isDiElectron && !(selectedEvent.isDiMuon || selectedEvent.isMuonElectron || selectedEvent.isMuonElectron)){
            FillerTrigger(selectedEvent,eff1Maps["trigDilepEE"],false,(selectedEvent.trigDiEle&&selectedEvent.trigDiEleMatch));
            }
         if(selectedEvent.isDiMuon && !(selectedEvent.isDiElectron || selectedEvent.isMuonElectron || selectedEvent.isMuonElectron)){
            FillerTrigger(selectedEvent,eff1Maps["trigDilepMM"],false,(selectedEvent.trigDiMu&&selectedEvent.trigDiMuMatch));
            }
         if((selectedEvent.isMuonElectron||selectedEvent.isElectronMuon) && !(selectedEvent.isDiElectron || selectedEvent.isDiMuon)){
            FillerTrigger(selectedEvent,eff1Maps["trigDilepEM"],false,(selectedEvent.trigMuEle&&selectedEvent.trigMuEleMatch));
            }
      }
   }
   ClearVariables();
   if(SelectEventTriggerStudies(TRIGDILEP_ptcuts)){
      if (selectedEvent.trigHt){//baselineTrigger
         if(selectedEvent.isDiElectron && !(selectedEvent.isDiMuon || selectedEvent.isMuonElectron || selectedEvent.isMuonElectron)){
            FillerTrigger(selectedEvent,eff1Maps["trigDilepEE_ptcuts"],false,(selectedEvent.trigDiEle&&selectedEvent.trigDiEleMatch));
         }
         if(selectedEvent.isDiMuon && !(selectedEvent.isDiElectron || selectedEvent.isMuonElectron || selectedEvent.isMuonElectron)){
            FillerTrigger(selectedEvent,eff1Maps["trigDilepMM_ptcuts"],false,(selectedEvent.trigDiMu&&selectedEvent.trigDiMuMatch));
         }
         if((selectedEvent.isMuonElectron||selectedEvent.isElectronMuon) && !(selectedEvent.isDiElectron || selectedEvent.isDiMuon)){
            FillerTrigger(selectedEvent,eff1Maps["trigDilepEM_ptcuts"],false,(selectedEvent.trigMuEle&&selectedEvent.trigMuEleMatch));
         }
      }
   }
}



void ntupleAnalyzerShrinked::FillHistograms(){
   ClearVariables();

   if(config_selectionsToProcessMap[DILEP]){
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
         if(selectedEvent.selLeptonSize==selectedEvent.matchedLeptonSize){

            if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)&&(selectedEvent.ETmiss<100.)){ //ONZ+MET<100
               if(selectedEvent.isDiMuon || selectedEvent.isDiElectron) Filler(selectedEvent,h1Maps["onZMet0100"],true);
               if (selectedEvent.isDiElectron){ 
                  Filler(selectedEvent,h1Maps["onZMet0100EE"],true);
               }else{
                  if (selectedEvent.isDiMuon) Filler(selectedEvent,h1Maps["onZMet0100MM"],true);
                  if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron) Filler(selectedEvent,h1Maps["onZMet0100EM"],true);
               }
            }
            
            if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)&&(selectedEvent.ETmiss>=100.)){
               if(selectedEvent.isDiMuon || selectedEvent.isDiElectron){
                  Filler(selectedEvent,h1Maps["onZMet100"],true);
                  }
               if (selectedEvent.isDiElectron){
                  Filler(selectedEvent,h1Maps["onZMet100EE"],true);
               }
               if (selectedEvent.isDiMuon){
                     Filler(selectedEvent,h1Maps["onZMet100MM"],true);
               }
               if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron){
                   Filler(selectedEvent,h1Maps["onZMet100EM"],true);
               }
            }
            
            if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)&&(selectedEvent.ETmiss>=150.)){ //ONZ+MET>150
               if(selectedEvent.isDiMuon || selectedEvent.isDiElectron){
                  Filler(selectedEvent,h1Maps["onZMet150"],true);
                  }
               if (selectedEvent.isDiElectron){
                  Filler(selectedEvent,h1Maps["onZMet150EE"],true);
               }
               if (selectedEvent.isDiMuon){
                     Filler(selectedEvent,h1Maps["onZMet150MM"],true);
               }
               if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron){
                   Filler(selectedEvent,h1Maps["onZMet150EM"],true);
               }
            }
            
            if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)&&(selectedEvent.ETmiss<150.)&&(selectedEvent.ETmiss>=100.)){ //ONZ+MET>100<150
               if(selectedEvent.isDiMuon || selectedEvent.isDiElectron){
                  Filler(selectedEvent,h1Maps["onZMet100150"],true);
                  }
               if (selectedEvent.isDiElectron){
                  Filler(selectedEvent,h1Maps["onZMet100150EE"],true);
               }
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
   
   //ClearVariables();
   //clearCutFlowMap();
   //if(config_selectionsToProcessMap[ControlRegionDY]){
      //if(SelectEvent(ONZ)){
         //if(selectedEvent.selLeptonSize==selectedEvent.matchedLeptonSize){
            //if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)&&(selectedEvent.ETmiss<100.)){ //ONZ+MET<100 ~ CR DY/Z(+gamma)
               //if(selectedEvent.isDiMuon || selectedEvent.isDiElectron) Filler(selectedEvent,h1Maps["CRDY"],true);
               //if (selectedEvent.isDiElectron){
                  //Filler(selectedEvent,h1Maps["CRDYEE"],true);
               //}else{                          
                  //if (selectedEvent.isDiMuon){Filler(selectedEvent,h1Maps["CRDYMM"],true);}
                  //if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron) Filler(selectedEvent,h1Maps["CRDYEM"],true);
               //}
            //}
         //}
      //}
   //}
   ClearVariables();
   clearCutFlowMap();
   if(config_selectionsToProcessMap[ControlRegionDY]){
      if(SelectEvent(ONZ)){
         if(selectedEvent.selLeptonSize==selectedEvent.matchedLeptonSize){
            if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)){//&&(selectedEvent.ETmiss<100.)){ //ONZ+MET<100 ~ CR DY/Z(+gamma)
               if(selectedEvent.ETmiss<100.){
                  if(selectedEvent.isDiMuon || selectedEvent.isDiElectron) Filler(selectedEvent,cr1Maps["CRDY"]["nom"],true,true);
                  if (selectedEvent.isDiElectron){
                     Filler(selectedEvent,cr1Maps["CRDYEE"]["nom"],true,true);
                  }else{                          
                     if (selectedEvent.isDiMuon){Filler(selectedEvent,cr1Maps["CRDYMM"]["nom"],true,true);}
                     if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron) Filler(selectedEvent,cr1Maps["CRDYEM"]["nom"],true,true);
                  }
                  
                  for(int i=0; i<nWeights; i++){
                     if(selectedEvent.isDiMuon || selectedEvent.isDiElectron) Filler(selectedEvent,cr1Maps["CRDY"][to_string(i)],true,true,i);
                     if (selectedEvent.isDiElectron){
                        Filler(selectedEvent,cr1Maps["CRDYEE"][to_string(i)],true,true,i);
                     }else{                          
                        if (selectedEvent.isDiMuon){Filler(selectedEvent,cr1Maps["CRDYMM"][to_string(i)],true,true,i);}
                        if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron) Filler(selectedEvent,cr1Maps["CRDYEM"][to_string(i)],true,true,i);
                     }
                  }
               }
               if((met_JESu->p.Pt()<100.)){ //ONZ+MET<100 ~ CR DY/Z(+gamma)
                  if(selectedEvent.isDiMuon || selectedEvent.isDiElectron) Filler(selectedEvent,cr1Maps["CRDY"]["JESu"],true,true,9999,"JESu");
                  if (selectedEvent.isDiElectron){
                     Filler(selectedEvent,cr1Maps["CRDYEE"]["JESu"],true,true,9999,"JESu");
                  }else{                          
                     if (selectedEvent.isDiMuon){Filler(selectedEvent,cr1Maps["CRDYMM"]["JESu"],true,true,9999,"JESu");}
                     if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron) Filler(selectedEvent,cr1Maps["CRDYEM"]["JESu"],true,true,9999,"JESu");
                  }
               }
               if((met_JESd->p.Pt()<100.)){ //ONZ+MET<100 ~ CR DY/Z(+gamma)
                  if(selectedEvent.isDiMuon || selectedEvent.isDiElectron) Filler(selectedEvent,cr1Maps["CRDY"]["JESd"],true,true,9999,"JESd");
                  if (selectedEvent.isDiElectron){
                     Filler(selectedEvent,cr1Maps["CRDYEE"]["JESd"],true,true,9999,"JESd");
                  }else{                          
                     if (selectedEvent.isDiMuon){Filler(selectedEvent,cr1Maps["CRDYMM"]["JESd"],true,true,9999,"JESd");}
                     if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron) Filler(selectedEvent,cr1Maps["CRDYEM"]["JESd"],true,true,9999,"JESd");
                  }
               }
               if((met_JERu->p.Pt()<100.)){ //ONZ+MET<100 ~ CR DY/Z(+gamma)
                  if(selectedEvent.isDiMuon || selectedEvent.isDiElectron) Filler(selectedEvent,cr1Maps["CRDY"]["JERu"],true,true,9999,"JERu");
                  if (selectedEvent.isDiElectron){
                     Filler(selectedEvent,cr1Maps["CRDYEE"]["JERu"],true,true,9999,"JERu");
                  }else{                          
                     if (selectedEvent.isDiMuon){Filler(selectedEvent,cr1Maps["CRDYMM"]["JERu"],true,true,9999,"JERu");}
                     if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron) Filler(selectedEvent,cr1Maps["CRDYEM"]["JERu"],true,true,9999,"JERu");
                  }
               }
               if((met_JESd->p.Pt()<100.)){ //ONZ+MET<100 ~ CR DY/Z(+gamma)
                  if(selectedEvent.isDiMuon || selectedEvent.isDiElectron) Filler(selectedEvent,cr1Maps["CRDY"]["JERd"],true,true,9999,"JERd");
                  if (selectedEvent.isDiElectron){
                     Filler(selectedEvent,cr1Maps["CRDYEE"]["JERd"],true,true,9999,"JERd");
                  }else{                          
                     if (selectedEvent.isDiMuon){Filler(selectedEvent,cr1Maps["CRDYMM"]["JERd"],true,true,9999,"JERd");}
                     if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron) Filler(selectedEvent,cr1Maps["CRDYEM"]["JERd"],true,true,9999,"JERd");
                  }
               }
            }
         }
      }
   }
   if(config_selectionsToProcessMap[ControlRegionWW]){
      if(SelectEvent(DILEP)){
         if(selectedEvent.selLeptonSize==selectedEvent.matchedLeptonSize){
            if ((selectedEvent.mll<81. || selectedEvent.mll>101.)&&(selectedEvent.ETmiss<100. && selectedEvent.ETmiss>20.)){ //ONZ+MET<100 ~ CR DY/Z(+gamma)
               if((selectedEvent.selJets.size()==0)&&((selectedEvent.l1+selectedEvent.l2).Pt()>45.)){
               
                  if(selectedEvent.isDiMuon || selectedEvent.isDiElectron || selectedEvent.isElectronMuon || selectedEvent.isMuonElectron){
                     //Filler(selectedEvent,h1Maps["CRWWLL"],false,false);
                     Filler(selectedEvent,h1Maps["CRWWLL"],false);
                   }
                  if(selectedEvent.isDiMuon || selectedEvent.isDiElectron){
                     //Filler(selectedEvent,h1Maps["CRWW"],false,false);
                     Filler(selectedEvent,h1Maps["CRWW"],false);
                   }
                  if (selectedEvent.isDiElectron){
                     //Filler(selectedEvent,h1Maps["CRWWEE"],false,false);
                     Filler(selectedEvent,h1Maps["CRWWEE"],false);
                  }                      
                  if (selectedEvent.isDiMuon){
                     //Filler(selectedEvent,h1Maps["CRWWMM"],false,false);
                     Filler(selectedEvent,h1Maps["CRWWMM"],false);
                  }
                  if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron){
                     //Filler(selectedEvent,h1Maps["CRWWEM"],false,false);
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
   //if(config_selectionsToProcessMap[ControlRegionTT]){
      //if(SelectEvent(SEL)){
         //if(selectedEvent.selLeptonSize==selectedEvent.matchedLeptonSize){
            //if ((selectedEvent.selPhotons.size()!=0)){ //Different Flavor + 1Photon ~ CR TT(+gamma)
               //if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron){
                  //Filler(selectedEvent,h1Maps["CRTTEM"],true);
                  //if(selectedEvent.selPhotons.at(0).p.Pt()>=80.){
                     //Filler(selectedEvent,h1Maps["CRTT80EM"],true);
                  //}else{
                     //Filler(selectedEvent,h1Maps["CRTT080EM"],true);
                  //}
               //}
            //}
         //}
      //}
   //}
   if(config_selectionsToProcessMap[ControlRegionTT]){
      if(SelectEvent(SEL)){
         if(selectedEvent.selLeptonSize==selectedEvent.matchedLeptonSize){
            if ((selectedEvent.selPhotons.size()!=0)){ //Different Flavor + 1Photon ~ CR TT(+gamma)
               if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron){
                  Filler(selectedEvent,cr1Maps["CRTTEM"]["nom"],true,true);
                  Filler(selectedEvent,cr1Maps["CRTTEM"]["JESu"],true,true,9999,"JESu");
                  Filler(selectedEvent,cr1Maps["CRTTEM"]["JESd"],true,true,9999,"JESd");
                  Filler(selectedEvent,cr1Maps["CRTTEM"]["JERu"],true,true,9999,"JERu");
                  Filler(selectedEvent,cr1Maps["CRTTEM"]["JERd"],true,true,9999,"JERd");
                  for(int i=0; i<nWeights; i++){
                     if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron) Filler(selectedEvent,cr1Maps["CRTTEM"][to_string(i)],true,true,i);
                  }
                  if(selectedEvent.selPhotons.at(0).p.Pt()>=80.){
                     Filler(selectedEvent,cr1Maps["CRTT80EM"]["nom"],true,true);
                     //Filler(selectedEvent,cr1Maps["CRTT80EM"]["JESup"],true,9999,"JESup");
                     //Filler(selectedEvent,cr1Maps["CRTT80EM"]["JESup"],true,9999,"JESdn");
                     //Filler(selectedEvent,cr1Maps["CRTT80EM"]["JERup"],true,9999,"JERup");
                     //Filler(selectedEvent,cr1Maps["CRTT80EM"]["JERup"],true,9999,"JERdn");
                  //for(int i=0; i<nWeights; i++){
                     //if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron) Filler(selectedEvent,cr1Maps["CRTTEM"][to_string(i)],true,i);
                  //}
                  }else{
                     Filler(selectedEvent,cr1Maps["CRTT080EM"]["nom"],true,true);
                     //Filler(selectedEvent,cr1Maps["CRTT080EM"]["JESup"],true,9999,"JESup");
                     //Filler(selectedEvent,cr1Maps["CRTT080EM"]["JESup"],true,9999,"JESdn");
                     //Filler(selectedEvent,cr1Maps["CRTT080EM"]["JERup"],true,9999,"JERup");
                     //Filler(selectedEvent,cr1Maps["CRTT080EM"]["JERup"],true,9999,"JERdn");
                  //}
               
                     //for(int i=0; i<nWeights; i++){
                        //if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron) Filler(selectedEvent,cr1Maps["CRTTEM"][to_string(i)],true,i);
                     //}
                  }
               }
            }
         }
      }
   }
   ClearVariables();
   clearCutFlowMap();
   if(config_selectionsToProcessMap[ControlRegionZZ]){
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
                        //FillerZZ(selectedEvent,h1Maps["CRZZMM"],false,lep1,lep2,lep3,lep4);
                        FillerZZ(selectedEvent,cr1Maps["CRZZMM"]["nom"],false,lep1,lep2,lep3,lep4,true);
                        FillerZZ(selectedEvent,cr1Maps["CRZZMM"]["JESu"],false,lep1,lep2,lep3,lep4,true,9999,"JESu");
                        FillerZZ(selectedEvent,cr1Maps["CRZZMM"]["JESd"],false,lep1,lep2,lep3,lep4,true,9999,"JESd");
                        FillerZZ(selectedEvent,cr1Maps["CRZZMM"]["JERu"],false,lep1,lep2,lep3,lep4,true,9999,"JERu");
                        FillerZZ(selectedEvent,cr1Maps["CRZZMM"]["JERd"],false,lep1,lep2,lep3,lep4,true,9999,"JERd");
                        for(int i=0; i<nWeights; i++){
                           FillerZZ(selectedEvent,cr1Maps["CRZZMM"][to_string(i)],false,lep1,lep2,lep3,lep4,true,i);
                        }
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
                        //FillerZZ(selectedEvent,h1Maps["CRZZEE"],false,lep1,lep2,lep3,lep4);
                        FillerZZ(selectedEvent,cr1Maps["CRZZEE"]["nom"],false,lep1,lep2,lep3,lep4,true);
                        FillerZZ(selectedEvent,cr1Maps["CRZZEE"]["JESu"],false,lep1,lep2,lep3,lep4,true,9999,"JESu");
                        FillerZZ(selectedEvent,cr1Maps["CRZZEE"]["JESd"],false,lep1,lep2,lep3,lep4,true,9999,"JESd");
                        FillerZZ(selectedEvent,cr1Maps["CRZZEE"]["JERu"],false,lep1,lep2,lep3,lep4,true,9999,"JERu");
                        FillerZZ(selectedEvent,cr1Maps["CRZZEE"]["JERd"],false,lep1,lep2,lep3,lep4,true,9999,"JERd");
                        for(int i=0; i<nWeights; i++){
                           FillerZZ(selectedEvent,cr1Maps["CRZZEE"][to_string(i)],false,lep1,lep2,lep3,lep4,true,i);
                        }
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
                              //FillerZZ(selectedEvent,h1Maps["CRZZMM"],false,lep1,lep2,lep3,lep4);
                              FillerZZ(selectedEvent,cr1Maps["CRZZMM"]["nom"],false,lep1,lep2,lep3,lep4,true);
                              FillerZZ(selectedEvent,cr1Maps["CRZZMM"]["JESu"],false,lep1,lep2,lep3,lep4,true,9999,"JESu");
                              FillerZZ(selectedEvent,cr1Maps["CRZZMM"]["JESd"],false,lep1,lep2,lep3,lep4,true,9999,"JESd");
                              FillerZZ(selectedEvent,cr1Maps["CRZZMM"]["JERu"],false,lep1,lep2,lep3,lep4,true,9999,"JERu");
                              FillerZZ(selectedEvent,cr1Maps["CRZZMM"]["JERd"],false,lep1,lep2,lep3,lep4,true,9999,"JERd");
                              for(int i=0; i<nWeights; i++){
                                 FillerZZ(selectedEvent,cr1Maps["CRZZMM"][to_string(i)],false,lep1,lep2,lep3,lep4,true,i);
                              }
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
                              //FillerZZ(selectedEvent,h1Maps["CRZZEE"],false,lep1,lep2,lep3,lep4);
                              FillerZZ(selectedEvent,cr1Maps["CRZZEE"]["nom"],false,lep1,lep2,lep3,lep4,true);
                              FillerZZ(selectedEvent,cr1Maps["CRZZEE"]["JESu"],false,lep1,lep2,lep3,lep4,true,9999,"JESu");
                              FillerZZ(selectedEvent,cr1Maps["CRZZEE"]["JESd"],false,lep1,lep2,lep3,lep4,true,9999,"JESd");
                              FillerZZ(selectedEvent,cr1Maps["CRZZEE"]["JERu"],false,lep1,lep2,lep3,lep4,true,9999,"JERu");
                              FillerZZ(selectedEvent,cr1Maps["CRZZEE"]["JERd"],false,lep1,lep2,lep3,lep4,true,9999,"JERd");
                              for(int i=0; i<nWeights; i++){
                                 FillerZZ(selectedEvent,cr1Maps["CRZZEE"][to_string(i)],false,lep1,lep2,lep3,lep4,true,i);
                              }
                           }
                        }
                     }
                  }
                  float mll1=(lep1.vec+lep2.vec).M();
                  float mll2=(lep3.vec+lep4.vec).M();
                  if(mll1<106. && mll1>76. && mll2<130. && mll2>50.){
                     //FillerZZ(selectedEvent,h1Maps["CRZZ"],false,lep1,lep2,lep3,lep4);
                     FillerZZ(selectedEvent,cr1Maps["CRZZ"]["nom"],false,lep1,lep2,lep3,lep4,true);
                     FillerZZ(selectedEvent,cr1Maps["CRZZ"]["JESu"],false,lep1,lep2,lep3,lep4,true,9999,"JESu");
                     FillerZZ(selectedEvent,cr1Maps["CRZZ"]["JESd"],false,lep1,lep2,lep3,lep4,true,9999,"JESd");
                     FillerZZ(selectedEvent,cr1Maps["CRZZ"]["JERu"],false,lep1,lep2,lep3,lep4,true,9999,"JERu");
                     FillerZZ(selectedEvent,cr1Maps["CRZZ"]["JERd"],false,lep1,lep2,lep3,lep4,true,9999,"JERd");
                     for(int i=0; i<nWeights; i++){
                        FillerZZ(selectedEvent,cr1Maps["CRZZ"][to_string(i)],false,lep1,lep2,lep3,lep4,true,i);
                     }
                  }
               }
            }
         }
      }
   }
   
ClearVariables();
   clearCutFlowMap();
   if(config_selectionsToProcessMap[ControlRegionWZ]){
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
                              //FillerWZ(selectedEvent,h1Maps["CRWZMM"],false,lep1,lep2,lep3);
                              FillerWZ(selectedEvent,cr1Maps["CRWZMM"]["nom"],false,lep1,lep2,lep3,true);
                              for(int i=0; i<nWeights; i++){
                                 FillerWZ(selectedEvent,cr1Maps["CRWZMM"][to_string(i)],false,lep1,lep2,lep3,true,i);
                              }
                           }
                           if(met_JESu->p.Pt()>70.){
                              FillerWZ(selectedEvent,cr1Maps["CRWZMM"]["JESu"],false,lep1,lep2,lep3,true,9999,"JESu");
                           }
                           if(met_JESd->p.Pt()>70.){
                              FillerWZ(selectedEvent,cr1Maps["CRWZMM"]["JESd"],false,lep1,lep2,lep3,true,9999,"JESd");
                           }
                           if(met_JERu->p.Pt()>70.){
                              FillerWZ(selectedEvent,cr1Maps["CRWZMM"]["JERu"],false,lep1,lep2,lep3,true,9999,"JERu");
                           }
                           if(met_JERd->p.Pt()>70.){
                              FillerWZ(selectedEvent,cr1Maps["CRWZMM"]["JERd"],false,lep1,lep2,lep3,true,9999,"JERd");
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
                           //if(selectedEvent.ETmiss>70.){
                              //FillerWZ(selectedEvent,h1Maps["CRWZEE"],false,lep1,lep2,lep3);
                           //}
                           if(selectedEvent.ETmiss>70.){
                              FillerWZ(selectedEvent,cr1Maps["CRWZEE"]["nom"],false,lep1,lep2,lep3,true);
                              for(int i=0; i<nWeights; i++){
                                 FillerWZ(selectedEvent,cr1Maps["CRWZEE"][to_string(i)],false,lep1,lep2,lep3,true,i);
                              }
                           }
                           if(met_JESu->p.Pt()>70.){
                              FillerWZ(selectedEvent,cr1Maps["CRWZEE"]["JESu"],false,lep1,lep2,lep3,true,9999,"JESu");
                           }
                           if(met_JESd->p.Pt()>70.){
                              FillerWZ(selectedEvent,cr1Maps["CRWZEE"]["JESd"],false,lep1,lep2,lep3,true,9999,"JESd");
                           }
                           if(met_JERu->p.Pt()>70.){
                              FillerWZ(selectedEvent,cr1Maps["CRWZEE"]["JERu"],false,lep1,lep2,lep3,true,9999,"JERu");
                           }
                           if(met_JERd->p.Pt()>70.){
                              FillerWZ(selectedEvent,cr1Maps["CRWZEE"]["JERd"],false,lep1,lep2,lep3,true,9999,"JERd");
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
                           //if(selectedEvent.ETmiss>70.){
                              //FillerWZ(selectedEvent,h1Maps["CRWZEE"],false,lep1,lep2,lep3);
                           //}
                           if(selectedEvent.ETmiss>70.){
                              FillerWZ(selectedEvent,cr1Maps["CRWZEE"]["nom"],false,lep1,lep2,lep3,true);
                              for(int i=0; i<nWeights; i++){
                                 FillerWZ(selectedEvent,cr1Maps["CRWZEE"][to_string(i)],false,lep1,lep2,lep3,true,i);
                              }
                           }
                           if(met_JESu->p.Pt()>70.){
                              FillerWZ(selectedEvent,cr1Maps["CRWZEE"]["JESu"],false,lep1,lep2,lep3,true,9999,"JESu");
                           }
                           if(met_JESd->p.Pt()>70.){
                              FillerWZ(selectedEvent,cr1Maps["CRWZEE"]["JESd"],false,lep1,lep2,lep3,true,9999,"JESd");
                           }
                           if(met_JERu->p.Pt()>70.){
                              FillerWZ(selectedEvent,cr1Maps["CRWZEE"]["JERu"],false,lep1,lep2,lep3,true,9999,"JERu");
                           }
                           if(met_JERd->p.Pt()>70.){
                              FillerWZ(selectedEvent,cr1Maps["CRWZEE"]["JERd"],false,lep1,lep2,lep3,true,9999,"JERd");
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
                           //if(selectedEvent.ETmiss>70.){
                              //FillerWZ(selectedEvent,h1Maps["CRWZMM"],false,lep1,lep2,lep3);
                           //}
                           if(selectedEvent.ETmiss>70.){
                              //FillerWZ(selectedEvent,h1Maps["CRWZMM"],false,lep1,lep2,lep3);
                              FillerWZ(selectedEvent,cr1Maps["CRWZMM"]["nom"],false,lep1,lep2,lep3,true);
                              for(int i=0; i<nWeights; i++){
                                 FillerWZ(selectedEvent,cr1Maps["CRWZMM"][to_string(i)],false,lep1,lep2,lep3,true,i);
                              }
                           }
                           if(met_JESu->p.Pt()>70.){
                              FillerWZ(selectedEvent,cr1Maps["CRWZMM"]["JESu"],false,lep1,lep2,lep3,true,9999,"JESu");
                           }
                           if(met_JESd->p.Pt()>70.){
                              FillerWZ(selectedEvent,cr1Maps["CRWZMM"]["JESd"],false,lep1,lep2,lep3,true,9999,"JESd");
                           }
                           if(met_JERu->p.Pt()>70.){
                              FillerWZ(selectedEvent,cr1Maps["CRWZMM"]["JERu"],false,lep1,lep2,lep3,true,9999,"JERu");
                           }
                           if(met_JERd->p.Pt()>70.){
                              FillerWZ(selectedEvent,cr1Maps["CRWZMM"]["JERd"],false,lep1,lep2,lep3,true,9999,"JERd");
                           }                        
                        }
                     }                  
                  }
                  
                               
                  
               float mll1=(lep1.vec+lep2.vec).M();
               if(mll1<106. && mll1>76.){
                  if( (lep1.vec+selectedEvent.ETmiss_vec).Mt() > 50. ){
                     //if(selectedEvent.ETmiss>70.){
                        //FillerWZ(selectedEvent,h1Maps["CRWZ"],false,lep1,lep2,lep3);
                     //}
                     if(selectedEvent.ETmiss>70.){
                        FillerWZ(selectedEvent,cr1Maps["CRWZ"]["nom"],false,lep1,lep2,lep3,true);
                        for(int i=0; i<nWeights; i++){
                           FillerWZ(selectedEvent,cr1Maps["CRWZ"][to_string(i)],false,lep1,lep2,lep3,true,i);
                        }
                     }
                     if(met_JESu->p.Pt()>70.){
                        FillerWZ(selectedEvent,cr1Maps["CRWZ"]["JESu"],false,lep1,lep2,lep3,true,9999,"JESu");
                     }
                     if(met_JESd->p.Pt()>70.){
                        FillerWZ(selectedEvent,cr1Maps["CRWZ"]["JESd"],false,lep1,lep2,lep3,true,9999,"JESd");
                     }
                     if(met_JERu->p.Pt()>70.){
                        FillerWZ(selectedEvent,cr1Maps["CRWZ"]["JERu"],false,lep1,lep2,lep3,true,9999,"JERu");
                     }
                     if(met_JERd->p.Pt()>70.){
                        FillerWZ(selectedEvent,cr1Maps["CRWZ"]["JERd"],false,lep1,lep2,lep3,true,9999,"JERd");
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
      if(SelectEvent(ONZ)){
            if(selectedEvent.selLeptonSize==selectedEvent.matchedLeptonSize){
            if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)){
               if(selectedEvent.isDiMuon || selectedEvent.isDiElectron){
                  Filler(selectedEvent,h1Maps["onZ"],true);
               }
               if (selectedEvent.isDiElectron){
                  Filler(selectedEvent,h1Maps["onZEE"],true);
               }else{
                  if (selectedEvent.isDiMuon){
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


void ntupleAnalyzerShrinked::Filler(selEvent& ev, map<Histograms1D,TH1F>& m,bool withPhoton,bool slimmed,int changePDF,string changeMET){
   
   float weight_TopPt=1.;
   float weight_LeptonPt=1.;
   float weight_nIsr=1.;
   float weight_EWKinoPt=1.;
   
   float weight_PDF=1.;
   
   if(changePDF>150){
      weight_PDF=1.;
   }else{
      weight_PDF=pdf_weights->at(changePDF);
   }
   
   
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
   
   tempWeight = tempWeight*weight_TopPt*weight_LeptonPt*weight_nIsr*weight_EWKinoPt*weight_PDF;
   
   m.at(WEIGHT_EWKINOPAIRPT).Fill(weight_EWKinoPt);
   m.at(WEIGHT_LEPTONPAIRPT).Fill(weight_LeptonPt);
   m.at(WEIGHT_NISR).Fill(weight_nIsr);
   m.at(WEIGHT_TOPPT).Fill(weight_TopPt);
   m.at(WEIGHT_PDF).Fill(weight_PDF);
   
   if(changeMET=="N"){
      m.at(ETMISS).Fill(ev.ETmiss, tempWeight);
   }
   else{
      if(changeMET=="JESu"){
         m.at(ETMISS).Fill(met_JESu->p.Pt(),tempWeight);
      }else{
         if(changeMET=="JESd"){
            m.at(ETMISS).Fill(met_JESd->p.Pt(),tempWeight);
         }else{
            if(changeMET=="JERu"){
               m.at(ETMISS).Fill(met_JERu->p.Pt(),tempWeight);
            }else{
               if(changeMET=="JERd"){
                  m.at(ETMISS).Fill(met_JERd->p.Pt(),tempWeight);
               }
            }
         }
      }
   }
   
   
   m.at(PT1).Fill(ev.pt1, tempWeight);
   m.at(PT2).Fill(ev.pt2, tempWeight);
   m.at(MLL).Fill(ev.mll,tempWeight);   
   m.at(NPHOTONS).Fill(ev.selPhotons.size(),tempWeight);
   m.at(HT).Fill(ev.calcHt,tempWeight);
   m.at(NJETS).Fill(ev.nselJets,tempWeight);
   m.at(ETA1).Fill(fabs(ev.eta1),tempWeight);
   m.at(ETA2).Fill(fabs(ev.eta2),tempWeight);
   m.at(PHI1).Fill(fabs(ev.phi2),tempWeight);
   m.at(PHI2).Fill(fabs(ev.phi2),tempWeight);

   //m.at(MT2).Fill(ev.MT2_val,tempWeight);

   
   int jetSize=ev.selJets.size();



   
   
   if(!slimmed){
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
      m.at(NVTX).Fill(*nGoodVertices,tempWeight);
      m.at(GENHT).Fill(*genHt,tempWeight);
      m.at(NBJETS).Fill(ev.nselBJets,tempWeight);
      m.at(DeltaEtaLL).Fill(fabs(ev.eta1-ev.eta2),tempWeight);
      m.at(DeltaPhiLL).Fill(fabs(ev.l1.DeltaPhi(ev.l2)),tempWeight);
      m.at(DeltaRLL).Fill(fabs(ev.deltaRll),tempWeight);
      m.at(DeltaEtaLL_neg).Fill(ev.eta1-ev.eta2,tempWeight);
      m.at(DeltaPhiLL_neg).Fill(ev.l1.DeltaPhi(ev.l2),tempWeight);
      m.at(DeltaRLL_neg).Fill(ev.deltaRll,tempWeight);
      m.at(ZPT).Fill((ev.l1+ev.l2).Pt(),tempWeight);
      m.at(MTLL).Fill((ev.l1+ev.l2).Mt(),tempWeight);
      m.at(ST).Fill(ev.pt1+ev.pt2,tempWeight);
      m.at(DeltaPhiLLMet).Fill(fabs((ev.l1+ev.l2).DeltaPhi(ev.ETmiss_vec)),tempWeight);
      m.at(DeltaEtaLLMet).Fill(fabs((ev.l1+ev.l2).Eta() - ev.ETmiss_vec.Eta()),tempWeight);
      m.at(DeltaRLLMet).Fill(fabs((ev.l1+ev.l2).DeltaR(ev.ETmiss_vec)),tempWeight);   
   }
   
   
   //if(ev.evtHasGenPhotonVeto) m.at(VetoCompare).Fill(1.,tempWeight);
   //if(!ev.evtHasGenPhotonVeto) m.at(VetoCompare).Fill(0.,tempWeight);
   
   if(withPhoton){
      m.at(PTG1).Fill(ev.selPhotons.at(0).p.Pt(),tempWeight);
      m.at(PHIG1).Fill(ev.selPhotons.at(0).p.Phi(),tempWeight);
      m.at(ETAG1).Fill(ev.selPhotons.at(0).p.Eta(),tempWeight);
      
      if(ev.selPhotons.at(0).matchedToPhoton){
         m.at(Fakes).Fill(0.,tempWeight);
      }else{
         if(ev.selPhotons.at(0).matchedToElectron){
            m.at(Fakes).Fill(1.,tempWeight);
         }else{
            if(ev.selPhotons.at(0).matchedToJet){
               m.at(Fakes).Fill(2.,tempWeight);
            }else{
               m.at(Fakes).Fill(3.,tempWeight);
            }
         }
      }
      
      if(!slimmed){
         m.at(SIGMAIETAIETAG1).Fill(ev.selPhotons.at(0).sigmaIetaIeta,tempWeight);
         m.at(SIGMAIPHIIPHIG1).Fill(ev.selPhotons.at(0).sigmaIphiIphi,tempWeight);
         m.at(R9).Fill(ev.selPhotons.at(0).r9,tempWeight);
         m.at(HOVERE).Fill(ev.selPhotons.at(0).hOverE,tempWeight);
         m.at(DELTARGL1).Fill(ev.selPhotons.at(0).deltaR1,tempWeight);
         m.at(DELTARGL2).Fill(ev.selPhotons.at(0).deltaR2,tempWeight);
         m.at(DeltaEtaLLG).Fill(fabs((ev.l1+ev.l2).Eta()-ev.selPhotons.at(0).vec.Eta()), tempWeight);
         m.at(DeltaPhiLLG).Fill(fabs((ev.l1+ev.l2).DeltaPhi(ev.selPhotons.at(0).vec)),tempWeight);
         m.at(DeltaRLLG).Fill(fabs((ev.l1+ev.l2).DeltaR(ev.selPhotons.at(0).vec)),tempWeight);
         m.at(DeltaPhiGMet).Fill(fabs((ev.ETmiss_vec).DeltaPhi(ev.selPhotons.at(0).vec)),tempWeight);
         m.at(DeltaRGMet).Fill(fabs((ev.ETmiss_vec).DeltaR(ev.selPhotons.at(0).vec)),tempWeight);
         m.at(MTLLG).Fill((ev.l1+ev.l2+ev.selPhotons.at(0).vec).Mt(),tempWeight);
         m.at(MTL1MET).Fill((ev.l1+ev.ETmiss_vec).Mt(),tempWeight);
         m.at(MTL2MET).Fill((ev.l2+ev.ETmiss_vec).Mt(),tempWeight);
         m.at(MTGMET).Fill((ev.selPhotons.at(0).vec+ev.ETmiss_vec).Mt(),tempWeight);
         m.at(MTLLMET).Fill((ev.l1+ev.l2+ev.ETmiss_vec).Mt(),tempWeight);
         m.at(MTLLGMET).Fill((ev.l1+ev.l2+ev.selPhotons.at(0).vec+ev.ETmiss_vec).Mt(),tempWeight);
         m.at(STG).Fill(ev.pt1+ev.pt2+ev.selPhotons.at(0).p.Pt(),tempWeight);
         m.at(STMET).Fill(ev.pt1+ev.pt2+ev.selPhotons.at(0).p.Pt()+ev.ETmiss,tempWeight);
         m.at(MLLG).Fill((ev.l1+ev.l2+ev.selPhotons.at(0).vec).M(),tempWeight);      
      }
      
      
      
      //m.at(PT_llg).Fill((ev.l1+ev.l2+ev.selPhotons.at(0).vec).Pt(),tempWeight);
      
      //if(FindGenPhotonMatch(ev.selPhotons.at(0))) m.at(genPhotonPT).Fill(GetGenPhotonMatch(ev.selPhotons.at(0)).p.Pt(),tempWeight);
      
      //if(ev.evtHasGenPhotonVeto) if(FindGenPhotonMatch(ev.selPhotons.at(0))) m.at(genPhotonPT_Veto).Fill(GetGenPhotonMatch(ev.selPhotons.at(0)).p.Pt(),tempWeight);
      //if(ev.evtHasGenPhotonVeto) m.at(PTG1_Veto).Fill(ev.selPhotons.at(0).p.Pt(),tempWeight);
      //if(!ev.evtHasGenPhotonVeto) if(FindGenPhotonMatch(ev.selPhotons.at(0))) m.at(genPhotonPT_NoVeto).Fill(GetGenPhotonMatch(ev.selPhotons.at(0)).p.Pt(),tempWeight);
      //if(!ev.evtHasGenPhotonVeto) m.at(PTG1_NoVeto).Fill(ev.selPhotons.at(0).p.Pt(),tempWeight);
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
void ntupleAnalyzerShrinked::FillerZZ(selEvent& ev, map<Histograms1D,TH1F>& m,bool withPhoton,selLepton& l1, selLepton& l2, selLepton& l3, selLepton& l4,bool slimmed,int changePDF,string changeMET){
   
   float weight_TopPt=1.;
   float weight_LeptonPt=1.;
   float weight_nIsr=1.;
   float weight_EWKinoPt=1.;
   
   
   float weight_PDF=1.;
   
   if(changePDF>150){
      weight_PDF=1.;
   }else{
      //weight_PDF=isData? 1. : pdf_weights->at(changePDF);
      weight_PDF=pdf_weights->at(changePDF);
      //weight_PDF=1.;
   }
   
   
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
   
   tempWeight = tempWeight*weight_TopPt*weight_LeptonPt*weight_nIsr*weight_EWKinoPt*weight_PDF;
   
   m.at(WEIGHT_EWKINOPAIRPT).Fill(weight_EWKinoPt);
   m.at(WEIGHT_LEPTONPAIRPT).Fill(weight_LeptonPt);
   m.at(WEIGHT_NISR).Fill(weight_nIsr);
   m.at(WEIGHT_TOPPT).Fill(weight_TopPt);
   m.at(WEIGHT_PDF).Fill(weight_PDF);

   
   if(changeMET=="N"){
      m.at(ETMISS).Fill(ev.ETmiss, tempWeight);
   }
   else{
      if(changeMET=="JESu"){
         m.at(ETMISS).Fill(met_JESu->p.Pt(),tempWeight);
      }else{
         if(changeMET=="JESd"){
            m.at(ETMISS).Fill(met_JESd->p.Pt(),tempWeight);
         }else{
            if(changeMET=="JERu"){
               m.at(ETMISS).Fill(met_JERu->p.Pt(),tempWeight);
            }else{
               if(changeMET=="JERd"){
                  m.at(ETMISS).Fill(met_JERd->p.Pt(),tempWeight);
               }
            }
         }
      }
   }
   
   
   //m.at(ETMISS).Fill(ev.ETmiss, tempWeight);
   m.at(PT1).Fill(l1.p.Pt(), tempWeight);
   m.at(PT2).Fill(l2.p.Pt(), tempWeight);
   m.at(PT3).Fill(l3.p.Pt(), tempWeight);
   m.at(PT4).Fill(l4.p.Pt(), tempWeight);
   m.at(MLL).Fill((l1.vec+l2.vec).M(),tempWeight);   
   m.at(MLL2).Fill((l3.vec+l4.vec).M(),tempWeight);   
   m.at(NPHOTONS).Fill(ev.selPhotons.size(),tempWeight);
   m.at(HT).Fill(ev.calcHt,tempWeight);
   m.at(NJETS).Fill(ev.selJets.size(),tempWeight);
   m.at(ETA1).Fill(fabs(l1.p.Eta()),tempWeight);
   m.at(ETA2).Fill(fabs(l2.p.Eta()),tempWeight);
   m.at(ETA3).Fill(fabs(l3.p.Eta()),tempWeight);
   m.at(ETA4).Fill(fabs(l4.p.Eta()),tempWeight);
   m.at(PHI1).Fill(fabs(l1.p.Phi()),tempWeight);
   m.at(PHI2).Fill(fabs(l2.p.Phi()),tempWeight);
   m.at(PHI3).Fill(fabs(l3.p.Phi()),tempWeight);
   m.at(PHI4).Fill(fabs(l4.p.Phi()),tempWeight);




   int jetSize=ev.selJets.size();









   if(!slimmed){
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
      m.at(NVTX).Fill(*nGoodVertices,tempWeight);
      m.at(MLLLL).Fill((l1.vec+l2.vec+l3.vec+l4.vec).M(),tempWeight);   
      m.at(GENHT).Fill(*genHt,tempWeight);
      m.at(ZPT).Fill((l1.vec+l2.vec).Pt(),tempWeight);
      m.at(ZPT2).Fill((l3.vec+l4.vec).Pt(),tempWeight);
      m.at(ST).Fill(ev.pt1+ev.pt2,tempWeight);
      m.at(MT2).Fill(ev.MT2_val,tempWeight);
      m.at(NElectrons).Fill(ev.selElectrons.size(),tempWeight);
      m.at(NMuons).Fill(ev.selMuons.size(),tempWeight);
   }

   
   if(withPhoton){
      m.at(PTG1).Fill(ev.selPhotons.at(0).p.Pt(),tempWeight);
      m.at(PHIG1).Fill(ev.selPhotons.at(0).p.Phi(),tempWeight);
      m.at(ETAG1).Fill(ev.selPhotons.at(0).p.Eta(),tempWeight);
      
      
      //if(FindGenPhotonMatch(ev.selPhotons.at(0))) m.at(genPhotonPT).Fill(GetGenPhotonMatch(ev.selPhotons.at(0)).p.Pt(),tempWeight);
      //
      //if(ev.evtHasGenPhotonVeto) if(FindGenPhotonMatch(ev.selPhotons.at(0))) m.at(genPhotonPT_Veto).Fill(GetGenPhotonMatch(ev.selPhotons.at(0)).p.Pt(),tempWeight);
      //if(ev.evtHasGenPhotonVeto) m.at(PTG1_Veto).Fill(ev.selPhotons.at(0).p.Pt(),tempWeight);
      //if(!ev.evtHasGenPhotonVeto) if(FindGenPhotonMatch(ev.selPhotons.at(0))) m.at(genPhotonPT_NoVeto).Fill(GetGenPhotonMatch(ev.selPhotons.at(0)).p.Pt(),tempWeight);
      //if(!ev.evtHasGenPhotonVeto) m.at(PTG1_NoVeto).Fill(ev.selPhotons.at(0).p.Pt(),tempWeight);
      
      
      if(!slimmed){
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
      }
      
   }
}

void ntupleAnalyzerShrinked::FillerWZ(selEvent& ev, map<Histograms1D,TH1F>& m,bool withPhoton,selLepton& l1, selLepton& l2, selLepton& l3,bool slimmed,int changePDF,string changeMET){
   
   float weight_TopPt=1.;
   float weight_LeptonPt=1.;
   float weight_nIsr=1.;
   float weight_EWKinoPt=1.;
   
   
   float weight_PDF=1.;
   
   if(changePDF>150){
      weight_PDF=1.;
   }else{
      weight_PDF=pdf_weights->at(changePDF);
      //weight_PDF=isData? 1. : pdf_weights->at(changePDF);
      //weight_PDF=1.;
   }
   
   
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
   
   tempWeight = tempWeight*weight_TopPt*weight_LeptonPt*weight_nIsr*weight_EWKinoPt*weight_PDF;
   
   m.at(WEIGHT_EWKINOPAIRPT).Fill(weight_EWKinoPt);
   m.at(WEIGHT_LEPTONPAIRPT).Fill(weight_LeptonPt);
   m.at(WEIGHT_NISR).Fill(weight_nIsr);
   m.at(WEIGHT_TOPPT).Fill(weight_TopPt);
   m.at(WEIGHT_PDF).Fill(weight_PDF);
   
   if(changeMET=="N"){
      m.at(ETMISS).Fill(ev.ETmiss, tempWeight);
   }
   else{
      if(changeMET=="JESu"){
         m.at(ETMISS).Fill(met_JESu->p.Pt(),tempWeight);
      }else{
         if(changeMET=="JESd"){
            m.at(ETMISS).Fill(met_JESd->p.Pt(),tempWeight);
         }else{
            if(changeMET=="JERu"){
               m.at(ETMISS).Fill(met_JERu->p.Pt(),tempWeight);
            }else{
               if(changeMET=="JERd"){
                  m.at(ETMISS).Fill(met_JERd->p.Pt(),tempWeight);
               }
            }
         }
      }
   }
   
   
   
   //m.at(ETMISS).Fill(ev.ETmiss, tempWeight);
   m.at(PT1).Fill(l1.p.Pt(), tempWeight);
   m.at(PT2).Fill(l2.p.Pt(), tempWeight);
   m.at(PT3).Fill(l3.p.Pt(), tempWeight);
   m.at(MLL).Fill((l1.vec+l2.vec).M(),tempWeight);   
   m.at(NPHOTONS).Fill(ev.selPhotons.size(),tempWeight);
   m.at(HT).Fill(ev.calcHt,tempWeight);
   m.at(NJETS).Fill(ev.selJets.size(),tempWeight);
   m.at(ETA1).Fill(fabs(l1.p.Eta()),tempWeight);
   m.at(ETA2).Fill(fabs(l2.p.Eta()),tempWeight);
   m.at(ETA3).Fill(fabs(l3.p.Eta()),tempWeight);
   m.at(PHI1).Fill(fabs(l1.p.Phi()),tempWeight);
   m.at(PHI2).Fill(fabs(l2.p.Phi()),tempWeight);
   m.at(PHI3).Fill(fabs(l3.p.Phi()),tempWeight);


   int jetSize=ev.selJets.size();








   if(!slimmed){
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
      m.at(NVTX).Fill(*nGoodVertices,tempWeight);
      m.at(GENHT).Fill(*genHt,tempWeight);
      m.at(ZPT).Fill((l1.vec+l2.vec).Pt(),tempWeight);
      m.at(ST).Fill(ev.pt1+ev.pt2,tempWeight);
      m.at(MT2).Fill(ev.MT2_val,tempWeight);

      m.at(DeltaEtaLL).Fill(fabs(l1.p.Eta()-l2.p.Eta()),tempWeight);
      m.at(DeltaPhiLL).Fill(fabs(l1.p.Phi()-l2.p.Phi()),tempWeight);
      m.at(DeltaRLL).Fill(fabs(l1.vec.DeltaR(l2.vec)),tempWeight);

      m.at(DeltaEtaLL_neg).Fill(l1.p.Eta()-l2.p.Eta(),tempWeight);
      m.at(DeltaPhiLL_neg).Fill(l1.p.Phi()-l2.p.Phi(),tempWeight);
      m.at(DeltaRLL_neg).Fill(l1.vec.DeltaR(l2.vec),tempWeight);   
      m.at(MTL3MET).Fill((l1.vec+ev.ETmiss_vec).Mt(),tempWeight);

      m.at(NElectrons).Fill(ev.selElectrons.size(),tempWeight);
      m.at(NMuons).Fill(ev.selMuons.size(),tempWeight);   
   }


   
   if(withPhoton){
      m.at(PTG1).Fill(ev.selPhotons.at(0).p.Pt(),tempWeight);
      m.at(PHIG1).Fill(ev.selPhotons.at(0).p.Phi(),tempWeight);
      m.at(ETAG1).Fill(ev.selPhotons.at(0).p.Eta(),tempWeight);
      
      
      
      if(!slimmed){
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
      
      }
      
      
      
      //if(FindGenPhotonMatch(ev.selPhotons.at(0))) m.at(genPhotonPT).Fill(GetGenPhotonMatch(ev.selPhotons.at(0)).p.Pt(),tempWeight);
      //
      //if(ev.evtHasGenPhotonVeto) if(FindGenPhotonMatch(ev.selPhotons.at(0))) m.at(genPhotonPT_Veto).Fill(GetGenPhotonMatch(ev.selPhotons.at(0)).p.Pt(),tempWeight);
      //if(ev.evtHasGenPhotonVeto) m.at(PTG1_Veto).Fill(ev.selPhotons.at(0).p.Pt(),tempWeight);
      //if(!ev.evtHasGenPhotonVeto) if(FindGenPhotonMatch(ev.selPhotons.at(0))) m.at(genPhotonPT_NoVeto).Fill(GetGenPhotonMatch(ev.selPhotons.at(0)).p.Pt(),tempWeight);
      //if(!ev.evtHasGenPhotonVeto) m.at(PTG1_NoVeto).Fill(ev.selPhotons.at(0).p.Pt(),tempWeight);
      
   }
}
void ntupleAnalyzerShrinked::Filler2D(selEvent& ev, map<Histograms2D,TH2F>& m,bool withPhoton){
   if(withPhoton){
      //m.at(PTG1).Fill(ev.selPhotons.at(0).p.Pt(),ev.totalWeight);
      m.at(ISRVFSR).Fill(ev.mll,(ev.l1+selectedEvent.l2+ev.selPhotons.at(0).vec).M(),ev.totalWeight);  
      m.at(PTGvsMLLG).Fill((ev.l1+selectedEvent.l2+ev.selPhotons.at(0).vec).M(),ev.selPhotons.at(0).p.Pt(),ev.totalWeight);  
   }
}

void ntupleAnalyzerShrinked::FillerTrigger(selEvent& ev, map<Histograms1D,TEfficiency>& m,bool withPhoton,bool TriggerBool){
   //m.at(ETMISS).Fill(TriggerBool,ev.ETmiss);
   //m.at(PT1).Fill(TriggerBool,ev.pt1);
   //m.at(PT2).Fill(TriggerBool,ev.pt2);
   //m.at(MLL).Fill(TriggerBool,ev.mll);   
   //m.at(NPHOTONS).Fill(TriggerBool,ev.selPhotons.size());
   //m.at(NVTX).Fill(TriggerBool,*nGoodVertices);
   //m.at(HT).Fill(TriggerBool,ev.calcHt);
   //m.at(GENHT).Fill(TriggerBool,*genHt);
   //m.at(NJETS).Fill(TriggerBool,ev.selJets.size());
   //m.at(ETA1).Fill(TriggerBool,fabs(ev.eta1));
   //m.at(ETA2).Fill(TriggerBool,fabs(ev.eta2));
   //m.at(PHI1).Fill(TriggerBool,fabs(ev.phi2));
   //m.at(PHI2).Fill(TriggerBool,fabs(ev.phi2));
   m.at(ETMISS).FillWeighted(TriggerBool,ev.totalWeight,ev.ETmiss);
   m.at(PT1).Fill(TriggerBool,ev.totalWeight,ev.pt1);
   m.at(PT2).Fill(TriggerBool,ev.totalWeight,ev.pt2);
   m.at(MLL).Fill(TriggerBool,ev.mll);   
   m.at(NPHOTONS).FillWeighted(TriggerBool,ev.totalWeight,ev.selPhotons.size());
   m.at(NVTX).FillWeighted(TriggerBool,ev.totalWeight,*nGoodVertices);
   m.at(HT).FillWeighted(TriggerBool,ev.totalWeight,ev.calcHt);
   m.at(GENHT).FillWeighted(TriggerBool,ev.totalWeight,*genHt);
   m.at(NJETS).FillWeighted(TriggerBool,ev.totalWeight,ev.selJets.size());
   m.at(ETA1).FillWeighted(TriggerBool,ev.totalWeight,fabs(ev.eta1));
   m.at(ETA2).FillWeighted(TriggerBool,ev.totalWeight,fabs(ev.eta2));
   m.at(PHI1).FillWeighted(TriggerBool,ev.totalWeight,fabs(ev.phi2));
   m.at(PHI2).FillWeighted(TriggerBool,ev.totalWeight,fabs(ev.phi2));
   if(withPhoton){
      m.at(PTG1).Fill(TriggerBool,ev.selPhotons.at(0).p.Pt());
      m.at(PHIG1).Fill(TriggerBool,fabs(ev.selPhotons.at(0).p.Phi()));
      m.at(ETAG1).Fill(TriggerBool,fabs(ev.selPhotons.at(0).p.Eta()));
      m.at(SIGMAIETAIETAG1).Fill(TriggerBool,ev.selPhotons.at(0).sigmaIetaIeta);
   }
}

//void ntupleAnalyzerShrinked::FillerSignal(selEvent& ev, map<Histograms1D,TH1F>& m, float divideFactor=1.){
void ntupleAnalyzerShrinked::FillerSignal(selEvent& ev, map<Histograms1D,TH1F>& m, float divideFactor){
   
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
   
   m.at(ETMISS).Fill(ev.ETmiss, tempWeight*1./(nGen/divideFactor));
}




void ntupleAnalyzerShrinked::InitAllHistos(){
   if (config_selectionsToProcessMap[DILEP]==true){
      h1Maps["dilepEE"]=InitHistograms(DILEP); 
      h1Maps["dilepMM"]=InitHistograms(DILEP);
      h1Maps["dilepEM"]=InitHistograms(DILEP);
      h1Maps["dilep"]=InitHistograms(DILEP);
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
   //if (config_selectionsToProcessMap[ControlRegionDY]==true){
      //h1Maps["CRDYEE"]=InitHistograms(ONZ);
      //h1Maps["CRDYMM"]=InitHistograms(ONZ); 
      //h1Maps["CRDYEM"]=InitHistograms(ONZ); 
      //h1Maps["CRDY"]=InitHistograms(ONZ); 
   //}
   if (config_selectionsToProcessMap[ControlRegionDY]==true){
      InitWeightHistos(cr1Maps,ONZ,"CRDY");
   }
   //if (config_selectionsToProcessMap[ControlRegionTT]==true){
      //h1Maps["CRTTEE"]=InitHistograms(ONZ);
      //h1Maps["CRTTMM"]=InitHistograms(ONZ); 
      //h1Maps["CRTTEM"]=InitHistograms(ONZ); 
      //h1Maps["CRTT"]=InitHistograms(ONZ); 
      //h1Maps["CRTT080EE"]=InitHistograms(ONZ);
      //h1Maps["CRTT080MM"]=InitHistograms(ONZ); 
      //h1Maps["CRTT080EM"]=InitHistograms(ONZ); 
      //h1Maps["CRTT080"]=InitHistograms(ONZ); 
      //h1Maps["CRTT80EE"]=InitHistograms(ONZ);
      //h1Maps["CRTT80MM"]=InitHistograms(ONZ); 
      //h1Maps["CRTT80EM"]=InitHistograms(ONZ); 
      //h1Maps["CRTT80"]=InitHistograms(ONZ); 
   //}
   if (config_selectionsToProcessMap[ControlRegionTT]==true){
      InitWeightHistos(cr1Maps,ONZ,"CRTT");
      InitWeightHistos(cr1Maps,ONZ,"CRTT080");
      InitWeightHistos(cr1Maps,ONZ,"CRTT80");
   }
   //if (config_selectionsToProcessMap[ControlRegionZZ]==true){
      //h1Maps["CRZZEE"]=InitHistograms(ControlRegionZZ);
      //h1Maps["CRZZMM"]=InitHistograms(ControlRegionZZ); 
      //h1Maps["CRZZEM"]=InitHistograms(ControlRegionZZ); 
      //h1Maps["CRZZ"]=InitHistograms(ControlRegionZZ); 
   //}
   if (config_selectionsToProcessMap[ControlRegionZZ]==true){
      InitWeightHistos(cr1Maps,ControlRegionZZ,"CRZZ");
   }
   //if (config_selectionsToProcessMap[ControlRegionWZ]==true){
      //h1Maps["CRWZEE"]=InitHistograms(ControlRegionZZ);
      //h1Maps["CRWZMM"]=InitHistograms(ControlRegionZZ); 
      //h1Maps["CRWZEM"]=InitHistograms(ControlRegionZZ); 
      //h1Maps["CRWZ"]=InitHistograms(ControlRegionZZ); 
   //}
   if (config_selectionsToProcessMap[ControlRegionWZ]==true){
      InitWeightHistos(cr1Maps,ControlRegionZZ,"CRWZ"); 
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
   //if (config_selectionsToProcessMap[SEL]==true){
      //h2Maps["selEE"]=Init2DHistograms(SEL);
      //h2Maps["selMM"]=Init2DHistograms(SEL); 
      //h2Maps["selEM"]=Init2DHistograms(SEL); 
      //h2Maps["sel"]=Init2DHistograms(SEL); 
   //}
   //if (config_selectionsToProcessMap[ONZ]==true){
      //h2Maps["onZEE"]=Init2DHistograms(ONZ);
      //h2Maps["onZMM"]=Init2DHistograms(ONZ); 
      //h2Maps["onZEM"]=Init2DHistograms(ONZ); 
      //h2Maps["onZ"]=Init2DHistograms(ONZ); 
   //}
   //if (config_selectionsToProcessMap[DILEP]==true){
      //h2Maps["dilepEE"]=Init2DHistograms(DILEP);
      //h2Maps["dilepMM"]=Init2DHistograms(DILEP); 
      //h2Maps["dilepEM"]=Init2DHistograms(DILEP); 
      //h2Maps["dilep"]=Init2DHistograms(DILEP); 
   //}
   //if (config_selectionsToProcessMap[EXO]==true){
      //h2Maps["exoEE"]=Init2DHistograms(EXO);
      //h2Maps["exoMM"]=Init2DHistograms(EXO); 
      //h2Maps["exoEM"]=Init2DHistograms(EXO); 
      //h2Maps["exo"]=Init2DHistograms(EXO); 
   //}
}

map<Histograms1D,TH1F> ntupleAnalyzerShrinked::InitCutFlowHistograms(const selectionType selection){
   map<Histograms1D,TH1F> cMap;
   cMap[CUTFLOW] = TH1F("","",5,0,5);
   cMap[CUTFLOW].Fill("triggered", 0);
   cMap[CUTFLOW].Fill("2leptons", 0);
   cMap[CUTFLOW].Fill("m50", 0);
   cMap[CUTFLOW].Fill("1photon", 0);
   cMap[CUTFLOW].Fill("Z", 0);
   return cMap;
}


map<Histograms1D,TH1F> ntupleAnalyzerShrinked::InitCutFlowHistograms_Fine(const selectionType selection){
   map<Histograms1D,TH1F> cMap;
   cMap[CUTFLOW_fine] = TH1F("","",18,0,18);
   cMap[CUTFLOW_fine].Fill("genZToLL",0);
   cMap[CUTFLOW_fine].Fill("triggered", 0);
   cMap[CUTFLOW_fine].Fill("triggeredMatched", 0);
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


void ntupleAnalyzerShrinked::InitCutFlowHistos(){
   c1Maps["cutFlow_onZEE"]=InitCutFlowHistograms(ONZ);
   c1Maps["cutFlow_onZMM"]=InitCutFlowHistograms(ONZ);
   c1Maps["cutFlow_onZEM"]=InitCutFlowHistograms(ONZ);
}


void ntupleAnalyzerShrinked::InitCutFlowHistos_Fine(){
   c1Maps["cutFlow_Fine_onZEE"]=InitCutFlowHistograms_Fine(ONZ);
   c1Maps["cutFlow_Fine_onZMM"]=InitCutFlowHistograms_Fine(ONZ);
   c1Maps["cutFlow_Fine_onZEM"]=InitCutFlowHistograms_Fine(ONZ);
}


map<Histograms1D,TH1F> ntupleAnalyzerShrinked::InitHistograms(const selectionType selection_){
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
   hMap[NBJETS] = TH1F("","n_{BJets}",20,0,20);
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
   hMap[WEIGHT_PDF]=TH1F("","weight",10000,0,10);


   hMap[DeltaEtaLL] = TH1F("", ";#Delta#Eta_{ll}", 24000, -12., 12.);
   hMap[DeltaPhiLL] = TH1F("", ";#Delta#Phi_{ll}", 24000, -12., 12.);
   hMap[DeltaRLL] = TH1F("", ";#DeltaR_{ll}", 24000, -12, 12.);

   hMap[DeltaEtaLL_neg] = TH1F("", ";#Delta#Eta_{ll}", 24000, -12., 12.);
   hMap[DeltaPhiLL_neg] = TH1F("", ";#Delta#Phi_{ll}", 24000, -12., 12.);
   hMap[DeltaRLL_neg] = TH1F("", ";#DeltaR_{ll}", 24000, -12., 12.);
         //
   if(selection_!=ControlRegionZZ){
   hMap[DeltaPhiLLMet] = TH1F("", "#Delta#Phi_{ll,MET}", 24000, -12., 12.);
   hMap[DeltaEtaLLMet] = TH1F("", ";#Delta#Eta_{ll,MET}", 24000, -12., 12.);
   hMap[DeltaRLLMet] = TH1F("", "#DeltaR_{ll,MET}", 24000, -12., 12.);
   hMap[MTLL] = TH1F("", ";m_{T}^{ll}", 5000, 0, 5000); 
   }

   if(selection_==ControlRegionZZ){
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
   
   if ((selection_==PHOTON)||(selection_==SEL)||(selection_==ONZ)){
      
      //hMap[FakeElectron] = TH1F("", "", 2, 0, 2);
      //hMap[FakeJet] = TH1F("", "", 2, 0, 2);
      //hMap[FakePhoton] = TH1F("", "", 2, 0, 2);
      hMap[Fakes] = TH1F("", "", 4, 0, 4);
      
      
      hMap[PTG1] = TH1F("", ";#it{p}_{T}^{#gamma 1} (GeV)", 5000, 0, 5000);
      hMap[ETAG1] = TH1F("", ";|#eta_{#gamma 1}|", 260, 0, 2.6);
      hMap[PHIG1] = TH1F("", ";|#phi_{#gamma 1}|", 350, 0, 3.5);
      hMap[SIGMAIETAIETAG1] = TH1F("", ";#sigma_{i#etai#eta}^{#gamma 1}", 400, 0, 0.04);
      hMap[SIGMAIPHIIPHIG1] = TH1F("", ";#sigma_{i#phii#phi}^{#gamma 1}", 2000, 0, 0.2);
      hMap[R9] = TH1F("", ";r9", 1500, 0, 1.5);
      hMap[HOVERE] = TH1F("", ";H/E", 1000, 0, 0.1);
      hMap[DELTARGL1] = TH1F("", ";#DeltaR_{l1,#gamma}", 24000, -12., 12.);
      hMap[DELTARGL2] = TH1F("", ";#DeltaR_{l2,#gamma}", 24000, -12., 12.);
      hMap[DeltaRLLG] = TH1F("", ";#DeltaR_{ll,#gamma}", 24000, -12., 12.);
      hMap[DeltaEtaLLG] = TH1F("", ";#Delta#Eta_{ll,#gamma}", 2400, -12., 12.);
      hMap[DeltaPhiLLG] = TH1F("", ";#Delta#Phi_{ll,#gamma}", 2400, -12., 12.);
      hMap[DeltaPhiGMet] = TH1F("", ";#Delta#Phi_{met,#gamma}", 2400, -12., 6.);
      hMap[DeltaRGMet] = TH1F("", ";#DeltaR_{met,#gamma}", 2400, -12., 6.);
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


map<Histograms1D,TH1F> ntupleAnalyzerShrinked::InitHistogramsSlimmed(const selectionType selection_){
   map<Histograms1D,TH1F> hMap;

   hMap[ETMISS] = TH1F("", ";#it{p}_{T}^{miss} (GeV)", 5000, 0, 5000);
   hMap[PT1] = TH1F("", ";#it{p}_{T}^{leading} (GeV)", 5000, 0, 5000);
   hMap[PT2] = TH1F("", ";#it{p}_{T}^{trailing} (GeV)", 5000, 0, 5000);
   hMap[MLL] = TH1F("", ";#it{m}_{ll} (GeV)", 5000, 0, 5000);
   hMap[NPHOTONS] = TH1F("","n_#gamma",10,0,10);
   //hMap[NVTX] = TH1F("","n_{Vtx}",60,0,60);
   hMap[HT] = TH1F("","#it{H}_{T} (GeV)",5000,0,5000);
   //hMap[GENHT] = TH1F("","#it{H}_{T}^{gen} (GeV)",5000,0,5000);
   hMap[NJETS] = TH1F("","n_{Jets}",20,0,20);
   //hMap[NBJETS] = TH1F("","n_{BJets}",20,0,20);
   hMap[ETA1] = TH1F("", ";|#eta_{trailing}|", 260, 0, 2.6);
   hMap[ETA2] = TH1F("", ";|#eta_{leading}|", 260, 0, 2.6);
   hMap[PHI1] = TH1F("", ";|#phi_{trailing}|", 350, 0, 3.5);
   hMap[PHI2] = TH1F("", ";|#phi_{leading}|", 350, 0, 3.5);
   //hMap[ZPT] = TH1F("", ";Z_{p_T}", 5000, 0, 5000);
   //hMap[ST] = TH1F("", ";S_T", 5000, 0, 5000.);
   //hMap[MT2] = TH1F("", ";M_{T2}", 50000, 0, 5000.);
   //hMap[VetoCompare] = TH1F("", "", 2, 0, 2);
   //hMap[NElectrons] = TH1F("","nElectrons",20,0,20);
   //hMap[NMuons] = TH1F("","nMuons",20,0,20);

   //hMap[JetPt1] = TH1F("", ";#it{p}_{T}^{Jet 1} (GeV)", 5000, 0, 5000);
   //hMap[JetPt2] = TH1F("", ";#it{p}_{T}^{Jet 2} (GeV)", 5000, 0, 5000);
   //hMap[JetPt3] = TH1F("", ";#it{p}_{T}^{Jet 3} (GeV)", 5000, 0, 5000);
   //hMap[JetPt4] = TH1F("", ";#it{p}_{T}^{Jet 4} (GeV)", 5000, 0, 5000);
   //hMap[JetPhi1] = TH1F("", ";|#phi_{Jet 1}|", 350, 0, 3.5);
   //hMap[JetPhi2] = TH1F("", ";|#phi_{Jet 2}|", 350, 0, 3.5);
   //hMap[JetPhi3] = TH1F("", ";|#phi_{Jet 3}|", 350, 0, 3.5);
   //hMap[JetPhi4] = TH1F("", ";|#phi_{Jet 4}|", 350, 0, 3.5);
   //hMap[JetEta1] = TH1F("", ";|#eta_{Jet 1}|", 260, 0, 2.6);
   //hMap[JetEta2] = TH1F("", ";|#eta_{Jet 2}|", 260, 0, 2.6);
   //hMap[JetEta3] = TH1F("", ";|#eta_{Jet 3}|", 260, 0, 2.6);
   //hMap[JetEta4] = TH1F("", ";|#eta_{Jet 4}|", 260, 0, 2.6);

   hMap[WEIGHT_NISR] = TH1F("","weight",10000,0,10);
   hMap[WEIGHT_TOPPT] = TH1F("","weight",10000,0,10);
   hMap[WEIGHT_EWKINOPAIRPT] = TH1F("","weight",10000,0,10);
   hMap[WEIGHT_LEPTONPAIRPT] = TH1F("","weight",10000,0,10);
   hMap[WEIGHT_PDF]=TH1F("","weight",10000,0,10);


   //hMap[DeltaEtaLL] = TH1F("", ";#Delta#Eta_{ll}", 24000, -12., 12.);
   //hMap[DeltaPhiLL] = TH1F("", ";#Delta#Phi_{ll}", 24000, -12., 12.);
   //hMap[DeltaRLL] = TH1F("", ";#DeltaR_{ll}", 24000, -12, 12.);

   //hMap[DeltaEtaLL_neg] = TH1F("", ";#Delta#Eta_{ll}", 24000, -12., 12.);
   //hMap[DeltaPhiLL_neg] = TH1F("", ";#Delta#Phi_{ll}", 24000, -12., 12.);
   //hMap[DeltaRLL_neg] = TH1F("", ";#DeltaR_{ll}", 24000, -12., 12.);
         //
   //if(selection_!=ControlRegionZZ){
   //hMap[DeltaPhiLLMet] = TH1F("", "#Delta#Phi_{ll,MET}", 24000, -12., 12.);
   //hMap[DeltaEtaLLMet] = TH1F("", ";#Delta#Eta_{ll,MET}", 24000, -12., 12.);
   //hMap[DeltaRLLMet] = TH1F("", "#DeltaR_{ll,MET}", 24000, -12., 12.);
   //hMap[MTLL] = TH1F("", ";m_{T}^{ll}", 5000, 0, 5000); 
   //}

   if(selection_==ControlRegionZZ){
      hMap[PT3] = TH1F("", ";#it{p}_{T}^{trailing} (GeV)", 5000, 0, 5000);
      hMap[PT4] = TH1F("", ";#it{p}_{T}^{trailing} (GeV)", 5000, 0, 5000);
      hMap[ETA3] = TH1F("", ";|#eta_{leading}|", 260, 0, 2.6);
      hMap[ETA4] = TH1F("", ";|#eta_{leading}|", 260, 0, 2.6);
      hMap[PHI3] = TH1F("", ";|#phi_{leading}|", 350, 0, 3.5);
      hMap[PHI4] = TH1F("", ";|#phi_{leading}|", 350, 0, 3.5);
      //hMap[ZPT2] = TH1F("", ";Z_{p_T}", 5000, 0, 5000);
      hMap[MLL2] = TH1F("", ";#it{m}_{ll} (GeV)", 5000, 0, 5000);
      //hMap[MTL3MET] = TH1F("", ";#it{m}_{T}^{l2,Met} (GeV)", 5000, 0, 5000);
      //hMap[MLLLL] = TH1F("", ";#it{m}_{llll} (GeV)", 5000, 0, 5000);
      
   }
   
   if ((selection_==PHOTON)||(selection_==SEL)||(selection_==ONZ)){
      hMap[PTG1] = TH1F("", ";#it{p}_{T}^{#gamma 1} (GeV)", 5000, 0, 5000);
      hMap[ETAG1] = TH1F("", ";|#eta_{#gamma 1}|", 260, 0, 2.6);
      hMap[PHIG1] = TH1F("", ";|#phi_{#gamma 1}|", 350, 0, 3.5);
      
      //hMap[FakeElectron] = TH1F("", "", 2, 0, 2);
      //hMap[FakeJet] = TH1F("", "", 2, 0, 2);
      //hMap[FakePhoton] = TH1F("", "", 2, 0, 2);
      hMap[Fakes] = TH1F("", "", 3, 0, 3);
      
      //hMap[SIGMAIETAIETAG1] = TH1F("", ";#sigma_{i#etai#eta}^{#gamma 1}", 400, 0, 0.04);
      //hMap[SIGMAIPHIIPHIG1] = TH1F("", ";#sigma_{i#phii#phi}^{#gamma 1}", 2000, 0, 0.2);
      //hMap[R9] = TH1F("", ";r9", 1500, 0, 1.5);
      //hMap[HOVERE] = TH1F("", ";H/E", 1000, 0, 0.1);
      //hMap[DELTARGL1] = TH1F("", ";#DeltaR_{l1,#gamma}", 24000, -12., 12.);
      //hMap[DELTARGL2] = TH1F("", ";#DeltaR_{l2,#gamma}", 24000, -12., 12.);
      //hMap[DeltaRLLG] = TH1F("", ";#DeltaR_{ll,#gamma}", 24000, -12., 12.);
      //hMap[DeltaEtaLLG] = TH1F("", ";#Delta#Eta_{ll,#gamma}", 2400, -12., 12.);
      //hMap[DeltaPhiLLG] = TH1F("", ";#Delta#Phi_{ll,#gamma}", 2400, -12., 12.);
      //hMap[DeltaPhiGMet] = TH1F("", ";#Delta#Phi_{met,#gamma}", 2400, -12., 6.);
      //hMap[DeltaRGMet] = TH1F("", ";#DeltaR_{met,#gamma}", 2400, -12., 6.);
      //hMap[STG] = TH1F("", ";S_T", 5000, 0, 5000.);
      //hMap[STMET] = TH1F("", ";S_T + #it{p}_{T}^{miss} (GeV)", 5000, 0, 5000.);   
      //hMap[MTLLG] = TH1F("", ";m_{T}^{ll#gamma}", 5000, 0, 5000);
      //hMap[MTL1MET] = TH1F("", ";m_{T}^{l1,met}", 5000, 0, 5000);
      //hMap[MTL2MET] = TH1F("", ";m_{T}^{l2,met}", 5000, 0, 5000);
      //hMap[MTGMET] = TH1F("", ";m_{T}^{#gamma,met}", 5000, 0, 5000);
      //hMap[MTLLMET] = TH1F("", ";m_{T}^{ll,met}", 5000, 0, 5000);
      //hMap[MTLLGMET] = TH1F("", ";m_{T}^{ll#gamma,met}", 5000, 0, 5000);
      //hMap[MOTHERID] = TH1F("", ";ID_{mother}", 200, 0, 200);
      //hMap[MLLG] = TH1F("", ";ID_{ll#gamma}", 5000, 0, 5000);
      //hMap[PT_llg] = TH1F("", ";p_{T}^{ll#gamma}", 5000, 0, 5000);
      //hMap[MZG_exo] = TH1F("", ";m_{Z#gamma}", 5000, 0, 5000);
      //hMap[gammaMotherID] = TH1F("", ";motherID_{#gamma}", 5000, 0, 5000);
      //hMap[genPhotonPT] = TH1F("", ";gen p_T^{#gamma,matched}", 5000, 0, 5000);
      //hMap[genPhotonPT_Veto] = TH1F("", ";gen p_T^{#gamma,matched,veto}", 5000, 0, 5000);
      //hMap[PTG1_Veto] = TH1F("", ";gen p_T^{#gamma,veto}", 5000, 0, 5000);
      //hMap[genPhotonPT_NoVeto] = TH1F("", ";gen p_T^{#gamma,matched,Noveto}", 5000, 0, 5000);
      //hMap[PTG1_NoVeto] = TH1F("", ";gen p_T^{#gamma,Noveto}", 5000, 0, 5000);
      //hMap[VetoCompare] = TH1F("", "", 2, 0, 2);
   }
   return hMap;
}

map<Histograms2D,TH2F> ntupleAnalyzerShrinked::Init2DHistograms(const selectionType selection_){
   map<Histograms2D,TH2F> h2Map;
   if ((selection_==PHOTON)||(selection_==SEL)||(selection_==ONZ)||(selection_==EXO) ){
      h2Map[PTGvsMLLG] = TH2F("",";#it{m}_{ll#gamma};#it{p}_{T}^{gamma}",1000,0,1000,1000,0,1000);   
      h2Map[ISRVFSR] = TH2F("",";#it{m}_{ll};#it{m}{ll#gamma}",200,0,200,500,0,500);   
   }
   return h2Map;
}





map<Histograms1D,TH1F> ntupleAnalyzerShrinked::InitSignalScanHistograms(const selectionType selection_){
   map<Histograms1D,TH1F> sMap;
   
   sMap[WEIGHT_NISR] = TH1F("","weight",10000,0,10);
   sMap[WEIGHT_TOPPT] = TH1F("","weight",10000,0,10);
   sMap[WEIGHT_EWKINOPAIRPT] = TH1F("","weight",10000,0,10);
   sMap[WEIGHT_LEPTONPAIRPT] = TH1F("","weight",10000,0,10);
   
   
   sMap[ETMISS] = TH1F("", ";#it{p}_{T}^{miss} (GeV)", 5000, 0, 5000);
   
   return sMap;
}

void ntupleAnalyzerShrinked::InitSignalScanHistos(string masspoint){
   s1Maps[masspoint]["sig"]=InitSignalScanHistograms(ONZ);
   s1Maps[masspoint]["sigEE"]=InitSignalScanHistograms(ONZ);
   s1Maps[masspoint]["sigMM"]=InitSignalScanHistograms(ONZ);
   s1Maps[masspoint]["sigEM"]=InitSignalScanHistograms(ONZ);
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
   }
}

void ntupleAnalyzerShrinked::InitWeightHistos(map<string,map<string,map<Histograms1D,TH1F>>>& map_,const selectionType selection_, string name_){
   //map_[name_+"EE"]["nom"]=InitHistograms(selection_);
   //map_[name_+"MM"]["nom"]=InitHistograms(selection_);
   //map_[name_+"EM"]["nom"]=InitHistograms(selection_);
   //map_[name_]["nom"]=InitHistograms(selection_);
   //map_[name_+"EE"]["JESu"]=InitHistograms(selection_);
   //map_[name_+"MM"]["JESu"]=InitHistograms(selection_);
   //map_[name_+"EM"]["JESu"]=InitHistograms(selection_);
   //map_[name_]["JESu"]=InitHistograms(selection_);
   //map_[name_+"EE"]["JESd"]=InitHistograms(selection_);
   //map_[name_+"MM"]["JESd"]=InitHistograms(selection_);
   //map_[name_+"EM"]["JESd"]=InitHistograms(selection_);
   //map_[name_]["JESd"]=InitHistograms(selection_);
   //map_[name_+"EE"]["JERu"]=InitHistograms(selection_);
   //map_[name_+"MM"]["JERu"]=InitHistograms(selection_);
   //map_[name_+"EM"]["JERu"]=InitHistograms(selection_);
   //map_[name_]["JERu"]=InitHistograms(selection_);
   //map_[name_+"EE"]["JERd"]=InitHistograms(selection_);
   //map_[name_+"MM"]["JERd"]=InitHistograms(selection_);
   //map_[name_+"EM"]["JERd"]=InitHistograms(selection_);
   //map_[name_]["JERd"]=InitHistograms(selection_);
   
   //for(int i=0; i<nWeights; i++){
      //map_[name_+"EE"][to_string(i)]=InitHistograms(selection_);
      //map_[name_+"MM"][to_string(i)]=InitHistograms(selection_);
      //map_[name_+"EM"][to_string(i)]=InitHistograms(selection_);
      //map_[name_][to_string(i)]=InitHistograms(selection_);
   //}
   map_[name_+"EE"]["nom"]=InitHistogramsSlimmed(selection_);
   map_[name_+"MM"]["nom"]=InitHistogramsSlimmed(selection_);
   map_[name_+"EM"]["nom"]=InitHistogramsSlimmed(selection_);
   map_[name_]["nom"]=InitHistogramsSlimmed(selection_);
   map_[name_+"EE"]["JESu"]=InitHistogramsSlimmed(selection_);
   map_[name_+"MM"]["JESu"]=InitHistogramsSlimmed(selection_);
   map_[name_+"EM"]["JESu"]=InitHistogramsSlimmed(selection_);
   map_[name_]["JESu"]=InitHistogramsSlimmed(selection_);
   map_[name_+"EE"]["JESd"]=InitHistogramsSlimmed(selection_);
   map_[name_+"MM"]["JESd"]=InitHistogramsSlimmed(selection_);
   map_[name_+"EM"]["JESd"]=InitHistogramsSlimmed(selection_);
   map_[name_]["JESd"]=InitHistogramsSlimmed(selection_);
   map_[name_+"EE"]["JERu"]=InitHistogramsSlimmed(selection_);
   map_[name_+"MM"]["JERu"]=InitHistogramsSlimmed(selection_);
   map_[name_+"EM"]["JERu"]=InitHistogramsSlimmed(selection_);
   map_[name_]["JERu"]=InitHistogramsSlimmed(selection_);
   map_[name_+"EE"]["JERd"]=InitHistogramsSlimmed(selection_);
   map_[name_+"MM"]["JERd"]=InitHistogramsSlimmed(selection_);
   map_[name_+"EM"]["JERd"]=InitHistogramsSlimmed(selection_);
   map_[name_]["JERd"]=InitHistogramsSlimmed(selection_);
   
   for(int i=0; i<nWeights; i++){
      map_[name_+"EE"][to_string(i)]=InitHistogramsSlimmed(selection_);
      map_[name_+"MM"][to_string(i)]=InitHistogramsSlimmed(selection_);
      map_[name_+"EM"][to_string(i)]=InitHistogramsSlimmed(selection_);
      map_[name_][to_string(i)]=InitHistogramsSlimmed(selection_);
   }
}

//void ntupleAnalyzerShrinked::InitHistosForWeights(map<string,map<string,map<Histograms1D,TH1F>>>& map_,const selectionType selection_, string name1_, string name2_){
   //map_[name1_][name2_]=InitHistograms(selection_);
//}

map<Histograms1D,TEfficiency> ntupleAnalyzerShrinked::InitTriggerStudies(const selectionType selection_){
    
   map<Histograms1D,TEfficiency> hMap;
    
   hMap[ETMISS] = TEfficiency("", ";#it{p}_{T}^{miss} (GeV)", 200, 0, 1000);
   hMap[ETMISS].SetUseWeightedEvents();
   hMap[PT1] = TEfficiency("", ";#it{p}_{T}^{leading} (GeV)", 1000, 0, 1000);
   hMap[PT1].SetUseWeightedEvents();
   hMap[PT2] = TEfficiency("", ";#it{p}_{T}^{trailing} (GeV)", 1000, 0, 1000);
   hMap[PT2].SetUseWeightedEvents();
   hMap[MLL] = TEfficiency("", ";#it{m}_{ll} (GeV)", 200, 0, 1000);
   hMap[MLL].SetUseWeightedEvents();
   hMap[NPHOTONS] = TEfficiency("","n_#gamma",10,0,10);
   hMap[NPHOTONS].SetUseWeightedEvents();
   hMap[NVTX] = TEfficiency("","n_{Vtx}",60,0,60);
   hMap[NVTX].SetUseWeightedEvents();
   hMap[HT] = TEfficiency("","#it{H}_{T} (GeV)",200,0,1000);
   hMap[HT].SetUseWeightedEvents();
   hMap[GENHT] = TEfficiency("","#it{H}_{T}^{gen} (GeV)",200,0,1000);
   hMap[GENHT].SetUseWeightedEvents();
   hMap[NJETS] = TEfficiency("","n_{Jets}",20,0,20);
   hMap[NJETS].SetUseWeightedEvents();
   hMap[ETA1] = TEfficiency("", ";|#eta_{leading}|", 260, 0, 2.6);
   hMap[ETA1].SetUseWeightedEvents();
   hMap[ETA2] = TEfficiency("", ";|#eta_{trailing}|", 260, 0, 2.6);
   hMap[ETA2].SetUseWeightedEvents();
   hMap[PHI1] = TEfficiency("", ";|#phi_{leading}|", 350, 0, 3.5);
   hMap[PHI1].SetUseWeightedEvents();
   hMap[PHI2] = TEfficiency("", ";|#phi_{trailing}|", 350, 0, 3.5);
   hMap[PHI2].SetUseWeightedEvents();

   if ((selection_==TRIGONZ)||(selection_==TRIGSEL)||(selection_==TRIGSEL_ptcuts)||(selection_==TRIGSEL_ptcuts)){
      hMap[PTG1] = TEfficiency("", ";#it{p}_{T}^{#gamma 1} (GeV)", 200, 0, 1000);
      hMap[PTG1].SetUseWeightedEvents();
      hMap[ETAG1] = TEfficiency("", ";|#eta_{#gamma 1}|", 260, 0, 2.6);
      hMap[ETAG1].SetUseWeightedEvents();
      hMap[PHIG1] = TEfficiency("", ";|#phi_{#gamma 1}|", 350, 0, 3.5);
      hMap[PHIG1].SetUseWeightedEvents();
      hMap[SIGMAIETAIETAG1] = TEfficiency("", ";#sigma_{i#etai#eta}^{#gamma 1}", 400, 0, 0.04);
      hMap[SIGMAIETAIETAG1].SetUseWeightedEvents();
   }

   return hMap;

}



void ntupleAnalyzerShrinked::InitTriggerStudiesHistos(){
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
   eff1Maps["trigSelEE_ptcuts"]=InitTriggerStudies(TRIGSEL_ptcuts);
   eff1Maps["trigSelMM_ptcuts"]=InitTriggerStudies(TRIGSEL_ptcuts);
   eff1Maps["trigSelEM_ptcuts"]=InitTriggerStudies(TRIGSEL_ptcuts);
   eff1Maps["trigOnZEE_ptcuts"]=InitTriggerStudies(TRIGONZ_ptcuts);
   eff1Maps["trigOnZMM_ptcuts"]=InitTriggerStudies(TRIGONZ_ptcuts);
   eff1Maps["trigOnZEM_ptcuts"]=InitTriggerStudies(TRIGONZ_ptcuts);

}


void ntupleAnalyzerShrinked::InitScaleFactorsFinal(){
   
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

   PhotonIDWeighter = Weighter("newScaleFactors/photon/scaleFactor_photon_ID.root","EGamma_SF2D"); // pt vs eta
   PhotonIDWeighter.fillOverflow2d();
   
   PhotonConversionWeighter = Weighter("newScaleFactors/photon/scaleFactor_photon_conv.root","Scaling_Factors_HasPix_R9 Inclusive"); // pt vs eta
   PhotonConversionWeighter.fillOverflow2d();
}




float ntupleAnalyzerShrinked::GetScaleFactorAndErrorFinal(float pt, float eta,bool isFastSim,bool isEle){
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


float ntupleAnalyzerShrinked::GetScaleFactorAndErrorPhotons(vector<selPhoton>& vecGamma){
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




bool ntupleAnalyzerShrinked::GenPhotonVeto(const int a){
   bool isGenPhoton=false;
   bool vetoPt130ZG=false;
   if(!isData){
      if(noPromptPhotons || isZGammaInclusive){
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
            }
         }
      }
   }
   return (isGenPhoton && noPromptPhotons)||(vetoPt130ZG);
}


bool ntupleAnalyzerShrinked::CheckParticles(){
   return ((muons->size()>=2.) || (electrons->size()>=2.) || ((muons->size()>=1)&&(electrons->size()>=1)) );      
}
bool ntupleAnalyzerShrinked::Check2Ele(){
   return (electrons->size()>=2.);
}
bool ntupleAnalyzerShrinked::Check2Mu(){
   return (muons->size()>=2.);
}
bool ntupleAnalyzerShrinked::CheckEMu(){
   return ((muons->size()>=1.) && (electrons->size()>=1.));
}





void ntupleAnalyzerShrinked::clearCutFlowMap(){
   decisionMapCutFlowFine[TRIGGERED]=false;
   decisionMapCutFlowFine[TRIGGEREDMATCHED]=false;
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
   decisionMapCutFlowFine_weight[TRIGGEREDMATCHED]=totalWeight;
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








bool ntupleAnalyzerShrinked::Cleaning() //return if event is valid (2 Leptons exist)
{
   if(!isSignal){
      selectedEvent.trigDiEle = *hlt_ele17_ele12_iso || *hlt_ele23_ele12_iso || *hlt_doubleEle33 || *hlt_doubleEle33_mw;
      selectedEvent.trigDiMu  = *hlt_mu17_mu8_iso || *hlt_mu17_tkMu8_iso || *hlt_mu17_mu8_iso_dz || *hlt_mu17_tkMu8_iso_dz
                  || *hlt_tkMu17_tkMu8_iso_dz || *hlt_mu27_tkMu8 || *hlt_mu30_tkMu11;
      selectedEvent.trigMuEle = *hlt_mu17_ele12_iso || *hlt_mu23_ele8_iso || *hlt_mu23_ele8_iso_dz || *hlt_mu23_ele12_iso
                  || *hlt_mu23_ele12_iso_dz || *hlt_mu8_ele17_iso || *hlt_mu8_ele23_iso
                  || *hlt_mu8_ele23_iso_dz || *hlt_mu12_ele23_iso || *hlt_mu12_ele23_iso_dz
                  || *hlt_mu30_ele30 || *hlt_mu33_ele33;
      selectedEvent.trigHt    = *hlt_ht200 || *hlt_ht250 || *hlt_ht300 || *hlt_ht350 || *hlt_ht400
                  || *hlt_ht475 || *hlt_ht600 || *hlt_ht650 || *hlt_ht800;
   }else{
      selectedEvent.trigDiEle = true;
      selectedEvent.trigDiMu = true;
      selectedEvent.trigMuEle = true;
      selectedEvent.trigHt = true;
   }
   

   if (isData){
      bool isEESelection = inputName.find("DoubleEG") != string::npos;
      bool isMuMuSelection = inputName.find("DoubleMuon") != string::npos;
      bool isEMuSelection = inputName.find("MuonEG") != string::npos;
      bool isHTSelection = inputName.find("JetHT")!= string::npos;
      auto missingET = *met;
      
      float pt1_ = 0.;
      float pt2_ = 0.;
      //float tol = 0.0001;
      float tol = 0.00005;
  
      std::string leptonFlavor1 = "N";
      std::string leptonFlavor2 = "N";
      
     for(auto& it : *electrons){
        if ((it.p.Pt()-tol) > pt1_){
           pt1_ = it.p.Pt();
           leptonFlavor1 = "E";
        }
      }
      for(auto& it : *muons){
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
        if ((it.p.Pt() < (pt1_-tol)) && ((it.p.Pt()-tol) > pt2_)){
           pt2_ = it.p.Pt();
           leptonFlavor2 = "M";
        }
      }
      if (leptonFlavor1 == "E" && leptonFlavor2 == "E") {
         if(isEESelection && Check2Ele() && selectedEvent.trigDiEle){ CalculateVariables(electrons->at(0),electrons->at(1),E); return true;}
      }
      else if (leptonFlavor1 == "M" && leptonFlavor2 == "M") {
         if(isMuMuSelection && Check2Mu() && selectedEvent.trigDiMu){ CalculateVariables(muons->at(0),muons->at(1),M); return true;}      
      }
      else if (leptonFlavor1 == "E" && leptonFlavor2 == "M") {
         if(isEMuSelection && CheckEMu() && selectedEvent.trigMuEle){ CalculateVariables(electrons->at(0),muons->at(0),EM); return true;}      
      }
      else if (leptonFlavor1 == "M" && leptonFlavor2 == "E") {
         if(isEMuSelection && CheckEMu() && selectedEvent.trigMuEle){ CalculateVariables(muons->at(0),electrons->at(0),ME); return true;}      
      }

      return false;
      }else{
         auto missingET = *met;

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
           if ((it.p.Pt() < (pt1_-tol)) && ((it.p.Pt()-tol) > pt2_)){
              pt2_ = it.p.Pt();
              leptonFlavor2 = "M";
           }
         }
         //cout<<leptonFlavor1<<leptonFlavor2<<endl;
         if (leptonFlavor1 == "E" && leptonFlavor2 == "E") {
            if(Check2Ele() && selectedEvent.trigDiEle){CalculateVariables(electrons->at(0),electrons->at(1),E); return true;}
         }
         else if (leptonFlavor1 == "M" && leptonFlavor2 == "M") {
            if(Check2Mu() && selectedEvent.trigDiMu){ CalculateVariables(muons->at(0),muons->at(1),M); return true;}
         }
         else if (leptonFlavor1 == "E" && leptonFlavor2 == "M") {
            if(CheckEMu() && selectedEvent.trigMuEle){ CalculateVariables(electrons->at(0),muons->at(0),EM); return true;}     
         }
         else if (leptonFlavor1 == "M" && leptonFlavor2 == "E") {
            if(CheckEMu() && selectedEvent.trigMuEle){ CalculateVariables(muons->at(0),electrons->at(0),ME); return true;}     
         }
         return false;
   }
   return false;
}


bool ntupleAnalyzerShrinked::CleaningTriggerStudies() //return if event is valid (2 Leptons exist)
{
   if(!isSignal){
      selectedEvent.trigDiEle = *hlt_ele17_ele12_iso || *hlt_ele23_ele12_iso || *hlt_doubleEle33 || *hlt_doubleEle33_mw;
      selectedEvent.trigDiMu  = *hlt_mu17_mu8_iso || *hlt_mu17_tkMu8_iso || *hlt_mu17_mu8_iso_dz || *hlt_mu17_tkMu8_iso_dz
                  || *hlt_tkMu17_tkMu8_iso_dz || *hlt_mu27_tkMu8 || *hlt_mu30_tkMu11;
      selectedEvent.trigMuEle = *hlt_mu17_ele12_iso || *hlt_mu23_ele8_iso || *hlt_mu23_ele8_iso_dz || *hlt_mu23_ele12_iso
                  || *hlt_mu23_ele12_iso_dz || *hlt_mu8_ele17_iso || *hlt_mu8_ele23_iso
                  || *hlt_mu8_ele23_iso_dz || *hlt_mu12_ele23_iso || *hlt_mu12_ele23_iso_dz
                  || *hlt_mu30_ele30 || *hlt_mu33_ele33;
      selectedEvent.trigHt    = *hlt_ht200 || *hlt_ht250 || *hlt_ht300 || *hlt_ht350 || *hlt_ht400
                  || *hlt_ht475 || *hlt_ht600 || *hlt_ht650 || *hlt_ht800;
   }else{
      selectedEvent.trigDiEle = true;
      selectedEvent.trigDiMu = true;
      selectedEvent.trigMuEle = false;
      selectedEvent.trigHt = true;
   }
   

   if (isData){
      //find Selection (Dataset)
      bool isEESelection = inputName.find("DoubleEG") != string::npos;
      bool isMuMuSelection = inputName.find("DoubleMuon") != string::npos;
      bool isEMuSelection = inputName.find("MuonEG") != string::npos;
      bool isHTSelection = inputName.find("JetHT")!= string::npos;
      auto missingET = *met;
        

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
        if ((it.p.Pt() < (pt1_-tol)) && ((it.p.Pt()-tol) > pt2_)){
           pt2_ = it.p.Pt();
           leptonFlavor2 = "M";
        }
      }
      

      
      if (leptonFlavor1 == "E" && leptonFlavor2 == "E") {
         if(isHTSelection && Check2Ele() && selectedEvent.trigHt){ CalculateVariables(electrons->at(0),electrons->at(1),E); return true;}
      }
      else if (leptonFlavor1 == "M" && leptonFlavor2 == "M") {
         if(isHTSelection && Check2Mu() && selectedEvent.trigHt){ CalculateVariables(muons->at(0),muons->at(1),M); return true;}      
      }
      else if (leptonFlavor1 == "E" && leptonFlavor2 == "M") {
         if(isHTSelection && CheckEMu() && selectedEvent.trigHt){ CalculateVariables(electrons->at(0),muons->at(0),EM); return true;}      
      }
      else if (leptonFlavor1 == "M" && leptonFlavor2 == "E") {
         if(isHTSelection && CheckEMu() && selectedEvent.trigHt){ CalculateVariables(muons->at(0),electrons->at(0),ME); return true;}      
      }
      return false;      

   }else{
      auto missingET = *met;
      
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
           if ((it.p.Pt() < (pt1_-tol)) && ((it.p.Pt()-tol) > pt2_)){
              pt2_ = it.p.Pt();
              leptonFlavor2 = "M";
           }
         }
         if (leptonFlavor1 == "E" && leptonFlavor2 == "E") {
            if(Check2Ele() && selectedEvent.trigHt){CalculateVariables(electrons->at(0),electrons->at(1),E); return true;}
         }
         else if (leptonFlavor1 == "M" && leptonFlavor2 == "M") {
            if(Check2Mu() && selectedEvent.trigHt){ CalculateVariables(muons->at(0),muons->at(1),M); return true;}
         }
         else if (leptonFlavor1 == "E" && leptonFlavor2 == "M") {
            if(CheckEMu() && selectedEvent.trigHt){ CalculateVariables(electrons->at(0),muons->at(0),EM); return true;}     
         }
         else if (leptonFlavor1 == "M" && leptonFlavor2 == "E") {
            if(CheckEMu() && selectedEvent.trigHt){ CalculateVariables(muons->at(0),electrons->at(0),ME); return true;}     
         }
      
   return false;
   }
   return false;
}




//ELECTRON
bool ntupleAnalyzerShrinked::testSelection(const tree::Electron& pa,selectionType selection ,bool leading){
   bool decision=false;
   if (selection==UNCUT){
      decision = pa.isPassConvVeto && pa.passImpactParameter && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && pa.isTightMVA;
   }
   if(selection==SEL || selection==DILEP || selection==ControlRegionZZ || selection==ControlRegionWZ || selection==SEL){
      decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (selectedEvent.deltaRll>0.1) && pa.isTightMVA;
   
   }
   if(selection==ONZ){
      decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (selectedEvent.deltaRll>0.1) && pa.isTightMVA;
      if (pa.isPassConvVeto && pa.isTightMVA){
         leading? decisionMapCutFlowFine[LEPTONIDPure_leading]=true : decisionMapCutFlowFine[LEPTONIDPure_trailing]=true;
         if (pa.passImpactParameter){
            leading? decisionMapCutFlowFine[LEPTONIDImpact_leading]=true : decisionMapCutFlowFine[LEPTONIDImpact_trailing]=true;
            if((fabs(pa.p.Eta())<2.4)){
               leading? decisionMapCutFlowFine[LEPTONIDEta_leading]=true : decisionMapCutFlowFine[LEPTONIDEta_trailing]=true;
               if(pa.miniIso<0.1){
                  leading? decisionMapCutFlowFine[LEPTONIDIso_leading]=true : decisionMapCutFlowFine[LEPTONIDIso_trailing]=true;
                  if(selectedEvent.deltaRll>0.1){
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
   if(selection==TRIGSEL || selection==TRIGDILEP || selection==TRIGONZ){
      decision =(*ht>200.) && pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>0.) : (pa.p.Pt()>0.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (selectedEvent.deltaRll>0.1) && pa.isTightMVA;
   }
   if(selection==TRIGSEL_ptcuts || selection==TRIGDILEP_ptcuts || selection==TRIGONZ_ptcuts){
      decision =(*ht>200.) && pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (selectedEvent.deltaRll>0.1) && pa.isTightMVA;
   }
   return decision;
}
//MUON
bool ntupleAnalyzerShrinked::testSelection(const tree::Muon& pa, selectionType selection_, bool leading){
   bool decision=false;
   if (selection_==UNCUT){
      decision = pa.passImpactParameter  && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (pa.isTight);
   }
   if (selection_==SEL || selection_==DILEP || selection_==ControlRegionZZ || selection_==ControlRegionWZ){
      decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.2) && (pa.isMedium) && (selectedEvent.deltaRll>0.1);
   }
   if (selection_ == ONZ){
      decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.2) && (pa.isMedium) && (selectedEvent.deltaRll>0.1);
      if((pa.isMedium)){
         leading? decisionMapCutFlowFine[LEPTONIDPure_leading]=true : decisionMapCutFlowFine[LEPTONIDPure_trailing]=true;
         if(pa.passImpactParameter){
            leading? decisionMapCutFlowFine[LEPTONIDImpact_leading]=true : decisionMapCutFlowFine[LEPTONIDImpact_trailing]=true;
            if((fabs(pa.p.Eta())<2.4) ){
               leading? decisionMapCutFlowFine[LEPTONIDEta_leading]=true : decisionMapCutFlowFine[LEPTONIDEta_trailing]=true;
               if(pa.miniIso<0.2){
                  leading? decisionMapCutFlowFine[LEPTONIDIso_leading]=true : decisionMapCutFlowFine[LEPTONIDIso_trailing]=true;
                  if((selectedEvent.deltaRll>0.1)){
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
   if(selection_==TRIGSEL || selection_==TRIGDILEP || selection_==TRIGONZ){
      decision = (*ht>200.) && pa.passImpactParameter && (leading ? (pa.p.Pt()>0.) : (pa.p.Pt()>0.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (selectedEvent.deltaRll>0.1);
   }
   if(selection_==TRIGSEL_ptcuts || selection_==TRIGDILEP_ptcuts || selection_==TRIGONZ_ptcuts){
      decision = (*ht>200.)&& pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (selectedEvent.deltaRll>0.1);
   }
   return decision;
}
//PHOTON
bool ntupleAnalyzerShrinked::testSelection(const selPhoton& pa, selectionType selection){
   bool decision=false;
   if(selection==SEL || selection==DILEP || selection==EGRegression || selection==ControlRegionZZ || selection==TRIGDILEP || selection==TRIGSEL || selection==TRIGONZ){
      decision = (pa.p.Pt()>25.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isLoose) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3); //study Delta R cut ;
   }
   if(selection==ONZ){
      decision = (pa.p.Pt()>25.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isLoose) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3); //study Delta R cut ;
      if((pa.isLoose)){
         if(selection==ONZ) decisionMapCutFlowFine[PHOTON1ID]=true;
         if(!(pa.hasPixelSeed)){
            if(selection==ONZ) decisionMapCutFlowFine[PHOTON1SEED]=true;
            if(fabs(pa.p.Eta()<1.4442)){
               if(selection==ONZ) decisionMapCutFlowFine[PHOTON1ETA]=true;
               if((pa.p.Pt()>25.)){
                  if(selection==ONZ) decisionMapCutFlowFine[PHOTON1PT]=true;
                  if((pa.deltaR1>0.3) && (pa.deltaR2>0.3)){
                     if(selection==ONZ) decisionMapCutFlowFine[PHOTON1DR]=true;
                  }
               }
            }
         }
      }
   }
   return decision;
}
//SELMUON
bool ntupleAnalyzerShrinked::testSelection(const selMuon& pa, selectionType selection_){
   bool decision=false;
   if (selection_==LooseLeptons){
      decision = pa.passImpactParameter && (pa.p.Pt()>20.) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.2) && (pa.isMedium);
   }
   return decision;
}
//SELELECTRON
bool ntupleAnalyzerShrinked::testSelection(const selElectron& pa, selectionType selection_){
   bool decision=false;
   if(selection_==LooseLeptons){
      decision = pa.isPassConvVeto && pa.passImpactParameter && (pa.p.Pt()>20.) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && pa.isTightMVA;
   }
   return decision;
}
//JET
bool ntupleAnalyzerShrinked::testSelection(const selJet& pa){
   bool decision=false;
   decision = (pa.p.Pt()>30.) && (fabs(pa.p.Eta())<2.4) && (pa.isLoose) && !(pa.hasElectronMatch) && !(pa.hasPhotonMatch) && !(pa.hasMuonMatch); //&& (pa.deltaR1>0.3) && (pa.deltaR2>0.3) //study Delta R cut ;
   return decision;
}

bool ntupleAnalyzerShrinked::matchLepton(const selElectron& pa){
   bool foundMatch=false;
   if((pa.deltaR1 < 0.05 ) || (pa.deltaR2 < 0.05 )){
      foundMatch=true;
   }else{
      foundMatch=false;
   }
   return foundMatch;
}
bool ntupleAnalyzerShrinked::matchLepton(const selMuon& pa){
   bool foundMatch=false;
   TLorentzVector dummyLepton(0.,0.,0.,0.);
   dummyLepton.SetPtEtaPhiM(pa.p.Pt(),pa.p.Eta(),pa.p.Phi(),0.51099895e-3);
   if((pa.deltaR1 < 0.05 ) || (pa.deltaR2 < 0.05 )){
      foundMatch=true;
   }else{
      foundMatch=false;
   }
   return foundMatch;
}


bool ntupleAnalyzerShrinked::matchRecoPhoton(const selPhoton& pho, const string collection){
   bool matched=false;
   if(collection=="J"){
      for(vector<tree::Particle>::iterator it = genJets->begin(); it != genJets->end();it++){
         if(pho.p.DeltaR(it->p)<0.1){
            matched=matched||true;
         }
      }
   }else{
      if(collection=="P"){
         for(vector<tree::GenParticle>::iterator it = genParticles->begin(); it != genParticles->end(); it++){
            if((abs(it->pdgId)==22)&&(it->statusID==1)){
               if(pho.p.DeltaR(it->p)<0.1) matched=matched||true;
            }
         }
      }else{
         if(collection=="E"){
            for(vector<tree::GenParticle>::iterator it = genParticles->begin(); it != genParticles->end(); it++){
               if((abs(it->pdgId)==11)&&(it->statusID==1)){
                  if(pho.p.DeltaR(it->p)<0.1) matched=matched||true;
               }
            }
         }
      }
   }
   return matched;
}



void ntupleAnalyzerShrinked::CalculateVariables(const tree::Lepton& lep1, const tree::Lepton& lep2 ,const particleType particle=DUMMYPARTICLE){
   
   //leptons
   if (particle==E){
      selectedEvent.isDiElectron=true;
   }else{
      selectedEvent.isDiElectron=false;
   }
   if (particle==M){
      selectedEvent.isDiMuon=true;
   }else{
      selectedEvent.isDiMuon=false;
   }
   if (particle==EM){
      selectedEvent.isElectronMuon=true;
   }else{
      selectedEvent.isElectronMuon=false;
   }
   if (particle==ME){
      selectedEvent.isMuonElectron=true;
   }else{
      selectedEvent.isMuonElectron=false;
   }


   if(!isSignal){
      selectedEvent.trigDiEle = *hlt_ele17_ele12_iso || *hlt_ele23_ele12_iso || *hlt_doubleEle33 || *hlt_doubleEle33_mw;
      selectedEvent.trigDiMu  = *hlt_mu17_mu8_iso || *hlt_mu17_tkMu8_iso || *hlt_mu17_mu8_iso_dz || *hlt_mu17_tkMu8_iso_dz
                  || *hlt_tkMu17_tkMu8_iso_dz || *hlt_mu27_tkMu8 || *hlt_mu30_tkMu11;
      selectedEvent.trigMuEle = *hlt_mu17_ele12_iso || *hlt_mu23_ele8_iso || *hlt_mu23_ele8_iso_dz || *hlt_mu23_ele12_iso
                  || *hlt_mu23_ele12_iso_dz || *hlt_mu8_ele17_iso || *hlt_mu8_ele23_iso
                  || *hlt_mu8_ele23_iso_dz || *hlt_mu12_ele23_iso || *hlt_mu12_ele23_iso_dz
                  || *hlt_mu30_ele30 || *hlt_mu33_ele33;
      selectedEvent.trigHt    = *hlt_ht200 || *hlt_ht250 || *hlt_ht300 || *hlt_ht350 || *hlt_ht400
                  || *hlt_ht475 || *hlt_ht600 || *hlt_ht650 || *hlt_ht800;
   }else{
      selectedEvent.trigDiEle = true;
      selectedEvent.trigDiMu = true;
      selectedEvent.trigMuEle = true;
      selectedEvent.trigHt = true;
   }

   //selectedEvent.trigMuEleMatch=selectedEvent.trigMuEle;
   //trigObjEleMatched=false
   if(!isSignal){
      if(selectedEvent.isDiElectron){
         if(selectedEvent.trigDiEle){
            bool matchedL1=false;
            bool matchedL2=false;
            for(vector<tree::Particle>::iterator it = trigObj_ele17_ele12->begin(); it != trigObj_ele17_ele12->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_ele23_ele12->begin(); it != trigObj_ele23_ele12->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_ele33_ele33->begin(); it != trigObj_ele33_ele33->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_ele33_ele33_mw->begin(); it != trigObj_ele33_ele33_mw->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            selectedEvent.trigDiEleMatch=(matchedL1 && matchedL2);
         }
      }
      
      if(selectedEvent.isDiMuon){
         if(selectedEvent.trigDiMu){
            bool matchedL1=false;
            bool matchedL2=false;
            for(vector<tree::Particle>::iterator it = trigObj_mu17_mu8->begin(); it != trigObj_mu17_mu8->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu17_mu8tk->begin(); it != trigObj_mu17_mu8tk->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu17_mu8_dz->begin(); it != trigObj_mu17_mu8_dz->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu17_mu8tk_dz->begin(); it != trigObj_mu17_mu8tk_dz->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu17tk_mu8tk_dz->begin(); it != trigObj_mu17tk_mu8tk_dz->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu27_mu8->begin(); it != trigObj_mu27_mu8->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu30_mu11->begin(); it != trigObj_mu30_mu11->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            selectedEvent.trigDiMuMatch=(matchedL1 && matchedL2);
         }
      }
      if(selectedEvent.isMuonElectron||selectedEvent.isElectronMuon){
         if(selectedEvent.trigMuEle){
            bool matchedL1=false;
            bool matchedL2=false;
            for(vector<tree::Particle>::iterator it = trigObj_mu17_ele12_eleLeg->begin(); it != trigObj_mu17_ele12_eleLeg->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu17_ele12_muLeg->begin(); it != trigObj_mu17_ele12_muLeg->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu23_ele8_eleLeg->begin(); it != trigObj_mu23_ele8_eleLeg->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu23_ele8_muLeg->begin(); it != trigObj_mu23_ele8_muLeg->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu23_ele8_dz->begin(); it != trigObj_mu23_ele8_dz->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu23_ele12_eleLeg->begin(); it != trigObj_mu23_ele12_eleLeg->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu23_ele12_muLeg->begin(); it != trigObj_mu23_ele12_muLeg->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu23_ele12_dz->begin(); it != trigObj_mu23_ele12_dz->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu8_ele17_eleLeg->begin(); it != trigObj_mu8_ele17_eleLeg->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu8_ele17_muLeg->begin(); it != trigObj_mu8_ele17_muLeg->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu8_ele23_eleLeg->begin(); it != trigObj_mu8_ele23_eleLeg->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu8_ele23_muLeg->begin(); it != trigObj_mu8_ele23_muLeg->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu8_ele23_dz->begin(); it != trigObj_mu8_ele23_dz->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu12_ele23_muLeg->begin(); it != trigObj_mu12_ele23_muLeg->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu12_ele23_eleLeg->begin(); it != trigObj_mu12_ele23_eleLeg->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu12_ele23_dz->begin(); it != trigObj_mu12_ele23_dz->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu30_ele30_eleLeg->begin(); it != trigObj_mu30_ele30_eleLeg->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu30_ele30_muLeg->begin(); it != trigObj_mu30_ele30_muLeg->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu33_ele33_eleLeg->begin(); it != trigObj_mu33_ele33_eleLeg->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            for(vector<tree::Particle>::iterator it = trigObj_mu33_ele33_muLeg->begin(); it != trigObj_mu33_ele33_muLeg->end();it++){
               if(lep1.p.DeltaR(it->p)<0.05){
                  matchedL1=matchedL1 || true;
               }
               if(lep2.p.DeltaR(it->p)<0.05){
                  matchedL2=matchedL2 || true;
               }
            }
            selectedEvent.trigMuEleMatch=(matchedL1 && matchedL2);
         }
      }
   }else{
      selectedEvent.trigDiEleMatch=true;
      selectedEvent.trigDiMuMatch=true;
      selectedEvent.trigMuEleMatch=true;
   }


   float leptonMass1 = (selectedEvent.isDiElectron||selectedEvent.isElectronMuon)?  0.51099895e-3 : 0.1056584;
   float leptonMass2 = (selectedEvent.isDiMuon||selectedEvent.isMuonElectron)?  0.1056584 : 0.51099895e-3;

   int charge1 = (int) lep1.charge;
   int charge2 = (int) lep2.charge;
   
   selectedEvent.l1.SetPtEtaPhiM(lep1.p.Pt(),lep1.p.Eta(),lep1.p.Phi(),leptonMass1);
   selectedEvent.l2.SetPtEtaPhiM(lep2.p.Pt(),lep2.p.Eta(),lep2.p.Phi(),leptonMass2);
   
   selectedEvent.pt1 = selectedEvent.l1.Pt();  
   selectedEvent.pt2 = selectedEvent.l2.Pt();  
   selectedEvent.phi1 = selectedEvent.l1.Phi();  
   selectedEvent.phi2 = selectedEvent.l2.Phi();  
   selectedEvent.eta1 = selectedEvent.l1.Eta();  
   selectedEvent.eta2 = selectedEvent.l2.Eta();

   selectedEvent.chargeProduct = charge1*charge2;
   selectedEvent.mll = (selectedEvent.l1+selectedEvent.l2).M();
   selectedEvent.miniIso1 = lep1.miniIso;
   selectedEvent.miniIso2 = lep2.miniIso;
   selectedEvent.deltaRll=selectedEvent.l1.DeltaR(selectedEvent.l2);
   //selectedEvent.ETmiss = mett.p.Pt();
   selectedEvent.ETmiss = met->p.Pt();
   selectedEvent.ETmiss_vec.SetPtEtaPhiM(met->p.Pt(),met->p.Eta(),met->p.Phi(),0.);;
   
   //double pa[3] = {lep1.M(),lep1.Px(),lep1.Py()};
   //double pb[3] = {lep2.M(),lep2.Px(),lep2.Py()};
   //double pmiss[3] = {0.,missing.Px(),missing.Py()};
   //fctMT2_.set_mn(0.);
   //fctMT2_.set_momenta(pa,pb,pmiss,*evtNo);
   //MT2_val = static_cast<float>(fctMT2_.get_mt2()); 
   
   selectedEvent.evtHasGenPhotonVeto = GenPhotonVeto(config_veto);
   
   for (vector<tree::Photon>::iterator it = photons->begin(); it != photons->end(); it++){
      auto g = *it;
      selPhoton gamma;
      gamma.setAll(g);
      gamma.vec.SetPtEtaPhiM(gamma.p.Pt(),gamma.p.Eta(),gamma.p.Phi(),0.);
      gamma.deltaR1=gamma.vec.DeltaR(selectedEvent.l1);
      gamma.deltaR2=gamma.vec.DeltaR(selectedEvent.l2);
      
      bool matchedToPhoton=matchRecoPhoton(gamma,"P");
      bool matchedToJet=matchRecoPhoton(gamma,"J");
      bool matchedToElectron=matchRecoPhoton(gamma,"E");
      
      gamma.matchedToPhoton=matchedToPhoton;
      gamma.matchedToJet=matchedToJet;
      gamma.matchedToElectron=matchedToElectron;
      //gamma.matchedToPhoton=false;
      //gamma.matchedToJet=false;
      //gamma.matchedToElectron=false;
      
      if (testSelection(gamma,SEL)){
         selectedEvent.selPhotons.push_back(gamma);
      }
   }
   for (vector<tree::Muon>::iterator it = muons->begin(); it != muons->end(); it++){
      auto g = *it;
      selMuon muon;
      muon.setAll(g);
      muon.vec.SetPtEtaPhiM(muon.p.Pt(),muon.p.Eta(),muon.p.Phi(),0.1056584);
      if (testSelection(muon,LooseLeptons)){
         muon.deltaR1 = muon.vec.DeltaR(selectedEvent.l1);
         muon.deltaR2 = muon.vec.DeltaR(selectedEvent.l2);
         muon.chargeInt = (int) g.charge;
         if(matchLepton(muon)){
            muon.matched=true;
            selectedEvent.matchedMuSize+=1;
            selectedEvent.matchedLeptonSize+=1;
         }
            selectedEvent.selMuons.push_back(muon);
      }
   }
   for (vector<tree::Electron>::iterator it = electrons->begin(); it != electrons->end(); it++){
      auto g = *it;
      selElectron ele;
      ele.setAll(g);
      ele.vec.SetPtEtaPhiM(ele.p.Pt(),ele.p.Eta(),ele.p.Phi(),0.51099895e-3);
      if (testSelection(ele,LooseLeptons)){
         ele.deltaR1 = ele.vec.DeltaR(selectedEvent.l1);
         ele.deltaR2 = ele.vec.DeltaR(selectedEvent.l2);
         ele.chargeInt = (int) g.charge;
         if(matchLepton(ele)){
            ele.matched=true;
            selectedEvent.matchedEleSize+=1;
            selectedEvent.matchedLeptonSize+=1;
         }
            selectedEvent.selElectrons.push_back(ele);
      }
   }
   
   selectedEvent.selMuonSize=(selectedEvent.selMuons.size());
   selectedEvent.selElectronSize=(selectedEvent.selElectrons.size());   
   selectedEvent.selLeptonSize=(selectedEvent.selMuonSize+selectedEvent.selElectronSize);
   
   for (vector<tree::Jet>::iterator it = jets->begin(); it != jets->end(); it++){
      auto j = *it;
      selJet jet;
      jet.setAll(j);
      jet.vec.SetPtEtaPhiM(jet.p.Pt(),jet.p.Eta(),jet.p.Phi(),0.);
      jet.deltaR1=jet.vec.DeltaR(selectedEvent.l1);
      jet.deltaR2=jet.vec.DeltaR(selectedEvent.l2);
      jet.bTag = (jet.bDiscriminator>0.8484); //MEDIUM WP include scale factors
      if (testSelection(jet)){
         selectedEvent.selJets.push_back(jet);
         selectedEvent.nselJets+=1;
         if(jet.bTag) selectedEvent.nselBJets+=1;
         selectedEvent.calcHt += jet.vec.Pt();
      }
   }   
}

void ntupleAnalyzerShrinked::ClearVariables(){

   selectedEvent.isDiElectron=false;
   selectedEvent.isDiMuon=false;
   selectedEvent.isMuonElectron=false;
   selectedEvent.isElectronMuon=false;
   float leptonMass1 = 0.;
   float leptonMass2 = 0.;
   
   selectedEvent.trigDiEle=false;
   selectedEvent.trigDiMu=false;
   selectedEvent.trigMuEle=false;
   selectedEvent.trigHt=false;
   
   int charge1 = 0;
   int charge2 = 0;
   
   selectedEvent.l1.SetPtEtaPhiM(0.,0.,0.,0.);
   selectedEvent.l2.SetPtEtaPhiM(0.,0.,0.,0.);

   selectedEvent.pt1 = 0.;  
   selectedEvent.pt2 = 0.;  
   selectedEvent.phi1 = 0.;  
   selectedEvent.phi2 = 0.;  
   selectedEvent.eta1 = 0.;  
   selectedEvent.eta2 = 0.;

   selectedEvent.miniIso1=1000.;
   selectedEvent.miniIso2=1000.;

   selectedEvent.chargeProduct = charge1*charge2;

   selectedEvent.mll = 0.;
   selectedEvent.deltaRll=0.;
   selectedEvent.ETmiss = 0.;
   selectedEvent.ETmiss_vec.SetPtEtaPhiM(0.,0.,0.,0.);
   
   //double pa[3] = {lep1.M(),lep1.Px(),lep1.Py()};
   //double pb[3] = {lep2.M(),lep2.Px(),lep2.Py()};
  //
   //double pmiss[3] = {0.,missing.Px(),missing.Py()};
   //fctMT2_.set_mn(0.);
   //fctMT2_.set_momenta(pa,pb,pmiss,*evtNo);
   //MT2_val = static_cast<float>(fctMT2_.get_mt2()); 
   
   selectedEvent.evtHasGenPhotonVeto = false;   
   
   selectedEvent.selPhotons.clear();
   selectedEvent.selMuons.clear();
   selectedEvent.selElectrons.clear();
   selectedEvent.selJets.clear();
   
   selectedEvent.matchedEleSize=0;
   selectedEvent.matchedMuSize=0;
   selectedEvent.matchedLeptonSize=0;
   
   selectedEvent.selMuonSize=0;
   selectedEvent.selElectronSize=0;   
   selectedEvent.selLeptonSize=0;
   
   selectedEvent.nselJets=0;
   selectedEvent.nselBJets=0;
   
   selectedEvent.calcHt=0.;
   
   }





bool ntupleAnalyzerShrinked::SelectEvent(selectionType selection){
   clearCutFlowMap();
   if (CheckParticles()){
      if (Cleaning()){//decide category in dilep events
      //if (cleaned){//decide category in dilep events
         //cout<<selectedEvent.mll<<endl;
         //cout<<selectedEvent.matchedEleSize<<selectedEvent.matchedMuSize<<selectedEvent.matchedLeptonSize<<endl;
         //cout<<selectedEvent.selMuonSize<<selectedEvent.selElectronSize<<selectedEvent.selPhotons.size()<<selectedEvent.selJets.size()<<endl;
         if(!isData){
            if(selectedEvent.isDiElectron){
               selectedEvent.totalWeight = totalWeight * GetScaleFactorAndErrorFinal(selectedEvent.pt1,selectedEvent.eta1,isSignal,true)*GetScaleFactorAndErrorFinal(selectedEvent.pt2,selectedEvent.eta2,isSignal,true);
            }else{
               if(selectedEvent.isDiMuon){
                  selectedEvent.totalWeight = totalWeight * GetScaleFactorAndErrorFinal(selectedEvent.pt1,selectedEvent.eta1,isSignal,false)*GetScaleFactorAndErrorFinal(selectedEvent.pt2,selectedEvent.eta2,isSignal,false);
               }
               else{
                  if(selectedEvent.isMuonElectron||selectedEvent.isElectronMuon){
                     selectedEvent.totalWeight = totalWeight * GetScaleFactorAndErrorFinal(muons->at(0).p.Pt(),muons->at(0).p.Eta(),isSignal,false)*GetScaleFactorAndErrorFinal(electrons->at(0).p.Pt(),electrons->at(0).p.Eta(),isSignal,true);
                  }
               }
            }
         }else{
            selectedEvent.totalWeight=totalWeight;
         }
         //fill cutflowInfo DiEle or DiMu?
         cutflowDiEle=selectedEvent.isDiElectron;
         cutflowDiMu=selectedEvent.isDiMuon;
         decisionMapCutFlowFine[DIMUON]=selectedEvent.isDiMuon;
         decisionMapCutFlowFine[DIELECTRON]=selectedEvent.isDiElectron;
         decisionMapCutFlowFine[EMUON]=selectedEvent.isMuonElectron||selectedEvent.isElectronMuon;
         
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
         

         
         if (selectedEvent.isDiElectron){
            //Fill cutflow triggered bool
            if(selection==ONZ) cutflowIsTriggered=true;
            if(selection==ONZ){
               if(isSignal){
                  decisionMapCutFlowFine[TRIGGERED]=((genNNToG==1)&&(genNNToZ==1)&&(genZToLL==2));
                  decisionMapCutFlowFine[TRIGGEREDMATCHED]=((genNNToG==1)&&(genNNToZ==1)&&(genZToLL==2));
               }else{
                  decisionMapCutFlowFine[TRIGGERED]=selectedEvent.trigDiEle;
                  decisionMapCutFlowFine[TRIGGEREDMATCHED]=selectedEvent.trigDiEleMatch;
               }
            }
            if(selection==ONZ) decisionMapCutFlowFine_weight[TRIGGERED]=totalWeight;
            if(selection==ONZ) decisionMapCutFlowFine_weight[TRIGGEREDMATCHED]=totalWeight;
            
            auto e1 = electrons->at(0);
            auto e2 = electrons->at(1);
            
            
            if(selectedEvent.trigDiEleMatch){
               if(!selectedEvent.evtHasGenPhotonVeto){
                  if(selection==ONZ) decisionMapCutFlowFine[GENVETO]=true;
                  if(selection==ONZ) decisionMapCutFlowFine_weight[GENVETO]=totalWeight;
         
                  if(testSelection(e1,selection,true) && testSelection(e2,selection,false) && (selectedEvent.chargeProduct < 0.) ){

                     //Fill cutflow 2 leptons bool
                     if(selection==ONZ) cutflow2Leptons=true;
                     if((selection==UNCUT)? true : (selectedEvent.mll>50.)){
                        //Fill cutflow mll>50 bool
                        if(selection==ONZ) cutflowMll50=true;
                        if(selection==ONZ) decisionMapCutFlowFine[M50]=true;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDPure_leading]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDPure_trailing]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDImpact_leading]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDImpact_trailing]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDEta_leading]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDEta_trailing]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDIso_leading]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDIso_trailing]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDDeltaR_leading]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDDeltaR_trailing]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONPT_leading]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONPT_trailing]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[M50]=selectedEvent.totalWeight;
                   
                        if((selection==ONZ)&&(photons->size()!=0)) decisionMapCutFlowFine[PHOTON1]=true;
                        if((selection==ONZ)&&(photons->size()!=0)) decisionMapCutFlowFine_weight[PHOTON1]=selectedEvent.totalWeight;
                        
                        if(!isData){
                              selectedEvent.totalWeight=selectedEvent.totalWeight*GetScaleFactorAndErrorPhotons(selectedEvent.selPhotons);
                        }
                        
                        if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1ID]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1SEED]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1ETA]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1PT]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1DR]=selectedEvent.totalWeight;
                        
                        if(selectedEvent.selPhotons.size()!=0){
                           if(selection==ONZ) cutflow1Photon=true;
                           if((selectedEvent.mll>81.) && (selectedEvent.mll<101.)){
                              if(selection==ONZ) cutflowOnZ=true;
                              if(selection==ONZ) decisionMapCutFlowFine[ZMASS]=true;
                              if(selection==ONZ) decisionMapCutFlowFine_weight[ZMASS]=selectedEvent.totalWeight;
                           }
                        }
                        return true;
                     }
                  }
               }
            }
         }
         if(selectedEvent.isDiMuon){

            if(selection==ONZ) decisionMapCutFlowFine_weight[genZLL]=totalWeight;
            
            
            if(selection==ONZ) cutflowIsTriggered=true;
            if(selection==ONZ){
               if(isSignal){
                  decisionMapCutFlowFine[TRIGGERED]=((genNNToG==1)&&(genNNToZ==1)&&(genZToLL==2));
                  decisionMapCutFlowFine[TRIGGEREDMATCHED]=((genNNToG==1)&&(genNNToZ==1)&&(genZToLL==2));
               }else{
                  decisionMapCutFlowFine[TRIGGERED]=selectedEvent.trigDiMu;
                  decisionMapCutFlowFine[TRIGGEREDMATCHED]=selectedEvent.trigDiMuMatch;
               }
            }
            if(selection==ONZ) decisionMapCutFlowFine_weight[TRIGGERED]=totalWeight;
            if(selection==ONZ) decisionMapCutFlowFine_weight[TRIGGEREDMATCHED]=totalWeight;

            auto m1 = muons->at(0); 
            auto m2 = muons->at(1);
            
            if(selectedEvent.trigDiMuMatch){
               if(!selectedEvent.evtHasGenPhotonVeto){
                  
                  if(selection==ONZ) decisionMapCutFlowFine[GENVETO]=true;
                  if(selection==ONZ) decisionMapCutFlowFine_weight[GENVETO]=totalWeight;


                  if(testSelection(m1,selection,true) && testSelection(m2,selection,false) && (selectedEvent.chargeProduct < 0.)){

                     if(selection==ONZ) cutflow2Leptons=true;
                     if((selection==UNCUT)? true : (selectedEvent.mll>50.)){
                        if(selection==ONZ) {cutflowMll50=true;}
                        if(selection==ONZ) decisionMapCutFlowFine[M50]=true;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDPure_leading]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDPure_trailing]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDImpact_leading]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDImpact_trailing]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDEta_leading]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDEta_trailing]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDIso_leading]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDIso_trailing]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDDeltaR_leading]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONIDDeltaR_trailing]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONPT_leading]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONPT_trailing]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[M50]=selectedEvent.totalWeight;

                        if(selection==ONZ) if(photons->size()!=0) decisionMapCutFlowFine[PHOTON1]=true;
                        if(selection==ONZ) if(photons->size()!=0) decisionMapCutFlowFine_weight[PHOTON1]=selectedEvent.totalWeight;
                        if(!isData){
                              selectedEvent.totalWeight=selectedEvent.totalWeight*GetScaleFactorAndErrorPhotons(selectedEvent.selPhotons);
                        }
                        
                        if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1ID]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1SEED]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1ETA]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1PT]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1DR]=selectedEvent.totalWeight;

                        if(selectedEvent.selPhotons.size()>0){
                           if(selection==ONZ) cutflow1Photon=true;
                           if(selectedEvent.mll<101. && selectedEvent.mll>81.){
                              if(selection==ONZ) cutflowOnZ=true;
                              if(selection==ONZ) decisionMapCutFlowFine[ZMASS]=true;
                              if(selection==ONZ) decisionMapCutFlowFine_weight[ZMASS]=selectedEvent.totalWeight;
                           }
                        }
                        return true;
                     }
                  }
               }
            }
         }
            
         if(selectedEvent.isElectronMuon || selectedEvent.isMuonElectron){
         
            if(selection==ONZ) cutflowIsTriggered=true;
            if(selection==ONZ){
               if(isSignal){
                  decisionMapCutFlowFine[TRIGGERED]=((genNNToG==1)&&(genNNToZ==1)&&(genZToLL==2));
                  decisionMapCutFlowFine[TRIGGEREDMATCHED]=((genNNToG==1)&&(genNNToZ==1)&&(genZToLL==2));
               }else{
                  decisionMapCutFlowFine[TRIGGERED]=selectedEvent.trigMuEle;
                  decisionMapCutFlowFine[TRIGGEREDMATCHED]=selectedEvent.trigMuEleMatch;
               }
            }
            if(selection==ONZ) decisionMapCutFlowFine_weight[TRIGGERED]=totalWeight;
            if(selection==ONZ) decisionMapCutFlowFine_weight[TRIGGEREDMATCHED]=totalWeight;

            auto e0 = electrons->at(0); 
            auto m0 = muons->at(0);
            
            if(selectedEvent.trigMuEleMatch){
               if(!selectedEvent.evtHasGenPhotonVeto){
                  
                  if(selection==ONZ) decisionMapCutFlowFine[GENVETO]=true;
                  if(selection==ONZ) decisionMapCutFlowFine_weight[GENVETO]=totalWeight;

                  if(testSelection(m0,selection,selectedEvent.isMuonElectron) && testSelection(e0,selection,selectedEvent.isElectronMuon) && (selectedEvent.chargeProduct < 0.)){
                     if(selection==ONZ) cutflow2Leptons=true;
                     if((selection==UNCUT)? true : (selectedEvent.mll>50.)){
                        if(selection==ONZ) cutflowMll50=true;
                        if(selection==ONZ) decisionMapCutFlowFine[M50]=true;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[M50]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONID_leading]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONID_trailing]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONPT_leading]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[LEPTONPT_trailing]=selectedEvent.totalWeight;


                        if(selection==ONZ) if(photons->size()!=0) decisionMapCutFlowFine[PHOTON1]=true;
                        if(selection==ONZ) if(photons->size()!=0) decisionMapCutFlowFine_weight[PHOTON1]=selectedEvent.totalWeight;
                        if(!isData){
                              selectedEvent.totalWeight=selectedEvent.totalWeight*GetScaleFactorAndErrorPhotons(selectedEvent.selPhotons);
                        }
                        
                        if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1ID]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1PT]=selectedEvent.totalWeight;
                        if(selection==ONZ) decisionMapCutFlowFine_weight[PHOTON1DR]=selectedEvent.totalWeight;
   ;
                           
                        if(selectedEvent.selPhotons.size()>0){
                           if(selection==ONZ) cutflow1Photon=true;
                           if(selectedEvent.mll<101. && selectedEvent.mll>81.){
                              if(selection==ONZ) cutflowOnZ=true;
                              if(selection==ONZ) decisionMapCutFlowFine[ZMASS]=true;
                              if(selection==ONZ) decisionMapCutFlowFine_weight[ZMASS]=selectedEvent.totalWeight;
                           }
                        }
                        return true;
                     }
                  }
               }
            }
         }
         return false;
      }else{return false;}
   }else{return false;}
  return false;
}

bool ntupleAnalyzerShrinked::SelectEventTriggerStudies(selectionType selection){
   if (CheckParticles()){//only Fill Histograms if nMu>=2 or nEle>=2, nGamma>=0
      if (CleaningTriggerStudies()){//decide in events with DiEle and DiMu where to contribute
         if(!isData){
            if(selectedEvent.isDiElectron){
               selectedEvent.totalWeight = totalWeight * GetScaleFactorAndErrorFinal(selectedEvent.pt1,selectedEvent.eta1,isSignal,true)*GetScaleFactorAndErrorFinal(selectedEvent.pt2,selectedEvent.eta2,isSignal,true);
            }else{
               if(selectedEvent.isDiMuon){
                  selectedEvent.totalWeight = totalWeight * GetScaleFactorAndErrorFinal(selectedEvent.pt1,selectedEvent.eta1,isSignal,false)*GetScaleFactorAndErrorFinal(selectedEvent.pt2,selectedEvent.eta2,isSignal,false);
               }
               else{
                  if(selectedEvent.isMuonElectron||selectedEvent.isElectronMuon){
                     selectedEvent.totalWeight = totalWeight * GetScaleFactorAndErrorFinal(muons->at(0).p.Pt(),muons->at(0).p.Eta(),isSignal,false)*GetScaleFactorAndErrorFinal(electrons->at(0).p.Pt(),electrons->at(0).p.Eta(),isSignal,true);
                  }
               }
            }
         }else{
            selectedEvent.totalWeight=totalWeight;
         }
         if (selectedEvent.isDiElectron){
            auto e1 = electrons->at(0);
            auto e2 = electrons->at(1);
            if(!selectedEvent.evtHasGenPhotonVeto){
               if(testSelection(e1,selection,true) && testSelection(e2,selection,false) && (selectedEvent.chargeProduct < 0.) && ((selection==UNCUT)? true : (selectedEvent.mll>50.)) ){
                  return true;
            }
         }
      }
      else{
         if (selectedEvent.isDiMuon){ 
         auto m1 = muons->at(0); 
         auto m2 = muons->at(1);
         if(!selectedEvent.evtHasGenPhotonVeto){
            if(testSelection(m1,selection,true) && testSelection(m2,selection,false) && (selectedEvent.chargeProduct < 0.) && ((selection==UNCUT)? true : (selectedEvent.mll>50.))){
               return true;
            }
         }
      }else{
            if (selectedEvent.isElectronMuon || selectedEvent.isMuonElectron){
               auto e0 = electrons->at(0); 
               auto m0 = muons->at(0);
               if(!selectedEvent.evtHasGenPhotonVeto){
                  if(testSelection(e0,selection,selectedEvent.isElectronMuon) && testSelection(m0,selection,selectedEvent.isMuonElectron) && (selectedEvent.chargeProduct < 0.) && ((selection==UNCUT)? true : (selectedEvent.mll>50.))){
                     return true;
                  }
               }
            }
         }
         }
      }else{return false;}
   }else{return false;}
  return false;
}



void ntupleAnalyzerShrinked::FillCutFlowHistograms(){
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
void ntupleAnalyzerShrinked::SetCutFlowHistogramsStatus(){
   cutflowIsTriggered=false;
   cutflow2Leptons=false;
   cutflowMll50=false;
   cutflow1Photon=false;
   cutflowOnZ=false;
}


void ntupleAnalyzerShrinked::FillCutFlowHistograms_Fine(){
   auto m=decisionMapCutFlowFine;
   auto mW=decisionMapCutFlowFine_weight;
   if(m[DIELECTRON] && !(m[DIMUON]||m[EMUON])){
            
      if (m[TRIGGERED]){
         c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("triggered",mW[TRIGGERED]);
         if(m[TRIGGEREDMATCHED]){
            c1Maps["cutFlow_Fine_onZEE"].at(CUTFLOW_fine).Fill("triggeredMatched",mW[TRIGGEREDMATCHED]);
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
      }
   }
   if(m[DIMUON] && !(m[DIELECTRON]||m[EMUON])){
      if (m[TRIGGERED]){
         c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("triggered",mW[TRIGGERED]);
         if(m[TRIGGEREDMATCHED]){
            c1Maps["cutFlow_Fine_onZMM"].at(CUTFLOW_fine).Fill("triggeredMatched",mW[TRIGGEREDMATCHED]);
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
   }
   if(m[EMUON] && !(m[DIELECTRON]||m[DIMUON])){
      if (m[TRIGGERED]){
         c1Maps["cutFlow_Fine_onZEM"].at(CUTFLOW_fine).Fill("triggered",mW[TRIGGERED]);
         if(m[TRIGGEREDMATCHED]){
            c1Maps["cutFlow_Fine_onZEM"].at(CUTFLOW_fine).Fill("triggeredMatched",mW[TRIGGEREDMATCHED]);
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
            h.second.Write(histoNames[h.first].c_str(), TObject::kWriteDelete);
         }
         file.cd();
      }
   }
}





void ntupleAnalyzerShrinked::Terminate()
{
   auto outputName = "output"+config_outputfolder+"/"+getOutputFilename(inputName);

   
   TFile file(outputName.c_str(), "RECREATE");
   save2File(h1Maps, file);
   save2File(eff1Maps, file);
   save2File(h2Maps, file);
   save2File(c1Maps, file);
   save2File(s1Maps, file);
   save2File(cr1Maps, file);
   if(! isTotalSignal){
      cutFlow.Write("hCutFlow");
   }
   file.Close();
   cout << "Created " << outputName << " in " << (time(NULL) - startTime)/60 << " min" << endl;
}






