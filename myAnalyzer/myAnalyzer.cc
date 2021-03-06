#include "myAnalyzer.h"

myAnalyzer::myAnalyzer() :

        selPhotons(fReader,"selPhotons"),
        selJets(fReader,"selJets"),
        selElectrons(fReader,"selElectrons"),
        selMuons(fReader,"selMuons"),

        chargeProduct(fReader,"chargeProduct"),
        mll(fReader,"mll"),
        miniIso1(fReader,"miniIso1"),
        miniIso2(fReader,"miniIso2"),
        deltaRll(fReader,"deltaRll"),

        pt1(fReader,"pt1"),
        pt2(fReader,"pt2"),
        phi1(fReader,"phi1"),
        phi2(fReader,"phi2"),
        eta1(fReader,"eta1"),
        eta2(fReader,"eta2"),
        l1(fReader,"l1"),
        l2(fReader,"l2"),

        evtHasGenPhotonVeto(fReader,"genPhotonVeto"),

        ETmiss(fReader,"ETmiss"),
        ETmiss_vec(fReader,"ETmiss_vec"),

        selMuonSize(fReader,"selMuonSize"),
        selElectronSize(fReader,"selElectronSize"),
        selLeptonSize(fReader,"selLeptonSize"),
        selPhotonSize(fReader,"selPhotonSize"),
        selJetSize(fReader,"selJetSize"),
        selBJetSize(fReader,"selBJetSize"),
        matchedEleSize(fReader,"matchedEleSize"),
        matchedLeptonSize(fReader,"matchedLeptonSize"),
        matchedMuSize(fReader,"matchedMuSize"),

        isDiElectron(fReader,"isDiElectron"),
        isDiMuon(fReader,"isDiMuon"),
        isMuonElectron(fReader,"isMuonElectron"),
        isElectronMuon(fReader,"isElectronMuon"),

        calcHt(fReader,"calcHt"),
        Mt2(fReader,"MT2"),

        lepSF_weight(fReader,"lepSF_weight"),
        lepSF_weightUp(fReader,"lepSF_weightUp"),
        lepSF_weightDown(fReader,"lepSF_weightDown"),
        photonSF_weight(fReader,"photonSF_weight"),
        photonSF_weightUp(fReader,"photonSF_weightUp"),
        photonSF_weightDown(fReader,"photonSF_weightDown"),
        topPt_weight(fReader,"topPt_weight"),
        isr_weight(fReader,"isr_weight"),
        isr_weightUp(fReader,"isr_weightUp"),
        isr_weightDown(fReader,"isr_weightDown"),
        ewk_weight(fReader,"isr_weight"),
        ewk_weightUp(fReader,"ewk_weightUp"),
        ewk_weightDown(fReader,"ewk_weightDown"),


        countNegCharge(fReader,"countNegCharge"),
        countPosCharge(fReader,"countPosCharge"),

        negElectrons(fReader,"negElectrons"),
        posElectrons(fReader,"posElectrons"),
        negMuons(fReader,"negMuons"),
        posMuons(fReader,"posMuons"),

        trigHt(fReader,"htTriggered"),
        trigMET(fReader,"metTriggered"),
        trigDiEle(fReader,"eeTriggered"),
        trigDiEleMatch(fReader,"eeTriggeredMatch"),
        trigDiMu(fReader,"mmTriggered"),
        trigDiMuMatch(fReader,"mmTriggeredMatch"),
        trigMuEle(fReader,"emTriggered"),
        trigMuEleMatch(fReader,"emTriggeredMatch"),

        photons(fReader, "photons"),
        jets(fReader, "jets"),
        electrons(fReader, "electrons"),
        muons(fReader, "muons"),
        genJets(fReader, "genJets"),
        genParticles(fReader, "genParticles"),
        intermediateGenParticles(fReader, "intermediateGenParticles"),
        met(fReader, "met"),
        metRaw(fReader, "metRaw"),
        met_JESu(fReader, "metJESu"),
        met_gen(fReader, "metgen"),
        met_JESd(fReader, "metJESd"),
        met_JERu(fReader, "metJERu"),
        met_JERd(fReader, "metJERd"),
        nGoodVertices(fReader, "nGoodVertices"),
        nTracksPV(fReader, "nTracksPV"),
        pu_weight(fReader, "pu_weight"),
        pu_weightUp(fReader, "pu_weightUp"),
        pu_weightDown(fReader, "pu_weightDown"),
        mc_weight(fReader, "mc_weight"),
        pdf_weights(fReader, "pdf_weights"),
        genHt(fReader, "genHt"),
        ht(fReader, "ht"),
        nTruePV(fReader, "nTruePV"),

        nISR(fReader, "nISR"),
        EWKinoPairPt(fReader, "EWKinoPairPt"),
        leptonPairPt(fReader, "leptonPairPt"),
        topPt1(fReader, "topPt1"),
        topPt2(fReader, "topPt2"),

        runNo(fReader, "runNo"),
        lumNo(fReader, "lumNo"),
        evtNo(fReader, "evtNo"),

        signal_m1(fReader,"signal_m1"),
        signal_m2(fReader,"signal_m2"),
        nBinos(fReader,"nBinos"),
        startTime(time(NULL))//,
        //rand()
{
}


void myAnalyzer::Init(TTree *tree)
{
        fReader.SetTree(tree);
        inputName = fReader.GetTree()->GetCurrentFile()->GetName();

        //cout<<inputName<<endl;

        isTotalSignal = (inputName.find("SMS-TChiNG_BF") != string::npos) || (inputName.find("GMSB_GravitinoLSP") != string::npos) || (inputName.find("SMS-T6ttZg_myTuple")!=string::npos) || (inputName.find("SMS-T5bbbbZg_myTuple")!=string::npos) || (inputName.find("GGM")!=string::npos);
        isData = inputName.find("Run201") != string::npos;
        isSignal = (inputName.find("SMS") != string::npos) || (inputName.find("GGM") != string::npos) || (inputName.find("GMSB") != string::npos);

        isTTInclusive = (inputName.find("TTTo2L2Nu") != string::npos)|| (inputName.find("TT_my") != string::npos);
        isZZ2L = (inputName.find("ZZTo2L") != string::npos);
        isZZ4L = (inputName.find("ZZTo4L") != string::npos);

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
        config_selectionsToProcessMap[ControlRegionDY] = propertyTree.get<bool>("selectionsToProcess.control");
        config_selectionsToProcessMap[ControlRegionTT] = propertyTree.get<bool>("selectionsToProcess.control");
        config_selectionsToProcessMap[ControlRegionZZ] = propertyTree.get<bool>("selectionsToProcess.control");
        config_selectionsToProcessMap[ControlRegionWZ] = propertyTree.get<bool>("selectionsToProcess.control");
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
        config_doTrigger = propertyTree.get<bool>("generalSettings.htstudies");
        config_doTriggerPure = propertyTree.get<bool>("generalSettings.htstudiespure");
        config_doTriggerMET = propertyTree.get<bool>("generalSettings.metstudies");

        config_doPdfSplit = propertyTree.get<bool>("generalSettings.pdfsplit");
        config_pdfStart = propertyTree.get<int>("generalSettings.pdfstart");
        config_pdfStep = propertyTree.get<int>("generalSettings.pdfstep");

        config_eeWeight = propertyTree.get<float>("generalSettings.eeweight");
        config_emWeight = propertyTree.get<float>("generalSettings.emweight");
        config_mmWeight = propertyTree.get<float>("generalSettings.mmweight");
        config_eeEff = propertyTree.get<float>("generalSettings.eeeff");
        config_emEff = propertyTree.get<float>("generalSettings.emeff");
        config_mmEff = propertyTree.get<float>("generalSettings.mmeff");

        if(isSignal) {
                config_doPdfSplit=false;
        }

        if(isTotalSignal) {
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
                config_doTrigger = false;

                config_doPdfSplit = false;
        }





        //map<unsigned short,unsigned short> massPointsForMaps;
        vector<SignalPoint> massPointsForMaps;
        massPointsForMaps.clear();

        if(!isTotalSignal) {
                //cutFlow = *((TH1F*)fReader.GetTree()->GetCurrentFile()->Get("TreeWriter/hCutFlow"));
                cutFlow = *((TH1F*)fReader.GetTree()->GetCurrentFile()->Get("hCutFlow"));
                //massPointsForMaps[0]=[0];
                sp_.first=0;
                sp_.second=0;
                massPointsForMaps.push_back(sp_);

                //if(!config_doTrigger)weightHisto=*((TH1F*)fReader.GetTree()->GetCurrentFile()->Get("weightHisto"));
                weightHisto=*((TH1F*)fReader.GetTree()->GetCurrentFile()->Get("weightHisto"));

        }else{
                TList* list = fReader.GetTree()->GetCurrentFile()->GetListOfKeys();
                if (!list) { printf("<E> No keys found in file\n"); exit(1); }
                TIter next(list);
                TKey* key;
                TObject* obj;

                while ( (key = ((TKey*)next())) ) {
                        obj = key->ReadObj();
                        if (    (strcmp(obj->IsA()->GetName(),"TProfile")!=0)
                                && (!obj->InheritsFrom("TH1"))
                                && (!obj->InheritsFrom("TH1"))
                                ) {
                                printf("<W> Object %s is not 1D or 2D histogram : "
                                       "will not be converted\n",obj->GetName());
                        }else{

                                string text = obj->GetName();
                                vector<string> results;
                                results.clear();

                                boost::split(results, text, boost::is_any_of("_"), boost::token_compress_on);

                                unsigned short toFillM1=0;
                                unsigned short toFillM2=0;

                                if((inputName.find("GGM") != string::npos)) {
                                        if((text.find("weightHisto") != string::npos)) {
                                                continue;
                                        }
                                }

                                cout<<text<<endl;
                                cout<<results.size()<<endl;

                                string tempWeightHistName="";

                                if(results.size()>2) {
                                        toFillM1=atoi(results.at(1).c_str());
                                        toFillM2=atoi(results.at(2).c_str());
                                        if((inputName.find("GGM") != string::npos)) {
                                                tempWeightHistName=tempWeightHistName+"_"+results.at(1).replace(0,2,"").c_str()+"_"+results.at(2).replace(0,2,"").c_str();
                                                cout<<results.at(1)<<"----"<<results.at(2)<<endl;
                                                toFillM1=atoi(results.at(1).replace(0,0,"").c_str());
                                                toFillM2=atoi(results.at(2).replace(0,0,"").c_str());
                                        }else{
                                                tempWeightHistName=tempWeightHistName+"_"+results.at(1).c_str()+"_"+results.at(2).c_str();
                                                toFillM1=atoi(results.at(1).c_str());
                                                toFillM2=atoi(results.at(2).c_str());
                                        }
                                }else{
                                        if(results.size()>1) {
                                                toFillM1=atoi(results.at(1).c_str());
                                                tempWeightHistName=tempWeightHistName+"_"+results.at(1).c_str();
                                        }
                                }

                                cout<<tempWeightHistName<<endl;
                                cout<<toFillM1<<toFillM2<<endl;
                                sp_.first=toFillM1;
                                sp_.second=toFillM2;
                                massPointsForMaps.push_back(sp_);

                                //if(!config_doTrigger)weightHistoMap["weightHisto"+tempWeightHistName]=*((TH1F*)fReader.GetTree()->GetCurrentFile()->Get(("weightHisto"+tempWeightHistName).c_str()));
                                if(!(config_doTrigger||config_doTriggerPure||config_doTriggerMET)) weightHistoMap["weightHisto"+tempWeightHistName]=*((TH1F*)fReader.GetTree()->GetCurrentFile()->Get(("weightHisto"+tempWeightHistName).c_str()));

                        }
                }



        }




        fReader.GetEntries(true);
        nEntries=fReader.GetTree()->GetEntries();




        setHistoNames();
        setFolderNames();




        if(config_doPdfSplit) {
                if(isSignal) {
                        nWeights=((inputName.find("GGM") != string::npos)||(inputName.find("GMSB") != string::npos)) ? 0 : 9;
                }else{
                        //nWeights=(isData||(inputName.find("HadronicDecays") != string::npos)||(inputName.find("GGM") != string::npos)||(inputName.find("GMSB") != string::npos))? 0 : 110;
                        nWeights=(isData||(inputName.find("HadronicDecays") != string::npos)||(inputName.find("GGM") != string::npos)||(inputName.find("GMSB") != string::npos)) ? 0 : config_pdfStep;
                }
        }else{
                if(isSignal) {
                        nWeights=((inputName.find("GGM") != string::npos)||(inputName.find("GMSB") != string::npos)) ? 0 : 9;
                        config_pdfStart=0;
                        config_pdfStep=((inputName.find("GGM") != string::npos)||(inputName.find("GMSB") != string::npos)) ? 0 : 10;
                }else{

                        //nWeights=(isData||(inputName.find("HadronicDecays") != string::npos)||(inputName.find("GGM") != string::npos)||(inputName.find("GMSB") != string::npos)) ? 0 : 110;
                        nWeights=(isData||(inputName.find("HadronicDecays") != string::npos)||(inputName.find("GGM") != string::npos)||(inputName.find("GMSB") != string::npos)) ? 0 : 10;
                        //nWeights=(isData||(inputName.find("HadronicDecays") != string::npos)||(inputName.find("GGM") != string::npos)||(inputName.find("GMSB") != string::npos))? 0 : 10;

                        config_pdfStart=0;
                        //config_pdfStep=110;
                        //config_pdfStep=(isData||(inputName.find("HadronicDecays") != string::npos)||(inputName.find("GGM") != string::npos)||(inputName.find("GMSB") != string::npos))? 0:110;
                        config_pdfStep=(isData||(inputName.find("HadronicDecays") != string::npos)||(inputName.find("GGM") != string::npos)||(inputName.find("GMSB") != string::npos)) ? 0 : 0;
                        //config_pdfStep=(isData||(inputName.find("HadronicDecays") != string::npos)||(inputName.find("GGM") != string::npos)||(inputName.find("GMSB") != string::npos))? 0:10;
                        //nWeights=(isData||(inputName.find("HadronicDecays") != string::npos)||(inputName.find("GGM") != string::npos)||(inputName.find("GMSB") != string::npos))? 0 : 0;
                        //nWeights=(isData||(inputName.find("HadronicDecays") != string::npos)||(inputName.find("GGM") != string::npos)||(inputName.find("GMSB") != string::npos))? 0 : 9;
                }
                //config_pdfStart=0;
                //config_pdfStep=110;
                //config_pdfStep=0;
        }



        //if(!config_doTrigger){
        if(!(config_doTrigger||config_doTriggerMET||config_doTriggerPure)) {
                if (config_docutflow) InitCutFlowHistos();
                if (config_docutflowfine) InitCutFlowHistos_Fine();
        }
        if(!isTotalSignal) {
                //if(config_doTrigger){
                if(config_doTrigger||config_doTriggerMET||config_doTriggerPure) {
                        InitTriggerStudiesHistos();
                }else{
                        InitAllHistos();
                        for(vector<SignalPoint>::iterator it = massPointsForMaps.begin(); it != massPointsForMaps.end(); it++) {
                                InitSignalScanHistos(*it);
                        }
                }
        }else{
                float counter=0.;
                float maxCounter=massPointsForMaps.size();
                float initProgress= counter/maxCounter;
                for(vector<SignalPoint>::iterator it = massPointsForMaps.begin(); it != massPointsForMaps.end(); it++) {
                        initProgress= counter/maxCounter;
                        InitSignalScanHistos(*it);
                        cout<<counter<<"/"<<maxCounter<<" ["<<initProgress*100.<<"%] in "<<(time(NULL) - startTime)/60 << " min" <<endl;
                        counter+=1.;
                }
        }
        //cout<<"done"<<endl;
}


void myAnalyzer::SlaveBegin(TTree *tree)
{
}



Bool_t myAnalyzer::Process(Long64_t entry){

        double tempPercentage = (double) entry/ (double)nEntries;
        if(!(abs(config_eventpercentage-100.)<0.1)) {
                if (tempPercentage>config_eventpercentage/100.) {
                        return kTRUE;
                }
        }



        fReader.SetLocalEntry(entry);

        if(((inputName.find("ZZTo2L") != string::npos))||((inputName.find("ZGToLL") != string::npos))) {
                if(pdf_weights->size()<110) {
                        return kTRUE;
                }
        }



        if(isSignal && !isTotalSignal) {
                sp_.first=0;
                sp_.second=0;
        }else{
                sp_.first=*signal_m1;
                sp_.second=*signal_m2;
        }

        totalWeight = *mc_weight * *pu_weight;

        zzKFactor=GetZZKFactor();



        puWeights[normalPU]=*pu_weight;
        puWeights[upPU]=*pu_weightUp;
        puWeights[downPU]=*pu_weightDown;
        puWeights[noPU]=1.;
        puWeights[noPUUp]=1.;
        puWeights[noPUDown]=1.;

        lepSfWeights[normalLEPSF]=*lepSF_weight;
        lepSfWeights[upLEPSF]=*lepSF_weightUp;
        lepSfWeights[downLEPSF]=*lepSF_weightDown;

        photonSfWeights[normalPHOTONSF]=*photonSF_weight;
        photonSfWeights[upPHOTONSF]=*photonSF_weightUp;
        photonSfWeights[downPHOTONSF]=*photonSF_weightDown;

        ewkWeights[normalEWK]=*ewk_weight;
        ewkWeights[upEWK]=*ewk_weightUp;
        ewkWeights[downEWK]=*ewk_weightDown;

        isrWeights[normalISR]=*isr_weight;
        isrWeights[upISR]=*isr_weightUp;
        isrWeights[downISR]=*isr_weightDown;


        long progress = tempPercentage*100.;
        //if(entry%100000==0){
        //if(entry%10000==0){
        if(entry%100==0) {

                std::cout<<"[";
                for(long i=0; i<100; i++)
                        if(i<progress)
                                std::cout<<'=';
                        else if(i==progress)
                                std::cout<<'>';
                        else
                                std::cout<<' ';
                //std::cout<<"] "<<progress<<" %"<<" "<<getOutputFilename(inputName)<<'\r';
                std::cout<<"] "<<tempPercentage*100.<<" %"<<" "<<getOutputFilename(inputName)<<'\r';
                std::cout.flush();

        }


        string cutFlowName = "hCutFlow";

        if(isTotalSignal) {
                if(inputName.find("TChi")!=string::npos) {
                        cutFlowName+="TChiNG";
                }else{
                        if(inputName.find("T5")!=string::npos) {
                                cutFlowName+="T5bbbbZg";
                        }else{
                                if(inputName.find("GGM")!=string::npos) {
                                        cutFlowName+="GGM";
                                }else{
                                        if(inputName.find("GMSB")!=string::npos) {
                                                cutFlowName+="GMSB";
                                        }else{
                                                if(inputName.find("T6")!=string::npos) {
                                                        cutFlowName+="T6ttZg";
                                                }
                                        }
                                }
                        }
                }
                if(inputName.find("GGM")!=string::npos) {
                        cutFlowName += "_M1"+to_string(*signal_m1);
                }else{
                        cutFlowName += "_"+to_string(*signal_m1);
                }
                if (*signal_m2) {
                        if(inputName.find("GGM")!=string::npos) {
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
                if(cutFlow.GetEntries()) {
                        nGen = cutFlow.GetBinContent(2);
                }
        }


        if (!isTotalSignal) {
                //if(!config_doTrigger){
                if(!(config_doTrigger||config_doTriggerMET||config_doTriggerPure)) {
                        SetCutFlowHistogramsStatus();
                }

        }

        if (!isTotalSignal) {
                //if(!config_doTrigger){
                if(!(config_doTrigger||config_doTriggerMET||config_doTriggerPure)) {
                        FillHistograms();
                }
        }
        //if(!config_doTrigger){
        if(!(config_doTrigger||config_doTriggerMET||config_doTriggerPure)) {
                if (config_docutflow) FillCutFlowHistograms();
                if (config_docutflowfine) FillCutFlowHistograms_Fine();
                if (!isTotalSignal) {
                        clearCutFlowMap();
                }
        }



        //if((inputName.find("JetHT")!= string::npos)||(!isData && !isTotalSignal)){
        if((inputName.find("MET")!= string::npos)||(inputName.find("JetHT")!= string::npos)||(!isData && !isTotalSignal)) {
                //if(config_doTrigger){
                if(config_doTrigger||config_doTriggerMET||config_doTriggerPure) {
                        FillTriggerStudies();
                }
        }

        if (config_dosignalscan) {
                //if(!config_doTrigger){
                if(!(config_doTrigger||config_doTriggerMET||config_doTriggerPure)) {
                        FillSignalHistograms();
                }
        }

        if (!isTotalSignal) {
                //if(!config_doTrigger){
                if(!(config_doTrigger||config_doTriggerMET||config_doTriggerPure)) {
                        clearCutFlowMap();
                }
        }

        return kTRUE;
}

//void myAnalyzer::FillHistograms2D(){
//if(SelectEvent(SEL)){
//if (selectedEvent.selPhotons.size()!=0){
//Filler2D(h2Maps["sel"],true);
//if (selectedEvent.isDiElectron){
//Filler2D(selectedEvent,h2Maps["selEE"],true);
//}else{
//if (selectedEvent.isDiMuon) Filler2D(selectedEvent,h2Maps["selMM"],true);
//}
//}
//}

//if(SelectEvent(ONZ)){
//if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>81. && selectedEvent.mll<101.)){
//Filler2D(selectedEvent,h2Maps["onZ"],true);
//if (selectedEvent.isDiElectron){
//Filler2D(selectedEvent,h2Maps["onZEE"],true);
//}else{
//if (selectedEvent.isDiMuon) Filler2D(selectedEvent,h2Maps["onZMM"],true);
//}
//}
//}
//if(SelectEvent(EXO)){
//if ((selectedEvent.selPhotons.size()!=0)&&(selectedEvent.mll>50. && selectedEvent.mll<130.)){
//Filler2D(selectedEvent,h2Maps["exo"],true);
//if (selectedEvent.isDiElectron){
//Filler2D(selectedEvent,h2Maps["exoEE"],true);
//}else{
//if (selectedEvent.isDiMuon) Filler2D(selectedEvent,h2Maps["exoMM"],true);
//}
//}
//}
//}

void myAnalyzer::FillSignalHistograms(){
        auto signalPoint=sp_;
        float TChiNG_weight=1.;
        bool useTChiNG=true;
        bool isTChiNG=false;
        // original: gg=4/16,zz=1/16,hh=1/16,gz=4/16,gh=4/16,zh=2/16
        if(inputName.find("TChiNG") != string::npos) {
                isTChiNG=true;
                //cout<<"here"<<endl;
                //cout<<(*intermediateGenParticles)[0].daughters.size()<<(*intermediateGenParticles)[1].daughters.size()<<endl;
                //cout<<(*intermediateGenParticles)[0].daughters[0].pdgId<<(*intermediateGenParticles)[0].daughters[1].pdgId<<endl;
                //cout<<(*intermediateGenParticles)[1].daughters[0].pdgId<<(*intermediateGenParticles)[1].daughters[1].pdgId<<endl;
                //GG case -> weight=1. 0.25to0.25
                bool oneOfThem=false;
                if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2) {
                        if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 22) {
                                //cout<<"buh"<<endl;
                                if (fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 22) {
                                        //cout<<"sdfdfsd"<<endl;
                                        useTChiNG=true;
                                        TChiNG_weight=1.;
                                        oneOfThem=oneOfThem||true;

                                }
                        }
                }
                //ZZ case ->weight=4. 1/16 to 0.25
                if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2) {
                        if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 23) {
                                if (fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 23) {
                                        useTChiNG=true;
                                        TChiNG_weight=4.;
                                        oneOfThem=oneOfThem||true;

                                }
                        }
                }
                //cout<<"test"<<endl;
                //ZG case -> weight = 2 0.25to0.5
                if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2) {
                        //cout<<"test 2"<<endl;
                        if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 23) {
                                //cout<<"test3"<<endl;
                                if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 22) {
                                        useTChiNG=true;
                                        TChiNG_weight=2.;
                                        oneOfThem=oneOfThem||true;
                                }
                        }
                }
                if(oneOfThem==false) {
                        useTChiNG=false;
                        TChiNG_weight=1.;
                }



        }
        bool usePoint=(isTChiNG&&useTChiNG) || (!isTChiNG);
        if(usePoint) {
                if(SelectEvent(ONZ)) {




                        if(*selLeptonSize==*matchedLeptonSize && *selLeptonSize==2 && (*selPhotonSize!=0)) {
                                //if((*selPhotonSize!=0)&&(*mll>81. && *mll<101.)&&(*ETmiss>=150.)){
                                //if((*selPhotonSize!=0)&&(*mll>81. && *mll<101.)){
                                //if((*selPhotonSize!=0)&&(*mll>81. && *mll<101.)&&(*Mt2>=100.)){
                                //if((*selPhotonSize!=0)&&(*mll>81. && *mll<101.)&&((*Mt2>=100.)||(*Mt2<100. && selPhotons->at(0).p.Pt()>80.))){
                                if((*selPhotonSize!=0)&&(*mll>81. && *mll<101.)) {
                                        float tempPhotonLeadingPt = selPhotons->at(0).p.Pt();

                                        //if(usePoint){
                                        //string signalPoint;
                                        //if(inputName.find("GGM")!=string::npos){
                                        //signalPoint = getSignalPointName(10,*signal_m1,*signal_m2);
                                        //}else{
                                        //if(inputName.find("GMSB")!=string::npos){
                                        //signalPoint = getSignalPointName(11,*signal_m1,*signal_m2);
                                        //}else{
                                        //signalPoint = getSignalPointName(*nBinos,*signal_m1,*signal_m2);
                                        //}
                                        //}
                                        //SignalPoint sp_;
                                        //if(!(s1Maps.count(signalPoint)>0)){
                                        //if(!(s1Maps.count(sp_)>0)){
                                        //InitSignalScanHistos(signalPoint);
                                        //InitSignalScanHistos(sp_);
                                        //}
                                        //auto signalPoint=sp_;
                                        if(*isDiMuon || *isDiElectron) {

                                                //if(*ETmiss>=150.){
                                                //if((*ETmiss>=150. && *Mt2>=100.)||(*ETmiss>=100. && selPhotons->at(0).p.Pt()>=80. && (*Mt2<100||*ETmiss<150.))){
                                                if((*ETmiss>=150. && *Mt2>=100.)) {
                                                        //cout<<isTChiNG<<useTChiNG<<TChiNG_weight<<endl;
                                                        if(isTChiNG && useTChiNG) {
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(nom),TChiNG_weight,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        }else{
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        }
                                                        //if(*Mt2>100.) FillerSignal(s1Maps.at(signalPoint).at(sigMt2).at(LL).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        //FillerSignal(s2Maps.at(signalPoint).at(sig).at(LL).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                }
                                                //if(met_JESu->p.Pt()>=150.){
                                                //if((met_JESu->p.Pt()>=150. && *Mt2>=100.)||(met_JESu->p.Pt()>=100. && selPhotons->at(0).p.Pt()>=80. && (*Mt2<100||met_JESu->p.Pt()<150.))){
                                                if((met_JESu->p.Pt()>=150. && *Mt2>=100.)) {
                                                        if(isTChiNG && useTChiNG) {
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(JESu),TChiNG_weight,9999,JESUP,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        }else{
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(JESu),1.,9999,JESUP,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        }
                                                        //if(*Mt2>100.) FillerSignal(s1Maps.at(signalPoint).at(sigMt2).at(LL).at(JESu),1.,9999,JESUP,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        //FillerSignal(s2Maps.at(signalPoint).at(sig).at(LL).at(JESu),1.,9999,JESUP,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                }
                                                //if(met_gen->p.Pt()>=150.){
                                                //if((met_gen->p.Pt()>=150. && *Mt2>=100.)||(met_gen->p.Pt()>=100. && selPhotons->at(0).p.Pt()>=80. && (*Mt2<100||met_gen->p.Pt()<150.))){
                                                if((met_gen->p.Pt()>=150. && *Mt2>=100.)) {
                                                        if(isTChiNG && useTChiNG) {
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(GENMET),TChiNG_weight,9999,METGEN,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        }else{
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(GENMET),1.,9999,METGEN,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        }
                                                        //if(*Mt2>100.) FillerSignal(s1Maps.at(signalPoint).at(sigMt2).at(LL).at(GENMET),1.,9999,METGEN,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        //FillerSignal(s2Maps.at(signalPoint).at(sig).at(LL).at(GENMET),1.,9999,METGEN,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                }
                                                //if(met_JESd->p.Pt()>=150.){
                                                //if((met_JESd->p.Pt()>=150. && *Mt2>=100.)||(met_JESd->p.Pt()>=100. && selPhotons->at(0).p.Pt()>=80. && (*Mt2<100||met_JESd->p.Pt()<150.))){
                                                if((met_JESd->p.Pt()>=150. && *Mt2>=100.)) {
                                                        if(isTChiNG && useTChiNG) {
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(JESd),TChiNG_weight,9999,JESDOWN,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        }else{
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(JESd),1.,9999,JESDOWN,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        }
                                                        //FillerSignal(s2Maps.at(signalPoint).at(sig).at(LL).at(JESd),1.,9999,JESDOWN,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        //if(*Mt2>100.) FillerSignal(s1Maps.at(signalPoint).at(sigMt2).at(LL).at(JESd),1.,9999,JESDOWN,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                }
                                                //if(met_JERu->p.Pt()>=150.){
                                                //if((met_JERu->p.Pt()>=150. && *Mt2>=100.)||(met_JERu->p.Pt()>=100. && selPhotons->at(0).p.Pt()>=80. && (*Mt2<100||met_JERu->p.Pt()<150.))){
                                                if((met_JERu->p.Pt()>=150. && *Mt2>=100.)) {
                                                        if(isTChiNG && useTChiNG) {
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(JERu),TChiNG_weight,9999,JERUP,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        }else{
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(JERu),1.,9999,JERUP,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        }
                                                        //if(*Mt2>100.) FillerSignal(s1Maps.at(signalPoint).at(sigMt2).at(LL).at(JERu),1.,9999,JERUP,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        //FillerSignal(s2Maps.at(signalPoint).at(sig).at(LL).at(JERu),1.,9999,JERUP,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                }
                                                //if(met_JERd->p.Pt()>=150.){
                                                //if((met_JERd->p.Pt()>=150. && *Mt2>=100.)||(met_JERd->p.Pt()>=100. && selPhotons->at(0).p.Pt()>=80. && (*Mt2<100||met_JERd->p.Pt()<150.))){
                                                if((met_JERd->p.Pt()>=150. && *Mt2>=100.)) {
                                                        if(isTChiNG && useTChiNG) {
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(JERd),TChiNG_weight,9999,JERDOWN,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        }else{
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(JERd),1.,9999,JERDOWN,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        }
                                                        //if(*Mt2>100.) FillerSignal(s1Maps.at(signalPoint).at(sigMt2).at(LL).at(JERd),1.,9999,JERDOWN,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        //FillerSignal(s2Maps.at(signalPoint).at(sig).at(LL).at(JERd),1.,9999,JERDOWN,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                }

                                                //if(*ETmiss>=150.){
                                                //if((*ETmiss>=150. && *Mt2>=100.)||(*ETmiss>=100. && selPhotons->at(0).p.Pt()>=80. && (*Mt2<100||*ETmiss<150.))){
                                                if((*ETmiss>=150. && *Mt2>=100.)) {
                                                        if(isTChiNG && useTChiNG) {
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(PUUP),TChiNG_weight,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        }else{
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(PUUP),1.,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        }
                                                        //if(*Mt2>100.) FillerSignal(s1Maps.at(signalPoint).at(sigMt2).at(LL).at(PUUP),1.,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        //FillerSignal(s2Maps.at(signalPoint).at(sig).at(LL).at(PUUP),1.,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        if(isTChiNG && useTChiNG) {
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(PUDOWN),TChiNG_weight,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(LEPSFUP),TChiNG_weight,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        }else{
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(PUDOWN),1.,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(LEPSFUP),1.,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        }
                                                        //if(*Mt2>100.) FillerSignal(s1Maps.at(signalPoint).at(sigMt2).at(LL).at(PUDOWN),1.,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        //if(*Mt2>100.) FillerSignal(s1Maps.at(signalPoint).at(sigMt2).at(LL).at(LEPSFUP),1.,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        //FillerSignal(s2Maps.at(signalPoint).at(sig).at(LL).at(PUDOWN),1.,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        //FillerSignal(s2Maps.at(signalPoint).at(sig).at(LL).at(LEPSFUP),1.,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);

                                                        //FillerSignal(s1Maps.at(signalPoint).at(sigMt2).at(LL).at(NOPU),1.,9999,normal,noPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        //FillerSignal(s1Maps.at(signalPoint).at(sigMt2).at(LL).at(NOPU),1.,9999,normal,noPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        //FillerSignal(s2Maps.at(signalPoint).at(sig).at(LL).at(NOPU),1.,9999,normal,noPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);

                                                        if (*nTruePV>=20) {
                                                                if(isTChiNG && useTChiNG) {
                                                                        FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(NOPUUP),TChiNG_weight,9999,normal,noPUUp,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                }else{
                                                                        FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(NOPUUP),1.,9999,normal,noPUUp,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                }
                                                                //if(*Mt2>100.) FillerSignal(s1Maps.at(signalPoint).at(sigMt2).at(LL).at(NOPUUP),1.,9999,normal,noPUUp,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(LL).at(NOPUUP),1.,9999,normal,noPUUp,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        }else{
                                                                if(isTChiNG && useTChiNG) {
                                                                        FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(NOPUDOWN),TChiNG_weight,9999,normal,noPUDown,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                }else{
                                                                        FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(NOPUDOWN),1.,9999,normal,noPUDown,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                }
                                                                //if(*Mt2>100.) FillerSignal(s1Maps.at(signalPoint).at(sigMt2).at(LL).at(NOPUDOWN),1.,9999,normal,noPUDown,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(LL).at(NOPUDOWN),1.,9999,normal,noPUDown,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        }
                                                        if(isTChiNG && useTChiNG) {
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(LEPSFDOWN),TChiNG_weight,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(PHOTONSFUP),TChiNG_weight,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                        }else{
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(LEPSFDOWN),1.,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(PHOTONSFUP),1.,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                        }
                                                        //if(*Mt2>100.) FillerSignal(s1Maps.at(signalPoint).at(sigMt2).at(LL).at(LEPSFDOWN),1.,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        //if(*Mt2>100.) FillerSignal(s1Maps.at(signalPoint).at(sigMt2).at(LL).at(PHOTONSFUP),1.,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                        //FillerSignal(s2Maps.at(signalPoint).at(sig).at(LL).at(LEPSFDOWN),1.,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        //FillerSignal(s2Maps.at(signalPoint).at(sig).at(LL).at(PHOTONSFUP),1.,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);

                                                        if(isTChiNG && useTChiNG) {
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(PHOTONSFDOWN),TChiNG_weight,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(ISRUP),TChiNG_weight,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                        }else{
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(PHOTONSFDOWN),1.,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(ISRUP),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                        }
                                                        //if(*Mt2>100.) FillerSignal(s1Maps.at(signalPoint).at(sigMt2).at(LL).at(PHOTONSFDOWN),1.,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                        //if(*Mt2>100.) FillerSignal(s1Maps.at(signalPoint).at(sigMt2).at(LL).at(ISRUP),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                        //FillerSignal(s2Maps.at(signalPoint).at(sig).at(LL).at(PHOTONSFDOWN),1.,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                        //FillerSignal(s2Maps.at(signalPoint).at(sig).at(LL).at(ISRUP),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);

                                                        if(isTChiNG && useTChiNG) {
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(ISRDOWN),TChiNG_weight,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(EWKUP),TChiNG_weight,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                        }else{
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(ISRDOWN),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(EWKUP),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                        }
                                                        //if(*Mt2>100.) FillerSignal(s1Maps.at(signalPoint).at(sigMt2).at(LL).at(ISRDOWN),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                        //if(*Mt2>100.) FillerSignal(s1Maps.at(signalPoint).at(sigMt2).at(LL).at(EWKUP),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                        //FillerSignal(s2Maps.at(signalPoint).at(sig).at(LL).at(ISRDOWN),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                        //FillerSignal(s2Maps.at(signalPoint).at(sig).at(LL).at(EWKUP),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);

                                                        if(isTChiNG && useTChiNG) {
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(EWKDOWN),TChiNG_weight,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                        }else{
                                                                FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(EWKDOWN),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                        }
                                                        //if(*Mt2>100.) FillerSignal(s1Maps.at(signalPoint).at(sigMt2).at(LL).at(EWKDOWN),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                        //FillerSignal(s2Maps.at(signalPoint).at(sig).at(LL).at(EWKDOWN),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);

                                                        //for(int i=0; i<nWeights; i++){
                                                        for(int i=config_pdfStart; i<config_pdfStart+config_pdfStep; i++) {
                                                                if(isTChiNG && useTChiNG) {
                                                                        FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(PDFNAMES[i]),TChiNG_weight,i,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                }else{
                                                                        FillerSignal(s1Maps.at(signalPoint).at(sig).at(LL).at(PDFNAMES[i]),1.,i,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                }
                                                                //if(*Mt2>100.) FillerSignal(s1Maps.at(signalPoint).at(sigMt2).at(LL).at(PDFNAMES[i]),1.,i,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(LL).at(PDFNAMES[i]),1.,i,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        }


                                                        if(tempPhotonLeadingPt>=80.) {
                                                                if(isTChiNG && useTChiNG) {
                                                                        FillerSignal(s1Maps.at(signalPoint).at(sig80).at(LL).at(nom),TChiNG_weight,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                }else{
                                                                        FillerSignal(s1Maps.at(signalPoint).at(sig80).at(LL).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                }
                                                                //if(*Mt2>100.) FillerSignal(s1Maps.at(signalPoint).at(sig80Mt2).at(LL).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                //FillerSignal(s2Maps.at(signalPoint).at(sig80).at(LL).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        }else{
                                                                if(isTChiNG && useTChiNG) {
                                                                        FillerSignal(s1Maps.at(signalPoint).at(sig080).at(LL).at(nom),TChiNG_weight,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                }else{
                                                                        FillerSignal(s1Maps.at(signalPoint).at(sig080).at(LL).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                }
                                                                //if(*Mt2>100.) FillerSignal(s1Maps.at(signalPoint).at(sig080Mt2).at(LL).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                //FillerSignal(s2Maps.at(signalPoint).at(sig080).at(LL).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        }

                                                }
                                                //if (*isDiElectron){
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(EE).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(EE).at(JESu),1.,9999,JESUP,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(EE).at(JESd),1.,9999,JESDOWN,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(EE).at(JERu),1.,9999,JERUP,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(EE).at(GENMET),1.,9999,METGEN,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(EE).at(PUUP),1.,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(EE).at(PUDOWN),1.,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(EE).at(NOPU),1.,9999,normal,noPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);

                                                //if (*nTruePV>=20){
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(EE).at(NOPUUP),1.,9999,normal,noPUUp,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //}else{
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(EE).at(NOPUDOWN),1.,9999,normal,noPUDown,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //}
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(EE).at(LEPSFUP),1.,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(EE).at(LEPSFUP),1.,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(EE).at(PHOTONSFUP),1.,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(EE).at(PHOTONSFDOWN),1.,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(EE).at(ISRUP),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(EE).at(ISRDOWN),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(EE).at(EWKUP),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(EE).at(EWKDOWN),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(EE).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(EE).at(JESu),1.,9999,JESUP,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(EE).at(JESd),1.,9999,JESDOWN,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(EE).at(JERu),1.,9999,JERUP,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(EE).at(JERd),1.,9999,JERDOWN,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(EE).at(GENMET),1.,9999,METGEN,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(EE).at(PUUP),1.,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(EE).at(PUDOWN),1.,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(EE).at(NOPU),1.,9999,normal,noPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);

                                                //if (*nTruePV>=20){
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(EE).at(NOPUUP),1.,9999,normal,noPUUp,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //}else{
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(EE).at(NOPUDOWN),1.,9999,normal,noPUDown,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //}
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(EE).at(LEPSFUP),1.,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(EE).at(LEPSFUP),1.,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(EE).at(PHOTONSFUP),1.,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(EE).at(PHOTONSFDOWN),1.,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(EE).at(ISRUP),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(EE).at(ISRDOWN),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(EE).at(EWKUP),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(EE).at(EWKDOWN),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                ///////for(int i=0; i<nWeights; i++){
                                                //for(int i=config_pdfStart; i<config_pdfStart+config_pdfStep; i++){
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(EE).at(PDFNAMES[i]),1.,i,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(EE).at(PDFNAMES[i]),1.,i,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //}

                                                //if(tempPhotonLeadingPt>=80.){
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig80).at(EE).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig80).at(EE).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //}else{
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig080).at(EE).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig080).at(EE).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //}
                                                //}
                                                //if (*isDiMuon){
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(MM).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(MM).at(JESu),1.,9999,JESUP,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(MM).at(JESd),1.,9999,JESDOWN,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(MM).at(JERu),1.,9999,JERUP,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(MM).at(JERd),1.,9999,JERDOWN,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(MM).at(GENMET),1.,9999,METGEN,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(MM).at(PUUP),1.,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(MM).at(PUDOWN),1.,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(MM).at(NOPU),1.,9999,normal,noPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);

                                                //if (*nTruePV>=20){
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(MM).at(NOPUUP),1.,9999,normal,noPUUp,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //}else{
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(MM).at(NOPUDOWN),1.,9999,normal,noPUDown,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //}
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(MM).at(LEPSFUP),1.,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(MM).at(LEPSFUP),1.,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(MM).at(PHOTONSFUP),1.,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(MM).at(PHOTONSFDOWN),1.,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(MM).at(ISRUP),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(MM).at(ISRDOWN),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(MM).at(EWKUP),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(MM).at(EWKDOWN),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(MM).at(JESu),1.,9999,JESUP,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(MM).at(JESd),1.,9999,JESDOWN,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(MM).at(JERu),1.,9999,JERUP,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(MM).at(JERd),1.,9999,JERDOWN,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(MM).at(GENMET),1.,9999,METGEN,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(MM).at(PUUP),1.,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(MM).at(PUDOWN),1.,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(MM).at(NOPU),1.,9999,normal,noPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);

                                                //if (*nTruePV>=20){
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(MM).at(NOPUUP),1.,9999,normal,noPUUp,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //}else{
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(MM).at(NOPUDOWN),1.,9999,normal,noPUDown,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //}
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(MM).at(LEPSFUP),1.,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(MM).at(LEPSFUP),1.,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(MM).at(PHOTONSFUP),1.,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(MM).at(PHOTONSFDOWN),1.,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(MM).at(ISRUP),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(MM).at(ISRDOWN),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(MM).at(EWKUP),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(MM).at(EWKDOWN),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                ///////for(int i=0; i<nWeights; i++){
                                                //for(int i=config_pdfStart; i<config_pdfStart+config_pdfStep; i++){
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig).at(MM).at(PDFNAMES[i]),1.,i,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig).at(MM).at(PDFNAMES[i]),1.,i,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //}

                                                //if(tempPhotonLeadingPt>=80.){
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig80).at(MM).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig80).at(MM).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //}else{
                                                //FillerSignal(s1Maps.at(signalPoint).at(sig080).at(MM).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerSignal(s2Maps.at(signalPoint).at(sig080).at(MM).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //}
                                                //}
                                        }
                                        //if(*isElectronMuon || *isMuonElectron){
                                        //FillerSignal(s1Maps.at(signalPoint).at(sig).at(EM).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                        //FillerSignal(s2Maps.at(signalPoint).at(sig).at(EM).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                        //}

                                        if(config_dosignalscanSplit) {

                                                if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2) {
                                                        if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 22) {
                                                                if (fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 22) {
                                                                        if(*isDiMuon || *isDiElectron) {
                                                                                FillerSignal(s1Maps.at(signalPoint).at(sig_gg).at(LL).at(nom),4.0);
                                                                        }
                                                                        if (*isDiElectron) {
                                                                                FillerSignal(s1Maps.at(signalPoint).at(sig_gg).at(EE).at(nom),4.0);
                                                                        }
                                                                        if (*isDiMuon) {
                                                                                FillerSignal(s1Maps.at(signalPoint).at(sig_gg).at(MM).at(nom),4.0);
                                                                        }
                                                                        if(*isElectronMuon || *isMuonElectron) {
                                                                                FillerSignal(s1Maps.at(signalPoint).at(sig_gg).at(EM).at(nom),4.0);
                                                                        }
                                                                }
                                                        }
                                                }
                                                if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2) {
                                                        if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 23) {
                                                                if (fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 23) {
                                                                        if(*isDiMuon || *isDiElectron) {
                                                                                FillerSignal(s1Maps.at(signalPoint).at(sig_zz).at(LL).at(nom),16.0);
                                                                        }
                                                                        if (*isDiElectron) {
                                                                                FillerSignal(s1Maps.at(signalPoint).at(sig_zz).at(EE).at(nom),16.0);
                                                                        }
                                                                        if (*isDiMuon) {
                                                                                FillerSignal(s1Maps.at(signalPoint).at(sig_zz).at(MM).at(nom),16.0);
                                                                        }
                                                                        if(*isElectronMuon || *isMuonElectron) {
                                                                                FillerSignal(s1Maps.at(signalPoint).at(sig_zz).at(EM).at(nom),16.0);
                                                                        }
                                                                }
                                                        }
                                                }
                                                if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2) {
                                                        if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 23) {
                                                                if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 22) {
                                                                        if(*isDiMuon || *isDiElectron) {
                                                                                FillerSignal(s1Maps.at(signalPoint).at(sig_gz).at(LL).at(nom),4.0);
                                                                        }
                                                                        if (*isDiElectron) {
                                                                                FillerSignal(s1Maps.at(signalPoint).at(sig_gz).at(EE).at(nom),4.0);
                                                                        }
                                                                        if (*isDiMuon) {
                                                                                FillerSignal(s1Maps.at(signalPoint).at(sig_gz).at(MM).at(nom),4.0);
                                                                        }
                                                                        if(*isElectronMuon || *isMuonElectron) {
                                                                                FillerSignal(s1Maps.at(signalPoint).at(sig_gz).at(EM).at(nom),4.0);

                                                                        }
                                                                }
                                                        }
                                                }
                                        }
                                }
                                if ((*selPhotonSize!=0)&&(*mll>81. && *mll<101.)) {//&&(selectedEvent.ETmiss<100.)){ //ONZ+MET<100 ~ CR DY/Z(+gamma)
                                        if(*ETmiss<100.) {
                                                if(*isDiMuon || *isDiElectron) {
                                                        FillerSignal(s1Maps.at(signalPoint).at(controlregionDY).at(LL).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                }
                                        }
                                }
                                //if ((*selPhotonSize!=0)&&(*mll>81. && *mll<101.)&&(*ETmiss>=100.)&&(*ETmiss<150.)){ //ONZ+MET<150>100 ~ VR
                                //if ((*selPhotonSize!=0)&&(*mll>81. && *mll<101.)&&(*ETmiss>=100.)&&((*ETmiss<150.)||(*ETmiss>=150. && *Mt2<100.))){ //ONZ+MET<150>100 ~ VR
                                if ((*selPhotonSize!=0)&&(*mll>81. && *mll<101.)&&(*ETmiss>=100.)&&((*ETmiss<150.)||(*Mt2<100.))) { //ONZ+MET<150>100 ~ VR
                                        //if(selPhotons->at(0).p.Pt()<80.){
                                        if(*isDiMuon || *isDiElectron) {
                                                FillerSignal(s1Maps.at(signalPoint).at(validationregion).at(LL).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                        }
                                        //}
                                }
                                if ((*selPhotonSize!=0)) { //Different Flavor + 1Photon ~ CR TT(+gamma)
                                        if (*isElectronMuon || *isMuonElectron) {
                                                FillerSignal(s1Maps.at(signalPoint).at(controlregionTT).at(EM).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                        }
                                }
                        }
                }
                if(SelectEvent(ControlRegionZZ)) {
                        if((*selLeptonSize==4)&&(*matchedLeptonSize==2)) {
                                if((*selMuonSize==4)||(*selElectronSize==4)||(*selMuonSize==2 && *selElectronSize==2)) {

                                        float ZMass = 91.1876;

                                        if(*countNegCharge==2 && *countPosCharge==2) {

                                                selLepton lep1,lep2,lep3,lep4;

                                                if(*selMuonSize==4) {

                                                        if((negMuons->size()==2)&&(posMuons->size()==2)) {
                                                                selMuon neg1 = negMuons->at(0);
                                                                selMuon neg2 = negMuons->at(1);
                                                                selMuon pos1 = posMuons->at(0);
                                                                selMuon pos2 = posMuons->at(1);

                                                                //there are 2 combinations n1p1/n2p2 - n1p2/n2p1
                                                                float n1p1 = abs((neg1.vec+pos1.vec).M()-ZMass);
                                                                float n1p2 = abs((neg1.vec+pos2.vec).M()-ZMass);
                                                                float n2p1 = abs((neg2.vec+pos1.vec).M()-ZMass);
                                                                float n2p2 = abs((neg2.vec+pos2.vec).M()-ZMass);

                                                                if((n1p1<n1p2)&&(n1p1<n2p1)&&(n1p1<n2p2)) {
                                                                        if(neg1.p.Pt()>pos1.p.Pt()) {
                                                                                lep1.setAll(neg1);
                                                                                lep2.setAll(pos1);
                                                                        }else{
                                                                                lep1.setAll(pos1);
                                                                                lep2.setAll(neg1);
                                                                        }
                                                                        if(neg2.p.Pt()>pos2.p.Pt()) {
                                                                                lep3.setAll(neg2);
                                                                                lep4.setAll(pos2);
                                                                        }else{
                                                                                lep3.setAll(pos2);
                                                                                lep4.setAll(neg2);
                                                                        }
                                                                }else{
                                                                        if((n1p2<n1p1)&&(n1p2<n2p1)&&(n1p2<n2p2)) {
                                                                                if(neg1.p.Pt()>pos2.p.Pt()) {
                                                                                        lep1.setAll(neg1);
                                                                                        lep2.setAll(pos2);
                                                                                }else{
                                                                                        lep1.setAll(pos2);
                                                                                        lep2.setAll(neg1);
                                                                                }
                                                                                if(neg2.p.Pt()>pos1.p.Pt()) {
                                                                                        lep3.setAll(neg2);
                                                                                        lep4.setAll(pos1);
                                                                                }else{
                                                                                        lep3.setAll(pos1);
                                                                                        lep4.setAll(neg2);
                                                                                }
                                                                        }else{
                                                                                if((n2p1<n1p1)&&(n2p1<n1p2)&&(n2p1<n2p2)) {
                                                                                        if(neg2.p.Pt()>pos1.p.Pt()) {
                                                                                                lep1.setAll(neg2);
                                                                                                lep2.setAll(pos1);
                                                                                        }else{
                                                                                                lep1.setAll(pos1);
                                                                                                lep2.setAll(neg2);
                                                                                        }
                                                                                        if(neg1.p.Pt()>pos2.p.Pt()) {
                                                                                                lep3.setAll(neg1);
                                                                                                lep4.setAll(pos2);
                                                                                        }else{
                                                                                                lep3.setAll(pos2);
                                                                                                lep4.setAll(neg1);
                                                                                        }
                                                                                }else{
                                                                                        if(neg2.p.Pt()>pos2.p.Pt()) {
                                                                                                lep1.setAll(neg2);
                                                                                                lep2.setAll(pos2);
                                                                                        }else{
                                                                                                lep1.setAll(pos2);
                                                                                                lep2.setAll(neg2);
                                                                                        }
                                                                                        if(neg1.p.Pt()>pos1.p.Pt()) {
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
                                                        if(mll1<106. && mll1>76. && mll2<130. && mll2>50.) {
                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][nom],false,lep1,lep2,lep3,lep4,true);
                                                                FillerSignal(s1Maps.at(signalPoint).at(controlregionZZ).at(LL).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        }
                                                }

                                                if(*selElectronSize==4) {

                                                        if((negElectrons->size()==2)&&(posElectrons->size()==2)) {
                                                                selElectron neg1 = negElectrons->at(0);
                                                                selElectron neg2 = negElectrons->at(1);
                                                                selElectron pos1 = posElectrons->at(0);
                                                                selElectron pos2 = posElectrons->at(1);

                                                                //there are 2 combinations n1p1/n2p2 - n1p2/n2p1
                                                                float n1p1 = abs((neg1.vec+pos1.vec).M()-ZMass);
                                                                float n1p2 = abs((neg1.vec+pos2.vec).M()-ZMass);
                                                                float n2p1 = abs((neg2.vec+pos1.vec).M()-ZMass);
                                                                float n2p2 = abs((neg2.vec+pos2.vec).M()-ZMass);


                                                                if((n1p1<n1p2)&&(n1p1<n2p1)&&(n1p1<n2p2)) {
                                                                        if(neg1.p.Pt()>pos1.p.Pt()) {
                                                                                lep1.setAll(neg1);
                                                                                lep2.setAll(pos1);
                                                                        }else{
                                                                                lep1.setAll(pos1);
                                                                                lep2.setAll(neg1);
                                                                        }
                                                                        if(neg2.p.Pt()>pos2.p.Pt()) {
                                                                                lep3.setAll(neg2);
                                                                                lep4.setAll(pos2);
                                                                        }else{
                                                                                lep3.setAll(pos2);
                                                                                lep4.setAll(neg2);
                                                                        }
                                                                }else{
                                                                        if((n1p2<n1p1)&&(n1p2<n2p1)&&(n1p2<n2p2)) {
                                                                                if(neg1.p.Pt()>pos2.p.Pt()) {
                                                                                        lep1.setAll(neg1);
                                                                                        lep2.setAll(pos2);
                                                                                }else{
                                                                                        lep1.setAll(pos2);
                                                                                        lep2.setAll(neg1);
                                                                                }
                                                                                if(neg2.p.Pt()>pos1.p.Pt()) {
                                                                                        lep3.setAll(neg2);
                                                                                        lep4.setAll(pos1);
                                                                                }else{
                                                                                        lep3.setAll(pos1);
                                                                                        lep4.setAll(neg2);
                                                                                }
                                                                        }else{
                                                                                if((n2p1<n1p1)&&(n2p1<n1p2)&&(n2p1<n2p2)) {
                                                                                        if(neg2.p.Pt()>pos1.p.Pt()) {
                                                                                                lep1.setAll(neg2);
                                                                                                lep2.setAll(pos1);
                                                                                        }else{
                                                                                                lep1.setAll(pos1);
                                                                                                lep2.setAll(neg2);
                                                                                        }
                                                                                        if(neg1.p.Pt()>pos2.p.Pt()) {
                                                                                                lep3.setAll(neg1);
                                                                                                lep4.setAll(pos2);
                                                                                        }else{
                                                                                                lep3.setAll(pos2);
                                                                                                lep4.setAll(neg1);
                                                                                        }
                                                                                }else{
                                                                                        if(neg2.p.Pt()>pos2.p.Pt()) {
                                                                                                lep1.setAll(neg2);
                                                                                                lep2.setAll(pos2);
                                                                                        }else{
                                                                                                lep1.setAll(pos2);
                                                                                                lep2.setAll(neg2);
                                                                                        }
                                                                                        if(neg1.p.Pt()>pos1.p.Pt()) {
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
                                                        if(mll1<106. && mll1>76. && mll2<130. && mll2>50.) {
                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][nom],false,lep1,lep2,lep3,lep4,true);
                                                                FillerSignal(s1Maps.at(signalPoint).at(controlregionZZ).at(LL).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        }

                                                }

                                                if(*selMuonSize==2 && *selElectronSize==2) {

                                                        if((negMuons->size()==1)&&(posMuons->size()==1)&&(negElectrons->size()==1)&&(posElectrons->size()==1)) {

                                                                selMuon m1 = negMuons->at(0);
                                                                selMuon m2 = posMuons->at(0);
                                                                selElectron e1 = negElectrons->at(0);
                                                                selElectron e2 = posElectrons->at(0);

                                                                float m1m2 = abs((m1.vec+m2.vec).M()-ZMass);
                                                                float e1e2 = abs((e1.vec+e2.vec).M()-ZMass);

                                                                if(m1m2<e1e2) {
                                                                        if(m1.p.Pt()>m2.p.Pt()) {
                                                                                lep1.setAll(m1);
                                                                                lep2.setAll(m2);
                                                                        }else{
                                                                                lep1.setAll(m2);
                                                                                lep2.setAll(m1);
                                                                        }
                                                                        if(e1.p.Pt()>e2.p.Pt()) {
                                                                                lep3.setAll(e1);
                                                                                lep4.setAll(e2);
                                                                        }else{
                                                                                lep3.setAll(e2);
                                                                                lep4.setAll(e1);
                                                                        }

                                                                        float mll1=(lep1.vec+lep2.vec).M();
                                                                        float mll2=(lep3.vec+lep4.vec).M();
                                                                        if(mll1<106. && mll1>76. && mll2<130. && mll2>50.) {
                                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][nom],false,lep1,lep2,lep3,lep4,true);
                                                                                FillerSignal(s1Maps.at(signalPoint).at(controlregionZZ).at(LL).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                        }

                                                                }else{
                                                                        if(m1.p.Pt()>m2.p.Pt()) {
                                                                                lep3.setAll(m1);
                                                                                lep4.setAll(m2);
                                                                        }else{
                                                                                lep3.setAll(m2);
                                                                                lep4.setAll(m1);
                                                                        }
                                                                        if(e1.p.Pt()>e2.p.Pt()) {
                                                                                lep1.setAll(e1);
                                                                                lep2.setAll(e2);
                                                                        }else{
                                                                                lep1.setAll(e2);
                                                                                lep2.setAll(e1);
                                                                        }

                                                                        float mll1=(lep1.vec+lep2.vec).M();
                                                                        float mll2=(lep3.vec+lep4.vec).M();
                                                                        if(mll1<106. && mll1>76. && mll2<130. && mll2>50.) {
                                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][nom],false,lep1,lep2,lep3,lep4,true);
                                                                                FillerSignal(s1Maps.at(signalPoint).at(controlregionZZ).at(LL).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                        }
                                                                }
                                                        }
                                                }
                                                //float mll1=(lep1.vec+lep2.vec).M();
                                                //float mll2=(lep3.vec+lep4.vec).M();
                                                //if(mll1<106. && mll1>76. && mll2<130. && mll2>50.){
                                                //FillerSignal(s1Maps.at(signalPoint).at(controlregionZZ).at(LL).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerZZ(cr1Maps[controlregionZZ][LL][nom],false,lep1,lep2,lep3,lep4,true);
                                                //}
                                        }
                                }
                        }
                }
                if(SelectEvent(ControlRegionZZ)) {
                        if((*selLeptonSize==3)&&(*matchedLeptonSize==2)) {
                                if((*selMuonSize==3)||(*selElectronSize==3)||(*selMuonSize==1 && *selElectronSize==2)||(*selMuonSize==2 && *selElectronSize==1)) {

                                        float ZMass = 91.1876;


                                        if((*countNegCharge==1 && *countPosCharge==2)||(*countNegCharge==2 && *countPosCharge==1)) {
                                                //find lepton combinations

                                                selLepton lep1,lep2,lep3;
                                                TLorentzVector temp_JESu(0.,0.,0.,0.);
                                                temp_JESu.SetPtEtaPhiM(met_JESu->p.Pt(),met_JESu->p.Eta(),met_JESu->p.Phi(),0.);
                                                TLorentzVector temp_JESd(0.,0.,0.,0.);
                                                temp_JESd.SetPtEtaPhiM(met_JESd->p.Pt(),met_JESd->p.Eta(),met_JESd->p.Phi(),0.);
                                                TLorentzVector temp_JERu(0.,0.,0.,0.);
                                                temp_JERu.SetPtEtaPhiM(met_JERu->p.Pt(),met_JERu->p.Eta(),met_JERu->p.Phi(),0.);
                                                TLorentzVector temp_JERd(0.,0.,0.,0.);
                                                temp_JERd.SetPtEtaPhiM(met_JERd->p.Pt(),met_JERd->p.Eta(),met_JERd->p.Phi(),0.);


                                                if(*selMuonSize==3) {

                                                        if((negMuons->size()==2)&&(posMuons->size()==1)) {

                                                                selMuon neg1 = negMuons->at(0);
                                                                selMuon neg2 = negMuons->at(1);
                                                                selMuon pos1 = posMuons->at(0);

                                                                //there are 2 combinations n1p1/n2 - n2p1/n1
                                                                float n1p1 = abs((neg1.vec+pos1.vec).M()-ZMass);
                                                                float n2p1 = abs((neg2.vec+pos1.vec).M()-ZMass);

                                                                if(n1p1<n2p1) {
                                                                        if(neg1.p.Pt()>pos1.p.Pt()) {
                                                                                lep1.setAll(neg1);
                                                                                lep2.setAll(pos1);
                                                                        }else{
                                                                                lep1.setAll(pos1);
                                                                                lep2.setAll(neg1);
                                                                        }
                                                                        lep3.setAll(neg2);
                                                                }else{
                                                                        if(n2p1<n1p1) {
                                                                                if(neg2.p.Pt()>pos1.p.Pt()) {
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
                                                        if((negMuons->size()==1)&&(posMuons->size()==2)) {

                                                                selMuon neg1 = negMuons->at(0);
                                                                selMuon pos1 = posMuons->at(0);
                                                                selMuon pos2 = posMuons->at(1);

                                                                //there are 2 combinations n1p1/p2 - n1p2/p1
                                                                float n1p1 = abs((neg1.vec+pos1.vec).M()-ZMass);
                                                                float n1p2 = abs((neg1.vec+pos2.vec).M()-ZMass);

                                                                if(n1p1<n1p2) {
                                                                        if(neg1.p.Pt()>pos1.p.Pt()) {
                                                                                lep1.setAll(neg1);
                                                                                lep2.setAll(pos1);
                                                                        }else{
                                                                                lep1.setAll(pos1);
                                                                                lep2.setAll(neg1);
                                                                        }
                                                                        lep3.setAll(pos2);
                                                                }else{
                                                                        if(n1p2<n1p1) {
                                                                                if(neg1.p.Pt()>pos2.p.Pt()) {
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
                                                        if(mll1<106. && mll1>76.) {
                                                                if( (lep1.vec+ETmiss_vec->vec).Mt() > 50. ) {
                                                                        if(*ETmiss>70.) {
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][nom],false,lep1,lep2,lep3,true);
                                                                                FillerSignal(s1Maps.at(signalPoint).at(controlregionWZ).at(LL).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                        }
                                                                }
                                                        }

                                                }


                                                if(*selElectronSize==3) {

                                                        if((negElectrons->size()==2)&&(posElectrons->size()==1)) {

                                                                selElectron neg1 = negElectrons->at(0);
                                                                selElectron neg2 = negElectrons->at(1);
                                                                selElectron pos1 = posElectrons->at(0);

                                                                //there are 2 combinations n1p1/n2 - n2p1/n1
                                                                float n1p1 = abs((neg1.vec+pos1.vec).M()-ZMass);
                                                                float n2p1 = abs((neg2.vec+pos1.vec).M()-ZMass);

                                                                if(n1p1<n2p1) {
                                                                        if(neg1.p.Pt()>pos1.p.Pt()) {
                                                                                lep1.setAll(neg1);
                                                                                lep2.setAll(pos1);
                                                                        }else{
                                                                                lep1.setAll(pos1);
                                                                                lep2.setAll(neg1);
                                                                        }
                                                                        lep3.setAll(neg2);
                                                                }else{
                                                                        if(n2p1<n1p1) {
                                                                                if(neg2.p.Pt()>pos1.p.Pt()) {
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
                                                        if((negElectrons->size()==1)&&(posElectrons->size()==2)) {

                                                                selElectron neg1 = negElectrons->at(0);
                                                                selElectron pos1 = posElectrons->at(0);
                                                                selElectron pos2 = posElectrons->at(1);

                                                                //there are 2 combinations n1p1/p2 - n1p2/p1
                                                                float n1p1 = abs((neg1.vec+pos1.vec).M()-ZMass);
                                                                float n1p2 = abs((neg1.vec+pos2.vec).M()-ZMass);

                                                                if(n1p1<n1p2) {
                                                                        if(neg1.p.Pt()>pos1.p.Pt()) {
                                                                                lep1.setAll(neg1);
                                                                                lep2.setAll(pos1);
                                                                        }else{
                                                                                lep1.setAll(pos1);
                                                                                lep2.setAll(neg1);
                                                                        }
                                                                        lep3.setAll(pos2);
                                                                }else{
                                                                        if(n1p2<n1p1) {
                                                                                if(neg1.p.Pt()>pos2.p.Pt()) {
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
                                                        if(mll1<106. && mll1>76.) {
                                                                if( (lep1.vec+ETmiss_vec->vec).Mt() > 50. ) {
                                                                        if(*ETmiss>70.) {
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][nom],false,lep1,lep2,lep3,true);
                                                                                FillerSignal(s1Maps.at(signalPoint).at(controlregionWZ).at(LL).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                        }
                                                                }
                                                        }
                                                }

                                                if((*selElectronSize==2)&&(*selMuonSize==1)) {

                                                        if(negElectrons->size()==1 && posElectrons->size()==1 && posMuons->size()==1) {
                                                                selElectron neg1 = negElectrons->at(0);
                                                                selElectron pos1 = posElectrons->at(0);
                                                                selMuon pos2 = posMuons->at(0);

                                                                if(neg1.p.Pt()>pos1.p.Pt()) {
                                                                        lep1.setAll(neg1);
                                                                        lep2.setAll(pos1);
                                                                }else{
                                                                        lep1.setAll(pos1);
                                                                        lep2.setAll(neg1);
                                                                }
                                                                lep3.setAll(pos2);
                                                        }
                                                        if(negElectrons->size()==1 && posElectrons->size()==1 && negMuons->size()==1) {
                                                                selElectron neg1 = negElectrons->at(0);
                                                                selElectron pos1 = posElectrons->at(0);
                                                                selMuon neg2 = negMuons->at(0);

                                                                if(neg1.p.Pt()>pos1.p.Pt()) {
                                                                        lep1.setAll(neg1);
                                                                        lep2.setAll(pos1);
                                                                }else{
                                                                        lep1.setAll(pos1);
                                                                        lep2.setAll(neg1);
                                                                }
                                                                lep3.setAll(neg2);
                                                        }

                                                        float mll1=(lep1.vec+lep2.vec).M();
                                                        if(mll1<106. && mll1>76.) {
                                                                if( (lep1.vec+ETmiss_vec->vec).Mt() > 50. ) {
                                                                        if(*ETmiss>70.) {
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][nom],false,lep1,lep2,lep3,true);
                                                                                FillerSignal(s1Maps.at(signalPoint).at(controlregionWZ).at(LL).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                        }
                                                                }
                                                        }

                                                }

                                                if((*selElectronSize==1)&&(*selMuonSize==2)) {

                                                        if(negMuons->size()==1 && posMuons->size()==1 && posElectrons->size()==1) {
                                                                selMuon neg1 = negMuons->at(0);
                                                                selMuon pos1 = posMuons->at(0);
                                                                selElectron pos2 = posElectrons->at(0);

                                                                if(neg1.p.Pt()>pos1.p.Pt()) {
                                                                        lep1.setAll(neg1);
                                                                        lep2.setAll(pos1);
                                                                }else{
                                                                        lep1.setAll(pos1);
                                                                        lep2.setAll(neg1);
                                                                }
                                                                lep3.setAll(pos2);
                                                        }
                                                        if(negMuons->size()==1 && posMuons->size()==1 && negElectrons->size()==1) {
                                                                selMuon neg1 = negMuons->at(0);
                                                                selMuon pos1 = posMuons->at(0);
                                                                selElectron neg2 = negElectrons->at(0);

                                                                if(neg1.p.Pt()>pos1.p.Pt()) {
                                                                        lep1.setAll(neg1);
                                                                        lep2.setAll(pos1);
                                                                }else{
                                                                        lep1.setAll(pos1);
                                                                        lep2.setAll(neg1);
                                                                }
                                                                lep3.setAll(neg2);
                                                        }

                                                        float mll1=(lep1.vec+lep2.vec).M();
                                                        if(mll1<106. && mll1>76.) {
                                                                if( (lep1.vec+ETmiss_vec->vec).Mt() > 50. ) {
                                                                        if(*ETmiss>70.) {
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][nom],false,lep1,lep2,lep3,true);
                                                                                FillerSignal(s1Maps.at(signalPoint).at(controlregionWZ).at(LL).at(nom),1.,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
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


void myAnalyzer::FillTriggerStudies(){
        if(SelectEventTriggerStudies(TRIGDILEP)) {
                if(config_doTrigger) {
                        if (*trigHt) {//baselineTrigger
                                if(*isDiElectron && !(*isDiMuon || *isMuonElectron || *isMuonElectron)) {
                                        FillerTrigger(eff1Maps[trigdilep][EE],false,(*trigDiEle&&*trigDiEleMatch));
                                }
                                if(*isDiMuon && !(*isDiElectron || *isMuonElectron || *isMuonElectron)) {
                                        FillerTrigger(eff1Maps[trigdilep][MM],false,(*trigDiMu&&*trigDiMuMatch));
                                }
                                if((*isMuonElectron||*isElectronMuon) && !(*isDiElectron || *isDiMuon)) {
                                        FillerTrigger(eff1Maps[trigdilep][EM],false,(*trigMuEle&&*trigMuEleMatch));
                                }
                        }
                }
                if(config_doTriggerPure) {
                        //if (*trigHt){//baselineTrigger
                        if(*isDiElectron && !(*isDiMuon || *isMuonElectron || *isMuonElectron)) {
                                FillerTrigger(eff1Maps[trigdilep][EE],false,(*trigDiEle&&*trigDiEleMatch));
                        }
                        if(*isDiMuon && !(*isDiElectron || *isMuonElectron || *isMuonElectron)) {
                                FillerTrigger(eff1Maps[trigdilep][MM],false,(*trigDiMu&&*trigDiMuMatch));
                        }
                        if((*isMuonElectron||*isElectronMuon) && !(*isDiElectron || *isDiMuon)) {
                                FillerTrigger(eff1Maps[trigdilep][EM],false,(*trigMuEle&&*trigMuEleMatch));
                        }
                        //}
                }
                if(config_doTriggerMET) {
                        if (*trigMET) {//baselineTrigger
                                if(*isDiElectron && !(*isDiMuon || *isMuonElectron || *isMuonElectron)) {
                                        FillerTrigger(eff1Maps[trigdilep][EE],false,(*trigDiEle&&*trigDiEleMatch));
                                }
                                if(*isDiMuon && !(*isDiElectron || *isMuonElectron || *isMuonElectron)) {
                                        FillerTrigger(eff1Maps[trigdilep][MM],false,(*trigDiMu&&*trigDiMuMatch));
                                }
                                if((*isMuonElectron||*isElectronMuon) && !(*isDiElectron || *isDiMuon)) {
                                        FillerTrigger(eff1Maps[trigdilep][EM],false,(*trigMuEle&&*trigMuEleMatch));
                                }
                        }
                }
        }

        if(SelectEventTriggerStudies(TRIGDILEP_ptcuts)) {
                if(config_doTrigger) {
                        if (*trigHt) {//baselineTrigger
                                if(*isDiElectron && !(*isDiMuon || *isMuonElectron || *isMuonElectron)) {
                                        FillerTrigger(eff1Maps[trigdilep_ptcuts][EE],false,(*trigDiEle&&*trigDiEleMatch));
                                }
                                if(*isDiMuon && !(*isDiElectron || *isMuonElectron || *isMuonElectron)) {
                                        FillerTrigger(eff1Maps[trigdilep_ptcuts][MM],false,(*trigDiMu&&*trigDiMuMatch));
                                }
                                if((*isMuonElectron||*isElectronMuon) && !(*isDiElectron || *isDiMuon)) {
                                        FillerTrigger(eff1Maps[trigdilep_ptcuts][EM],false,(*trigMuEle&&*trigMuEleMatch));
                                }
                        }
                }
                if(config_doTriggerPure) {
                        //if (*trigHt){//baselineTrigger
                        if(*isDiElectron && !(*isDiMuon || *isMuonElectron || *isMuonElectron)) {
                                FillerTrigger(eff1Maps[trigdilep_ptcuts][EE],false,(*trigDiEle&&*trigDiEleMatch));
                        }
                        if(*isDiMuon && !(*isDiElectron || *isMuonElectron || *isMuonElectron)) {
                                FillerTrigger(eff1Maps[trigdilep_ptcuts][MM],false,(*trigDiMu&&*trigDiMuMatch));
                        }
                        if((*isMuonElectron||*isElectronMuon) && !(*isDiElectron || *isDiMuon)) {
                                FillerTrigger(eff1Maps[trigdilep_ptcuts][EM],false,(*trigMuEle&&*trigMuEleMatch));
                        }
                        //}
                }
                if(config_doTriggerMET) {
                        if (*trigMET) {//baselineTrigger
                                if(*isDiElectron && !(*isDiMuon || *isMuonElectron || *isMuonElectron)) {
                                        FillerTrigger(eff1Maps[trigdilep_ptcuts][EE],false,(*trigDiEle&&*trigDiEleMatch));
                                }
                                if(*isDiMuon && !(*isDiElectron || *isMuonElectron || *isMuonElectron)) {
                                        FillerTrigger(eff1Maps[trigdilep_ptcuts][MM],false,(*trigDiMu&&*trigDiMuMatch));
                                }
                                if((*isMuonElectron||*isElectronMuon) && !(*isDiElectron || *isDiMuon)) {
                                        FillerTrigger(eff1Maps[trigdilep_ptcuts][EM],false,(*trigMuEle&&*trigMuEleMatch));
                                }
                        }
                }
        }
}



void myAnalyzer::FillHistograms(){

        if(config_selectionsToProcessMap[DILEP]) {
                if(SelectEvent(DILEP)) {
                        if(*selLeptonSize==*matchedLeptonSize && *selLeptonSize==2) {
                                if((*isDiMuon || *isDiElectron) && (!(*isElectronMuon || *isMuonElectron))) {
                                        Filler(h1Maps[dilep][LL],false);
                                }
                                if (*isDiElectron && (!*isDiMuon) && (!(*isElectronMuon || *isMuonElectron))) {
                                        Filler(h1Maps[dilep][EE],false);
                                }else{
                                        if (*isDiMuon && !(*isMuonElectron||*isElectronMuon)) {Filler(h1Maps[dilep][MM],false);}
                                        if (*isElectronMuon || *isMuonElectron) Filler(h1Maps[dilep][EM],false);
                                }
                        }
                }
        }

        if(config_selectionsToProcessMap[SEL]) {
                if(SelectEvent(SEL)) {
                        if(*selLeptonSize==*matchedLeptonSize && *selLeptonSize==2) {
                                if (*selPhotonSize!=0) {
                                        if(*isDiMuon || *isDiElectron) {
                                                Filler(h1Maps[sel][LL],true);
                                        }
                                        if (*isDiElectron) {
                                                Filler(h1Maps[sel][EE],true);
                                        }else{
                                                if (*isDiMuon) {
                                                        Filler(h1Maps[sel][MM],true);
                                                }
                                                if (*isElectronMuon || *isMuonElectron) {
                                                        Filler(h1Maps[sel][EM],true);
                                                }
                                        }
                                }
                        }
                }
        }




        if(config_selectionsToProcessMap[ONZMET]) {

                if(SelectEvent(ONZ)) {
                        if(*selLeptonSize==*matchedLeptonSize && *selLeptonSize==2) {

                                if ((*selPhotonSize!=0)&&(*mll>81. && *mll<101.)&&(*ETmiss<100.)) { //ONZ+MET<100
                                        if(*isDiMuon || *isDiElectron) {
                                                Filler(h1Maps[onzmet0100][LL],true);
                                        }
                                        if(*isDiElectron) {
                                                Filler(h1Maps[onzmet0100][EE],true);
                                        }else{
                                                if (*isDiMuon) {
                                                        Filler(h1Maps[onzmet0100][MM],true);
                                                }
                                                if (*isElectronMuon || *isMuonElectron) {
                                                        Filler(h1Maps[onzmet0100][EM],true);
                                                }
                                        }
                                }

                                if ((*selPhotonSize!=0)&&(*mll>81. && *mll<101.)&&(*ETmiss>=150.)) { //ONZ+MET>150
                                        if(*isDiMuon || *isDiElectron) {
                                                Filler(h1Maps[onzmet150][LL],true);
                                                //Filler2D(h2Maps[onzmet150][LL],true);
                                                if (*isDiElectron) {
                                                        Filler(h1Maps[onzmet150][EE],true);
                                                        //Filler2D(h2Maps[onzmet150][EE],true);
                                                }
                                                if (*isDiMuon) {
                                                        Filler(h1Maps[onzmet150][MM],true);
                                                        //Filler2D(h2Maps[onzmet150][MM],true);
                                                }
                                        }
                                        if (*isElectronMuon || *isMuonElectron) {
                                                Filler(h1Maps[onzmet150][EM],true);
                                                //Filler2D(h2Maps[onzmet150][EM],true);
                                        }
                                }

                                if ((*selPhotonSize!=0)&&(*mll>81. && *mll<101.)&&(*ETmiss<150.)&&(*ETmiss>=100.)) { //ONZ+MET>100<150
                                        if(*isDiMuon || *isDiElectron) {
                                                Filler(h1Maps[onzmet100150][LL],true);
                                                if (*isDiElectron) {
                                                        Filler(h1Maps[onzmet100150][EE],true);
                                                }
                                                if (*isDiMuon) {
                                                        Filler(h1Maps[onzmet100150][MM],true);
                                                }
                                        }
                                        if (*isElectronMuon || *isMuonElectron) {
                                                Filler(h1Maps[onzmet100150][EM],true);

                                        }
                                }
                        }
                }
        }


        clearCutFlowMap();
        if(config_selectionsToProcessMap[ControlRegionDY]) {
                if(SelectEvent(ONZ)) {
                        if(*selLeptonSize==*matchedLeptonSize && *selLeptonSize==2) {
                                if ((*selPhotonSize!=0)&&(*mll>81. && *mll<101.)) {//&&(selectedEvent.ETmiss<100.)){ //ONZ+MET<100 ~ CR DY/Z(+gamma)
                                        if(*ETmiss<100.) {
                                                if(*isDiMuon || *isDiElectron) {
                                                        Filler(cr1Maps[controlregionDY][LL][nom],true,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        Filler(cr1Maps[controlregionDY][LL][ISRUP],true,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                        Filler(cr1Maps[controlregionDY][LL][ISRDOWN],true,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                        Filler(cr1Maps[controlregionDY][LL][EWKUP],true,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                        Filler(cr1Maps[controlregionDY][LL][EWKDOWN],true,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                        Filler(cr1Maps[controlregionDY][LL][LEPSFUP],true,true,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        Filler(cr1Maps[controlregionDY][LL][LEPSFDOWN],true,true,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        Filler(cr1Maps[controlregionDY][LL][PHOTONSFUP],true,true,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                        Filler(cr1Maps[controlregionDY][LL][PHOTONSFDOWN],true,true,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                        Filler(cr1Maps[controlregionDY][LL][PUUP],true,true,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        Filler(cr1Maps[controlregionDY][LL][PUDOWN],true,true,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                }
                                                //if (*isDiElectron){
                                                //Filler(cr1Maps[controlregionDY][EE][nom],true,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //Filler(cr1Maps[controlregionDY][EE][ISRUP],true,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                //Filler(cr1Maps[controlregionDY][EE][ISRDOWN],true,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                //Filler(cr1Maps[controlregionDY][EE][EWKUP],true,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                //Filler(cr1Maps[controlregionDY][EE][EWKDOWN],true,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                //Filler(cr1Maps[controlregionDY][EE][LEPSFUP],true,true,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //Filler(cr1Maps[controlregionDY][EE][LEPSFDOWN],true,true,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //Filler(cr1Maps[controlregionDY][EE][PHOTONSFUP],true,true,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                //Filler(cr1Maps[controlregionDY][EE][PHOTONSFDOWN],true,true,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                //Filler(cr1Maps[controlregionDY][EE][PUUP],true,true,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //Filler(cr1Maps[controlregionDY][EE][PUDOWN],true,true,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //}
                                                else{
                                                        //if (*isDiMuon){
                                                        //Filler(cr1Maps[controlregionDY][MM][nom],true,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        //Filler(cr1Maps[controlregionDY][MM][ISRUP],true,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                        //Filler(cr1Maps[controlregionDY][MM][ISRDOWN],true,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                        //Filler(cr1Maps[controlregionDY][MM][EWKUP],true,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                        //Filler(cr1Maps[controlregionDY][MM][EWKDOWN],true,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                        //Filler(cr1Maps[controlregionDY][MM][LEPSFUP],true,true,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        //Filler(cr1Maps[controlregionDY][MM][LEPSFDOWN],true,true,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        //Filler(cr1Maps[controlregionDY][MM][PHOTONSFUP],true,true,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                        //Filler(cr1Maps[controlregionDY][MM][PHOTONSFDOWN],true,true,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                        //Filler(cr1Maps[controlregionDY][MM][PUUP],true,true,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        //Filler(cr1Maps[controlregionDY][MM][PUDOWN],true,true,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                        //}
                                                        //if (*isElectronMuon || *isMuonElectron){
                                                        //Filler(cr1Maps[controlregionDY][EM][nom],true,true);
                                                        //}
                                                }

                                                //for(int i=0; i<nWeights; i++){
                                                for(int i=config_pdfStart; i<(config_pdfStart+config_pdfStep); i++) {
                                                        if(*isDiMuon || *isDiElectron) {
                                                                Filler(cr1Maps[controlregionDY][LL][PDFNAMES[i]],true,true,i);
                                                        }
                                                        //if (*isDiElectron){
                                                        //Filler(cr1Maps[controlregionDY][EE][PDFNAMES[i]],true,true,i);
                                                        //}else{
                                                        //if (*isDiMuon){
                                                        //Filler(cr1Maps[controlregionDY][MM][PDFNAMES[i]],true,true,i);
                                                        //}
                                                        //if (*isElectronMuon || *isMuonElectron){
                                                        //Filler(cr1Maps[controlregionDY][EM][PDFNAMES[i]],true,true,i);
                                                        //}
                                                        //}
                                                }
                                        }
                                        if((met_JESu->p.Pt()<100.)) { //ONZ+MET<100 ~ CR DY/Z(+gamma)
                                                if(*isDiMuon || *isDiElectron) {
                                                        Filler(cr1Maps[controlregionDY][LL][JESu],true,true,9999,JESUP);
                                                }
                                                //if (*isDiElectron){
                                                //Filler(cr1Maps[controlregionDY][EE][JESu],true,true,9999,JESUP);
                                                //}else{
                                                //if (*isDiMuon){
                                                //Filler(cr1Maps[controlregionDY][MM][JESu],true,true,9999,JESUP);
                                                //}
                                                //if (*isElectronMuon || *isMuonElectron){
                                                //Filler(cr1Maps[controlregionDY][EM][JESu],true,true,9999,JESUP);
                                                //}
                                                //}
                                        }
                                        if((met_JESd->p.Pt()<100.)) { //ONZ+MET<100 ~ CR DY/Z(+gamma)
                                                if(*isDiMuon || *isDiElectron) {
                                                        Filler(cr1Maps[controlregionDY][LL][JESd],true,true,9999,JESDOWN);
                                                }
                                                //if (*isDiElectron){
                                                //Filler(cr1Maps[controlregionDY][EE][JESd],true,true,9999,JESDOWN);
                                                //}else{
                                                //if (*isDiMuon){
                                                //Filler(cr1Maps[controlregionDY][MM][JESd],true,true,9999,JESDOWN);
                                                //}
                                                //if (*isElectronMuon || *isMuonElectron){
                                                //Filler(cr1Maps[controlregionDY][EM][JESd],true,true,9999,JESDOWN);
                                                //}
                                                //}
                                        }
                                        if((met_JERu->p.Pt()<100.)) { //ONZ+MET<100 ~ CR DY/Z(+gamma)
                                                if(*isDiMuon || *isDiElectron) {
                                                        Filler(cr1Maps[controlregionDY][LL][JERu],true,true,9999,JERUP);
                                                }
                                                //if (*isDiElectron){
                                                //Filler(cr1Maps[controlregionDY][EE][JERu],true,true,9999,JERUP);
                                                //}else{
                                                //if (*isDiMuon){
                                                //Filler(cr1Maps[controlregionDY][MM][JERu],true,true,9999,JERUP);
                                                //}
                                                //if (*isElectronMuon || *isMuonElectron){
                                                //Filler(cr1Maps[controlregionDY][EM][JERu],true,true,9999,JERUP);
                                                //}
                                                //}
                                        }
                                        if((met_JESd->p.Pt()<100.)) { //ONZ+MET<100 ~ CR DY/Z(+gamma)
                                                if(*isDiMuon || *isDiElectron) {
                                                        Filler(cr1Maps[controlregionDY][LL][JERd],true,true,9999,JERDOWN);
                                                }
                                                //if(*isDiElectron){
                                                //Filler(cr1Maps[controlregionDY][EE][JERd],true,true,9999,JERDOWN);
                                                //}else{
                                                //if (*isDiMuon){
                                                //Filler(cr1Maps[controlregionDY][MM][JERd],true,true,9999,JERDOWN);
                                                //}
                                                //if (*isElectronMuon || *isMuonElectron){
                                                //Filler(cr1Maps[controlregionDY][EM][JERd],true,true,9999,JERDOWN);
                                                //}
                                                //}
                                        }
                                }
                        }
                }
        }

        clearCutFlowMap();
        if(config_selectionsToProcessMap[ValidationRegion]) {
                if(SelectEvent(ONZ)) {
                        if(*selLeptonSize==*matchedLeptonSize && *selLeptonSize==2) {
                                //if ((*selPhotonSize!=0)&&(*mll>81. && *mll<101.)&&(*ETmiss>=100.)&&(*ETmiss<150.)){ //ONZ+MET<150>100 ~ VR
                                //if ((*selPhotonSize!=0)&&(*mll>81. && *mll<101.)&&(*ETmiss>=100.)&&((*ETmiss<150.)||(*ETmiss>=150. && *Mt2<100.))){ //ONZ+MET<150>100 ~ VR
                                if ((*selPhotonSize!=0)&&(*mll>81. && *mll<101.)&&(*ETmiss>=100.)&&((*ETmiss<150.)||(*Mt2<100.))) { //ONZ+MET<150>100 ~ VR
                                        //if(selPhotons->at(0).p.Pt()<80.){
                                        if(*isDiMuon || *isDiElectron) {
                                                Filler(h1Maps[validationregion][LL],true);
                                                if(selPhotons->at(0).p.Pt()>=80.) {
                                                        Filler(h1Maps[validationregion80][LL],true);
                                                }else{
                                                        Filler(h1Maps[validationregion080][LL],true);
                                                }
                                        }
                                        //if (*isDiElectron){
                                        //Filler(h1Maps[validationregion][EE],true);
                                        //if(selPhotons->at(0).p.Pt()>=80.){
                                        //Filler(h1Maps[validationregion80][EE],true);
                                        //}else{
                                        //Filler(h1Maps[validationregion080][EE],true);
                                        //}
                                        //}
                                        //if (*isDiMuon){
                                        //Filler(h1Maps[validationregion][MM],true);
                                        //if(selPhotons->at(0).p.Pt()>=80.){
                                        //Filler(h1Maps[validationregion80][MM],true);
                                        //}else{
                                        //Filler(h1Maps[validationregion080][MM],true);
                                        //}
                                        //}
                                        //if (*isElectronMuon || *isMuonElectron){
                                        //Filler(h1Maps[validationregion][EM],true);
                                        //if(selPhotons->at(0).p.Pt()>=80.){
                                        //Filler(h1Maps[validationregion80][EM],true);
                                        //}else{
                                        //Filler(h1Maps[validationregion080][EM],true);
                                        //}
                                        //}
                                        //}
                                }
                        }
                }
        }

        clearCutFlowMap();

        if(config_selectionsToProcessMap[ControlRegionTT]) {
                if(SelectEvent(SEL)) {
                        if(*selLeptonSize==*matchedLeptonSize && *selLeptonSize==2) {
                                if ((*selPhotonSize!=0)) { //Different Flavor + 1Photon ~ CR TT(+gamma)
                                        if (*isElectronMuon || *isMuonElectron) {
                                                Filler(cr1Maps[controlregionTT][EM][nom],true,true);
                                                Filler(cr1Maps[controlregionTT][EM][ISRUP],true,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                Filler(cr1Maps[controlregionTT][EM][ISRDOWN],true,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                Filler(cr1Maps[controlregionTT][EM][EWKUP],true,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                Filler(cr1Maps[controlregionTT][EM][EWKDOWN],true,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                Filler(cr1Maps[controlregionTT][EM][JESu],true,true,9999,JESUP);
                                                Filler(cr1Maps[controlregionTT][EM][JESd],true,true,9999,JESDOWN);
                                                Filler(cr1Maps[controlregionTT][EM][JERu],true,true,9999,JERUP);
                                                Filler(cr1Maps[controlregionTT][EM][JERd],true,true,9999,JERDOWN);
                                                Filler(cr1Maps[controlregionTT][EM][LEPSFUP],true,true,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                Filler(cr1Maps[controlregionTT][EM][LEPSFDOWN],true,true,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                Filler(cr1Maps[controlregionTT][EM][PHOTONSFUP],true,true,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                Filler(cr1Maps[controlregionTT][EM][PHOTONSFDOWN],true,true,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                Filler(cr1Maps[controlregionTT][EM][PUUP],true,true,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                Filler(cr1Maps[controlregionTT][EM][PUDOWN],true,true,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //for(int i=0; i<nWeights; i++){
                                                for(int i=config_pdfStart; i<config_pdfStart+config_pdfStep; i++) {
                                                        Filler(cr1Maps[controlregionTT][EM][PDFNAMES[i]],true,true,i);
                                                }
                                                if(selPhotons->at(0).p.Pt()>=80.) {
                                                        Filler(cr1Maps[controlregionTT80][EM][nom],true,true);
                                                }else{
                                                        Filler(cr1Maps[controlregionTT080][EM][nom],true,true);
                                                }
                                        }
                                }
                        }
                }
        }

        clearCutFlowMap();
        if(config_selectionsToProcessMap[ControlRegionZZ]) {
                if(SelectEvent(ControlRegionZZ)) {
                        if((*selLeptonSize==4)&&(*matchedLeptonSize==2)) {
                                if((*selMuonSize==4)||(*selElectronSize==4)||(*selMuonSize==2 && *selElectronSize==2)) {


                                        float ZMass = 91.1876;

                                        if(*countNegCharge==2 && *countPosCharge==2) {

                                                selLepton lep1,lep2,lep3,lep4;

                                                if(*selMuonSize==4) {

                                                        if((negMuons->size()==2)&&(posMuons->size()==2)) {
                                                                selMuon neg1 = negMuons->at(0);
                                                                selMuon neg2 = negMuons->at(1);
                                                                selMuon pos1 = posMuons->at(0);
                                                                selMuon pos2 = posMuons->at(1);

                                                                //there are 2 combinations n1p1/n2p2 - n1p2/n2p1
                                                                float n1p1 = abs((neg1.vec+pos1.vec).M()-ZMass);
                                                                float n1p2 = abs((neg1.vec+pos2.vec).M()-ZMass);
                                                                float n2p1 = abs((neg2.vec+pos1.vec).M()-ZMass);
                                                                float n2p2 = abs((neg2.vec+pos2.vec).M()-ZMass);

                                                                if((n1p1<n1p2)&&(n1p1<n2p1)&&(n1p1<n2p2)) {
                                                                        if(neg1.p.Pt()>pos1.p.Pt()) {
                                                                                lep1.setAll(neg1);
                                                                                lep2.setAll(pos1);
                                                                        }else{
                                                                                lep1.setAll(pos1);
                                                                                lep2.setAll(neg1);
                                                                        }
                                                                        if(neg2.p.Pt()>pos2.p.Pt()) {
                                                                                lep3.setAll(neg2);
                                                                                lep4.setAll(pos2);
                                                                        }else{
                                                                                lep3.setAll(pos2);
                                                                                lep4.setAll(neg2);
                                                                        }
                                                                }else{
                                                                        if((n1p2<n1p1)&&(n1p2<n2p1)&&(n1p2<n2p2)) {
                                                                                if(neg1.p.Pt()>pos2.p.Pt()) {
                                                                                        lep1.setAll(neg1);
                                                                                        lep2.setAll(pos2);
                                                                                }else{
                                                                                        lep1.setAll(pos2);
                                                                                        lep2.setAll(neg1);
                                                                                }
                                                                                if(neg2.p.Pt()>pos1.p.Pt()) {
                                                                                        lep3.setAll(neg2);
                                                                                        lep4.setAll(pos1);
                                                                                }else{
                                                                                        lep3.setAll(pos1);
                                                                                        lep4.setAll(neg2);
                                                                                }
                                                                        }else{
                                                                                if((n2p1<n1p1)&&(n2p1<n1p2)&&(n2p1<n2p2)) {
                                                                                        if(neg2.p.Pt()>pos1.p.Pt()) {
                                                                                                lep1.setAll(neg2);
                                                                                                lep2.setAll(pos1);
                                                                                        }else{
                                                                                                lep1.setAll(pos1);
                                                                                                lep2.setAll(neg2);
                                                                                        }
                                                                                        if(neg1.p.Pt()>pos2.p.Pt()) {
                                                                                                lep3.setAll(neg1);
                                                                                                lep4.setAll(pos2);
                                                                                        }else{
                                                                                                lep3.setAll(pos2);
                                                                                                lep4.setAll(neg1);
                                                                                        }
                                                                                }else{
                                                                                        if(neg2.p.Pt()>pos2.p.Pt()) {
                                                                                                lep1.setAll(neg2);
                                                                                                lep2.setAll(pos2);
                                                                                        }else{
                                                                                                lep1.setAll(pos2);
                                                                                                lep2.setAll(neg2);
                                                                                        }
                                                                                        if(neg1.p.Pt()>pos1.p.Pt()) {
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
                                                        if(mll1<106. && mll1>76. && mll2<130. && mll2>50.) {
                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][nom],false,lep1,lep2,lep3,lep4,true);
                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][JESu],false,lep1,lep2,lep3,lep4,true,9999,JESUP);
                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][JESd],false,lep1,lep2,lep3,lep4,true,9999,JESDOWN);
                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][JERu],false,lep1,lep2,lep3,lep4,true,9999,JERUP);
                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][JERd],false,lep1,lep2,lep3,lep4,true,9999,JERDOWN);
                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][LEPSFUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][LEPSFDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][PHOTONSFUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][PHOTONSFDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][PUUP],false,lep1,lep2,lep3,lep4,true,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][PUDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][ISRUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][ISRDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][EWKUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][EWKDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][nom],false,lep1,lep2,lep3,lep4,true);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][JESu],false,lep1,lep2,lep3,lep4,true,9999,JESUP);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][JESd],false,lep1,lep2,lep3,lep4,true,9999,JESDOWN);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][JERu],false,lep1,lep2,lep3,lep4,true,9999,JERUP);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][JERd],false,lep1,lep2,lep3,lep4,true,9999,JERDOWN);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][LEPSFUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][LEPSFDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][PHOTONSFUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][PHOTONSFDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][PUUP],false,lep1,lep2,lep3,lep4,true,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][PUDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][ISRUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][ISRDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][EWKUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][EWKDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                                //for(int i=0; i<nWeights; i++){
                                                                for(int i=config_pdfStart; i<config_pdfStart+config_pdfStep; i++) {
                                                                        //FillerZZ(cr1Maps[controlregionZZ][MM][PDFNAMES[i]],false,lep1,lep2,lep3,lep4,true,i);
                                                                        FillerZZ(cr1Maps[controlregionZZ][LL][PDFNAMES[i]],false,lep1,lep2,lep3,lep4,true,i);
                                                                }
                                                        }
                                                }

                                                if(*selElectronSize==4) {

                                                        if((negElectrons->size()==2)&&(posElectrons->size()==2)) {
                                                                selElectron neg1 = negElectrons->at(0);
                                                                selElectron neg2 = negElectrons->at(1);
                                                                selElectron pos1 = posElectrons->at(0);
                                                                selElectron pos2 = posElectrons->at(1);

                                                                //there are 2 combinations n1p1/n2p2 - n1p2/n2p1
                                                                float n1p1 = abs((neg1.vec+pos1.vec).M()-ZMass);
                                                                float n1p2 = abs((neg1.vec+pos2.vec).M()-ZMass);
                                                                float n2p1 = abs((neg2.vec+pos1.vec).M()-ZMass);
                                                                float n2p2 = abs((neg2.vec+pos2.vec).M()-ZMass);


                                                                if((n1p1<n1p2)&&(n1p1<n2p1)&&(n1p1<n2p2)) {
                                                                        if(neg1.p.Pt()>pos1.p.Pt()) {
                                                                                lep1.setAll(neg1);
                                                                                lep2.setAll(pos1);
                                                                        }else{
                                                                                lep1.setAll(pos1);
                                                                                lep2.setAll(neg1);
                                                                        }
                                                                        if(neg2.p.Pt()>pos2.p.Pt()) {
                                                                                lep3.setAll(neg2);
                                                                                lep4.setAll(pos2);
                                                                        }else{
                                                                                lep3.setAll(pos2);
                                                                                lep4.setAll(neg2);
                                                                        }
                                                                }else{
                                                                        if((n1p2<n1p1)&&(n1p2<n2p1)&&(n1p2<n2p2)) {
                                                                                if(neg1.p.Pt()>pos2.p.Pt()) {
                                                                                        lep1.setAll(neg1);
                                                                                        lep2.setAll(pos2);
                                                                                }else{
                                                                                        lep1.setAll(pos2);
                                                                                        lep2.setAll(neg1);
                                                                                }
                                                                                if(neg2.p.Pt()>pos1.p.Pt()) {
                                                                                        lep3.setAll(neg2);
                                                                                        lep4.setAll(pos1);
                                                                                }else{
                                                                                        lep3.setAll(pos1);
                                                                                        lep4.setAll(neg2);
                                                                                }
                                                                        }else{
                                                                                if((n2p1<n1p1)&&(n2p1<n1p2)&&(n2p1<n2p2)) {
                                                                                        if(neg2.p.Pt()>pos1.p.Pt()) {
                                                                                                lep1.setAll(neg2);
                                                                                                lep2.setAll(pos1);
                                                                                        }else{
                                                                                                lep1.setAll(pos1);
                                                                                                lep2.setAll(neg2);
                                                                                        }
                                                                                        if(neg1.p.Pt()>pos2.p.Pt()) {
                                                                                                lep3.setAll(neg1);
                                                                                                lep4.setAll(pos2);
                                                                                        }else{
                                                                                                lep3.setAll(pos2);
                                                                                                lep4.setAll(neg1);
                                                                                        }
                                                                                }else{
                                                                                        if(neg2.p.Pt()>pos2.p.Pt()) {
                                                                                                lep1.setAll(neg2);
                                                                                                lep2.setAll(pos2);
                                                                                        }else{
                                                                                                lep1.setAll(pos2);
                                                                                                lep2.setAll(neg2);
                                                                                        }
                                                                                        if(neg1.p.Pt()>pos1.p.Pt()) {
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
                                                        if(mll1<106. && mll1>76. && mll2<130. && mll2>50.) {
                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][nom],false,lep1,lep2,lep3,lep4,true);
                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][JESu],false,lep1,lep2,lep3,lep4,true,9999,JESUP);
                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][JESd],false,lep1,lep2,lep3,lep4,true,9999,JESDOWN);
                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][JERu],false,lep1,lep2,lep3,lep4,true,9999,JERUP);
                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][JERd],false,lep1,lep2,lep3,lep4,true,9999,JERDOWN);
                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][LEPSFUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][LEPSFDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][PHOTONSFUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][PHOTONSFDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][PUUP],false,lep1,lep2,lep3,lep4,true,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][PUDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][ISRUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][ISRDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][EWKUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][EWKDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][nom],false,lep1,lep2,lep3,lep4,true);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][JESu],false,lep1,lep2,lep3,lep4,true,9999,JESUP);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][JESd],false,lep1,lep2,lep3,lep4,true,9999,JESDOWN);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][JERu],false,lep1,lep2,lep3,lep4,true,9999,JERUP);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][JERd],false,lep1,lep2,lep3,lep4,true,9999,JERDOWN);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][LEPSFUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][LEPSFDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][PHOTONSFUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][PHOTONSFDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][PUUP],false,lep1,lep2,lep3,lep4,true,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][PUDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][ISRUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][ISRDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][EWKUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                                FillerZZ(cr1Maps[controlregionZZ][LL][EWKDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                                //for(int i=0; i<nWeights; i++){
                                                                for(int i=config_pdfStart; i<config_pdfStart+config_pdfStep; i++) {
                                                                        //FillerZZ(cr1Maps[controlregionZZ][EE][PDFNAMES[i]],false,lep1,lep2,lep3,lep4,true,i);
                                                                        FillerZZ(cr1Maps[controlregionZZ][LL][PDFNAMES[i]],false,lep1,lep2,lep3,lep4,true,i);
                                                                }
                                                        }

                                                }

                                                if(*selMuonSize==2 && *selElectronSize==2) {

                                                        if((negMuons->size()==1)&&(posMuons->size()==1)&&(negElectrons->size()==1)&&(posElectrons->size()==1)) {

                                                                selMuon m1 = negMuons->at(0);
                                                                selMuon m2 = posMuons->at(0);
                                                                selElectron e1 = negElectrons->at(0);
                                                                selElectron e2 = posElectrons->at(0);

                                                                float m1m2 = abs((m1.vec+m2.vec).M()-ZMass);
                                                                float e1e2 = abs((e1.vec+e2.vec).M()-ZMass);

                                                                if(m1m2<e1e2) {
                                                                        if(m1.p.Pt()>m2.p.Pt()) {
                                                                                lep1.setAll(m1);
                                                                                lep2.setAll(m2);
                                                                        }else{
                                                                                lep1.setAll(m2);
                                                                                lep2.setAll(m1);
                                                                        }
                                                                        if(e1.p.Pt()>e2.p.Pt()) {
                                                                                lep3.setAll(e1);
                                                                                lep4.setAll(e2);
                                                                        }else{
                                                                                lep3.setAll(e2);
                                                                                lep4.setAll(e1);
                                                                        }

                                                                        float mll1=(lep1.vec+lep2.vec).M();
                                                                        float mll2=(lep3.vec+lep4.vec).M();
                                                                        if(mll1<106. && mll1>76. && mll2<130. && mll2>50.) {
                                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][nom],false,lep1,lep2,lep3,lep4,true);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][JESu],false,lep1,lep2,lep3,lep4,true,9999,JESUP);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][JESd],false,lep1,lep2,lep3,lep4,true,9999,JESDOWN);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][JERu],false,lep1,lep2,lep3,lep4,true,9999,JERDOWN);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][JERd],false,lep1,lep2,lep3,lep4,true,9999,JERDOWN);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][LEPSFUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][LEPSFDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][PHOTONSFUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][PHOTONSFDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][PUUP],false,lep1,lep2,lep3,lep4,true,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][PUDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][ISRUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][ISRDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][EWKUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][MM][EWKDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][nom],false,lep1,lep2,lep3,lep4,true);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][JESu],false,lep1,lep2,lep3,lep4,true,9999,JESUP);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][JESd],false,lep1,lep2,lep3,lep4,true,9999,JESDOWN);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][JERu],false,lep1,lep2,lep3,lep4,true,9999,JERDOWN);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][JERd],false,lep1,lep2,lep3,lep4,true,9999,JERDOWN);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][LEPSFUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][LEPSFDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][PHOTONSFUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][PHOTONSFDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][PUUP],false,lep1,lep2,lep3,lep4,true,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][PUDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][ISRUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][ISRDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][EWKUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][EWKDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                                                //for(int i=0; i<nWeights; i++){
                                                                                for(int i=config_pdfStart; i<config_pdfStart+config_pdfStep; i++) {
                                                                                        //FillerZZ(cr1Maps[controlregionZZ][MM][PDFNAMES[i]],false,lep1,lep2,lep3,lep4,true,i);
                                                                                        FillerZZ(cr1Maps[controlregionZZ][LL][PDFNAMES[i]],false,lep1,lep2,lep3,lep4,true,i);
                                                                                }
                                                                        }

                                                                }else{
                                                                        if(m1.p.Pt()>m2.p.Pt()) {
                                                                                lep3.setAll(m1);
                                                                                lep4.setAll(m2);
                                                                        }else{
                                                                                lep3.setAll(m2);
                                                                                lep4.setAll(m1);
                                                                        }
                                                                        if(e1.p.Pt()>e2.p.Pt()) {
                                                                                lep1.setAll(e1);
                                                                                lep2.setAll(e2);
                                                                        }else{
                                                                                lep1.setAll(e2);
                                                                                lep2.setAll(e1);
                                                                        }

                                                                        float mll1=(lep1.vec+lep2.vec).M();
                                                                        float mll2=(lep3.vec+lep4.vec).M();
                                                                        if(mll1<106. && mll1>76. && mll2<130. && mll2>50.) {
                                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][nom],false,lep1,lep2,lep3,lep4,true);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][JESu],false,lep1,lep2,lep3,lep4,true,9999,JESUP);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][JESd],false,lep1,lep2,lep3,lep4,true,9999,JESDOWN);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][JERu],false,lep1,lep2,lep3,lep4,true,9999,JERUP);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][JERd],false,lep1,lep2,lep3,lep4,true,9999,JERDOWN);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][LEPSFUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][LEPSFDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][PHOTONSFUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][PHOTONSFDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][PUUP],false,lep1,lep2,lep3,lep4,true,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][PUDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][ISRUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][ISRDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][EWKUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                                                //FillerZZ(cr1Maps[controlregionZZ][EE][EWKDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                                                //for(int i=0; i<nWeights; i++){
                                                                                for(int i=config_pdfStart; i<config_pdfStart+config_pdfStep; i++) {
                                                                                        //FillerZZ(cr1Maps[controlregionZZ][EE][PDFNAMES[i]],false,lep1,lep2,lep3,lep4,true,i);
                                                                                        FillerZZ(cr1Maps[controlregionZZ][LL][PDFNAMES[i]],false,lep1,lep2,lep3,lep4,true,i);
                                                                                }
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][nom],false,lep1,lep2,lep3,lep4,true);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][JESu],false,lep1,lep2,lep3,lep4,true,9999,JESUP);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][JESd],false,lep1,lep2,lep3,lep4,true,9999,JESDOWN);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][JERu],false,lep1,lep2,lep3,lep4,true,9999,JERUP);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][JERd],false,lep1,lep2,lep3,lep4,true,9999,JERDOWN);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][LEPSFUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][LEPSFDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][PHOTONSFUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][PHOTONSFDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][PUUP],false,lep1,lep2,lep3,lep4,true,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][PUDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][ISRUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][ISRDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][EWKUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                                                FillerZZ(cr1Maps[controlregionZZ][LL][EWKDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                                        }
                                                                }
                                                        }
                                                }
                                                //float mll1=(lep1.vec+lep2.vec).M();
                                                //float mll2=(lep3.vec+lep4.vec).M();
                                                //if(mll1<106. && mll1>76. && mll2<130. && mll2>50.){
                                                //FillerZZ(cr1Maps[controlregionZZ][LL][nom],false,lep1,lep2,lep3,lep4,true);
                                                //FillerZZ(cr1Maps[controlregionZZ][LL][JESu],false,lep1,lep2,lep3,lep4,true,9999,JESUP);
                                                //FillerZZ(cr1Maps[controlregionZZ][LL][JESd],false,lep1,lep2,lep3,lep4,true,9999,JESDOWN);
                                                //FillerZZ(cr1Maps[controlregionZZ][LL][JERu],false,lep1,lep2,lep3,lep4,true,9999,JERUP);
                                                //FillerZZ(cr1Maps[controlregionZZ][LL][JERd],false,lep1,lep2,lep3,lep4,true,9999,JERDOWN);
                                                //FillerZZ(cr1Maps[controlregionZZ][LL][LEPSFUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerZZ(cr1Maps[controlregionZZ][LL][LEPSFDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerZZ(cr1Maps[controlregionZZ][LL][PHOTONSFUP],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                //FillerZZ(cr1Maps[controlregionZZ][LL][PHOTONSFDOWN],false,lep1,lep2,lep3,lep4,true,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                //for(int i=0; i<nWeights; i++){
                                                //FillerZZ(cr1Maps[controlregionZZ][LL][PDFNAMES[i]],false,lep1,lep2,lep3,lep4,true,i);
                                                //}
                                                //}
                                        }
                                }
                        }
                }
        }


        clearCutFlowMap();
        if(config_selectionsToProcessMap[ControlRegionWZ]) {
                if(SelectEvent(ControlRegionZZ)) {
                        if((*selLeptonSize==3)&&(*matchedLeptonSize==2)) {
                                if((*selMuonSize==3)||(*selElectronSize==3)||(*selMuonSize==1 && *selElectronSize==2)||(*selMuonSize==2 && *selElectronSize==1)) {

                                        float ZMass = 91.1876;


                                        if((*countNegCharge==1 && *countPosCharge==2)||(*countNegCharge==2 && *countPosCharge==1)) {
                                                //find lepton combinations

                                                selLepton lep1,lep2,lep3;
                                                TLorentzVector temp_JESu(0.,0.,0.,0.);
                                                temp_JESu.SetPtEtaPhiM(met_JESu->p.Pt(),met_JESu->p.Eta(),met_JESu->p.Phi(),0.);
                                                TLorentzVector temp_JESd(0.,0.,0.,0.);
                                                temp_JESd.SetPtEtaPhiM(met_JESd->p.Pt(),met_JESd->p.Eta(),met_JESd->p.Phi(),0.);
                                                TLorentzVector temp_JERu(0.,0.,0.,0.);
                                                temp_JERu.SetPtEtaPhiM(met_JERu->p.Pt(),met_JERu->p.Eta(),met_JERu->p.Phi(),0.);
                                                TLorentzVector temp_JERd(0.,0.,0.,0.);
                                                temp_JERd.SetPtEtaPhiM(met_JERd->p.Pt(),met_JERd->p.Eta(),met_JERd->p.Phi(),0.);


                                                if(*selMuonSize==3) {

                                                        if((negMuons->size()==2)&&(posMuons->size()==1)) {

                                                                selMuon neg1 = negMuons->at(0);
                                                                selMuon neg2 = negMuons->at(1);
                                                                selMuon pos1 = posMuons->at(0);

                                                                //there are 2 combinations n1p1/n2 - n2p1/n1
                                                                float n1p1 = abs((neg1.vec+pos1.vec).M()-ZMass);
                                                                float n2p1 = abs((neg2.vec+pos1.vec).M()-ZMass);

                                                                if(n1p1<n2p1) {
                                                                        if(neg1.p.Pt()>pos1.p.Pt()) {
                                                                                lep1.setAll(neg1);
                                                                                lep2.setAll(pos1);
                                                                        }else{
                                                                                lep1.setAll(pos1);
                                                                                lep2.setAll(neg1);
                                                                        }
                                                                        lep3.setAll(neg2);
                                                                }else{
                                                                        if(n2p1<n1p1) {
                                                                                if(neg2.p.Pt()>pos1.p.Pt()) {
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
                                                        if((negMuons->size()==1)&&(posMuons->size()==2)) {

                                                                selMuon neg1 = negMuons->at(0);
                                                                selMuon pos1 = posMuons->at(0);
                                                                selMuon pos2 = posMuons->at(1);

                                                                //there are 2 combinations n1p1/p2 - n1p2/p1
                                                                float n1p1 = abs((neg1.vec+pos1.vec).M()-ZMass);
                                                                float n1p2 = abs((neg1.vec+pos2.vec).M()-ZMass);

                                                                if(n1p1<n1p2) {
                                                                        if(neg1.p.Pt()>pos1.p.Pt()) {
                                                                                lep1.setAll(neg1);
                                                                                lep2.setAll(pos1);
                                                                        }else{
                                                                                lep1.setAll(pos1);
                                                                                lep2.setAll(neg1);
                                                                        }
                                                                        lep3.setAll(pos2);
                                                                }else{
                                                                        if(n1p2<n1p1) {
                                                                                if(neg1.p.Pt()>pos2.p.Pt()) {
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
                                                        if(mll1<106. && mll1>76.) {
                                                                if( (lep1.vec+ETmiss_vec->vec).Mt() > 50. ) {
                                                                        if(*ETmiss>70.) {
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][nom],false,lep1,lep2,lep3,true);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][LEPSFUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][LEPSFDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][PHOTONSFUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][PHOTONSFDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][PUUP],false,lep1,lep2,lep3,true,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][PUDOWN],false,lep1,lep2,lep3,true,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][ISRUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][ISRDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][EWKUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][EWKDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][nom],false,lep1,lep2,lep3,true);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][LEPSFUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][LEPSFDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][PHOTONSFUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][PHOTONSFDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][PUUP],false,lep1,lep2,lep3,true,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][PUDOWN],false,lep1,lep2,lep3,true,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][ISRUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][ISRDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][EWKUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][EWKDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                                                //for(int i=0; i<nWeights; i++){
                                                                                for(int i=config_pdfStart; i<config_pdfStart+config_pdfStep; i++) {
                                                                                        //FillerWZ(cr1Maps[controlregionWZ][MM][PDFNAMES[i]],false,lep1,lep2,lep3,true,i);
                                                                                        FillerWZ(cr1Maps[controlregionWZ][LL][PDFNAMES[i]],false,lep1,lep2,lep3,true,i);
                                                                                }
                                                                        }
                                                                        if(met_JESu->p.Pt()>70.) {
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][JESu],false,lep1,lep2,lep3,true,9999,JESUP);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][JESu],false,lep1,lep2,lep3,true,9999,JESUP);
                                                                        }
                                                                        if(met_JESd->p.Pt()>70.) {
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][JESd],false,lep1,lep2,lep3,true,9999,JESDOWN);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][JESd],false,lep1,lep2,lep3,true,9999,JESDOWN);
                                                                        }
                                                                        if(met_JERu->p.Pt()>70.) {
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][JERu],false,lep1,lep2,lep3,true,9999,JERUP);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][JERu],false,lep1,lep2,lep3,true,9999,JERUP);
                                                                        }
                                                                        if(met_JERd->p.Pt()>70.) {
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][JERd],false,lep1,lep2,lep3,true,9999,JERDOWN);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][JERd],false,lep1,lep2,lep3,true,9999,JERDOWN);
                                                                        }
                                                                }
                                                        }

                                                }


                                                if(*selElectronSize==3) {

                                                        if((negElectrons->size()==2)&&(posElectrons->size()==1)) {

                                                                selElectron neg1 = negElectrons->at(0);
                                                                selElectron neg2 = negElectrons->at(1);
                                                                selElectron pos1 = posElectrons->at(0);

                                                                //there are 2 combinations n1p1/n2 - n2p1/n1
                                                                float n1p1 = abs((neg1.vec+pos1.vec).M()-ZMass);
                                                                float n2p1 = abs((neg2.vec+pos1.vec).M()-ZMass);

                                                                if(n1p1<n2p1) {
                                                                        if(neg1.p.Pt()>pos1.p.Pt()) {
                                                                                lep1.setAll(neg1);
                                                                                lep2.setAll(pos1);
                                                                        }else{
                                                                                lep1.setAll(pos1);
                                                                                lep2.setAll(neg1);
                                                                        }
                                                                        lep3.setAll(neg2);
                                                                }else{
                                                                        if(n2p1<n1p1) {
                                                                                if(neg2.p.Pt()>pos1.p.Pt()) {
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
                                                        if((negElectrons->size()==1)&&(posElectrons->size()==2)) {

                                                                selElectron neg1 = negElectrons->at(0);
                                                                selElectron pos1 = posElectrons->at(0);
                                                                selElectron pos2 = posElectrons->at(1);

                                                                //there are 2 combinations n1p1/p2 - n1p2/p1
                                                                float n1p1 = abs((neg1.vec+pos1.vec).M()-ZMass);
                                                                float n1p2 = abs((neg1.vec+pos2.vec).M()-ZMass);

                                                                if(n1p1<n1p2) {
                                                                        if(neg1.p.Pt()>pos1.p.Pt()) {
                                                                                lep1.setAll(neg1);
                                                                                lep2.setAll(pos1);
                                                                        }else{
                                                                                lep1.setAll(pos1);
                                                                                lep2.setAll(neg1);
                                                                        }
                                                                        lep3.setAll(pos2);
                                                                }else{
                                                                        if(n1p2<n1p1) {
                                                                                if(neg1.p.Pt()>pos2.p.Pt()) {
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
                                                        if(mll1<106. && mll1>76.) {
                                                                if( (lep1.vec+ETmiss_vec->vec).Mt() > 50. ) {
                                                                        if(*ETmiss>70.) {
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][nom],false,lep1,lep2,lep3,true);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][LEPSFUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][LEPSFDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][PHOTONSFUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][PHOTONSFDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][PUUP],false,lep1,lep2,lep3,true,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][PUDOWN],false,lep1,lep2,lep3,true,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][ISRUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][ISRDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][EWKUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][EWKDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][nom],false,lep1,lep2,lep3,true);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][LEPSFUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][LEPSFDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][PHOTONSFUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][PHOTONSFDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][PUUP],false,lep1,lep2,lep3,true,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][PUDOWN],false,lep1,lep2,lep3,true,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][ISRUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][ISRDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][EWKUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][EWKDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                                                //for(int i=0; i<nWeights; i++){
                                                                                for(int i=config_pdfStart; i<config_pdfStart+config_pdfStep; i++) {
                                                                                        //FillerWZ(cr1Maps[controlregionWZ][EE][PDFNAMES[i]],false,lep1,lep2,lep3,true,i);
                                                                                        FillerWZ(cr1Maps[controlregionWZ][LL][PDFNAMES[i]],false,lep1,lep2,lep3,true,i);
                                                                                }
                                                                        }
                                                                }
                                                        }
                                                        if( (lep1.vec+temp_JESu).Mt() > 50. ) {
                                                                if(met_JESu->p.Pt()>70.) {
                                                                        //FillerWZ(cr1Maps[controlregionWZ][EE][JESu],false,lep1,lep2,lep3,true,9999,JESUP);
                                                                        FillerWZ(cr1Maps[controlregionWZ][LL][JESu],false,lep1,lep2,lep3,true,9999,JESUP);
                                                                }
                                                        }
                                                        if( (lep1.vec+temp_JESd).Mt() > 50. ) {
                                                                if(met_JESd->p.Pt()>70.) {
                                                                        //FillerWZ(cr1Maps[controlregionWZ][EE][JESd],false,lep1,lep2,lep3,true,9999,JESDOWN);
                                                                        FillerWZ(cr1Maps[controlregionWZ][LL][JESd],false,lep1,lep2,lep3,true,9999,JESDOWN);
                                                                }
                                                        }
                                                        if( (lep1.vec+temp_JERu).Mt() > 50. ) {
                                                                if(met_JERu->p.Pt()>70.) {
                                                                        //FillerWZ(cr1Maps[controlregionWZ][EE][JERu],false,lep1,lep2,lep3,true,9999,JERUP);
                                                                        FillerWZ(cr1Maps[controlregionWZ][LL][JERu],false,lep1,lep2,lep3,true,9999,JERUP);
                                                                }
                                                        }
                                                        if( (lep1.vec+temp_JERd).Mt() > 50. ) {
                                                                if(met_JERd->p.Pt()>70.) {
                                                                        //FillerWZ(cr1Maps[controlregionWZ][EE][JERd],false,lep1,lep2,lep3,true,9999,JERDOWN);
                                                                        FillerWZ(cr1Maps[controlregionWZ][LL][JERd],false,lep1,lep2,lep3,true,9999,JERDOWN);
                                                                }
                                                        }


                                                }

                                                if((*selElectronSize==2)&&(*selMuonSize==1)) {

                                                        if(negElectrons->size()==1 && posElectrons->size()==1 && posMuons->size()==1) {
                                                                selElectron neg1 = negElectrons->at(0);
                                                                selElectron pos1 = posElectrons->at(0);
                                                                selMuon pos2 = posMuons->at(0);

                                                                if(neg1.p.Pt()>pos1.p.Pt()) {
                                                                        lep1.setAll(neg1);
                                                                        lep2.setAll(pos1);
                                                                }else{
                                                                        lep1.setAll(pos1);
                                                                        lep2.setAll(neg1);
                                                                }
                                                                lep3.setAll(pos2);
                                                        }
                                                        if(negElectrons->size()==1 && posElectrons->size()==1 && negMuons->size()==1) {
                                                                selElectron neg1 = negElectrons->at(0);
                                                                selElectron pos1 = posElectrons->at(0);
                                                                selMuon neg2 = negMuons->at(0);

                                                                if(neg1.p.Pt()>pos1.p.Pt()) {
                                                                        lep1.setAll(neg1);
                                                                        lep2.setAll(pos1);
                                                                }else{
                                                                        lep1.setAll(pos1);
                                                                        lep2.setAll(neg1);
                                                                }
                                                                lep3.setAll(neg2);
                                                        }

                                                        float mll1=(lep1.vec+lep2.vec).M();
                                                        if(mll1<106. && mll1>76.) {
                                                                if( (lep1.vec+ETmiss_vec->vec).Mt() > 50. ) {
                                                                        if(*ETmiss>70.) {
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][nom],false,lep1,lep2,lep3,true);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][LEPSFUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][LEPSFDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][PHOTONSFUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][PHOTONSFDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][PUUP],false,lep1,lep2,lep3,true,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][PUDOWN],false,lep1,lep2,lep3,true,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][ISRUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][ISRDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][EWKUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][EWKDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][nom],false,lep1,lep2,lep3,true);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][LEPSFUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][LEPSFDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][PHOTONSFUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][PHOTONSFDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][PUUP],false,lep1,lep2,lep3,true,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][PUDOWN],false,lep1,lep2,lep3,true,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][ISRUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][ISRDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][EWKUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][EWKDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                                                //for(int i=0; i<nWeights; i++){
                                                                                for(int i=config_pdfStart; i<config_pdfStart+config_pdfStep; i++) {
                                                                                        //FillerWZ(cr1Maps[controlregionWZ][EE][PDFNAMES[i]],false,lep1,lep2,lep3,true,i);
                                                                                        FillerWZ(cr1Maps[controlregionWZ][LL][PDFNAMES[i]],false,lep1,lep2,lep3,true,i);
                                                                                }
                                                                        }
                                                                }
                                                                if( (lep1.vec+temp_JESu).Mt() > 50. ) {
                                                                        if(met_JESu->p.Pt()>70.) {
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][JESu],false,lep1,lep2,lep3,true,9999,JESUP);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][JESu],false,lep1,lep2,lep3,true,9999,JESUP);
                                                                        }
                                                                }
                                                                if( (lep1.vec+temp_JESd).Mt() > 50. ) {
                                                                        if(met_JESd->p.Pt()>70.) {
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][JESd],false,lep1,lep2,lep3,true,9999,JESDOWN);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][JESd],false,lep1,lep2,lep3,true,9999,JESDOWN);
                                                                        }
                                                                }
                                                                if( (lep1.vec+temp_JERu).Mt() > 50. ) {
                                                                        if(met_JERu->p.Pt()>70.) {
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][JERu],false,lep1,lep2,lep3,true,9999,JERUP);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][JERu],false,lep1,lep2,lep3,true,9999,JERUP);
                                                                        }
                                                                }
                                                                if( (lep1.vec+temp_JERd).Mt() > 50. ) {
                                                                        if(met_JERd->p.Pt()>70.) {
                                                                                //FillerWZ(cr1Maps[controlregionWZ][EE][JERd],false,lep1,lep2,lep3,true,9999,JERDOWN);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][JERd],false,lep1,lep2,lep3,true,9999,JERDOWN);
                                                                        }
                                                                }
                                                        }

                                                }

                                                if((*selElectronSize==1)&&(*selMuonSize==2)) {

                                                        if(negMuons->size()==1 && posMuons->size()==1 && posElectrons->size()==1) {
                                                                selMuon neg1 = negMuons->at(0);
                                                                selMuon pos1 = posMuons->at(0);
                                                                selElectron pos2 = posElectrons->at(0);

                                                                if(neg1.p.Pt()>pos1.p.Pt()) {
                                                                        lep1.setAll(neg1);
                                                                        lep2.setAll(pos1);
                                                                }else{
                                                                        lep1.setAll(pos1);
                                                                        lep2.setAll(neg1);
                                                                }
                                                                lep3.setAll(pos2);
                                                        }
                                                        if(negMuons->size()==1 && posMuons->size()==1 && negElectrons->size()==1) {
                                                                selMuon neg1 = negMuons->at(0);
                                                                selMuon pos1 = posMuons->at(0);
                                                                selElectron neg2 = negElectrons->at(0);

                                                                if(neg1.p.Pt()>pos1.p.Pt()) {
                                                                        lep1.setAll(neg1);
                                                                        lep2.setAll(pos1);
                                                                }else{
                                                                        lep1.setAll(pos1);
                                                                        lep2.setAll(neg1);
                                                                }
                                                                lep3.setAll(neg2);
                                                        }

                                                        float mll1=(lep1.vec+lep2.vec).M();
                                                        if(mll1<106. && mll1>76.) {
                                                                if( (lep1.vec+ETmiss_vec->vec).Mt() > 50. ) {
                                                                        if(*ETmiss>70.) {
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][nom],false,lep1,lep2,lep3,true);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][LEPSFUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][LEPSFDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][PHOTONSFUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][PHOTONSFDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][PUUP],false,lep1,lep2,lep3,true,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][PUDOWN],false,lep1,lep2,lep3,true,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][ISRUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][ISRDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][EWKUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][EWKDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][nom],false,lep1,lep2,lep3,true);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][LEPSFUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][LEPSFDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][PHOTONSFUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][PHOTONSFDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][PUUP],false,lep1,lep2,lep3,true,9999,normal,upPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][PUDOWN],false,lep1,lep2,lep3,true,9999,normal,downPU,normalLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][ISRUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,upISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][ISRDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,downISR,normalEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][EWKUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,upEWK);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][EWKDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,normalPHOTONSF,normalISR,downEWK);
                                                                                //for(int i=0; i<nWeights; i++){
                                                                                for(int i=config_pdfStart; i<config_pdfStart+config_pdfStep; i++) {
                                                                                        //FillerWZ(cr1Maps[controlregionWZ][MM][PDFNAMES[i]],false,lep1,lep2,lep3,true,i);
                                                                                        FillerWZ(cr1Maps[controlregionWZ][LL][PDFNAMES[i]],false,lep1,lep2,lep3,true,i);
                                                                                }
                                                                        }
                                                                }
                                                                if( (lep1.vec+temp_JESu).Mt() > 50. ) {
                                                                        if(met_JESu->p.Pt()>70.) {
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][JESu],false,lep1,lep2,lep3,true,9999,JESUP);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][JESu],false,lep1,lep2,lep3,true,9999,JESUP);
                                                                        }
                                                                }
                                                                if( (lep1.vec+temp_JESd).Mt() > 50. ) {
                                                                        if(met_JESd->p.Pt()>70.) {
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][JESd],false,lep1,lep2,lep3,true,9999,JESDOWN);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][JESd],false,lep1,lep2,lep3,true,9999,JESDOWN);
                                                                        }
                                                                }
                                                                if( (lep1.vec+temp_JERu).Mt() > 50. ) {
                                                                        if(met_JERu->p.Pt()>70.) {
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][JERu],false,lep1,lep2,lep3,true,9999,JERUP);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][JERu],false,lep1,lep2,lep3,true,9999,JERUP);
                                                                        }
                                                                }
                                                                if( (lep1.vec+temp_JERd).Mt() > 50. ) {
                                                                        if(met_JERd->p.Pt()>70.) {
                                                                                //FillerWZ(cr1Maps[controlregionWZ][MM][JERd],false,lep1,lep2,lep3,true,9999,JERDOWN);
                                                                                FillerWZ(cr1Maps[controlregionWZ][LL][JERd],false,lep1,lep2,lep3,true,9999,JERDOWN);
                                                                        }
                                                                }
                                                        }
                                                }



                                                //float mll1=(lep1.vec+lep2.vec).M();
                                                //if(mll1<106. && mll1>76.){
                                                //if( (lep1.vec+ETmiss_vec->vec).Mt() > 50. ){
                                                //if(*ETmiss>70.){
                                                //FillerWZ(cr1Maps[controlregionWZ][LL][nom],false,lep1,lep2,lep3,true);
                                                //FillerWZ(cr1Maps[controlregionWZ][LL][LEPSFUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,upLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerWZ(cr1Maps[controlregionWZ][LL][LEPSFDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,downLEPSF,normalPHOTONSF,normalISR,normalEWK);
                                                //FillerWZ(cr1Maps[controlregionWZ][LL][PHOTONSFUP],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,upPHOTONSF,normalISR,normalEWK);
                                                //FillerWZ(cr1Maps[controlregionWZ][LL][PHOTONSFDOWN],false,lep1,lep2,lep3,true,9999,normal,normalPU,normalLEPSF,downPHOTONSF,normalISR,normalEWK);
                                                //for(int i=0; i<nWeights; i++){
                                                //FillerWZ(cr1Maps[controlregionWZ][LL][PDFNAMES[i]],false,lep1,lep2,lep3,true,i);
                                                //}
                                                //}
                                                //}
                                                //if( (lep1.vec+temp_JESu).Mt() > 50. ){
                                                //if(met_JESu->p.Pt()>70.){
                                                //FillerWZ(cr1Maps[controlregionWZ][LL][JESu],false,lep1,lep2,lep3,true,9999,JESUP);
                                                //}
                                                //}
                                                //if( (lep1.vec+temp_JESd).Mt() > 50. ){
                                                //if(met_JESd->p.Pt()>70.){
                                                //FillerWZ(cr1Maps[controlregionWZ][LL][JESd],false,lep1,lep2,lep3,true,9999,JESDOWN);
                                                //}
                                                //}
                                                //if( (lep1.vec+temp_JERu).Mt() > 50. ){
                                                //if(met_JERu->p.Pt()>70.){
                                                //FillerWZ(cr1Maps[controlregionWZ][LL][JERu],false,lep1,lep2,lep3,true,9999,JERUP);
                                                //}
                                                //}
                                                //if( (lep1.vec+temp_JERd).Mt() > 50. ){
                                                //if(met_JERd->p.Pt()>70.){
                                                //FillerWZ(cr1Maps[controlregionWZ][LL][JERd],false,lep1,lep2,lep3,true,9999,JERDOWN);
                                                //}
                                                //}
                                                //}
                                        }
                                }

                        }
                }
        }


        clearCutFlowMap();
        if(config_selectionsToProcessMap[ONZ]) {
                if(SelectEvent(ONZ)) {
                        if(*selLeptonSize==*matchedLeptonSize && *selLeptonSize==2) {
                                if ((*selPhotonSize!=0)&&(*mll>81. && *mll<101.)) {
                                        if(*isDiMuon || *isDiElectron) {
                                                Filler(h1Maps[onz][LL],true);
                                                Filler2D(h2Maps[onz][LL],true);
                                        }
                                        if (*isDiElectron) {
                                                Filler(h1Maps[onz][EE],true);
                                        }
                                        //}else{
                                        if (*isDiMuon) {
                                                Filler(h1Maps[onz][MM],true);
                                        }
                                        if (*isElectronMuon || *isMuonElectron) {
                                                Filler(h1Maps[onz][EM],true);
                                        }

                                        //}
                                }
                        }
                }
        }

}


void myAnalyzer::Filler(map<Histograms1D,TH1F>& m,bool withPhoton,bool slimmed,int changePDF,changemet changeMET,changepu changePU,changeLEPSF changeLepSF,changePHOTONSF changePhotonSF,changeISR changeisr,changeEWK changeewk){

        //tempWeight = tempWeight*weight_LepSF *weight_PhotonSF *weight_TopPt*weight_nIsr*weight_EWKinoPt*weight_PDF;
        float tempWeight=*mc_weight;
        float weight_PDF=1.;
        if(changePDF>150) {
                weight_PDF=1.;
        }else{
                weight_PDF=pdf_weights->at(changePDF);
        }
        //tempWeight = tempWeight*puWeights[(int)changePU]*lepSfWeights[(int)changeLepSF] *photonSfWeights[(int)changePhotonSF] **topPt_weight*isrWeights[(int)changeisr]*ewkWeights[(int)changeewk]*weight_PDF;
        tempWeight = tempWeight*puWeights[(int)changePU]*lepSfWeights[(int)changeLepSF] *photonSfWeights[(int)changePhotonSF] **topPt_weight*isrWeights[(int)changeisr]*ewkWeights[(int)changeewk]*weight_PDF*zzKFactor;
        if(!isData) {
                if(!isSignal) {
                        if(*isDiElectron) tempWeight=tempWeight*config_eeWeight;
                        if(*isMuonElectron||*isElectronMuon) tempWeight=tempWeight*config_emWeight;
                        if(*isDiMuon) tempWeight=tempWeight*config_mmWeight;
                }else{
                        if(*isDiElectron) tempWeight=tempWeight*config_eeEff;
                        if(*isMuonElectron||*isElectronMuon) tempWeight=tempWeight*config_emEff;
                        if(*isDiMuon) tempWeight=tempWeight*config_mmEff;
                }
        }

        float TChiNG_weight=1.;
        bool useTChiNG=true;
        bool isTChiNG=false;
        if(inputName.find("TChiNG") != string::npos) {
                isTChiNG=true;
                //cout<<"here"<<endl;
                //cout<<(*intermediateGenParticles)[0].daughters.size()<<(*intermediateGenParticles)[1].daughters.size()<<endl;
                //cout<<(*intermediateGenParticles)[0].daughters[0].pdgId<<(*intermediateGenParticles)[0].daughters[1].pdgId<<endl;
                //cout<<(*intermediateGenParticles)[1].daughters[0].pdgId<<(*intermediateGenParticles)[1].daughters[1].pdgId<<endl;
                //GG case -> weight=1. 0.25to0.25
                bool oneOfThem=false;
                if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2) {
                        if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 22) {
                                //cout<<"buh"<<endl;
                                if (fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 22) {
                                        //cout<<"sdfdfsd"<<endl;
                                        useTChiNG=true;
                                        TChiNG_weight=1.;
                                        oneOfThem=oneOfThem||true;

                                }
                        }
                }
                //ZZ case ->weight=4. 1/16 to 0.25
                if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2) {
                        if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 23) {
                                if (fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 23) {
                                        useTChiNG=true;
                                        TChiNG_weight=4.;
                                        oneOfThem=oneOfThem||true;

                                }
                        }
                }
                //cout<<"test"<<endl;
                //ZG case -> weight = 2 0.25to0.5
                if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2) {
                        //cout<<"test 2"<<endl;
                        if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 23) {
                                //cout<<"test3"<<endl;
                                if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 22) {
                                        useTChiNG=true;
                                        TChiNG_weight=2.;
                                        oneOfThem=oneOfThem||true;
                                }
                        }
                }
                if(oneOfThem==false) {
                        useTChiNG=false;
                        TChiNG_weight=1.;
                }



        }

        bool usePoint=(isTChiNG&&useTChiNG) || (!isTChiNG);
        if(usePoint) {
                if(isTChiNG) {
                        tempWeight=tempWeight*TChiNG_weight;
                }

                //m.at(WEIGHT_EWKINOPAIRPT).Fill(ewkWeights[(int)changeewk]);
                //m.at(WEIGHT_NISR).Fill(isrWeights[(int)changeisr]);
                //m.at(WEIGHT_TOPPT).Fill(*topPt_weight);
                //m.at(WEIGHT_PDF).Fill(weight_PDF);

                if(changeMET==normal) {
                        m.at(ETMISS).Fill(*ETmiss, tempWeight);
                }
                else{
                        if(changeMET==JESUP) {
                                m.at(ETMISS).Fill(met_JESu->p.Pt(),tempWeight);
                        }else{
                                if(changeMET==JESDOWN) {
                                        m.at(ETMISS).Fill(met_JESd->p.Pt(),tempWeight);
                                }else{
                                        if(changeMET==JERUP) {
                                                m.at(ETMISS).Fill(met_JERu->p.Pt(),tempWeight);
                                        }else{
                                                if(changeMET==JERDOWN) {
                                                        m.at(ETMISS).Fill(met_JERd->p.Pt(),tempWeight);
                                                }
                                        }
                                }
                        }
                }


                m.at(PT1).Fill(*pt1, tempWeight);
                m.at(PT2).Fill(*pt2, tempWeight);
                m.at(MLL).Fill(*mll,tempWeight);
                m.at(NPHOTONS).Fill(*selPhotonSize,tempWeight);
                m.at(HT).Fill(*calcHt,tempWeight);
                m.at(NJETS).Fill(*selJetSize,tempWeight);
                m.at(ETA1).Fill(fabs(*eta1),tempWeight);
                m.at(ETA2).Fill(fabs(*eta2),tempWeight);
                m.at(PHI1).Fill(fabs(*phi2),tempWeight);
                m.at(PHI2).Fill(fabs(*phi2),tempWeight);

                m.at(DeltaPhiLL).Fill(fabs(l1->vec.DeltaPhi(l2->vec)),tempWeight);

                m.at(MT2).Fill(*Mt2,tempWeight);

                if(!slimmed) {
                        if(*selJetSize>0) {
                                m.at(JetPt1).Fill(selJets->at(0).p.Pt(),tempWeight);
                                m.at(JetPhi1).Fill(selJets->at(0).p.Phi(),tempWeight);
                                m.at(JetEta1).Fill(selJets->at(0).p.Eta(),tempWeight);
                                if(*selJetSize>1) {
                                        m.at(JetPt2).Fill(selJets->at(1).p.Pt(),tempWeight);
                                        m.at(JetPhi2).Fill(selJets->at(1).p.Phi(),tempWeight);
                                        m.at(JetEta2).Fill(selJets->at(1).p.Eta(),tempWeight);
                                        if(*selJetSize>2) {
                                                m.at(JetPt3).Fill(selJets->at(2).p.Pt(),tempWeight);
                                                m.at(JetPhi3).Fill(selJets->at(2).p.Phi(),tempWeight);
                                                m.at(JetEta3).Fill(selJets->at(2).p.Eta(),tempWeight);
                                                if(*selJetSize>3) {
                                                        m.at(JetPt4).Fill(selJets->at(3).p.Pt(),tempWeight);
                                                        m.at(JetPhi4).Fill(selJets->at(3).p.Phi(),tempWeight);
                                                        m.at(JetEta4).Fill(selJets->at(3).p.Eta(),tempWeight);
                                                }
                                        }
                                }
                        }

                        m.at(NElectrons).Fill(*selElectronSize,tempWeight);
                        m.at(NMuons).Fill(*selMuonSize,tempWeight);
                        //m.at(NVTX).Fill(*nGoodVertices,tempWeight);
                        m.at(NVTX).Fill(*nTruePV,tempWeight);
                        m.at(GENHT).Fill(*genHt,tempWeight);
                        m.at(NBJETS).Fill(*selBJetSize,tempWeight);
                        m.at(DeltaEtaLL).Fill(fabs(*eta1-*eta2),tempWeight);
                        //m.at(DeltaPhiLL).Fill(fabs(l1->vec.DeltaPhi(l2->vec)),tempWeight);
                        m.at(DeltaRLL).Fill(fabs(*deltaRll),tempWeight);
                        m.at(DeltaEtaLL_neg).Fill(*eta1-*eta2,tempWeight);
                        m.at(DeltaPhiLL_neg).Fill(l1->vec.DeltaPhi(l2->vec),tempWeight);
                        m.at(DeltaRLL_neg).Fill(*deltaRll,tempWeight);
                        m.at(ZPT).Fill((l1->vec+l2->vec).Pt(),tempWeight);
                        m.at(MTLL).Fill((l1->vec+l2->vec).Mt(),tempWeight);
                        m.at(ST).Fill(*pt1+*pt2,tempWeight);
                        m.at(DeltaPhiLLMet).Fill(fabs((l1->vec+l2->vec).DeltaPhi(ETmiss_vec->vec)),tempWeight);
                        m.at(DeltaEtaLLMet).Fill(fabs((l1->vec+l2->vec).Eta() - ETmiss_vec->vec.Eta()),tempWeight);
                        m.at(DeltaRLLMet).Fill(fabs((l1->vec+l2->vec).DeltaR(ETmiss_vec->vec)),tempWeight);
                }


                //if(ev.evtHasGenPhotonVeto) m.at(VetoCompare).Fill(1.,tempWeight);
                //if(!ev.evtHasGenPhotonVeto) m.at(VetoCompare).Fill(0.,tempWeight);

                if(withPhoton) {
                        m.at(PTG1).Fill(selPhotons->at(0).p.Pt(),tempWeight);
                        m.at(PHIG1).Fill(selPhotons->at(0).p.Phi(),tempWeight);
                        m.at(ETAG1).Fill(selPhotons->at(0).p.Eta(),tempWeight);

                        if(selPhotons->at(0).matchedToPhoton) {
                                m.at(Fakes).Fill(1.,tempWeight);
                        }
                        if(selPhotons->at(0).matchedToElectron) {
                                m.at(Fakes).Fill(2.,tempWeight);
                        }
                        if(selPhotons->at(0).matchedToJet) {
                                m.at(Fakes).Fill(3.,tempWeight);
                        }
                        if(!((selPhotons->at(0).matchedToPhoton)||(selPhotons->at(0).matchedToElectron)||(selPhotons->at(0).matchedToJet))) {
                                m.at(Fakes).Fill(4.,tempWeight);
                        }
                        m.at(Fakes).Fill(0.,tempWeight);

                        //m.at(DELTARGL1).Fill(selPhotons->at(0).vec.DeltaR(l1->vec),tempWeight);
                        //m.at(DELTARGL2).Fill(selPhotons->at(0).vec.DeltaR(l2->vec),tempWeight);
                        m.at(DeltaPhiL1G).Fill(selPhotons->at(0).vec.DeltaPhi(l1->vec),tempWeight);
                        m.at(DeltaPhiL2G).Fill(selPhotons->at(0).vec.DeltaPhi(l2->vec),tempWeight);
                        m.at(DeltaPhiLLG).Fill(fabs((l1->vec+l2->vec).DeltaPhi(selPhotons->at(0).vec)),tempWeight);

                        m.at(ML1G).Fill((l1->vec+selPhotons->at(0).vec).M(),tempWeight);
                        m.at(ML2G).Fill((l2->vec+selPhotons->at(0).vec).M(),tempWeight);
                        m.at(MLLG).Fill((l1->vec+l2->vec+selPhotons->at(0).vec).M(),tempWeight);


                        if(!slimmed) {
                                //m.at(SIGMAIETAIETAG1).Fill(selPhotons->at(0).sigmaIetaIeta,tempWeight);
                                //m.at(SIGMAIPHIIPHIG1).Fill(selPhotons->at(0).sigmaIphiIphi,tempWeight);
                                //m.at(R9).Fill(selPhotons->at(0).r9,tempWeight);
                                //m.at(HOVERE).Fill(selPhotons->at(0).hOverE,tempWeight);
                                //m.at(DELTARGL1).Fill(selPhotons->at(0).deltaR1,tempWeight);
                                //m.at(DELTARGL2).Fill(selPhotons->at(0).deltaR2,tempWeight);
                                m.at(DeltaEtaLLG).Fill(fabs((l1->vec+l2->vec).Eta()-selPhotons->at(0).vec.Eta()), tempWeight);
                                //m.at(DeltaPhiLLG).Fill(fabs((l1->vec+l2->vec).DeltaPhi(selPhotons->at(0).vec)),tempWeight);
                                m.at(DeltaRLLG).Fill(fabs((l1->vec+l2->vec).DeltaR(selPhotons->at(0).vec)),tempWeight);
                                m.at(DeltaPhiGMet).Fill(fabs((ETmiss_vec->vec).DeltaPhi(selPhotons->at(0).vec)),tempWeight);
                                m.at(DeltaRGMet).Fill(fabs((ETmiss_vec->vec).DeltaR(selPhotons->at(0).vec)),tempWeight);
                                m.at(MTLLG).Fill((l1->vec+l2->vec+selPhotons->at(0).vec).Mt(),tempWeight);
                                m.at(MTL1MET).Fill((l1->vec+ETmiss_vec->vec).Mt(),tempWeight);
                                m.at(MTL2MET).Fill((l2->vec+ETmiss_vec->vec).Mt(),tempWeight);
                                m.at(MTGMET).Fill((selPhotons->at(0).vec+ETmiss_vec->vec).Mt(),tempWeight);
                                m.at(MTLLMET).Fill((l1->vec+l2->vec+ETmiss_vec->vec).Mt(),tempWeight);
                                m.at(MTLLGMET).Fill((l1->vec+l2->vec+selPhotons->at(0).vec+ETmiss_vec->vec).Mt(),tempWeight);
                                m.at(STG).Fill(*pt1+*pt2+selPhotons->at(0).p.Pt(),tempWeight);
                                m.at(STMET).Fill(*pt1+*pt2+selPhotons->at(0).p.Pt()+*ETmiss,tempWeight);
                                //m.at(MLLG).Fill((l1->vec+l2->vec+selPhotons->at(0).vec).M(),tempWeight);
                        }



                        //m.at(PT_llg).Fill((l1->vec+l2->vec+selPhotons->at(0).vec).Pt(),tempWeight);

                        //if(FindGenPhotonMatch(ev.selPhotons->at(0))) m.at(genPhotonPT).Fill(GetGenPhotonMatch(ev.selPhotons->at(0)).p.Pt(),tempWeight);

                        //if(ev.evtHasGenPhotonVeto) if(FindGenPhotonMatch(ev.selPhotons->at(0))) m.at(genPhotonPT_Veto).Fill(GetGenPhotonMatch(ev.selPhotons->at(0)).p.Pt(),tempWeight);
                        //if(ev.evtHasGenPhotonVeto) m.at(PTG1_Veto).Fill(ev.selPhotons->at(0).p.Pt(),tempWeight);
                        //if(!ev.evtHasGenPhotonVeto) if(FindGenPhotonMatch(ev.selPhotons->at(0))) m.at(genPhotonPT_NoVeto).Fill(GetGenPhotonMatch(ev.selPhotons->at(0)).p.Pt(),tempWeight);
                        //if(!ev.evtHasGenPhotonVeto) m.at(PTG1_NoVeto).Fill(ev.selPhotons->at(0).p.Pt(),tempWeight);
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
}
void myAnalyzer::FillerZZ(map<Histograms1D,TH1F>& m,bool withPhoton,selLepton& le1, selLepton& le2, selLepton& le3, selLepton& le4,bool slimmed,int changePDF,changemet changeMET,changepu changePU,changeLEPSF changeLepSF,changePHOTONSF changePhotonSF,changeISR changeisr,changeEWK changeewk){

        float tempWeight=*mc_weight;
        float weight_PDF=1.;
        if(changePDF>150) {
                weight_PDF=1.;
        }else{
                weight_PDF=pdf_weights->at(changePDF);
        }
        //tempWeight = tempWeight*puWeights[(int)changePU]*lepSfWeights[(int)changeLepSF] *photonSfWeights[(int)changePhotonSF] **topPt_weight*isrWeights[(int)changeisr]*ewkWeights[(int)changeewk]*weight_PDF;
        tempWeight = tempWeight*puWeights[(int)changePU]*lepSfWeights[(int)changeLepSF] *photonSfWeights[(int)changePhotonSF] **topPt_weight*isrWeights[(int)changeisr]*ewkWeights[(int)changeewk]*weight_PDF*zzKFactor;
        if(!isData) {
                if(!isSignal) {
                        if(*isDiElectron) tempWeight=tempWeight*config_eeWeight;
                        if(*isMuonElectron||*isElectronMuon) tempWeight=tempWeight*config_emWeight;
                        if(*isDiMuon) tempWeight=tempWeight*config_mmWeight;
                }else{
                        if(*isDiElectron) tempWeight=tempWeight*config_eeEff;
                        if(*isMuonElectron||*isElectronMuon) tempWeight=tempWeight*config_emEff;
                        if(*isDiMuon) tempWeight=tempWeight*config_mmEff;
                }
        }



        float TChiNG_weight=1.;
        bool useTChiNG=true;
        bool isTChiNG=false;
        if(inputName.find("TChiNG") != string::npos) {
                isTChiNG=true;
                //cout<<"here"<<endl;
                //cout<<(*intermediateGenParticles)[0].daughters.size()<<(*intermediateGenParticles)[1].daughters.size()<<endl;
                //cout<<(*intermediateGenParticles)[0].daughters[0].pdgId<<(*intermediateGenParticles)[0].daughters[1].pdgId<<endl;
                //cout<<(*intermediateGenParticles)[1].daughters[0].pdgId<<(*intermediateGenParticles)[1].daughters[1].pdgId<<endl;
                //GG case -> weight=1. 0.25to0.25
                bool oneOfThem=false;
                if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2) {
                        if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 22) {
                                //cout<<"buh"<<endl;
                                if (fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 22) {
                                        //cout<<"sdfdfsd"<<endl;
                                        useTChiNG=true;
                                        TChiNG_weight=1.;
                                        oneOfThem=oneOfThem||true;

                                }
                        }
                }
                //ZZ case ->weight=4. 1/16 to 0.25
                if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2) {
                        if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 23) {
                                if (fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 23) {
                                        useTChiNG=true;
                                        TChiNG_weight=4.;
                                        oneOfThem=oneOfThem||true;

                                }
                        }
                }
                //cout<<"test"<<endl;
                //ZG case -> weight = 2 0.25to0.5
                if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2) {
                        //cout<<"test 2"<<endl;
                        if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 23) {
                                //cout<<"test3"<<endl;
                                if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 22) {
                                        useTChiNG=true;
                                        TChiNG_weight=2.;
                                        oneOfThem=oneOfThem||true;
                                }
                        }
                }
                if(oneOfThem==false) {
                        useTChiNG=false;
                        TChiNG_weight=1.;
                }



        }

        bool usePoint=(isTChiNG&&useTChiNG) || (!isTChiNG);
        if(usePoint) {
                if(isTChiNG) {
                        tempWeight=tempWeight*TChiNG_weight;
                }
                //m.at(WEIGHT_EWKINOPAIRPT).Fill(ewkWeights[(int)changeewk]);
                //m.at(WEIGHT_NISR).Fill(isrWeights[(int)changeisr]);
                //m.at(WEIGHT_TOPPT).Fill(*topPt_weight);
                //m.at(WEIGHT_PDF).Fill(weight_PDF);


                if(changeMET==normal) {
                        m.at(ETMISS).Fill(*ETmiss, tempWeight);
                }
                else{
                        if(changeMET==JESUP) {
                                m.at(ETMISS).Fill(met_JESu->p.Pt(),tempWeight);
                        }else{
                                if(changeMET==JESDOWN) {
                                        m.at(ETMISS).Fill(met_JESd->p.Pt(),tempWeight);
                                }else{
                                        if(changeMET==JERUP) {
                                                m.at(ETMISS).Fill(met_JERu->p.Pt(),tempWeight);
                                        }else{
                                                if(changeMET==JERDOWN) {
                                                        m.at(ETMISS).Fill(met_JERd->p.Pt(),tempWeight);
                                                }
                                        }
                                }
                        }
                }


                m.at(PT1).Fill(le1.p.Pt(), tempWeight);
                m.at(PT2).Fill(le2.p.Pt(), tempWeight);
                m.at(PT3).Fill(le3.p.Pt(), tempWeight);
                m.at(PT4).Fill(le4.p.Pt(), tempWeight);
                m.at(MLL).Fill((le1.vec+le2.vec).M(),tempWeight);
                m.at(MLL2).Fill((le3.vec+le4.vec).M(),tempWeight);
                m.at(NPHOTONS).Fill(*selPhotonSize,tempWeight);
                m.at(HT).Fill(*calcHt,tempWeight);
                m.at(NJETS).Fill(*selJetSize,tempWeight);
                m.at(ETA1).Fill(fabs(le1.p.Eta()),tempWeight);
                m.at(ETA2).Fill(fabs(le2.p.Eta()),tempWeight);
                m.at(ETA3).Fill(fabs(le3.p.Eta()),tempWeight);
                m.at(ETA4).Fill(fabs(le4.p.Eta()),tempWeight);
                m.at(PHI1).Fill(fabs(le1.p.Phi()),tempWeight);
                m.at(PHI2).Fill(fabs(le2.p.Phi()),tempWeight);
                m.at(PHI3).Fill(fabs(le3.p.Phi()),tempWeight);
                m.at(PHI4).Fill(fabs(le4.p.Phi()),tempWeight);


                m.at(MT2).Fill(*Mt2,tempWeight);


                if(!slimmed) {
                        if(*selJetSize>0) {
                                m.at(JetPt1).Fill(selJets->at(0).p.Pt(),tempWeight);
                                m.at(JetPhi1).Fill(selJets->at(0).p.Phi(),tempWeight);
                                m.at(JetEta1).Fill(selJets->at(0).p.Eta(),tempWeight);
                                if(*selJetSize>1) {
                                        m.at(JetPt2).Fill(selJets->at(1).p.Pt(),tempWeight);
                                        m.at(JetPhi2).Fill(selJets->at(1).p.Phi(),tempWeight);
                                        m.at(JetEta2).Fill(selJets->at(1).p.Eta(),tempWeight);
                                        if(*selJetSize>2) {
                                                m.at(JetPt3).Fill(selJets->at(2).p.Pt(),tempWeight);
                                                m.at(JetPhi3).Fill(selJets->at(2).p.Phi(),tempWeight);
                                                m.at(JetEta3).Fill(selJets->at(2).p.Eta(),tempWeight);
                                                if(*selJetSize>3) {
                                                        m.at(JetPt4).Fill(selJets->at(3).p.Pt(),tempWeight);
                                                        m.at(JetPhi4).Fill(selJets->at(3).p.Phi(),tempWeight);
                                                        m.at(JetEta4).Fill(selJets->at(3).p.Eta(),tempWeight);
                                                }
                                        }
                                }
                        }
                        m.at(NVTX).Fill(*nGoodVertices,tempWeight);
                        m.at(MLLLL).Fill((le1.vec+le2.vec+le3.vec+le4.vec).M(),tempWeight);
                        m.at(GENHT).Fill(*genHt,tempWeight);
                        m.at(ZPT).Fill((le1.vec+le2.vec).Pt(),tempWeight);
                        m.at(ZPT2).Fill((le3.vec+le4.vec).Pt(),tempWeight);
                        m.at(ST).Fill(le1.p.Pt()+le2.p.Pt(),tempWeight);
                        m.at(NElectrons).Fill(*selElectronSize,tempWeight);
                        m.at(NMuons).Fill(*selMuonSize,tempWeight);
                }


                if(withPhoton) {
                        m.at(PTG1).Fill(selPhotons->at(0).p.Pt(),tempWeight);
                        m.at(PHIG1).Fill(selPhotons->at(0).p.Phi(),tempWeight);
                        m.at(ETAG1).Fill(selPhotons->at(0).p.Eta(),tempWeight);


                        if(selPhotons->at(0).matchedToPhoton) {
                                m.at(Fakes).Fill(1.,tempWeight);
                        }
                        if(selPhotons->at(0).matchedToElectron) {
                                m.at(Fakes).Fill(2.,tempWeight);
                        }
                        if(selPhotons->at(0).matchedToJet) {
                                m.at(Fakes).Fill(3.,tempWeight);
                        }
                        if(!((selPhotons->at(0).matchedToPhoton)||(selPhotons->at(0).matchedToElectron)||(selPhotons->at(0).matchedToJet))) {
                                m.at(Fakes).Fill(4.,tempWeight);
                        }
                        m.at(Fakes).Fill(0.,tempWeight);



                        //m.at(DELTARGL1).Fill(selPhotons->at(0).vec.DeltaR(l1->vec),tempWeight);
                        //m.at(DELTARGL2).Fill(selPhotons->at(0).vec.DeltaR(l2->vec),tempWeight);
                        m.at(DeltaPhiL1G).Fill(selPhotons->at(0).vec.DeltaPhi(l1->vec),tempWeight);
                        m.at(DeltaPhiL2G).Fill(selPhotons->at(0).vec.DeltaPhi(l2->vec),tempWeight);
                        m.at(DeltaPhiLLG).Fill(fabs((l1->vec+l2->vec).DeltaPhi(selPhotons->at(0).vec)),tempWeight);


                        //if(FindGenPhotonMatch(ev.selPhotons->at(0))) m.at(genPhotonPT).Fill(GetGenPhotonMatch(ev.selPhotons->at(0)).p.Pt(),tempWeight);
                        //
                        //if(ev.evtHasGenPhotonVeto) if(FindGenPhotonMatch(ev.selPhotons->at(0))) m.at(genPhotonPT_Veto).Fill(GetGenPhotonMatch(ev.selPhotons->at(0)).p.Pt(),tempWeight);
                        //if(ev.evtHasGenPhotonVeto) m.at(PTG1_Veto).Fill(ev.selPhotons->at(0).p.Pt(),tempWeight);
                        //if(!ev.evtHasGenPhotonVeto) if(FindGenPhotonMatch(ev.selPhotons->at(0))) m.at(genPhotonPT_NoVeto).Fill(GetGenPhotonMatch(ev.selPhotons->at(0)).p.Pt(),tempWeight);
                        //if(!ev.evtHasGenPhotonVeto) m.at(PTG1_NoVeto).Fill(ev.selPhotons->at(0).p.Pt(),tempWeight);


                        if(!slimmed) {
                                //m.at(SIGMAIETAIETAG1).Fill(selPhotons->at(0).sigmaIetaIeta,tempWeight);
                                //m.at(SIGMAIPHIIPHIG1).Fill(selPhotons->at(0).sigmaIphiIphi,tempWeight);
                                //m.at(R9).Fill(selPhotons->at(0).r9,tempWeight);
                                //m.at(HOVERE).Fill(selPhotons->at(0).hOverE,tempWeight);
                                m.at(DELTARGL1).Fill(selPhotons->at(0).deltaR1,tempWeight);
                                m.at(DELTARGL2).Fill(selPhotons->at(0).deltaR2,tempWeight);
                                m.at(DeltaEtaLLG).Fill(fabs((le1.vec+le2.vec).Eta()-selPhotons->at(0).vec.Eta()), tempWeight);
                                //m.at(DeltaPhiLLG).Fill(fabs((le1.vec+le2.vec).Phi()-selPhotons->at(0).vec.Phi()),tempWeight);
                                m.at(DeltaRLLG).Fill(fabs((le1.vec+le2.vec).DeltaR(selPhotons->at(0).vec)),tempWeight);
                                m.at(MTLLG).Fill((le1.vec+le2.vec+selPhotons->at(0).vec).Mt(),tempWeight);
                                m.at(MTL1MET).Fill((le1.vec+ETmiss_vec->vec).Mt(),tempWeight);
                                m.at(MTL2MET).Fill((le2.vec+ETmiss_vec->vec).Mt(),tempWeight);
                                m.at(MTGMET).Fill((selPhotons->at(0).vec+ETmiss_vec->vec).Mt(),tempWeight);
                                m.at(MTLLMET).Fill((le1.vec+le2.vec+ETmiss_vec->vec).Mt(),tempWeight);
                                m.at(MTLLGMET).Fill((le1.vec+le2.vec+selPhotons->at(0).vec+ETmiss_vec->vec).Mt(),tempWeight);
                                m.at(STG).Fill(le1.vec.Pt()+le2.vec.Pt()+selPhotons->at(0).p.Pt(),tempWeight);
                                m.at(STMET).Fill(le1.vec.Pt()+le2.vec.Pt()+selPhotons->at(0).p.Pt()+*ETmiss,tempWeight);
                                m.at(MLLG).Fill((le1.vec+le2.vec+selPhotons->at(0).vec).M(),tempWeight);
                        }

                }
        }
}

void myAnalyzer::FillerWZ(map<Histograms1D,TH1F>& m,bool withPhoton,selLepton& le1, selLepton& le2, selLepton& le3,bool slimmed,int changePDF,changemet changeMET,changepu changePU,changeLEPSF changeLepSF,changePHOTONSF changePhotonSF,changeISR changeisr,changeEWK changeewk){

        float tempWeight=*mc_weight;
        float weight_PDF=1.;
        if(changePDF>150) {
                weight_PDF=1.;
        }else{
                weight_PDF=pdf_weights->at(changePDF);
        }
        //tempWeight = tempWeight*puWeights[(int)changePU]*lepSfWeights[(int)changeLepSF] *photonSfWeights[(int)changePhotonSF] **topPt_weight*isrWeights[(int)changeisr]*ewkWeights[(int)changeewk]*weight_PDF;
        tempWeight = tempWeight*puWeights[(int)changePU]*lepSfWeights[(int)changeLepSF] *photonSfWeights[(int)changePhotonSF] **topPt_weight*isrWeights[(int)changeisr]*ewkWeights[(int)changeewk]*weight_PDF*zzKFactor;
        if(!isData) {
                if(!isSignal) {
                        if(*isDiElectron) tempWeight=tempWeight*config_eeWeight;
                        if(*isMuonElectron||*isElectronMuon) tempWeight=tempWeight*config_emWeight;
                        if(*isDiMuon) tempWeight=tempWeight*config_mmWeight;
                }else{
                        if(*isDiElectron) tempWeight=tempWeight*config_eeEff;
                        if(*isMuonElectron||*isElectronMuon) tempWeight=tempWeight*config_emEff;
                        if(*isDiMuon) tempWeight=tempWeight*config_mmEff;
                }
        }


        float TChiNG_weight=1.;
        bool useTChiNG=true;
        bool isTChiNG=false;
        if(inputName.find("TChiNG") != string::npos) {
                isTChiNG=true;
                //cout<<"here"<<endl;
                //cout<<(*intermediateGenParticles)[0].daughters.size()<<(*intermediateGenParticles)[1].daughters.size()<<endl;
                //cout<<(*intermediateGenParticles)[0].daughters[0].pdgId<<(*intermediateGenParticles)[0].daughters[1].pdgId<<endl;
                //cout<<(*intermediateGenParticles)[1].daughters[0].pdgId<<(*intermediateGenParticles)[1].daughters[1].pdgId<<endl;
                //GG case -> weight=1. 0.25to0.25
                bool oneOfThem=false;
                if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2) {
                        if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 22) {
                                //cout<<"buh"<<endl;
                                if (fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 22) {
                                        //cout<<"sdfdfsd"<<endl;
                                        useTChiNG=true;
                                        TChiNG_weight=1.;
                                        oneOfThem=oneOfThem||true;

                                }
                        }
                }
                //ZZ case ->weight=4. 1/16 to 0.25
                if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2) {
                        if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 23) {
                                if (fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 23) {
                                        useTChiNG=true;
                                        TChiNG_weight=4.;
                                        oneOfThem=oneOfThem||true;

                                }
                        }
                }
                //cout<<"test"<<endl;
                //ZG case -> weight = 2 0.25to0.5
                if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2) {
                        //cout<<"test 2"<<endl;
                        if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 23) {
                                //cout<<"test3"<<endl;
                                if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 22) {
                                        useTChiNG=true;
                                        TChiNG_weight=2.;
                                        oneOfThem=oneOfThem||true;
                                }
                        }
                }
                if(oneOfThem==false) {
                        useTChiNG=false;
                        TChiNG_weight=1.;
                }



        }

        bool usePoint=(isTChiNG&&useTChiNG) || (!isTChiNG);
        if(usePoint) {
                if(isTChiNG) {
                        tempWeight=tempWeight*TChiNG_weight;
                }



                //m.at(WEIGHT_EWKINOPAIRPT).Fill(ewkWeights[(int)changeewk]);
                //m.at(WEIGHT_NISR).Fill(isrWeights[(int)changeisr]);
                //m.at(WEIGHT_TOPPT).Fill(*topPt_weight);
                //m.at(WEIGHT_PDF).Fill(weight_PDF);

                if(changeMET==normal) {
                        m.at(ETMISS).Fill(*ETmiss, tempWeight);
                }
                else{
                        if(changeMET==JESUP) {
                                m.at(ETMISS).Fill(met_JESu->p.Pt(),tempWeight);
                        }else{
                                if(changeMET==JESDOWN) {
                                        m.at(ETMISS).Fill(met_JESd->p.Pt(),tempWeight);
                                }else{
                                        if(changeMET==JERUP) {
                                                m.at(ETMISS).Fill(met_JERu->p.Pt(),tempWeight);
                                        }else{
                                                if(changeMET==JERDOWN) {
                                                        m.at(ETMISS).Fill(met_JERd->p.Pt(),tempWeight);
                                                }
                                        }
                                }
                        }
                }

                m.at(PT1).Fill(le1.p.Pt(), tempWeight);
                m.at(PT2).Fill(le2.p.Pt(), tempWeight);
                m.at(PT3).Fill(le3.p.Pt(), tempWeight);
                m.at(MLL).Fill((le1.vec+le2.vec).M(),tempWeight);
                m.at(NPHOTONS).Fill(*selPhotonSize,tempWeight);
                m.at(HT).Fill(*calcHt,tempWeight);
                m.at(NJETS).Fill(*selJetSize,tempWeight);
                m.at(ETA1).Fill(fabs(le1.p.Eta()),tempWeight);
                m.at(ETA2).Fill(fabs(le2.p.Eta()),tempWeight);
                m.at(ETA3).Fill(fabs(le3.p.Eta()),tempWeight);
                m.at(PHI1).Fill(fabs(le1.p.Phi()),tempWeight);
                m.at(PHI2).Fill(fabs(le2.p.Phi()),tempWeight);
                m.at(PHI3).Fill(fabs(le3.p.Phi()),tempWeight);

                m.at(MT2).Fill(*Mt2,tempWeight);


                if(!slimmed) {
                        if(*selJetSize>0) {
                                m.at(JetPt1).Fill(selJets->at(0).p.Pt(),tempWeight);
                                m.at(JetPhi1).Fill(selJets->at(0).p.Phi(),tempWeight);
                                m.at(JetEta1).Fill(selJets->at(0).p.Eta(),tempWeight);
                                if(*selJetSize>1) {
                                        m.at(JetPt2).Fill(selJets->at(1).p.Pt(),tempWeight);
                                        m.at(JetPhi2).Fill(selJets->at(1).p.Phi(),tempWeight);
                                        m.at(JetEta2).Fill(selJets->at(1).p.Eta(),tempWeight);
                                        if(*selJetSize>2) {
                                                m.at(JetPt3).Fill(selJets->at(2).p.Pt(),tempWeight);
                                                m.at(JetPhi3).Fill(selJets->at(2).p.Phi(),tempWeight);
                                                m.at(JetEta3).Fill(selJets->at(2).p.Eta(),tempWeight);
                                                if(*selJetSize>3) {
                                                        m.at(JetPt4).Fill(selJets->at(3).p.Pt(),tempWeight);
                                                        m.at(JetPhi4).Fill(selJets->at(3).p.Phi(),tempWeight);
                                                        m.at(JetEta4).Fill(selJets->at(3).p.Eta(),tempWeight);
                                                }
                                        }
                                }
                        }
                        m.at(NVTX).Fill(*nGoodVertices,tempWeight);
                        m.at(GENHT).Fill(*genHt,tempWeight);
                        m.at(ZPT).Fill((le1.vec+le2.vec).Pt(),tempWeight);
                        m.at(ST).Fill(le1.vec.Pt()+le2.vec.Pt(),tempWeight);

                        m.at(DeltaEtaLL).Fill(fabs(le1.p.Eta()-le2.p.Eta()),tempWeight);
                        m.at(DeltaPhiLL).Fill(fabs(le1.p.Phi()-le2.p.Phi()),tempWeight);
                        m.at(DeltaRLL).Fill(fabs(le1.vec.DeltaR(le2.vec)),tempWeight);

                        m.at(DeltaEtaLL_neg).Fill(le1.p.Eta()-le2.p.Eta(),tempWeight);
                        m.at(DeltaPhiLL_neg).Fill(le1.p.Phi()-le2.p.Phi(),tempWeight);
                        m.at(DeltaRLL_neg).Fill(le1.vec.DeltaR(le2.vec),tempWeight);
                        m.at(MTL3MET).Fill((le1.vec+ETmiss_vec->vec).Mt(),tempWeight);

                        m.at(NElectrons).Fill(*selElectronSize,tempWeight);
                        m.at(NMuons).Fill(*selMuonSize,tempWeight);
                }



                if(withPhoton) {
                        m.at(PTG1).Fill(selPhotons->at(0).p.Pt(),tempWeight);
                        m.at(PHIG1).Fill(selPhotons->at(0).p.Phi(),tempWeight);
                        m.at(ETAG1).Fill(selPhotons->at(0).p.Eta(),tempWeight);

                        if(selPhotons->at(0).matchedToPhoton) {
                                m.at(Fakes).Fill(1.,tempWeight);
                        }
                        if(selPhotons->at(0).matchedToElectron) {
                                m.at(Fakes).Fill(2.,tempWeight);
                        }
                        if(selPhotons->at(0).matchedToJet) {
                                m.at(Fakes).Fill(3.,tempWeight);
                        }
                        if(!((selPhotons->at(0).matchedToPhoton)||(selPhotons->at(0).matchedToElectron)||(selPhotons->at(0).matchedToJet))) {
                                m.at(Fakes).Fill(4.,tempWeight);
                        }
                        m.at(Fakes).Fill(0.,tempWeight);


                        //m.at(DELTARGL1).Fill(selPhotons->at(0).vec.DeltaR(l1->vec),tempWeight);
                        //m.at(DELTARGL2).Fill(selPhotons->at(0).vec.DeltaR(l2->vec),tempWeight);
                        m.at(DeltaPhiL1G).Fill(selPhotons->at(0).vec.DeltaPhi(l1->vec),tempWeight);
                        m.at(DeltaPhiL2G).Fill(selPhotons->at(0).vec.DeltaPhi(l2->vec),tempWeight);
                        m.at(DeltaPhiLLG).Fill(fabs((l1->vec+l2->vec).DeltaPhi(selPhotons->at(0).vec)),tempWeight);


                        if(!slimmed) {
                                //m.at(SIGMAIETAIETAG1).Fill(selPhotons->at(0).sigmaIetaIeta,tempWeight);
                                //m.at(SIGMAIPHIIPHIG1).Fill(selPhotons->at(0).sigmaIphiIphi,tempWeight);
                                //m.at(R9).Fill(selPhotons->at(0).r9,tempWeight);
                                //m.at(HOVERE).Fill(selPhotons->at(0).hOverE,tempWeight);
                                m.at(DELTARGL1).Fill(selPhotons->at(0).deltaR1,tempWeight);
                                m.at(DELTARGL2).Fill(selPhotons->at(0).deltaR2,tempWeight);
                                m.at(DeltaEtaLLG).Fill(fabs((le1.vec+le2.vec).Eta()-selPhotons->at(0).vec.Eta()), tempWeight);
                                //m.at(DeltaPhiLLG).Fill(fabs((le1.vec+le2.vec).Phi()-selPhotons->at(0).vec.Phi()),tempWeight);
                                m.at(DeltaRLLG).Fill(fabs((le1.vec+le2.vec).DeltaR(selPhotons->at(0).vec)),tempWeight);
                                m.at(MTLLG).Fill((le1.vec+le2.vec+selPhotons->at(0).vec).Mt(),tempWeight);
                                m.at(MTL1MET).Fill((le1.vec+ETmiss_vec->vec).Mt(),tempWeight);
                                m.at(MTL2MET).Fill((le2.vec+ETmiss_vec->vec).Mt(),tempWeight);
                                m.at(MTGMET).Fill((selPhotons->at(0).vec+ETmiss_vec->vec).Mt(),tempWeight);
                                m.at(MTLLMET).Fill((le1.vec+le2.vec+ETmiss_vec->vec).Mt(),tempWeight);
                                m.at(MTLLGMET).Fill((le1.vec+le2.vec+selPhotons->at(0).vec+ETmiss_vec->vec).Mt(),tempWeight);
                                m.at(STG).Fill(le1.vec.Pt()+le2.vec.Pt()+selPhotons->at(0).p.Pt(),tempWeight);
                                m.at(STMET).Fill(le1.vec.Pt()+le2.vec.Pt()+selPhotons->at(0).p.Pt()+*ETmiss,tempWeight);
                                m.at(MLLG).Fill((le1.vec+le2.vec+selPhotons->at(0).vec).M(),tempWeight);

                        }

                        //if(FindGenPhotonMatch(ev.selPhotons->at(0))) m.at(genPhotonPT).Fill(GetGenPhotonMatch(ev.selPhotons->at(0)).p.Pt(),tempWeight);
                        //
                        //if(ev.evtHasGenPhotonVeto) if(FindGenPhotonMatch(ev.selPhotons->at(0))) m.at(genPhotonPT_Veto).Fill(GetGenPhotonMatch(ev.selPhotons->at(0)).p.Pt(),tempWeight);
                        //if(ev.evtHasGenPhotonVeto) m.at(PTG1_Veto).Fill(ev.selPhotons->at(0).p.Pt(),tempWeight);
                        //if(!ev.evtHasGenPhotonVeto) if(FindGenPhotonMatch(ev.selPhotons->at(0))) m.at(genPhotonPT_NoVeto).Fill(GetGenPhotonMatch(ev.selPhotons->at(0)).p.Pt(),tempWeight);
                        //if(!ev.evtHasGenPhotonVeto) m.at(PTG1_NoVeto).Fill(ev.selPhotons->at(0).p.Pt(),tempWeight);

                }
        }
}
void myAnalyzer::Filler2D(map<Histograms2D,TH2F>& m,bool withPhoton,bool slimmed,int changePDF,changemet changeMET,changepu changePU,changeLEPSF changeLepSF,changePHOTONSF changePhotonSF,changeISR changeisr,changeEWK changeewk){


        float tempWeight=*mc_weight;
        float weight_PDF=1.;
        if(changePDF>150) {
                weight_PDF=1.;
        }else{
                weight_PDF=pdf_weights->at(changePDF);
        }
        //tempWeight = tempWeight*puWeights[(int)changePU]*lepSfWeights[(int)changeLepSF] *photonSfWeights[(int)changePhotonSF] **topPt_weight*isrWeights[(int)changeisr]*ewkWeights[(int)changeewk]*weight_PDF;
        tempWeight = tempWeight*puWeights[(int)changePU]*lepSfWeights[(int)changeLepSF] *photonSfWeights[(int)changePhotonSF] **topPt_weight*isrWeights[(int)changeisr]*ewkWeights[(int)changeewk]*weight_PDF*zzKFactor;
        if(!isData) {
                if(!isSignal) {
                        if(*isDiElectron) tempWeight=tempWeight*config_eeWeight;
                        if(*isMuonElectron||*isElectronMuon) tempWeight=tempWeight*config_emWeight;
                        if(*isDiMuon) tempWeight=tempWeight*config_mmWeight;
                }else{
                        if(*isDiElectron) tempWeight=tempWeight*config_eeEff;
                        if(*isMuonElectron||*isElectronMuon) tempWeight=tempWeight*config_emEff;
                        if(*isDiMuon) tempWeight=tempWeight*config_mmEff;
                }
        }

        //if(changeMET==normal){
        //m.at(ETMISS).Fill(*ETmiss, tempWeight);
        //}

        //m.at(PT1).Fill(*pt1, tempWeight);


        if(withPhoton) {
                //m.at(PTG1).Fill(ev.selPhotons->at(0).p.Pt(),ev.totalWeight);
                //m.at(ISRVFSR).Fill(*mll,(l1->vec+l2->vec+selPhotons->at(0).vec).M(),totalWeight);
                //m.at(PTGvsMLLG).Fill((l1->vec+l2->vec+selPhotons->at(0).vec).M(),selPhotons->at(0).p.Pt(),totalWeight);
                //m.at(MetZpt).Fill(*ETmiss,(l1->vec+l2->vec).Pt(),tempWeight);
                //m.at(MetDeltaPhiLL).Fill(*ETmiss,fabs(l1->vec.DeltaPhi(l2->vec)),tempWeight);
                //m.at(MetMllg).Fill(*ETmiss,(l1->vec+l2->vec+selPhotons->at(0).vec).M(),tempWeight);
                m.at(MetMt2).Fill(*ETmiss,*Mt2,tempWeight);
        }
}

void myAnalyzer::FillerTrigger(map<Histograms1D,TEfficiency>& m,bool withPhoton,bool TriggerBool){
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
        m.at(ETMISS).Fill(TriggerBool,*ETmiss);
        m.at(PT1).Fill(TriggerBool,*pt1);
        m.at(PT2).Fill(TriggerBool,*pt2);
        m.at(MLL).Fill(TriggerBool,*mll);
        //m.at(NPHOTONS).Fill(TriggerBool,*selPhotons.size());
        m.at(NVTX).Fill(TriggerBool,*nGoodVertices);
        //m.at(HT).Fill(TriggerBool,*calcHt);
        m.at(HT).Fill(TriggerBool,*ht);
        //m.at(GENHT).Fill(TriggerBool,*genHt);
        //m.at(NJETS).Fill(TriggerBool,ev.selJets.size());
        //m.at(ETA1).Fill(TriggerBool,fabs(ev.eta1));
        //m.at(ETA2).Fill(TriggerBool,fabs(ev.eta2));
        //m.at(PHI1).Fill(TriggerBool,fabs(ev.phi2));
        //m.at(PHI2).Fill(TriggerBool,fabs(ev.phi2));
}

//void myAnalyzer::FillerSignal(selEvent& ev, map<Histograms1D,TH1F>& m, float divideFactor=1.){
void myAnalyzer::FillerSignal(map<Histograms1D,TH1F>& m, float divideFactor,int changePDF,changemet changeMET,changepu changePU,changeLEPSF changeLepSF,changePHOTONSF changePhotonSF,changeISR changeisr,changeEWK changeewk){
        //float tempWeight=*mc_weight * *pu_weight;
        float tempWeight=*mc_weight;


        float weight_PDF=1.;
        if(changePDF>150) {
                weight_PDF=1.;
        }else{
                weight_PDF=pdf_weights->at(changePDF);
        }
        //tempWeight = tempWeight*puWeights[(int)changePU]*lepSfWeights[(int)changeLepSF] *photonSfWeights[(int)changePhotonSF] **topPt_weight*isrWeights[(int)changeisr]*ewkWeights[(int)changeewk]*weight_PDF;
        tempWeight = tempWeight*puWeights[(int)changePU]*lepSfWeights[(int)changeLepSF] *photonSfWeights[(int)changePhotonSF] **topPt_weight*isrWeights[(int)changeisr]*ewkWeights[(int)changeewk]*weight_PDF*zzKFactor;
        if(!isData) {
                if(!isSignal) {
                        if(*isDiElectron) tempWeight=tempWeight*config_eeWeight;
                        if(*isMuonElectron||*isElectronMuon) tempWeight=tempWeight*config_emWeight;
                        if(*isDiMuon) tempWeight=tempWeight*config_mmWeight;
                }else{
                        if(*isDiElectron) tempWeight=tempWeight*config_eeEff;
                        if(*isMuonElectron||*isElectronMuon) tempWeight=tempWeight*config_emEff;
                        if(*isDiMuon) tempWeight=tempWeight*config_mmEff;
                }
        }
        //m.at(WEIGHT_EWKINOPAIRPT).Fill(ewkWeights[(int)changeewk]);
        //m.at(WEIGHT_NISR).Fill(isrWeights[(int)changeisr]);
        //m.at(WEIGHT_TOPPT).Fill(*topPt_weight);
        //m.at(WEIGHT_PDF).Fill(weight_PDF);

        //m.at(MT2).Fill(*Mt2, tempWeight*1./(nGen/divideFactor));


        if(changeMET==normal) {
                m.at(ETMISS).Fill(*ETmiss, tempWeight*1./(nGen/divideFactor));
                m.at(MT2).Fill(*Mt2, tempWeight*1./(nGen/divideFactor));
                //m.at(ETMISSRAW).Fill(*ETmiss, *mc_weight);
                if(*mc_weight>0) {
                        m.at(ETMISSRAW).Fill(*ETmiss, *mc_weight);
                }
        }
        else{
                if(changeMET==JESUP) {
                        m.at(ETMISS).Fill(met_JESu->p.Pt(), tempWeight*1./(nGen/divideFactor));
                }else{
                        if(changeMET==JESDOWN) {
                                m.at(ETMISS).Fill(met_JESd->p.Pt(), tempWeight*1./(nGen/divideFactor));
                        }else{
                                if(changeMET==JERUP) {
                                        m.at(ETMISS).Fill(met_JERu->p.Pt(), tempWeight*1./(nGen/divideFactor));
                                }else{
                                        if(changeMET==JERDOWN) {
                                                m.at(ETMISS).Fill(met_JERd->p.Pt(), tempWeight*1./(nGen/divideFactor));
                                        }else{
                                                if(changeMET==METGEN) {
                                                        m.at(ETMISS).Fill(met_gen->p.Pt(), tempWeight*1./(nGen/divideFactor));
                                                }
                                        }
                                }
                        }
                }
        }
}
void myAnalyzer::FillerSignal(map<Histograms2D,TH2F>& m, float divideFactor,int changePDF,changemet changeMET,changepu changePU,changeLEPSF changeLepSF,changePHOTONSF changePhotonSF,changeISR changeisr,changeEWK changeewk){
        //float tempWeight=*mc_weight * *pu_weight;
        float tempWeight=*mc_weight;


        float weight_PDF=1.;
        if(changePDF>150) {
                weight_PDF=1.;
        }else{
                weight_PDF=pdf_weights->at(changePDF);
        }
        //tempWeight = tempWeight*puWeights[(int)changePU]*lepSfWeights[(int)changeLepSF] *photonSfWeights[(int)changePhotonSF] **topPt_weight*isrWeights[(int)changeisr]*ewkWeights[(int)changeewk]*weight_PDF;
        tempWeight = tempWeight*puWeights[(int)changePU]*lepSfWeights[(int)changeLepSF] *photonSfWeights[(int)changePhotonSF] **topPt_weight*isrWeights[(int)changeisr]*ewkWeights[(int)changeewk]*weight_PDF*zzKFactor;
        if(!isData) {
                if(!isSignal) {
                        if(*isDiElectron) tempWeight=tempWeight*config_eeWeight;
                        if(*isMuonElectron||*isElectronMuon) tempWeight=tempWeight*config_emWeight;
                        if(*isDiMuon) tempWeight=tempWeight*config_mmWeight;
                }else{
                        if(*isDiElectron) tempWeight=tempWeight*config_eeEff;
                        if(*isMuonElectron||*isElectronMuon) tempWeight=tempWeight*config_emEff;
                        if(*isDiMuon) tempWeight=tempWeight*config_mmEff;
                }
        }
        //m.at(WEIGHT_EWKINOPAIRPT).Fill(ewkWeights[(int)changeewk]);
        //m.at(WEIGHT_NISR).Fill(isrWeights[(int)changeisr]);
        //m.at(WEIGHT_TOPPT).Fill(*topPt_weight);
        //m.at(WEIGHT_PDF).Fill(weight_PDF);

        if(changeMET==normal) {
                //m.at(ETMISS).Fill(*ETmiss, tempWeight*1./(nGen/divideFactor));
                //m.at(MetZpt).Fill(*ETmiss,(l1->vec+l2->vec).Pt(),tempWeight);
                m.at(MetDeltaPhiLL).Fill(*ETmiss,fabs(l1->vec.DeltaPhi(l2->vec)),tempWeight);
                //m.at(MetMllg).Fill(*ETmiss,(l1->vec+l2->vec+selPhotons->at(0).vec).M(),tempWeight);
                m.at(MetMt2).Fill(*ETmiss,*Mt2,tempWeight);
        }
        else{
                if(changeMET==JESUP) {
                        //m.at(ETMISS).Fill(met_JESu->p.Pt(), tempWeight*1./(nGen/divideFactor));
                        //m.at(MetZpt).Fill(met_JESu->p.Pt(),(l1->vec+l2->vec).Pt(),tempWeight);
                        m.at(MetDeltaPhiLL).Fill(met_JESu->p.Pt(),fabs(l1->vec.DeltaPhi(l2->vec)),tempWeight);
                        //m.at(MetMllg).Fill(met_JESu->p.Pt(),(l1->vec+l2->vec+selPhotons->at(0).vec).M(),tempWeight);
                        m.at(MetMt2).Fill(met_JESu->p.Pt(),*Mt2,tempWeight);
                }else{
                        if(changeMET==JESDOWN) {
                                //m.at(ETMISS).Fill(met_JESd->p.Pt(), tempWeight*1./(nGen/divideFactor));
                                //m.at(MetZpt).Fill(met_JESd->p.Pt(),(l1->vec+l2->vec).Pt(),tempWeight);
                                m.at(MetDeltaPhiLL).Fill(met_JESd->p.Pt(),fabs(l1->vec.DeltaPhi(l2->vec)),tempWeight);
                                //m.at(MetMllg).Fill(met_JESd->p.Pt(),(l1->vec+l2->vec+selPhotons->at(0).vec).M(),tempWeight);
                                m.at(MetMt2).Fill(met_JESd->p.Pt(),*Mt2,tempWeight);
                        }else{
                                if(changeMET==JERUP) {
                                        //m.at(ETMISS).Fill(met_JERu->p.Pt(), tempWeight*1./(nGen/divideFactor));
                                        //m.at(MetZpt).Fill(met_JERu->p.Pt(),(l1->vec+l2->vec).Pt(),tempWeight);
                                        m.at(MetDeltaPhiLL).Fill(met_JERu->p.Pt(),fabs(l1->vec.DeltaPhi(l2->vec)),tempWeight);
                                        //m.at(MetMllg).Fill(met_JERu->p.Pt(),(l1->vec+l2->vec+selPhotons->at(0).vec).M(),tempWeight);
                                        m.at(MetMt2).Fill(met_JERu->p.Pt(),*Mt2,tempWeight);
                                }else{
                                        if(changeMET==JERDOWN) {
                                                //m.at(ETMISS).Fill(met_JERd->p.Pt(), tempWeight*1./(nGen/divideFactor));
                                                //m.at(MetZpt).Fill(met_JERd->p.Pt(),(l1->vec+l2->vec).Pt(),tempWeight);
                                                m.at(MetDeltaPhiLL).Fill(met_JERd->p.Pt(),fabs(l1->vec.DeltaPhi(l2->vec)),tempWeight);
                                                //m.at(MetMllg).Fill(met_JERd->p.Pt(),(l1->vec+l2->vec+selPhotons->at(0).vec).M(),tempWeight);
                                                m.at(MetMt2).Fill(met_JERd->p.Pt(),*Mt2,tempWeight);
                                        }
                                        else{
                                                if(changeMET==JERDOWN) {
                                                        //m.at(ETMISS).Fill(met_JERd->p.Pt(), tempWeight*1./(nGen/divideFactor));
                                                        //m.at(MetZpt).Fill(met_JERd->p.Pt(),(l1->vec+l2->vec).Pt(),tempWeight);
                                                        m.at(MetDeltaPhiLL).Fill(met_JERd->p.Pt(),fabs(l1->vec.DeltaPhi(l2->vec)),tempWeight);
                                                        //m.at(MetMllg).Fill(met_JERd->p.Pt(),(l1->vec+l2->vec+selPhotons->at(0).vec).M(),tempWeight);
                                                        m.at(MetMt2).Fill(met_JERd->p.Pt(),*Mt2,tempWeight);
                                                }else{
                                                        if(changeMET==METGEN) {
                                                                //m.at(MetZpt).Fill(met_gen->p.Pt(),(l1->vec+l2->vec).Pt(),tempWeight);
                                                                m.at(MetDeltaPhiLL).Fill(met_gen->p.Pt(),fabs(l1->vec.DeltaPhi(l2->vec)),tempWeight);
                                                                //m.at(MetMllg).Fill(met_gen->p.Pt(),(l1->vec+l2->vec+selPhotons->at(0).vec).M(),tempWeight);
                                                                m.at(MetMt2).Fill(met_gen->p.Pt(),*Mt2,tempWeight);
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }
}




void myAnalyzer::InitAllHistos(){
        if (config_selectionsToProcessMap[DILEP]==true) {
                h1Maps[dilep][EE]=InitHistograms(DILEP);
                h1Maps[dilep][MM]=InitHistograms(DILEP);
                h1Maps[dilep][EM]=InitHistograms(DILEP);
                h1Maps[dilep][LL]=InitHistograms(DILEP);

        }
        if (config_selectionsToProcessMap[SEL]==true) {
                h1Maps[sel][EE]=InitHistograms(SEL);
                h1Maps[sel][MM]=InitHistograms(SEL);
                h1Maps[sel][EM]=InitHistograms(SEL);
                h1Maps[sel][LL]=InitHistograms(SEL);
        }
        if (config_selectionsToProcessMap[ONZ]==true) {
                h1Maps[onz][EE]=InitHistograms(ONZ);
                h1Maps[onz][MM]=InitHistograms(ONZ);
                h1Maps[onz][EM]=InitHistograms(ONZ);
                h1Maps[onz][LL]=InitHistograms(ONZ);
                h2Maps[onz][LL]=Init2DHistograms(ONZ);
        }
        if (config_selectionsToProcessMap[ONZMET]==true) {
                h1Maps[onzmet0100][EE]=InitHistograms(ONZ); //<100
                h1Maps[onzmet0100][MM]=InitHistograms(ONZ);
                h1Maps[onzmet0100][EM]=InitHistograms(ONZ);
                h1Maps[onzmet0100][LL]=InitHistograms(ONZ);

                h1Maps[onzmet150][EE]=InitHistograms(ONZ); //>150
                h1Maps[onzmet150][MM]=InitHistograms(ONZ);
                h1Maps[onzmet150][EM]=InitHistograms(ONZ);
                h1Maps[onzmet150][LL]=InitHistograms(ONZ);
                //h2Maps[onzmet150][EM]=Init2DHistograms(ONZ);
                //h2Maps[onzmet150][EE]=Init2DHistograms(ONZ);
                //h2Maps[onzmet150][MM]=Init2DHistograms(ONZ);
                //h2Maps[onzmet150][LL]=Init2DHistograms(ONZ);

                h1Maps[onzmet100150][EE]=InitHistograms(ONZ); //>100<150
                h1Maps[onzmet100150][MM]=InitHistograms(ONZ);
                h1Maps[onzmet100150][EM]=InitHistograms(ONZ);
                h1Maps[onzmet100150][LL]=InitHistograms(ONZ);
        }
        if (config_selectionsToProcessMap[ControlRegionDY]==true) {
                InitWeightHistos(cr1Maps,ONZ,controlregionDY);
        }

        if (config_selectionsToProcessMap[ControlRegionTT]==true) {
                InitWeightHistos(cr1Maps,ONZ,controlregionTT);
                InitWeightHistos(cr1Maps,ONZ,controlregionTT080);
                InitWeightHistos(cr1Maps,ONZ,controlregionTT80);
        }
        if (config_selectionsToProcessMap[ControlRegionZZ]==true) {
                InitWeightHistos(cr1Maps,ControlRegionZZ,controlregionZZ);
        }
        if (config_selectionsToProcessMap[ControlRegionWZ]==true) {
                InitWeightHistos(cr1Maps,ControlRegionZZ,controlregionWZ);
        }
        if (config_selectionsToProcessMap[ValidationRegion]==true) {
                h1Maps[validationregion][EE]=InitHistograms(ONZ);
                h1Maps[validationregion][MM]=InitHistograms(ONZ);
                h1Maps[validationregion][EM]=InitHistograms(ONZ);
                h1Maps[validationregion][LL]=InitHistograms(ONZ);
                h1Maps[validationregion080][EE]=InitHistograms(ONZ);
                h1Maps[validationregion080][MM]=InitHistograms(ONZ);
                h1Maps[validationregion080][EM]=InitHistograms(ONZ);
                h1Maps[validationregion080][LL]=InitHistograms(ONZ);
                h1Maps[validationregion80][EE]=InitHistograms(ONZ);
                h1Maps[validationregion80][MM]=InitHistograms(ONZ);
                h1Maps[validationregion80][EM]=InitHistograms(ONZ);
                h1Maps[validationregion80][LL]=InitHistograms(ONZ);
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
}

map<Histograms1D,TH1F> myAnalyzer::InitCutFlowHistograms(const selectionType selection){
        map<Histograms1D,TH1F> cMap;
        cMap[CUTFLOW] = TH1F("","",5,0,5);
        cMap[CUTFLOW].Fill("triggered", 0);
        cMap[CUTFLOW].Fill("2leptons", 0);
        cMap[CUTFLOW].Fill("m50", 0);
        cMap[CUTFLOW].Fill("1photon", 0);
        cMap[CUTFLOW].Fill("Z", 0);
        return cMap;
}


map<Histograms1D,TH1F> myAnalyzer::InitCutFlowHistograms_Fine(const selectionType selection){
        map<Histograms1D,TH1F> cMap;
        //cMap[CUTFLOW_fine] = TH1F("","",18,0,18);
        cMap[CUTFLOW_fine] = TH1F("","",21,0,21);
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
        cMap[CUTFLOW_fine].Fill("ex2Lep", 0);
        cMap[CUTFLOW_fine].Fill("MT2", 0);
        cMap[CUTFLOW_fine].Fill("Met", 0);
        return cMap;
}


void myAnalyzer::InitCutFlowHistos(){
        //c1Maps["cutFlow_onZEE"]=InitCutFlowHistograms(ONZ);
        //c1Maps["cutFlow_onZMM"]=InitCutFlowHistograms(ONZ);
        //c1Maps["cutFlow_onZEM"]=InitCutFlowHistograms(ONZ);
        c1Maps[cutFlow_onZEE]=InitCutFlowHistograms(ONZ);
        c1Maps[cutFlow_onZMM]=InitCutFlowHistograms(ONZ);
        c1Maps[cutFlow_onZEM]=InitCutFlowHistograms(ONZ);
}


void myAnalyzer::InitCutFlowHistos_Fine(){
        //c1Maps["cutFlow_Fine_onZEE"]=InitCutFlowHistograms_Fine(ONZ);
        //c1Maps["cutFlow_Fine_onZMM"]=InitCutFlowHistograms_Fine(ONZ);
        //c1Maps["cutFlow_Fine_onZEM"]=InitCutFlowHistograms_Fine(ONZ);
        c1Maps[cutFlow_Fine_onZEE]=InitCutFlowHistograms_Fine(ONZ);
        c1Maps[cutFlow_Fine_onZMM]=InitCutFlowHistograms_Fine(ONZ);
        c1Maps[cutFlow_Fine_onZEM]=InitCutFlowHistograms_Fine(ONZ);
}


map<Histograms1D,TH1F> myAnalyzer::InitHistograms(const selectionType selection_){
        map<Histograms1D,TH1F> hMap;

        hMap[ETMISS] = TH1F(emptyLabelPtr, ";#p_{T}^{miss} (GeV)", 5000, 0, 5000);
        hMap[PT1] = TH1F(emptyLabelPtr, ";#p_{T}^{leading} (GeV)", 5000, 0, 5000);
        hMap[PT2] = TH1F(emptyLabelPtr, ";#p_{T}^{trailing} (GeV)", 5000, 0, 5000);
        hMap[MLL] = TH1F(emptyLabelPtr, ";#m_{ll} (GeV)", 5000, 0, 5000);
        hMap[NPHOTONS] = TH1F(emptyLabelPtr,"n_#gamma",10,0,10);
        hMap[NVTX] = TH1F(emptyLabelPtr,"n_{Vtx}",60,0,60);
        hMap[HT] = TH1F(emptyLabelPtr,"#H_{T} (GeV)",5000,0,5000);
        hMap[GENHT] = TH1F(emptyLabelPtr,"#H_{T}^{gen} (GeV)",5000,0,5000);
        hMap[NJETS] = TH1F(emptyLabelPtr,"n_{Jets}",20,0,20);
        hMap[NBJETS] = TH1F(emptyLabelPtr,"n_{BJets}",20,0,20);
        hMap[ETA1] = TH1F(emptyLabelPtr, ";|#eta_{trailing}|", 260, 0, 2.6);
        hMap[ETA2] = TH1F(emptyLabelPtr, ";|#eta_{leading}|", 260, 0, 2.6);
        hMap[PHI1] = TH1F(emptyLabelPtr, ";|#phi_{trailing}|", 350, 0, 3.5);
        hMap[PHI2] = TH1F(emptyLabelPtr, ";|#phi_{leading}|", 350, 0, 3.5);
        hMap[ZPT] = TH1F(emptyLabelPtr, ";Z_{p_T}", 5000, 0, 5000);
        hMap[ST] = TH1F(emptyLabelPtr, ";S_T", 5000, 0, 5000.);
        hMap[MT2] = TH1F(emptyLabelPtr, ";M_{T2}", 50000, 0, 5000.);
        hMap[VetoCompare] = TH1F(emptyLabelPtr, emptyLabelPtr, 2, 0, 2);
        hMap[NElectrons] = TH1F(emptyLabelPtr,"nElectrons",20,0,20);
        hMap[NMuons] = TH1F(emptyLabelPtr,"nMuons",20,0,20);

        hMap[JetPt1] = TH1F(emptyLabelPtr, ";#it{p}_{T}^{Jet 1} (GeV)", 5000, 0, 5000);
        hMap[JetPt2] = TH1F(emptyLabelPtr, ";#it{p}_{T}^{Jet 2} (GeV)", 5000, 0, 5000);
        hMap[JetPt3] = TH1F(emptyLabelPtr, ";#it{p}_{T}^{Jet 3} (GeV)", 5000, 0, 5000);
        hMap[JetPt4] = TH1F(emptyLabelPtr, ";#it{p}_{T}^{Jet 4} (GeV)", 5000, 0, 5000);
        hMap[JetPhi1] = TH1F(emptyLabelPtr, ";|#phi_{Jet 1}|", 350, 0, 3.5);
        hMap[JetPhi2] = TH1F(emptyLabelPtr, ";|#phi_{Jet 2}|", 350, 0, 3.5);
        hMap[JetPhi3] = TH1F(emptyLabelPtr, ";|#phi_{Jet 3}|", 350, 0, 3.5);
        hMap[JetPhi4] = TH1F(emptyLabelPtr, ";|#phi_{Jet 4}|", 350, 0, 3.5);
        hMap[JetEta1] = TH1F(emptyLabelPtr, ";|#eta_{Jet 1}|", 260, 0, 2.6);
        hMap[JetEta2] = TH1F(emptyLabelPtr, ";|#eta_{Jet 2}|", 260, 0, 2.6);
        hMap[JetEta3] = TH1F(emptyLabelPtr, ";|#eta_{Jet 3}|", 260, 0, 2.6);
        hMap[JetEta4] = TH1F(emptyLabelPtr, ";|#eta_{Jet 4}|", 260, 0, 2.6);

        //hMap[WEIGHT_NISR] = TH1F(emptyLabelPtr,"weight",10000,0,10);
        //hMap[WEIGHT_TOPPT] = TH1F(emptyLabelPtr,"weight",10000,0,10);
        //hMap[WEIGHT_EWKINOPAIRPT] = TH1F(emptyLabelPtr,"weight",10000,0,10);
        //////hMap[WEIGHT_LEPTONPAIRPT] = TH1F(emptyLabelPtr,"weight",10000,0,10);
        //hMap[WEIGHT_PDF]=TH1F(emptyLabelPtr,"weight",10000,0,10);


        hMap[DeltaEtaLL] = TH1F(emptyLabelPtr, ";#Delta#Eta_{ll}", 24000, -12., 12.);
        hMap[DeltaPhiLL] = TH1F(emptyLabelPtr, ";#Delta#Phi_{ll}", 24000, -12., 12.);
        hMap[DeltaRLL] = TH1F(emptyLabelPtr, ";#DeltaR_{ll}", 24000, -12, 12.);

        hMap[DeltaEtaLL_neg] = TH1F(emptyLabelPtr, ";#Delta#Eta_{ll}", 24000, -12., 12.);
        hMap[DeltaPhiLL_neg] = TH1F(emptyLabelPtr, ";#Delta#Phi_{ll}", 24000, -12., 12.);
        hMap[DeltaRLL_neg] = TH1F(emptyLabelPtr, ";#DeltaR_{ll}", 24000, -12., 12.);
        //
        if(selection_!=ControlRegionZZ) {
                hMap[DeltaPhiLLMet] = TH1F(emptyLabelPtr, "#Delta#Phi_{ll,MET}", 24000, -12., 12.);
                hMap[DeltaEtaLLMet] = TH1F(emptyLabelPtr, ";#Delta#Eta_{ll,MET}", 24000, -12., 12.);
                hMap[DeltaRLLMet] = TH1F(emptyLabelPtr, "#DeltaR_{ll,MET}", 24000, -12., 12.);
                hMap[MTLL] = TH1F(emptyLabelPtr, ";m_{T}^{ll}", 5000, 0, 5000);
        }

        if(selection_==ControlRegionZZ) {
                hMap[PT3] = TH1F(emptyLabelPtr, ";#p_{T}^{trailing} (GeV)", 5000, 0, 5000);
                hMap[PT4] = TH1F(emptyLabelPtr, ";#p_{T}^{trailing} (GeV)", 5000, 0, 5000);
                hMap[ETA3] = TH1F(emptyLabelPtr, ";|#eta_{leading}|", 260, 0, 2.6);
                hMap[ETA4] = TH1F(emptyLabelPtr, ";|#eta_{leading}|", 260, 0, 2.6);
                hMap[PHI3] = TH1F(emptyLabelPtr, ";|#phi_{leading}|", 350, 0, 3.5);
                hMap[PHI4] = TH1F(emptyLabelPtr, ";|#phi_{leading}|", 350, 0, 3.5);
                hMap[ZPT2] = TH1F(emptyLabelPtr, ";Z_{p_T}", 5000, 0, 5000);
                hMap[MLL2] = TH1F(emptyLabelPtr, ";#m_{ll} (GeV)", 5000, 0, 5000);
                hMap[MTL3MET] = TH1F(emptyLabelPtr, ";#m_{T}^{l2,Met} (GeV)", 5000, 0, 5000);
                hMap[MLLLL] = TH1F(emptyLabelPtr, ";#m_{llll} (GeV)", 5000, 0, 5000);

        }

        if ((selection_==PHOTON)||(selection_==SEL)||(selection_==ONZ)) {

                hMap[Fakes] = TH1F(emptyLabelPtr, emptyLabelPtr, 5, 0, 5);


                hMap[DELTARGL1] = TH1F(emptyLabelPtr, ";#DeltaR_{l1,#gamma}", 24000, -12., 12.);
                hMap[DELTARGL2] = TH1F(emptyLabelPtr, ";#DeltaR_{l2,#gamma}", 24000, -12., 12.);
                hMap[DeltaPhiL1G] = TH1F(emptyLabelPtr, ";#DeltaPhi_{l1,#gamma}", 24000, -12., 12.);
                hMap[DeltaPhiL2G] = TH1F(emptyLabelPtr, ";#DeltaPhi_{l2,#gamma}", 24000, -12., 12.);
                hMap[DeltaPhiLLG] = TH1F(emptyLabelPtr, ";#Delta#Phi_{ll,#gamma}", 2400, -12., 12.);

                hMap[ML1G] = TH1F(emptyLabelPtr, ";m_{l1,#gamma}", 5000, 0., 5000.);
                hMap[ML2G] = TH1F(emptyLabelPtr, ";m_{l2,#gamma}", 5000, 0., 5000.);

                hMap[PTG1] = TH1F(emptyLabelPtr, ";#p_{T}^{#gamma} (GeV)", 5000, 0, 5000);
                hMap[ETAG1] = TH1F(emptyLabelPtr, ";|#eta_{#gamma}|", 260, 0, 2.6);
                hMap[PHIG1] = TH1F(emptyLabelPtr, ";|#phi_{#gamma}|", 350, 0, 3.5);
                //hMap[SIGMAIETAIETAG1] = TH1F(emptyLabelPtr, ";#sigma_{i#etai#eta}^{#gamma}", 400, 0, 0.04);
                //hMap[SIGMAIPHIIPHIG1] = TH1F(emptyLabelPtr, ";#sigma_{i#phii#phi}^{#gamma}", 2000, 0, 0.2);
                //hMap[R9] = TH1F(emptyLabelPtr, ";r9", 1500, 0, 1.5);
                //hMap[HOVERE] = TH1F(emptyLabelPtr, ";H/E", 1000, 0, 0.1);
                //hMap[DELTARGL1] = TH1F(emptyLabelPtr, ";#DeltaR_{l1,#gamma}", 24000, -12., 12.);
                //hMap[DELTARGL2] = TH1F(emptyLabelPtr, ";#DeltaR_{l2,#gamma}", 24000, -12., 12.);
                hMap[DeltaRLLG] = TH1F(emptyLabelPtr, ";#DeltaR_{ll,#gamma}", 24000, -12., 12.);
                hMap[DeltaEtaLLG] = TH1F(emptyLabelPtr, ";#Delta#Eta_{ll,#gamma}", 2400, -12., 12.);
                //hMap[DeltaPhiLLG] = TH1F(emptyLabelPtr, ";#Delta#Phi_{ll,#gamma}", 2400, -12., 12.);
                hMap[DeltaPhiGMet] = TH1F(emptyLabelPtr, ";#Delta#Phi_{met,#gamma}", 2400, -12., 6.);
                hMap[DeltaRGMet] = TH1F(emptyLabelPtr, ";#DeltaR_{met,#gamma}", 2400, -12., 6.);
                hMap[STG] = TH1F(emptyLabelPtr, ";S_T", 5000, 0, 5000.);
                hMap[STMET] = TH1F(emptyLabelPtr, ";S_T + #it{p}_{T}^{miss} (GeV)", 5000, 0, 5000.);
                hMap[MTLLG] = TH1F(emptyLabelPtr, ";m_{T}^{ll#gamma}", 5000, 0, 5000);
                hMap[MTL1MET] = TH1F(emptyLabelPtr, ";m_{T}^{l1,met}", 5000, 0, 5000);
                hMap[MTL2MET] = TH1F(emptyLabelPtr, ";m_{T}^{l2,met}", 5000, 0, 5000);
                hMap[MTGMET] = TH1F(emptyLabelPtr, ";m_{T}^{#gamma,met}", 5000, 0, 5000);
                hMap[MTLLMET] = TH1F(emptyLabelPtr, ";m_{T}^{ll,met}", 5000, 0, 5000);
                hMap[MTLLGMET] = TH1F(emptyLabelPtr, ";m_{T}^{ll#gamma,met}", 5000, 0, 5000);
                hMap[MOTHERID] = TH1F(emptyLabelPtr, ";ID_{mother}", 200, 0, 200);
                hMap[MLLG] = TH1F(emptyLabelPtr, ";m_{ll#gamma}", 5000, 0, 5000);
                hMap[PT_llg] = TH1F(emptyLabelPtr, ";p_{T}^{ll#gamma}", 5000, 0, 5000);
                hMap[MZG_exo] = TH1F(emptyLabelPtr, ";m_{Z#gamma}", 5000, 0, 5000);
                hMap[gammaMotherID] = TH1F(emptyLabelPtr, ";motherID_{#gamma}", 5000, 0, 5000);
                hMap[genPhotonPT] = TH1F(emptyLabelPtr, ";gen p_T^{#gamma,matched}", 5000, 0, 5000);
                hMap[genPhotonPT_Veto] = TH1F(emptyLabelPtr, ";gen p_T^{#gamma,matched,veto}", 5000, 0, 5000);
                hMap[PTG1_Veto] = TH1F(emptyLabelPtr, ";gen p_T^{#gamma,veto}", 5000, 0, 5000);
                hMap[genPhotonPT_NoVeto] = TH1F(emptyLabelPtr, ";gen p_T^{#gamma,matched,Noveto}", 5000, 0, 5000);
                hMap[PTG1_NoVeto] = TH1F(emptyLabelPtr, ";gen p_T^{#gamma,Noveto}", 5000, 0, 5000);
                //hMap[VetoCompare] = TH1F(emptyLabelPtr, emptyLabelPtr, 2, 0, 2);
        }
        return hMap;
}


map<Histograms1D,TH1F> myAnalyzer::InitHistogramsSlimmed(const selectionType selection_){
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
        hMap[MT2] = TH1F("", ";M_{T2}", 50000, 0, 5000.);
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

        //hMap[WEIGHT_NISR] = TH1F("","weight",10000,0,10);
        //hMap[WEIGHT_TOPPT] = TH1F("","weight",10000,0,10);
        //hMap[WEIGHT_EWKINOPAIRPT] = TH1F("","weight",10000,0,10);
        //hMap[WEIGHT_PDF]=TH1F("","weight",10000,0,10);

        //hMap[DELTARGL1] = TH1F(emptyLabelPtr, ";#DeltaR_{l1,#gamma}", 24000, -12., 12.);
        //hMap[DELTARGL2] = TH1F(emptyLabelPtr, ";#DeltaR_{l2,#gamma}", 24000, -12., 12.);
        //hMap[DeltaPhiL1G] = TH1F(emptyLabelPtr, ";#DeltaPhi_{l1,#gamma}", 24000, -12., 12.);
        //hMap[DeltaPhiL2G] = TH1F(emptyLabelPtr, ";#DeltaPhi_{l2,#gamma}", 24000, -12., 12.);
        //hMap[DeltaPhiLLG] = TH1F(emptyLabelPtr, ";#Delta#Phi_{ll,#gamma}", 2400, -12., 12.);


        //hMap[DeltaEtaLL] = TH1F("", ";#Delta#Eta_{ll}", 24000, -12., 12.);
        hMap[DeltaPhiLL] = TH1F("", ";#Delta#Phi_{ll}", 24000, -12., 12.);
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

        if(selection_==ControlRegionZZ) {
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

        if ((selection_==PHOTON)||(selection_==SEL)||(selection_==ONZ)) {
                //hMap[PTG1] = TH1F("", ";#it{p}_{T}^{#gamma 1} (GeV)", 5000, 0, 5000);
                //hMap[ETAG1] = TH1F("", ";|#eta_{#gamma 1}|", 260, 0, 2.6);
                //hMap[PHIG1] = TH1F("", ";|#phi_{#gamma 1}|", 350, 0, 3.5);

                //hMap[Fakes] = TH1F("", "", 5, 0, 5);


                //hMap[ML1G] = TH1F(emptyLabelPtr, ";m_{l1,#gamma}", 5000, 0., 5000.);
                //hMap[ML2G] = TH1F(emptyLabelPtr, ";m_{l2,#gamma}", 5000, 0., 5000.);
                //hMap[MLLG] = TH1F(emptyLabelPtr, ";m_{ll#gamma}", 5000, 0, 5000);

                //hMap[DELTARGL1] = TH1F(emptyLabelPtr, ";#DeltaR_{l1,#gamma}", 24000, -12., 12.);
                //hMap[DELTARGL2] = TH1F(emptyLabelPtr, ";#DeltaR_{l2,#gamma}", 24000, -12., 12.);
                //hMap[DeltaPhiL1G] = TH1F(emptyLabelPtr, ";#DeltaPhi_{l1,#gamma}", 24000, -12., 12.);
                //hMap[DeltaPhiL2G] = TH1F(emptyLabelPtr, ";#DeltaPhi_{l2,#gamma}", 24000, -12., 12.);
                //hMap[DeltaPhiLLG] = TH1F(emptyLabelPtr, ";#Delta#Phi_{ll,#gamma}", 2400, -12., 12.);

                hMap[Fakes] = TH1F(emptyLabelPtr, emptyLabelPtr, 5, 0, 5);


                //hMap[DELTARGL1] = TH1F(emptyLabelPtr, ";#DeltaR_{l1,#gamma}", 24000, -12., 12.);
                //hMap[DELTARGL2] = TH1F(emptyLabelPtr, ";#DeltaR_{l2,#gamma}", 24000, -12., 12.);
                hMap[DeltaPhiL1G] = TH1F(emptyLabelPtr, ";#DeltaPhi_{l1,#gamma}", 24000, -12., 12.);
                hMap[DeltaPhiL2G] = TH1F(emptyLabelPtr, ";#DeltaPhi_{l2,#gamma}", 24000, -12., 12.);
                hMap[DeltaPhiLLG] = TH1F(emptyLabelPtr, ";#Delta#Phi_{ll,#gamma}", 2400, -12., 12.);

                hMap[ML1G] = TH1F(emptyLabelPtr, ";m_{l1,#gamma}", 5000, 0., 5000.);
                hMap[ML2G] = TH1F(emptyLabelPtr, ";m_{l2,#gamma}", 5000, 0., 5000.);

                hMap[PTG1] = TH1F(emptyLabelPtr, ";#it{p}_{T}^{#gamma 1} (GeV)", 5000, 0, 5000);
                hMap[ETAG1] = TH1F(emptyLabelPtr, ";|#eta_{#gamma 1}|", 260, 0, 2.6);
                hMap[PHIG1] = TH1F(emptyLabelPtr, ";|#phi_{#gamma 1}|", 350, 0, 3.5);
                //hMap[SIGMAIETAIETAG1] = TH1F(emptyLabelPtr, ";#sigma_{i#etai#eta}^{#gamma 1}", 400, 0, 0.04);
                //hMap[SIGMAIPHIIPHIG1] = TH1F(emptyLabelPtr, ";#sigma_{i#phii#phi}^{#gamma 1}", 2000, 0, 0.2);
                //hMap[R9] = TH1F(emptyLabelPtr, ";r9", 1500, 0, 1.5);
                //hMap[HOVERE] = TH1F(emptyLabelPtr, ";H/E", 1000, 0, 0.1);
                //hMap[DELTARGL1] = TH1F(emptyLabelPtr, ";#DeltaR_{l1,#gamma}", 24000, -12., 12.);
                //hMap[DELTARGL2] = TH1F(emptyLabelPtr, ";#DeltaR_{l2,#gamma}", 24000, -12., 12.);
                //hMap[DeltaRLLG] = TH1F(emptyLabelPtr, ";#DeltaR_{ll,#gamma}", 24000, -12., 12.);
                //hMap[DeltaEtaLLG] = TH1F(emptyLabelPtr, ";#Delta#Eta_{ll,#gamma}", 2400, -12., 12.);
                //hMap[DeltaPhiLLG] = TH1F(emptyLabelPtr, ";#Delta#Phi_{ll,#gamma}", 2400, -12., 12.);
                //hMap[DeltaPhiGMet] = TH1F(emptyLabelPtr, ";#Delta#Phi_{met,#gamma}", 2400, -12., 6.);
                //hMap[DeltaRGMet] = TH1F(emptyLabelPtr, ";#DeltaR_{met,#gamma}", 2400, -12., 6.);
                //hMap[STG] = TH1F(emptyLabelPtr, ";S_T", 5000, 0, 5000.);
                //hMap[STMET] = TH1F(emptyLabelPtr, ";S_T + #it{p}_{T}^{miss} (GeV)", 5000, 0, 5000.);
                //hMap[MTLLG] = TH1F(emptyLabelPtr, ";m_{T}^{ll#gamma}", 5000, 0, 5000);
                //hMap[MTL1MET] = TH1F(emptyLabelPtr, ";m_{T}^{l1,met}", 5000, 0, 5000);
                //hMap[MTL2MET] = TH1F(emptyLabelPtr, ";m_{T}^{l2,met}", 5000, 0, 5000);
                //hMap[MTGMET] = TH1F(emptyLabelPtr, ";m_{T}^{#gamma,met}", 5000, 0, 5000);
                //hMap[MTLLMET] = TH1F(emptyLabelPtr, ";m_{T}^{ll,met}", 5000, 0, 5000);
                //hMap[MTLLGMET] = TH1F(emptyLabelPtr, ";m_{T}^{ll#gamma,met}", 5000, 0, 5000);
                //hMap[MOTHERID] = TH1F(emptyLabelPtr, ";ID_{mother}", 200, 0, 200);
                hMap[MLLG] = TH1F(emptyLabelPtr, ";m_{ll#gamma}", 5000, 0, 5000);
                //hMap[PT_llg] = TH1F(emptyLabelPtr, ";p_{T}^{ll#gamma}", 5000, 0, 5000);
                //hMap[MZG_exo] = TH1F(emptyLabelPtr, ";m_{Z#gamma}", 5000, 0, 5000);
                //hMap[gammaMotherID] = TH1F(emptyLabelPtr, ";motherID_{#gamma}", 5000, 0, 5000);
                //hMap[genPhotonPT] = TH1F(emptyLabelPtr, ";gen p_T^{#gamma,matched}", 5000, 0, 5000);
                //hMap[genPhotonPT_Veto] = TH1F(emptyLabelPtr, ";gen p_T^{#gamma,matched,veto}", 5000, 0, 5000);
                //hMap[PTG1_Veto] = TH1F(emptyLabelPtr, ";gen p_T^{#gamma,veto}", 5000, 0, 5000);
                //hMap[genPhotonPT_NoVeto] = TH1F(emptyLabelPtr, ";gen p_T^{#gamma,matched,Noveto}", 5000, 0, 5000);
                //hMap[PTG1_NoVeto] = TH1F(emptyLabelPtr, ";gen p_T^{#gamma,Noveto}", 5000, 0, 5000);



        }
        return hMap;
}

map<Histograms2D,TH2F> myAnalyzer::Init2DHistograms(const selectionType selection_){
        map<Histograms2D,TH2F> h2Map;
        if ((selection_==PHOTON)||(selection_==SEL)||(selection_==ONZ)) {
                //h2Map[PTGvsMLLG] = TH2F("",";#it{m}_{ll#gamma};#it{p}_{T}^{gamma}",1000,0,1000,1000,0,1000);
                //h2Map[ISRVFSR] = TH2F("",";#it{m}_{ll};#it{m}{ll#gamma}",200,0,200,500,0,500);
                //h2Map[MetZpt] = TH2F(emptyLabelPtr,";met;zpt",200,0,1000,200,0,1000);
                //h2Map[MetDeltaPhiLL] = TH2F(emptyLabelPtr,";met;DeltaPhiLL",200,0,1000,400,0,4);
                //h2Map[MetMllg] = TH2F(emptyLabelPtr,";met;#it{m}{ll#gamma}",200,0,1000,200,0,1000);
                h2Map[MetMt2] = TH2F(emptyLabelPtr,";met;#it{m}{T2}",200,0,1000,400,0,2000);
        }
        return h2Map;
}





map<Histograms1D,TH1F> myAnalyzer::InitSignalScanHistograms(const selectionType selection_){
        map<Histograms1D,TH1F> sMap;

        //sMap[WEIGHT_NISR] = TH1F("","weight",10000,0,10);
        //sMap[WEIGHT_TOPPT] = TH1F("","weight",10000,0,10);
        //sMap[WEIGHT_EWKINOPAIRPT] = TH1F("","weight",10000,0,10);
        //sMap[WEIGHT_PDF] = TH1F("","weight",10000,0,10);
        //sMap[WEIGHT_NISR] = TH1F("","weight",100,0,10);
        //sMap[WEIGHT_TOPPT] = TH1F("","weight",100,0,10);
        //sMap[WEIGHT_EWKINOPAIRPT] = TH1F("","weight",10000,0,10);
        //sMap[WEIGHT_PDF] = TH1F("","weight",100,0,10);

        //sMap[WEIGHT_NISR] = TH1F(emptyLabelPtr,weightLabelPtr,1000,0,10);
        //sMap[WEIGHT_TOPPT] = TH1F(emptyLabelPtr,weightLabelPtr,1000,0,10);
        //sMap[WEIGHT_EWKINOPAIRPT] = TH1F(emptyLabelPtr,weightLabelPtr,1000,0,10);
        //sMap[WEIGHT_PDF] = TH1F(emptyLabelPtr,weightLabelPtr,1000,0,10);

        sMap[ETMISS] = TH1F(emptyLabelPtr, metLabelPtr, 500, 0, 5000);
        sMap[ETMISSRAW] = TH1F(emptyLabelPtr, metLabelPtr, 500, 0, 5000);
        sMap[MT2] = TH1F(emptyLabelPtr, metLabelPtr, 100, 0, 2000);

        return sMap;
}
map<Histograms2D,TH2F> myAnalyzer::InitSignalScanHistograms2D(const selectionType selection_){
        map<Histograms2D,TH2F> sMap;

        //sMap[WEIGHT_NISR] = TH1F("","weight",10000,0,10);
        //sMap[WEIGHT_TOPPT] = TH1F("","weight",10000,0,10);
        //sMap[WEIGHT_EWKINOPAIRPT] = TH1F("","weight",10000,0,10);
        //sMap[WEIGHT_PDF] = TH1F("","weight",10000,0,10);
        //sMap[WEIGHT_NISR] = TH1F("","weight",100,0,10);
        //sMap[WEIGHT_TOPPT] = TH1F("","weight",100,0,10);
        //sMap[WEIGHT_EWKINOPAIRPT] = TH1F("","weight",10000,0,10);
        //sMap[WEIGHT_PDF] = TH1F("","weight",100,0,10);
        //sMap[WEIGHT_NISR] = TH1F(emptyLabelPtr,weightLabelPtr,10000,0,10);
        //sMap[WEIGHT_TOPPT] = TH1F(emptyLabelPtr,weightLabelPtr,10000,0,10);
        //sMap[WEIGHT_EWKINOPAIRPT] = TH1F(emptyLabelPtr,weightLabelPtr,10000,0,10);
        //sMap[WEIGHT_PDF] = TH1F(emptyLabelPtr,weightLabelPtr,10000,0,10);

        //sMap[MetZpt] = TH2F(emptyLabelPtr,";met;zpt",200,0,1000,200,0,1000);
        //sMap[MetDeltaPhiLL] = TH2F(emptyLabelPtr,";met;DeltaPhiLL",200,0,1000,400,0,4);
        //sMap[MetMllg] = TH2F(emptyLabelPtr,";met;#it{m}{ll#gamma}",200,0,1000,200,0,1000);
        //sMap[MetMt2] = TH2F(emptyLabelPtr,";met;#it{m}{T2}",200,0,1000,400,0,2000);
        //sMap[MetZpt] = TH2F(emptyLabelPtr,";met;zpt",200,0,2000,100,0,1000);
        sMap[MetDeltaPhiLL] = TH2F(emptyLabelPtr,";met;DeltaPhiLL",200,0,2000,40,0,4);
        //sMap[MetMllg] = TH2F(emptyLabelPtr,";met;#it{m}{ll#gamma}",200,0,2000,200,0,1000);
        sMap[MetMt2] = TH2F(emptyLabelPtr,";met;#it{m}{T2}",200,0,2000,200,0,2000);

        return sMap;
}

//void myAnalyzer::InitSignalScanHistos(string masspoint){
void myAnalyzer::InitSignalScanHistos(SignalPoint masspoint){
        s1Maps[masspoint][sig][LL][nom]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sigMt2][LL][nom]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EE][nom]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][MM][nom]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EM][nom]=InitSignalScanHistograms(ONZ);

        //s2Maps[masspoint][sig][LL][nom]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sigMt2][LL][nom]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EE][nom]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][MM][nom]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EM][nom]=InitSignalScanHistograms2D(ONZ);

        s1Maps[masspoint][validationregion][LL][nom]=InitSignalScanHistograms(ONZ);
        s1Maps[masspoint][controlregionTT][EM][nom]=InitSignalScanHistograms(ONZ);
        s1Maps[masspoint][controlregionDY][LL][nom]=InitSignalScanHistograms(ONZ);
        s1Maps[masspoint][controlregionWZ][LL][nom]=InitSignalScanHistograms(ControlRegionZZ);
        s1Maps[masspoint][controlregionZZ][LL][nom]=InitSignalScanHistograms(ControlRegionZZ);
        //s1Maps[masspoint][validationregion][LL][nom]=InitSignalScanHistograms(ONZ);


        s1Maps[masspoint][sig080][LL][nom]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig080Mt2][LL][nom]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig080][EE][nom]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig080][MM][nom]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig080][EM][nom]=InitSignalScanHistograms(ONZ);
        s1Maps[masspoint][sig80][LL][nom]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig80Mt2][LL][nom]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig80][EE][nom]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig80][MM][nom]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig80][EM][nom]=InitSignalScanHistograms(ONZ);

        //s2Maps[masspoint][sig080][LL][nom]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig080Mt2][LL][nom]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig080][EE][nom]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig080][MM][nom]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig080][EM][nom]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig80][LL][nom]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig80Mt2][LL][nom]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig80][EE][nom]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig80][MM][nom]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig80][EM][nom]=InitSignalScanHistograms2D(ONZ);

        s1Maps[masspoint][sig][LL][JESu]=InitSignalScanHistograms(ONZ);
        s1Maps[masspoint][sig][LL][JESd]=InitSignalScanHistograms(ONZ);
        s1Maps[masspoint][sig][LL][JERu]=InitSignalScanHistograms(ONZ);
        s1Maps[masspoint][sig][LL][JERd]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sigMt2][LL][JESu]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sigMt2][LL][JESd]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sigMt2][LL][JERu]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sigMt2][LL][JERd]=InitSignalScanHistograms(ONZ);

        //s2Maps[masspoint][sig][LL][JESu]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][LL][JESd]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][LL][JERu]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][LL][JERd]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sigMt2][LL][JESu]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sigMt2][LL][JESd]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sigMt2][LL][JERu]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sigMt2][LL][JERd]=InitSignalScanHistograms2D(ONZ);

        //s1Maps[masspoint][sig][EE][JESu]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EE][JESd]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EE][JERu]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EE][JERd]=InitSignalScanHistograms(ONZ);

        //s2Maps[masspoint][sig][EE][JESu]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EE][JESd]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EE][JERu]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EE][JERd]=InitSignalScanHistograms2D(ONZ);

        //s1Maps[masspoint][sig][MM][JESu]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][MM][JESd]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][MM][JERu]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][MM][JERd]=InitSignalScanHistograms(ONZ);

        //s2Maps[masspoint][sig][MM][JESu]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][MM][JESd]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][MM][JERu]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][MM][JERd]=InitSignalScanHistograms2D(ONZ);

        //s1Maps[masspoint][sig][EM][JESu]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EM][JESd]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EM][JERu]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EM][JERd]=InitSignalScanHistograms(ONZ);

        //s2Maps[masspoint][sig][EM][JESu]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EM][JESd]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EM][JERu]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EM][JERd]=InitSignalScanHistograms2D(ONZ);

        //s2Maps[masspoint][sig][LL][GENMET]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EE][GENMET]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][MM][GENMET]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EM][GENMET]=InitSignalScanHistograms2D(ONZ);

        s1Maps[masspoint][sig][LL][GENMET]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sigMt2][LL][GENMET]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EE][GENMET]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][MM][GENMET]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EM][GENMET]=InitSignalScanHistograms(ONZ);

        s1Maps[masspoint][sig][LL][LEPSFUP]=InitSignalScanHistograms(ONZ);
        s1Maps[masspoint][sig][LL][LEPSFDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sigMt2][LL][LEPSFUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sigMt2][LL][LEPSFDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EE][LEPSFUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EE][LEPSFDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EM][LEPSFUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EM][LEPSFDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][MM][LEPSFUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][MM][LEPSFDOWN]=InitSignalScanHistograms(ONZ);

        //s2Maps[masspoint][sig][LL][LEPSFUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][LL][LEPSFDOWN]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EE][LEPSFUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EE][LEPSFDOWN]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EM][LEPSFUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EM][LEPSFDOWN]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][MM][LEPSFUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][MM][LEPSFDOWN]=InitSignalScanHistograms2D(ONZ);

        s1Maps[masspoint][sig][LL][PHOTONSFUP]=InitSignalScanHistograms(ONZ);
        s1Maps[masspoint][sig][LL][PHOTONSFDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sigMt2][LL][PHOTONSFUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sigMt2][LL][PHOTONSFDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EE][PHOTONSFUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EE][PHOTONSFDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][MM][PHOTONSFUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][MM][PHOTONSFDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EM][PHOTONSFUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EM][PHOTONSFDOWN]=InitSignalScanHistograms(ONZ);

        //s2Maps[masspoint][sig][LL][PHOTONSFUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][LL][PHOTONSFDOWN]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EE][PHOTONSFUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EE][PHOTONSFDOWN]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][MM][PHOTONSFUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][MM][PHOTONSFDOWN]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EM][PHOTONSFUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EM][PHOTONSFDOWN]=InitSignalScanHistograms2D(ONZ);

        s1Maps[masspoint][sig][LL][PUUP]=InitSignalScanHistograms(ONZ);
        s1Maps[masspoint][sig][LL][PUDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sigMt2][LL][PUUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sigMt2][LL][PUDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EE][PUUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EE][PUDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][MM][PUUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][MM][PUDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EM][PUUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EM][PUDOWN]=InitSignalScanHistograms(ONZ);

        //s2Maps[masspoint][sig][LL][PUUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][LL][PUDOWN]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EE][PUUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EE][PUDOWN]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][MM][PUUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][MM][PUDOWN]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EM][PUUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EM][PUDOWN]=InitSignalScanHistograms2D(ONZ);

        s1Maps[masspoint][sig][LL][NOPUUP]=InitSignalScanHistograms(ONZ);
        s1Maps[masspoint][sig][LL][NOPUDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sigMt2][LL][NOPUUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sigMt2][LL][NOPUDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EE][NOPUUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EE][NOPUDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][MM][NOPUUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][MM][NOPUDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EM][NOPUUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EM][NOPUDOWN]=InitSignalScanHistograms(ONZ);

        //s2Maps[masspoint][sig][LL][NOPUUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][LL][NOPUDOWN]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EE][NOPUUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EE][NOPUDOWN]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][MM][NOPUUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][MM][NOPUDOWN]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EM][NOPUUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EM][NOPUDOWN]=InitSignalScanHistograms2D(ONZ);

        s1Maps[masspoint][sig][LL][NOPU]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sigMt2][LL][NOPU]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EE][NOPU]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][MM][NOPU]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EM][NOPU]=InitSignalScanHistograms(ONZ);

        //s2Maps[masspoint][sig][LL][NOPU]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EE][NOPU]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][MM][NOPU]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EM][NOPU]=InitSignalScanHistograms2D(ONZ);

        s1Maps[masspoint][sig][LL][ISRUP]=InitSignalScanHistograms(ONZ);
        s1Maps[masspoint][sig][LL][ISRDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sigMt2][LL][ISRUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sigMt2][LL][ISRDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EE][ISRUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EE][ISRDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][MM][ISRUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][MM][ISRDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EM][ISRUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EM][ISRDOWN]=InitSignalScanHistograms(ONZ);

        //s2Maps[masspoint][sig][LL][ISRUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][LL][ISRDOWN]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EE][ISRUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EE][ISRDOWN]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][MM][ISRUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][MM][ISRDOWN]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EM][ISRUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EM][ISRDOWN]=InitSignalScanHistograms2D(ONZ);

        s1Maps[masspoint][sig][LL][EWKUP]=InitSignalScanHistograms(ONZ);
        s1Maps[masspoint][sig][LL][EWKDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sigMt2][LL][EWKUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sigMt2][LL][EWKDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EE][EWKUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EE][EWKDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][MM][EWKUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][MM][EWKDOWN]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EM][EWKUP]=InitSignalScanHistograms(ONZ);
        //s1Maps[masspoint][sig][EM][EWKDOWN]=InitSignalScanHistograms(ONZ);

        //s2Maps[masspoint][sig][LL][EWKUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][LL][EWKDOWN]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EE][EWKUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EE][EWKDOWN]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][MM][EWKUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][MM][EWKDOWN]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EM][EWKUP]=InitSignalScanHistograms2D(ONZ);
        //s2Maps[masspoint][sig][EM][EWKDOWN]=InitSignalScanHistograms2D(ONZ);

        //for(int i=0; i<nWeights; i++){
        for(int i=config_pdfStart; i<config_pdfStart+config_pdfStep; i++) {
                //s1Maps[masspoint][sig][EE][PDFNAMES[i]]=InitSignalScanHistograms(ONZ);
                //s1Maps[masspoint][sig][MM][PDFNAMES[i]]=InitSignalScanHistograms(ONZ);
                //s1Maps[masspoint][sig][EM][PDFNAMES[i]]=InitSignalScanHistograms(ONZ);
                s1Maps[masspoint][sig][LL][PDFNAMES[i]]=InitSignalScanHistograms(ONZ);
                //s1Maps[masspoint][sigMt2][LL][PDFNAMES[i]]=InitSignalScanHistograms(ONZ);

                //s2Maps[masspoint][sig][EE][PDFNAMES[i]]=InitSignalScanHistograms2D(ONZ);
                //s2Maps[masspoint][sig][MM][PDFNAMES[i]]=InitSignalScanHistograms2D(ONZ);
                //s2Maps[masspoint][sig][EM][PDFNAMES[i]]=InitSignalScanHistograms2D(ONZ);
                //s2Maps[masspoint][sig][LL][PDFNAMES[i]]=InitSignalScanHistograms2D(ONZ);
        }

        if(config_dosignalscanSplit) {
                s1Maps[masspoint][sig_gg][LL][nom]=InitSignalScanHistograms(ONZ);
                s1Maps[masspoint][sig_gg][EE][nom]=InitSignalScanHistograms(ONZ);
                s1Maps[masspoint][sig_gg][MM][nom]=InitSignalScanHistograms(ONZ);
                s1Maps[masspoint][sig_gg][EM][nom]=InitSignalScanHistograms(ONZ);
                s1Maps[masspoint][sig_zz][LL][nom]=InitSignalScanHistograms(ONZ);
                s1Maps[masspoint][sig_zz][EE][nom]=InitSignalScanHistograms(ONZ);
                s1Maps[masspoint][sig_zz][MM][nom]=InitSignalScanHistograms(ONZ);
                s1Maps[masspoint][sig_zz][EM][nom]=InitSignalScanHistograms(ONZ);
                s1Maps[masspoint][sig_gz][LL][nom]=InitSignalScanHistograms(ONZ);
                s1Maps[masspoint][sig_gz][EE][nom]=InitSignalScanHistograms(ONZ);
                s1Maps[masspoint][sig_gz][MM][nom]=InitSignalScanHistograms(ONZ);
                s1Maps[masspoint][sig_gz][EM][nom]=InitSignalScanHistograms(ONZ);
        }
}

//void myAnalyzer::InitWeightHistos(map<string,map<string,map<Histograms1D,TH1F>>>& map_,const selectionType selection_, string name_){
void myAnalyzer::InitWeightHistos(map<selectionFolderName,map<selectionFolderName,map<selectionFolderName,map<Histograms1D,TH1F> > > >& map_,const selectionType selection_, selectionFolderName name_){
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
        //map_[name_][EE][nom]=InitHistogramsSlimmed(selection_);
        //map_[name_][MM][nom]=InitHistogramsSlimmed(selection_);
        map_[name_][EM][nom]=InitHistogramsSlimmed(selection_);
        map_[name_][LL][nom]=InitHistogramsSlimmed(selection_);

        //map_[name_][EE][JESu]=InitHistogramsSlimmed(selection_);
        //map_[name_][MM][JESu]=InitHistogramsSlimmed(selection_);
        map_[name_][EM][JESu]=InitHistogramsSlimmed(selection_);
        map_[name_][LL][JESu]=InitHistogramsSlimmed(selection_);

        //map_[name_][EE][JESd]=InitHistogramsSlimmed(selection_);
        //map_[name_][MM][JESd]=InitHistogramsSlimmed(selection_);
        map_[name_][EM][JESd]=InitHistogramsSlimmed(selection_);
        map_[name_][LL][JESd]=InitHistogramsSlimmed(selection_);

        //map_[name_][EE][JERu]=InitHistogramsSlimmed(selection_);
        //map_[name_][MM][JERu]=InitHistogramsSlimmed(selection_);
        map_[name_][EM][JERu]=InitHistogramsSlimmed(selection_);
        map_[name_][LL][JERu]=InitHistogramsSlimmed(selection_);

        //map_[name_][EE][JERd]=InitHistogramsSlimmed(selection_);
        //map_[name_][MM][JERd]=InitHistogramsSlimmed(selection_);
        map_[name_][EM][JERd]=InitHistogramsSlimmed(selection_);
        map_[name_][LL][JERd]=InitHistogramsSlimmed(selection_);

        //map_[name_][EE][GENMET]=InitHistogramsSlimmed(selection_);
        //map_[name_][MM][GENMET]=InitHistogramsSlimmed(selection_);
        //map_[name_][EM][GENMET]=InitHistogramsSlimmed(selection_);
        //map_[name_][LL][GENMET]=InitHistogramsSlimmed(selection_);

        //map_[name_][EE][LEPSFUP]=InitHistogramsSlimmed(selection_);
        map_[name_][EM][LEPSFUP]=InitHistogramsSlimmed(selection_);
        //map_[name_][MM][LEPSFUP]=InitHistogramsSlimmed(selection_);
        map_[name_][LL][LEPSFUP]=InitHistogramsSlimmed(selection_);

        //map_[name_][EE][LEPSFDOWN]=InitHistogramsSlimmed(selection_);
        map_[name_][EM][LEPSFDOWN]=InitHistogramsSlimmed(selection_);
        //map_[name_][MM][LEPSFDOWN]=InitHistogramsSlimmed(selection_);
        map_[name_][LL][LEPSFDOWN]=InitHistogramsSlimmed(selection_);

        //map_[name_][EE][PHOTONSFUP]=InitHistogramsSlimmed(selection_);
        map_[name_][EM][PHOTONSFUP]=InitHistogramsSlimmed(selection_);
        //map_[name_][MM][PHOTONSFUP]=InitHistogramsSlimmed(selection_);
        map_[name_][LL][PHOTONSFUP]=InitHistogramsSlimmed(selection_);

        //map_[name_][EE][PHOTONSFDOWN]=InitHistogramsSlimmed(selection_);
        map_[name_][EM][PHOTONSFDOWN]=InitHistogramsSlimmed(selection_);
        //map_[name_][MM][PHOTONSFDOWN]=InitHistogramsSlimmed(selection_);
        map_[name_][LL][PHOTONSFDOWN]=InitHistogramsSlimmed(selection_);

        //map_[name_][EE][PUUP]=InitHistogramsSlimmed(selection_);
        map_[name_][EM][PUUP]=InitHistogramsSlimmed(selection_);
        //map_[name_][MM][PUUP]=InitHistogramsSlimmed(selection_);
        map_[name_][LL][PUUP]=InitHistogramsSlimmed(selection_);
        //map_[name_][EE][NOPUUP]=InitHistogramsSlimmed(selection_);
        //map_[name_][EM][NOPUUP]=InitHistogramsSlimmed(selection_);
        //map_[name_][MM][NOPUUP]=InitHistogramsSlimmed(selection_);
        //map_[name_][LL][NOPUUP]=InitHistogramsSlimmed(selection_);
        //map_[name_][EE][NOPU]=InitHistogramsSlimmed(selection_);
        //map_[name_][EM][NOPU]=InitHistogramsSlimmed(selection_);
        //map_[name_][MM][NOPU]=InitHistogramsSlimmed(selection_);
        //map_[name_][LL][NOPU]=InitHistogramsSlimmed(selection_);

        //map_[name_][EE][NOPUDOWN]=InitHistogramsSlimmed(selection_);
        //map_[name_][EM][NOPUDOWN]=InitHistogramsSlimmed(selection_);
        //map_[name_][MM][NOPUDOWN]=InitHistogramsSlimmed(selection_);
        //map_[name_][LL][NOPUDOWN]=InitHistogramsSlimmed(selection_);
        //map_[name_][EE][PUDOWN]=InitHistogramsSlimmed(selection_);
        map_[name_][EM][PUDOWN]=InitHistogramsSlimmed(selection_);
        //map_[name_][MM][PUDOWN]=InitHistogramsSlimmed(selection_);
        map_[name_][LL][PUDOWN]=InitHistogramsSlimmed(selection_);

        //map_[name_][EE][ISRUP]=InitHistogramsSlimmed(selection_);
        map_[name_][EM][ISRUP]=InitHistogramsSlimmed(selection_);
        //map_[name_][MM][ISRUP]=InitHistogramsSlimmed(selection_);
        map_[name_][LL][ISRUP]=InitHistogramsSlimmed(selection_);

        //map_[name_][EE][ISRDOWN]=InitHistogramsSlimmed(selection_);
        map_[name_][EM][ISRDOWN]=InitHistogramsSlimmed(selection_);
        //map_[name_][MM][ISRDOWN]=InitHistogramsSlimmed(selection_);
        map_[name_][LL][ISRDOWN]=InitHistogramsSlimmed(selection_);

        //map_[name_][EE][EWKUP]=InitHistogramsSlimmed(selection_);
        map_[name_][EM][EWKUP]=InitHistogramsSlimmed(selection_);
        //map_[name_][MM][EWKUP]=InitHistogramsSlimmed(selection_);
        map_[name_][LL][EWKUP]=InitHistogramsSlimmed(selection_);

        //map_[name_][EE][EWKDOWN]=InitHistogramsSlimmed(selection_);
        map_[name_][EM][EWKDOWN]=InitHistogramsSlimmed(selection_);
        //map_[name_][MM][EWKDOWN]=InitHistogramsSlimmed(selection_);
        map_[name_][LL][EWKDOWN]=InitHistogramsSlimmed(selection_);


        //for(int i=0; i<nWeights; i++){
        for(int i=config_pdfStart; i<config_pdfStart+config_pdfStep; i++) {
                //map_[name_+"EE"][to_string(i)]=InitHistogramsSlimmed(selection_);
                //map_[name_+"MM"][to_string(i)]=InitHistogramsSlimmed(selection_);
                //map_[name_+"EM"][to_string(i)]=InitHistogramsSlimmed(selection_);
                //map_[name_][to_string(i)]=InitHistogramsSlimmed(selection_);
                //map_[name_][EE][PDFNAMES[i]]=InitHistogramsSlimmed(selection_);
                //map_[name_][MM][PDFNAMES[i]]=InitHistogramsSlimmed(selection_);
                map_[name_][EM][PDFNAMES[i]]=InitHistogramsSlimmed(selection_);
                map_[name_][LL][PDFNAMES[i]]=InitHistogramsSlimmed(selection_);
        }
}



map<Histograms1D,TEfficiency> myAnalyzer::InitTriggerStudies(const selectionType selection_){

        map<Histograms1D,TEfficiency> hMap;

        hMap[ETMISS] = TEfficiency(emptyLabelPtr, ";#it{p}_{T}^{miss} (GeV)", 200, 0, 1000);
        hMap[PT1] = TEfficiency(emptyLabelPtr, ";#it{p}_{T}^{leading} (GeV)", 1000, 0, 1000);
        hMap[PT2] = TEfficiency(emptyLabelPtr, ";#it{p}_{T}^{trailing} (GeV)", 1000, 0, 1000);
        hMap[MLL] = TEfficiency(emptyLabelPtr, ";#it{m}_{ll} (GeV)", 200, 0, 1000);
        hMap[NPHOTONS] = TEfficiency(emptyLabelPtr,"n_#gamma",10,0,10);
        hMap[NVTX] = TEfficiency(emptyLabelPtr,"n_{Vtx}",60,0,60);
        hMap[HT] = TEfficiency(emptyLabelPtr,"#it{H}_{T} (GeV)",200,0,1000);
        hMap[GENHT] = TEfficiency(emptyLabelPtr,"#it{H}_{T}^{gen} (GeV)",200,0,1000);
        hMap[NJETS] = TEfficiency(emptyLabelPtr,"n_{Jets}",20,0,20);
        hMap[ETA1] = TEfficiency(emptyLabelPtr, ";|#eta_{leading}|", 260, 0, 2.6);
        hMap[ETA2] = TEfficiency(emptyLabelPtr, ";|#eta_{trailing}|", 260, 0, 2.6);
        hMap[PHI1] = TEfficiency(emptyLabelPtr, ";|#phi_{leading}|", 350, 0, 3.5);
        hMap[PHI2] = TEfficiency(emptyLabelPtr, ";|#phi_{trailing}|", 350, 0, 3.5);

        if ((selection_==TRIGONZ)||(selection_==TRIGSEL)||(selection_==TRIGSEL_ptcuts)||(selection_==TRIGSEL_ptcuts)) {
                hMap[PTG1] = TEfficiency(emptyLabelPtr, ";#it{p}_{T}^{#gamma 1} (GeV)", 200, 0, 1000);
                hMap[ETAG1] = TEfficiency(emptyLabelPtr, ";|#eta_{#gamma 1}|", 260, 0, 2.6);
                hMap[PHIG1] = TEfficiency(emptyLabelPtr, ";|#phi_{#gamma 1}|", 350, 0, 3.5);
                hMap[SIGMAIETAIETAG1] = TEfficiency(emptyLabelPtr, ";#sigma_{i#etai#eta}^{#gamma 1}", 400, 0, 0.04);
        }

        return hMap;

}



void myAnalyzer::InitTriggerStudiesHistos(){
        eff1Maps[trigdilep][EE]=InitTriggerStudies(TRIGDILEP);
        eff1Maps[trigdilep][MM]=InitTriggerStudies(TRIGDILEP);
        eff1Maps[trigdilep][EM]=InitTriggerStudies(TRIGDILEP);
        eff1Maps[trigdilep_ptcuts][EE]=InitTriggerStudies(TRIGDILEP_ptcuts);
        eff1Maps[trigdilep_ptcuts][MM]=InitTriggerStudies(TRIGDILEP_ptcuts);
        eff1Maps[trigdilep_ptcuts][EM]=InitTriggerStudies(TRIGDILEP_ptcuts);
}


















bool myAnalyzer::CheckParticles(){
        return ((muons->size()>=2.) || (electrons->size()>=2.) || ((muons->size()>=1)&&(electrons->size()>=1)) );
}
bool myAnalyzer::Check2Ele(){
        return (electrons->size()>=2.);
}
bool myAnalyzer::Check2Mu(){
        return (muons->size()>=2.);
}
bool myAnalyzer::CheckEMu(){
        return ((muons->size()>=1.) && (electrons->size()>=1.));
}





void myAnalyzer::clearCutFlowMap(){
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
        decisionMapCutFlowFine[ex2LEP]=false;
        decisionMapCutFlowFine[MT2Cut]=false;
        decisionMapCutFlowFine[MetCut]=false;

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
        decisionMapCutFlowFine_weight[ex2LEP]=totalWeight;
        decisionMapCutFlowFine_weight[MT2Cut]=totalWeight;
        decisionMapCutFlowFine_weight[MetCut]=totalWeight;
}







//ELECTRON
//bool myAnalyzer::testSelection(const tree::Electron& pa,selectionType selection ,bool leading){
bool myAnalyzer::testSelection(const selElectron& pa,selectionType selection,bool leading){
        bool decision=false;
        if (selection==UNCUT) {
                //decision = pa.isPassConvVeto && pa.passImpactParameter && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && pa.isTightMVA;
                decision = pa.isPassConvVeto && pa.passImpactParameter && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && pa.isMedium;
        }
        if(selection==SEL || selection==DILEP || selection==ControlRegionZZ || selection==ControlRegionWZ || selection==SEL) {
                //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (*deltaRll>0.1) && pa.isTightMVA;
                decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (*deltaRll>0.1) && pa.isMedium;

        }
        if(selection==ONZ) {
                //decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (*deltaRll>0.1) && pa.isTightMVA;
                decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (*deltaRll>0.1) && pa.isMedium;
                //if (pa.isPassConvVeto && pa.isTightMVA){
                if (pa.isPassConvVeto && pa.isMedium) {
                        leading ? decisionMapCutFlowFine[LEPTONIDPure_leading]=true : decisionMapCutFlowFine[LEPTONIDPure_trailing]=true;
                        if (pa.passImpactParameter) {
                                leading ? decisionMapCutFlowFine[LEPTONIDImpact_leading]=true : decisionMapCutFlowFine[LEPTONIDImpact_trailing]=true;
                                if((fabs(pa.p.Eta())<2.4)) {
                                        leading ? decisionMapCutFlowFine[LEPTONIDEta_leading]=true : decisionMapCutFlowFine[LEPTONIDEta_trailing]=true;
                                        if(pa.miniIso<0.1) {
                                                leading ? decisionMapCutFlowFine[LEPTONIDIso_leading]=true : decisionMapCutFlowFine[LEPTONIDIso_trailing]=true;
                                                if(*deltaRll>0.1) {
                                                        leading ? decisionMapCutFlowFine[LEPTONIDDeltaR_leading]=true : decisionMapCutFlowFine[LEPTONIDDeltaR_trailing]=true;
                                                        if((leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.))) {
                                                                leading ? decisionMapCutFlowFine[LEPTONPT_leading]=true : decisionMapCutFlowFine[LEPTONPT_trailing]=true;
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }
        if(selection==TRIGSEL || selection==TRIGDILEP || selection==TRIGONZ) {
                //decision =(*ht>200.) && pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>0.) : (pa.p.Pt()>0.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (*deltaRll>0.1) && pa.isTightMVA;
                decision =(*ht>200.) && pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>0.) : (pa.p.Pt()>0.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (*deltaRll>0.1) && pa.isMedium;
                //decision =pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>0.) : (pa.p.Pt()>0.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (*deltaRll>0.1) && pa.isMedium;
                //decision =(*ETmiss>150.)&& pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>0.) : (pa.p.Pt()>0.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (*deltaRll>0.1) && pa.isMedium;
        }
        if(selection==TRIGSEL_ptcuts || selection==TRIGDILEP_ptcuts || selection==TRIGONZ_ptcuts) {
                //decision =(*ht>200.) && pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (*deltaRll>0.1) && pa.isTightMVA;
                decision =(*ht>200.) && pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (*deltaRll>0.1) && pa.isMedium;
                //decision =pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (*deltaRll>0.1) && pa.isMedium;
                decision = (*ETmiss>150.)&& pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (*deltaRll>0.1) && pa.isMedium;
                decision =pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>35.) : (pa.p.Pt()>35.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (*deltaRll>0.1) && pa.isMedium;
        }
        return decision;
}
//MUON
//bool myAnalyzer::testSelection(const tree::Muon& pa, selectionType selection_, bool leading){
bool myAnalyzer::testSelection(const selMuon& pa, selectionType selection_, bool leading){
        bool decision=false;
        if (selection_==UNCUT) {
                //decision = pa.passImpactParameter  && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (pa.isTight);
                decision = pa.passImpactParameter  && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (pa.isMedium);
        }
        if (selection_==SEL || selection_==DILEP || selection_==ControlRegionZZ || selection_==ControlRegionWZ) {
                decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.2) && (pa.isMedium) && (*deltaRll>0.1);
        }
        if (selection_ == ONZ) {
                decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.2) && (pa.isMedium) && (*deltaRll>0.1);
                if((pa.isMedium)) {
                        leading ? decisionMapCutFlowFine[LEPTONIDPure_leading]=true : decisionMapCutFlowFine[LEPTONIDPure_trailing]=true;
                        if(pa.passImpactParameter) {
                                leading ? decisionMapCutFlowFine[LEPTONIDImpact_leading]=true : decisionMapCutFlowFine[LEPTONIDImpact_trailing]=true;
                                if((fabs(pa.p.Eta())<2.4) ) {
                                        leading ? decisionMapCutFlowFine[LEPTONIDEta_leading]=true : decisionMapCutFlowFine[LEPTONIDEta_trailing]=true;
                                        if(pa.miniIso<0.2) {
                                                leading ? decisionMapCutFlowFine[LEPTONIDIso_leading]=true : decisionMapCutFlowFine[LEPTONIDIso_trailing]=true;
                                                if((*deltaRll>0.1)) {
                                                        leading ? decisionMapCutFlowFine[LEPTONIDDeltaR_leading]=true : decisionMapCutFlowFine[LEPTONIDDeltaR_trailing]=true;
                                                        if((leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.))) {
                                                                leading ? decisionMapCutFlowFine[LEPTONPT_leading]=true : decisionMapCutFlowFine[LEPTONPT_trailing]=true;
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }
        if(selection_==TRIGSEL || selection_==TRIGDILEP || selection_==TRIGONZ) {
                //decision = (*ht>200.) && pa.passImpactParameter && (leading ? (pa.p.Pt()>0.) : (pa.p.Pt()>0.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (*deltaRll>0.1);
                decision = (*ht>200.) && pa.passImpactParameter && (leading ? (pa.p.Pt()>0.) : (pa.p.Pt()>0.)) && (fabs(pa.p.Eta())<2.4)  && (pa.miniIso<0.2) && (pa.isMedium) && (*deltaRll>0.1);
                //decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>0.) : (pa.p.Pt()>0.)) && (fabs(pa.p.Eta())<2.4)  && (pa.miniIso<0.2) && (pa.isMedium) && (*deltaRll>0.1);
                //decision = (*ETmiss>150.)&& pa.passImpactParameter && (leading ? (pa.p.Pt()>0.) : (pa.p.Pt()>0.)) && (fabs(pa.p.Eta())<2.4)  && (pa.miniIso<0.2) && (pa.isMedium) && (*deltaRll>0.1);
        }
        if(selection_==TRIGSEL_ptcuts || selection_==TRIGDILEP_ptcuts || selection_==TRIGONZ_ptcuts) {
                //decision = (*ht>200.)&& pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (*deltaRll>0.1);
                decision = (*ht>200.)&& pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4)  && (pa.miniIso<0.2) && (pa.isMedium) && (*deltaRll>0.1);
                //decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4)  && (pa.miniIso<0.2) && (pa.isMedium) && (*deltaRll>0.1);
                //decision = (*ETmiss>150.)&& pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4)  && (pa.miniIso<0.2) && (pa.isMedium) && (*deltaRll>0.1);
                //decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>35.) : (pa.p.Pt()>35.)) && (fabs(pa.p.Eta())<2.4)  && (pa.miniIso<0.2) && (pa.isMedium) && (*deltaRll>0.1);
        }
        return decision;
}
//PHOTON
bool myAnalyzer::testSelection(const selPhoton& pa, selectionType selection){
        bool decision=false;
        if(selection==SEL || selection==DILEP || selection==EGRegression || selection==ControlRegionZZ || selection==TRIGDILEP || selection==TRIGSEL || selection==TRIGONZ) {
                decision = (pa.p.Pt()>20.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isLoose) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3); //study Delta R cut ;
                //decision = (pa.p.Pt()>20.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<2.4) && (pa.isLoose) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3); //study Delta R cut ;
        }
        if(selection==ONZ) {
                //decision = (pa.p.Pt()>20.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<2.4) && (pa.isLoose) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3); //study Delta R cut ;
                decision = (pa.p.Pt()>20.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isLoose) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3); //study Delta R cut ;
                if((pa.isLoose)) {
                        if(selection==ONZ) decisionMapCutFlowFine[PHOTON1ID]=true;
                        if(!(pa.hasPixelSeed)) {
                                if(selection==ONZ) decisionMapCutFlowFine[PHOTON1SEED]=true;
                                if(fabs(pa.p.Eta()<1.4442)) {
                                        //if(fabs(pa.p.Eta()<2.4)){
                                        if(selection==ONZ) decisionMapCutFlowFine[PHOTON1ETA]=true;
                                        //if((pa.p.Pt()>25.)){
                                        if((pa.p.Pt()>20.)) {
                                                if(selection==ONZ) decisionMapCutFlowFine[PHOTON1PT]=true;
                                                if((pa.deltaR1>0.3) && (pa.deltaR2>0.3)) {
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
bool myAnalyzer::testSelection(const selMuon& pa, selectionType selection_){
        bool decision=false;
        if (selection_==LooseLeptons) {
                decision = pa.passImpactParameter && (pa.p.Pt()>20.) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.2) && (pa.isMedium);
        }
        return decision;
}
//SELELECTRON
//bool myAnalyzer::testSelection(const selElectron& pa, selectionType selection_){
//bool decision=false;
//if(selection_==LooseLeptons){
//decision = pa.isPassConvVeto && pa.passImpactParameter && (pa.p.Pt()>20.) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && pa.isTightMVA;
//decision = pa.isPassConvVeto && pa.passImpactParameter && (pa.p.Pt()>20.) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && pa.isMedium;
//}
//return decision;
//}
//JET
bool myAnalyzer::testSelection(const selJet& pa){
        bool decision=false;
        decision = (pa.p.Pt()>30.) && (fabs(pa.p.Eta())<2.4) && (pa.isLoose) && !(pa.hasElectronMatch) && !(pa.hasPhotonMatch) && !(pa.hasMuonMatch); //&& (pa.deltaR1>0.3) && (pa.deltaR2>0.3) //study Delta R cut ;
        return decision;
}






bool myAnalyzer::GenPhotonTTVeto(){
        //bool isGenPhoton=false;
        bool isGenPhoton2=false;
        //bool deltaRGreater02=true;
        float minDeltaR=10.;
        if(!isData) {
                if(isTTInclusive) {
                        for (vector<tree::GenParticle>::iterator it = genParticles->begin(); it != genParticles->end(); it++) {
                                minDeltaR=10.;
                                //if((it->pdgId==22)&&(it->statusID==1)){
                                if((abs(it->pdgId)==22)&&(it->statusID==1)) {
                                        if(it->p.Pt()>13.) {
                                                if(it->p.Eta()<3.0) {
                                                        if(((abs(it->motherId)>0)&&(abs(it->motherId)<7))  || ((abs(it->motherId)>10)&&(abs(it->motherId)<19)) || (abs(it->motherId)==9)|| (abs(it->motherId)==21)|| (abs(it->motherId)==24)|| (abs(it->motherId)==23)) {
                                                                for (vector<tree::GenParticle>::iterator it2 = genParticles->begin(); it2 != genParticles->end(); it2++) {
                                                                        if((it->p.DeltaR(it2->p)<minDeltaR)&&(it->p.DeltaR(it2->p)>0.001)) {
                                                                                minDeltaR=it->p.DeltaR(it2->p);
                                                                        }
                                                                }
                                                                if(minDeltaR>0.2) {
                                                                        isGenPhoton2=isGenPhoton2||true;
                                                                }
                                                        }
                                                }
                                        }
                                }
                                //if((it->pdgId==22)&&(it->statusID==1)){
                                //if(((it->isPromptFinalState))){
                                //isGenPhoton=isGenPhoton||true;
                                //}
                                //}
                        }
                }
        }
        //if((isGenPhoton+isGenPhoton2)>0){
        //if((isGenPhoton2)){
        //cout<<isGenPhoton<<" "<<isGenPhoton2<<endl;
        //}
        //return (isGenPhoton && isTTInclusive);
        return (isGenPhoton2 && isTTInclusive);
        //return (false);
}


float myAnalyzer::GetZZKFactor(){
        //float Zmass=91.1876;
        float k=0.0;
        if(!isData) {
                if(isZZ2L||isZZ4L) {
                        if(intermediateGenParticles->size()==2) {
                                auto Z1=intermediateGenParticles->at(0);
                                auto Z2=intermediateGenParticles->at(1);
                                int finalState=0;
                                //float PtZ1=Z1.p.Pt();
                                //float PtZ2=Z2.p.Pt();
                                float GENpTZZ=(Z1.p+Z2.p).Pt();
                                int flavor1=abs(Z1.daughters.at(0).pdgId);
                                int flavor2=abs(Z2.daughters.at(0).pdgId);
                                if (flavor1==flavor2) {
                                        finalState=1;
                                }else{
                                        finalState=2;
                                }




                                if (finalState==1) {
                                        //cout<<"1 FS"<<endl;
                                        k+=0.64155491983*(abs(GENpTZZ)>0.0&&abs(GENpTZZ)<=5.0);
                                        k+=1.09985240531*(abs(GENpTZZ)>5.0&&abs(GENpTZZ)<=10.0);
                                        k+=1.29390628654*(abs(GENpTZZ)>10.0&&abs(GENpTZZ)<=15.0);
                                        k+=1.37859998571*(abs(GENpTZZ)>15.0&&abs(GENpTZZ)<=20.0);
                                        k+=1.42430263312*(abs(GENpTZZ)>20.0&&abs(GENpTZZ)<=25.0);
                                        k+=1.45038493266*(abs(GENpTZZ)>25.0&&abs(GENpTZZ)<=30.0);
                                        k+=1.47015377651*(abs(GENpTZZ)>30.0&&abs(GENpTZZ)<=35.0);
                                        k+=1.48828685748*(abs(GENpTZZ)>35.0&&abs(GENpTZZ)<=40.0);
                                        k+=1.50573440448*(abs(GENpTZZ)>40.0&&abs(GENpTZZ)<=45.0);
                                        k+=1.50211655928*(abs(GENpTZZ)>45.0&&abs(GENpTZZ)<=50.0);
                                        k+=1.50918720827*(abs(GENpTZZ)>50.0&&abs(GENpTZZ)<=55.0);
                                        k+=1.52463089491*(abs(GENpTZZ)>55.0&&abs(GENpTZZ)<=60.0);
                                        k+=1.52400838378*(abs(GENpTZZ)>60.0&&abs(GENpTZZ)<=65.0);
                                        k+=1.52418067701*(abs(GENpTZZ)>65.0&&abs(GENpTZZ)<=70.0);
                                        k+=1.55424382578*(abs(GENpTZZ)>70.0&&abs(GENpTZZ)<=75.0);
                                        k+=1.52544284222*(abs(GENpTZZ)>75.0&&abs(GENpTZZ)<=80.0);
                                        k+=1.57896384602*(abs(GENpTZZ)>80.0&&abs(GENpTZZ)<=85.0);
                                        k+=1.53034682567*(abs(GENpTZZ)>85.0&&abs(GENpTZZ)<=90.0);
                                        k+=1.56147329708*(abs(GENpTZZ)>90.0&&abs(GENpTZZ)<=95.0);
                                        k+=1.54468169268*(abs(GENpTZZ)>95.0&&abs(GENpTZZ)<=100.0);
                                        k+=1.57222952415*(abs(GENpTZZ)>100.0);
                                }

                                if (finalState==2) {
                                        //cout<<"2 FS"<<endl;
                                        k+=0.743602533303*(abs(GENpTZZ)>0.0&&abs(GENpTZZ)<=5.0);
                                        k+=1.14789453219*(abs(GENpTZZ)>5.0&&abs(GENpTZZ)<=10.0);
                                        k+=1.33815867892*(abs(GENpTZZ)>10.0&&abs(GENpTZZ)<=15.0);
                                        k+=1.41420044104*(abs(GENpTZZ)>15.0&&abs(GENpTZZ)<=20.0);
                                        k+=1.45511318916*(abs(GENpTZZ)>20.0&&abs(GENpTZZ)<=25.0);
                                        k+=1.47569225244*(abs(GENpTZZ)>25.0&&abs(GENpTZZ)<=30.0);
                                        k+=1.49053003693*(abs(GENpTZZ)>30.0&&abs(GENpTZZ)<=35.0);
                                        k+=1.50622827695*(abs(GENpTZZ)>35.0&&abs(GENpTZZ)<=40.0);
                                        k+=1.50328889799*(abs(GENpTZZ)>40.0&&abs(GENpTZZ)<=45.0);
                                        k+=1.52186945281*(abs(GENpTZZ)>45.0&&abs(GENpTZZ)<=50.0);
                                        k+=1.52043468754*(abs(GENpTZZ)>50.0&&abs(GENpTZZ)<=55.0);
                                        k+=1.53977869986*(abs(GENpTZZ)>55.0&&abs(GENpTZZ)<=60.0);
                                        k+=1.53491994434*(abs(GENpTZZ)>60.0&&abs(GENpTZZ)<=65.0);
                                        k+=1.51772882172*(abs(GENpTZZ)>65.0&&abs(GENpTZZ)<=70.0);
                                        k+=1.54494489131*(abs(GENpTZZ)>70.0&&abs(GENpTZZ)<=75.0);
                                        k+=1.57762411697*(abs(GENpTZZ)>75.0&&abs(GENpTZZ)<=80.0);
                                        k+=1.55078339014*(abs(GENpTZZ)>80.0&&abs(GENpTZZ)<=85.0);
                                        k+=1.57078191891*(abs(GENpTZZ)>85.0&&abs(GENpTZZ)<=90.0);
                                        k+=1.56162666568*(abs(GENpTZZ)>90.0&&abs(GENpTZZ)<=95.0);
                                        k+=1.54183774627*(abs(GENpTZZ)>95.0&&abs(GENpTZZ)<=100.0);
                                        k+=1.58485762205*(abs(GENpTZZ)>100.0);
                                }



                        }
                }
        }
        float kToReturn=1.;
        if (k==0.0) {
                kToReturn= 1.1;
        }else{
                kToReturn= k; // if something goes wrong return inclusive k-factor
        }

        if(isZZ2L||isZZ4L) {
                return kToReturn;
        }else{
                return 1.;
        }

}









bool myAnalyzer::SelectEvent(selectionType selection){
        clearCutFlowMap();
        if(*selLeptonSize>1) {
                //tempWeight = tempWeight*puWeights[(int)changePU]*lepSfWeights[(int)changeLepSF] *photonSfWeights[(int)changePhotonSF] **topPt_weight*isrWeights[(int)changeisr]*ewkWeights[(int)changeewk]*weight_PDF;

                //float totalWeightCalc=totalWeight;
                float totalWeightCalc=*mc_weight**pu_weight**ewk_weight**isr_weight**topPt_weight;
                //float totalWeightCalc=totalWeight**pu_weight**ewk_weight**isr_weight;
                //fill cutflowInfo DiEle or DiMu?
                if((selection==ONZ) && (!isTotalSignal)) {
                        cutflowDiEle=*isDiElectron;
                        cutflowDiMu=*isDiMuon;
                        decisionMapCutFlowFine[DIMUON]=*isDiMuon;
                        decisionMapCutFlowFine[DIELECTRON]=*isDiElectron;
                        decisionMapCutFlowFine[EMUON]=*isMuonElectron||*isElectronMuon;
                }
                int genNNToZ=0;
                int genNNToG=0;
                //int genNNToH=0;
                int genZToLL=0;

                if((selection==ONZ) && (!isTotalSignal)) {
                        //get "triggered" signal events for cutflow
                        for(auto& intermediate: *intermediateGenParticles) {
                                for(auto& daughter: intermediate.daughters) {
                                        if((abs(daughter.pdgId)==23)&&(intermediate.pdgId>1000021)) {
                                                genNNToZ=genNNToZ+1;
                                        }
                                        if((abs(daughter.pdgId)==22)&&(intermediate.pdgId>1000021)) {
                                                genNNToG=genNNToG+1;
                                        }
                                        //if((abs(daughter.pdgId)==25)&&(intermediate.pdgId>1000021)){
                                        //genNNToH=genNNToH+1;
                                        //}
                                        if(((abs(daughter.pdgId)==11)||(abs(daughter.pdgId)==13))&&(intermediate.pdgId==23)) {
                                                genZToLL=genZToLL+1;
                                        }
                                }
                        }
                }
                if(isSignal) {
                        if(*isDiElectron) totalWeightCalc=totalWeightCalc*config_eeEff;
                        if(*isMuonElectron||*isElectronMuon) totalWeightCalc=totalWeightCalc*config_emEff;
                        if(*isDiMuon) totalWeightCalc=totalWeightCalc*config_mmEff;
                }else{
                        if(!isData) {
                                if(*isDiElectron) totalWeightCalc=totalWeightCalc*config_eeWeight;
                                if(*isMuonElectron||*isElectronMuon) totalWeightCalc=totalWeightCalc*config_emWeight;
                                if(*isDiMuon) totalWeightCalc=totalWeightCalc*config_mmWeight;
                        }
                }

                if ((*isDiElectron)&&(*selElectronSize>1)) {
                        //Fill cutflow triggered bool
                        if((selection==ONZ) && (!isTotalSignal)) {
                                cutflowIsTriggered=true;
                                if(isSignal) {
                                        //decisionMapCutFlowFine[TRIGGERED]=((genNNToG==1)&&(genNNToZ==1)&&(genZToLL==2));
                                        //decisionMapCutFlowFine[TRIGGEREDMATCHED]=((genNNToG==1)&&(genNNToZ==1)&&(genZToLL==2));
                                        decisionMapCutFlowFine[TRIGGERED]=true;
                                        decisionMapCutFlowFine[TRIGGEREDMATCHED]=true;
                                }else{

                                        decisionMapCutFlowFine[TRIGGERED]=*trigDiEle;
                                        decisionMapCutFlowFine[TRIGGEREDMATCHED]=*trigDiEleMatch;
                                }
                                //decisionMapCutFlowFine_weight[TRIGGERED]=totalWeight;
                                decisionMapCutFlowFine_weight[TRIGGERED]=totalWeightCalc;
                                decisionMapCutFlowFine_weight[TRIGGEREDMATCHED]=totalWeightCalc;
                        }

                        auto e1 = selElectrons->at(0);
                        auto e2 = selElectrons->at(1);


                        if(*trigDiEleMatch) {
                                if((!*evtHasGenPhotonVeto)&&(!GenPhotonTTVeto())) {
                                        if((selection==ONZ) && (!isTotalSignal)) {
                                                decisionMapCutFlowFine[GENVETO]=true;
                                                //decisionMapCutFlowFine_weight[GENVETO]=totalWeight;
                                                decisionMapCutFlowFine_weight[GENVETO]=totalWeightCalc;
                                        }

                                        if(testSelection(e1,selection,true) && testSelection(e2,selection,false) && (*chargeProduct < 0.) ) {


                                                //totalWeightCalc=totalWeight* *lepSF_weight;
                                                //totalWeightCalc=totalWeight**pu_weight**ewk_weight**isr_weight* *lepSF_weight;
                                                //totalWeightCalc=totalWeightCalc* *lepSF_weight;

                                                //Fill cutflow 2 leptons bool
                                                if((selection==ONZ) && (!isTotalSignal)) cutflow2Leptons=true;
                                                if((selection==UNCUT) ? true : (*mll>50.)) {
                                                        //Fill cutflow mll>50 bool
                                                        if((selection==ONZ) && (!isTotalSignal)) {
                                                                cutflowMll50=true;
                                                                decisionMapCutFlowFine[M50]=true;
                                                                decisionMapCutFlowFine_weight[LEPTONIDPure_leading]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDPure_trailing]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDImpact_leading]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDImpact_trailing]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDEta_leading]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDEta_trailing]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDIso_leading]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDIso_trailing]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDDeltaR_leading]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDDeltaR_trailing]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONPT_leading]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONPT_trailing]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[M50]=totalWeightCalc* *lepSF_weight;

                                                                if((photons->size()!=0)) {
                                                                        decisionMapCutFlowFine[PHOTON1]=true;
                                                                        //decisionMapCutFlowFine_weight[PHOTON1]=totalWeightCalc* *lepSF_weight;
                                                                        decisionMapCutFlowFine_weight[PHOTON1]=totalWeightCalc;
                                                                }
                                                        }
                                                        //if(!isData){
                                                        //totalWeightCalc=totalWeightCalc* *lepSF_weight* *photonSF_weight;
                                                        //}
                                                        if((selection==ONZ) && (!isTotalSignal)) {
                                                                if(*selPhotonSize!=0) {
                                                                        bool bla=testSelection(selPhotons->at(0),selection);
                                                                }
                                                                decisionMapCutFlowFine_weight[PHOTON1ID]=totalWeightCalc* *lepSF_weight* *photonSF_weight;
                                                                decisionMapCutFlowFine_weight[PHOTON1SEED]=totalWeightCalc* *lepSF_weight* *photonSF_weight;
                                                                decisionMapCutFlowFine_weight[PHOTON1ETA]=totalWeightCalc* *lepSF_weight* *photonSF_weight;
                                                                decisionMapCutFlowFine_weight[PHOTON1PT]=totalWeightCalc* *lepSF_weight* *photonSF_weight;
                                                                decisionMapCutFlowFine_weight[PHOTON1DR]=totalWeightCalc* *lepSF_weight* *photonSF_weight;

                                                                if(*selPhotonSize!=0) {
                                                                        if(selection==ONZ) cutflow1Photon=true;
                                                                        if((*mll>81.) && (*mll<101.)) {
                                                                                cutflowOnZ=true;
                                                                                decisionMapCutFlowFine[ZMASS]=true;
                                                                                decisionMapCutFlowFine_weight[ZMASS]=totalWeightCalc* *lepSF_weight* *photonSF_weight;
                                                                                if((*selLeptonSize==*matchedLeptonSize && *selLeptonSize==2)) {
                                                                                        decisionMapCutFlowFine[ex2LEP]=true;
                                                                                        decisionMapCutFlowFine_weight[ex2LEP]=totalWeightCalc* *lepSF_weight* *photonSF_weight;
                                                                                        if(!isData) {
                                                                                                if(*Mt2>100.) {
                                                                                                        decisionMapCutFlowFine[MT2Cut]=true;
                                                                                                        decisionMapCutFlowFine_weight[MT2Cut]=totalWeightCalc* *lepSF_weight* *photonSF_weight;
                                                                                                        if(*ETmiss>150.) {
                                                                                                                decisionMapCutFlowFine[MetCut]=true;
                                                                                                                decisionMapCutFlowFine_weight[MetCut]=totalWeightCalc* *lepSF_weight* *photonSF_weight;
                                                                                                        }
                                                                                                }
                                                                                        }
                                                                                }
                                                                        }
                                                                }
                                                        }
                                                        return true;
                                                }
                                        }
                                }
                        }
                }
                if((*isDiMuon)&&(*selMuonSize>1)) {

                        if((selection==ONZ) && (!isTotalSignal)) if(selection==ONZ) decisionMapCutFlowFine_weight[genZLL]=totalWeight;


                        if((selection==ONZ) && (!isTotalSignal)) {
                                cutflowIsTriggered=true;
                                if(isSignal) {
                                        //decisionMapCutFlowFine[TRIGGERED]=((genNNToG==1)&&(genNNToZ==1)&&(genZToLL==2));
                                        //decisionMapCutFlowFine[TRIGGEREDMATCHED]=((genNNToG==1)&&(genNNToZ==1)&&(genZToLL==2));
                                        decisionMapCutFlowFine[TRIGGERED]=true;
                                        decisionMapCutFlowFine[TRIGGEREDMATCHED]=true;
                                }else{
                                        decisionMapCutFlowFine[TRIGGERED]=*trigDiMu;
                                        decisionMapCutFlowFine[TRIGGEREDMATCHED]=*trigDiMuMatch;
                                }
                                //decisionMapCutFlowFine_weight[TRIGGERED]=totalWeight;
                                //decisionMapCutFlowFine_weight[TRIGGEREDMATCHED]=totalWeight;
                                //decisionMapCutFlowFine_weight[TRIGGERED]=totalWeight**pu_weight**ewk_weight**isr_weight;
                                //decisionMapCutFlowFine_weight[TRIGGEREDMATCHED]=totalWeight**pu_weight**ewk_weight**isr_weight;
                                decisionMapCutFlowFine_weight[TRIGGERED]=totalWeightCalc;
                                decisionMapCutFlowFine_weight[TRIGGEREDMATCHED]=totalWeightCalc;
                        }
                        auto m1 = selMuons->at(0);
                        auto m2 = selMuons->at(1);

                        if(*trigDiMuMatch) {
                                if((!*evtHasGenPhotonVeto)&&(!GenPhotonTTVeto())) {

                                        if((selection==ONZ) && (!isTotalSignal)) {
                                                decisionMapCutFlowFine[GENVETO]=true;
                                                //decisionMapCutFlowFine_weight[GENVETO]=totalWeight;
                                                //decisionMapCutFlowFine_weight[GENVETO]=totalWeight**pu_weight**ewk_weight**isr_weight;
                                                decisionMapCutFlowFine_weight[GENVETO]=totalWeightCalc;
                                        }


                                        if(testSelection(m1,selection,true) && testSelection(m2,selection,false) && (*chargeProduct < 0.)) {

                                                //totalWeightCalc=totalWeight* *lepSF_weight;
                                                //totalWeightCalc=totalWeight**pu_weight**ewk_weight**isr_weight* *lepSF_weight;

                                                if((selection==ONZ) && (!isTotalSignal)) cutflow2Leptons=true;
                                                if((selection==UNCUT) ? true : (*mll>50.)) {
                                                        if((selection==ONZ) && (!isTotalSignal)) {
                                                                cutflowMll50=true;
                                                                decisionMapCutFlowFine[M50]=true;
                                                                decisionMapCutFlowFine_weight[LEPTONIDPure_leading]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDPure_trailing]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDImpact_leading]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDImpact_trailing]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDEta_leading]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDEta_trailing]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDIso_leading]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDIso_trailing]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDDeltaR_leading]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDDeltaR_trailing]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONPT_leading]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONPT_trailing]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[M50]=totalWeightCalc* *lepSF_weight;

                                                                if(photons->size()!=0) {
                                                                        decisionMapCutFlowFine[PHOTON1]=true;
                                                                        decisionMapCutFlowFine_weight[PHOTON1]=totalWeightCalc* *lepSF_weight;
                                                                }
                                                        }
                                                        //if(!isData){
                                                        //totalWeightCalc=totalWeightCalc* *photonSF_weight;
                                                        //}
                                                        if((selection==ONZ) && (!isTotalSignal)) {
                                                                if(*selPhotonSize!=0) {
                                                                        bool bla=testSelection(selPhotons->at(0),selection);
                                                                }
                                                                decisionMapCutFlowFine_weight[PHOTON1ID]=totalWeightCalc* *lepSF_weight * *photonSF_weight;
                                                                decisionMapCutFlowFine_weight[PHOTON1SEED]=totalWeightCalc* *lepSF_weight * *photonSF_weight;
                                                                decisionMapCutFlowFine_weight[PHOTON1ETA]=totalWeightCalc* *lepSF_weight * *photonSF_weight;
                                                                decisionMapCutFlowFine_weight[PHOTON1PT]=totalWeightCalc* *lepSF_weight * *photonSF_weight;
                                                                decisionMapCutFlowFine_weight[PHOTON1DR]=totalWeightCalc* *lepSF_weight * *photonSF_weight;

                                                                if(*selPhotonSize>0) {
                                                                        cutflow1Photon=true;
                                                                        if(*mll<101. && *mll>81.) {
                                                                                cutflowOnZ=true;
                                                                                decisionMapCutFlowFine[ZMASS]=true;
                                                                                decisionMapCutFlowFine_weight[ZMASS]=totalWeightCalc* *lepSF_weight * *photonSF_weight;
                                                                                if((*selLeptonSize==*matchedLeptonSize && *selLeptonSize==2)) {
                                                                                        decisionMapCutFlowFine[ex2LEP]=true;
                                                                                        decisionMapCutFlowFine_weight[ex2LEP]=totalWeightCalc* *lepSF_weight* *photonSF_weight;
                                                                                        if(!isData) {
                                                                                                if(*Mt2>=100.) {
                                                                                                        decisionMapCutFlowFine[MT2Cut]=true;
                                                                                                        decisionMapCutFlowFine_weight[MT2Cut]=totalWeightCalc* *lepSF_weight* *photonSF_weight;
                                                                                                        if(*ETmiss>=150.) {
                                                                                                                decisionMapCutFlowFine[MetCut]=true;
                                                                                                                decisionMapCutFlowFine_weight[MetCut]=totalWeightCalc* *lepSF_weight* *photonSF_weight;
                                                                                                        }
                                                                                                }
                                                                                        }
                                                                                }
                                                                        }
                                                                }
                                                        }
                                                        return true;
                                                }
                                        }
                                }
                        }
                }

                if((*isElectronMuon || *isMuonElectron)&&(*selElectronSize>0 && *selMuonSize>0)) {

                        if((selection==ONZ) && (!isTotalSignal)) {
                                cutflowIsTriggered=true;
                                if(isSignal) {
                                        //decisionMapCutFlowFine[TRIGGERED]=((genNNToG==1)&&(genNNToZ==1)&&(genZToLL==2));
                                        //decisionMapCutFlowFine[TRIGGEREDMATCHED]=((genNNToG==1)&&(genNNToZ==1)&&(genZToLL==2));
                                        decisionMapCutFlowFine[TRIGGERED]=true;
                                        decisionMapCutFlowFine[TRIGGEREDMATCHED]=true;
                                }else{
                                        decisionMapCutFlowFine[TRIGGERED]=*trigMuEle;
                                        decisionMapCutFlowFine[TRIGGEREDMATCHED]=*trigMuEleMatch;
                                }
                                //decisionMapCutFlowFine_weight[TRIGGERED]=totalWeight;
                                //decisionMapCutFlowFine_weight[TRIGGEREDMATCHED]=totalWeight;
                                //decisionMapCutFlowFine_weight[TRIGGERED]=totalWeight**pu_weight**ewk_weight**isr_weight;
                                //decisionMapCutFlowFine_weight[TRIGGEREDMATCHED]=totalWeight**pu_weight**ewk_weight**isr_weight;
                                decisionMapCutFlowFine_weight[TRIGGERED]=totalWeightCalc;
                                decisionMapCutFlowFine_weight[TRIGGEREDMATCHED]=totalWeightCalc;
                        }

                        auto e0 = selElectrons->at(0);
                        auto m0 = selMuons->at(0);

                        if(*trigMuEleMatch) {
                                if((!*evtHasGenPhotonVeto)&&(!GenPhotonTTVeto())) {

                                        if((selection==ONZ) && (!isTotalSignal)) {
                                                decisionMapCutFlowFine[GENVETO]=true;
                                                //decisionMapCutFlowFine_weight[GENVETO]=totalWeight;
                                                //decisionMapCutFlowFine_weight[GENVETO]=totalWeight**pu_weight**ewk_weight**isr_weight;
                                                decisionMapCutFlowFine_weight[GENVETO]=totalWeightCalc;
                                        }

                                        if(testSelection(m0,selection,*isMuonElectron) && testSelection(e0,selection,*isElectronMuon) && (*chargeProduct < 0.)) {

                                                //totalWeightCalc=totalWeight* *lepSF_weight;
                                                //totalWeightCalc=totalWeight**pu_weight**ewk_weight**isr_weight* *lepSF_weight;

                                                if((selection==ONZ) && (!isTotalSignal)) cutflow2Leptons=true;
                                                if((selection==UNCUT) ? true : (*mll>50.)) {
                                                        if((selection==ONZ) && (!isTotalSignal)) {
                                                                cutflowMll50=true;
                                                                decisionMapCutFlowFine[M50]=true;
                                                                decisionMapCutFlowFine_weight[LEPTONIDPure_leading]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDPure_trailing]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDImpact_leading]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDImpact_trailing]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDEta_leading]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDEta_trailing]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDIso_leading]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDIso_trailing]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDDeltaR_leading]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONIDDeltaR_trailing]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONPT_leading]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[LEPTONPT_trailing]=totalWeightCalc* *lepSF_weight;
                                                                decisionMapCutFlowFine_weight[M50]=totalWeightCalc* *lepSF_weight;

                                                                if(photons->size()!=0) {
                                                                        decisionMapCutFlowFine[PHOTON1]=true;
                                                                        decisionMapCutFlowFine_weight[PHOTON1]=totalWeightCalc* *lepSF_weight;
                                                                }
                                                        }
                                                        //if(!isData){
                                                        //totalWeightCalc=totalWeightCalc* *photonSF_weight;
                                                        //}
                                                        if((selection==ONZ) && (!isTotalSignal)) {
                                                                if(*selPhotonSize!=0) {
                                                                        bool bla=testSelection(selPhotons->at(0),selection);
                                                                }
                                                                decisionMapCutFlowFine_weight[PHOTON1ID]=totalWeightCalc* *lepSF_weight* *photonSF_weight;
                                                                decisionMapCutFlowFine_weight[PHOTON1SEED]=totalWeightCalc* *lepSF_weight* *photonSF_weight;
                                                                decisionMapCutFlowFine_weight[PHOTON1ETA]=totalWeightCalc* *lepSF_weight* *photonSF_weight;
                                                                decisionMapCutFlowFine_weight[PHOTON1PT]=totalWeightCalc* *lepSF_weight* *photonSF_weight;
                                                                decisionMapCutFlowFine_weight[PHOTON1DR]=totalWeightCalc* *lepSF_weight* *photonSF_weight;

                                                                if(*selPhotonSize>0) {
                                                                        cutflow1Photon=true;
                                                                        if(*mll<101. && *mll>81.) {
                                                                                cutflowOnZ=true;
                                                                                decisionMapCutFlowFine[ZMASS]=true;
                                                                                decisionMapCutFlowFine_weight[ZMASS]=totalWeightCalc* *lepSF_weight* *photonSF_weight;
                                                                                if((*selLeptonSize==*matchedLeptonSize && *selLeptonSize==2)) {
                                                                                        decisionMapCutFlowFine[ex2LEP]=true;
                                                                                        decisionMapCutFlowFine_weight[ex2LEP]=totalWeightCalc* *lepSF_weight* *photonSF_weight;
                                                                                        if(!isData) {
                                                                                                if(*Mt2>100.) {
                                                                                                        decisionMapCutFlowFine[MT2Cut]=true;
                                                                                                        decisionMapCutFlowFine_weight[MT2Cut]=totalWeightCalc* *lepSF_weight* *photonSF_weight;
                                                                                                        if(*ETmiss>150.) {
                                                                                                                decisionMapCutFlowFine[MetCut]=true;
                                                                                                                decisionMapCutFlowFine_weight[MetCut]=totalWeightCalc* *lepSF_weight* *photonSF_weight;
                                                                                                        }
                                                                                                }
                                                                                        }
                                                                                }
                                                                        }
                                                                }
                                                        }
                                                        return true;
                                                }
                                        }
                                }
                        }
                }
        }else{return false;}
        return false;
}

bool myAnalyzer::SelectEventTriggerStudies(selectionType selection){
        if ((*isDiElectron)&&(*selElectronSize>1)) {
                auto e1 = selElectrons->at(0);
                auto e2 = selElectrons->at(1);
                //if(! *evtHasGenPhotonVeto){
                if((!*evtHasGenPhotonVeto)&&(!GenPhotonTTVeto())) {
                        if(testSelection(e1,selection,true) && testSelection(e2,selection,false) && (*chargeProduct < 0.) && ((selection==UNCUT) ? true : (*mll>50.)) ) {
                                //if(testSelection(e1,selection,true) && testSelection(e2,selection,false) && (*chargeProduct < 0.) && ((selection==UNCUT)? true : (*mll>20.)) ){
                                return true;
                        }
                }
        }
        else{
                if ((*isDiMuon)&&(*selMuonSize>1)) {
                        auto m1 = selMuons->at(0);
                        auto m2 = selMuons->at(1);
                        //if(! *evtHasGenPhotonVeto){
                        if((!*evtHasGenPhotonVeto)&&(!GenPhotonTTVeto())) {
                                if(testSelection(m1,selection,true) && testSelection(m2,selection,false) && (*chargeProduct < 0.) && ((selection==UNCUT) ? true : (*mll>50.))) {
                                        //if(testSelection(m1,selection,true) && testSelection(m2,selection,false) && (*chargeProduct < 0.) && ((selection==UNCUT)? true : (*mll>20.))){
                                        return true;
                                }
                        }
                }else{
                        if ((*isElectronMuon || *isMuonElectron)&&(*selElectronSize>0 && *selMuonSize>0)) {
                                auto e0 = selElectrons->at(0);
                                auto m0 = selMuons->at(0);
                                //if(! *evtHasGenPhotonVeto){
                                if((!*evtHasGenPhotonVeto)&&(!GenPhotonTTVeto())) {
                                        if(testSelection(e0,selection,*isElectronMuon) && testSelection(m0,selection,*isMuonElectron) && (*chargeProduct < 0.) && ((selection==UNCUT) ? true : (*mll>50.))) {
                                                //if(testSelection(e0,selection,*isElectronMuon) && testSelection(m0,selection,*isMuonElectron) && (*chargeProduct < 0.) && ((selection==UNCUT)? true : (*mll>20.))){
                                                return true;
                                        }
                                }
                        }
                }
        }
        return false;
}



void myAnalyzer::FillCutFlowHistograms(){
        if(cutflowDiEle && !(cutflowDiMu)) {
                if (cutflowIsTriggered) c1Maps[cutFlow_onZEE].at(CUTFLOW).Fill("triggered",totalWeight);
                if (cutflow2Leptons) c1Maps[cutFlow_onZEE].at(CUTFLOW).Fill("2leptons",totalWeight);
                if (cutflowMll50) c1Maps[cutFlow_onZEE].at(CUTFLOW).Fill("m50",totalWeight);
                if (cutflow1Photon) c1Maps[cutFlow_onZEE].at(CUTFLOW).Fill("1photon",totalWeight);
                if (cutflowOnZ) c1Maps[cutFlow_onZEE].at(CUTFLOW).Fill("Z",totalWeight);
        }
        if(cutflowDiMu && !(cutflowDiEle)) {
                if (cutflowIsTriggered) c1Maps[cutFlow_onZMM].at(CUTFLOW).Fill("triggered",totalWeight);
                if (cutflow2Leptons) c1Maps[cutFlow_onZMM].at(CUTFLOW).Fill("2leptons",totalWeight);
                if (cutflowMll50) c1Maps[cutFlow_onZMM].at(CUTFLOW).Fill("m50",totalWeight);
                if (cutflow1Photon) c1Maps[cutFlow_onZMM].at(CUTFLOW).Fill("1photon",totalWeight);
                if (cutflowOnZ) c1Maps[cutFlow_onZMM].at(CUTFLOW).Fill("Z",totalWeight);
        }

}
void myAnalyzer::SetCutFlowHistogramsStatus(){
        cutflowIsTriggered=false;
        cutflow2Leptons=false;
        cutflowMll50=false;
        cutflow1Photon=false;
        cutflowOnZ=false;
}


void myAnalyzer::FillCutFlowHistograms_Fine(){
        auto m=decisionMapCutFlowFine;
        auto mW=decisionMapCutFlowFine_weight;
        float TChiNG_weight=1.;
        bool useTChiNG=true;
        bool isTChiNG=false;
        float tempWeight=1.;
        if(inputName.find("TChiNG") != string::npos) {
                isTChiNG=true;
                //cout<<"here"<<endl;
                //cout<<(*intermediateGenParticles)[0].daughters.size()<<(*intermediateGenParticles)[1].daughters.size()<<endl;
                //cout<<(*intermediateGenParticles)[0].daughters[0].pdgId<<(*intermediateGenParticles)[0].daughters[1].pdgId<<endl;
                //cout<<(*intermediateGenParticles)[1].daughters[0].pdgId<<(*intermediateGenParticles)[1].daughters[1].pdgId<<endl;
                //GG case -> weight=1. 0.25to0.25
                bool oneOfThem=false;
                if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2) {
                        if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 22) {
                                //cout<<"buh"<<endl;
                                if (fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 22) {
                                        //cout<<"sdfdfsd"<<endl;
                                        useTChiNG=true;
                                        TChiNG_weight=1.;
                                        oneOfThem=oneOfThem||true;

                                }
                        }
                }
                //ZZ case ->weight=4. 1/16 to 0.25
                if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2) {
                        if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 23) {
                                if (fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 23) {
                                        useTChiNG=true;
                                        TChiNG_weight=4.;
                                        oneOfThem=oneOfThem||true;

                                }
                        }
                }
                //cout<<"test"<<endl;
                //ZG case -> weight = 2 0.25to0.5
                if ((*intermediateGenParticles)[0].daughters.size()==2 && (*intermediateGenParticles)[1].daughters.size()==2) {
                        //cout<<"test 2"<<endl;
                        if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 23 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 23) {
                                //cout<<"test3"<<endl;
                                if (fabs((*intermediateGenParticles)[0].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[0].daughters[1].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[0].pdgId) == 22 || fabs((*intermediateGenParticles)[1].daughters[1].pdgId) == 22) {
                                        useTChiNG=true;
                                        TChiNG_weight=2.;
                                        oneOfThem=oneOfThem||true;
                                }
                        }
                }
                if(oneOfThem==false) {
                        useTChiNG=false;
                        TChiNG_weight=1.;
                }



        }

        //cout<<isTChiNG<<useTChiNG<<TChiNG_weight<<endl;

        bool usePoint=(isTChiNG && useTChiNG) || (!isTChiNG);
        if(usePoint) {
                tempWeight=1.;
                if(isTChiNG) {
                        tempWeight=tempWeight*TChiNG_weight;
                }
                //cout<<tempWeight<<endl;
                //cout<<isTChiNG<<useTChiNG<<TChiNG_weight<<tempWeight<<endl;
                if(m[DIELECTRON] && !(m[DIMUON]||m[EMUON])) {

                        if (m[TRIGGERED]) {
                                c1Maps[cutFlow_Fine_onZEE].at(CUTFLOW_fine).Fill("triggered",mW[TRIGGERED]*tempWeight);
                                if(m[TRIGGEREDMATCHED]) {
                                        c1Maps[cutFlow_Fine_onZEE].at(CUTFLOW_fine).Fill("triggeredMatched",mW[TRIGGEREDMATCHED]*tempWeight);
                                        if(m[GENVETO]) {
                                                c1Maps[cutFlow_Fine_onZEE].at(CUTFLOW_fine).Fill("GenPhotonVeto",mW[GENVETO]*tempWeight);
                                                if(m[LEPTONIDPure_leading] && m[LEPTONIDPure_trailing]) {
                                                        c1Maps[cutFlow_Fine_onZEE].at(CUTFLOW_fine).Fill("2LeptonIDPure",mW[LEPTONIDPure_leading]*tempWeight);
                                                        if(m[LEPTONIDImpact_leading] && m[LEPTONIDImpact_trailing]) {
                                                                c1Maps[cutFlow_Fine_onZEE].at(CUTFLOW_fine).Fill("2LeptonIDImpact",mW[LEPTONIDImpact_leading]*tempWeight);
                                                                if(m[LEPTONIDEta_leading] && m[LEPTONIDEta_trailing]) {
                                                                        c1Maps[cutFlow_Fine_onZEE].at(CUTFLOW_fine).Fill("2LeptonIDEta",mW[LEPTONIDEta_leading]*tempWeight);
                                                                        if(m[LEPTONIDIso_leading] && m[LEPTONIDIso_trailing]) {
                                                                                c1Maps[cutFlow_Fine_onZEE].at(CUTFLOW_fine).Fill("2LeptonIDIso",mW[LEPTONIDIso_leading]*tempWeight);
                                                                                if(m[LEPTONIDDeltaR_leading] && m[LEPTONIDDeltaR_trailing]) {
                                                                                        c1Maps[cutFlow_Fine_onZEE].at(CUTFLOW_fine).Fill("2LeptonIDDeltaR",mW[LEPTONIDDeltaR_leading]*tempWeight);
                                                                                        if(m[LEPTONPT_leading] && m[LEPTONPT_trailing]) {
                                                                                                c1Maps[cutFlow_Fine_onZEE].at(CUTFLOW_fine).Fill("2LeptonPT",mW[LEPTONPT_leading]*tempWeight);
                                                                                                if(m[M50]) {
                                                                                                        c1Maps[cutFlow_Fine_onZEE].at(CUTFLOW_fine).Fill("m50",mW[M50]*tempWeight);
                                                                                                        if(m[PHOTON1]) {
                                                                                                                c1Maps[cutFlow_Fine_onZEE].at(CUTFLOW_fine).Fill("1Photon",mW[PHOTON1]*tempWeight);
                                                                                                                if(m[PHOTON1ID]) {
                                                                                                                        c1Maps[cutFlow_Fine_onZEE].at(CUTFLOW_fine).Fill("1PhotonID",mW[PHOTON1ID]*tempWeight);
                                                                                                                        if(m[PHOTON1SEED]) {
                                                                                                                                c1Maps[cutFlow_Fine_onZEE].at(CUTFLOW_fine).Fill("1PhotonSeed",mW[PHOTON1SEED]*tempWeight);
                                                                                                                                if(m[PHOTON1ETA]) {
                                                                                                                                        c1Maps[cutFlow_Fine_onZEE].at(CUTFLOW_fine).Fill("1PhotonEta",mW[PHOTON1ETA]*tempWeight);
                                                                                                                                        if(m[PHOTON1PT]) {
                                                                                                                                                c1Maps[cutFlow_Fine_onZEE].at(CUTFLOW_fine).Fill("1PhotonPT",mW[PHOTON1PT]*tempWeight);
                                                                                                                                                if(m[PHOTON1DR]) {
                                                                                                                                                        c1Maps[cutFlow_Fine_onZEE].at(CUTFLOW_fine).Fill("1PhotonDeltaR",mW[PHOTON1DR]*tempWeight);
                                                                                                                                                        if(m[ZMASS]) {
                                                                                                                                                                c1Maps[cutFlow_Fine_onZEE].at(CUTFLOW_fine).Fill("Z",mW[ZMASS]*tempWeight);
                                                                                                                                                                if(m[ex2LEP]) {
                                                                                                                                                                        c1Maps[cutFlow_Fine_onZEE].at(CUTFLOW_fine).Fill("ex2Lep",mW[ex2LEP]*tempWeight);
                                                                                                                                                                        if(m[MT2Cut]) {
                                                                                                                                                                                c1Maps[cutFlow_Fine_onZEE].at(CUTFLOW_fine).Fill("MT2",mW[MT2Cut]*tempWeight);
                                                                                                                                                                                if(m[MetCut]) {
                                                                                                                                                                                        c1Maps[cutFlow_Fine_onZEE].at(CUTFLOW_fine).Fill("Met",mW[MetCut]*tempWeight);
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
                        }
                }
                if(m[DIMUON] && !(m[DIELECTRON]||m[EMUON])) {
                        if (m[TRIGGERED]) {
                                c1Maps[cutFlow_Fine_onZMM].at(CUTFLOW_fine).Fill("triggered",mW[TRIGGERED]*tempWeight);
                                if(m[TRIGGEREDMATCHED]) {
                                        c1Maps[cutFlow_Fine_onZMM].at(CUTFLOW_fine).Fill("triggeredMatched",mW[TRIGGEREDMATCHED]*tempWeight);
                                        if(m[GENVETO]) {
                                                c1Maps[cutFlow_Fine_onZMM].at(CUTFLOW_fine).Fill("GenPhotonVeto",mW[GENVETO]*tempWeight);
                                                if(m[LEPTONIDPure_leading] && m[LEPTONIDPure_trailing]) {
                                                        c1Maps[cutFlow_Fine_onZMM].at(CUTFLOW_fine).Fill("2LeptonIDPure",mW[LEPTONIDPure_leading]*tempWeight);
                                                        if(m[LEPTONIDImpact_leading] && m[LEPTONIDImpact_trailing]) {
                                                                c1Maps[cutFlow_Fine_onZMM].at(CUTFLOW_fine).Fill("2LeptonIDImpact",mW[LEPTONIDImpact_leading]*tempWeight);
                                                                if(m[LEPTONIDEta_leading] && m[LEPTONIDEta_trailing]) {
                                                                        c1Maps[cutFlow_Fine_onZMM].at(CUTFLOW_fine).Fill("2LeptonIDEta",mW[LEPTONIDEta_leading]*tempWeight);
                                                                        if(m[LEPTONIDIso_leading] && m[LEPTONIDIso_trailing]) {
                                                                                c1Maps[cutFlow_Fine_onZMM].at(CUTFLOW_fine).Fill("2LeptonIDIso",mW[LEPTONIDIso_leading]*tempWeight);
                                                                                if(m[LEPTONIDDeltaR_leading] && m[LEPTONIDDeltaR_trailing]) {
                                                                                        c1Maps[cutFlow_Fine_onZMM].at(CUTFLOW_fine).Fill("2LeptonIDDeltaR",mW[LEPTONIDDeltaR_leading]*tempWeight);
                                                                                        if(m[LEPTONPT_leading] && m[LEPTONPT_trailing]) {
                                                                                                c1Maps[cutFlow_Fine_onZMM].at(CUTFLOW_fine).Fill("2LeptonPT",mW[LEPTONPT_leading]*tempWeight);
                                                                                                if(m[M50]) {
                                                                                                        c1Maps[cutFlow_Fine_onZMM].at(CUTFLOW_fine).Fill("m50",mW[M50]*tempWeight);
                                                                                                        if(m[PHOTON1]) {
                                                                                                                c1Maps[cutFlow_Fine_onZMM].at(CUTFLOW_fine).Fill("1Photon",mW[PHOTON1]*tempWeight);
                                                                                                                if(m[PHOTON1ID]) {
                                                                                                                        c1Maps[cutFlow_Fine_onZMM].at(CUTFLOW_fine).Fill("1PhotonID",mW[PHOTON1ID]*tempWeight);
                                                                                                                        if(m[PHOTON1SEED]) {
                                                                                                                                c1Maps[cutFlow_Fine_onZMM].at(CUTFLOW_fine).Fill("1PhotonSeed",mW[PHOTON1SEED]*tempWeight);
                                                                                                                                if(m[PHOTON1ETA]) {
                                                                                                                                        c1Maps[cutFlow_Fine_onZMM].at(CUTFLOW_fine).Fill("1PhotonEta",mW[PHOTON1ETA]*tempWeight);
                                                                                                                                        if(m[PHOTON1PT]) {
                                                                                                                                                c1Maps[cutFlow_Fine_onZMM].at(CUTFLOW_fine).Fill("1PhotonPT",mW[PHOTON1PT]*tempWeight);
                                                                                                                                                if(m[PHOTON1DR]) {
                                                                                                                                                        c1Maps[cutFlow_Fine_onZMM].at(CUTFLOW_fine).Fill("1PhotonDeltaR",mW[PHOTON1DR]*tempWeight);
                                                                                                                                                        if(m[ZMASS]) {
                                                                                                                                                                c1Maps[cutFlow_Fine_onZMM].at(CUTFLOW_fine).Fill("Z",mW[ZMASS]*tempWeight);
                                                                                                                                                                if(m[ex2LEP]) {
                                                                                                                                                                        c1Maps[cutFlow_Fine_onZMM].at(CUTFLOW_fine).Fill("ex2Lep",mW[ex2LEP]*tempWeight);
                                                                                                                                                                        if(m[MT2Cut]) {
                                                                                                                                                                                c1Maps[cutFlow_Fine_onZMM].at(CUTFLOW_fine).Fill("MT2",mW[MT2Cut]*tempWeight);
                                                                                                                                                                                if(m[MetCut]) {
                                                                                                                                                                                        c1Maps[cutFlow_Fine_onZMM].at(CUTFLOW_fine).Fill("Met",mW[MetCut]*tempWeight);
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
                        }
                }
                if(m[EMUON] && !(m[DIELECTRON]||m[DIMUON])) {
                        if (m[TRIGGERED]) {
                                c1Maps[cutFlow_Fine_onZEM].at(CUTFLOW_fine).Fill("triggered",mW[TRIGGERED]*tempWeight);
                                if(m[TRIGGEREDMATCHED]) {
                                        c1Maps[cutFlow_Fine_onZEM].at(CUTFLOW_fine).Fill("triggeredMatched",mW[TRIGGEREDMATCHED]*tempWeight);
                                        if(m[GENVETO]) {
                                                c1Maps[cutFlow_Fine_onZEM].at(CUTFLOW_fine).Fill("GenPhotonVeto",mW[GENVETO]*tempWeight);
                                                if(m[LEPTONIDPure_leading] && m[LEPTONIDPure_trailing]) {
                                                        c1Maps[cutFlow_Fine_onZEM].at(CUTFLOW_fine).Fill("2LeptonIDPure",mW[LEPTONIDPure_leading]*tempWeight);
                                                        if(m[LEPTONIDImpact_leading] && m[LEPTONIDImpact_trailing]) {
                                                                c1Maps[cutFlow_Fine_onZEM].at(CUTFLOW_fine).Fill("2LeptonIDImpact",mW[LEPTONIDImpact_leading]*tempWeight);
                                                                if(m[LEPTONIDEta_leading] && m[LEPTONIDEta_trailing]) {
                                                                        c1Maps[cutFlow_Fine_onZEM].at(CUTFLOW_fine).Fill("2LeptonIDEta",mW[LEPTONIDEta_leading]*tempWeight);
                                                                        if(m[LEPTONIDIso_leading] && m[LEPTONIDIso_trailing]) {
                                                                                c1Maps[cutFlow_Fine_onZEM].at(CUTFLOW_fine).Fill("2LeptonIDIso",mW[LEPTONIDIso_leading]*tempWeight);
                                                                                if(m[LEPTONIDDeltaR_leading] && m[LEPTONIDDeltaR_trailing]) {
                                                                                        c1Maps[cutFlow_Fine_onZEM].at(CUTFLOW_fine).Fill("2LeptonIDDeltaR",mW[LEPTONIDDeltaR_leading]*tempWeight);
                                                                                        if(m[LEPTONPT_leading] && m[LEPTONPT_trailing]) {
                                                                                                c1Maps[cutFlow_Fine_onZEM].at(CUTFLOW_fine).Fill("2LeptonPT",mW[LEPTONPT_leading]*tempWeight);
                                                                                                if(m[M50]) {
                                                                                                        c1Maps[cutFlow_Fine_onZEM].at(CUTFLOW_fine).Fill("m50",mW[M50]*tempWeight);
                                                                                                        if(m[PHOTON1]) {
                                                                                                                c1Maps[cutFlow_Fine_onZEM].at(CUTFLOW_fine).Fill("1Photon",mW[PHOTON1]*tempWeight);
                                                                                                                if(m[PHOTON1ID]) {
                                                                                                                        c1Maps[cutFlow_Fine_onZEM].at(CUTFLOW_fine).Fill("1PhotonID",mW[PHOTON1ID]*tempWeight);
                                                                                                                        if(m[PHOTON1SEED]) {
                                                                                                                                c1Maps[cutFlow_Fine_onZEM].at(CUTFLOW_fine).Fill("1PhotonSeed",mW[PHOTON1SEED]*tempWeight);
                                                                                                                                if(m[PHOTON1ETA]) {
                                                                                                                                        c1Maps[cutFlow_Fine_onZEM].at(CUTFLOW_fine).Fill("1PhotonEta",mW[PHOTON1ETA]*tempWeight);
                                                                                                                                        if(m[PHOTON1PT]) {
                                                                                                                                                c1Maps[cutFlow_Fine_onZEM].at(CUTFLOW_fine).Fill("1PhotonPT",mW[PHOTON1PT]*tempWeight);
                                                                                                                                                if(m[PHOTON1DR]) {
                                                                                                                                                        c1Maps[cutFlow_Fine_onZEM].at(CUTFLOW_fine).Fill("1PhotonDeltaR",mW[PHOTON1DR]*tempWeight);
                                                                                                                                                        if(m[ZMASS]) {
                                                                                                                                                                c1Maps[cutFlow_Fine_onZEM].at(CUTFLOW_fine).Fill("Z",mW[ZMASS]*tempWeight);
                                                                                                                                                                if(m[ex2LEP]) {
                                                                                                                                                                        c1Maps[cutFlow_Fine_onZEM].at(CUTFLOW_fine).Fill("ex2Lep",mW[ex2LEP]*tempWeight);
                                                                                                                                                                        if(m[MT2Cut]) {
                                                                                                                                                                                c1Maps[cutFlow_Fine_onZEM].at(CUTFLOW_fine).Fill("MT2",mW[MT2Cut]*tempWeight);
                                                                                                                                                                                if(m[MetCut]) {
                                                                                                                                                                                        c1Maps[cutFlow_Fine_onZEM].at(CUTFLOW_fine).Fill("Met",mW[MetCut]*tempWeight);
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
                        }
                }
        }
}

template<typename T>
//void save2File(const map<string,T>& hMaps, TFile& file)
void save2File(map<string,T>& hMaps, TFile& file)
{
        for (auto& hMapIt : hMaps) {

                //if (!file.Get(hMapIt.first.c_str())) {
                //file.mkdir(hMapIt.first.c_str());
                //}
                //file.cd(hMapIt.first.c_str());
                //for (auto& h : hMapIt.second) {
                string saveName = hMapIt.first;
                //saveName=saveName.replace(0,11,"");
                //cout<<"saving "<<saveName<<endl;
                //hMapIt.second.Write(hMapIt.first.replace(0,10,"").c_str(), TObject::kWriteDelete);
                //cout<<"a"<<endl;
                hMapIt.second.Write(saveName.c_str(), TObject::kWriteDelete);
                //cout<<"b"<<endl;
                hMapIt.second.SetDirectory(0);
                //cout<<"c"<<endl;
                //}
                //file.cd();
        }
}

//template<typename T>
//void save2File(const map<string,map<Histograms1D,T>>& hMaps, TFile& file)
//{
//for (auto& hMapIt : hMaps) {
//if (!file.Get(hMapIt.first.c_str())) {
//file.mkdir(hMapIt.first.c_str());
//}
//file.cd(hMapIt.first.c_str());
//for (auto& h : hMapIt.second) {
//h.second.Write(histoNames[h.first].c_str(), TObject::kWriteDelete);
//}
//file.cd();
//}
//}
template<typename T>
void save2File(const map<cutFlowMapName,map<Histograms1D,T> >& hMaps, TFile& file)
{
        for (auto& hMapIt : hMaps) {
                if (!file.Get(cutFlowMapNameString[hMapIt.first].c_str())) {
                        file.mkdir(cutFlowMapNameString[hMapIt.first].c_str());
                }
                file.cd(cutFlowMapNameString[hMapIt.first].c_str());
                for (auto& h : hMapIt.second) {
                        h.second.Write(histoNames[h.first].c_str(), TObject::kWriteDelete);
                }
                file.cd();
        }
}
template<typename T>
void save2File(const map<selectionFolderName,map<Histograms1D,T> >& hMaps, TFile& file)
{
        for (auto& hMapIt : hMaps) {
                if (!file.Get(selectionFolderNameString[hMapIt.first].c_str())) {
                        file.mkdir(selectionFolderNameString[hMapIt.first].c_str());
                }
                file.cd(selectionFolderNameString[hMapIt.first].c_str());
                for (auto& h : hMapIt.second) {
                        h.second.Write(histoNames[h.first].c_str(), TObject::kWriteDelete);
                }
                file.cd();
        }
}
template<typename T>
void save2File(const map<selectionFolderName,map<Histograms2D,T> >& hMaps, TFile& file)
{
        for (auto& hMapIt : hMaps) {
                if (!file.Get(selectionFolderNameString[hMapIt.first].c_str())) {
                        file.mkdir(selectionFolderNameString[hMapIt.first].c_str());
                }
                file.cd(selectionFolderNameString[hMapIt.first].c_str());
                for (auto& h : hMapIt.second) {
                        h.second.Write(histoNames2D[h.first].c_str(), TObject::kWriteDelete);
                }
                file.cd();
        }
}
template<typename T>
void save2File(const map<string,map<selectionFolderName,map<Histograms1D,T> > >& hMaps, TFile& file)
{
        for (auto& hMapIt : hMaps) {
                if (!file.Get(hMapIt.first.c_str())) {
                        file.mkdir(hMapIt.first.c_str());
                }
                file.cd(hMapIt.first.c_str());
                for(auto& hMapIt2 : hMapIt.second) {
                        if (!file.Get(selectionFolderNameString[hMapIt2.first].c_str())) {
                                file.mkdir((hMapIt.first+"/"+selectionFolderNameString[hMapIt2.first]).c_str());
                        }
                        file.cd((selectionFolderNameString[hMapIt.first]+"/"+selectionFolderNameString[hMapIt2.first]).c_str());
                        for (auto h: hMapIt2.second) {
                                h.second.Write(histoNames[h.first].c_str(), TObject::kWriteDelete);
                        }
                        file.cd();
                }
        }
}
template<typename T>
void save2File(const map<selectionFolderName,map<selectionFolderName,map<Histograms1D,T> > >& hMaps, TFile& file)
{
        for (auto& hMapIt : hMaps) {
                if (!file.Get(selectionFolderNameString[hMapIt.first].c_str())) {
                        file.mkdir(selectionFolderNameString[hMapIt.first].c_str());
                }
                file.cd(selectionFolderNameString[hMapIt.first].c_str());
                for(auto& hMapIt2 : hMapIt.second) {
                        if (!file.Get(selectionFolderNameString[hMapIt2.first].c_str())) {
                                file.mkdir((selectionFolderNameString[hMapIt.first]+"/"+selectionFolderNameString[hMapIt2.first]).c_str());
                        }
                        file.cd((selectionFolderNameString[hMapIt.first]+"/"+selectionFolderNameString[hMapIt2.first]).c_str());
                        for (auto h: hMapIt2.second) {
                                h.second.Write(histoNames[h.first].c_str(), TObject::kWriteDelete);
                        }
                        file.cd();
                }
        }
}
template<typename T>
void save2File(const map<selectionFolderName,map<selectionFolderName,map<Histograms2D,T> > >& hMaps, TFile& file)
{
        for (auto& hMapIt : hMaps) {
                if (!file.Get(selectionFolderNameString[hMapIt.first].c_str())) {
                        file.mkdir(selectionFolderNameString[hMapIt.first].c_str());
                }
                file.cd(selectionFolderNameString[hMapIt.first].c_str());
                for(auto& hMapIt2 : hMapIt.second) {
                        if (!file.Get(selectionFolderNameString[hMapIt2.first].c_str())) {
                                file.mkdir((selectionFolderNameString[hMapIt.first]+"/"+selectionFolderNameString[hMapIt2.first]).c_str());
                        }
                        file.cd((selectionFolderNameString[hMapIt.first]+"/"+selectionFolderNameString[hMapIt2.first]).c_str());
                        for (auto h: hMapIt2.second) {
                                h.second.Write(histoNames2D[h.first].c_str(), TObject::kWriteDelete);
                        }
                        file.cd();
                }
        }
}
//template<typename T>
void save2File(const map<string,map<selectionFolderName,map<selectionFolderName,map<Histograms1D,TH1F> > > >& hMaps, TFile& file)
{
        for (auto& hMapIt : hMaps) {
                if (!file.Get(hMapIt.first.c_str())) {
                        file.mkdir(hMapIt.first.c_str());
                }
                file.cd(hMapIt.first.c_str());
                for(auto& hMapIt2 : hMapIt.second) {
                        if (!file.Get(selectionFolderNameString[hMapIt2.first].c_str())) {
                                file.mkdir((hMapIt.first+"/"+selectionFolderNameString[hMapIt2.first]).c_str());
                        }
                        file.cd((hMapIt.first+"/"+selectionFolderNameString[hMapIt2.first]).c_str());
                        for(auto& hMapIt3 : hMapIt2.second) {
                                if (!file.Get(selectionFolderNameString[hMapIt3.first].c_str())) {
                                        file.mkdir((hMapIt.first+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]).c_str());
                                }
                                file.cd((hMapIt.first+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]).c_str());
                                for (auto h: hMapIt3.second) {
                                        h.second.Write(histoNames[h.first].c_str(), TObject::kWriteDelete);
                                }
                                file.cd();
                        }
                }
        }
}
void save2File(const map<string,map<selectionFolderName,map<selectionFolderName,map<selectionFolderName,map<Histograms1D,TH1F> > > > >& hMaps, TFile& file)
{
        for (auto& hMapIt : hMaps) {
                if (!file.Get(hMapIt.first.c_str())) {
                        file.mkdir(hMapIt.first.c_str());
                }
                file.cd(hMapIt.first.c_str());
                for(auto& hMapIt2 : hMapIt.second) {
                        if (!file.Get(selectionFolderNameString[hMapIt2.first].c_str())) {
                                file.mkdir((hMapIt.first+"/"+selectionFolderNameString[hMapIt2.first]).c_str());
                        }
                        file.cd((hMapIt.first+"/"+selectionFolderNameString[hMapIt2.first]).c_str());
                        for(auto& hMapIt3 : hMapIt2.second) {
                                if (!file.Get(selectionFolderNameString[hMapIt3.first].c_str())) {
                                        file.mkdir((hMapIt.first+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]).c_str());
                                }
                                file.cd((hMapIt.first+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]).c_str());
                                for(auto& hMapIt4 : hMapIt3.second) {
                                        if (!file.Get(selectionFolderNameString[hMapIt4.first].c_str())) {
                                                file.mkdir((hMapIt.first+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]+"/"+selectionFolderNameString[hMapIt4.first]).c_str());
                                        }
                                        file.cd((hMapIt.first+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]+"/"+selectionFolderNameString[hMapIt4.first]).c_str());
                                        for (auto h: hMapIt4.second) {
                                                h.second.Write(histoNames[h.first].c_str(), TObject::kWriteDelete);
                                        }
                                        file.cd();
                                }
                        }
                }
        }
}
void save2File(const map<SignalPoint,map<selectionFolderName,map<selectionFolderName,map<selectionFolderName,map<Histograms1D,TH1F> > > > >& hMaps, TFile& file,int nbinos, string nameOfFile)
{
        for (auto& hMapIt : hMaps) {
                //if (!file.Get(hMapIt.first.c_str())) {
                //file.mkdir(hMapIt.first.c_str());
                //}
                string inputStringName = nameOfFile;
                string dirSignalPoint;
                if(inputStringName.find("GGM")!=string::npos) {
                        dirSignalPoint = getSignalPointNameShort(10);
                }else{
                        if(inputStringName.find("GMSB")!=string::npos) {
                                dirSignalPoint = getSignalPointNameShort(11);
                        }else{
                                dirSignalPoint = getSignalPointNameShort(nbinos);
                        }
                }
                auto dir = (dirSignalPoint+"_"+to_string(hMapIt.first.first)+"_"+to_string(hMapIt.first.second));
                if (!file.Get(dir.c_str())) {
                        file.mkdir(dir.c_str());
                }
                //file.cd(hMapIt.first.c_str());
                file.cd(dir.c_str());
                for(auto& hMapIt2 : hMapIt.second) {
                        //cout<<selectionFolderNameString[hMapIt2.first]<<endl;
                        //if(selectionFolderNameString[hMapIt2.first]=="VR"){
                        //file.ls();
                        //}
                        //if (!file.Get(selectionFolderNameString[hMapIt2.first].c_str())) {
                        if (!file.Get((dir+"/"+selectionFolderNameString[hMapIt2.first]).c_str())) {
                                //cout<<"here"<<endl;
                                //file.mkdir((hMapIt.first+"/"+selectionFolderNameString[hMapIt2.first]).c_str());
                                file.mkdir((dir+"/"+selectionFolderNameString[hMapIt2.first]).c_str());
                                //cout<<"created"<<endl;
                        }
                        //file.cd((hMapIt.first+"/"+selectionFolderNameString[hMapIt2.first]).c_str());
                        file.cd((dir+"/"+selectionFolderNameString[hMapIt2.first]).c_str());
                        for(auto& hMapIt3 : hMapIt2.second) {
                                if (!file.Get(selectionFolderNameString[hMapIt3.first].c_str())) {
                                        //file.mkdir((hMapIt.first+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]).c_str());
                                        file.mkdir((dir+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]).c_str());
                                }
                                //file.cd((hMapIt.first+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]).c_str());
                                file.cd((dir+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]).c_str());
                                for(auto& hMapIt4 : hMapIt3.second) {
                                        if (!file.Get(selectionFolderNameString[hMapIt4.first].c_str())) {
                                                //file.mkdir((hMapIt.first+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]+"/"+selectionFolderNameString[hMapIt4.first]).c_str());
                                                file.mkdir((dir+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]+"/"+selectionFolderNameString[hMapIt4.first]).c_str());
                                        }
                                        //file.cd((hMapIt.first+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]+"/"+selectionFolderNameString[hMapIt4.first]).c_str());
                                        file.cd((dir+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]+"/"+selectionFolderNameString[hMapIt4.first]).c_str());
                                        for (auto h: hMapIt4.second) {
                                                h.second.Write(histoNames[h.first].c_str(), TObject::kWriteDelete);
                                        }
                                        file.cd();
                                }
                        }
                }
        }
}
void save2File(const map<SignalPoint,map<selectionFolderName,map<selectionFolderName,map<selectionFolderName,map<Histograms2D,TH2F> > > > >& hMaps, TFile& file,int nbinos, string nameOfFile)
{
        for (auto& hMapIt : hMaps) {
                //if (!file.Get(hMapIt.first.c_str())) {
                //file.mkdir(hMapIt.first.c_str());
                //}
                string inputStringName = nameOfFile;
                string dirSignalPoint;
                if(inputStringName.find("GGM")!=string::npos) {
                        dirSignalPoint = getSignalPointNameShort(10);
                }else{
                        if(inputStringName.find("GMSB")!=string::npos) {
                                dirSignalPoint = getSignalPointNameShort(11);
                        }else{
                                dirSignalPoint = getSignalPointNameShort(nbinos);
                        }
                }
                auto dir = (dirSignalPoint+"_"+to_string(hMapIt.first.first)+"_"+to_string(hMapIt.first.second));
                if (!file.Get(dir.c_str())) {
                        file.mkdir(dir.c_str());
                }
                //file.cd(hMapIt.first.c_str());
                file.cd(dir.c_str());
                for(auto& hMapIt2 : hMapIt.second) {
                        //cout<<selectionFolderNameString[hMapIt2.first]<<endl;
                        //if(selectionFolderNameString[hMapIt2.first]=="VR"){
                        //file.ls();
                        //}
                        //if (!file.Get(selectionFolderNameString[hMapIt2.first].c_str())) {
                        if (!file.Get((dir+"/"+selectionFolderNameString[hMapIt2.first]).c_str())) {
                                //cout<<"here"<<endl;
                                //file.mkdir((hMapIt.first+"/"+selectionFolderNameString[hMapIt2.first]).c_str());
                                file.mkdir((dir+"/"+selectionFolderNameString[hMapIt2.first]).c_str());
                                //cout<<"created"<<endl;
                        }
                        //file.cd((hMapIt.first+"/"+selectionFolderNameString[hMapIt2.first]).c_str());
                        file.cd((dir+"/"+selectionFolderNameString[hMapIt2.first]).c_str());
                        for(auto& hMapIt3 : hMapIt2.second) {
                                if (!file.Get(selectionFolderNameString[hMapIt3.first].c_str())) {
                                        //file.mkdir((hMapIt.first+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]).c_str());
                                        file.mkdir((dir+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]).c_str());
                                }
                                //file.cd((hMapIt.first+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]).c_str());
                                file.cd((dir+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]).c_str());
                                for(auto& hMapIt4 : hMapIt3.second) {
                                        if (!file.Get(selectionFolderNameString[hMapIt4.first].c_str())) {
                                                //file.mkdir((hMapIt.first+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]+"/"+selectionFolderNameString[hMapIt4.first]).c_str());
                                                file.mkdir((dir+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]+"/"+selectionFolderNameString[hMapIt4.first]).c_str());
                                        }
                                        //file.cd((hMapIt.first+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]+"/"+selectionFolderNameString[hMapIt4.first]).c_str());
                                        file.cd((dir+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]+"/"+selectionFolderNameString[hMapIt4.first]).c_str());
                                        for (auto h: hMapIt4.second) {
                                                h.second.Write(histoNames2D[h.first].c_str(), TObject::kWriteDelete);
                                        }
                                        file.cd();
                                }
                        }
                }
        }
}
//template<typename T>
void save2File(const map<selectionFolderName,map<selectionFolderName,map<selectionFolderName,map<Histograms1D,TH1F> > > >& hMaps, TFile& file)
{
        for (auto& hMapIt : hMaps) {
                if (!file.Get(selectionFolderNameString[hMapIt.first].c_str())) {
                        file.mkdir(selectionFolderNameString[hMapIt.first].c_str());
                }
                file.cd(selectionFolderNameString[hMapIt.first].c_str());
                for(auto& hMapIt2 : hMapIt.second) {
                        if (!file.Get(selectionFolderNameString[hMapIt2.first].c_str())) {
                                file.mkdir((selectionFolderNameString[hMapIt.first]+"/"+selectionFolderNameString[hMapIt2.first]).c_str());
                        }
                        file.cd((selectionFolderNameString[hMapIt.first]+"/"+selectionFolderNameString[hMapIt2.first]).c_str());
                        for(auto& hMapIt3 : hMapIt2.second) {
                                if (!file.Get(selectionFolderNameString[hMapIt3.first].c_str())) {
                                        file.mkdir((selectionFolderNameString[hMapIt.first]+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]).c_str());
                                }
                                file.cd((selectionFolderNameString[hMapIt.first]+"/"+selectionFolderNameString[hMapIt2.first]+"/"+selectionFolderNameString[hMapIt3.first]).c_str());
                                for (auto h: hMapIt3.second) {
                                        h.second.Write(histoNames[h.first].c_str(), TObject::kWriteDelete);
                                }
                                file.cd();
                        }
                }
        }
}

//template<typename T>
void save2File(const map<selectionFolderName,map<selectionFolderName,map<string,map<Histograms1D,TH1F> > > >& hMaps, TFile& file)
{
        for (auto& hMapIt : hMaps) {
                if (!file.Get(selectionFolderNameString[hMapIt.first].c_str())) {
                        file.mkdir(selectionFolderNameString[hMapIt.first].c_str());
                }
                file.cd(selectionFolderNameString[hMapIt.first].c_str());
                for(auto& hMapIt2 : hMapIt.second) {
                        if (!file.Get(selectionFolderNameString[hMapIt2.first].c_str())) {
                                file.mkdir((selectionFolderNameString[hMapIt.first]+"/"+selectionFolderNameString[hMapIt2.first]).c_str());
                        }
                        file.cd((selectionFolderNameString[hMapIt.first]+"/"+selectionFolderNameString[hMapIt2.first]).c_str());
                        for(auto& hMapIt3 : hMapIt2.second) {
                                if (!file.Get(hMapIt3.first.c_str())) {
                                        file.mkdir((selectionFolderNameString[hMapIt.first]+"/"+selectionFolderNameString[hMapIt2.first]+"/"+hMapIt3.first).c_str());
                                }
                                file.cd((selectionFolderNameString[hMapIt.first]+"/"+selectionFolderNameString[hMapIt2.first]+"/"+hMapIt3.first).c_str());
                                for (auto h: hMapIt3.second) {
                                        h.second.Write(histoNames[h.first].c_str(), TObject::kWriteDelete);
                                }
                                file.cd();
                        }
                }
        }
}






void myAnalyzer::Terminate()
{
        auto outputName = "output"+config_outputfolder+"/"+getOutputFilename(inputName);
        //if(!config_doTrigger){
        if(!(config_doTrigger||config_doTriggerMET||config_doTriggerPure)) {
                outputName = "output"+config_outputfolder+"/"+getOutputFilename(inputName);
        }else{
                if(config_doTrigger) {
                        outputName = "output"+config_outputfolder+"/ht/"+getOutputFilename(inputName);
                }
                if(config_doTriggerPure) {
                        outputName = "output"+config_outputfolder+"/htPure/"+getOutputFilename(inputName);
                }
                if(config_doTriggerMET) {
                        outputName = "output"+config_outputfolder+"/met/"+getOutputFilename(inputName);
                }
        }

        if(config_doPdfSplit) {
                outputName = outputName + "_" +to_string(config_pdfStart)+"-"+to_string(config_pdfStep-1);
        }

        TFile file(outputName.c_str(), "RECREATE");
        //cout<<"1"<<endl;
        save2File(h1Maps, file);
        //cout<<"2"<<endl;
        save2File(eff1Maps, file);
        //cout<<"3"<<endl;
        save2File(h2Maps, file);
        //cout<<"4"<<endl;
        save2File(c1Maps, file);
        //save2File(s1Maps, file);
        //cout<<"5"<<endl;
        save2File(s1Maps, file,*nBinos,inputName);
        save2File(s2Maps, file,*nBinos,inputName);
        //cout<<"6"<<endl;
        save2File(cr1Maps, file);
        //cout<<"7"<<endl;
        if(!isTotalSignal) {
                cutFlow.Write("hCutFlow");
                //if(!config_doTrigger)weightHisto.Write("weightHisto");
                weightHisto.Write("weightHisto");
        }else{
                //if(!config_doTrigger)save2File(weightHistoMap,file);
                if(!(config_doTrigger||config_doTriggerMET||config_doTriggerPure)) save2File(weightHistoMap,file);
        }
        file.Close();
        cout << "Created " << outputName << " in " << (time(NULL) - startTime)/60 << " min" << endl;
}
