#include "treeSplitter.h"

treeSplitter::treeSplitter() :
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
        met_gen(fReader, "met_gen"),
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

        hlt_met110(fReader, "HLT_PFMET110_PFMHT110_IDTight_v"),
        // HLT_PFMET110_PFMHT110_IDTight_v
        // HLT_PFMET110_PFMHT110_IDTight_v
        hlt_met120(fReader, "HLT_PFMET120_PFMHT120_IDTight_v"),
        hlt_met170_Noise(fReader, "HLT_PFMET170_NoiseCleaned_v"),
        hlt_met170_HBHE(fReader, "HLT_PFMET170_HBHECleaned_v"),
        hlt_met170_Jet(fReader, "HLT_PFMET170_JetIdCleaned_v"),
        hlt_met170_Not(fReader, "HLT_PFMET170_NotCleaned_v"),
        hlt_met300(fReader, "HLT_PFMET300_v"),
        hlt_met400(fReader, "HLT_PFMET400_v"),
        hlt_met500(fReader, "HLT_PFMET500_v"),
        hlt_met600(fReader, "HLT_PFMET600_v"),



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


void treeSplitter::Init(TTree *tree)
{
        fReader.SetTree(tree);
        inputName = fReader.GetTree()->GetCurrentFile()->GetName();

        isTotalSignal = (inputName.find("SMS-TChiNG_BF") != string::npos) || (inputName.find("GMSB_GravitinoLSP") != string::npos) || (inputName.find("SMS-T6ttZg_nTuple")!=string::npos) || (inputName.find("SMS-T5bbbbZg_nTuple")!=string::npos) || (inputName.find("GGM")!=string::npos);
        isData = inputName.find("Run201") != string::npos;
        isSignal = (inputName.find("SMS") != string::npos) || (inputName.find("GGM") != string::npos) || (inputName.find("GMSB") != string::npos);

        boost::property_tree::ini_parser::read_ini("example.ini", propertyTree);

        //doWeights_TopPt = (inputName.find("TTTo") != string::npos) || (inputName.find("TTGamma") != string::npos) || (inputName.find("ST_") != string::npos);
        //doWeights_TopPt = (inputName.find("TTTo") != string::npos) || (inputName.find("TTGamma") != string::npos);
        doWeights_TopPt = (inputName.find("TTTo") != string::npos);
        //doWeights_TopPt = false;
        doWeights_EWKinoPairPt = (inputName.find("TChi") != string::npos);
        doWeights_LeptonPairPt = false;
        //doWeights_nISR = (inputName.find("T5") != string::npos) || (inputName.find("TTGamma") != string::npos);
        doWeights_nISR = (inputName.find("T5") != string::npos);

        config_veto = propertyTree.get<int>("generalSettings.veto");
        config_eventpercentage = propertyTree.get<float>("generalSettings.eventpercentage");
        config_outputfolder = propertyTree.get<string>("generalSettings.outputfolder");
        config_doHt = propertyTree.get<bool>("generalSettings.ht");
        config_doHtPure = propertyTree.get<bool>("generalSettings.htpure");
        config_doMET = propertyTree.get<bool>("generalSettings.met");
        //cout<<config_doHt<<endl;


        if(!isTotalSignal) {
                cutFlow = *((TH1F*)fReader.GetTree()->GetCurrentFile()->Get("TreeWriter/hCutFlow"));

                weightHisto = TH1F("weightHisto","weightHisto",120,0,120);
                weightHisto.Fill("none",0);
                weightHisto.Fill("pu", 0);
                weightHisto.Fill("pu_mc", 0);
                weightHisto.Fill("pu_mc_nisr", 0);
                weightHisto.Fill("pu_mc_nisrUp", 0);
                weightHisto.Fill("pu_mc_nisrDn", 0);
                weightHisto.Fill("pu_mc_toppt", 0);
                weightHisto.Fill("pu_mc_ewk", 0);
                weightHisto.Fill("pu_mc_ewkUp", 0);
                weightHisto.Fill("pu_mc_ewkDn", 0);
                for(int i=0; i<(110); i++) {
                        string tempPDFstring="pu_mc_pdf_";
                        tempPDFstring=tempPDFstring+to_string(i);
                        weightHisto.Fill(tempPDFstring.c_str(), 0);
                }
        }else{
                //TList* list = fReader.GetTree()->GetCurrentFile()->GetListOfKeys() ;
                //TList* list = fReader.GetTree()->GetCurrentFile()->cd("TreeWriter")->GetListOfKeys() ;
                TList* list = fReader.GetTree()->GetCurrentFile()->GetDirectory("TreeWriter")->GetListOfKeys();
                //auto tempFile=fReader.GetTree()->GetCurrentFile();
                //tempFile->cd("TreeWriter");
                //TList* list = GetListOfKeys() ;
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

                                //cout<<text<<endl;

                                boost::split(results, text, boost::is_any_of("_"), boost::token_compress_on);

                                //unsigned short toFillM1=0;
                                //unsigned short toFillM2=0;

                                //string tempWeightHistName="weightHisto"
                                string tempWeightHistName="";

                                //cout<<results.size()<<endl;

                                if(results.size()>2) {
                                        //toFillM1=atoi(results.at(1).c_str());
                                        //toFillM2=atoi(results.at(2).c_str());
                                        if((inputName.find("GGM") != string::npos)) {
                                                tempWeightHistName=tempWeightHistName+"_"+results.at(1).replace(0,2,"").c_str()+"_"+results.at(2).replace(0,2,"").c_str();
                                        }else{
                                                tempWeightHistName=tempWeightHistName+"_"+results.at(1).c_str()+"_"+results.at(2).c_str();
                                        }
                                        //cout<<tempWeightHistName<<endl;
                                }else{
                                        if(results.size()>1) {
                                                //toFillM1=atoi(results.at(1).c_str());
                                                tempWeightHistName=tempWeightHistName+"_"+results.at(1).c_str();
                                        }
                                }
                                //cout<<tempWeightHistName<<endl;
                                weightHistoMap[tempWeightHistName] = TH1F(("weightHisto"+tempWeightHistName).c_str(),("weightHisto"+tempWeightHistName).c_str(),120,0,120);
                                weightHistoMap[tempWeightHistName].Fill("none",0);
                                weightHistoMap[tempWeightHistName].Fill("pu", 0);
                                weightHistoMap[tempWeightHistName].Fill("pu_mc", 0);
                                weightHistoMap[tempWeightHistName].Fill("pu_mc_nisr", 0);
                                weightHistoMap[tempWeightHistName].Fill("pu_mc_nisrUp", 0);
                                weightHistoMap[tempWeightHistName].Fill("pu_mc_nisrDn", 0);
                                weightHistoMap[tempWeightHistName].Fill("pu_mc_toppt", 0);
                                weightHistoMap[tempWeightHistName].Fill("pu_mc_ewk", 0);
                                weightHistoMap[tempWeightHistName].Fill("pu_mc_ewkUp", 0);
                                weightHistoMap[tempWeightHistName].Fill("pu_mc_ewkDn", 0);
                                for(int i=0; i<(110); i++) {
                                        //weightHistoMap[tempWeightHistName].Fill("pu_pdf"+to_string(i), 0);
                                        string tempPDFstring="pu_mc_pdf_";
                                        tempPDFstring=tempPDFstring+to_string(i);
                                        weightHistoMap[tempWeightHistName].Fill(tempPDFstring.c_str(), 0);
                                }
                                //weightHisto = TH1F("weightHisto","weightHisto",120,0,120);
                                //weightHisto.Fill("none",0);
                                //weightHisto.Fill("pu", 0);
                                //weightHisto.Fill("pu_mc", 0);
                                //weightHisto.Fill("pu_mc_nisr", 0);
                                //weightHisto.Fill("pu_mc_nisrUp", 0);
                                //weightHisto.Fill("pu_mc_nisrDn", 0);
                                //weightHisto.Fill("pu_mc_toppt", 0);
                                //weightHisto.Fill("pu_mc_ewk", 0);
                                //weightHisto.Fill("pu_mc_ewkUp", 0);
                                //weightHisto.Fill("pu_mc_ewkDn", 0);
                                //for(int i=0; i<(110); i++){
                                //string tempPDFstring="pu_mc_pdf_";
                                //tempPDFstring=tempPDFstring+to_string(i);
                                //weightHisto.Fill(tempPDFstring.c_str(), 0);
                                //}

                                //weightHistoMap[]
                                //sp_.first=toFillM1;
                                //sp_.second=toFillM2;
                                //massPointsForMaps.push_back(sp_);
                        }
                }
        }
        fReader.GetEntries(true);

        nEntries=fReader.GetTree()->GetEntries();


        string pileUpHistotoUseUp;
        string pileUpHistotoUseDown;
        if((inputName.find("GMSB") != string::npos)) {
                pileUpHistotoUseUp = "pileupWeightUp_mix_2016_25ns_Moriond17MC_PoissonOOTPU";
                pileUpHistotoUseDown = "pileupWeightDown_mix_2016_25ns_Moriond17MC_PoissonOOTPU";
        }else{
                if((inputName.find("GGM") != string::npos)) {
                        pileUpHistotoUseUp = "pileupWeightUp_mix_2016_25ns_Moriond17MC_PoissonOOTPU";
                        pileUpHistotoUseDown = "pileupWeightDown_mix_2016_25ns_Moriond17MC_PoissonOOTPU";
                }else{
                        if((inputName.find("SMS") != string::npos)) {
                                pileUpHistotoUseUp = "pileupWeightUp_mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU";
                                pileUpHistotoUseDown = "pileupWeightDown_mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU";
                        }else{
                                pileUpHistotoUseUp = "pileupWeightUp_mix_2016_25ns_Moriond17MC_PoissonOOTPU";
                                pileUpHistotoUseDown = "pileupWeightDown_mix_2016_25ns_Moriond17MC_PoissonOOTPU";
                        }
                }
        }

        puWeighterUp= Weighter("/home/home4/institut_1b/swuchterl/cmssw/TreeWriter_16/CMSSW_8_0_26_patch2/src/TreeWriter/PUreweighting/data/puWeights.root", pileUpHistotoUseUp);
        puWeighterDown= Weighter("/home/home4/institut_1b/swuchterl/cmssw/TreeWriter_16/CMSSW_8_0_26_patch2/src/TreeWriter/PUreweighting/data/puWeights.root", pileUpHistotoUseDown);




        if(isSignal) {
                nWeights=9;
        }else{
                nWeights=isData ? 0 : 110;
        }

        noPromptPhotons = inputName.find("DYJets") != string::npos;
        isZGammaInclusive = (inputName.find("ZGTo2LG_ext") != string::npos) || (inputName.find("ZGTo2LG_Total") != string::npos);

        //InitScaleFactorsFinal();
        InitScaleFactorsFinalMedium();




        if(config_doHt || config_doMET||config_doHtPure) {
                InitHTTree();
        }else{
                InitTree();
        }
}






void treeSplitter::InitTree(){


        outputFilenameCompressed = getOutputFilename(inputName,"myTuple");

        //cout<<config_outputfolder<<endl;

        auto saveStr= config_outputfolder+"/"+outputFilenameCompressed;

        //TFile newfile(saveStr.c_str(),"RECREATE");

        //TFile* treeFile = new TFile(saveStr.c_str(),"RECREATE");
        treeFile = new TFile(saveStr.c_str(),"RECREATE");

        eeTree = new TTree("Tree","Data");


        eeTree->Branch("selPhotons",&selectedEvent.selPhotons);
        eeTree->Branch("selJets",&selectedEvent.selJets);
        eeTree->Branch("selElectrons",&selectedEvent.selElectrons);
        eeTree->Branch("selMuons",&selectedEvent.selMuons);

        eeTree->Branch("chargeProduct",&selectedEvent.chargeProduct);
        eeTree->Branch("mll",&selectedEvent.mll);
        eeTree->Branch("miniIso1",&selectedEvent.miniIso1);
        eeTree->Branch("miniIso2",&selectedEvent.miniIso2);
        eeTree->Branch("deltaRll",&selectedEvent.deltaRll);

        eeTree->Branch("pt1",&selectedEvent.pt1);
        eeTree->Branch("pt2",&selectedEvent.pt2);
        eeTree->Branch("phi1",&selectedEvent.phi1);
        eeTree->Branch("phi2",&selectedEvent.phi2);
        eeTree->Branch("eta1",&selectedEvent.eta1);
        eeTree->Branch("eta2",&selectedEvent.eta2);
        eeTree->Branch("l1",&selectedEvent.l1);
        eeTree->Branch("l2",&selectedEvent.l2);

        eeTree->Branch("negMuons",&temp_negMuons);
        eeTree->Branch("posMuons",&temp_posMuons);
        eeTree->Branch("negElectrons",&temp_negElectrons);
        eeTree->Branch("posElectrons",&temp_posElectrons);

        eeTree->Branch("countNegCharge",&temp_countNegCharge);
        eeTree->Branch("countPosCharge",&temp_countPosCharge);

        eeTree->Branch("MT2",&selectedEvent.MT2_val);

        eeTree->Branch("genPhotonVeto",&selectedEvent.evtHasGenPhotonVeto);

        //eeTree->Branch("genZToLLCounter",&genZToLLCounter);

        eeTree->Branch("ETmiss",&selectedEvent.ETmiss);
        eeTree->Branch("ETmiss_vec",&selectedEvent.ETmiss_vec);

        eeTree->Branch("selMuonSize",&selectedEvent.selMuonSize);
        eeTree->Branch("selElectronSize",&selectedEvent.selElectronSize);
        eeTree->Branch("selLeptonSize",&selectedEvent.selLeptonSize);
        eeTree->Branch("selPhotonSize",&selectedEvent.selPhotonSize);
        eeTree->Branch("selJetSize",&selectedEvent.nselJets);
        eeTree->Branch("selBJetSize",&selectedEvent.nselBJets);

        eeTree->Branch("matchedEleSize",&selectedEvent.matchedEleSize);
        eeTree->Branch("matchedMuSize",&selectedEvent.matchedMuSize);
        eeTree->Branch("matchedLeptonSize",&selectedEvent.matchedLeptonSize);

        eeTree->Branch("isDiElectron",&selectedEvent.isDiElectron);
        eeTree->Branch("isDiMuon",&selectedEvent.isDiMuon);
        eeTree->Branch("isMuonElectron",&selectedEvent.isMuonElectron);
        eeTree->Branch("isElectronMuon",&selectedEvent.isElectronMuon);

        eeTree->Branch("calcHt",&selectedEvent.calcHt);

        eeTree->Branch("lepSF_weight",&selectedEvent.lepSF_weight);
        eeTree->Branch("lepSF_weightUp",&selectedEvent.lepSF_weightUp);
        eeTree->Branch("lepSF_weightDown",&selectedEvent.lepSF_weightDown);

        eeTree->Branch("photonSF_weight",&selectedEvent.photonSF_weight);
        eeTree->Branch("photonSF_weightUp",&selectedEvent.photonSF_weightUp);
        eeTree->Branch("photonSF_weightDown",&selectedEvent.photonSF_weightDown);

        eeTree->Branch("topPt_weight",&selectedEvent.topWeight);
        //eeTree->Branch("topPt_weightUp",&selectedEvent.topWeightUp);
        //eeTree->Branch("topPt_weightDown",&selectedEvent.topWeightDown);

        eeTree->Branch("isr_weight",&selectedEvent.isrWeight);
        eeTree->Branch("isr_weightUp",&selectedEvent.isrWeightUp);
        eeTree->Branch("isr_weightDown",&selectedEvent.isrWeightDown);

        eeTree->Branch("ewk_weight",&selectedEvent.ewkWeight);
        eeTree->Branch("ewk_weightUp",&selectedEvent.ewkWeightUp);
        eeTree->Branch("ewk_weightDown",&selectedEvent.ewkWeightDown);


        eeTree->Branch("electrons",&temp_electrons);
        eeTree->Branch("muons",&temp_muons);
        eeTree->Branch("photons",&temp_photons);
        eeTree->Branch("jets",&temp_jets);
        //eeTree->Branch("genJets",&*genJets);
        eeTree->Branch("genJets",&temp_genJets);
        eeTree->Branch("genParticles",&temp_genParticles);
        eeTree->Branch("intermediateGenParticles",&temp_intermediateGenParticles);
        eeTree->Branch("met",&temp_met);
        eeTree->Branch("metRaw",&temp_metRaw);
        eeTree->Branch("metJESu",&temp_met_JESu);
        eeTree->Branch("metgen",&temp_met_gen);
        eeTree->Branch("metJESd",&temp_met_JESd);
        eeTree->Branch("metJERu",&temp_met_JERu);
        eeTree->Branch("metJERd",&temp_met_JERd);
        eeTree->Branch("nGoodVertices",&temp_nGoodVertices);
        eeTree->Branch("nTracksPV",&temp_nTracksPV);

        eeTree->Branch("pu_weight",&temp_pu_weight);
        eeTree->Branch("pu_weightUp",&temp_pu_weightUp);
        eeTree->Branch("pu_weightDown",&temp_pu_weightDown);
        eeTree->Branch("mc_weight",&temp_mc_weight);
        eeTree->Branch("pdf_weights",&temp_pdf_weights);

        eeTree->Branch("genHt",&temp_genHt);
        eeTree->Branch("ht",&temp_ht);
        eeTree->Branch("nTruePV",&temp_nTruePV);

        eeTree->Branch("nISR",&temp_nISR);
        eeTree->Branch("EWKinoPairPt",&temp_EWKinoPairPt);
        eeTree->Branch("leptonPairPt",&temp_leptonPairPt);
        eeTree->Branch("topPt1",&temp_topPt1);
        eeTree->Branch("topPt2",&temp_topPt2);

        eeTree->Branch("runNo",&temp_runNo);
        eeTree->Branch("lumNo",&temp_lumNo);
        eeTree->Branch("evtNo",&temp_evtNo);

        eeTree->Branch("htTriggered",&selectedEvent.trigHt);
        eeTree->Branch("eeTriggered",&selectedEvent.trigDiEle);
        eeTree->Branch("eeTriggeredMatch",&selectedEvent.trigDiEleMatch);
        eeTree->Branch("mmTriggered",&selectedEvent.trigDiMu);
        eeTree->Branch("mmTriggeredMatch",&selectedEvent.trigDiMuMatch);
        eeTree->Branch("emTriggered",&selectedEvent.trigMuEle);
        eeTree->Branch("emTriggeredMatch",&selectedEvent.trigMuEleMatch);
        eeTree->Branch("metTriggered",&trigMET);

        eeTree->Branch("signal_m1",&temp_signal_m1);
        eeTree->Branch("signal_m2",&temp_signal_m2);
        eeTree->Branch("nBinos",&temp_nBinos);


        //treeFile->Add(eeTree)

}
void treeSplitter::InitHTTree(){

        outputFilenameCompressed = getOutputFilename(inputName,"myTuple");

        //cout<<config_outputfolder<<endl;

        auto saveStr= config_outputfolder+"/ht/"+outputFilenameCompressed;

        if(config_doHt) {
                saveStr= config_outputfolder+"/ht/"+outputFilenameCompressed;
        }
        if(config_doHtPure) {
                saveStr= config_outputfolder+"/htPure/"+outputFilenameCompressed;
        }
        if(config_doMET) {
                saveStr= config_outputfolder+"/met/"+outputFilenameCompressed;
        }

        //TFile newfile(saveStr.c_str(),"RECREATE");

        //TFile* treeFile = new TFile(saveStr.c_str(),"RECREATE");
        treeFile = new TFile(saveStr.c_str(),"RECREATE");

        htTree = new TTree("Tree","Data");


        htTree->Branch("selPhotons",&selectedEvent.selPhotons);
        htTree->Branch("selJets",&selectedEvent.selJets);
        htTree->Branch("selElectrons",&selectedEvent.selElectrons);
        htTree->Branch("selMuons",&selectedEvent.selMuons);

        htTree->Branch("chargeProduct",&selectedEvent.chargeProduct);
        htTree->Branch("mll",&selectedEvent.mll);
        htTree->Branch("miniIso1",&selectedEvent.miniIso1);
        htTree->Branch("miniIso2",&selectedEvent.miniIso2);
        htTree->Branch("deltaRll",&selectedEvent.deltaRll);

        htTree->Branch("pt1",&selectedEvent.pt1);
        htTree->Branch("pt2",&selectedEvent.pt2);
        htTree->Branch("phi1",&selectedEvent.phi1);
        htTree->Branch("phi2",&selectedEvent.phi2);
        htTree->Branch("eta1",&selectedEvent.eta1);
        htTree->Branch("eta2",&selectedEvent.eta2);
        htTree->Branch("l1",&selectedEvent.l1);
        htTree->Branch("l2",&selectedEvent.l2);

        htTree->Branch("negMuons",&temp_negMuons);
        htTree->Branch("posMuons",&temp_posMuons);
        htTree->Branch("negElectrons",&temp_negElectrons);
        htTree->Branch("posElectrons",&temp_posElectrons);

        htTree->Branch("countNegCharge",&temp_countNegCharge);
        htTree->Branch("countPosCharge",&temp_countPosCharge);

        htTree->Branch("MT2",&selectedEvent.MT2_val);

        htTree->Branch("genPhotonVeto",&selectedEvent.evtHasGenPhotonVeto);

        //htTree->Branch("genZToLLCounter",&genZToLLCounter);

        htTree->Branch("ETmiss",&selectedEvent.ETmiss);
        htTree->Branch("ETmiss_vec",&selectedEvent.ETmiss_vec);

        htTree->Branch("selMuonSize",&selectedEvent.selMuonSize);
        htTree->Branch("selElectronSize",&selectedEvent.selElectronSize);
        htTree->Branch("selLeptonSize",&selectedEvent.selLeptonSize);
        htTree->Branch("selPhotonSize",&selectedEvent.selPhotonSize);
        htTree->Branch("selJetSize",&selectedEvent.nselJets);
        htTree->Branch("selBJetSize",&selectedEvent.nselBJets);

        htTree->Branch("matchedEleSize",&selectedEvent.matchedEleSize);
        htTree->Branch("matchedMuSize",&selectedEvent.matchedMuSize);
        htTree->Branch("matchedLeptonSize",&selectedEvent.matchedLeptonSize);

        htTree->Branch("isDiElectron",&selectedEvent.isDiElectron);
        htTree->Branch("isDiMuon",&selectedEvent.isDiMuon);
        htTree->Branch("isMuonElectron",&selectedEvent.isMuonElectron);
        htTree->Branch("isElectronMuon",&selectedEvent.isElectronMuon);

        htTree->Branch("calcHt",&selectedEvent.calcHt);

        htTree->Branch("lepSF_weight",&selectedEvent.lepSF_weight);
        htTree->Branch("lepSF_weightUp",&selectedEvent.lepSF_weightUp);
        htTree->Branch("lepSF_weightDown",&selectedEvent.lepSF_weightDown);

        htTree->Branch("photonSF_weight",&selectedEvent.photonSF_weight);
        htTree->Branch("photonSF_weightUp",&selectedEvent.photonSF_weightUp);
        htTree->Branch("photonSF_weightDown",&selectedEvent.photonSF_weightDown);

        htTree->Branch("topPt_weight",&selectedEvent.topWeight);
        //htTree->Branch("topPt_weightUp",&selectedEvent.topWeightUp);
        //htTree->Branch("topPt_weightDown",&selectedEvent.topWeightDown);

        htTree->Branch("isr_weight",&selectedEvent.isrWeight);
        htTree->Branch("isr_weightUp",&selectedEvent.isrWeightUp);
        htTree->Branch("isr_weightDown",&selectedEvent.isrWeightDown);

        htTree->Branch("ewk_weight",&selectedEvent.ewkWeight);
        htTree->Branch("ewk_weightUp",&selectedEvent.ewkWeightUp);
        htTree->Branch("ewk_weightDown",&selectedEvent.ewkWeightDown);


        htTree->Branch("electrons",&temp_electrons);
        htTree->Branch("muons",&temp_muons);
        htTree->Branch("photons",&temp_photons);
        htTree->Branch("jets",&temp_jets);
        //htTree->Branch("genJets",&*genJets);
        htTree->Branch("genJets",&temp_genJets);
        htTree->Branch("genParticles",&temp_genParticles);
        htTree->Branch("intermediateGenParticles",&temp_intermediateGenParticles);
        htTree->Branch("met",&temp_met);
        htTree->Branch("metRaw",&temp_metRaw);
        htTree->Branch("metJESu",&temp_met_JESu);
        htTree->Branch("metgen",&temp_met_gen);
        htTree->Branch("metJESd",&temp_met_JESd);
        htTree->Branch("metJERu",&temp_met_JERu);
        htTree->Branch("metJERd",&temp_met_JERd);
        htTree->Branch("nGoodVertices",&temp_nGoodVertices);
        htTree->Branch("nTracksPV",&temp_nTracksPV);

        htTree->Branch("pu_weight",&temp_pu_weight);
        htTree->Branch("pu_weightUp",&temp_pu_weightUp);
        htTree->Branch("pu_weightDown",&temp_pu_weightDown);
        htTree->Branch("mc_weight",&temp_mc_weight);
        htTree->Branch("pdf_weights",&temp_pdf_weights);

        htTree->Branch("genHt",&temp_genHt);
        htTree->Branch("ht",&temp_ht);
        htTree->Branch("nTruePV",&temp_nTruePV);

        htTree->Branch("nISR",&temp_nISR);
        htTree->Branch("EWKinoPairPt",&temp_EWKinoPairPt);
        htTree->Branch("leptonPairPt",&temp_leptonPairPt);
        htTree->Branch("topPt1",&temp_topPt1);
        htTree->Branch("topPt2",&temp_topPt2);

        htTree->Branch("runNo",&temp_runNo);
        htTree->Branch("lumNo",&temp_lumNo);
        htTree->Branch("evtNo",&temp_evtNo);

        htTree->Branch("htTriggered",&selectedEvent.trigHt);
        htTree->Branch("metTriggered",&trigMET);
        htTree->Branch("eeTriggered",&selectedEvent.trigDiEle);
        htTree->Branch("eeTriggeredMatch",&selectedEvent.trigDiEleMatch);
        htTree->Branch("mmTriggered",&selectedEvent.trigDiMu);
        htTree->Branch("mmTriggeredMatch",&selectedEvent.trigDiMuMatch);
        htTree->Branch("emTriggered",&selectedEvent.trigMuEle);
        htTree->Branch("emTriggeredMatch",&selectedEvent.trigMuEleMatch);

        htTree->Branch("signal_m1",&temp_signal_m1);
        htTree->Branch("signal_m2",&temp_signal_m2);
        htTree->Branch("nBinos",&temp_nBinos);

}






void treeSplitter::SlaveBegin(TTree *tree)
{
}



Bool_t treeSplitter::Process(Long64_t entry){

        float tempPercentage = (float) entry/ (float)nEntries;
        if(!(abs(config_eventpercentage-100.)<0.1)) {
                if (tempPercentage>config_eventpercentage/100.) {
                        return kTRUE;
                }
        }



        fReader.SetLocalEntry(entry);


        // genHist=TH1F("gen", ";p_{T}^{#gamma} (GeV)", 5000, 0, 5000);
        // recoHist=TH1F("reco", ";p_{T}^{#gamma} (GeV)", 5000, 0, 5000);

        // float countGen=0;
        // float countReco=0;

        // cout<<"----"<<endl;
        for (vector<tree::IntermediateGenParticle>::iterator it = intermediateGenParticles->begin(); it != intermediateGenParticles->end(); it++) {
                for (vector<tree::GenParticle>::iterator itD = it->daughters.begin(); itD != it->daughters.end(); itD++) {
                        if((abs(it->pdgId)>10000)&&(abs(itD->pdgId)==22)&&(itD->p.Pt()>20)) {
                                // cout<<"gen:"<<itD->p.Pt()<<endl;
                                countGen+=1;
                                genHist.Fill(itD->p.Pt());
                                for (vector<tree::Photon>::iterator it2 = photons->begin(); it2 !=photons->end(); it2++) {
                                        // if((itD->p.DeltaR(it2->p)<0.3)) {
                                        // if((it2->p.Pt()>20)&&(!(it2->hasPixelSeed))&& (fabs(it2->p.Eta())<1.4442) && (it2->isLoose)) {
                                        if((it2->p.Pt()>20)&&(!(it2->hasPixelSeed))&& (it2->isLoose)) {
                                                // if((abs(itD->p.Pt()-it2->p.Pt())/itD->p.Pt())<1.) {
                                                // cout<<"reco:"<<it2->p.Pt()<<endl;
                                                countReco+=1;
                                                recoHist.Fill(it2->p.Pt());
                                        }

                                        // }
                                }
                        }
                }
        }
        // cout<<countReco/countGen<<endl;
        // cout<<"----"<<endl;


        trigMET=false;

        //totalWeight = *mc_weight * *pu_weight;

        int progress = tempPercentage*100.;
        if(entry%100000==0) {

                std::cout<<"[";
                for(int i=0; i<100; i++)
                        if(i<progress)
                                std::cout<<'=';
                        else if(i==progress)
                                std::cout<<'>';
                        else
                                std::cout<<' ';
                std::cout<<"] "<<progress<<" %"<<" "<<getOutputFilename(inputName,"myTuple")<<'\r';
                std::cout.flush();

        }


        string cutFlowName = "TreeWriter/hCutFlow";

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

                //cout<<cutFlowName<<endl;

                //if(!cutFlowMap[cutFlowName]){
                if(!(cutFlowMap.count(cutFlowName)>0)) {
                        //cout<<"adding "<<cutFlowName<<" Histogram"<<endl;
                        cutFlowMap[cutFlowName]=cutFlow;
                }

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


        ClearVariables();
        if(config_doHt||config_doMET||config_doHtPure) {
                if(CheckParticles()) {
                        if(CleaningTriggerStudies()) {

                                if(!isTotalSignal) {
                                        weightHisto.Fill("none",1.);
                                        weightHisto.Fill("pu",1.**pu_weight);
                                        weightHisto.Fill("pu_mc",1.**pu_weight**mc_weight);
                                        weightHisto.Fill("pu_mc_nisr",1.**pu_weight**mc_weight*selectedEvent.isrWeight);
                                        weightHisto.Fill("pu_mc_nisrUp",1.**pu_weight**mc_weight*selectedEvent.isrWeightUp);
                                        weightHisto.Fill("pu_mc_nisrDn",1.**pu_weight**mc_weight*selectedEvent.isrWeightDown);
                                        weightHisto.Fill("pu_mc_toppt",1.**pu_weight**mc_weight*selectedEvent.topWeight);
                                        weightHisto.Fill("pu_mc_ewk",1.**pu_weight**mc_weight*selectedEvent.ewkWeight);
                                        weightHisto.Fill("pu_mc_ewkUp",1.**pu_weight**mc_weight*selectedEvent.ewkWeightUp);
                                        weightHisto.Fill("pu_mc_ewkDn",1.**pu_weight**mc_weight*selectedEvent.ewkWeightDown);
                                        //for(int i=6; i<(117); i++){
                                        for(unsigned int i=0; i<pdf_weights->size(); i++) {
                                                string tempPDFstring="pu_mc_pdf_";
                                                tempPDFstring=tempPDFstring+to_string(i);
                                                weightHisto.Fill(tempPDFstring.c_str(), 1.**pu_weight**mc_weight*pdf_weights->at(i));
                                        }
                                }else{
                                        string tempWeightHistName="";
                                        tempWeightHistName=tempWeightHistName+"_"+to_string(*signal_m1);
                                        if(*signal_m2>0) {
                                                tempWeightHistName=tempWeightHistName+"_"+to_string(*signal_m2);
                                        }
                                        weightHistoMap[tempWeightHistName].Fill("none",1.);
                                        weightHistoMap[tempWeightHistName].Fill("pu",1.**pu_weight);
                                        weightHistoMap[tempWeightHistName].Fill("pu_mc",1.**pu_weight**mc_weight);
                                        weightHistoMap[tempWeightHistName].Fill("pu_mc_nisr",1.**pu_weight**mc_weight*selectedEvent.isrWeight);
                                        weightHistoMap[tempWeightHistName].Fill("pu_mc_nisrUp",1.**pu_weight**mc_weight*selectedEvent.isrWeightUp);
                                        weightHistoMap[tempWeightHistName].Fill("pu_mc_nisrDn",1.**pu_weight**mc_weight*selectedEvent.isrWeightDown);
                                        weightHistoMap[tempWeightHistName].Fill("pu_mc_toppt",1.**pu_weight**mc_weight*selectedEvent.topWeight);
                                        weightHistoMap[tempWeightHistName].Fill("pu_mc_ewk",1.**pu_weight**mc_weight*selectedEvent.ewkWeight);
                                        weightHistoMap[tempWeightHistName].Fill("pu_mc_ewkUp",1.**pu_weight**mc_weight*selectedEvent.ewkWeightUp);
                                        weightHistoMap[tempWeightHistName].Fill("pu_mc_ewkDn",1.**pu_weight**mc_weight*selectedEvent.ewkWeightDown);
                                        for(unsigned int i=0; i<pdf_weights->size(); i++) {
                                                string tempPDFstring="pu_mc_pdf_";
                                                tempPDFstring=tempPDFstring+to_string(i);
                                                weightHistoMap[tempWeightHistName].Fill(tempPDFstring.c_str(), 1.**pu_weight**mc_weight*pdf_weights->at(i));
                                        }
                                }



                                temp_electrons=*electrons;
                                temp_muons=*muons;
                                temp_photons=*photons;
                                temp_jets=*jets;

                                temp_genJets=*genJets;
                                temp_genParticles = *genParticles;
                                temp_intermediateGenParticles = *intermediateGenParticles;
                                temp_met = *met;
                                temp_metRaw = *metRaw;
                                temp_met_JESu = *met_JESu;
                                temp_met_gen = *met_gen;
                                temp_met_JESd = *met_JESd;
                                temp_met_JERu = *met_JERu;
                                temp_met_JERd = *met_JERd;
                                temp_nGoodVertices = *nGoodVertices;
                                temp_nTracksPV=*nTracksPV;

                                temp_pu_weight = *pu_weight;
                                temp_pu_weightUp = puWeighterUp.getWeight(*nTruePV);
                                temp_pu_weightDown = puWeighterDown.getWeight(*nTruePV);
                                temp_mc_weight = *mc_weight;
                                temp_pdf_weights = *pdf_weights;

                                temp_genHt = *genHt;
                                temp_ht = *ht;
                                temp_nTruePV=*nTruePV;

                                temp_nISR = *nISR;
                                temp_EWKinoPairPt = *EWKinoPairPt;
                                temp_leptonPairPt=*leptonPairPt;
                                temp_topPt1=*topPt1;
                                temp_topPt2=*topPt2;

                                temp_runNo=*runNo;
                                temp_lumNo=*lumNo;
                                temp_evtNo=*evtNo;

                                temp_nBinos=*nBinos;
                                temp_signal_m1=*signal_m1;
                                temp_signal_m2=*signal_m2;

                                htTree->Fill();
                        }
                }
        }else{
                float temp_isrWeight=1.;
                float temp_isrWeightUp=1.;
                float temp_isrWeightDown=1.;
                if(doWeights_nISR) {
                        temp_isrWeight = isrReweighting(*nISR,false);
                        temp_isrWeightUp= temp_isrWeight+isrReweighting(*nISR,true);
                        temp_isrWeightDown=temp_isrWeight-isrReweighting(*nISR,true);
                }
                float temp_ewkWeight=1.;
                float temp_ewkWeightUp=1.;
                float temp_ewkWeightDown=1.;
                if(doWeights_EWKinoPairPt) {
                        temp_ewkWeight = isrReweightingEWK(*EWKinoPairPt,false);
                        temp_ewkWeightUp= temp_ewkWeight + abs(1.-isrReweightingEWK(*EWKinoPairPt,false));
                        temp_ewkWeightDown=temp_ewkWeight - abs(1.-isrReweightingEWK(*EWKinoPairPt,false));
                }
                float temp_topWeight=1.;
                if(doWeights_TopPt) {
                        float temp_topWeight=topPtReweighting(*topPt1,*topPt2);
                }

                if(!isTotalSignal) {
                        weightHisto.Fill("none",1.);
                        weightHisto.Fill("pu",1.**pu_weight);
                        weightHisto.Fill("pu_mc",1.**pu_weight**mc_weight);
                        weightHisto.Fill("pu_mc_nisr",1.**pu_weight**mc_weight*temp_isrWeight);
                        weightHisto.Fill("pu_mc_nisrUp",1.**pu_weight**mc_weight*temp_isrWeightUp);
                        weightHisto.Fill("pu_mc_nisrDn",1.**pu_weight**mc_weight*temp_isrWeightDown);
                        //weightHisto.Fill("pu_mc_toppt",1.**pu_weight**mc_weight*selectedEvent.topWeight);
                        weightHisto.Fill("pu_mc_toppt",1.**pu_weight**mc_weight*temp_topWeight);
                        weightHisto.Fill("pu_mc_ewk",1.**pu_weight**mc_weight*temp_ewkWeight);
                        weightHisto.Fill("pu_mc_ewkUp",1.**pu_weight**mc_weight*temp_ewkWeightUp);
                        weightHisto.Fill("pu_mc_ewkDn",1.**pu_weight**mc_weight*temp_ewkWeightDown);
                        for(unsigned int i=0; i<pdf_weights->size(); i++) {
                                string tempPDFstring="pu_mc_pdf_";
                                tempPDFstring=tempPDFstring+to_string(i);
                                weightHisto.Fill(tempPDFstring.c_str(), 1.**pu_weight**mc_weight*pdf_weights->at(i));
                        }
                }else{
                        string tempWeightHistName="";
                        tempWeightHistName=tempWeightHistName+"_"+to_string(*signal_m1);
                        if(*signal_m2>0) {
                                tempWeightHistName=tempWeightHistName+"_"+to_string(*signal_m2);
                        }
                        weightHistoMap[tempWeightHistName].Fill("none",1.);
                        weightHistoMap[tempWeightHistName].Fill("pu",1.**pu_weight);
                        weightHistoMap[tempWeightHistName].Fill("pu_mc",1.**pu_weight**mc_weight);
                        weightHistoMap[tempWeightHistName].Fill("pu_mc_nisr",1.**pu_weight**mc_weight*temp_isrWeight);
                        weightHistoMap[tempWeightHistName].Fill("pu_mc_nisrUp",1.**pu_weight**mc_weight*temp_isrWeightUp);
                        weightHistoMap[tempWeightHistName].Fill("pu_mc_nisrDn",1.**pu_weight**mc_weight*temp_isrWeightDown);
                        //weightHistoMap[tempWeightHistName].Fill("pu_mc_toppt",1.**pu_weight**mc_weight*topWeight);
                        weightHistoMap[tempWeightHistName].Fill("pu_mc_toppt",1.**pu_weight**mc_weight*temp_topWeight);
                        weightHistoMap[tempWeightHistName].Fill("pu_mc_ewk",1.**pu_weight**mc_weight*temp_ewkWeight);
                        weightHistoMap[tempWeightHistName].Fill("pu_mc_ewkUp",1.**pu_weight**mc_weight*temp_ewkWeightUp);
                        weightHistoMap[tempWeightHistName].Fill("pu_mc_ewkDn",1.**pu_weight**mc_weight*temp_ewkWeightDown);
                        for(unsigned int i=0; i<pdf_weights->size(); i++) {
                                string tempPDFstring="pu_mc_pdf_";
                                tempPDFstring=tempPDFstring+to_string(i);
                                weightHistoMap[tempWeightHistName].Fill(tempPDFstring.c_str(), 1.**pu_weight**mc_weight*pdf_weights->at(i));
                        }
                }



                if(CheckParticles()) {
                        if(Cleaning()) {

                                //if(!isTotalSignal){
                                //weightHisto.Fill("none",1.);
                                //weightHisto.Fill("pu",1.**pu_weight);
                                //weightHisto.Fill("pu_mc",1.**pu_weight**mc_weight);
                                //weightHisto.Fill("pu_mc_nisr",1.**pu_weight**mc_weight*selectedEvent.isrWeight);
                                //weightHisto.Fill("pu_mc_nisrUp",1.**pu_weight**mc_weight*selectedEvent.isrWeightUp);
                                //weightHisto.Fill("pu_mc_nisrDn",1.**pu_weight**mc_weight*selectedEvent.isrWeightDown);
                                //weightHisto.Fill("pu_mc_toppt",1.**pu_weight**mc_weight*selectedEvent.topWeight);
                                //weightHisto.Fill("pu_mc_ewk",1.**pu_weight**mc_weight*selectedEvent.ewkWeight);
                                //weightHisto.Fill("pu_mc_ewkUp",1.**pu_weight**mc_weight*selectedEvent.ewkWeightUp);
                                //weightHisto.Fill("pu_mc_ewkDn",1.**pu_weight**mc_weight*selectedEvent.ewkWeightDown);
                                //for(unsigned int i=0; i<pdf_weights->size(); i++){
                                //string tempPDFstring="pu_mc_pdf_";
                                //tempPDFstring=tempPDFstring+to_string(i);
                                //weightHisto.Fill(tempPDFstring.c_str(), 1.**pu_weight**mc_weight*pdf_weights->at(i));
                                //}
                                //}else{
                                //string tempWeightHistName="";
                                //tempWeightHistName=tempWeightHistName+"_"+to_string(*signal_m1);
                                //if(*signal_m2>0){
                                //tempWeightHistName=tempWeightHistName+"_"+to_string(*signal_m2);
                                //}
                                //weightHistoMap[tempWeightHistName].Fill("none",1.);
                                //weightHistoMap[tempWeightHistName].Fill("pu",1.**pu_weight);
                                //weightHistoMap[tempWeightHistName].Fill("pu_mc",1.**pu_weight**mc_weight);
                                //weightHistoMap[tempWeightHistName].Fill("pu_mc_nisr",1.**pu_weight**mc_weight*selectedEvent.isrWeight);
                                //weightHistoMap[tempWeightHistName].Fill("pu_mc_nisrUp",1.**pu_weight**mc_weight*selectedEvent.isrWeightUp);
                                //weightHistoMap[tempWeightHistName].Fill("pu_mc_nisrDn",1.**pu_weight**mc_weight*selectedEvent.isrWeightDown);
                                //weightHistoMap[tempWeightHistName].Fill("pu_mc_toppt",1.**pu_weight**mc_weight*selectedEvent.topWeight);
                                //weightHistoMap[tempWeightHistName].Fill("pu_mc_ewk",1.**pu_weight**mc_weight*selectedEvent.ewkWeight);
                                //weightHistoMap[tempWeightHistName].Fill("pu_mc_ewkUp",1.**pu_weight**mc_weight*selectedEvent.ewkWeightUp);
                                //weightHistoMap[tempWeightHistName].Fill("pu_mc_ewkDn",1.**pu_weight**mc_weight*selectedEvent.ewkWeightDown);
                                //for(unsigned int i=0; i<pdf_weights->size(); i++){
                                //string tempPDFstring="pu_mc_pdf_";
                                //tempPDFstring=tempPDFstring+to_string(i);
                                //weightHistoMap[tempWeightHistName].Fill(tempPDFstring.c_str(), 1.**pu_weight**mc_weight*pdf_weights->at(i));
                                //}
                                //}


                                temp_electrons=*electrons;
                                temp_muons=*muons;
                                temp_photons=*photons;
                                temp_jets=*jets;

                                temp_genJets=*genJets;
                                temp_genParticles = *genParticles;
                                temp_intermediateGenParticles = *intermediateGenParticles;
                                temp_met = *met;
                                temp_metRaw = *metRaw;
                                temp_met_JESu = *met_JESu;
                                temp_met_gen = *met_gen;
                                temp_met_JESd = *met_JESd;
                                temp_met_JERu = *met_JERu;
                                temp_met_JERd = *met_JERd;
                                temp_nGoodVertices = *nGoodVertices;
                                temp_nTracksPV=*nTracksPV;

                                temp_pu_weight = *pu_weight;
                                temp_pu_weightUp = puWeighterUp.getWeight(*nTruePV);
                                temp_pu_weightDown = puWeighterDown.getWeight(*nTruePV);
                                temp_mc_weight = *mc_weight;
                                temp_pdf_weights = *pdf_weights;

                                temp_genHt = *genHt;
                                temp_ht = *ht;
                                temp_nTruePV=*nTruePV;

                                temp_nISR = *nISR;
                                temp_EWKinoPairPt = *EWKinoPairPt;
                                temp_leptonPairPt=*leptonPairPt;
                                temp_topPt1=*topPt1;
                                temp_topPt2=*topPt2;

                                temp_runNo=*runNo;
                                temp_lumNo=*lumNo;
                                temp_evtNo=*evtNo;

                                temp_nBinos=*nBinos;
                                temp_signal_m1=*signal_m1;
                                temp_signal_m2=*signal_m2;

                                eeTree->Fill();
                        }
                }
        }






        return kTRUE;
}























































void treeSplitter::InitScaleFactorsFinal(){

        //DiEleWeighterID=Weighter("scaleFactors/Alternative/ScaleFactorElectronID_MVATight.root","EGamma_SF2D");// pt vs |eta|
        //DiEleWeighterID.fillOverflow2d();
        //DiEleWeighterTrack=Weighter("scaleFactors/Alternative/egammaEffi.txt_EGM2D.root","EGamma_SF2D"); //pt vs |eta|
        //DiEleWeighterTrack.fillOverflow2d();
        DiEleWeighterID=Weighter("newScaleFactors/electron/MVATight/ScaleFactorElectronID_MVATight.root","EGamma_SF2D");// pt vs |eta|
        DiEleWeighterID.fillOverflow2d();
        DiEleWeighterTrack=Weighter("newScaleFactors/electron/MVATight/egammaEffi.txt_EGM2D.root","EGamma_SF2D"); //pt vs |eta|
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
void treeSplitter::InitScaleFactorsFinalMedium(){

        DiEleWeighterID=Weighter("newScaleFactors/electron/mediumCut/mediumCut.root","EGamma_SF2D");// pt(y) vs |eta|(x) EGAMMAPOG
        DiEleWeighterID.fillOverflow2d();
        DiEleWeighterTrack=Weighter("newScaleFactors/electron/mediumCut/egammaEffi.txt_EGM2D.root","EGamma_SF2D"); // pt(y) vs |eta|(x) EGAMMAPOG
        DiEleWeighterTrack.fillOverflow2d();

        //DiEleWeighterID=Weighter("newScaleFactors/electron/mediumSUSY/scaleFactors.root","GsfElectronToCutBasedSpring15M");// eta(y) vs pt(x)  SUSY PAG
        //DiEleWeighterID.fillOverflow2d();
        //DiEleWeighterIso=Weighter("newScaleFactors/electron/mediumSUSY/scaleFactors.root","MVAVLooseElectronToMini");// eta(y) vs pt(x)  SUSY PAG
        //DiEleWeighterIso.fillOverflow2d();
        //DiEleWeighterConv=Weighter("newScaleFactors/electron/mediumSUSY/scaleFactors.root","MVATightElectronToConvVetoIHit0");// eta(y) vs pt(x)  SUSY PAG
        //DiEleWeighterConv.fillOverflow2d();
        //DiEleWeighterTrack=Weighter("newScaleFactors/electron/mediumCut/egammaEffi.txt_EGM2D.root","EGamma_SF2D"); //pt vs |eta| EGAMMAPOG
        //DiEleWeighterTrack.fillOverflow2d();

        DiMuWeighterID = Weighter("newScaleFactors/muon/scaleFactor_muon_ID.root","SF");// |eta| vs pt
        DiMuWeighterID.fillOverflow2d();
        DiMuWeighterIso = Weighter("newScaleFactors/muon/scaleFactor_muon_iso.root","SF");// |eta| vs pt
        DiMuWeighterIso.fillOverflow2d();
        DiMuWeighterIP2D = Weighter("newScaleFactors/muon/scaleFactor_muon_IP2D.root","SF");// |eta| vs pt
        DiMuWeighterIP2D.fillOverflow2d();
        DiMuWeighterSIP3D = Weighter("newScaleFactors/muon/scaleFactor_muon_SIP3D.root","SF");// |eta| vs pt
        DiMuWeighterSIP3D.fillOverflow2d();
        //DiMuWeighterTrack = Weighter("newScaleFactors/muon/TrackScaleFactorsMuons.root","muonTrackScaleFactorEtaHisto");// |eta|
        //DiMuWeighterTrack.fillOverflow2d();

        PhotonIDWeighter = Weighter("newScaleFactors/photon/scaleFactor_photon_ID.root","EGamma_SF2D"); // pt vs eta
        PhotonIDWeighter.fillOverflow2d();

        PhotonConversionWeighter = Weighter("newScaleFactors/photon/scaleFactor_photon_conv.root","Scaling_Factors_HasPix_R9 Inclusive"); // pt vs eta
        PhotonConversionWeighter.fillOverflow2d();
}




float treeSplitter::GetScaleFactorAndErrorFinal(float pt, float eta,bool isFastSim,bool isEle,bool error){
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

        if(isEle) {
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


                DiEleSF = DiEleSF * DiEleWeighterID.getWeight(absEta,pt); // pt(y) vs |eta|(x) EGAMMAPOG
                //DiEleSF = DiEleSF * DiEleWeighterID.getWeight(pt,absEta);     // eta(y) vs pt(x)  SUSY PAG
                //DiEleSF = DiEleSF * DiEleWeighterIso.getWeight(pt,absEta);   // eta(y) vs pt(x)  SUSY PAG
                //DiEleSF = DiEleSF * DiEleWeighterConv.getWeight(pt,absEta);  // eta(y) vs pt(x)  SUSY PAG
                DiEleSF = DiEleSF * DiEleWeighterTrack.getWeight(absEta,pt);

                DiEleSFErr = DiEleWeighterID.getError(absEta,pt);
                DiEleSFErr = pow((pow(DiEleSFErr,2.) + pow(DiEleWeighterTrack.getError(absEta,pt),2.)),0.5);
                //DiEleSFErr = pow((pow(DiEleSFErr,2.) + pow(DiEleWeighterIso.getError(pt,absEta),2.)),0.5);
                //DiEleSFErr = pow((pow(DiEleSFErr,2.) + pow(DiEleWeighterConv.getError(pt,absEta),2.)),0.5);
                DiEleSFErr = pow((pow(DiEleSFErr,2.) + pow((0.01*DiEleSF*(pt<20. && pt>75.)),2.)),0.5); //additional 1% uncertainty for all (normal just pt<20 and >75)

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
                //DiMuSF = DiMuSF * DiMuWeighterTrack.getWeight(absEta);

                DiMuSFErr = 0.03;

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
        float SFError=1.;
        //if(isFastSim){
        //SF = isEle ?  FastSimDiEleSF : FastSimDiMuSF;
        //}else{
        SF = isEle ?  DiEleSF : DiMuSF;
        SFError=isEle ? DiEleSFErr : DiMuSFErr;
        //}
        if(error) {
                return SFError;
        }else{
                return SF;
        }
}


float treeSplitter::GetScaleFactorAndErrorPhotons(vector<selPhoton>& vecGamma,bool error){
        float totalSF = 1.;
        float totalSFError = 0.;
        for (vector<selPhoton>::iterator it = vecGamma.begin(); it != vecGamma.end(); it++) {
                float tempPt=it->p.Pt();
                float tempEta=it->p.Eta();
                //0.9938 +- 0.0119 https://twiki.cern.ch/twiki/pub/CMS/EgammaIDRecipesRun2/EleVetoScalingFactors_Moriond17.pdf
                totalSF = totalSF * PhotonIDWeighter.getWeight(fabs(tempEta),tempPt)*PhotonConversionWeighter.getWeight(fabs(tempEta),tempPt);
                totalSFError = pow(pow(totalSFError,2.) + pow(PhotonIDWeighter.getError(fabs(tempEta),tempPt),2.) + pow(PhotonConversionWeighter.getError(fabs(tempEta),tempPt),2.),0.5);
        }
        if(error) {
                return totalSFError;
        }else{
                return totalSF;
        }
}




bool treeSplitter::GenPhotonVeto(const int a){
        bool isGenPhoton=false;
        bool vetoPt130ZG=false;
        if(!isData) {
                if(noPromptPhotons || isZGammaInclusive) {
                        for (vector<tree::GenParticle>::iterator it = genParticles->begin(); it != genParticles->end(); it++) {


                                if(it->pdgId==22) {

                                        if(isZGammaInclusive) {
                                                if(it->p.Pt()>=130.) {
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
                                                isGenPhoton=isGenPhoton || it->isPrompt;
                                                break;
                                        case 4:
                                                isGenPhoton=isGenPhoton || !(it->promptStatus==DIRECTPROMPT);
                                                break;
                                        case 5:
                                                isGenPhoton=isGenPhoton || ((it->statusID==1) && (it->isPrompt));
                                                break;
                                        case 6:
                                                isGenPhoton=isGenPhoton || ((it->isPromptFinalState));
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


bool treeSplitter::CheckParticles(){
        return ((muons->size()>=2.) || (electrons->size()>=2.) || ((muons->size()>=1)&&(electrons->size()>=1)) );
}
bool treeSplitter::Check2Ele(){
        return (electrons->size()>=2.);
}
bool treeSplitter::Check2Mu(){
        return (muons->size()>=2.);
}
bool treeSplitter::CheckEMu(){
        return ((muons->size()>=1.) && (electrons->size()>=1.));
}













bool treeSplitter::Cleaning() //return if event is valid (2 Leptons exist)
{
        if(!isSignal) {
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


        if (isData) {
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

                for(auto& it : *electrons) {
                        if ((it.p.Pt()-tol) > pt1_) {
                                pt1_ = it.p.Pt();
                                leptonFlavor1 = "E";
                        }
                }
                for(auto& it : *muons) {
                        if ((it.p.Pt()-tol) > pt1_) {
                                pt1_ = it.p.Pt();
                                leptonFlavor1 = "M";
                        }
                }

                for(auto& it : *electrons) {
                        if ((it.p.Pt() < (pt1_-tol)) && ((it.p.Pt()-tol) > pt2_)) {
                                pt2_ = it.p.Pt();
                                leptonFlavor2 = "E";
                        }
                }
                for(auto& it : *muons) {
                        if ((it.p.Pt() < (pt1_-tol)) && ((it.p.Pt()-tol) > pt2_)) {
                                pt2_ = it.p.Pt();
                                leptonFlavor2 = "M";
                        }
                }
                if (leptonFlavor1 == "E" && leptonFlavor2 == "E") {
                        if(isEESelection && Check2Ele() && selectedEvent.trigDiEle) { CalculateVariables(electrons->at(0),electrons->at(1),E); return true;}
                }
                else if (leptonFlavor1 == "M" && leptonFlavor2 == "M") {
                        if(isMuMuSelection && Check2Mu() && selectedEvent.trigDiMu) { CalculateVariables(muons->at(0),muons->at(1),M); return true;}
                }
                else if (leptonFlavor1 == "E" && leptonFlavor2 == "M") {
                        if(isEMuSelection && CheckEMu() && selectedEvent.trigMuEle) { CalculateVariables(electrons->at(0),muons->at(0),EM); return true;}
                }
                else if (leptonFlavor1 == "M" && leptonFlavor2 == "E") {
                        if(isEMuSelection && CheckEMu() && selectedEvent.trigMuEle) { CalculateVariables(muons->at(0),electrons->at(0),ME); return true;}
                }

                return false;
        }else{
                auto missingET = *met;

                float pt1_ = 0.;
                float pt2_ = 0.;
                float tol = 0.0001;

                std::string leptonFlavor1 = "N";
                std::string leptonFlavor2 = "N";

                for(auto& it : *electrons) {
                        if ((it.p.Pt()-tol) > pt1_) {
                                pt1_ = it.p.Pt();
                                leptonFlavor1 = "E";
                        }
                }
                for(auto& it : *muons) {
                        if ((it.p.Pt()-tol) > pt1_) {
                                pt1_ = it.p.Pt();
                                leptonFlavor1 = "M";
                        }
                }
                for(auto& it : *electrons) {
                        if ((it.p.Pt() < (pt1_-tol)) && ((it.p.Pt()-tol) > pt2_)) {
                                pt2_ = it.p.Pt();
                                leptonFlavor2 = "E";
                        }
                }
                for(auto& it : *muons) {
                        if ((it.p.Pt() < (pt1_-tol)) && ((it.p.Pt()-tol) > pt2_)) {
                                pt2_ = it.p.Pt();
                                leptonFlavor2 = "M";
                        }
                }
                //cout<<leptonFlavor1<<leptonFlavor2<<endl;
                if (leptonFlavor1 == "E" && leptonFlavor2 == "E") {
                        if(Check2Ele() && selectedEvent.trigDiEle) {CalculateVariables(electrons->at(0),electrons->at(1),E); return true;}
                }
                else if (leptonFlavor1 == "M" && leptonFlavor2 == "M") {
                        if(Check2Mu() && selectedEvent.trigDiMu) { CalculateVariables(muons->at(0),muons->at(1),M); return true;}
                }
                else if (leptonFlavor1 == "E" && leptonFlavor2 == "M") {
                        if(CheckEMu() && selectedEvent.trigMuEle) { CalculateVariables(electrons->at(0),muons->at(0),EM); return true;}
                }
                else if (leptonFlavor1 == "M" && leptonFlavor2 == "E") {
                        if(CheckEMu() && selectedEvent.trigMuEle) { CalculateVariables(muons->at(0),electrons->at(0),ME); return true;}
                }
                return false;
        }
        return false;
}


bool treeSplitter::CleaningTriggerStudies() //return if event is valid (2 Leptons exist)
{
        if(!isSignal) {
                selectedEvent.trigDiEle = *hlt_ele17_ele12_iso || *hlt_ele23_ele12_iso || *hlt_doubleEle33 || *hlt_doubleEle33_mw;
                selectedEvent.trigDiMu  = *hlt_mu17_mu8_iso || *hlt_mu17_tkMu8_iso || *hlt_mu17_mu8_iso_dz || *hlt_mu17_tkMu8_iso_dz
                                          || *hlt_tkMu17_tkMu8_iso_dz || *hlt_mu27_tkMu8 || *hlt_mu30_tkMu11;
                selectedEvent.trigMuEle = *hlt_mu17_ele12_iso || *hlt_mu23_ele8_iso || *hlt_mu23_ele8_iso_dz || *hlt_mu23_ele12_iso
                                          || *hlt_mu23_ele12_iso_dz || *hlt_mu8_ele17_iso || *hlt_mu8_ele23_iso
                                          || *hlt_mu8_ele23_iso_dz || *hlt_mu12_ele23_iso || *hlt_mu12_ele23_iso_dz
                                          || *hlt_mu30_ele30 || *hlt_mu33_ele33;
                selectedEvent.trigHt    = *hlt_ht200 || *hlt_ht250 || *hlt_ht300 || *hlt_ht350 || *hlt_ht400
                                          || *hlt_ht475 || *hlt_ht600 || *hlt_ht650 || *hlt_ht800;
                trigMET    = *hlt_met110 || *hlt_met120 || *hlt_met170_HBHE || *hlt_met170_Jet || *hlt_met170_Noise
                             || *hlt_met170_Not || *hlt_met300 || *hlt_met400 || *hlt_met500 || *hlt_met600;
        }else{
                selectedEvent.trigDiEle = true;
                selectedEvent.trigDiMu = true;
                selectedEvent.trigMuEle = true;
                selectedEvent.trigHt = true;
                trigMET = true;
        }


        if (isData) {
                //find Selection (Dataset)
                bool isEESelection = inputName.find("DoubleEG") != string::npos;
                bool isMuMuSelection = inputName.find("DoubleMuon") != string::npos;
                bool isEMuSelection = inputName.find("MuonEG") != string::npos;
                bool isHTSelection = inputName.find("JetHT")!= string::npos;
                bool isMETSelection = inputName.find("MET")!= string::npos;
                auto missingET = *met;


                float pt1_ = 0.;
                float pt2_ = 0.;
                float tol = 0.0001;

                std::string leptonFlavor1 = "N";
                std::string leptonFlavor2 = "N";

                for(auto& it : *electrons) {
                        if ((it.p.Pt()-tol) > pt1_) {
                                pt1_ = it.p.Pt();
                                leptonFlavor1 = "E";
                        }
                }
                for(auto& it : *muons) {
                        if ((it.p.Pt()-tol) > pt1_) {
                                pt1_ = it.p.Pt();
                                leptonFlavor1 = "M";
                        }
                }
                for(auto& it : *electrons) {
                        if ((it.p.Pt() < (pt1_-tol)) && ((it.p.Pt()-tol) > pt2_)) {
                                pt2_ = it.p.Pt();
                                leptonFlavor2 = "E";
                        }
                }
                for(auto& it : *muons) {
                        if ((it.p.Pt() < (pt1_-tol)) && ((it.p.Pt()-tol) > pt2_)) {
                                pt2_ = it.p.Pt();
                                leptonFlavor2 = "M";
                        }
                }



                if (leptonFlavor1 == "E" && leptonFlavor2 == "E") {
                        if(config_doHt) {
                                if(isHTSelection && Check2Ele() && selectedEvent.trigHt) { CalculateVariables(electrons->at(0),electrons->at(1),E); return true;}
                        }else{
                                if(isMETSelection && Check2Ele() && trigMET) { CalculateVariables(electrons->at(0),electrons->at(1),E); return true;}
                        }
                }
                else if (leptonFlavor1 == "M" && leptonFlavor2 == "M") {
                        if(config_doHt) {
                                if(isHTSelection && Check2Mu() && selectedEvent.trigHt) { CalculateVariables(muons->at(0),muons->at(1),M); return true;}
                        }else{
                                if(isMETSelection && Check2Mu() && trigMET) { CalculateVariables(muons->at(0),muons->at(1),M); return true;}
                        }
                }
                else if (leptonFlavor1 == "E" && leptonFlavor2 == "M") {
                        if(config_doHt) {
                                if(isHTSelection && CheckEMu() && selectedEvent.trigHt) { CalculateVariables(electrons->at(0),muons->at(0),EM); return true;}
                        }else{
                                if(isMETSelection && CheckEMu() && trigMET) { CalculateVariables(electrons->at(0),muons->at(0),EM); return true;}
                        }
                }
                else if (leptonFlavor1 == "M" && leptonFlavor2 == "E") {
                        if(config_doHt) {
                                if(isHTSelection && CheckEMu() && selectedEvent.trigHt) { CalculateVariables(muons->at(0),electrons->at(0),ME); return true;}
                        }else{
                                if(isMETSelection && CheckEMu() && trigMET) { CalculateVariables(muons->at(0),electrons->at(0),ME); return true;}
                        }
                }
                return false;

        }else{
                auto missingET = *met;

                float pt1_ = 0.;
                float pt2_ = 0.;
                float tol = 0.0001;

                std::string leptonFlavor1 = "N";
                std::string leptonFlavor2 = "N";

                for(auto& it : *electrons) {
                        if ((it.p.Pt()-tol) > pt1_) {
                                pt1_ = it.p.Pt();
                                leptonFlavor1 = "E";
                        }
                }
                for(auto& it : *muons) {
                        if ((it.p.Pt()-tol) > pt1_) {
                                pt1_ = it.p.Pt();
                                leptonFlavor1 = "M";
                        }
                }
                for(auto& it : *electrons) {
                        if ((it.p.Pt() < (pt1_-tol)) && ((it.p.Pt()-tol) > pt2_)) {
                                pt2_ = it.p.Pt();
                                leptonFlavor2 = "E";
                        }
                }
                for(auto& it : *muons) {
                        if ((it.p.Pt() < (pt1_-tol)) && ((it.p.Pt()-tol) > pt2_)) {
                                pt2_ = it.p.Pt();
                                leptonFlavor2 = "M";
                        }
                }
                if (leptonFlavor1 == "E" && leptonFlavor2 == "E") {
                        if(config_doHtPure) {
                                if(Check2Ele()) {CalculateVariables(electrons->at(0),electrons->at(1),E); return true;}
                        }else{
                                if(config_doMET) {
                                        if(Check2Ele() && trigMET) {CalculateVariables(electrons->at(0),electrons->at(1),E); return true;}
                                }else{
                                        if(Check2Ele() && selectedEvent.trigHt) {CalculateVariables(electrons->at(0),electrons->at(1),E); return true;}
                                }
                        }
                }
                else if (leptonFlavor1 == "M" && leptonFlavor2 == "M") {
                        if(config_doHtPure) {
                                if(Check2Mu()) { CalculateVariables(muons->at(0),muons->at(1),M); return true;}
                        }else{
                                if(config_doMET) {
                                        if(Check2Mu() && trigMET) { CalculateVariables(muons->at(0),muons->at(1),M); return true;}
                                }else{
                                        if(Check2Mu() && selectedEvent.trigHt) { CalculateVariables(muons->at(0),muons->at(1),M); return true;}
                                }
                        }
                }
                else if (leptonFlavor1 == "E" && leptonFlavor2 == "M") {
                        if(config_doHtPure) {
                                if(CheckEMu()) { CalculateVariables(electrons->at(0),muons->at(0),EM); return true;}
                        }else{
                                if(config_doMET) {
                                        if(CheckEMu() && trigMET) { CalculateVariables(electrons->at(0),muons->at(0),EM); return true;}
                                }else{
                                        if(CheckEMu() && selectedEvent.trigHt) { CalculateVariables(electrons->at(0),muons->at(0),EM); return true;}
                                }
                        }

                }
                else if (leptonFlavor1 == "M" && leptonFlavor2 == "E") {
                        if(config_doHtPure) {
                                if(CheckEMu()) { CalculateVariables(muons->at(0),electrons->at(0),ME); return true;}
                        }else{
                                if(config_doMET) {
                                        if(CheckEMu() && trigMET) { CalculateVariables(muons->at(0),electrons->at(0),ME); return true;}
                                }else{
                                        if(CheckEMu() && selectedEvent.trigHt) { CalculateVariables(muons->at(0),electrons->at(0),ME); return true;}
                                }
                        }
                }

                return false;
        }
        return false;
}




//ELECTRON
//bool treeSplitter::testSelection(const tree::Electron& pa,selectionType selection ,bool leading){
//bool decision=false;
//if (selection==UNCUT){
//decision = pa.isPassConvVeto && pa.passImpactParameter && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && pa.isTightMVA;
//}
//if(selection==SEL || selection==DILEP || selection==ControlRegionZZ || selection==ControlRegionWZ || selection==SEL){
//decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (selectedEvent.deltaRll>0.1) && pa.isTightMVA;

//}
//if(selection==ONZ){
//decision = pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (selectedEvent.deltaRll>0.1) && pa.isTightMVA;
//}
//if(selection==TRIGSEL || selection==TRIGDILEP || selection==TRIGONZ){
//decision =(*ht>200.) && pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>0.) : (pa.p.Pt()>0.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (selectedEvent.deltaRll>0.1) && pa.isTightMVA;
//}
//if(selection==TRIGSEL_ptcuts || selection==TRIGDILEP_ptcuts || selection==TRIGONZ_ptcuts){
//decision =(*ht>200.) && pa.isPassConvVeto && pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.1) && (selectedEvent.deltaRll>0.1) && pa.isTightMVA;
//}
//return decision;
//}
//MUON
//bool treeSplitter::testSelection(const tree::Muon& pa, selectionType selection_, bool leading){
//bool decision=false;
//if (selection_==UNCUT){
//decision = pa.passImpactParameter  && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && (pa.isTight);
//}
//if (selection_==SEL || selection_==DILEP || selection_==ControlRegionZZ || selection_==ControlRegionWZ){
//decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.2) && (pa.isMedium) && (selectedEvent.deltaRll>0.1);
//}
//if (selection_ == ONZ){
//decision = pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.2) && (pa.isMedium) && (selectedEvent.deltaRll>0.1);
//}
//if(selection_==TRIGSEL || selection_==TRIGDILEP || selection_==TRIGONZ){
//decision = (*ht>200.) && pa.passImpactParameter && (leading ? (pa.p.Pt()>0.) : (pa.p.Pt()>0.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (selectedEvent.deltaRll>0.1);
//}
//if(selection_==TRIGSEL_ptcuts || selection_==TRIGDILEP_ptcuts || selection_==TRIGONZ_ptcuts){
//decision = (*ht>200.)&& pa.passImpactParameter && (leading ? (pa.p.Pt()>25.) : (pa.p.Pt()>20.)) && (fabs(pa.p.Eta())<2.4) && ((fabs(pa.p.Eta())<1.4)|| (fabs(pa.p.Eta()))>1.6) && (pa.miniIso<0.2) && (pa.isMedium) && (selectedEvent.deltaRll>0.1);
//}
//return decision;
//}
//PHOTON
bool treeSplitter::testSelection(const selPhoton& pa, selectionType selection){
        bool decision=false;
        if(selection==SEL || selection==DILEP || selection==EGRegression || selection==ControlRegionZZ || selection==TRIGDILEP || selection==TRIGSEL || selection==TRIGONZ) {
                //decision = (pa.p.Pt()>25.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isLoose) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3); //study Delta R cut ;
                decision = (pa.p.Pt()>20.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isLoose) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3); //study Delta R cut ;
                //decision = (pa.p.Pt()>20.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<2.4) && (pa.isLoose) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3); //study Delta R cut ;
                //decision = (pa.p.Pt()>10.) && !(pa.hasPixelSeed)  && (pa.isLoose) ; //study Delta R cut ;
        }
        if(selection==ONZ) {
                //decision = (pa.p.Pt()>25.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isLoose) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3); //study Delta R cut ;
                decision = (pa.p.Pt()>20.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<1.4442) && (pa.isLoose) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3); //study Delta R cut ;
                //decision = (pa.p.Pt()>20.) && !(pa.hasPixelSeed) && (fabs(pa.p.Eta())<2.4) && (pa.isLoose) && (pa.deltaR1>0.3) && (pa.deltaR2>0.3); //study Delta R cut ;
        }
        return decision;
}
//SELMUON
bool treeSplitter::testSelection(const selMuon& pa, selectionType selection_){
        bool decision=false;
        if (selection_==LooseLeptons) {
                decision = pa.passImpactParameter && (pa.p.Pt()>20.) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.2) && (pa.isMedium);
                //decision = pa.passImpactParameter && (pa.p.Pt()>10.) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.2) && (pa.isMedium);
        }
        if (selection_==LooseLeptonsTrigger) {
                //decision = pa.passImpactParameter && (pa.p.Pt()>20.) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.2) && (pa.isMedium);
                decision = pa.passImpactParameter && (pa.p.Pt()>10.) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.2) && (pa.isMedium);
        }
        return decision;
}
//SELELECTRON
bool treeSplitter::testSelection(const selElectron& pa, selectionType selection_){
        bool decision=false;
        if(selection_==LooseLeptons) {
                //decision = pa.isPassConvVeto && pa.passImpactParameter && (pa.p.Pt()>20.) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && pa.isTightMVA;
                decision = pa.isPassConvVeto && pa.passImpactParameter && (pa.p.Pt()>20.) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && pa.isMedium;
                //decision = pa.isPassConvVeto && pa.passImpactParameter && (pa.p.Pt()>10.) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && pa.isMedium;
        }
        if(selection_==LooseLeptonsTrigger) {
                //decision = pa.isPassConvVeto && pa.passImpactParameter && (pa.p.Pt()>20.) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && pa.isTightMVA;
                //decision = pa.isPassConvVeto && pa.passImpactParameter && (pa.p.Pt()>20.) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && pa.isMedium;
                decision = pa.isPassConvVeto && pa.passImpactParameter && (pa.p.Pt()>10.) && (fabs(pa.p.Eta())<2.4) && (pa.miniIso<0.1) && pa.isMedium;
        }
        return decision;
}
//JET
bool treeSplitter::testSelection(const selJet& pa){
        bool decision=false;
        decision = (pa.p.Pt()>30.) && (fabs(pa.p.Eta())<2.4) && (pa.isLoose) && !(pa.hasElectronMatch) && !(pa.hasPhotonMatch) && !(pa.hasMuonMatch); //&& (pa.deltaR1>0.3) && (pa.deltaR2>0.3) //study Delta R cut ;
        return decision;
}

bool treeSplitter::matchLepton(const selElectron& pa){
        bool foundMatch=false;
        if((pa.deltaR1 < 0.05 ) || (pa.deltaR2 < 0.05 )) {
                foundMatch=true;
        }else{
                foundMatch=false;
        }
        return foundMatch;
}
bool treeSplitter::matchLepton(const selMuon& pa){
        bool foundMatch=false;
        //TLorentzVector dummyLepton(0.,0.,0.,0.);
        //dummyLepton.SetPtEtaPhiM(pa.p.Pt(),pa.p.Eta(),pa.p.Phi(),0.51099895e-3);
        if((pa.deltaR1 < 0.05 ) || (pa.deltaR2 < 0.05 )) {
                foundMatch=true;
        }else{
                foundMatch=false;
        }
        return foundMatch;
}


bool treeSplitter::matchRecoPhoton(const selPhoton& pho, const string collection){
        bool matched=false;
        if(collection=="J") {
                //for(vector<tree::Particle>::iterator it = genJets->begin(); it != genJets->end();it++){
                //if(pho.p.DeltaR(it->p)<0.1){
                //matched=matched||true;
                //}
                //}
                for(vector<tree::GenParticle>::iterator it = genParticles->begin(); it != genParticles->end(); it++) {
                        if((abs(it->pdgId)<7)||(abs(it->pdgId)==21)||(abs(it->pdgId)==9))
                                if(pho.p.DeltaR(it->p)<0.1) {
                                        //if((abs(pho.p.Pt()-it->p.Pt())/pho.p.Pt())<0.5){
                                        //if((abs(pho.p.Pt()-it->p.Pt())/pho.p.Pt())<3.){
                                        if((abs(pho.p.Pt()-it->p.Pt())/it->p.Pt())<3.) {
                                                matched=matched||true;
                                        }
                                }
                }
        }else{
                if(collection=="P") {
                        for(vector<tree::GenParticle>::iterator it2 = genParticles->begin(); it2 != genParticles->end(); it2++) {
                                if((abs(it2->pdgId)==22)&&(it2->statusID==1)) {
                                        if(pho.p.DeltaR(it2->p)<0.1) {
                                                //if((abs(pho.p.Pt()-it2->p.Pt())/pho.p.Pt())<0.5){
                                                if((abs(pho.p.Pt()-it2->p.Pt())/it2->p.Pt())<0.5) {
                                                        matched=matched||true;
                                                }
                                        }
                                }
                        }
                }else{
                        if(collection=="E") {
                                for(vector<tree::GenParticle>::iterator it3 = genParticles->begin(); it3 != genParticles->end(); it3++) {
                                        if((abs(it3->pdgId)==11)&&(it3->statusID==1)) {
                                                if(pho.p.DeltaR(it3->p)<0.1) {
                                                        //if((abs(pho.p.Pt()-it3->p.Pt())/pho.p.Pt())<0.5){
                                                        if((abs(pho.p.Pt()-it3->p.Pt())/it3->p.Pt())<0.5) {
                                                                matched=matched||true;
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }
        return matched;
}



void treeSplitter::CalculateVariables(const tree::Lepton& lep1, const tree::Lepton& lep2,const particleType particle=DUMMYPARTICLE){

        //leptons
        if (particle==E) {
                selectedEvent.isDiElectron=true;
        }else{
                selectedEvent.isDiElectron=false;
        }
        if (particle==M) {
                selectedEvent.isDiMuon=true;
        }else{
                selectedEvent.isDiMuon=false;
        }
        if (particle==EM) {
                selectedEvent.isElectronMuon=true;
        }else{
                selectedEvent.isElectronMuon=false;
        }
        if (particle==ME) {
                selectedEvent.isMuonElectron=true;
        }else{
                selectedEvent.isMuonElectron=false;
        }


        if(!isSignal) {
                selectedEvent.trigDiEle = *hlt_ele17_ele12_iso || *hlt_ele23_ele12_iso || *hlt_doubleEle33 || *hlt_doubleEle33_mw;
                selectedEvent.trigDiMu  = *hlt_mu17_mu8_iso || *hlt_mu17_tkMu8_iso || *hlt_mu17_mu8_iso_dz || *hlt_mu17_tkMu8_iso_dz
                                          || *hlt_tkMu17_tkMu8_iso_dz || *hlt_mu27_tkMu8 || *hlt_mu30_tkMu11;
                selectedEvent.trigMuEle = *hlt_mu17_ele12_iso || *hlt_mu23_ele8_iso || *hlt_mu23_ele8_iso_dz || *hlt_mu23_ele12_iso
                                          || *hlt_mu23_ele12_iso_dz || *hlt_mu8_ele17_iso || *hlt_mu8_ele23_iso
                                          || *hlt_mu8_ele23_iso_dz || *hlt_mu12_ele23_iso || *hlt_mu12_ele23_iso_dz
                                          || *hlt_mu30_ele30 || *hlt_mu33_ele33;
                selectedEvent.trigHt    = *hlt_ht200 || *hlt_ht250 || *hlt_ht300 || *hlt_ht350 || *hlt_ht400
                                          || *hlt_ht475 || *hlt_ht600 || *hlt_ht650 || *hlt_ht800;
                trigMET    = *hlt_met110 || *hlt_met120 || *hlt_met170_HBHE || *hlt_met170_Jet || *hlt_met170_Noise
                             || *hlt_met170_Not || *hlt_met300 || *hlt_met400 || *hlt_met500 || *hlt_met600;
        }else{
                selectedEvent.trigDiEle = true;
                selectedEvent.trigDiMu = true;
                selectedEvent.trigMuEle = true;
                selectedEvent.trigHt = true;
                trigMET = true;
        }

        //selectedEvent.trigMuEleMatch=selectedEvent.trigMuEle;
        //trigObjEleMatched=false
        if(!isSignal) {
                if(selectedEvent.isDiElectron) {
                        if(selectedEvent.trigDiEle) {
                                bool matchedL1=false;
                                bool matchedL2=false;
                                for(vector<tree::Particle>::iterator it = trigObj_ele17_ele12->begin(); it != trigObj_ele17_ele12->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_ele23_ele12->begin(); it != trigObj_ele23_ele12->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_ele33_ele33->begin(); it != trigObj_ele33_ele33->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_ele33_ele33_mw->begin(); it != trigObj_ele33_ele33_mw->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                selectedEvent.trigDiEleMatch=(matchedL1 && matchedL2);
                        }
                }

                if(selectedEvent.isDiMuon) {
                        if(selectedEvent.trigDiMu) {
                                bool matchedL1=false;
                                bool matchedL2=false;
                                for(vector<tree::Particle>::iterator it = trigObj_mu17_mu8->begin(); it != trigObj_mu17_mu8->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu17_mu8tk->begin(); it != trigObj_mu17_mu8tk->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu17_mu8_dz->begin(); it != trigObj_mu17_mu8_dz->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu17_mu8tk_dz->begin(); it != trigObj_mu17_mu8tk_dz->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu17tk_mu8tk_dz->begin(); it != trigObj_mu17tk_mu8tk_dz->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu27_mu8->begin(); it != trigObj_mu27_mu8->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu30_mu11->begin(); it != trigObj_mu30_mu11->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                selectedEvent.trigDiMuMatch=(matchedL1 && matchedL2);
                        }
                }
                if(selectedEvent.isMuonElectron||selectedEvent.isElectronMuon) {
                        if(selectedEvent.trigMuEle) {
                                bool matchedL1=false;
                                bool matchedL2=false;
                                for(vector<tree::Particle>::iterator it = trigObj_mu17_ele12_eleLeg->begin(); it != trigObj_mu17_ele12_eleLeg->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu17_ele12_muLeg->begin(); it != trigObj_mu17_ele12_muLeg->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu23_ele8_eleLeg->begin(); it != trigObj_mu23_ele8_eleLeg->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu23_ele8_muLeg->begin(); it != trigObj_mu23_ele8_muLeg->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu23_ele8_dz->begin(); it != trigObj_mu23_ele8_dz->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu23_ele12_eleLeg->begin(); it != trigObj_mu23_ele12_eleLeg->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu23_ele12_muLeg->begin(); it != trigObj_mu23_ele12_muLeg->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu23_ele12_dz->begin(); it != trigObj_mu23_ele12_dz->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu8_ele17_eleLeg->begin(); it != trigObj_mu8_ele17_eleLeg->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu8_ele17_muLeg->begin(); it != trigObj_mu8_ele17_muLeg->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu8_ele23_eleLeg->begin(); it != trigObj_mu8_ele23_eleLeg->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu8_ele23_muLeg->begin(); it != trigObj_mu8_ele23_muLeg->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu8_ele23_dz->begin(); it != trigObj_mu8_ele23_dz->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu12_ele23_muLeg->begin(); it != trigObj_mu12_ele23_muLeg->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu12_ele23_eleLeg->begin(); it != trigObj_mu12_ele23_eleLeg->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu12_ele23_dz->begin(); it != trigObj_mu12_ele23_dz->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu30_ele30_eleLeg->begin(); it != trigObj_mu30_ele30_eleLeg->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu30_ele30_muLeg->begin(); it != trigObj_mu30_ele30_muLeg->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu33_ele33_eleLeg->begin(); it != trigObj_mu33_ele33_eleLeg->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
                                                matchedL2=matchedL2 || true;
                                        }
                                }
                                for(vector<tree::Particle>::iterator it = trigObj_mu33_ele33_muLeg->begin(); it != trigObj_mu33_ele33_muLeg->end(); it++) {
                                        if(lep1.p.DeltaR(it->p)<0.05) {
                                                matchedL1=matchedL1 || true;
                                        }
                                        if(lep2.p.DeltaR(it->p)<0.05) {
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


        float leptonMass1 = (selectedEvent.isDiElectron||selectedEvent.isElectronMuon) ?  0.51099895e-3 : 0.1056584;
        float leptonMass2 = (selectedEvent.isDiMuon||selectedEvent.isMuonElectron) ?  0.1056584 : 0.51099895e-3;

        int charge1 = (int) lep1.charge;
        int charge2 = (int) lep2.charge;

        selectedEvent.l1.vec.SetPtEtaPhiM(lep1.p.Pt(),lep1.p.Eta(),lep1.p.Phi(),leptonMass1);
        selectedEvent.l2.vec.SetPtEtaPhiM(lep2.p.Pt(),lep2.p.Eta(),lep2.p.Phi(),leptonMass2);

        selectedEvent.pt1 = selectedEvent.l1.vec.Pt();
        selectedEvent.pt2 = selectedEvent.l2.vec.Pt();
        selectedEvent.phi1 = selectedEvent.l1.vec.Phi();
        selectedEvent.phi2 = selectedEvent.l2.vec.Phi();
        selectedEvent.eta1 = selectedEvent.l1.vec.Eta();
        selectedEvent.eta2 = selectedEvent.l2.vec.Eta();

        selectedEvent.chargeProduct = charge1*charge2;
        selectedEvent.mll = (selectedEvent.l1.vec+selectedEvent.l2.vec).M();
        selectedEvent.miniIso1 = lep1.miniIso;
        selectedEvent.miniIso2 = lep2.miniIso;
        selectedEvent.deltaRll=selectedEvent.l1.vec.DeltaR(selectedEvent.l2.vec);
        //selectedEvent.ETmiss = mett.p.Pt();
        selectedEvent.ETmiss = met->p.Pt();
        //selectedEvent.ETmiss_vec.SetPtEtaPhiM(met->p.Pt(),met->p.Eta(),met->p.Phi(),0.);;
        selectedEvent.ETmiss_vec.vec.SetPtEtaPhiM(met->p.Pt(),met->p.Eta(),met->p.Phi(),0.);;

        //double pa[3] = {lep1.M(),lep1.Px(),lep1.Py()};
        //double pb[3] = {lep2.M(),lep2.Px(),lep2.Py()};
        //double pmiss[3] = {0.,missing.Px(),missing.Py()};
        //fctMT2_.set_mn(0.);
        //fctMT2_.set_momenta(pa,pb,pmiss,*evtNo);
        //MT2_val = static_cast<float>(fctMT2_.get_mt2());



        selectedEvent.evtHasGenPhotonVeto = GenPhotonVeto(config_veto);

        for (vector<tree::Photon>::iterator it = photons->begin(); it != photons->end(); it++) {
                auto g = *it;
                selPhoton gamma;
                gamma.setAll(g);
                gamma.vec.SetPtEtaPhiM(gamma.p.Pt(),gamma.p.Eta(),gamma.p.Phi(),0.);
                gamma.deltaR1=gamma.vec.DeltaR(selectedEvent.l1.vec);
                gamma.deltaR2=gamma.vec.DeltaR(selectedEvent.l2.vec);

                bool matchedToPhoton=false;
                bool matchedToJet=false;
                bool matchedToElectron=false;
                matchedToPhoton=matchRecoPhoton(gamma,"P");
                matchedToJet=matchRecoPhoton(gamma,"J");
                matchedToElectron=matchRecoPhoton(gamma,"E");

                gamma.matchedToPhoton=matchedToPhoton;
                gamma.matchedToJet=matchedToJet;
                gamma.matchedToElectron=matchedToElectron;

                gamma.scaleFactor=PhotonIDWeighter.getWeight(fabs(gamma.p.Eta()),gamma.p.Pt())*PhotonConversionWeighter.getWeight(fabs(gamma.p.Eta()),gamma.p.Pt());
                gamma.scaleFactorErr=PhotonIDWeighter.getError(fabs(gamma.p.Eta()),gamma.p.Pt())*PhotonConversionWeighter.getWeight(fabs(gamma.p.Eta()),gamma.p.Pt());


                if (testSelection(gamma,SEL)) {
                        selectedEvent.selPhotons.push_back(gamma);
                }
        }
        for (vector<tree::Muon>::iterator it = muons->begin(); it != muons->end(); it++) {
                auto g = *it;
                selMuon muon;
                muon.setAll(g);
                muon.vec.SetPtEtaPhiM(muon.p.Pt(),muon.p.Eta(),muon.p.Phi(),0.1056584);
                //if (testSelection(muon,config_doHt? LooseLeptonsTrigger : LooseLeptons)){
                if (testSelection(muon,(config_doHt||config_doMET||config_doHtPure) ? LooseLeptonsTrigger : LooseLeptons)) {
                        muon.deltaR1 = muon.vec.DeltaR(selectedEvent.l1.vec);
                        muon.deltaR2 = muon.vec.DeltaR(selectedEvent.l2.vec);
                        muon.chargeInt = (int) g.charge;

                        muon.scaleFactor= GetScaleFactorAndErrorFinal(muon.p.Pt(),muon.p.Eta(),false,false);
                        muon.scaleFactorErr= GetScaleFactorAndErrorFinal(muon.p.Pt(),muon.p.Eta(),false,false,true);

                        if(matchLepton(muon)) {
                                muon.matched=true;
                                selectedEvent.matchedMuSize+=1;
                                selectedEvent.matchedLeptonSize+=1;
                        }
                        selectedEvent.selMuons.push_back(muon);
                }
        }
        for (vector<tree::Electron>::iterator it = electrons->begin(); it != electrons->end(); it++) {
                auto g = *it;
                selElectron ele;
                ele.setAll(g);
                ele.vec.SetPtEtaPhiM(ele.p.Pt(),ele.p.Eta(),ele.p.Phi(),0.51099895e-3);
                //if (testSelection(ele,LooseLeptons)){
                if (testSelection(ele,(config_doHt||config_doMET||config_doHtPure) ? LooseLeptonsTrigger : LooseLeptons)) {
                        ele.deltaR1 = ele.vec.DeltaR(selectedEvent.l1.vec);
                        ele.deltaR2 = ele.vec.DeltaR(selectedEvent.l2.vec);
                        ele.chargeInt = (int) g.charge;

                        ele.scaleFactor=GetScaleFactorAndErrorFinal(ele.p.Pt(),ele.p.Eta(),false,true);
                        ele.scaleFactorErr=GetScaleFactorAndErrorFinal(ele.p.Pt(),ele.p.Eta(),false,true,true);

                        if(matchLepton(ele)) {
                                ele.matched=true;
                                selectedEvent.matchedEleSize+=1;
                                selectedEvent.matchedLeptonSize+=1;
                        }
                        selectedEvent.selElectrons.push_back(ele);
                }
        }

        selectedEvent.selPhotonSize=(selectedEvent.selPhotons.size());
        selectedEvent.selMuonSize=(selectedEvent.selMuons.size());
        selectedEvent.selElectronSize=(selectedEvent.selElectrons.size());
        selectedEvent.selLeptonSize=(selectedEvent.selMuonSize+selectedEvent.selElectronSize);

        for (vector<tree::Jet>::iterator it = jets->begin(); it != jets->end(); it++) {
                auto j = *it;
                selJet jet;
                jet.setAll(j);
                jet.vec.SetPtEtaPhiM(jet.p.Pt(),jet.p.Eta(),jet.p.Phi(),0.);
                jet.deltaR1=jet.vec.DeltaR(selectedEvent.l1.vec);
                jet.deltaR2=jet.vec.DeltaR(selectedEvent.l2.vec);
                jet.bTag = (jet.bDiscriminator>0.8484); //MEDIUM WP include scale factors
                if (testSelection(jet)) {
                        selectedEvent.selJets.push_back(jet);
                        selectedEvent.nselJets+=1;
                        if(jet.bTag) selectedEvent.nselBJets+=1;
                        selectedEvent.calcHt += jet.vec.Pt();
                }
        }

        float tempTotalWeightLep=1.;
        float tempTotalWeightLepErr=0.;
        float tempTotalWeightPhoton=1.;
        float tempTotalWeightPhotonErr=0.;

        for (vector<selElectron>::iterator it = selectedEvent.selElectrons.begin(); it != selectedEvent.selElectrons.end(); it++) {
                tempTotalWeightLep=tempTotalWeightLep*it->scaleFactor;
                tempTotalWeightLepErr=pow(pow(tempTotalWeightLepErr,2.)+pow(it->scaleFactorErr,2.),0.5);
        }
        for (vector<selMuon>::iterator it = selectedEvent.selMuons.begin(); it != selectedEvent.selMuons.end(); it++) {
                tempTotalWeightLep=tempTotalWeightLep*it->scaleFactor;
                tempTotalWeightLepErr=pow(pow(tempTotalWeightLepErr,2.)+pow(it->scaleFactorErr,2.),0.5);
        }
        for (vector<selPhoton>::iterator it = selectedEvent.selPhotons.begin(); it != selectedEvent.selPhotons.end(); it++) {
                tempTotalWeightPhoton=tempTotalWeightPhoton*it->scaleFactor;
                tempTotalWeightPhotonErr=pow(pow(tempTotalWeightPhotonErr,2.)+pow(it->scaleFactorErr,2.),0.5);
        }

        if(!isData) {
                selectedEvent.lepSF_weight=tempTotalWeightLep;
                selectedEvent.lepSF_weightUp=tempTotalWeightLep*(1.+tempTotalWeightLepErr);
                selectedEvent.lepSF_weightDown=tempTotalWeightLep*(1.-tempTotalWeightLepErr);

                selectedEvent.photonSF_weight = tempTotalWeightPhoton;
                selectedEvent.photonSF_weightUp=tempTotalWeightPhoton*(1.+tempTotalWeightPhotonErr);
                selectedEvent.photonSF_weightDown=tempTotalWeightPhoton*(1.-tempTotalWeightPhotonErr);
        }else{
                selectedEvent.lepSF_weight=1.;
                selectedEvent.lepSF_weightUp=1.;
                selectedEvent.lepSF_weightDown=1.;

                selectedEvent.photonSF_weight = 1.;
                selectedEvent.photonSF_weightUp=1.;
                selectedEvent.photonSF_weightDown=1.;
        }

        selectedEvent.topWeight=1.;
        if(doWeights_TopPt) {
                selectedEvent.topWeight = topPtReweighting(*topPt1,*topPt2);
                //selectedEvent.topWeightUp;
                //selectedEvent.topWeightDown;
        }
        selectedEvent.isrWeight=1.;
        selectedEvent.isrWeightUp=1.;
        selectedEvent.isrWeightDown=1.;
        if(doWeights_nISR) {
                selectedEvent.isrWeight = isrReweighting(*nISR,false);
                selectedEvent.isrWeightUp= selectedEvent.isrWeight+isrReweighting(*nISR,true);
                selectedEvent.isrWeightDown=selectedEvent.isrWeight-isrReweighting(*nISR,true);
        }
        selectedEvent.ewkWeight=1.;
        selectedEvent.ewkWeightUp=1.;
        selectedEvent.ewkWeightDown=1.;
        if(doWeights_EWKinoPairPt) {
                selectedEvent.ewkWeight = isrReweightingEWK(*EWKinoPairPt,false);
                selectedEvent.ewkWeightUp= selectedEvent.ewkWeight + abs(1.-isrReweightingEWK(*EWKinoPairPt,false));
                selectedEvent.ewkWeightDown=selectedEvent.ewkWeight - abs(1.-isrReweightingEWK(*EWKinoPairPt,false));
        }


        for (vector<selElectron>::iterator it = selectedEvent.selElectrons.begin(); it != selectedEvent.selElectrons.end(); it++) {
                if (it->chargeInt <0) {
                        temp_countNegCharge+=1;
                        temp_negElectrons.push_back(*it);
                }else{
                        temp_countPosCharge+=1;
                        temp_posElectrons.push_back(*it);
                }
        }
        for (vector<selMuon>::iterator it = selectedEvent.selMuons.begin(); it != selectedEvent.selMuons.end(); it++) {
                if (it->chargeInt <0) {
                        temp_countNegCharge+=1;
                        temp_negMuons.push_back(*it);
                }else{
                        temp_countPosCharge+=1;
                        temp_posMuons.push_back(*it);
                }
        }


        double pa[3] = {(selectedEvent.l1.vec+selectedEvent.l2.vec).M(),(selectedEvent.l1.vec+selectedEvent.l2.vec).Px(),(selectedEvent.l1.vec+selectedEvent.l2.vec).Py()};
        double pb[3] = {0.,0.,0.};

        if(selectedEvent.selPhotonSize>0) {
                pb[0] = {selectedEvent.selPhotons.at(0).vec.M()};
                pb[1] = {selectedEvent.selPhotons.at(0).vec.Px()};
                pb[2] = {selectedEvent.selPhotons.at(0).vec.Py()};
        }

        //cout<<pb[0]<<" "<<pb[1]<<" "<<pb[2]<<endl;

        double pmiss[3] = {0.,selectedEvent.ETmiss_vec.vec.Px(),selectedEvent.ETmiss_vec.vec.Py()};
        fctMT2_.set_mn(0.);
        fctMT2_.set_momenta(pa,pb,pmiss,*evtNo);
        selectedEvent.MT2_val = static_cast<float>(fctMT2_.get_mt2());




}

void treeSplitter::ClearVariables(){

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

        selectedEvent.trigDiEleMatch=false;
        selectedEvent.trigDiMuMatch=false;
        selectedEvent.trigMuEleMatch=false;

        trigMET=false;

        int charge1 = 0;
        int charge2 = 0;

        selectedEvent.l1.vec.SetPtEtaPhiM(0.,0.,0.,0.);
        selectedEvent.l2.vec.SetPtEtaPhiM(0.,0.,0.,0.);

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
        selectedEvent.ETmiss_vec.vec.SetPtEtaPhiM(0.,0.,0.,0.);

        //double pa[3] = {selectedEvent.l1.M(),selectedEvent.l1.Px(),selectedEvent.l1.Py()};
        //double pb[3] = {selectedEvent.l2.M(),selectedEvent.l2.Px(),selectedEvent.l2.Py()};

        //double pmiss[3] = {0.,selectedEvent.ETmiss_vec.Px(),selectedEvent.ETmiss_vec.Py()};
        //fctMT2_.set_mn(0.);
        //fctMT2_.set_momenta(pa,pb,pmiss,*evtNo);
        //selectedEvent.MT2_val = static_cast<float>(fctMT2_.get_mt2());

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
        selectedEvent.selPhotonSize=0;

        selectedEvent.nselJets=0;
        selectedEvent.nselBJets=0;

        selectedEvent.calcHt=0.;

        selectedEvent.lepSF_weight=1.;
        selectedEvent.lepSF_weightUp=1.;
        selectedEvent.lepSF_weightDown=1.;

        selectedEvent.photonSF_weight = 1.;
        selectedEvent.photonSF_weightUp=1.;
        selectedEvent.photonSF_weightDown=1.;


        selectedEvent.topWeight = 1.;



        selectedEvent.isrWeight = 1.;
        selectedEvent.isrWeightUp= 1.;
        selectedEvent.isrWeightDown=1.;


        selectedEvent.ewkWeight = 1.;
        selectedEvent.ewkWeightUp= 1.;
        selectedEvent.ewkWeightDown=1.;


        temp_countNegCharge=0;
        temp_countPosCharge=0;

        temp_posElectrons.clear();
        temp_posMuons.clear();
        temp_negElectrons.clear();
        temp_negMuons.clear();

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
//template<typename T>
//void save2File(const map<string,map<Histograms2D,T>>& hMaps, TFile& file)
//{
//for (auto& hMapIt : hMaps) {
//if (!file.Get(hMapIt.first.c_str())) {
//file.mkdir(hMapIt.first.c_str());
//}
//file.cd(hMapIt.first.c_str());
//for (auto& h : hMapIt.second) {
//h.second.Write(histoNames2D[h.first].c_str(), TObject::kWriteDelete);
//}
//file.cd();
//}
//}
//template<typename T>
//void save2File(const map<string,map<string,map<Histograms1D,T>>>& hMaps, TFile& file)
//{
//for (auto& hMapIt : hMaps) {
//if (!file.Get(hMapIt.first.c_str())) {
//file.mkdir(hMapIt.first.c_str());
//}
//file.cd(hMapIt.first.c_str());
//for(auto& hMapIt2 : hMapIt.second){
//if (!file.Get(hMapIt2.first.c_str())) {
//file.mkdir((hMapIt.first+"/"+hMapIt2.first).c_str());
//}
//file.cd((hMapIt.first+"/"+hMapIt2.first).c_str());
//for (auto h: hMapIt2.second){
//h.second.Write(histoNames[h.first].c_str(), TObject::kWriteDelete);
//}
//file.cd();
//}
//}
//}

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
                saveName=saveName.replace(0,11,"");
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


void treeSplitter::SaveTree(){

        outputFilenameCompressed = getOutputFilename(inputName,"myTuple");

        auto saveStr= config_outputfolder+"/"+outputFilenameCompressed;

        if(!isTotalSignal) {
                cutFlow.Write();
                cutFlow.SetDirectory(0);

                //cout<<weightHisto.Integral()<<endl;

                weightHisto.Write();
                weightHisto.SetDirectory(0);

        }else{
                save2File(cutFlowMap,*treeFile);
                save2File(weightHistoMap,*treeFile);
                cutFlow.SetDirectory(0);
        }

        eeTree->Write();
        eeTree->SetDirectory(0);
        treeFile->Write();
        treeFile->Close();

        cout << "Created " << saveStr << " in " << (time(NULL) - startTime)/60 << " min" << endl;
}
void treeSplitter::SaveHTTree(){

        outputFilenameCompressed = getOutputFilename(inputName,"myTuple");

        auto saveStr= config_outputfolder+"/ht/"+outputFilenameCompressed;

        if(config_doHt) {
                saveStr= config_outputfolder+"/ht/"+outputFilenameCompressed;
        }
        if(config_doHtPure) {
                saveStr= config_outputfolder+"/htPure/"+outputFilenameCompressed;
        }
        if(config_doMET) {
                saveStr= config_outputfolder+"/met/"+outputFilenameCompressed;
        }


        if(!isTotalSignal) {
                cutFlow.Write();
                cutFlow.SetDirectory(0);

                //cout<<weightHisto.Integral()<<endl;

                weightHisto.Write();
                weightHisto.SetDirectory(0);

        }else{
                save2File(cutFlowMap,*treeFile);
                cutFlow.SetDirectory(0);
        }


        htTree->Write();
        htTree->SetDirectory(0);
        treeFile->Write();
        treeFile->Close();

        cout << "Created " << saveStr << " in " << (time(NULL) - startTime)/60 << " min" << endl;
}






void treeSplitter::Terminate()
{
        //auto outputName = "output"+config_outputfolder+"/"+getOutputFilename(inputName);

        //TFile file(outputName.c_str(), "RECREATE");
        //if(! isTotalSignal){
        //cutFlow.Write("hCutFlow");
        //}
        //file.Close();
        cout<<countReco/countGen<<endl;

        recoHist.Write();
        recoHist.SetDirectory(0);
        genHist.Write();
        genHist.SetDirectory(0);

        if(config_doHt||config_doMET||config_doHtPure) {
                SaveHTTree();
        }else{
                SaveTree();
        }
        //save2File(cutFlowMap,*treeFile);
        //SaveHTTree();

        //cout << "Created " << outputName << " in " << (time(NULL) - startTime)/60 << " min" << endl;
}
