import ROOT
import copy
import auxiliary as aux
import os.path
import style


#path = "../minimal/output/"
#path = "../myAnalyzer/output/"
#path = "../myAnalyzer/output_medium/"
path = "../myAnalyzer/output_mediumPOG/"
#path = "../myAnalyzer/output_mediumPOG/ht/"

#path = "../minimal/output_signalScan/"

#path = "../minimal/output_uncorrected/"
#path = "../minimal/output_EGRegression/"

#path = "../minimal/output_mllWeight/"

#path = "../minimal/output_noVeto/"

#path = "../minimal/output_FSRVeto/"
#path = "../minimal/output_FSR+HardVeto/"
#path = "../minimal/output_01/"

class Dataset:
    names = []
    files = []
    xsecs = []
    ngens = []
    lname = []

    color = None
    label = ""

    def __radd__( self, dset ):
        if not dset: return self

        out = copy.deepcopy(self)

        out.names.extend( dset.names )
        out.files.extend( dset.files )
        out.xsecs.extend( dset.xsecs )
        out.ngens.extend( dset.ngens )
        out.label += " + "+dset.label
        out.lname.extend( dset.lname )

        if not self.color:
            out.color = dset.color
        return out

    def __init__( self, n, xsec=-1, col=ROOT.kBlack, fullname="", ngen=-1 ):
        fname = path+n+"_hists.root"
        if xsec == -1: xsec = aux.getXsecFromName( n )
        if ngen == -1 and os.path.isfile(fname): ngen = aux.getNgen( fname )
        self.names = [ n ]
        self.files = [ fname ]
        self.xsecs = [ xsec ]
        self.ngens = [ ngen ]
        self.color = col
        self.label = n
        self.lname = [ fullname ]

    def mcm( self ):
        for fullname in self.lname:
            print "https://cms-pdmv.cern.ch/mcm/requests?page=0&dataset_name="+fullname


    def __str__( self ):
        return "Dataset: " + self.label + "\ncolor: " + str(self.color) + \
            "\nnames: "+", ".join( self.names ) + \
            "\nfiles: "+", ".join( self.files ) + \
            "\nxsecs: "+", ".join( str(i) for i in self.xsecs ) + \
            "\nngens: "+", ".join( str(i) for i in self.ngens )

    def getNgenFromFile( self ):
        for ngen, file in zip(self.ngens, self.files):
            print file, ngen, aux.getNgen( file )


    def getHist( self, name ):
        h0 = None
        for i in range( len(self.files) ):
            h = aux.getFromFile( self.files[i], name )
            if isinstance( h, ROOT.TH1 ):
                if self.xsecs[i]:
                    if style.additionalPoissonUncertainty:
                        aux.addPoissonUncertainty(h)
                    h.Scale( aux.intLumi * self.xsecs[i] / self.ngens[i] )
                h.SetLineColor( self.color )
                h.SetMarkerColor( self.color )
            if h0: h0.Add( h )
            else: h0 = h
        return h0

    def getLatexTableHeader( self ):
        return "\\begin{tabular}{l|r|r}\nPrimary Dataset & cross section (pb) & effective Luminosity (/fb) \\\\\\hline"

    def getLatexTableLine( self ):
        # full samplename & xsec [pb] & effective Luminosity [/fb]
        out = ""
        for lname, xsec, ngen in zip(self.lname,self.xsecs,self.ngens):
            out += "{} & {:.3f} & {:.2f} \\\\\n".format( lname.replace("_","\\_"),xsec,0.001*ngen/xsec )
        return out

###############################################################################
# Data
###############################################################################
dataHt = Dataset("JetHT_Run2016B-03Feb2017_ver2-v2", 0, ROOT.kBlack ) \
    + Dataset("JetHT_Run2016C-03Feb2017-v1", 0, ROOT.kBlack ) \
    + Dataset("JetHT_Run2016D-03Feb2017-v1", 0, ROOT.kBlack ) \
    + Dataset("JetHT_Run2016E-03Feb2017-v1", 0, ROOT.kBlack ) \
    + Dataset("JetHT_Run2016F-03Feb2017-v1", 0, ROOT.kBlack ) \
    + Dataset("JetHT_Run2016G-03Feb2017-v1", 0, ROOT.kBlack ) \
    + Dataset("JetHT_Run2016H-03Feb2017_ver2-v1", 0, ROOT.kBlack ) \
    + Dataset("JetHT_Run2016H-03Feb2017_ver3-v1", 0, ROOT.kBlack )
dataHt.label = "Data (JetHt)"

dataDoubleMuon = Dataset("DoubleMuon_Run2016B-03Feb2017_ver2-v2",0,ROOT.kBlack) \
    + Dataset("DoubleMuon_Run2016C-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleMuon_Run2016D-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleMuon_Run2016E-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleMuon_Run2016F-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleMuon_Run2016G-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleMuon_Run2016H-03Feb2017_ver2-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleMuon_Run2016H-03Feb2017_ver3-v1", 0, ROOT.kBlack) 
dataDoubleMuon.label = "Data (#mu#mu)"

dataDoubleEG = Dataset("DoubleEG_Run2016B-03Feb2017_ver2-v2",0,ROOT.kBlack) \
    + Dataset("DoubleEG_Run2016C-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleEG_Run2016D-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleEG_Run2016E-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleEG_Run2016F-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleEG_Run2016G-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleEG_Run2016H-03Feb2017_ver3-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleEG_Run2016H-03Feb2017_ver2-v1", 0, ROOT.kBlack) 
dataDoubleEG.label = "Data (ee)"

dataMuonEG = Dataset("MuonEG_Run2016B-03Feb2017_ver2-v2",0,ROOT.kBlack) \
    + Dataset("MuonEG_Run2016C-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("MuonEG_Run2016D-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("MuonEG_Run2016E-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("MuonEG_Run2016F-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("MuonEG_Run2016G-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("MuonEG_Run2016H-03Feb2017_ver3-v1", 0, ROOT.kBlack) \
    + Dataset("MuonEG_Run2016H-03Feb2017_ver2-v1", 0, ROOT.kBlack) 
dataMuonEG.label = "Data (e#mu)"

dataDoubleSF = Dataset("DoubleEG_Run2016B-03Feb2017_ver2-v2",0,ROOT.kBlack) \
    + Dataset("DoubleEG_Run2016C-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleEG_Run2016D-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleEG_Run2016E-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleEG_Run2016F-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleEG_Run2016G-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleEG_Run2016H-03Feb2017_ver3-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleEG_Run2016H-03Feb2017_ver2-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleMuon_Run2016B-03Feb2017_ver2-v2",0,ROOT.kBlack) \
    + Dataset("DoubleMuon_Run2016C-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleMuon_Run2016D-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleMuon_Run2016E-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleMuon_Run2016F-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleMuon_Run2016G-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleMuon_Run2016H-03Feb2017_ver3-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleMuon_Run2016H-03Feb2017_ver2-v1", 0, ROOT.kBlack)  
dataDoubleSF.label = "Data (ee+#mu#mu)"

dataLL = Dataset("DoubleEG_Run2016B-03Feb2017_ver2-v2",0,ROOT.kBlack) \
    + Dataset("DoubleEG_Run2016C-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleEG_Run2016D-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleEG_Run2016E-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleEG_Run2016F-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleEG_Run2016G-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleEG_Run2016H-03Feb2017_ver3-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleEG_Run2016H-03Feb2017_ver2-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleMuon_Run2016B-03Feb2017_ver2-v2",0,ROOT.kBlack) \
    + Dataset("DoubleMuon_Run2016C-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleMuon_Run2016D-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleMuon_Run2016E-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleMuon_Run2016F-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleMuon_Run2016G-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleMuon_Run2016H-03Feb2017_ver3-v1", 0, ROOT.kBlack) \
    + Dataset("DoubleMuon_Run2016H-03Feb2017_ver2-v1", 0, ROOT.kBlack)  \
    + Dataset("MuonEG_Run2016B-03Feb2017_ver2-v2",0,ROOT.kBlack) \
    + Dataset("MuonEG_Run2016C-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("MuonEG_Run2016D-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("MuonEG_Run2016E-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("MuonEG_Run2016F-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("MuonEG_Run2016G-03Feb2017-v1", 0, ROOT.kBlack) \
    + Dataset("MuonEG_Run2016H-03Feb2017_ver3-v1", 0, ROOT.kBlack) \
    + Dataset("MuonEG_Run2016H-03Feb2017_ver2-v1", 0, ROOT.kBlack) 
dataLL.label = "Data (ee+#mu#mu+e#mu)"



###############################################################################
# Simulation
###############################################################################

###############################################################################
## WWGamma
wwgamma = Dataset("WWG-amcatnlo_ext", 0.2147, ROOT.kCyan+3)
wwgamma.label = "WW#gamma"
###############################################################################
## WW
ww = Dataset("WWTo2L2Nu", 12.178 , ROOT.kCyan-2)
ww.label = "WW"
###############################################################################
## WZ
#wz = Dataset("WZTo3LNu_ext", 4.42965, ROOT.kMagenta+4)
wz = Dataset("WZTo3LNu_Total", 4.42965, ROOT.kAzure-5)
wz.label = "WZ"

###############################################################################
## ZZ
#zz = Dataset("ZZTo2L2Nu_powheg", 0.5644, ROOT.kYellow) \
#zz =    Dataset("ZZTo2L2Nu_powheg_ext1", 0.5644, ROOT.kYellow)
zz =    Dataset("ZZTo2L2Nu_powheg_Total", 0.5644, ROOT.kYellow)
zz.label = "ZZ(#rightarrowll#nu#nu)"
###############################################################################
## ZZ(4l)
#zz4l =    Dataset("ZZTo4L_powheg_ext1", 1.212 , ROOT.kYellow+2)
#zz4l =    Dataset("ZZTo4L_powheg", 1.212 , ROOT.kYellow+2)
zz4l =    Dataset("ZZTo4L_powheg_Total", 1.212 , ROOT.kOrange-2)
zz4l.label = "ZZ(#rightarrow4l)"

###############################################################################
## ZGamma
#zgamma = Dataset("ZGTo2LG", 117.864 , ROOT.kGreen)\
    #+Dataset("ZGTo2LG_ext", 117.864 , ROOT.kGreen) \
    #+Dataset("ZGTo2LG_PtG-130", 0.1404 , ROOT.kGreen)
#zgamma = Dataset("ZGTo2LG_ext", (117.864-0.1404) , ROOT.kGreen-3) \
#zgamma = Dataset("ZGTo2LG_Total", (117.864-0.1404) , ROOT.kGreen-3) \
#zgamma = Dataset("ZGTo2LG_Total", (117.864) , ROOT.kGreen-3) \
zgamma = Dataset("ZGTo2LG_Total", (117.864-0.1404) , ROOT.kGreen-3) \
    +Dataset("ZGTo2LG_PtG-130", 0.1404 , ROOT.kGreen-3) #\
#zgamma = Dataset("ZGTo2LG_Total", (117.864) , ROOT.kGreen-3) #\
    #+Dataset("ZGTo2LG_PtG-130", 0.1404 , ROOT.kGreen-3) #\
#zgamma = Dataset("ZGTo2LG", 117.864 , ROOT.kGreen) #\
zgamma.label = "Z#gamma"

###############################################################################
## WZGamma
wzgamma = Dataset("WZG-amcatnlo", 0.04123 , ROOT.kBlue+3)
wzgamma.label = "WZ#gamma"

###############################################################################
## TTGamma_dilept
## NLO cross section calculated "by hand" from sample
ttgamma = Dataset("TTGamma_Dilept-amcatnlo", 0.6352 , ROOT.kRed+1)
#ttgamma = Dataset("TTGamma_Dilept-amcatnlo", 0.7056 , ROOT.kRed+1)
#ttgamma = Dataset("TTGamma_Dilept-amcatnlo", 0.899 , ROOT.kRed+1)
#ttgamma = Dataset("TTGamma_Dilept-amcatnlo", 0.899*1.6 , ROOT.kRed+1)
#ttgamma = Dataset("TTGamma_Dilept-amcatnlo", 0.6343 , ROOT.kRed+1)
ttgamma.label = "t#bar{t}#gamma"

###############################################################################
## DY+jets
#DYjets = Dataset("DYJetsToLL_M-50_ext", 1921.8*3.  , ROOT.kBlue-7)
#DYjets = Dataset("DYJetsToLL_M-50_ext", 6024.  , ROOT.kBlue-7)
#DYjets = Dataset("DYJetsToLL_M-50-amcatnloFXFX_ext", 1921.8*3.  , ROOT.kBlue-7) \
    #+Dataset("DYJetsToLL_M-10to50-amcatnloFXFX_ext",18610.,ROOT.kBlue-7)
DYjets = Dataset("DYJetsToLL_M-50-amcatnloFXFX_ext", 6024.  , ROOT.kGreen+3)
DYjets.label = "DY+jets"

DYjetsNLO = Dataset("DYJetsToLL_M-50-amcatnloFXFX_ext", 1921.8*3.  , ROOT.kGreen+3)#\
    #+Dataset("DYJetsToLL_M-10to50-amcatnloFXFX_ext",18610.,ROOT.kBlue-7)
#DYjetsNLO.label = "DY+jets (NLO)"
#DYjetsNLO.label = "DY+jets"
DYjetsNLO.label = "Drell-Yan/Z"
DYjetsLO = Dataset("DYJetsToLL_M-50-madgraphMLM_ext", 1921.8*3.  , ROOT.kGreen+3)#\
    #+Dataset("DYJetsToLL_M-10to50-madgraphMLM",18610.,ROOT.kBlue-7)
#DYjetsLO = Dataset("DYJetsToLL_M-50-madgraphMLM_ext", 6024.  , ROOT.kBlue-7)
DYjetsLO.label = "DY+jets (LO)"

###############################################################################
## ttbar
tt = Dataset("TTTo2L2Nu_TuneCUETP8M2_ttHtranche3", 87.31  , ROOT.kOrange+8)
tt.label = "t#bar{t}"
###############################################################################
## WJets
#wjets = Dataset("WJetsToLNu-madgraphMLM_ext", 61526.7  , ROOT.kGreen+3)#LO
#wjets = Dataset("WJetsToLNu-amcatnloFXFX_ext", 61526.7  , ROOT.kGreen+3)#NLO
wjets = Dataset("WJetsToLNu-amcatnloFXFX_Total", 61526.7  , ROOT.kBlue-9)#NLO
wjets.label = "W+jets"
###############################################################################
## WGamma
wgamma = Dataset("WGToLNuG-amcatnloFXFX_ext", 489.  , ROOT.kRed+3)#NLO
wgamma.label = "W+#gamma"

#
singletop = Dataset("ST_s-channel_4f_leptonDecays-amcatnlo",3.36,ROOT.kOrange+4)\
    +Dataset("ST_t-channel_antitop_4f_inclusiveDecaysV2",80.95,ROOT.kOrange+4)\
    +Dataset("ST_t-channel_top_4f_inclusiveDecaysV2",136.02,ROOT.kOrange+4)\
    +Dataset("ST_tW_antitop_5f_NoFullyHadronicDecays_ext",11.7,ROOT.kOrange+4)\
    +Dataset("ST_tW_top_5f_NoFullyHadronicDecays_ext",11.7,ROOT.kOrange+4)\
    +Dataset("ST_tWll_5f_LO-MadGraph",0.01123,ROOT.kOrange+4)
singletop.label= "single t"



#totalMC= Dataset("WWG-amcatnlo_ext", 0.2147, ROOT.kBlack)\
        #+Dataset("ZZTo2L2Nu_powheg_ext1", 0.5644, ROOT.kBlack)\
        #+Dataset("ZGTo2LG_ext", 117.864 , ROOT.kBlack)\
        #+Dataset("WZG-amcatnlo", 0.04123 , ROOT.kBlack)\
        #+Dataset("TTGamma_Dilept-amcatnlo", 0.6352 , ROOT.kBlack)\
        #+Dataset("DYJetsToLL_M-50_ext", 1921.8*3.  , ROOT.kBlack)\
        #+Dataset("TTTo2L2Nu_TuneCUETP8M2_ttHtranche3", 87.31  , ROOT.kOrange+8)\
totalMC= Dataset("WWG-amcatnlo_ext", 0.2147, ROOT.kBlack)\
        +Dataset("ZZTo2L2Nu_powheg_ext1", 0.5644, ROOT.kBlack)\
        +Dataset("ZGTo2LG_ext", 117.864 , ROOT.kBlack)\
        +Dataset("WZG-amcatnlo", 0.04123 , ROOT.kBlack)\
        +Dataset("TTGamma_Dilept-amcatnlo", 0.6352 , ROOT.kBlack)\
        +Dataset("DYJetsToLL_M-50_ext", 1921.8*3.  , ROOT.kBlack)\
        +Dataset("TTTo2L2Nu_TuneCUETP8M2_ttHtranche3", 87.31  , ROOT.kOrange+8)\


###############################################################################
# ##Signal samples
t5wg_1600_100 = Dataset("SMS-T5Wg_1600_100", 0.00810078, ROOT.kRed, "")
t5wg_1600_100.label = "T5Wg 1600 100"
t5wg_1600_1500 = Dataset("SMS-T5Wg_1600_1500", 0.00810078, ROOT.kRed+4, "")
t5wg_1600_1500.label = "T5Wg 1600 1500"
t5wg_2000_100 = Dataset("SMS-T5Wg_2000_100", 0.000981077, ROOT.kRed+4, "")
t5wg_2000_100.label = "T5Wg 2000 100"
t5wg_1750_1700 = Dataset("SMS-T5Wg_1750_1700", 0.00359842, ROOT.kRed+4, "")
t5wg_1750_1700 = Dataset("SMS-T5Wg_1750_1700", 0.00359842, ROOT.kRed+4, "")
t5wg_1600_800 = Dataset("SMS-T5Wg_1600_800", 0.000981077, ROOT.kRed+4, "")

t6gg_1750_1650 = Dataset("SMS-T6gg_1750_1650", 0.000646271, ROOT.kRed+4, "")
t6gg_1750_1650.label = "T6gg 1750 1650"
t6gg_1300_600 = Dataset("SMS-T6gg_1300_600", 0.0086557, ROOT.kRed+4, "")
t6gg_1100_600 = Dataset("SMS-T6gg_1100_600", 0.0313372, ROOT.kRed+4, "")

tchiwg_700 = Dataset("SMS-TChiWG_700", 9.51032/1000, ROOT.kRed+4, "")


tching_1200 = Dataset("SMS-TChiNG_1200",(0.196044+0.415851)/1000.,ROOT.kRed-2,"")
#tching_400 = Dataset("SMS-TChiNG_400",58.6311/1000. *4.,ROOT.kRed-3,"")
tching_400 = Dataset("SMS-TChiNG_400",(58.6311+121.013)/1000.,ROOT.kRed-3,"")
tching_600 = Dataset("SMS-TChiNG_600",(9.49913+20.1372)/1000.,ROOT.kRed-4,"")

#t5bbbbzg_1800_1700 = Dataset("SMS-T5bbbbZg_1800_1700",0.00276133*2.,ROOT.kBlue-4,"")
#t5bbbbzg_1800_400 = Dataset("SMS-T5bbbbZg_1800_400",0.00276133*2.,ROOT.kBlue+3,"")
#t5bbbbzg_1800_600 = Dataset("SMS-T5bbbbZg_1800_600",0.00276133*2.,ROOT.kGreen+1,"")
t5bbbbzg_1800_1700 = Dataset("SMS-T5bbbbZg_1800_1700",0.00276133,ROOT.kBlue-4,"")
t5bbbbzg_1800_400 = Dataset("SMS-T5bbbbZg_1800_400",0.00276133,ROOT.kBlue+3,"")
t5bbbbzg_1800_600 = Dataset("SMS-T5bbbbZg_1800_600",0.00276133,ROOT.kGreen+1,"")
t5bbbbzg_1500_1400 = Dataset("SMS-T5bbbbZg_1500_1400",0.0141903,ROOT.kBlue-4,"")
t5bbbbzg_1500_400 = Dataset("SMS-T5bbbbZg_1500_400",0.0141903,ROOT.kBlue+3,"")
t5bbbbzg_1500_600 = Dataset("SMS-T5bbbbZg_1500_600",0.0141903,ROOT.kGreen+1,"")
t5ttttzg_1800_400 = Dataset("SMS-T5ttttZg_1800_400",0.00276133,ROOT.kYellow-1,"")
t5ttttzg_1800_600 = Dataset("SMS-T5ttttZg_1800_600",0.00276133,ROOT.kYellow+1,"")
t6ttZg_600_300 = Dataset("SMS-T6ttZg_600_300",0.174599,ROOT.kYellow-1,"")
t6ttZg_900_500 = Dataset("SMS-T6ttZg_900_500",0.0128895,ROOT.kYellow+1,"")
gmsb_240_230 = Dataset("SMS-GMSB_240_230",645.0941/1000.,ROOT.kYellow+1,"")
gmsb_290_205 = Dataset("SMS-GMSB_240_230",1350.5800/1000.,ROOT.kYellow-1,"")

#ggm_m1550_m2750 = Dataset("SMS-GGM_M1550_M2750",0.0473802983761, ROOT.kOrange,"")
#ggm_m11400_m31000 = Dataset("SMS-GGM_M11400_M31000",0.00150699994992, ROOT.kYellow,"")
#ggm_m11200_m21000 = Dataset("SMS-GGM_M11200_M21000",0.00939011014998, ROOT.kYellow,"")
#ggm_m1800_m2600 = Dataset("SMS-GGM_M11200_M21000",0.166700005531, ROOT.kGreen+2,"")


import collections
class SampleCollection(collections.MutableMapping):
    """Dictionary used to store signal datasets.
    The first time a dataset is requested, it is created"""

    colors = [ROOT.kGreen+4, ROOT.kGreen+1] + [ ROOT.kGreen-i for i in range(5) ] + range(1000)

    def __init__(self, *args, **kwargs):
        self.store = dict()
        self.update(dict(*args, **kwargs))  # use the free update to set keys

    def __getitem__(self, key):
        if key not in self.store:
            self.store[key] = Dataset( key, col=self.colors[len(self)] )
            self.store[key].label = aux.getSignalLabel( key )
        return self.store[key]

    def __setitem__(self, key, value):
        self.store[key] = value

    def __delitem__(self, key):
        del self.store[key]

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)

signal = SampleCollection()

if __name__ == "__main__":
    # print information
    print gjets.getLatexTableHeader()
    for i in gjets_dr, qcd, ttjets_ht, ttg, wjets, wg, znunu, zg: print i.getLatexTableLine(),
    print "\\end{tabular}"
