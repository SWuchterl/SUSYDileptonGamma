import ROOT
import copy
import auxiliary as aux
import os.path
import style


path = "../minimal/output/"
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
dataDoubleSF.label = "Data (ll_{SF}})"


###############################################################################
# Simulation
###############################################################################

###############################################################################
## WWGamma
wwgamma = Dataset("WWG-amcatnlo_ext", 0.2147, ROOT.kCyan-7)
wwgamma.label = "WW#gamma"

###############################################################################
## ZZ
#zz = Dataset("ZZTo2L2Nu_powheg", 0.5644, ROOT.kYellow) \
zz =    Dataset("ZZTo2L2Nu_powheg_ext1", 0.5644, ROOT.kYellow)
zz.label = "ZZ"

###############################################################################
## ZGamma
#zgamma = Dataset("ZGTo2LG", 117.864 , ROOT.kGreen)\
    #+Dataset("ZGTo2LG_ext", 117.864 , ROOT.kGreen) \
    #+Dataset("ZGTo2LG_PtG-130", 0.1404 , ROOT.kGreen)
zgamma = Dataset("ZGTo2LG_ext", 117.864 , ROOT.kGreen-3) #\
    #+Dataset("ZGTo2LG_PtG-130", 0.1404 , ROOT.kGreen-3) #\
#zgamma = Dataset("ZGTo2LG", 117.864 , ROOT.kGreen) #\
zgamma.label = "Z#gamma"

###############################################################################
## WZGamma
wzgamma = Dataset("WZG-amcatnlo", 0.04123 , ROOT.kMagenta+2)
wzgamma.label = "WZ#gamma"

###############################################################################
## TTGamma_dilept
## LO cross section calculated "by hand" from sample
ttgamma = Dataset("TTGamma_Dilept-amcatnlo", 0.6352 , ROOT.kRed+1)
ttgamma.label = "tt#gamma"

###############################################################################
## DY+jets
DYjets = Dataset("DYJetsToLL_M-50_ext", 1921.8*3.  , ROOT.kBlue-7)
DYjets.label = "DY+jets"

###############################################################################
## ttbar
tt = Dataset("TTTo2L2Nu_TuneCUETP8M2_ttHtranche3", 87.31  , ROOT.kOrange+8)
tt.label = "t#bar{t}"
###############################################################################
## WJets
wjets = Dataset("WJetsToLNu-madgraphMLM_ext", 61526.7  , ROOT.kGreen+3)#LO
#wjets = Dataset("WJetsToLNu-madgraphMLM_ext_hists", 61526.7  , ROOT.kGreen+3)#NLO
wjets.label = "W+jets"


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


tching_900 = Dataset("SMS-TChiNG_900",)

t5bbbbzg_1800_1700 = Dataset("SMS-T5bbbbZg_1800_1700",0.00276133,ROOT.kRed+4,"")
t5bbbbzg_1800_400 = Dataset("SMS-T5bbbbZg_1800_400",0.00276133,ROOT.kBlue+2,"")


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
