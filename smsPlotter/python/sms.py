from array import *

#commonZmin = 4e-4
#commonZmax = 2e-2
#commonZmin = 1e-2

# commonZmin = 5e-3  # used
commonZmin = 1.8e-2  # used
# commonZmin = 5e-10  # used
commonZmax = 4.5e-2

squarkZmin = commonZmin
squarkZmax = commonZmax
gluinoZmin = commonZmin
gluinoZmax = commonZmax

gmsbZmin = 5e-3
gmsbZmax = 5e-1


class sms():

    def __init__(self, modelname):
        if modelname.find("T1tttt") != -1:
            self.T1tttt()
        if modelname.find("T5ttttDM175") != -1:
            self.T5ttttDM175()
        if modelname.find("T1bbbb") != -1:
            self.T1bbbb()
        if modelname.find("T1qqqq") != -1:
            self.T1qqqq()
        if modelname.find("T5gg") != -1:
            self.T5gg()
        if modelname.find("T5Wg") != -1:
            self.T5Wg()
        if modelname.find("T5bbbbZg") != -1:
            self.T5Zg()
        if modelname.find("T6ttZg") != -1:
            self.T6Zg()
        if modelname.find("T6gg") != -1:
            self.T6gg()
        if modelname.find("T6Wg") != -1:
            self.T6Wg()
        if modelname.find("GGM") != -1:
            self.GGM()
        if modelname.find("GMSB") != -1:
            self.GMSB()
        #if modelname.find("GMSB") != -1: self.GGM()
        if modelname.find("M2") != -1:
            self.GGMM1M2()
        if modelname.find("M3") != -1:
            self.GGMM1M3()
        if modelname.find("TChiNG") != -1:
            self.TChiNG_BR()

    def T6gg(self):
        # model name
        self.modelname = "T6gg"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label = "pp #rightarrow #tilde{q}#tilde{q}, #tilde{q} #rightarrow q%s, %s #rightarrow #gamma#tilde{G}" % (
            lsp_s, lsp_s)
        self.label2 = ""
        # scan range to plot
        self.Xmin = 1100.
        self.Xmax = 2050.
        self.Ymin = 0.
        self.Ymax = 2700.
        self.Zmin = squarkZmin
        self.Zmax = squarkZmax
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{q}}}} (GeV)"
        # LSP
        self.LSP = "m#kern[0.1]{_{" + lsp_s + "}} (GeV)"
        # diagonal lines
        self.diagOn = True

    def T6Wg(self):
        # model name
        self.modelname = "T6Wg"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label = "pp #rightarrow #tilde{q}#tilde{q}, #tilde{q} #rightarrow q%s, %s #rightarrow #gamma/W^{#pm}#tilde{G}" % (
            lsp_s, lsp_s)
        self.label2 = ""
        # scan range to plot
        self.Xmin = 1100.
        self.Xmax = 2050.
        self.Ymin = 0.
        self.Ymax = 2700.
        self.Zmin = squarkZmin
        self.Zmax = squarkZmax
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{q}}}} (GeV)"
        # LSP
        self.LSP = "m#kern[0.1]{_{" + lsp_s + "}} (GeV)"
        # diagonal lines
        self.diagOn = True

    def T5gg(self):
        # model name
        self.modelname = "T5gg"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label = "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow qq%s, %s #rightarrow #gamma#tilde{G}" % (
            lsp_s, lsp_s)
        self.label2 = ""
        # scan range to plot
        self.Xmin = 1400.
        self.Xmax = 2500.
        self.Ymin = 0.
        self.Ymax = 3300.
        self.Zmin = gluinoZmin
        self.Zmax = gluinoZmax
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV)"
        # LSP
        self.LSP = "m#kern[0.1]{_{" + lsp_s + "}} (GeV)"
        # diagonal lines
        self.diagOn = True

    def T5Wg(self):
        # model name
        self.modelname = "T5Wg"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label = "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow qq%s, %s #rightarrow #gamma/W^{#pm}#tilde{G}" % (
            lsp_s, lsp_s)
        self.label2 = ""
        # scan range to plot
        self.Xmin = 1400.
        self.Xmax = 2500.
        self.Ymin = 0.
        self.Ymax = 3300.
        self.Zmin = gluinoZmin
        self.Zmax = gluinoZmax
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV)"
        # LSP
        self.LSP = "m#kern[0.1]{_{" + lsp_s + "}} (GeV)"
        # diagonal lines
        self.diagOn = True

    def T5Zg(self):
        # model name
        self.modelname = "T5bbbbZg"
        # decay chain
        # lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label = "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow bb%s, %s #rightarrow #gamma/Z#tilde{G}" % (
            lsp_s, lsp_s)
        self.label2 = ""
        # scan range to plot
        #self.Xmin = 800.
        # self.Xmin = 800.
        # self.Xmin = 0.
        self.Xmin = 1000.  # yep
        self.Xmax = 2500.
        # self.Ymin = 0.
        self.Ymin = 100.
        self.Ymax = 3000.
        #self.Zmin = gluinoZmin
        self.Zmin = gluinoZmin
        self.Zmax = gluinoZmax
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} (GeV)"
        # LSP
        # self.LSP = "m#kern[0.1]{_{" + lsp_s + "}} (GeV)"
        self.LSP = "m#kern[0.1]{_{NLSP}} (GeV)"
        # diagonal lines
        self.diagOn = True

    def T6Zg(self):
        # model name
        self.modelname = "T5ttZg"
        # decay chain
        # lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label = "pp #rightarrow #tilde{t}#tilde{t}, #tilde{t} #rightarrow t%s, %s #rightarrow #gamma/Z#tilde{G}" % (
            lsp_s, lsp_s)
        self.label2 = ""
        # scan range to plot
        #self.Xmin = 800.
        # self.Xmin = 800.
        self.Xmin = 200.
        self.Xmax = 1500.
        #self.Ymin = 0.
        self.Ymin = 0.
        self.Ymax = 1500.
        #self.Zmin = gluinoZmin
        self.Zmin = gluinoZmin
        self.Zmax = gluinoZmax
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{t}}}} (GeV)"
        # LSP
        self.LSP = "m#kern[0.1]{_{" + lsp_s + "}} (GeV)"
        # diagonal lines
        self.diagOn = True

    def GGM(self):
        # model name
        self.modelname = "GGM"
        # decay chain
        self.label = "GGM"
        self.label2 = ""
        # scan range to plot
        self.Xmin = 205 - 12.5
        self.Xmax = 840 + 12.5
        self.Ymin = 390 - 12.5
        # self.Ymax = 1015+12.5
        self.Ymax = 900
        self.Zmin = 0.001
        self.Zmax = 2.
        # produce sparticle
        self.sParticle = "m_{#tilde{B}} [GeV]"
        # LSP
        self.LSP = "m_{#tilde{W}} [GeV]"
        # diagonal lines
        self.diagOn = True

    def GGMM1M2(self):
        # model name
        self.modelname = "GGM_M1_M2"
        # decay chain
        self.label = "GGM M1M2"
        self.label2 = ""
        # scan range to plot
        self.Xmin = 200 - 10
        self.Xmax = 1500 + 10
        self.Ymin = 200 - 10
        # self.Ymax = 1015+12.5
        self.Ymax = 1900 + 10
        self.Zmin = 0.001
        self.Zmax = 2.
        # produce sparticle
        self.sParticle = "m_{#tilde{B}} [GeV]"
        # LSP
        self.LSP = "m_{#tilde{W}} [GeV]"
        # diagonal lines
        self.diagOn = False

    def GMSB(self):
        # model name
        self.modelname = "GMSB"
        # decay chain
        self.label = "GMSB electroweak production"
        self.label2 = ""
        # scan range to plot
        # self.Xmin = 200 - 10
        # self.Xmax = 1015 + 10
        # self.Xmin = 205 - 12.5
        self.Xmin = 200
        self.Xmax = 840 + 12.5
        self.Ymin = 200 - 10
        # self.Ymax = 1015+12.5
        # self.Ymax = 900
        self.Ymax = 1000
        # self.Ymax = 1300 + 10
        # self.Zmin = 0.001
        # self.Zmax = 2.
        self.Zmin = gmsbZmin
        self.Zmax = gmsbZmax
        # produce sparticle
        self.sParticle = "m_{#tilde{B}} (GeV)"
        # LSP
        self.LSP = "m_{#tilde{W}} (GeV)"
        # diagonal lines
        self.diagOn = True

    def GGMM1M3(self):
        # model name
        self.modelname = "GGM_M1_M3"
        # decay chain
        self.label = "GGM M1M3"
        self.label2 = ""
        # scan range to plot
        self.Xmin = 50 - 10
        self.Xmax = 1500 + 10
        self.Ymin = 1000 - 10
        # self.Ymax = 1015+12.5
        self.Ymax = 2500 + 10
        self.Zmin = 0.001
        self.Zmax = 2.
        # produce sparticle
        self.sParticle = "m_{#tilde{B}} [GeV]"
        # LSP
        self.LSP = "m_{#tilde{W}} [GeV]"
        # diagonal lines
        self.diagOn = False

    def TChiNG_BR(self):
        print "working"
        # model name
        self.modelname = "TChiNG_BR_scaled"
        # decay chain
        self.label = "TChiNG_BR_scaled"
        self.label2 = ""
        # scan range to plot
        self.Xmin = 300 - 10
        self.Xmax = 1300 + 10
        self.Ymin = 0
        # self.Ymax = 1015+12.5
        self.Ymax = 100 + 50
        self.Zmin = 0.001
        self.Zmax = 2.
        # produce sparticle
        self.sParticle = "m_{NLSP} [GeV]"
        # LSP
        self.LSP = "BR(#rightarrow #gamma)"
        # diagonal lines
        self.diagOn = False

    def T1tttt(self):
        # model name
        self.modelname = "T1tttt"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label = "pp #rightarrow #tilde{g} #tilde{g}, #tilde{g} #rightarrow t #bar{t} " + lsp_s
        self.label2 = ""
        # scan range to plot
        self.Xmin = 600.
        self.Xmax = 2000.
        self.Ymin = 0.
        self.Ymax = 1800.
        self.Zmin = 0.001
        self.Zmax = 2.
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} [GeV]"
        # LSP
        self.LSP = "m#kern[0.1]{_{" + lsp_s + "}} [GeV]"
        # turn off diagonal lines
        self.diagOn = False

    def T5ttttDM175(self):
        # model name
        self.modelname = "T5ttttDM175"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label = "pp #rightarrow #tilde{g} #tilde{g},  #tilde{g} #rightarrow #tilde{t}_{1} t,  #tilde{t}_{1} #rightarrow #bar{t} " + lsp_s
        self.label2 = "m_{#tilde{t}_{1}} - m_{#tilde{#chi}^{0}_{1}} = 175 GeV"
        # scan range to plot
        self.Xmin = 600.
        self.Xmax = 1700.
        self.Ymin = 0.
        self.Ymax = 1800.
        self.Zmin = 0.001
        self.Zmax = 2.
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} [GeV]"
        # LSP
        self.LSP = "m#kern[0.1]{_{" + lsp_s + "}} [GeV]"
        # turn off diagonal lines
        self.diagOn = False

    def T1bbbb(self):
        # model name
        self.modelname = "T1bbbb"
        # decay chain
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        self.label = "pp #rightarrow #tilde{g} #tilde{g}, #tilde{g} #rightarrow b #bar{b} " + lsp_s
        self.label2 = ""
        # plot boundary. The top 1/4 of the y axis is taken by the legend
        self.Xmin = 600.
        self.Xmax = 1950.
        self.Ymin = 0.
        self.Ymax = 1800.
        self.Zmin = 0.001
        self.Zmax = 2.
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} [GeV]"
        # LSP
        self.LSP = "m#kern[0.1]{_{" + lsp_s + "}} [GeV]"
        # turn off diagonal lines
        self.diagOn = False

    def T1qqqq(self):
        # model name
        self.modelname = "T1qqqq"
        # decay chain
        self.label = "pp #rightarrow #tilde{g} #tilde{g}, #tilde{g} #rightarrow q #bar{q} #tilde{#chi}^{0}_{1}"
        self.label2 = ""
        # plot boundary. The top 1/4 of the y axis is taken by the legend
        self.Xmin = 600.
        self.Xmax = 1950.
        self.Ymin = 0.
        self.Ymax = 1600.
        self.Zmin = 0.001
        self.Zmax = 2.
        # produce sparticle
        self.sParticle = "m#kern[0.1]{_{#lower[-0.12]{#tilde{g}}}} [GeV]"
        # LSP
        self.LSP = "m#kern[0.1]{_{" + lsp_s + "}} [GeV]"
        # turn off diagonal lines
        self.diagOn = False
