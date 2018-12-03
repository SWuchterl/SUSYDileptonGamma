#from dataMC import labels,frange

import ROOT
from ROOT import *
from array import array
from include import *
import numpy as np
import pickle as pkl
#import CR_tt,CR_DY,CR_WZ,CR_ZZ

import style
style.defaultStyle()

pklZZ = pkl.load(open("plots_CR_zz/factors/CRZZ.pkl", "rb"))
pklDY = pkl.load(open("plots_CR_dy/factors/CRDY.pkl", "rb"))
pklTT = pkl.load(open("plots_CR_tt/factors/CRTT.pkl", "rb"))
#pklTTsmall = pkl.load( open( "plots_CR_tt/factors/CRTT_small.pkl", "rb" ) )
#pklTTbig = pkl.load( open( "plots_CR_tt/factors/CRTT_big.pkl", "rb" ) )
pklWZ = pkl.load(open("plots_CR_wz/factors/CRWZ.pkl", "rb"))
pklZZchi = pkl.load(open("plots_CR/chi/ZZ_chi.pkl", "rb"))
pklZZchiBig = pkl.load(open("plots_CR/chi/ZZ_chi_big.pkl", "rb"))
pklZZchiSmall = pkl.load(open("plots_CR/chi/ZZ_chi_small.pkl", "rb"))
#pklDYEEchi = pkl.load( open( "plots_CR/chi/DY_chi_EE.pkl", "rb" ) )
#pklDYMMchi = pkl.load( open( "plots_CR/chi/DY_chi_MM.pkl", "rb" ) )
#pklDYLLchi = pkl.load( open( "plots_CR/chi/DY_chi_LL.pkl", "rb" ) )
pklTTchi = pkl.load(open("plots_CR/chi/TT_chi.pkl", "rb"))
pklTTchiBig = pkl.load(open("plots_CR/chi/TT_chi_big.pkl", "rb"))
pklTTchiSmall = pkl.load(open("plots_CR/chi/TT_chi_small.pkl", "rb"))
pklWZchi = pkl.load(open("plots_CR/chi/WZ_chi.pkl", "rb"))
pklWZchiBig = pkl.load(open("plots_CR/chi/WZ_chi_big.pkl", "rb"))
pklWZchiSmall = pkl.load(open("plots_CR/chi/WZ_chi_small.pkl", "rb"))
pklDYchiLL = pkl.load(open("plots_CR/chi/DY_chi_LL.pkl", "rb"))
pklDYchiLLBig = pkl.load(open("plots_CR/chi/DY_chi_LL_big.pkl", "rb"))
pklDYchiLLSmall = pkl.load(open("plots_CR/chi/DY_chi_LL_small.pkl", "rb"))
#pklDYchiEE = pkl.load( open( "plots_CR/chi/DY_chi_EE.pkl", "rb" ) )
#pklDYchiEEBig = pkl.load( open( "plots_CR/chi/DY_chi_EE_big.pkl", "rb" ) )
#pklDYchiEESmall = pkl.load( open( "plots_CR/chi/DY_chi_EE_small.pkl", "rb" ) )
#pklDYchiMM = pkl.load( open( "plots_CR/chi/DY_chi_MM.pkl", "rb" ) )
#pklDYchiMMBig = pkl.load( open( "plots_CR/chi/DY_chi_MM_big.pkl", "rb" ) )
#pklDYchiMMSmall = pkl.load( open( "plots_CR/chi/DY_chi_MM_small.pkl", "rb" ) )

nameDict = {
    "phi1": "#phi_{leading}",
    "phi2": "#phi_{trailing}",
    "eta1": "#eta_{leading}",
    "eta2": "#eta_{trailing}",
    "pt1": "p_{T}^{leading}",
    "pt2": "p_{T}^{trailing}",
    "pt_g1": "p_{T}^{#gamma}",
    "met": "p_{T}^{miss}",
    "m_ll": "m_{ll}",
    "m_ll2": "m_{l_{3}l_{4}}",
    "ht": "H_{T}",
    "n_jets": "N_{Jets}",
}


def drawZZ():
    valueZZInt = pklZZ["LL"]["eta1"][0]
    errZZInt = pklZZ["LL"]["eta1"][1]
    arZZChi = []
    arZZChiUp = []
    arZZChiDn = []
    arZZChiName = []

    for comb in pklZZchiBig:
        for variable in pklZZchiBig[comb]:
            arZZChi.append(pklZZchiBig[comb][variable]["value"])
            arZZChiUp.append(pklZZchiBig[comb][variable]["erUp"])
            arZZChiDn.append(pklZZchiBig[comb][variable]["erDown"])
            arZZChiName.append(variable)
    for comb in pklZZchi:
        for variable in pklZZchi[comb]:
            arZZChi.append(pklZZchi[comb][variable]["value"])
            arZZChiUp.append(pklZZchi[comb][variable]["erUp"])
            arZZChiDn.append(pklZZchi[comb][variable]["erDown"])
            arZZChiName.append(variable)
    for comb in pklZZchiSmall:
        for variable in pklZZchiSmall[comb]:
            arZZChi.append(pklZZchiSmall[comb][variable]["value"])
            arZZChiUp.append(pklZZchiSmall[comb][variable]["erUp"])
            arZZChiDn.append(pklZZchiSmall[comb][variable]["erDown"])
            arZZChiName.append(variable)

    # print arZZChi
    # print arZZChiUp
    # print arZZChiDn
    # print arZZChiName

    c = TCanvas("canvas", "", 800, 1000)

    x = [valueZZInt]
    ex = [errZZInt]
    y = [len(arZZChi) / 2.]
    eyU = [round(len(arZZChi) / 2.) + 1]
    eyD = [len(arZZChi) / 2.]
    ge = TGraphAsymmErrors()
    for i in range(len(x)):
        ge.SetPoint(i, x[i], y[i])
        ge.SetPointError(i, ex[i], ex[i], eyD[i], eyU[i])
    # ge.SetFillStyle(3003)
    # ge.SetFillStyle(3004)
    ge.SetLineColor(ROOT.kWhite)
    # ge.SetFillColor(ROOT.kBlue)
    ge.SetFillColorAlpha(ROOT.kBlue, 0.3)
    ge.Draw("same a2")

    gr = TGraphAsymmErrors()

    for i in range(len(arZZChi)):
        y = float(i + 1)
        gr.SetPoint(i, arZZChi[i], y)
        gr.SetPointError(i, arZZChiDn[i], arZZChiUp[i], 0., 0.)

    gr.SetTitle("; Scale Factor #alpha_{ZZ}; ")
    ge.SetTitle("; Scale Factor #alpha_{ZZ}; ")
    # gr.SetMarkerColor(4)
    gr.SetMarkerColor(kRed)
    gr.SetMarkerStyle(21)
    gr.Draw("sameP")

    line = TLine(valueZZInt, 0, valueZZInt, len(arZZChi) + 1)
    line.SetLineColor(kBlue)
    line.SetLineWidth(2)
    line.Draw("SAME")

    mean = np.mean(arZZChi)
    lineMEAN = TLine(mean, 0, mean, len(arZZChi) + 1)
    lineMEAN.SetLineColor(kRed)
    lineMEAN.SetLineWidth(2)
    lineMEAN.Draw("SAME")

    ge.GetXaxis().SetLimits(0.6, 1.5)
    ge.GetYaxis().SetLabelOffset(1)
    ge.GetYaxis().SetNdivisions(len(arZZChi) * 2)

    #t = TText()
    t = TLatex()
    t.SetTextAlign(32)
    # t.SetTextAlign(42)
    # t.SetTextAngle(30)
    # t.SetNDC()
    t.SetTextSize(0.02)
    # t.SetTextFont(72)
    for i in range(len(arZZChi)):
        t.DrawLatex(0.58, (i + 1), nameDict[arZZChiName[i]])
        # t.DrawText(0.8,i+1,arTTChiName[i])

    smallT = TText()
    smallT.SetTextAngle(90)
    smallT.SetTextSize(0.025)
    smallT.SetTextFont(42)
    # smallT.SetNDC()
    smallT.DrawText(0.48, 3. * len(arZZChi) / 3. - 6, "double binning")
    smallT.DrawText(0.48, 2 * len(arZZChi) / 3. - 6, "normal binning")
    smallT.DrawText(0.48, len(arZZChi) / 3. - 6, "half binning")

    cutLine = TLine(0.6, 2 * len(arZZChi) / 3. + 0.5,
                    1.5, 2 * len(arZZChi) / 3. + 0.5)
    cutLine.SetLineColor(kBlack)
    cutLine.SetLineWidth(2)
    cutLine.Draw("SAME")

    cutLine2 = TLine(0.6, len(arZZChi) / 3. + 0.5,
                     1.5, len(arZZChi) / 3. + 0.5)
    cutLine2.SetLineColor(kBlack)
    cutLine2.SetLineWidth(2)
    cutLine2.Draw("SAME")

    l = ROOT.TLatex(
        # 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
        0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Private Work}}")
    l.SetNDC()
    # l.Draw()

    lum = ROOT.TLatex(.62, .95, "%.1f fb^{-1} (%s TeV)" %
                      (aux.intLumi / 1000., aux.Label.cmsEnergy))
    lum.SetNDC()
    lum.Draw()

    l4 = ROOT.TLatex(0.3, .89, "#scale[0.76]{#font[52]{ZZ Control Region}}")
    l4.SetNDC()
    l4.Draw()

    leg = TLegend(0.18, 0.45, 0.49, 0.55)
    leg.AddEntry(line, "#alpha (int. method)", "l")
    leg.AddEntry(ge, "stat. error (int. method)", "f")
    leg.AddEntry(gr, "#alpha (#chi^{2} method)", "lep")
    leg.AddEntry(lineMEAN, "mean (#chi^{2} method)", "l")
    leg.SetTextSize(0.03)
    # leg.SetFillStyle(0)
    leg.Draw("same")

    c.SetGridy()
    c.Update()

    c.SaveAs('plots_CR/chi/ZZ_Compare.pdf')


def drawDYLL():
    valueDYInt = pklDY["LL"]["eta1"][0]
    errDYInt = pklDY["LL"]["eta1"][1]
    arDYChi = []
    arDYChiUp = []
    arDYChiDn = []
    arDYChiName = []

    for comb in pklDYchiLLBig:
        for variable in pklDYchiLLBig[comb]:
            arDYChi.append(pklDYchiLLBig[comb][variable]["value"])
            arDYChiUp.append(pklDYchiLLBig[comb][variable]["erUp"])
            arDYChiDn.append(pklDYchiLLBig[comb][variable]["erDown"])
            arDYChiName.append(variable)
    for comb in pklDYchiLL:
        for variable in pklDYchiLL[comb]:
            print pklDYchiLL[comb][variable]
            arDYChi.append(pklDYchiLL[comb][variable]["value"])
            arDYChiUp.append(pklDYchiLL[comb][variable]["erUp"])
            arDYChiDn.append(pklDYchiLL[comb][variable]["erDown"])
            arDYChiName.append(variable)
    for comb in pklDYchiLLSmall:
        for variable in pklDYchiLLSmall[comb]:
            arDYChi.append(pklDYchiLLSmall[comb][variable]["value"])
            arDYChiUp.append(pklDYchiLLSmall[comb][variable]["erUp"])
            arDYChiDn.append(pklDYchiLLSmall[comb][variable]["erDown"])
            arDYChiName.append(variable)

    #print arDYChi
    # print arDYChiUp
    # print arDYChiDn
    # print arDYChiName

    c = TCanvas("canvas", "", 800, 1000)

    x = [valueDYInt]
    ex = [errDYInt]
    y = [len(arDYChi) / 2.]
    eyU = [round(len(arDYChi) / 2.) + 1]
    eyD = [len(arDYChi) / 2.]
    ge = TGraphAsymmErrors()
    for i in range(len(x)):
        ge.SetPoint(i, x[i], y[i])
        ge.SetPointError(i, ex[i], ex[i], eyD[i], eyU[i])
    # ge.SetFillStyle(3003)
    # ge.SetFillStyle(3004)
    # ge.SetFillStyle(3244)
    # ge.SetFillColor(ROOT.kBlue)
    ge.SetFillColorAlpha(ROOT.kBlue, 0.3)
    ge.SetLineColor(ROOT.kWhite)
    ge.Draw("same a2")

    gr = TGraphAsymmErrors()

    for i in range(len(arDYChi)):
        y = float(i + 1)
        gr.SetPoint(i, arDYChi[i], y)
        gr.SetPointError(i, arDYChiDn[i], arDYChiUp[i], 0., 0.)

    gr.SetTitle("; Scale Factor #alpha_{DY/Z(#gamma)}; ")
    ge.SetTitle("; Scale Factor #alpha_{DY/Z(#gamma)}; ")
    # gr.SetMarkerColor(4)
    gr.SetMarkerColor(kRed)
    gr.SetMarkerStyle(21)
    gr.Draw("sameP")

    line = TLine(valueDYInt, 0, valueDYInt, len(arDYChi) + 1)
    line.SetLineColor(kBlue)
    line.SetLineWidth(2)
    line.Draw("SAME")

    mean = np.mean(arDYChi)
    lineMEAN = TLine(mean, 0, mean, len(arDYChi) + 1)
    lineMEAN.SetLineColor(kRed)
    lineMEAN.SetLineWidth(2)
    lineMEAN.Draw("SAME")

    # ge.GetXaxis().SetLimits(1.05,1.12)
    ge.GetXaxis().SetLimits(1.02, 1.09)
    ge.GetYaxis().SetLabelOffset(1)
    ge.GetYaxis().SetNdivisions(len(arDYChi) * 2)

    #t = TText()
    t = TLatex()
    t.SetTextAlign(32)
    t.SetTextSize(0.02)
    # t.SetTextFont(72)
    for i in range(len(arDYChi)):
        t.DrawLatex(1.018, (i + 1), nameDict[arDYChiName[i]])
        # t.DrawText(0.8,i+1,arTTChiName[i])

    smallT = TText()
    smallT.SetTextAngle(90)
    smallT.SetTextSize(0.025)
    smallT.SetTextFont(42)
    smallT.DrawText(1.01, 3. * len(arDYChi) / 3. - 6, "double binning")
    smallT.DrawText(1.01, 2 * len(arDYChi) / 3. - 6, "normal binning")
    smallT.DrawText(1.01, len(arDYChi) / 3. - 6, "half binning")

    cutLine = TLine(1.02, 2 * len(arDYChi) / 3. + 0.5,
                    1.09, 2 * len(arDYChi) / 3. + 0.5)
    cutLine.SetLineColor(kBlack)
    cutLine.SetLineWidth(2)
    cutLine.Draw("SAME")

    cutLine2 = TLine(1.02, len(arDYChi) / 3. + 0.5,
                     1.09, len(arDYChi) / 3. + 0.5)
    cutLine2.SetLineColor(kBlack)
    cutLine2.SetLineWidth(2)
    cutLine2.Draw("SAME")

    l = ROOT.TLatex(
        # 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
        0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Private Work}}")
    l.SetNDC()
    # l.Draw()

    lum = ROOT.TLatex(.62, .95, "%.1f fb^{-1} (%s TeV)" %
                      (aux.intLumi / 1000., aux.Label.cmsEnergy))
    lum.SetNDC()
    lum.Draw()

    l4 = ROOT.TLatex(
        0.25, .89, "#scale[0.76]{#font[52]{DY/Z(#gamma) Control Region}}")
    l4.SetNDC()
    # l4.SetFillColor(0)
    l4.Draw()

    leg = TLegend(0.18, 0.45, 0.49, 0.55)
    leg.AddEntry(line, "#alpha (int. method)", "l")
    leg.AddEntry(ge, "stat. error (int. method)", "f")
    leg.AddEntry(gr, "#alpha (#chi^{2} method)", "lep")
    leg.AddEntry(lineMEAN, "mean (#chi^{2} method)", "l")
    leg.SetTextSize(0.03)
    # leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.Draw("same")

    c.SetGridy()
    c.Update()

    c.SaveAs('plots_CR/chi/DY_CompareLL.pdf')
# def drawDYEE():
    # valueDYInt=pklDY["EE"]["eta1"][0]
    # errDYInt=pklDY["EE"]["eta1"][1]
    # arDYChi=[]
    # arDYChiUp=[]
    # arDYChiDn=[]
    # arDYChiName=[]

    # for comb in pklDYchiEEBig:
    # for variable in pklDYchiEEBig[comb]:
    # arDYChi.append(pklDYchiEEBig[comb][variable]["value"])
    # arDYChiUp.append(pklDYchiEEBig[comb][variable]["erUp"])
    # arDYChiDn.append(pklDYchiEEBig[comb][variable]["erDown"])
    # arDYChiName.append(variable)
    # for comb in pklDYchiEE:
    # for variable in pklDYchiEE[comb]:
    # arDYChi.append(pklDYchiEE[comb][variable]["value"])
    # arDYChiUp.append(pklDYchiEE[comb][variable]["erUp"])
    # arDYChiDn.append(pklDYchiEE[comb][variable]["erDown"])
    # arDYChiName.append(variable)
    # for comb in pklDYchiEESmall:
    # for variable in pklDYchiEESmall[comb]:
    # arDYChi.append(pklDYchiEESmall[comb][variable]["value"])
    # arDYChiUp.append(pklDYchiEESmall[comb][variable]["erUp"])
    # arDYChiDn.append(pklDYchiEESmall[comb][variable]["erDown"])
    # arDYChiName.append(variable)

    #print arDYChi
    #print arDYChiUp
    #print arDYChiDn
    #print arDYChiName

    #c = TCanvas("canvas","",800,1200)

    # x=[valueDYInt]
    # ex=[errDYInt]
    # y=[len(arDYChi)/2.]
    # eyU=[round(len(arDYChi)/2.)+1]
    # eyD=[len(arDYChi)/2.]
    # ge=TGraphAsymmErrors()
    # for i in range(len(x)):
    # ge.SetPoint(i,x[i],y[i])
    # ge.SetPointError(i,ex[i],ex[i],eyD[i],eyU[i])
    # ge.SetFillStyle(3003)
    # ge.SetFillColor(ROOT.kBlue)
    #ge.Draw("same a2")

    #gr = TGraphAsymmErrors()

    # for i in range(len(arDYChi)):
    #y= float(i+1)
    # gr.SetPoint(i,arDYChi[i],y)
    # gr.SetPointError(i,arDYChiDn[i],arDYChiUp[i],0.,0.)

    # gr.SetTitle("; Scale Factor #alpha; ")
    # ge.SetTitle("; Scale Factor #alpha; ")
    # gr.SetMarkerColor(kRed)
    # gr.SetMarkerStyle(21)
    # gr.Draw("sameP")

    # line=TLine(valueDYInt,0,valueDYInt,len(arDYChi)+1)
    # line.SetLineColor(kBlue)
    # line.SetLineWidth(2)
    # line.Draw("SAME")

    # mean=np.mean(arDYChi)
    # lineMEAN=TLine(mean,0,mean,len(arDYChi)+1)
    # lineMEAN.SetLineColor(kRed)
    # lineMEAN.SetLineWidth(2)
    # lineMEAN.Draw("SAME")

    # ge.GetXaxis().SetLimits(1.02,1.12)
    # ge.GetYaxis().SetLabelOffset(1)
    # ge.GetYaxis().SetNdivisions(len(arDYChi))

    #t = TLatex()
    # t.SetTextAlign(32)
    # t.SetTextSize(0.02)
    # for i in range(len(arDYChi)):
    # t.DrawLatex(1.015,(i+1),nameDict[arDYChiName[i]])

    # smallT=TText()
    # smallT.SetTextAngle(90)
    # smallT.SetTextSize(0.025)
    # smallT.SetTextFont(42)
    #smallT.DrawText(1.005,3.*len(arDYChi)/3. -6,"NBins * 2")
    #smallT.DrawText(1.005,2*len(arDYChi)/3. -6,"normal binning")
    #smallT.DrawText(1.005,len(arDYChi)/3. -6,"NBins / 2")

    #cutLine=TLine(1.005,2*len(arDYChi)/3. +0.5,1.75,2*len(arDYChi)/3. +0.5)
    # cutLine.SetLineColor(kBlack)
    # cutLine.SetLineWidth(2)
    # cutLine.Draw("SAME")

    #cutLine2=TLine(1.005,len(arDYChi)/3. +0.5,1.75,len(arDYChi)/3. +0.5)
    # cutLine2.SetLineColor(kBlack)
    # cutLine2.SetLineWidth(2)
    # cutLine2.Draw("SAME")

    # l = ROOT.TLatex( 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
    # l.SetNDC()
    # l.Draw()

    #lum = ROOT.TLatex( .62, .95, "%.1f fb^{-1} (%s TeV)"%(aux.intLumi/1000., aux.Label.cmsEnergy) )
    # lum.SetNDC()
    # lum.Draw()

    # l4 = ROOT.TLatex( 0.17, .9, "#scale[0.66]{#font[52]{DY Control Region}}")
    # l4.SetNDC()
    # l4.Draw()

    # leg=TLegend(0.18,0.45,0.4,0.55)
    # leg.AddEntry(line,"#alpha from int. method","l")
    #leg.AddEntry(ge,"stat. error from int. method","f")
    # leg.AddEntry(gr,"#alpha from #chi^{2} method","lep")
    # leg.AddEntry(lineMEAN,"mean from #chi^{2} method","l")
    # leg.SetTextSize(0.02)
    # leg.Draw("same")

    # c.SetGridy()
    # c.Update()

    # c.SaveAs('plots_CR/chi/DY_CompareEE.pdf')
# def drawDYMM():
    # valueDYInt=pklDY["MM"]["eta1"][0]
    # errDYInt=pklDY["MM"]["eta1"][1]
    # arDYChi=[]
    # arDYChiUp=[]
    # arDYChiDn=[]
    # arDYChiName=[]

    # for comb in pklDYchiMMBig:
    # for variable in pklDYchiMMBig[comb]:
    # arDYChi.append(pklDYchiMMBig[comb][variable]["value"])
    # arDYChiUp.append(pklDYchiMMBig[comb][variable]["erUp"])
    # arDYChiDn.append(pklDYchiMMBig[comb][variable]["erDown"])
    # arDYChiName.append(variable)
    # for comb in pklDYchiMM:
    # for variable in pklDYchiMM[comb]:
    # arDYChi.append(pklDYchiMM[comb][variable]["value"])
    # arDYChiUp.append(pklDYchiMM[comb][variable]["erUp"])
    # arDYChiDn.append(pklDYchiMM[comb][variable]["erDown"])
    # arDYChiName.append(variable)
    # for comb in pklDYchiMMSmall:
    # for variable in pklDYchiMMSmall[comb]:
    # arDYChi.append(pklDYchiMMSmall[comb][variable]["value"])
    # arDYChiUp.append(pklDYchiMMSmall[comb][variable]["erUp"])
    # arDYChiDn.append(pklDYchiMMSmall[comb][variable]["erDown"])
    # arDYChiName.append(variable)

    #print arDYChi
    #print arDYChiUp
    #print arDYChiDn
    #print arDYChiName

    #c = TCanvas("canvas","",800,1200)

    # x=[valueDYInt]
    # ex=[errDYInt]
    # y=[len(arDYChi)/2.]
    # eyU=[round(len(arDYChi)/2.)+1]
    # eyD=[len(arDYChi)/2.]
    # ge=TGraphAsymmErrors()
    # for i in range(len(x)):
    # ge.SetPoint(i,x[i],y[i])
    # ge.SetPointError(i,ex[i],ex[i],eyD[i],eyU[i])
    # ge.SetFillStyle(3003)
    # ge.SetFillColor(ROOT.kBlue)
    #ge.Draw("same a2")

    #gr = TGraphAsymmErrors()

    # for i in range(len(arDYChi)):
    #y= float(i+1)
    # gr.SetPoint(i,arDYChi[i],y)
    # gr.SetPointError(i,arDYChiDn[i],arDYChiUp[i],0.,0.)

    # gr.SetTitle("; Scale Factor #alpha; ")
    # ge.SetTitle("; Scale Factor #alpha; ")
    # gr.SetMarkerColor(kRed)
    # gr.SetMarkerStyle(21)
    # gr.Draw("sameP")

    # line=TLine(valueDYInt,0,valueDYInt,len(arDYChi)+1)
    # line.SetLineColor(kBlue)
    # line.SetLineWidth(2)
    # line.Draw("SAME")

    # mean=np.mean(arDYChi)
    # lineMEAN=TLine(mean,0,mean,len(arDYChi)+1)
    # lineMEAN.SetLineColor(kRed)
    # lineMEAN.SetLineWidth(2)
    # lineMEAN.Draw("SAME")

    # ge.GetXaxis().SetLimits(1.05,1.12)
    # ge.GetYaxis().SetLabelOffset(1)
    # ge.GetYaxis().SetNdivisions(len(arDYChi))

    #t = TLatex()
    # t.SetTextAlign(32)
    # t.SetTextSize(0.02)
    # for i in range(len(arDYChi)):
    # t.DrawLatex(1.047,(i+1),nameDict[arDYChiName[i]])

    # smallT=TText()
    # smallT.SetTextAngle(90)
    # smallT.SetTextSize(0.025)
    # smallT.SetTextFont(42)
    #smallT.DrawText(1.04,3.*len(arDYChi)/3. -6,"NBins * 2")
    #smallT.DrawText(1.04,2*len(arDYChi)/3. -6,"normal binning")
    #smallT.DrawText(1.04,len(arDYChi)/3. -6,"NBins / 2")

    #cutLine=TLine(1.04,2*len(arDYChi)/3. +0.5,1.75,2*len(arDYChi)/3. +0.5)
    # cutLine.SetLineColor(kBlack)
    # cutLine.SetLineWidth(2)
    # cutLine.Draw("SAME")

    #cutLine2=TLine(1.04,len(arDYChi)/3. +0.5,1.75,len(arDYChi)/3. +0.5)
    # cutLine2.SetLineColor(kBlack)
    # cutLine2.SetLineWidth(2)
    # cutLine2.Draw("SAME")

    # l = ROOT.TLatex( 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
    # l.SetNDC()
    # l.Draw()

    #lum = ROOT.TLatex( .62, .95, "%.1f fb^{-1} (%s TeV)"%(aux.intLumi/1000., aux.Label.cmsEnergy) )
    # lum.SetNDC()
    # lum.Draw()

    # l4 = ROOT.TLatex( 0.17, .9, "#scale[0.66]{#font[52]{DY Control Region}}")
    # l4.SetNDC()
    # l4.Draw()

    # leg=TLegend(0.18,0.45,0.4,0.55)
    # leg.AddEntry(line,"#alpha from int. method","l")
    #leg.AddEntry(ge,"stat. error from int. method","f")
    # leg.AddEntry(gr,"#alpha from #chi^{2} method","lep")
    # leg.AddEntry(lineMEAN,"mean from #chi^{2} method","l")
    # leg.SetTextSize(0.02)
    # leg.Draw("same")

    # c.SetGridy()
    # c.Update()

    # c.SaveAs('plots_CR/chi/DY_CompareMM.pdf')


def drawTT():
    valueTTInt = pklTT["EM"]["eta1"][0]
    errTTInt = pklTT["EM"]["eta1"][1]
    arTTChi = []
    arTTChiUp = []
    arTTChiDn = []
    arTTChiName = []

    for comb in pklTTchiBig:
        for variable in pklTTchiBig[comb]:
            arTTChi.append(pklTTchiBig[comb][variable]["value"])
            arTTChiUp.append(pklTTchiBig[comb][variable]["erUp"])
            arTTChiDn.append(pklTTchiBig[comb][variable]["erDown"])
            arTTChiName.append(variable)
    for comb in pklTTchi:
        for variable in pklTTchi[comb]:
            arTTChi.append(pklTTchi[comb][variable]["value"])
            arTTChiUp.append(pklTTchi[comb][variable]["erUp"])
            arTTChiDn.append(pklTTchi[comb][variable]["erDown"])
            arTTChiName.append(variable)

    for comb in pklTTchiSmall:
        for variable in pklTTchiSmall[comb]:
            arTTChi.append(pklTTchiSmall[comb][variable]["value"])
            arTTChiUp.append(pklTTchiSmall[comb][variable]["erUp"])
            arTTChiDn.append(pklTTchiSmall[comb][variable]["erDown"])
            arTTChiName.append(variable)

    #print arTTChi
    #print arTTChiUp
    #print arTTChiDn
    #print arTTChiName

    c = TCanvas("canvas", "", 800, 1000)

    x = [valueTTInt]
    ex = [errTTInt]
    y = [len(arTTChi) / 2.]
    eyU = [round(len(arTTChi) / 2.) + 1]
    eyD = [len(arTTChi) / 2.]
    ge = TGraphAsymmErrors()
    for i in range(len(x)):
        ge.SetPoint(i, x[i], y[i])
        ge.SetPointError(i, ex[i], ex[i], eyD[i], eyU[i])
    # ge.SetFillStyle(3004)
    # ge.SetFillStyle(3001)
    # ge.SetFillColor(ROOT.kBlue)
    ge.SetFillColorAlpha(ROOT.kBlue, 0.3)
    ge.SetLineColor(ROOT.kWhite)
    ge.Draw("same a2")

    gr = TGraphAsymmErrors()

    for i in range(len(arTTChi)):
        y = float(i + 1)
        gr.SetPoint(i, arTTChi[i], y)
        gr.SetPointError(i, arTTChiDn[i], arTTChiUp[i], 0., 0.)

    gr.SetTitle("; Scale Factor #alpha_{t#bar{t}(#gamma)}; ")
    ge.SetTitle("; Scale Factor #alpha_{t#bar{t}(#gamma)}; ")
    # gr.SetMarkerColor(4)
    gr.SetMarkerColor(kRed)
    gr.SetMarkerStyle(21)
    gr.Draw("sameP")

    line = TLine(valueTTInt, 0, valueTTInt, len(arTTChi) + 1)
    line.SetLineColor(kBlue)
    line.SetLineWidth(2)
    line.Draw("SAME")

    mean = np.mean(arTTChi)
    lineMEAN = TLine(mean, 0, mean, len(arTTChi) + 1)
    lineMEAN.SetLineColor(kRed)
    lineMEAN.SetLineWidth(2)
    lineMEAN.Draw("SAME")

    # ge.GetXaxis().SetLimits(0.6,1.1)
    ge.GetXaxis().SetLimits(0.67, 0.88)
    ge.GetYaxis().SetLabelOffset(1)
    ge.GetYaxis().SetNdivisions(len(arTTChi) * 2)

    #t = TText()
    t = TLatex()
    t.SetTextAlign(32)
    t.SetTextSize(0.02)
    # t.SetTextFont(72)
    for i in range(len(arTTChi)):
        # t.DrawLatex(0.79,(i+1),nameDict[arTTChiName[i]])
        t.DrawLatex(0.665, (i + 1), nameDict[arTTChiName[i]])
        # t.DrawText(0.8,i+1,arTTChiName[i])

    smallT = TText()
    smallT.SetTextAngle(90)
    smallT.SetTextSize(0.025)
    smallT.SetTextFont(42)
    smallT.DrawText(0.64, 3. * len(arTTChi) / 3. - 6, "double binning")
    smallT.DrawText(0.64, 2 * len(arTTChi) / 3. - 6, "normal binning")
    smallT.DrawText(0.64, len(arTTChi) / 3. - 6, "half binning")

    cutLine = TLine(0.67, 2 * len(arTTChi) / 3. + 0.5,
                    0.88, 2 * len(arTTChi) / 3. + 0.5)
    cutLine.SetLineColor(kBlack)
    cutLine.SetLineWidth(2)
    cutLine.Draw("SAME")

    cutLine2 = TLine(0.67, len(arTTChi) / 3. + 0.5,
                     0.88, len(arTTChi) / 3. + 0.5)
    cutLine2.SetLineColor(kBlack)
    cutLine2.SetLineWidth(2)
    cutLine2.Draw("SAME")

    l = ROOT.TLatex(
        # 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
        0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Private Work}}")
    l.SetNDC()
    # l.Draw()

    lum = ROOT.TLatex(.62, .95, "%.1f fb^{-1} (%s TeV)" %
                      (aux.intLumi / 1000., aux.Label.cmsEnergy))
    lum.SetNDC()
    lum.Draw()

    l4 = ROOT.TLatex(
        0.25, .89, "#scale[0.76]{#font[52]{t#bar{t}(#gamma) Control Region}}")
    l4.SetNDC()
    l4.Draw()

    leg = TLegend(0.18, 0.45, 0.45, 0.55)
    leg.AddEntry(line, "#alpha (int. method)", "l")
    leg.AddEntry(ge, "stat. error (int. method)", "f")
    leg.AddEntry(gr, "#alpha (#chi^{2} method)", "lep")
    leg.AddEntry(lineMEAN, "mean (#chi^{2} method)", "l")
    leg.SetTextSize(0.03)
    # leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.Draw("same")

    c.SetGridy()
    c.Update()

    c.SaveAs('plots_CR/chi/TT_Compare.pdf')


def drawWZ():
    valueWZInt = pklWZ["LL"]["eta1"][0]
    errWZInt = pklWZ["LL"]["eta1"][1]
    arWZChi = []
    arWZChiUp = []
    arWZChiDn = []
    arWZChiName = []

    for comb in pklWZchiBig:
        for variable in pklWZchiBig[comb]:
            arWZChi.append(pklWZchiBig[comb][variable]["value"])
            arWZChiUp.append(pklWZchiBig[comb][variable]["erUp"])
            arWZChiDn.append(pklWZchiBig[comb][variable]["erDown"])
            arWZChiName.append(variable)
    for comb in pklWZchi:
        for variable in pklWZchi[comb]:
            arWZChi.append(pklWZchi[comb][variable]["value"])
            arWZChiUp.append(pklWZchi[comb][variable]["erUp"])
            arWZChiDn.append(pklWZchi[comb][variable]["erDown"])
            arWZChiName.append(variable)
    for comb in pklWZchiSmall:
        for variable in pklWZchiSmall[comb]:
            arWZChi.append(pklWZchiSmall[comb][variable]["value"])
            arWZChiUp.append(pklWZchiSmall[comb][variable]["erUp"])
            arWZChiDn.append(pklWZchiSmall[comb][variable]["erDown"])
            arWZChiName.append(variable)

    # print arWZChi
    # print arWZChiUp
    # print arWZChiDn
    # print arWZChiName

    c = TCanvas("canvas", "", 800, 1000)

    x = [valueWZInt]
    ex = [errWZInt]
    y = [len(arWZChi) / 2.]
    eyU = [round(len(arWZChi) / 2.) + 1]
    eyD = [len(arWZChi) / 2.]
    ge = TGraphAsymmErrors()
    for i in range(len(x)):
        ge.SetPoint(i, x[i], y[i])
        ge.SetPointError(i, ex[i], ex[i], eyD[i], eyU[i])
    # ge.SetFillStyle(3004)
    ge.SetFillColorAlpha(ROOT.kBlue, 0.3)
    # ge.SetFillColor(ROOT.kBlue)
    ge.SetLineColor(ROOT.kWhite)
    ge.Draw("same a2")

    gr = TGraphAsymmErrors()

    for i in range(len(arWZChi)):
        y = float(i + 1)
        gr.SetPoint(i, arWZChi[i], y)
        gr.SetPointError(i, arWZChiDn[i], arWZChiUp[i], 0., 0.)

    gr.SetTitle("; Scale Factor #alpha_{WZ}; ")
    ge.SetTitle("; Scale Factor #alpha_{WZ}; ")
    # gr.SetMarkerColor(4)
    gr.SetMarkerColor(kRed)
    gr.SetMarkerStyle(21)
    gr.Draw("sameP")

    line = TLine(valueWZInt, 0, valueWZInt, len(arWZChi) + 1)
    line.SetLineColor(kBlue)
    line.SetLineWidth(2)
    line.Draw("SAME")

    mean = np.mean(arWZChi)
    lineMEAN = TLine(mean, 0, mean, len(arWZChi) + 1)
    lineMEAN.SetLineColor(kRed)
    lineMEAN.SetLineWidth(2)
    lineMEAN.Draw("SAME")

    ge.GetXaxis().SetLimits(0.9, 1.3)
    ge.GetYaxis().SetLabelOffset(1)
    ge.GetYaxis().SetNdivisions(len(arWZChi) * 2)

    #t = TText()
    # t.SetTextAlign(32)
    # t.SetTextSize(0.035)
    # t.SetTextFont(72)
    # for i in range(len(arWZChi)):
    # t.DrawText(0.98,i+1,arWZChiName[i])

    #t = TText()
    t = TLatex()
    t.SetTextAlign(32)
    t.SetTextSize(0.02)
    # t.SetTextFont(72)
    for i in range(len(arWZChi)):
        t.DrawLatex(0.88, (i + 1), nameDict[arWZChiName[i]])
        # t.DrawText(0.8,i+1,arTTChiName[i])

    smallT = TText()
    smallT.SetTextAngle(90)
    smallT.SetTextSize(0.025)
    smallT.SetTextFont(42)
    smallT.DrawText(0.84, 3. * len(arWZChi) / 3. - 6, "double binning")
    smallT.DrawText(0.84, 2 * len(arWZChi) / 3. - 6, "normal binning")
    smallT.DrawText(0.84, len(arWZChi) / 3. - 6, "half binning")

    cutLine = TLine(0.9, 2 * len(arWZChi) / 3. + 0.5,
                    1.4, 2 * len(arWZChi) / 3. + 0.5)
    cutLine.SetLineColor(kBlack)
    cutLine.SetLineWidth(2)
    cutLine.Draw("SAME")

    cutLine2 = TLine(0.9, len(arWZChi) / 3. + 0.5,
                     1.4, len(arWZChi) / 3. + 0.5)
    cutLine2.SetLineColor(kBlack)
    cutLine2.SetLineWidth(2)
    cutLine2.Draw("SAME")

    l = ROOT.TLatex(
        # 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
        0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Private Work}}")
    l.SetNDC()
    # l.Draw()

    lum = ROOT.TLatex(.62, .95, "%.1f fb^{-1} (%s TeV)" %
                      (aux.intLumi / 1000., aux.Label.cmsEnergy))
    lum.SetNDC()
    lum.Draw()

    l4 = ROOT.TLatex(0.25, .89, "#scale[0.76]{#font[52]{WZ Control Region}}")
    l4.SetNDC()
    l4.Draw()

    leg = TLegend(0.18, 0.45, 0.45, 0.55)
    leg.AddEntry(line, "#alpha (int. method)", "l")
    leg.AddEntry(ge, "stat. error (int. method)", "f")
    leg.AddEntry(gr, "#alpha (#chi^{2} method)", "lep")
    leg.AddEntry(lineMEAN, "mean (#chi^{2} method)", "l")
    leg.SetTextSize(0.03)
    # leg.SetFillStyle(0)
    leg.Draw("same")

    c.SetGridy()
    c.Update()

    c.SaveAs('plots_CR/chi/WZ_Compare.pdf')


def main():
    drawZZ()
    drawTT()
    drawWZ()
    drawDYLL()
    # drawDYEE()
    # drawDYMM()


if __name__ == "__main__":
    main()
