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

pklZZ = pkl.load( open( "plots_CR_zz/factors/CRZZ.pkl", "rb" ) )
pklDY = pkl.load( open( "plots_CR_dy/factors/CRDY.pkl", "rb" ) )
pklTT = pkl.load( open( "plots_CR_tt/factors/CRTT.pkl", "rb" ) )
#pklTTsmall = pkl.load( open( "plots_CR_tt/factors/CRTT_small.pkl", "rb" ) )
#pklTTbig = pkl.load( open( "plots_CR_tt/factors/CRTT_big.pkl", "rb" ) )
pklWZ = pkl.load( open( "plots_CR_wz/factors/CRWZ.pkl", "rb" ) )

pklZZchi = pkl.load( open( "plots_CR/chi/ZZ_chi.pkl", "rb" ) )
pklZZchiBig = pkl.load( open( "plots_CR/chi/ZZ_chi_big.pkl", "rb" ) )
pklZZchiSmall = pkl.load( open( "plots_CR/chi/ZZ_chi_small.pkl", "rb" ) )

#pklDYEEchi = pkl.load( open( "plots_CR/chi/DY_chi_EE.pkl", "rb" ) )
#pklDYMMchi = pkl.load( open( "plots_CR/chi/DY_chi_MM.pkl", "rb" ) )
#pklDYLLchi = pkl.load( open( "plots_CR/chi/DY_chi_LL.pkl", "rb" ) )

pklTTchiJESu = pkl.load( open( "plots_CR/chi/TT_chi_JESu.pkl", "rb" ) )
pklTTchiJESd = pkl.load( open( "plots_CR/chi/TT_chi_JESd.pkl", "rb" ) )
pklTTchiJERu = pkl.load( open( "plots_CR/chi/TT_chi_JERu.pkl", "rb" ) )
pklTTchiJERd = pkl.load( open( "plots_CR/chi/TT_chi_JERd.pkl", "rb" ) )
pklTTchilepSFu = pkl.load( open( "plots_CR/chi/TT_chi_lepSFu.pkl", "rb" ) )
pklTTchilepSFd = pkl.load( open( "plots_CR/chi/TT_chi_lepSFd.pkl", "rb" ) )
pklTTchiphotonSFu = pkl.load( open( "plots_CR/chi/TT_chi_photonSFu.pkl", "rb" ) )
pklTTchiphotonSFd = pkl.load( open( "plots_CR/chi/TT_chi_photonSFd.pkl", "rb" ) )

pklZZchiJESu = pkl.load( open( "plots_CR/chi/ZZ_chi_JESu.pkl", "rb" ) )
pklZZchiJESd = pkl.load( open( "plots_CR/chi/ZZ_chi_JESd.pkl", "rb" ) )
pklZZchiJERu = pkl.load( open( "plots_CR/chi/ZZ_chi_JERu.pkl", "rb" ) )
pklZZchiJERd = pkl.load( open( "plots_CR/chi/ZZ_chi_JERd.pkl", "rb" ) )
pklZZchilepSFu = pkl.load( open( "plots_CR/chi/ZZ_chi_lepSFu.pkl", "rb" ) )
pklZZchilepSFd = pkl.load( open( "plots_CR/chi/ZZ_chi_lepSFd.pkl", "rb" ) )
pklZZchiphotonSFu = pkl.load( open( "plots_CR/chi/ZZ_chi_photonSFu.pkl", "rb" ) )
pklZZchiphotonSFd = pkl.load( open( "plots_CR/chi/ZZ_chi_photonSFd.pkl", "rb" ) )


pklWZchiJESu = pkl.load( open( "plots_CR/chi/WZ_chi_JESu.pkl", "rb" ) )
pklWZchiJESd = pkl.load( open( "plots_CR/chi/WZ_chi_JESd.pkl", "rb" ) )
pklWZchiJERu = pkl.load( open( "plots_CR/chi/WZ_chi_JERu.pkl", "rb" ) )
pklWZchiJERd = pkl.load( open( "plots_CR/chi/WZ_chi_JERd.pkl", "rb" ) )
pklWZchilepSFu = pkl.load( open( "plots_CR/chi/WZ_chi_lepSFu.pkl", "rb" ) )
pklWZchilepSFd = pkl.load( open( "plots_CR/chi/WZ_chi_lepSFd.pkl", "rb" ) )
pklWZchiphotonSFu = pkl.load( open( "plots_CR/chi/WZ_chi_photonSFu.pkl", "rb" ) )
pklWZchiphotonSFd = pkl.load( open( "plots_CR/chi/WZ_chi_photonSFd.pkl", "rb" ) )


pklDYchiMMJESu = pkl.load( open( "plots_CR/chi/DY_chi_MM_JESu.pkl", "rb" ) )
pklDYchiMMJESd = pkl.load( open( "plots_CR/chi/DY_chi_MM_JESd.pkl", "rb" ) )
pklDYchiMMJERu = pkl.load( open( "plots_CR/chi/DY_chi_MM_JERu.pkl", "rb" ) )
pklDYchiMMJERd = pkl.load( open( "plots_CR/chi/DY_chi_MM_JERd.pkl", "rb" ) )
pklDYchiMMlepSFu = pkl.load( open( "plots_CR/chi/DY_chi_MM_lepSFu.pkl", "rb" ) )
pklDYchiMMlepSFd = pkl.load( open( "plots_CR/chi/DY_chi_MM_lepSFd.pkl", "rb" ) )
pklDYchiMMphotonSFu = pkl.load( open( "plots_CR/chi/DY_chi_MM_photonSFu.pkl", "rb" ) )
pklDYchiMMphotonSFd = pkl.load( open( "plots_CR/chi/DY_chi_MM_photonSFd.pkl", "rb" ) )
pklDYchiLLJESu = pkl.load( open( "plots_CR/chi/DY_chi_LL_JESu.pkl", "rb" ) )
pklDYchiLLJESd = pkl.load( open( "plots_CR/chi/DY_chi_LL_JESd.pkl", "rb" ) )
pklDYchiLLJERu = pkl.load( open( "plots_CR/chi/DY_chi_LL_JERu.pkl", "rb" ) )
pklDYchiLLJERd = pkl.load( open( "plots_CR/chi/DY_chi_LL_JERd.pkl", "rb" ) )
pklDYchiLLlepSFu = pkl.load( open( "plots_CR/chi/DY_chi_LL_lepSFu.pkl", "rb" ) )
pklDYchiLLlepSFd = pkl.load( open( "plots_CR/chi/DY_chi_LL_lepSFd.pkl", "rb" ) )
pklDYchiLLphotonSFu = pkl.load( open( "plots_CR/chi/DY_chi_LL_photonSFu.pkl", "rb" ) )
pklDYchiLLphotonSFd = pkl.load( open( "plots_CR/chi/DY_chi_LL_photonSFd.pkl", "rb" ) )
pklDYchiEEJESu = pkl.load( open( "plots_CR/chi/DY_chi_EE_JESu.pkl", "rb" ) )
pklDYchiEEJESd = pkl.load( open( "plots_CR/chi/DY_chi_EE_JESd.pkl", "rb" ) )
pklDYchiEEJERu = pkl.load( open( "plots_CR/chi/DY_chi_EE_JERu.pkl", "rb" ) )
pklDYchiEEJERd = pkl.load( open( "plots_CR/chi/DY_chi_EE_JERd.pkl", "rb" ) )
pklDYchiEElepSFu = pkl.load( open( "plots_CR/chi/DY_chi_EE_lepSFu.pkl", "rb" ) )
pklDYchiEElepSFd = pkl.load( open( "plots_CR/chi/DY_chi_EE_lepSFd.pkl", "rb" ) )
pklDYchiEEphotonSFu = pkl.load( open( "plots_CR/chi/DY_chi_EE_photonSFu.pkl", "rb" ) )
pklDYchiEEphotonSFd = pkl.load( open( "plots_CR/chi/DY_chi_EE_photonSFd.pkl", "rb" ) )

pklWZchi = pkl.load( open( "plots_CR/chi/WZ_chi.pkl", "rb" ) )
pklWZchiBig = pkl.load( open( "plots_CR/chi/WZ_chi_big.pkl", "rb" ) )
pklWZchiSmall = pkl.load( open( "plots_CR/chi/WZ_chi_small.pkl", "rb" ) )

pklDYchiLL = pkl.load( open( "plots_CR/chi/DY_chi_LL.pkl", "rb" ) )
pklDYchiLLBig = pkl.load( open( "plots_CR/chi/DY_chi_LL_big.pkl", "rb" ) )
pklDYchiLLSmall = pkl.load( open( "plots_CR/chi/DY_chi_LL_small.pkl", "rb" ) )
pklDYchiEE = pkl.load( open( "plots_CR/chi/DY_chi_EE.pkl", "rb" ) )
pklDYchiEEBig = pkl.load( open( "plots_CR/chi/DY_chi_EE_big.pkl", "rb" ) )
pklDYchiEESmall = pkl.load( open( "plots_CR/chi/DY_chi_EE_small.pkl", "rb" ) )
pklDYchiMM = pkl.load( open( "plots_CR/chi/DY_chi_MM.pkl", "rb" ) )
pklDYchiMMBig = pkl.load( open( "plots_CR/chi/DY_chi_MM_big.pkl", "rb" ) )
pklDYchiMMSmall = pkl.load( open( "plots_CR/chi/DY_chi_MM_small.pkl", "rb" ) )

nameDict={
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
    valueZZInt=pklZZ["LL"]["eta1"][0]
    errZZInt=pklZZ["LL"]["eta1"][1]
    arZZChi=[]
    arZZChiUp=[]
    arZZChiDn=[]
    arZZChiName=[]

    for comb in pklZZchiJESu:
        for variable in pklZZchiJESu[comb]:
            arZZChi.append(pklZZchiJESu[comb][variable]["value"])
            arZZChiUp.append(pklZZchiJESu[comb][variable]["erUp"])
            arZZChiDn.append(pklZZchiJESu[comb][variable]["erDown"])
            arZZChiName.append(variable)
    for comb in pklZZchiJESd:
        for variable in pklZZchiJESd[comb]:
            arZZChi.append(pklZZchiJESd[comb][variable]["value"])
            arZZChiUp.append(pklZZchiJESd[comb][variable]["erUp"])
            arZZChiDn.append(pklZZchiJESd[comb][variable]["erDown"])
            arZZChiName.append(variable)
    for comb in pklZZchiJERu:
        for variable in pklZZchiJERu[comb]:
            arZZChi.append(pklZZchiJERu[comb][variable]["value"])
            arZZChiUp.append(pklZZchiJERu[comb][variable]["erUp"])
            arZZChiDn.append(pklZZchiJERu[comb][variable]["erDown"])
            arZZChiName.append(variable)
    for comb in pklZZchiJERd:
        for variable in pklZZchiJERd[comb]:
            arZZChi.append(pklZZchiJERd[comb][variable]["value"])
            arZZChiUp.append(pklZZchiJERd[comb][variable]["erUp"])
            arZZChiDn.append(pklZZchiJERd[comb][variable]["erDown"])
            arZZChiName.append(variable)

    #print arZZChi
    #print arZZChiUp
    #print arZZChiDn
    #print arZZChiName

    c = TCanvas("canvas","",800,1200)

    x=[valueZZInt]
    ex=[errZZInt]
    y=[len(arZZChi)/2.]
    eyU=[round(len(arZZChi)/2.)+1]
    eyD=[len(arZZChi)/2.]
    ge=TGraphAsymmErrors()
    for i in range(len(x)):
        ge.SetPoint(i,x[i],y[i])
        ge.SetPointError(i,ex[i],ex[i],eyD[i],eyU[i])
    ge.SetFillStyle(3003)
    ge.SetFillColor(ROOT.kBlue)
    ge.Draw("same a2")


    gr = TGraphAsymmErrors()

    for i in range(len(arZZChi)):
        y= float(i+1)
        gr.SetPoint(i,arZZChi[i],y)
        gr.SetPointError(i,arZZChiDn[i],arZZChiUp[i],0.,0.)

    gr.SetTitle("; Scale Factor #alpha; ")
    ge.SetTitle("; Scale Factor #alpha; ")
    gr.SetMarkerColor(4)
    gr.SetMarkerStyle(21)
    gr.Draw("sameP")

    line=TLine(valueZZInt,0,valueZZInt,len(arZZChi)+1)
    line.SetLineColor(kBlue)
    line.SetLineWidth(2)
    line.Draw("SAME")


    mean=np.mean(arZZChi)
    lineMEAN=TLine(mean,0,mean,len(arZZChi)+1)
    lineMEAN.SetLineColor(kRed)
    lineMEAN.SetLineWidth(2)
    lineMEAN.Draw("SAME")

    ge.GetXaxis().SetLimits(1.05,1.75)
    ge.GetYaxis().SetLabelOffset(1)
    ge.GetYaxis().SetNdivisions(len(arZZChi))

    #t = TText()
    t = TLatex()
    t.SetTextAlign(32)
    t.SetTextSize(0.02)
    #t.SetTextFont(72)
    for i in range(len(arZZChi)):
        t.DrawLatex(1.04,(i+1),nameDict[arZZChiName[i]])
        #t.DrawText(0.8,i+1,arTTChiName[i])
    
    smallT=TText()
    smallT.SetTextAngle(90)
    smallT.SetTextSize(0.025)
    smallT.SetTextFont(42)
    smallT.DrawText(0.98,4.*len(arZZChi)/4. -6,"JESu")
    smallT.DrawText(0.98,3.*len(arZZChi)/4. -6,"JESd")
    smallT.DrawText(0.98,2.*len(arZZChi)/4. -6,"JERu")
    smallT.DrawText(0.98,len(arZZChi)/4. -6,"JERd")
    
    cutLine=TLine(0.9,3*len(arZZChi)/4. +0.5,1.75,3*len(arZZChi)/4. +0.5)
    cutLine.SetLineColor(kBlack)
    cutLine.SetLineWidth(2)
    cutLine.Draw("SAME")
    
    cutLine2=TLine(0.9,2*len(arZZChi)/4. +0.5,1.75,2*len(arZZChi)/4. +0.5)
    cutLine2.SetLineColor(kBlack)
    cutLine2.SetLineWidth(2)
    cutLine2.Draw("SAME")
    
    cutLine2=TLine(0.9,len(arZZChi)/4. +0.5,1.75,len(arZZChi)/4. +0.5)
    cutLine2.SetLineColor(kBlack)
    cutLine2.SetLineWidth(2)
    cutLine2.Draw("SAME")

    l = ROOT.TLatex( 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
    l.SetNDC()
    l.Draw()

    lum = ROOT.TLatex( .62, .95, "%.1f fb^{-1} (%s TeV)"%(aux.intLumi/1000., aux.Label.cmsEnergy) )
    lum.SetNDC()
    lum.Draw()
    
    l4 = ROOT.TLatex( 0.17, .9, "#scale[0.66]{#font[52]{ZZ Control Region}}")
    l4.SetNDC()
    l4.Draw()
    
    
    leg=TLegend(0.18,0.45,0.4,0.55)
    leg.AddEntry(line,"#alpha from int. method","l")
    leg.AddEntry(ge,"stat. error from int. method","f")
    leg.AddEntry(gr,"#alpha from #chi^{2} method","lep")
    leg.AddEntry(lineMEAN,"mean from #chi^{2} method","l")
    leg.SetTextSize(0.02)
    #leg.SetFillStyle(0)
    leg.Draw("same")

    c.SetGridy()
    c.Update()

    c.SaveAs('plots_CR/chi/ZZ_Compare_JETsyst.pdf')
    
def drawDYLL():
    valueDYInt=pklDY["LL"]["eta1"][0]
    errDYInt=pklDY["LL"]["eta1"][1]
    arDYChi=[]
    arDYChiUp=[]
    arDYChiDn=[]
    arDYChiName=[]

    for comb in pklDYchiLLJESu:
        for variable in pklDYchiLLJESu[comb]:
            arDYChi.append(pklDYchiLLJESu[comb][variable]["value"])
            arDYChiUp.append(pklDYchiLLJESu[comb][variable]["erUp"])
            arDYChiDn.append(pklDYchiLLJESu[comb][variable]["erDown"])
            arDYChiName.append(variable)
    for comb in pklDYchiLLJESd:
        for variable in pklDYchiLLJESd[comb]:
            arDYChi.append(pklDYchiLLJESd[comb][variable]["value"])
            arDYChiUp.append(pklDYchiLLJESd[comb][variable]["erUp"])
            arDYChiDn.append(pklDYchiLLJESd[comb][variable]["erDown"])
            arDYChiName.append(variable)
    for comb in pklDYchiLLJERu:
        for variable in pklDYchiLLJERu[comb]:
            arDYChi.append(pklDYchiLLJERu[comb][variable]["value"])
            arDYChiUp.append(pklDYchiLLJERu[comb][variable]["erUp"])
            arDYChiDn.append(pklDYchiLLJERu[comb][variable]["erDown"])
            arDYChiName.append(variable)
    for comb in pklDYchiLLJERd:
        for variable in pklDYchiLLJERd[comb]:
            arDYChi.append(pklDYchiLLJERd[comb][variable]["value"])
            arDYChiUp.append(pklDYchiLLJERd[comb][variable]["erUp"])
            arDYChiDn.append(pklDYchiLLJERd[comb][variable]["erDown"])
            arDYChiName.append(variable)

    #print arDYChi
    #print arDYChiUp
    #print arDYChiDn
    #print arDYChiName

    c = TCanvas("canvas","",800,1200)

    x=[valueDYInt]
    ex=[errDYInt]
    y=[len(arDYChi)/2.]
    eyU=[round(len(arDYChi)/2.)+1]
    eyD=[len(arDYChi)/2.]
    ge=TGraphAsymmErrors()
    for i in range(len(x)):
        ge.SetPoint(i,x[i],y[i])
        ge.SetPointError(i,ex[i],ex[i],eyD[i],eyU[i])
    ge.SetFillStyle(3003)
    ge.SetFillColor(ROOT.kBlue)
    ge.Draw("same a2")


    gr = TGraphAsymmErrors()

    for i in range(len(arDYChi)):
        y= float(i+1)
        gr.SetPoint(i,arDYChi[i],y)
        gr.SetPointError(i,arDYChiDn[i],arDYChiUp[i],0.,0.)

    gr.SetTitle("; Scale Factor #alpha; ")
    ge.SetTitle("; Scale Factor #alpha; ")
    gr.SetMarkerColor(4)
    gr.SetMarkerStyle(21)
    gr.Draw("sameP")

    line=TLine(valueDYInt,0,valueDYInt,len(arDYChi)+1)
    line.SetLineColor(kBlue)
    line.SetLineWidth(2)
    line.Draw("SAME")


    mean=np.mean(arDYChi)
    lineMEAN=TLine(mean,0,mean,len(arDYChi)+1)
    lineMEAN.SetLineColor(kRed)
    lineMEAN.SetLineWidth(2)
    lineMEAN.Draw("SAME")

    ge.GetXaxis().SetLimits(1.05,1.12)
    ge.GetYaxis().SetLabelOffset(1)
    ge.GetYaxis().SetNdivisions(len(arDYChi))

    #t = TText()
    t = TLatex()
    t.SetTextAlign(32)
    t.SetTextSize(0.02)
    #t.SetTextFont(72)
    for i in range(len(arDYChi)):
        t.DrawLatex(1.047,(i+1),nameDict[arDYChiName[i]])
        #t.DrawText(0.8,i+1,arTTChiName[i])
    
    smallT=TText()
    smallT.SetTextAngle(90)
    smallT.SetTextSize(0.025)
    smallT.SetTextFont(42)
    smallT.DrawText(1.04,4.*len(arDYChi)/4. -6,"JESu")
    smallT.DrawText(1.04,3*len(arDYChi)/4. -6,"JESd")
    smallT.DrawText(1.04,2*len(arDYChi)/4. -6,"JERu")
    smallT.DrawText(1.04,len(arDYChi)/4. -6,"JERd")
    
    cutLine=TLine(1.04,3*len(arDYChi)/4. +0.5,1.12,3*len(arDYChi)/4. +0.5)
    cutLine.SetLineColor(kBlack)
    cutLine.SetLineWidth(2)
    cutLine.Draw("SAME")
    
    cutLine2=TLine(1.04,2*len(arDYChi)/4. +0.5,1.12,2*len(arDYChi)/4. +0.5)
    cutLine2.SetLineColor(kBlack)
    cutLine2.SetLineWidth(2)
    cutLine2.Draw("SAME")
    
    cutLine2=TLine(1.04,len(arDYChi)/4. +0.5,1.12,len(arDYChi)/4. +0.5)
    cutLine2.SetLineColor(kBlack)
    cutLine2.SetLineWidth(2)
    cutLine2.Draw("SAME")

    l = ROOT.TLatex( 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
    l.SetNDC()
    l.Draw()

    lum = ROOT.TLatex( .62, .95, "%.1f fb^{-1} (%s TeV)"%(aux.intLumi/1000., aux.Label.cmsEnergy) )
    lum.SetNDC()
    lum.Draw()
    
    l4 = ROOT.TLatex( 0.17, .9, "#scale[0.66]{#font[52]{DY Control Region}}")
    l4.SetNDC()
    l4.Draw()
    
    
    leg=TLegend(0.18,0.45,0.4,0.55)
    leg.AddEntry(line,"#alpha from int. method","l")
    leg.AddEntry(ge,"stat. error from int. method","f")
    leg.AddEntry(gr,"#alpha from #chi^{2} method","lep")
    leg.AddEntry(lineMEAN,"mean from #chi^{2} method","l")
    leg.SetTextSize(0.02)
    #leg.SetFillStyle(0)
    leg.Draw("same")

    c.SetGridy()
    c.Update()

    c.SaveAs('plots_CR/chi/DY_CompareLL_JETsyst.pdf')
def drawDYEE():
    valueDYInt=pklDY["EE"]["eta1"][0]
    errDYInt=pklDY["EE"]["eta1"][1]
    arDYChi=[]
    arDYChiUp=[]
    arDYChiDn=[]
    arDYChiName=[]

    for comb in pklDYchiEEJESu:
        for variable in pklDYchiEEJESu[comb]:
            arDYChi.append(pklDYchiEEJESu[comb][variable]["value"])
            arDYChiUp.append(pklDYchiEEJESu[comb][variable]["erUp"])
            arDYChiDn.append(pklDYchiEEJESu[comb][variable]["erDown"])
            arDYChiName.append(variable)
    for comb in pklDYchiEEJESd:
        for variable in pklDYchiEEJESd[comb]:
            arDYChi.append(pklDYchiEEJESd[comb][variable]["value"])
            arDYChiUp.append(pklDYchiEEJESd[comb][variable]["erUp"])
            arDYChiDn.append(pklDYchiEEJESd[comb][variable]["erDown"])
            arDYChiName.append(variable)
    for comb in pklDYchiEEJERu:
        for variable in pklDYchiEEJERu[comb]:
            arDYChi.append(pklDYchiEEJERu[comb][variable]["value"])
            arDYChiUp.append(pklDYchiEEJERu[comb][variable]["erUp"])
            arDYChiDn.append(pklDYchiEEJERu[comb][variable]["erDown"])
            arDYChiName.append(variable)
    for comb in pklDYchiEEJERd:
        for variable in pklDYchiEEJERd[comb]:
            arDYChi.append(pklDYchiEEJERd[comb][variable]["value"])
            arDYChiUp.append(pklDYchiEEJERd[comb][variable]["erUp"])
            arDYChiDn.append(pklDYchiEEJERd[comb][variable]["erDown"])
            arDYChiName.append(variable)

    #print arDYChi
    #print arDYChiUp
    #print arDYChiDn
    #print arDYChiName

    c = TCanvas("canvas","",800,1200)

    x=[valueDYInt]
    ex=[errDYInt]
    y=[len(arDYChi)/2.]
    eyU=[round(len(arDYChi)/2.)+1]
    eyD=[len(arDYChi)/2.]
    ge=TGraphAsymmErrors()
    for i in range(len(x)):
        ge.SetPoint(i,x[i],y[i])
        ge.SetPointError(i,ex[i],ex[i],eyD[i],eyU[i])
    ge.SetFillStyle(3003)
    ge.SetFillColor(ROOT.kBlue)
    ge.Draw("same a2")


    gr = TGraphAsymmErrors()

    for i in range(len(arDYChi)):
        y= float(i+1)
        gr.SetPoint(i,arDYChi[i],y)
        gr.SetPointError(i,arDYChiDn[i],arDYChiUp[i],0.,0.)

    gr.SetTitle("; Scale Factor #alpha; ")
    ge.SetTitle("; Scale Factor #alpha; ")
    gr.SetMarkerColor(4)
    gr.SetMarkerStyle(21)
    gr.Draw("sameP")

    line=TLine(valueDYInt,0,valueDYInt,len(arDYChi)+1)
    line.SetLineColor(kBlue)
    line.SetLineWidth(2)
    line.Draw("SAME")


    mean=np.mean(arDYChi)
    lineMEAN=TLine(mean,0,mean,len(arDYChi)+1)
    lineMEAN.SetLineColor(kRed)
    lineMEAN.SetLineWidth(2)
    lineMEAN.Draw("SAME")

    ge.GetXaxis().SetLimits(1.02,1.12)
    ge.GetYaxis().SetLabelOffset(1)
    ge.GetYaxis().SetNdivisions(len(arDYChi))

    #t = TText()
    t = TLatex()
    t.SetTextAlign(32)
    t.SetTextSize(0.02)
    #t.SetTextFont(72)
    for i in range(len(arDYChi)):
        t.DrawLatex(1.015,(i+1),nameDict[arDYChiName[i]])
        #t.DrawText(0.8,i+1,arTTChiName[i])
    
    smallT=TText()
    smallT.SetTextAngle(90)
    smallT.SetTextSize(0.025)
    smallT.SetTextFont(42)
    smallT.DrawText(1.005,4.*len(arDYChi)/4. -6,"JESu")
    smallT.DrawText(1.005,3*len(arDYChi)/4. -6,"JESd")
    smallT.DrawText(1.005,2*len(arDYChi)/4. -6,"JERu")
    smallT.DrawText(1.005,len(arDYChi)/4. -6,"JERd")
    
    cutLine=TLine(1.005,3*len(arDYChi)/4. +0.5,1.75,3*len(arDYChi)/4. +0.5)
    cutLine.SetLineColor(kBlack)
    cutLine.SetLineWidth(2)
    cutLine.Draw("SAME")
    
    cutLine2=TLine(1.005,2*len(arDYChi)/4. +0.5,1.75,2*len(arDYChi)/4. +0.5)
    cutLine2.SetLineColor(kBlack)
    cutLine2.SetLineWidth(2)
    cutLine2.Draw("SAME")
    
    cutLine2=TLine(1.005,len(arDYChi)/4. +0.5,1.75,len(arDYChi)/4. +0.5)
    cutLine2.SetLineColor(kBlack)
    cutLine2.SetLineWidth(2)
    cutLine2.Draw("SAME")

    l = ROOT.TLatex( 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
    l.SetNDC()
    l.Draw()

    lum = ROOT.TLatex( .62, .95, "%.1f fb^{-1} (%s TeV)"%(aux.intLumi/1000., aux.Label.cmsEnergy) )
    lum.SetNDC()
    lum.Draw()
    
    l4 = ROOT.TLatex( 0.17, .9, "#scale[0.66]{#font[52]{DY Control Region}}")
    l4.SetNDC()
    l4.Draw()
    
    
    leg=TLegend(0.18,0.45,0.4,0.55)
    leg.AddEntry(line,"#alpha from int. method","l")
    leg.AddEntry(ge,"stat. error from int. method","f")
    leg.AddEntry(gr,"#alpha from #chi^{2} method","lep")
    leg.AddEntry(lineMEAN,"mean from #chi^{2} method","l")
    leg.SetTextSize(0.02)
    #leg.SetFillStyle(0)
    leg.Draw("same")

    c.SetGridy()
    c.Update()

    c.SaveAs('plots_CR/chi/DY_CompareEE_JETsyst.pdf')
def drawDYMM():
    valueDYInt=pklDY["MM"]["eta1"][0]
    errDYInt=pklDY["MM"]["eta1"][1]
    arDYChi=[]
    arDYChiUp=[]
    arDYChiDn=[]
    arDYChiName=[]

    for comb in pklDYchiMMJESu:
        for variable in pklDYchiMMJESu[comb]:
            arDYChi.append(pklDYchiMMJESu[comb][variable]["value"])
            arDYChiUp.append(pklDYchiMMJESu[comb][variable]["erUp"])
            arDYChiDn.append(pklDYchiMMJESu[comb][variable]["erDown"])
            arDYChiName.append(variable)
    for comb in pklDYchiMMJESd:
        for variable in pklDYchiMMJESd[comb]:
            arDYChi.append(pklDYchiMMJESd[comb][variable]["value"])
            arDYChiUp.append(pklDYchiMMJESd[comb][variable]["erUp"])
            arDYChiDn.append(pklDYchiMMJESd[comb][variable]["erDown"])
            arDYChiName.append(variable)
    for comb in pklDYchiMMJERu:
        for variable in pklDYchiMMJERu[comb]:
            arDYChi.append(pklDYchiMMJERu[comb][variable]["value"])
            arDYChiUp.append(pklDYchiMMJERu[comb][variable]["erUp"])
            arDYChiDn.append(pklDYchiMMJERu[comb][variable]["erDown"])
            arDYChiName.append(variable)
    for comb in pklDYchiMMJERd:
        for variable in pklDYchiMMJERd[comb]:
            arDYChi.append(pklDYchiMMJERd[comb][variable]["value"])
            arDYChiUp.append(pklDYchiMMJERd[comb][variable]["erUp"])
            arDYChiDn.append(pklDYchiMMJERd[comb][variable]["erDown"])
            arDYChiName.append(variable)

    #print arDYChi
    #print arDYChiUp
    #print arDYChiDn
    #print arDYChiName

    c = TCanvas("canvas","",800,1200)

    x=[valueDYInt]
    ex=[errDYInt]
    y=[len(arDYChi)/2.]
    eyU=[round(len(arDYChi)/2.)+1]
    eyD=[len(arDYChi)/2.]
    ge=TGraphAsymmErrors()
    for i in range(len(x)):
        ge.SetPoint(i,x[i],y[i])
        ge.SetPointError(i,ex[i],ex[i],eyD[i],eyU[i])
    ge.SetFillStyle(3003)
    ge.SetFillColor(ROOT.kBlue)
    ge.Draw("same a2")


    gr = TGraphAsymmErrors()

    for i in range(len(arDYChi)):
        y= float(i+1)
        gr.SetPoint(i,arDYChi[i],y)
        gr.SetPointError(i,arDYChiDn[i],arDYChiUp[i],0.,0.)

    gr.SetTitle("; Scale Factor #alpha; ")
    ge.SetTitle("; Scale Factor #alpha; ")
    gr.SetMarkerColor(4)
    gr.SetMarkerStyle(21)
    gr.Draw("sameP")

    line=TLine(valueDYInt,0,valueDYInt,len(arDYChi)+1)
    line.SetLineColor(kBlue)
    line.SetLineWidth(2)
    line.Draw("SAME")


    mean=np.mean(arDYChi)
    lineMEAN=TLine(mean,0,mean,len(arDYChi)+1)
    lineMEAN.SetLineColor(kRed)
    lineMEAN.SetLineWidth(2)
    lineMEAN.Draw("SAME")

    ge.GetXaxis().SetLimits(1.05,1.12)
    ge.GetYaxis().SetLabelOffset(1)
    ge.GetYaxis().SetNdivisions(len(arDYChi))

    #t = TText()
    t = TLatex()
    t.SetTextAlign(32)
    t.SetTextSize(0.02)
    #t.SetTextFont(72)
    for i in range(len(arDYChi)):
        t.DrawLatex(1.047,(i+1),nameDict[arDYChiName[i]])
        #t.DrawText(0.8,i+1,arTTChiName[i])
    
    smallT=TText()
    smallT.SetTextAngle(90)
    smallT.SetTextSize(0.025)
    smallT.SetTextFont(42)
    smallT.DrawText(1.04,4.*len(arDYChi)/4. -6,"NBins * 2")
    smallT.DrawText(1.04,3*len(arDYChi)/4. -6,"normal binning")
    smallT.DrawText(1.04,2*len(arDYChi)/4. -6,"NBins / 2")
    smallT.DrawText(1.04,len(arDYChi)/4. -6,"NBins / 2")
    
    cutLine=TLine(1.04,3*len(arDYChi)/4. +0.5,1.75,3*len(arDYChi)/4. +0.5)
    cutLine.SetLineColor(kBlack)
    cutLine.SetLineWidth(2)
    cutLine.Draw("SAME")
    
    cutLine2=TLine(1.04,2*len(arDYChi)/4. +0.5,1.75,2*len(arDYChi)/4. +0.5)
    cutLine2.SetLineColor(kBlack)
    cutLine2.SetLineWidth(2)
    cutLine2.Draw("SAME")
    
    cutLine2=TLine(1.04,len(arDYChi)/4. +0.5,1.75,len(arDYChi)/4. +0.5)
    cutLine2.SetLineColor(kBlack)
    cutLine2.SetLineWidth(2)
    cutLine2.Draw("SAME")

    l = ROOT.TLatex( 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
    l.SetNDC()
    l.Draw()

    lum = ROOT.TLatex( .62, .95, "%.1f fb^{-1} (%s TeV)"%(aux.intLumi/1000., aux.Label.cmsEnergy) )
    lum.SetNDC()
    lum.Draw()
    
    l4 = ROOT.TLatex( 0.17, .9, "#scale[0.66]{#font[52]{DY Control Region}}")
    l4.SetNDC()
    l4.Draw()
    
    
    leg=TLegend(0.18,0.45,0.4,0.55)
    leg.AddEntry(line,"#alpha from int. method","l")
    leg.AddEntry(ge,"stat. error from int. method","f")
    leg.AddEntry(gr,"#alpha from #chi^{2} method","lep")
    leg.AddEntry(lineMEAN,"mean from #chi^{2} method","l")
    leg.SetTextSize(0.02)
    #leg.SetFillStyle(0)
    leg.Draw("same")

    c.SetGridy()
    c.Update()

    c.SaveAs('plots_CR/chi/DY_CompareMM_JETsyst.pdf')

def drawTT():
    valueTTInt=pklTT["EM"]["eta1"][0]
    errTTInt=pklTT["EM"]["eta1"][1]
    arTTChi=[]
    arTTChiUp=[]
    arTTChiDn=[]
    arTTChiName=[]


    for comb in pklTTchiJESu:
        for variable in pklTTchiJESu[comb]:
            arTTChi.append(pklTTchiJESu[comb][variable]["value"])
            arTTChiUp.append(pklTTchiJESu[comb][variable]["erUp"])
            arTTChiDn.append(pklTTchiJESu[comb][variable]["erDown"])
            arTTChiName.append(variable)
    for comb in pklTTchiJESd:
        for variable in pklTTchiJESd[comb]:
            arTTChi.append(pklTTchiJESd[comb][variable]["value"])
            arTTChiUp.append(pklTTchiJESd[comb][variable]["erUp"])
            arTTChiDn.append(pklTTchiJESd[comb][variable]["erDown"])
            arTTChiName.append(variable)
            
    for comb in pklTTchiJERu:
        for variable in pklTTchiJERu[comb]:
            arTTChi.append(pklTTchiJERu[comb][variable]["value"])
            arTTChiUp.append(pklTTchiJERu[comb][variable]["erUp"])
            arTTChiDn.append(pklTTchiJERu[comb][variable]["erDown"])
            arTTChiName.append(variable)
            
    for comb in pklTTchiJERd:
        for variable in pklTTchiJERd[comb]:
            arTTChi.append(pklTTchiJERd[comb][variable]["value"])
            arTTChiUp.append(pklTTchiJERd[comb][variable]["erUp"])
            arTTChiDn.append(pklTTchiJERd[comb][variable]["erDown"])
            arTTChiName.append(variable)
            
    for comb in pklTTchilepSFu:
        for variable in pklTTchilepSFu[comb]:
            arTTChi.append(pklTTchilepSFu[comb][variable]["value"])
            arTTChiUp.append(pklTTchilepSFu[comb][variable]["erUp"])
            arTTChiDn.append(pklTTchilepSFu[comb][variable]["erDown"])
            arTTChiName.append(variable)
            
    for comb in pklTTchilepSFd:
        for variable in pklTTchilepSFd[comb]:
            arTTChi.append(pklTTchilepSFd[comb][variable]["value"])
            arTTChiUp.append(pklTTchilepSFd[comb][variable]["erUp"])
            arTTChiDn.append(pklTTchilepSFd[comb][variable]["erDown"])
            arTTChiName.append(variable)
            
    for comb in pklTTchiphotonSFu:
        for variable in pklTTchiphotonSFu[comb]:
            arTTChi.append(pklTTchiphotonSFu[comb][variable]["value"])
            arTTChiUp.append(pklTTchiphotonSFu[comb][variable]["erUp"])
            arTTChiDn.append(pklTTchiphotonSFu[comb][variable]["erDown"])
            arTTChiName.append(variable)
            
    for comb in pklTTchiphotonSFd:
        for variable in pklTTchiphotonSFd[comb]:
            arTTChi.append(pklTTchiphotonSFd[comb][variable]["value"])
            arTTChiUp.append(pklTTchiphotonSFd[comb][variable]["erUp"])
            arTTChiDn.append(pklTTchiphotonSFd[comb][variable]["erDown"])
            arTTChiName.append(variable)

    #print arTTChi
    #print arTTChiUp
    #print arTTChiDn
    #print arTTChiName

    c = TCanvas("canvas","",800,1500)

    x=[valueTTInt]
    ex=[errTTInt]
    y=[len(arTTChi)/2.]
    eyU=[round(len(arTTChi)/2.)+1]
    eyD=[len(arTTChi)/2.]
    ge=TGraphAsymmErrors()
    for i in range(len(x)):
        ge.SetPoint(i,x[i],y[i])
        ge.SetPointError(i,ex[i],ex[i],eyD[i],eyU[i])
    #ge.SetFillStyle(3003)
    ge.SetFillStyle(3004)
    ge.SetFillColor(ROOT.kBlue)
    ge.Draw("same a2")


    gr = TGraphAsymmErrors()

    for i in range(len(arTTChi)):
        y= float(i+1)
        gr.SetPoint(i,arTTChi[i],y)
        gr.SetPointError(i,arTTChiDn[i],arTTChiUp[i],0.,0.)

    gr.SetTitle("; Scale Factor #alpha; ")
    ge.SetTitle("; Scale Factor #alpha; ")
    #gr.SetMarkerColor(4)
    gr.SetMarkerColor(kRed)
    gr.SetMarkerStyle(21)
    gr.Draw("sameP")

    line=TLine(valueTTInt,0,valueTTInt,len(arTTChi)+1)
    line.SetLineColor(kBlue)
    line.SetLineWidth(2)
    line.Draw("SAME")


    mean=np.mean(arTTChi)
    lineMEAN=TLine(mean,0,mean,len(arTTChi)+1)
    lineMEAN.SetLineColor(kRed)
    lineMEAN.SetLineWidth(2)
    lineMEAN.Draw("SAME")

    ge.GetXaxis().SetLimits(0.8,1.1)
    ge.GetYaxis().SetLabelOffset(1)
    #ge.GetYaxis().SetNdivisions(len(arTTChi)*2)
    ge.GetYaxis().SetNdivisions(len(arTTChi)*4)

    #t = TText()
    t = TLatex()
    t.SetTextAlign(32)
    t.SetTextSize(0.015)
    #t.SetTextFont(72)
    for i in range(len(arTTChi)):
        t.DrawLatex(0.79,(i+1),nameDict[arTTChiName[i]])
        #t.DrawText(0.8,i+1,arTTChiName[i])
    
    smallT=TText()
    smallT.SetTextAngle(90)
    smallT.SetTextSize(0.02)
    smallT.SetTextFont(42)
    
    smallT.DrawText(0.76,1.*len(arTTChi)/8. -6,"JESu")
    smallT.DrawText(0.76,2.*len(arTTChi)/8. -6,"JESd")
    smallT.DrawText(0.76,3.*len(arTTChi)/8. -6,"JERu")
    smallT.DrawText(0.76,4.*len(arTTChi)/8. -6,"JERd")
    smallT.DrawText(0.76,5.*len(arTTChi)/8. -6,"lepSFu")
    smallT.DrawText(0.76,6*len(arTTChi)/8. -6,"lepSFd")
    smallT.DrawText(0.76,7*len(arTTChi)/8. -6,"photonSFu")
    smallT.DrawText(0.76,8*len(arTTChi)/8. -6,"photonSFd")
    
    cutLine=TLine(0.75,3.*len(arTTChi)/8. +0.5,1.1,3.*len(arTTChi)/8. +0.5)
    cutLine.SetLineColor(kBlack)
    cutLine.SetLineWidth(1)
    cutLine.Draw("SAME")
    
    cutLine2=TLine(0.75,2.*len(arTTChi)/8. +0.5,1.1,2.*len(arTTChi)/8. +0.5)
    cutLine2.SetLineColor(kBlack)
    cutLine2.SetLineWidth(1)
    cutLine2.Draw("SAME")
    
    cutLine3=TLine(0.75,len(arTTChi)/8. +0.5,1.1,len(arTTChi)/8. +0.5)
    cutLine3.SetLineColor(kBlack)
    cutLine3.SetLineWidth(1)
    cutLine3.Draw("SAME")
    cutLine4=TLine(0.75,4.*len(arTTChi)/8. +0.5,1.1,4.*len(arTTChi)/8. +0.5)
    cutLine4.SetLineColor(kBlack)
    cutLine4.SetLineWidth(1)
    cutLine4.Draw("SAME")
    cutLine5=TLine(0.75,5.*len(arTTChi)/8. +0.5,1.1,5.*len(arTTChi)/8. +0.5)
    cutLine5.SetLineColor(kBlack)
    cutLine5.SetLineWidth(1)
    cutLine5.Draw("SAME")
    cutLine6=TLine(0.75,6.*len(arTTChi)/8. +0.5,1.1,6.*len(arTTChi)/8. +0.5)
    cutLine6.SetLineColor(kBlack)
    cutLine6.SetLineWidth(1)
    cutLine6.Draw("SAME")
    cutLine7=TLine(0.75,7.*len(arTTChi)/8. +0.5,1.1,7.*len(arTTChi)/8. +0.5)
    cutLine7.SetLineColor(kBlack)
    cutLine7.SetLineWidth(1)
    cutLine7.Draw("SAME")

    l = ROOT.TLatex( 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
    l.SetNDC()
    l.Draw()

    lum = ROOT.TLatex( .62, .95, "%.1f fb^{-1} (%s TeV)"%(aux.intLumi/1000., aux.Label.cmsEnergy) )
    lum.SetNDC()
    lum.Draw()
    
    l4 = ROOT.TLatex( 0.17, .9, "#scale[0.66]{#font[52]{t#bar{t}(+#gamma) Control Region}}")
    l4.SetNDC()
    l4.Draw()
    
    
    leg=TLegend(0.18,0.45,0.4,0.55)
    leg.AddEntry(line,"#alpha from int. method","l")
    leg.AddEntry(ge,"stat. error from int. method","f")
    leg.AddEntry(gr,"#alpha from #chi^{2} method","lep")
    leg.AddEntry(lineMEAN,"mean from #chi^{2} method","l")
    leg.SetTextSize(0.02)
    #leg.SetFillStyle(0)
    leg.Draw("same")
    
    c.SetGridy()
    c.Update()

    c.SaveAs('plots_CR/chi/TT_Compare_JETsyst.pdf')
    
def drawWZ():
    valueWZInt=pklWZ["LL"]["eta1"][0]
    errWZInt=pklWZ["LL"]["eta1"][1]
    arWZChi=[]
    arWZChiUp=[]
    arWZChiDn=[]
    arWZChiName=[]

    for comb in pklWZchiJESu:
        for variable in pklWZchiJESu[comb]:
            arWZChi.append(pklWZchiJESu[comb][variable]["value"])
            arWZChiUp.append(pklWZchiJESu[comb][variable]["erUp"])
            arWZChiDn.append(pklWZchiJESu[comb][variable]["erDown"])
            arWZChiName.append(variable)
    for comb in pklWZchiJESd:
        for variable in pklWZchiJESd[comb]:
            arWZChi.append(pklWZchiJESd[comb][variable]["value"])
            arWZChiUp.append(pklWZchiJESd[comb][variable]["erUp"])
            arWZChiDn.append(pklWZchiJESd[comb][variable]["erDown"])
            arWZChiName.append(variable)
    for comb in pklWZchiJERu:
        for variable in pklWZchiJERu[comb]:
            arWZChi.append(pklWZchiJERu[comb][variable]["value"])
            arWZChiUp.append(pklWZchiJERu[comb][variable]["erUp"])
            arWZChiDn.append(pklWZchiJERu[comb][variable]["erDown"])
            arWZChiName.append(variable)
    for comb in pklWZchiJERd:
        for variable in pklWZchiJERd[comb]:
            arWZChi.append(pklWZchiJERd[comb][variable]["value"])
            arWZChiUp.append(pklWZchiJERd[comb][variable]["erUp"])
            arWZChiDn.append(pklWZchiJERd[comb][variable]["erDown"])
            arWZChiName.append(variable)

    #print arWZChi
    #print arWZChiUp
    #print arWZChiDn
    #print arWZChiName

    c = TCanvas("canvas","",800,1200)

    x=[valueWZInt]
    ex=[errWZInt]
    y=[len(arWZChi)/2.]
    eyU=[round(len(arWZChi)/2.)+1]
    eyD=[len(arWZChi)/2.]
    ge=TGraphAsymmErrors()
    for i in range(len(x)):
        ge.SetPoint(i,x[i],y[i])
        ge.SetPointError(i,ex[i],ex[i],eyD[i],eyU[i])
    ge.SetFillStyle(3003)
    ge.SetFillColor(ROOT.kBlue)
    ge.Draw("same a2")
    
    
    gr = TGraphAsymmErrors()

    for i in range(len(arWZChi)):
        y= float(i+1)
        gr.SetPoint(i,arWZChi[i],y)
        gr.SetPointError(i,arWZChiDn[i],arWZChiUp[i],0.,0.)

    gr.SetTitle("; Scale Factor #alpha; ")
    ge.SetTitle("; Scale Factor #alpha; ")
    #gr.SetMarkerColor(4)
    gr.SetMarkerColor(kRed)
    gr.SetMarkerStyle(21)
    gr.Draw("sameP")

    line=TLine(valueWZInt,0,valueWZInt,len(arWZChi)+1)
    line.SetLineColor(kBlue)
    line.SetLineWidth(2)
    line.Draw("SAME")


    mean=np.mean(arWZChi)
    lineMEAN=TLine(mean,0,mean,len(arWZChi)+1)
    lineMEAN.SetLineColor(kRed)
    lineMEAN.SetLineWidth(2)
    lineMEAN.Draw("SAME")

    ge.GetXaxis().SetLimits(1.0,1.4)
    ge.GetYaxis().SetLabelOffset(1)
    ge.GetYaxis().SetNdivisions(len(arWZChi)*2)

    #t = TText()
    #t.SetTextAlign(32)
    #t.SetTextSize(0.035)
    #t.SetTextFont(72)
    #for i in range(len(arWZChi)):
        #t.DrawText(0.98,i+1,arWZChiName[i])

    #t = TText()
    t = TLatex()
    t.SetTextAlign(32)
    t.SetTextSize(0.02)
    #t.SetTextFont(72)
    for i in range(len(arWZChi)):
        t.DrawLatex(0.98,(i+1),nameDict[arWZChiName[i]])
        #t.DrawText(0.8,i+1,arTTChiName[i])
    
    smallT=TText()
    smallT.SetTextAngle(90)
    smallT.SetTextSize(0.025)
    smallT.SetTextFont(42)
    smallT.DrawText(0.94,4.*len(arWZChi)/4. -6,"JESu")
    smallT.DrawText(0.94,3*len(arWZChi)/4. -6,"JESd")
    smallT.DrawText(0.94,2*len(arWZChi)/4. -6,"JERu")
    smallT.DrawText(0.94,len(arWZChi)/4. -6,"JERd")
    
    cutLine=TLine(0.95,3*len(arWZChi)/4. +0.5,1.4,3*len(arWZChi)/4. +0.5)
    cutLine.SetLineColor(kBlack)
    cutLine.SetLineWidth(2)
    cutLine.Draw("SAME")
    
    cutLine2=TLine(0.95,2*len(arWZChi)/4. +0.5,1.4,2*len(arWZChi)/4. +0.5)
    cutLine2.SetLineColor(kBlack)
    cutLine2.SetLineWidth(2)
    cutLine2.Draw("SAME")
    
    cutLine2=TLine(0.95,len(arWZChi)/4. +0.5,1.4,len(arWZChi)/4. +0.5)
    cutLine2.SetLineColor(kBlack)
    cutLine2.SetLineWidth(2)
    cutLine2.Draw("SAME")

    l = ROOT.TLatex( 0.17, .95, "#font[61]{CMS} #scale[0.76]{#font[52]{Work in Progress}}")
    l.SetNDC()
    l.Draw()

    lum = ROOT.TLatex( .62, .95, "%.1f fb^{-1} (%s TeV)"%(aux.intLumi/1000., aux.Label.cmsEnergy) )
    lum.SetNDC()
    lum.Draw()
    
    l4 = ROOT.TLatex( 0.17, .9, "#scale[0.66]{#font[52]{WZ Control Region}}")
    l4.SetNDC()
    l4.Draw()
    
    
    leg=TLegend(0.18,0.45,0.4,0.55)
    leg.AddEntry(line,"#alpha from int. method","l")
    leg.AddEntry(ge,"stat. error from int. method","f")
    leg.AddEntry(gr,"#alpha from #chi^{2} method","lep")
    leg.AddEntry(lineMEAN,"mean from #chi^{2} method","l")
    leg.SetTextSize(0.02)
    #leg.SetFillStyle(0)
    leg.Draw("same")



    c.SetGridy()
    c.Update()

    c.SaveAs('plots_CR/chi/WZ_Compare_JETsyst.pdf')


def main():
    drawZZ()
    drawTT()
    drawWZ()
    drawDYLL()
    drawDYEE()
    drawDYMM()

if __name__=="__main__":
    main()
