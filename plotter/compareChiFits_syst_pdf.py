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
pklWZ = pkl.load( open( "plots_CR_wz/factors/CRWZ.pkl", "rb" ) )

pklZZchi = pkl.load( open( "plots_CR/chi/ZZ_chi_pdf.pkl", "rb" ) )

pklWZchi = pkl.load( open( "plots_CR/chi/WZ_chi_pdf.pkl", "rb" ) )

pklTTchi = pkl.load( open( "plots_CR/chi/TT_chi_pdf.pkl", "rb" ) )

pklDYchiLL = pkl.load( open( "plots_CR/chi/DY_chi_LL_pdf.pkl", "rb" ) )
pklDYchiEE = pkl.load( open( "plots_CR/chi/DY_chi_EE_pdf.pkl", "rb" ) )
pklDYchiMM = pkl.load( open( "plots_CR/chi/DY_chi_MM_pdf.pkl", "rb" ) )

pklWZchi = pkl.load( open( "plots_CR/chi/WZ_chi_pdf.pkl", "rb" ) )

#nameDict={
    #"phi1": "#phi_{leading}",
    #"phi2": "#phi_{trailing}",
    #"eta1": "#eta_{leading}",
    #"eta2": "#eta_{trailing}",
    #"pt1": "p_{T}^{leading}",
    #"pt2": "p_{T}^{trailing}",
    #"pt_g1": "p_{T}^{#gamma}",
    #"met": "p_{T}^{miss}",
    #"m_ll": "m_{ll}",
    #"m_ll2": "m_{l_{3}l_{4}}",
    #"ht": "H_{T}",
    #"n_jets": "N_{Jets}",
#}


def drawZZ():
    valueZZInt=pklZZ["LL"]["eta1"][0]
    errZZInt=pklZZ["LL"]["eta1"][1]
    arZZChi=[]
    arZZChiUp=[]
    arZZChiDn=[]
    arZZChiName=[]

    for comb in pklZZchi:
        for iPDF in pklZZchi[comb]["pt1"]:
            arZZChi.append(pklZZchi[comb]["pt1"][iPDF]["value"])
            arZZChiUp.append(pklZZchi[comb]["pt1"][iPDF]["erUp"])
            arZZChiDn.append(pklZZchi[comb]["pt1"][iPDF]["erDown"])
            arZZChiName.append(iPDF)

    print arZZChi
    print arZZChiUp
    print arZZChiDn
    print arZZChiName

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
        t.DrawLatex(1.04,(i+1),str(arZZChiName[i]))
        #t.DrawText(0.8,i+1,arTTChiName[i])
    
    #smallT=TText()
    #smallT.SetTextAngle(90)
    #smallT.SetTextSize(0.025)
    #smallT.SetTextFont(42)
    #smallT.DrawText(0.98,4.*len(arZZChi)/4. -6,"JESu")
    #smallT.DrawText(0.98,3.*len(arZZChi)/4. -6,"JESd")
    #smallT.DrawText(0.98,2.*len(arZZChi)/4. -6,"JERu")
    #smallT.DrawText(0.98,len(arZZChi)/4. -6,"JERd")
    
    #cutLine=TLine(0.9,3*len(arZZChi)/4. +0.5,1.75,3*len(arZZChi)/4. +0.5)
    #cutLine.SetLineColor(kBlack)
    #cutLine.SetLineWidth(2)
    #cutLine.Draw("SAME")
    
    #cutLine2=TLine(0.9,2*len(arZZChi)/4. +0.5,1.75,2*len(arZZChi)/4. +0.5)
    #cutLine2.SetLineColor(kBlack)
    #cutLine2.SetLineWidth(2)
    #cutLine2.Draw("SAME")
    
    #cutLine2=TLine(0.9,len(arZZChi)/4. +0.5,1.75,len(arZZChi)/4. +0.5)
    #cutLine2.SetLineColor(kBlack)
    #cutLine2.SetLineWidth(2)
    #cutLine2.Draw("SAME")

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

    c.SaveAs('plots_CR/chi/ZZ_Compare_PDFsyst.pdf')
    
def drawDYLL():
    valueDYInt=pklDY["LL"]["eta1"][0]
    errDYInt=pklDY["LL"]["eta1"][1]
    arDYChi=[]
    arDYChiUp=[]
    arDYChiDn=[]
    arDYChiName=[]

    for comb in pklDYchiLL:
        for iPDF in pklDYchiLL[comb]["pt1"]:
            arDYChi.append(pklDYchiLL[comb]["pt1"][iPDF]["value"])
            arDYChiUp.append(pklDYchiLL[comb]["pt1"][iPDF]["erUp"])
            arDYChiDn.append(pklDYchiLL[comb]["pt1"][iPDF]["erDown"])
            arDYChiName.append(iPDF)

    print arDYChi
    print arDYChiUp
    print arDYChiDn
    print arDYChiName

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
        t.DrawLatex(1.047,(i+1),str(arDYChiName[i]))
        #t.DrawText(0.8,i+1,arTTChiName[i])
    
    #smallT=TText()
    #smallT.SetTextAngle(90)
    #smallT.SetTextSize(0.025)
    #smallT.SetTextFont(42)
    #smallT.DrawText(1.04,4.*len(arDYChi)/4. -6,"JESu")
    #smallT.DrawText(1.04,3*len(arDYChi)/4. -6,"JESd")
    #smallT.DrawText(1.04,2*len(arDYChi)/4. -6,"JERu")
    #smallT.DrawText(1.04,len(arDYChi)/4. -6,"JERd")
    
    #cutLine=TLine(1.04,3*len(arDYChi)/4. +0.5,1.12,3*len(arDYChi)/4. +0.5)
    #cutLine.SetLineColor(kBlack)
    #cutLine.SetLineWidth(2)
    #cutLine.Draw("SAME")
    
    #cutLine2=TLine(1.04,2*len(arDYChi)/4. +0.5,1.12,2*len(arDYChi)/4. +0.5)
    #cutLine2.SetLineColor(kBlack)
    #cutLine2.SetLineWidth(2)
    #cutLine2.Draw("SAME")
    
    #cutLine2=TLine(1.04,len(arDYChi)/4. +0.5,1.12,len(arDYChi)/4. +0.5)
    #cutLine2.SetLineColor(kBlack)
    #cutLine2.SetLineWidth(2)
    #cutLine2.Draw("SAME")

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

    c.SaveAs('plots_CR/chi/DY_CompareLL_PDFsyst.pdf')
def drawDYEE():
    valueDYInt=pklDY["EE"]["eta1"][0]
    errDYInt=pklDY["EE"]["eta1"][1]
    arDYChi=[]
    arDYChiUp=[]
    arDYChiDn=[]
    arDYChiName=[]

    for comb in pklDYchiEE:
        for iPDF in pklDYchiEE[comb]["pt1"]:
            arDYChi.append(pklDYchiEE[comb]["pt1"][iPDF]["value"])
            arDYChiUp.append(pklDYchiEE[comb]["pt1"][iPDF]["erUp"])
            arDYChiDn.append(pklDYchiEE[comb]["pt1"][iPDF]["erDown"])
            arDYChiName.append(iPDF)


    print arDYChi
    print arDYChiUp
    print arDYChiDn
    print arDYChiName

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
        t.DrawLatex(1.015,(i+1),str(arDYChiName[i]))
        #t.DrawText(0.8,i+1,arTTChiName[i])
    
    #smallT=TText()
    #smallT.SetTextAngle(90)
    #smallT.SetTextSize(0.025)
    #smallT.SetTextFont(42)
    #smallT.DrawText(1.005,4.*len(arDYChi)/4. -6,"JESu")
    #smallT.DrawText(1.005,3*len(arDYChi)/4. -6,"JESd")
    #smallT.DrawText(1.005,2*len(arDYChi)/4. -6,"JERu")
    #smallT.DrawText(1.005,len(arDYChi)/4. -6,"JERd")
    
    #cutLine=TLine(1.005,3*len(arDYChi)/4. +0.5,1.75,3*len(arDYChi)/4. +0.5)
    #cutLine.SetLineColor(kBlack)
    #cutLine.SetLineWidth(2)
    #cutLine.Draw("SAME")
    
    #cutLine2=TLine(1.005,2*len(arDYChi)/4. +0.5,1.75,2*len(arDYChi)/4. +0.5)
    #cutLine2.SetLineColor(kBlack)
    #cutLine2.SetLineWidth(2)
    #cutLine2.Draw("SAME")
    
    #cutLine2=TLine(1.005,len(arDYChi)/4. +0.5,1.75,len(arDYChi)/4. +0.5)
    #cutLine2.SetLineColor(kBlack)
    #cutLine2.SetLineWidth(2)
    #cutLine2.Draw("SAME")

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

    c.SaveAs('plots_CR/chi/DY_CompareEE_PDFsyst.pdf')
def drawDYMM():
    valueDYInt=pklDY["MM"]["eta1"][0]
    errDYInt=pklDY["MM"]["eta1"][1]
    arDYChi=[]
    arDYChiUp=[]
    arDYChiDn=[]
    arDYChiName=[]

    for comb in pklDYchiMM:
        for iPDF in pklDYchiMM[comb]["pt1"]:
            arDYChi.append(pklDYchiMM[comb]["pt1"][iPDF]["value"])
            arDYChiUp.append(pklDYchiMM[comb]["pt1"][iPDF]["erUp"])
            arDYChiDn.append(pklDYchiMM[comb]["pt1"][iPDF]["erDown"])
            arDYChiName.append(iPDF)

    print arDYChi
    print arDYChiUp
    print arDYChiDn
    print arDYChiName

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
        t.DrawLatex(1.047,(i+1),str(arDYChiName[i]))
        #t.DrawText(0.8,i+1,arTTChiName[i])
    
    #smallT=TText()
    #smallT.SetTextAngle(90)
    #smallT.SetTextSize(0.025)
    #smallT.SetTextFont(42)
    #smallT.DrawText(1.04,4.*len(arDYChi)/4. -6,"NBins * 2")
    #smallT.DrawText(1.04,3*len(arDYChi)/4. -6,"normal binning")
    #smallT.DrawText(1.04,2*len(arDYChi)/4. -6,"NBins / 2")
    #smallT.DrawText(1.04,len(arDYChi)/4. -6,"NBins / 2")
    
    #cutLine=TLine(1.04,3*len(arDYChi)/4. +0.5,1.75,3*len(arDYChi)/4. +0.5)
    #cutLine.SetLineColor(kBlack)
    #cutLine.SetLineWidth(2)
    #cutLine.Draw("SAME")
    
    #cutLine2=TLine(1.04,2*len(arDYChi)/4. +0.5,1.75,2*len(arDYChi)/4. +0.5)
    #cutLine2.SetLineColor(kBlack)
    #cutLine2.SetLineWidth(2)
    #cutLine2.Draw("SAME")
    
    #cutLine2=TLine(1.04,len(arDYChi)/4. +0.5,1.75,len(arDYChi)/4. +0.5)
    #cutLine2.SetLineColor(kBlack)
    #cutLine2.SetLineWidth(2)
    #cutLine2.Draw("SAME")

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

    c.SaveAs('plots_CR/chi/DY_CompareMM_PDFsyst.pdf')

def drawTT():
    valueTTInt=pklTT["EM"]["eta1"][0]
    errTTInt=pklTT["EM"]["eta1"][1]
    arTTChi=[]
    arTTChiUp=[]
    arTTChiDn=[]
    arTTChiName=[]


    for comb in pklTTchi:
        for iPDF in pklTTchi[comb]["pt1"]:
            arTTChi.append(pklTTchi[comb]["pt1"][iPDF]["value"])
            arTTChiUp.append(pklTTchi[comb]["pt1"][iPDF]["erUp"])
            arTTChiDn.append(pklTTchi[comb]["pt1"][iPDF]["erDown"])
            arTTChiName.append(iPDF)

    #print arTTChi
    #print arTTChiUp
    #print arTTChiDn
    #print arTTChiName

    c = TCanvas("canvas","",800,1200)

    x=[valueTTInt]
    ex=[errTTInt]
    y=[len(arTTChi)/2.]
    eyU=[round(len(arTTChi)/2.)+1]
    eyD=[len(arTTChi)/2.]
    ge=TGraphAsymmErrors()
    for i in range(len(x)):
        ge.SetPoint(i,x[i],y[i])
        ge.SetPointError(i,ex[i],ex[i],eyD[i],eyU[i])
    ge.SetFillStyle(3003)
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
    ge.GetYaxis().SetNdivisions(len(arTTChi)*2)

    #t = TText()
    t = TLatex()
    t.SetTextAlign(32)
    t.SetTextSize(0.02)
    #t.SetTextFont(72)
    for i in range(len(arTTChi)):
        t.DrawLatex(0.79,(i+1),str(arTTChiName[i]))
        #t.DrawText(0.8,i+1,arTTChiName[i])
    
    #smallT=TText()
    #smallT.SetTextAngle(90)
    #smallT.SetTextSize(0.025)
    #smallT.SetTextFont(42)
    #smallT.DrawText(0.76,4.*len(arTTChi)/4. -6,"JESu")
    #smallT.DrawText(0.76,3*len(arTTChi)/4. -6,"JESd")
    #smallT.DrawText(0.76,2*len(arTTChi)/4. -6,"JERu")
    #smallT.DrawText(0.76,len(arTTChi)/4. -6,"JERd")
    
    #cutLine=TLine(0.75,3*len(arTTChi)/4. +0.5,1.1,3*len(arTTChi)/4. +0.5)
    #cutLine.SetLineColor(kBlack)
    #cutLine.SetLineWidth(2)
    #cutLine.Draw("SAME")
    
    #cutLine2=TLine(0.75,2*len(arTTChi)/4. +0.5,1.1,2*len(arTTChi)/4. +0.5)
    #cutLine2.SetLineColor(kBlack)
    #cutLine2.SetLineWidth(2)
    #cutLine2.Draw("SAME")
    
    #cutLine2=TLine(0.75,len(arTTChi)/4. +0.5,1.1,len(arTTChi)/4. +0.5)
    #cutLine2.SetLineColor(kBlack)
    #cutLine2.SetLineWidth(2)
    #cutLine2.Draw("SAME")

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

    c.SaveAs('plots_CR/chi/TT_Compare_PDFsyst.pdf')
    
def drawWZ():
    valueWZInt=pklWZ["LL"]["eta1"][0]
    errWZInt=pklWZ["LL"]["eta1"][1]
    arWZChi=[]
    arWZChiUp=[]
    arWZChiDn=[]
    arWZChiName=[]

    for comb in pklWZchi:
        for iPDF in pklWZchi[comb]["pt1"]:
            arWZChi.append(pklWZchi[comb]["pt1"][iPDF]["value"])
            arWZChiUp.append(pklWZchi[comb]["pt1"][iPDF]["erUp"])
            arWZChiDn.append(pklWZchi[comb]["pt1"][iPDF]["erDown"])
            arWZChiName.append(iPDF)


    print arWZChi
    print arWZChiUp
    print arWZChiDn
    print arWZChiName

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
        t.DrawLatex(0.98,(i+1),str(arWZChiName[i]))
        #t.DrawText(0.8,i+1,arTTChiName[i])
    
    #smallT=TText()
    #smallT.SetTextAngle(90)
    #smallT.SetTextSize(0.025)
    #smallT.SetTextFont(42)
    #smallT.DrawText(0.94,4.*len(arWZChi)/4. -6,"JESu")
    #smallT.DrawText(0.94,3*len(arWZChi)/4. -6,"JESd")
    #smallT.DrawText(0.94,2*len(arWZChi)/4. -6,"JERu")
    #smallT.DrawText(0.94,len(arWZChi)/4. -6,"JERd")
    
    #cutLine=TLine(0.95,3*len(arWZChi)/4. +0.5,1.4,3*len(arWZChi)/4. +0.5)
    #cutLine.SetLineColor(kBlack)
    #cutLine.SetLineWidth(2)
    #cutLine.Draw("SAME")
    
    #cutLine2=TLine(0.95,2*len(arWZChi)/4. +0.5,1.4,2*len(arWZChi)/4. +0.5)
    #cutLine2.SetLineColor(kBlack)
    #cutLine2.SetLineWidth(2)
    #cutLine2.Draw("SAME")
    
    #cutLine2=TLine(0.95,len(arWZChi)/4. +0.5,1.4,len(arWZChi)/4. +0.5)
    #cutLine2.SetLineColor(kBlack)
    #cutLine2.SetLineWidth(2)
    #cutLine2.Draw("SAME")

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

    c.SaveAs('plots_CR/chi/WZ_Compare_PDFsyst.pdf')


def main():
    drawZZ()
    drawTT()
    drawWZ()
    drawDYLL()
    drawDYEE()
    drawDYMM()

if __name__=="__main__":
    main()
