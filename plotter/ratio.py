import ROOT
import auxiliary as aux
from math import sqrt
import auxiliary as aux
import style


def clearXaxisCurrentPad():
    # Delete label and title of all histograms in the current pad
    for ding in ROOT.gPad.GetListOfPrimitives():
        if isinstance( ding, ROOT.TH1 ) or isinstance( ding, ROOT.THStack ):
            xaxis = ding.GetXaxis()
            xaxis.SetLabelSize(0)
            xaxis.SetLabelColor(0)
            xaxis.SetLabelOffset(1000)
            xaxis.SetTitle("")
            xaxis.SetTitleColor(0)
            xaxis.SetTitleSize(0)

def createBottomPad( r=.2 ):
    # r: the ratio in which the pad is splitted
    ROOT.gPad.SetBottomMargin( r + (1-r)*ROOT.gPad.GetBottomMargin() - r*ROOT.gPad.GetTopMargin() )
    rPad = ROOT.TPad( "rPad", "ratio", 0, 0, 1, 1 )
    rPad.SetTopMargin( (1-r) - (1-r)*rPad.GetBottomMargin() + r*rPad.GetTopMargin() )
    rPad.SetFillStyle(3955)
    rPad.Draw()
    rPad.cd()
    rPad.SetLogy(0)
    ROOT.SetOwnership(rPad, False)
    return rPad

def convertToTH1( profile ):
    h = profile.ProjectionX( aux.randomName() )
    h.SetLineColor( profile.GetLineColor() )
    return h

class Ratio:
    def __init__( self, title, numerator, denominator, sysHisto=None ):

        # convcert TProfiles to histograms
        if isinstance( numerator, ROOT.TProfile ): numerator = convertToTH1( numerator )
        if isinstance( denominator, ROOT.TProfile ): denominator = convertToTH1( denominator )

        self.title = title
        self.numerator = numerator
        self.denominator = denominator.Clone( aux.randomName() )
        self.sysHisto = sysHisto
        self.ratio = numerator.Clone( aux.randomName() )
        self.ratioStat = denominator.Clone( aux.randomName() )
        self.ratioSys = sysHisto.Clone( aux.randomName() ) if sysHisto else denominator.Clone(aux.randomName())
        self.totalUncert = denominator.Clone( aux.randomName() )
        self.allowUnsymmetricYaxis = False

        self.ratio.drawOption_ = "e0"

        # Set ratio properties
        for hist in [ self.ratio, self.ratioSys, self.ratioStat, self.totalUncert ]:
            hist.GetYaxis().SetNdivisions(2, 5, 2)
            hist.SetTitleOffset(1.2, "Y")
            hist.SetYTitle( self.title )

        aux.drawOpt(self.totalUncert, "totUnc")
        aux.drawOpt(self.ratioSys, "sysUnc")
        aux.drawOpt(self.ratioStat, "statUnc")
        #aux.drawOpt(self.ratioStat, "totUnc")

    def calculateRatio( self ):
        for bin in range(self.denominator.GetNbinsX()+2): self.denominator.SetBinError(bin,0)
        self.ratio.Divide( self.denominator )
        self.ratioGraph = ROOT.TGraphAsymmErrors(self.ratio)
        for bin in range(self.ratio.GetNbinsX()+1):
            den = self.denominator.GetBinContent(bin)
            if den:
                self.ratioGraph.SetPointEYhigh(bin-1, self.numerator.GetBinErrorUp(bin)/den)
                self.ratioGraph.SetPointEYlow(bin-1, self.numerator.GetBinErrorLow(bin)/den)
                if aux.integerContent(self.numerator, True) and style.divideByBinWidth:
                    bw = self.numerator.GetBinWidth(bin)
                    entries = int(round(self.numerator.GetBinContent(bin)*bw))
                    edn, eup = aux.getPoissonUnc(entries)
                    self.ratioGraph.SetPointEYhigh(bin-1, eup/den/bw)
                    self.ratioGraph.SetPointEYlow(bin-1, edn/den/bw)

        self.ratioStat.Divide( self.denominator )
        if self.sysHisto:
            self.ratioSys.Divide( self.denominator )
            for bin in range(self.denominator.GetNbinsX()+2):
                self.totalUncert.SetBinContent(bin, 1)
                self.totalUncert.SetBinError(bin, sqrt(self.ratioSys.GetBinError(bin)**2+self.ratioStat.GetBinError(bin)**2))


    def getYrange( self ):
        # If no minimum or maximum is specified, choose a minimum from 0 to .5
        # and a maximum from 1.5 to 50
        yMin = 0
        minimum = self.ratio.GetBinContent(self.ratio.GetMaximumBin())
        for bin in range( self.ratio.GetNbinsX()+2 ):
            minInBin = self.ratio.GetBinContent(bin)
            if minInBin < minimum and minInBin > 0:
                minimum = minInBin
        yMin = minimum*.95

        from math import ceil
        yMax = min( max(1.5, ceil(self.ratio.GetBinContent(self.ratio.GetMaximumBin()))), 50 )
        yMax = self.ratio.GetBinContent(self.ratio.GetMaximumBin())*1.05


        yValues = [ self.ratio.GetBinContent(bin) for bin in range( self.ratio.GetNbinsX()+2 ) ]
        yValues = filter( lambda a: a != 0, yValues )
        if self.allowUnsymmetricYaxis:
            return yMin, yMax
        else:
            yValuesAbsDiff = [ abs(x-1) for x in yValues ]
            yValuesAbsDiff.sort()
            maxYDiff = yValuesAbsDiff[-1]
            y = 0.5
            if maxYDiff < 0.05:
                y = 0.05
            if maxYDiff < 0.01:
                y = 0.01

            return 1-y, 1+y


    def draw( self, yMin=.9, yMax=1.1, stack=None, onlyTotal=False ):
        self.calculateRatio()

        #yMin, yMax = self.getYrange()
        for hist in [ self.ratio, self.ratioSys, self.ratioStat, self.totalUncert ]:
            hist.SetMinimum( yMin )
            hist.SetMaximum( yMax )

        clearXaxisCurrentPad()
        p = createBottomPad()

        if stack: aux.drawContributions(stack, yMin, yMax, self.title)

        self.ratioStat.Draw("e x0" + "same" if stack else "")
        #self.ratioStat.Draw("same e2")
        if self.sysHisto:
            if not onlyTotal: self.ratioSys.Draw("same e2")
            self.totalUncert.Draw("same e2")
        #self.ratio.Draw("same "+self.ratio.drawOption_)
        self.ratioGraph.Draw("same pz0")

        if yMin < 1 and yMax > 1:
            oneLine = ROOT.TLine()
            oneLine.SetLineStyle(2)
            axis = self.ratio.GetXaxis()
            oneLine.DrawLine( axis.GetBinLowEdge(axis.GetFirst()), 1.0, axis.GetBinLowEdge(1+axis.GetLast()), 1.0 )
