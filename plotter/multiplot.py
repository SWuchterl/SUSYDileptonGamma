import ROOT
import auxiliary as aux

class Multiplot:
    def __init__( self ):
        self.hists = []
        self.histsToStack = []

        # todo: impmlement setter and getter
        self.minimum = None
        self.maximum = None

        self.leg = ROOT.TLegend(.56,.59,.94,.915)
        self.leg.SetFillColor( ROOT.kWhite )
        self.leg.SetFillStyle(0)

    def add( self, h, label="" ):
        h.SetName(label)
        self.hists.append( h )

    def addStack( self, h, label="" ):
        h.SetName( label )
        self.histsToStack.append( h )

    def getMinimum( self ):
        return min( [ h.GetMinimum(0) for h in self.hists+self.histsToStack if not isinstance( h, ROOT.THStack ) and not isinstance( h, ROOT.TGraph ) ] )

    def getMaximum( self ):
        return max( [ h.GetMaximum() for h in self.hists ] )

    def getStack(self):
        stacks = [h for h in self.hists if isinstance(h, ROOT.THStack)]
        return stacks[0] if stacks else None

    def stackHists( self ):
        if not self.histsToStack:
            return
        stack = ROOT.THStack()
        stack.SetTitle( ";%s;%s"%(self.histsToStack[0].GetXaxis().GetTitle(),self.histsToStack[0].GetYaxis().GetTitle()) )
        stack.drawOption_ = "hist"
        for h in self.histsToStack:
            h.SetFillColor( h.GetLineColor() )
            h.SetLineColor( ROOT.kBlack )
            stack.Add( h )

        self.hists.append( stack )

    def sortStackByIntegral( self ):
        self.histsToStack = sorted( self.histsToStack, key=lambda x: x.Integral(0,-1) )


    def Draw( self ):
        if not self.hists and not self.histsToStack:
            return False
        self.stackHists()

        minimum = 1e-3
        maximum = 1.1*self.getMaximum()

        if self.maximum != None:
            maximum = self.maximum
        if self.minimum != None:
            minimum = self.minimum

        # Fill legend:
        # Data first
        for h in self.hists:
            if isinstance( h, ROOT.THStack ): continue
            if not hasattr( h, "drawOption_" ): h.drawOption_ = ""
            if aux.dataLikeName(h.GetName()):
                self.leg.AddEntry( h, h.GetName(), "pel" )

        # Stacked histograms
        for h in self.histsToStack[-1::-1]:
            h.SetLineColor(0)
            self.leg.AddEntry( h, h.GetName(), "f" )

        # Other histograms
        for h in self.hists:
            if not h.GetName(): continue
            if isinstance( h, ROOT.THStack ): continue
            if aux.dataLikeName(h.GetName()): continue

            if "p" in h.drawOption_:
                self.leg.AddEntry( h, h.GetName(), "ep" )
            elif "e2" in h.drawOption_:
                h.SetLineColor(0)
                self.leg.AddEntry( h, h.GetName(), "epf" )
            else:
                self.leg.AddEntry( h, h.GetName(), "l" )

        # change the order for drawing
        self.hists.reverse()
        for ih, h in enumerate(self.hists):
            import style
            if aux.dataLikeName(h.GetName()) and aux.integerContent(h, True) and style.divideByBinWidth:
                gr = ROOT.TGraphAsymmErrors(h)
                aux.saveStuff.append(gr)
                aux.drawOpt(gr, "data")
                for p in range(gr.GetN()):
                    bw = h.GetBinWidth(p+1)
                    entries = int(round(gr.GetY()[p]*bw))
                    edn, eup = aux.getPoissonUnc(entries)
                    gr.SetPointEYlow(p, edn/bw)
                    gr.SetPointEYhigh(p, eup/bw)
                gr.Draw(gr.drawOption_+"same")
            else:
                if not ih:
                    h.SetMinimum(minimum)
                    h.SetMaximum(maximum)
                else:
                    h.drawOption_ += "same"
                h.Draw(h.drawOption_)

        self.leg.Draw()

        return True

    def draw(self): # simple alias
        self.Draw()
