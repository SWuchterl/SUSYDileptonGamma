import ROOT

# own style options:
divideByBinWidth = False
minimumOne = True
additionalPoissonUncertainty = False

def defaultStyle():
    st = ROOT.TStyle("defaultStyle", "Knut's owns style" )
    st.SetCanvasColor( ROOT.kWhite )
    st.SetCanvasBorderMode( 0 )
    st.SetFrameBorderMode(0)
    st.SetCanvasDefH(600)
    st.SetCanvasDefW(600)

    st.SetPadTickX( 1 )
    st.SetPadTickY( 1 )

    st.SetPadColor( ROOT.kWhite )

    # Margins:
    st.SetPadTopMargin(0.06)
    st.SetPadBottomMargin(0.12)
    st.SetPadLeftMargin(0.16)
    st.SetPadRightMargin(0.04)

    st.SetTitleFillColor( ROOT.kWhite )
    st.SetTitleBorderSize( 0 )

    st.SetTitleOffset( 1.1, "x" )
    st.SetTitleOffset( 1.6, "y" )

    st.SetStatBorderSize(1)
    st.SetStatColor(0)

    st.SetLegendBorderSize(0)
    st.SetLegendFillColor( ROOT.kWhite )
    st.SetLegendFont( st.GetLabelFont() )
    #st.SetLegendTextSize( st.GetLabelSize() ) not in current ROOT version

    st.SetOptStat(0)

    textSize = 0.05
    st.SetLabelSize( textSize, "xyz" )
    st.SetTitleSize( textSize, "xyz" )

    st.SetTextFont( st.GetLabelFont() )
    st.SetTextSize( st.GetLabelSize() )

    st.SetNdivisions( 505, "xyz" )


    st.SetPalette(56)
    st.SetNumberContours(999)

    #st.SetErrorX(0)


    st.cd()
    return st

def style2d():
    st = defaultStyle()
    st.SetPadRightMargin(0.19)
    st.SetTitleOffset( 1.35, "z" )
    return st


def setPaletteRWB():
    # Sets the current palette to red -> white -> blue
    from array import array
    steps = array('d', [0.0, 0.5, 1.0])
    red   = array('d', [1.0, 1.0, 0.0])
    green = array('d', [0.0, 1.0, 0.0])
    blue  = array('d', [0.0, 1.0, 1.0])
    ROOT.TColor.CreateGradientColorTable(len(steps), steps, red, green, blue, ROOT.gStyle.GetNumberContours())

def setPaletteBWR():
    # Sets the current palette to blue -> white -> red
    from array import array
    steps = array('d', [0.0, 0.5, 1.0])
    red   = array('d', [0.0, 1.0, 1.0])
    green = array('d', [0.0, 1.0, 0.0])
    blue  = array('d', [1.0, 1.0, 0.0])
    ROOT.TColor.CreateGradientColorTable(len(steps), steps, red, green, blue, ROOT.gStyle.GetNumberContours())

defaultStyle()

# not style, but similar
ROOT.gROOT.SetBatch()
ROOT.TH1.SetDefaultSumw2()
ROOT.gROOT.ForceStyle()
