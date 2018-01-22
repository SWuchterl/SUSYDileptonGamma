import ROOT
# Usage of the corporate design colors of RWTH Aachen University
# RGB values copied from here:
# http://www9.rwth-aachen.de/global/show_document.asp?id=aaaaaaaaaadrbpl#page=27

# definition of the colors:
colors = {
    "blue": (0,84,159),
    "blue75": (64,127,183),
    "blue50": (142,186,229),
    "blue25": (199,221,242),
    "blue10": (232,241,250),
    "black": (0,0,0),
    "black75": (100,101,103),
    "black50": (156,158,159),
    "black25": (207,209,210),
    "black10": (236,237,237),
    "magenta": (227,0,102),
    "magenta75": (233,96,136),
    "magenta50": (241,158,177),
    "magenta25": (249,210,218),
    "magenta10": (253,238,240),
    "yellow": (255,237,0),
    "yellow75": (255,240,85),
    "yellow50": (255,245,155),
    "yellow25": (255,250,209),
    "yellow10": (255,253,238),
    "petrol": (0,97,101),
    "petrol75": (45,127,131),
    "petrol50": (125,164,167),
    "petrol25": (191,208,209),
    "petrol10": (230,236,236),
    "cyan": (0,152,161),
    "cyan75": (0,177,183),
    "cyan50": (137,204,207),
    "cyan25": (202,231,231),
    "cyan10": (235,246,246),
    "green": (87,171,39),
    "green75": (141,192,96),
    "green50": (184,214,152),
    "green25": (221,235,206),
    "green10": (242,247,236),
    "orange": (246,168,0),
    "orange75": (250,190,80),
    "orange50": (253,212,143),
    "orange25": (254,234,201),
    "orange10": (255,247,234),
    "red": (204,7,30),
    "red75": (216,92,65),
    "red50": (230,150,121),
    "red25": (243,205,187),
    "red10": (250,235,227),
    "bordeaux": (161,16,53),
    "bordeaux75": (182,82,86),
    "bordeaux50": (205,139,135),
    "bordeaux25": (229,197,192),
    "bordeaux10": (245,232,229),
    "violett": (97,33,88),
    "violett75": (131,78,117),
    "violett50": (168,133,158),
    "violett25": (210,192,205),
    "violett10": (237,229,234),
    "lila": (122,111,172),
    "lila75": (155,145,193),
    "lila50": (188,181,215),
    "lila25": (222,218,235),
    "lila10": (242,240,247),
    "mayGreen": (189,205,0),
    "mayGreen75": (208,217,92),
    "mayGreen50": (224,230,154),
    "mayGreen25": (240,243,208),
    "mayGreen10": (249,250,237),
    # From http://colorbrewer2.org/#type=diverging&scheme=RdYlBu&n=5
    "myRed": (215,25,28),
    "myOrange": (253,174,97),
    "myYellow": (255,255,191),
    "myLightBlue": (171,217,233),
    "myBlue": (44,123,182),
}

# This class will contain all colors as members
class rwth:
    colors = [] # colors are stored here to avoid garbage collector

startColor = 5001 # Starting point for ROOT internal colors naming scheme

for i, (name,(r,g,b)) in enumerate(colors.iteritems()):
    exec "rwth.{} = {:d}".format(name, startColor+i )
    rwth.colors.append( ROOT.TColor(startColor+i, r/255., g/255., b/255.) )

# Usage: histogram.SetLineColor( rwth.blue50 )

if __name__ == "__main__":
    import re
    unmodifiedKeys = [ i for i in colors if not re.match(".*\d\d",i)]
    t = ROOT.TLatex()
    c = ROOT.TCanvas("","",100,100)
    for ik, k in enumerate(unmodifiedKeys):
        exec "c1 = rwth.%s%s"%(k,"")
        exec "c2 = rwth.%s%s"%(k,"75")
        exec "c3 = rwth.%s%s"%(k,"50")
        exec "c4 = rwth.%s%s"%(k,"25")
        y = .1+0.8*ik/len(unmodifiedKeys)
        t.DrawLatex(.1,y,k)
        t.DrawLatex(.3,y,"#lower[.3]{#scale[2]{#color[%s]{#bullet} #color[%s]{#bullet} #color[%s]{#bullet} #color[%s]{#bullet}}}"%(c1,c2,c3,c4))
    ROOT.gPad.SaveAs("plots/colorsOvierview.pdf")
