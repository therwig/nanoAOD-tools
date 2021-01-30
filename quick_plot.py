import ROOT
ROOT.gStyle.SetOptStat(0)
colz=[ROOT.kBlack,
      ROOT.TColor.GetColor("#00B5E2"), # light blue
      ROOT.TColor.GetColor("#FED141"), # yellow (darker EAAA00)
      ROOT.TColor.GetColor("#AF272F"), # red
      ROOT.TColor.GetColor("#78BE21"), # green 78BE21 4C8C2C
      ROOT.TColor.GetColor("#F68D2E"), # orange F68D2E CB6014
      ROOT.TColor.GetColor("#004C97"), # blue
]

def defaultLegend( coords=(.60,.72, 0.90,.95) ):
    leg = ROOT.TLegend(*coords)
    leg.SetTextFont(42)
    leg.SetHeader("")
    leg.SetNColumns(1)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    return leg

f = ROOT.TFile("histOut.root","read")
c = ROOT.TCanvas("c","",400,300)
plotDir="plots"
ROOT.gSystem.mkdir(plotDir)

#
# plot two histograms on the same canvas
#
h1 = f.Get("plots/h_ele_pt")
h2 = f.Get("plots/h_ele_pt_match")

h1.SetLineWidth(2)
h2.SetLineWidth(2)
h1.SetLineColor(colz[0])
h2.SetLineColor(colz[1])
h1.Draw()
h2.Draw("same")

leg = defaultLegend()
leg.AddEntry(h1, "Reco electrons","l")
leg.AddEntry(h2, "Matched reco electrons","l")
leg.Draw()

c.SaveAs(plotDir+"/reco_pt.pdf")



#
# plot 2d histograms
#
h = f.Get("plots/h_ele_pt_eta")
h.Draw("colz")
c.SaveAs(plotDir+"/reco_pt_eta.pdf")



#
# plot efficiency
#
c.Clear()
eff = f.Get("plots/gen_ele_efficiency")
eff.SetLineWidth(2)
eff.SetLineColor(colz[0])
eff.Draw("AP")

leg = defaultLegend( coords=(.12,.80, 0.40,.90) )
leg.AddEntry(eff, "Reco ele efficiency","pel")
leg.Draw()

c.SaveAs(plotDir+"/eff_reco_match.pdf")
