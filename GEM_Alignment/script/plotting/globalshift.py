import ROOT, os
from make_plot import *

#f = ROOT.TFile("../condor/run342218_jun30/out_run342218_jun30.root")
#f = ROOT.TFile("/afs/cern.ch/user/t/toakhter/public/GEMAl2018Cv4.root")
#f = ROOT.TFile("/afs/cern.ch/work/d/daebi/analyser/CMSSW_11_3_0_pre5/src/GEMCSCBendingAnalyzer/MuonAnalyser/condor/MWGR4_condor/run342218_jul12/out_run342218_jul14.root")
#f = ROOT.TFile("out_ZeroBias_testbeam_nov12.root")
f = ROOT.TFile("out_ME11ana_CRAFT2022_ExpressCosmics_ME11Iter2_Align.root")

#event = f.Get("analyser/MuonData")
#event_ME11Seg = f.Get("analyser/ME11Seg_Prop")
#event_CSC = f.Get("analyser/CSC_Prop")
event_ME11ana = f.Get("ME11ana/Inner_Prop")


RdPhi_Low = -1
RdPhi_High = 1

h1 = plot1Dprofile("RdPhi:prop_location[3]", "ME11 RdPhi Profile CRAFT2022 R-1", event_ME11ana, [36, 1, 37, 10, RdPhi_Low, RdPhi_High], "Chamber", "RdPhi", "abs(RdPhi) < 20 && prop_location[0] == -1", "RdPhi", True)
c1 = ROOT.TCanvas("", "", 800, 600)
hist = ROOT.TH1D("h1", "ME11 RdPhi Profile CRAFT2022 R-1", 36, -(3.14159265)/36., 2*3.14159265 - (3.14159265)/36.)
for i in range(1, 37):
  hist.SetBinContent(i, h1.GetBinContent(i))

hist.GetYaxis().SetRangeUser(RdPhi_Low, RdPhi_High)
hist.GetXaxis().SetTitle("Phi")
hist.GetYaxis().SetTitle("RdPhi [cm]")
hist.Sumw2()
hist.Draw("Hist")
f1 = ROOT.TF1("f1", "[0]*sin(x) + [1]*cos(x) + [2]", 0, 2*3.14159265)
f1.SetParameters(.1, .1, 0)
hist.Fit("f1")
f1.Draw("same")
hist.SetStats(False)
latex = ROOT.TLatex()
latex.SetTextAlign(12)
latex.SetTextSize(0.04)
latex.SetTextFont(61)
latex.DrawLatex(2, 0.6, "{zer} sin(x) + {one} cos(x) + {two}".format(zer = round(f1.GetParameter(0), 4), one = round(f1.GetParameter(1), 4), two = round(f1.GetParameter(2), 4)))
c1.SaveAs("ME11_R-1_GlobalShift.png")



h1 = plot1Dprofile("RdPhi:prop_location[3]", "ME11 RdPhi Profile CRAFT2022 R1", event_ME11ana, [36, 1, 37, 10, RdPhi_Low, RdPhi_High], "Chamber", "RdPhi", "abs(RdPhi) < 20 && prop_location[0] == 1", "RdPhi", True)

c1 = ROOT.TCanvas("", "", 800, 600)
hist = ROOT.TH1D("h1", "ME11 RdPhi Profile CRAFT2022 R1", 36, -(3.14159265)/36., 2*3.14159265 - (3.14159265)/36.)
for i in range(1, 37):
  hist.SetBinContent(i, h1.GetBinContent(i))

hist.GetYaxis().SetRangeUser(RdPhi_Low, RdPhi_High)
hist.GetXaxis().SetTitle("Phi")
hist.GetYaxis().SetTitle("RdPhi [cm]")
hist.Sumw2()
hist.Draw("Hist")
f1 = ROOT.TF1("f1", "[0]*sin(x) + [1]*cos(x) + [2]", 0, 2*3.14159265)
f1.SetParameters(1, 1, 0)
hist.Fit("f1")
f1.Draw("same")
hist.SetStats(False)
latex = ROOT.TLatex()
latex.SetTextAlign(12)
latex.SetTextSize(0.04)
latex.SetTextFont(61)
latex.DrawLatex(2, 0.6, "{zer} sin(x) + {one} cos(x) + {two}".format(zer = round(f1.GetParameter(0), 4), one = round(f1.GetParameter(1), 4), two = round(f1.GetParameter(2), 4)))
c1.SaveAs("ME11_R1_GlobalShift.png")
