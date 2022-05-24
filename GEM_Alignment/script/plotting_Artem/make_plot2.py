import ROOT, array, math

imyFile = ROOT.TFile.Open("DTPlots/DT_tbma_cosmics.root", "READ")
#myFile = ROOT.TFile("CRUZET_march15.ROOT")
myTree = myFile.Get("DT_tbma/Inner_Prop")

# myHist = ROOT.TH1D("rdphi", "rdphi CRUZET R1", 100, -2, 2)
# myHist.Sumw2()

# which_residual = input("enter type of the residual: ")
# which_station = input("enter station: ")

pi = math.pi

#print("There are ", myTree.GetEntries(), " entries")
#for i in myTree:
#    if (i.RdPhi_Corrected < 100):
#         myHist.Fill(i.RdPhi_Corrected)

# myTree.Project("rdphi", "RdPhi_Corrected", "RdPhi_Corrected < 100 && prop_location[0] == 1")
canvas = ROOT.TCanvas("canvas","", 1200, 800)
canvas.SetWindowSize(500, 500)
ablue = array.array("d", [1,1,0])
ared = array.array("d", [0,1,1])
agreen = array.array("d", [0,1,0])
astop = array.array("d", [0,.5,1])
myPalette = []
fi = ROOT.TColor.CreateGradientColorTable(3, astop, ared, agreen, ablue, 100)
for x in range(100):
  myPalette.append(fi+x)
ROOT.gStyle.SetPalette(100, array.array("i", myPalette))
# my2dprofile.SetStats(0)
canvas.SetRightMargin(0.15)
canvas.SetWindowSize(500, 500)
# myHist.GetYaxis().SetTitle("events")
# myHist.GetXaxis().SetTitle("rdphi")
# myHist.Draw("h")
# canvas.SaveAs("AL2_hist_front.png")

#2d profile
for which_residual in ["dx", "dy"]:
    for which_station in [1, 2, 3, 4]:
        my2dprofile = ROOT.TProfile2D(which_residual + " vs global phi vs global z, staion" + str(which_station), which_residual + " vs global phi vs global z, staion" + str(which_station),  200, -pi, pi, 200, -700, 700, -5, 5)
        my2dprofile.GetXaxis().SetTitle("global z")
        my2dprofile.GetYaxis().SetTitle("global phi")
        my2dprofile.GetZaxis().SetTitle(which_residual)
        my2dprofile.Sumw2()
        myTree.Project(which_residual + " vs global phi vs global z, staion" + str(which_station), which_residual + ":prop_GP[2]:prop_global_phi", which_residual + " < 100 && prop_location[1] == " + str(which_station) )
        my2dprofile.Draw("colz")
        canvas.SaveAs("DTPlots/stations/" + str(which_residual) + "_Station" + str(which_station) + ".png")

# myHist.SetDirectory(0)
myFile.Close()
