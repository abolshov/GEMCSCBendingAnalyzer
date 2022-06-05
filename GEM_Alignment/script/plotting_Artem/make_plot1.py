import ROOT, array, math

myFile = ROOT.TFile.Open("DTPlots/DT_tbma_collisionMC_allRes.root", "READ")
#myFile = ROOT.TFile("CRUZET_march15.ROOT")
myTree = myFile.Get("DT_tbma/Inner_Prop_ChamberLevel")

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

canvas.SetRightMargin(0.15)
canvas.SetWindowSize(500, 500)
# myHist.GetYaxis().SetTitle("events")
# myHist.GetXaxis().SetTitle("rdphi")
# myHist.Draw("h")
# canvas.SaveAs("AL2_hist_front.png")

#2d profile

for which_residual in ["res_dx", "res_dy"]:
    for which_wheel in [-2, -1, 0, 1, 2]:
        for which_station in [1, 2, 3, 4]:
            for which_sector in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]:
                name = which_residual + " vs local y vs local x , wheel " + str(which_wheel) + ", station" + str(which_station) + ", sector " + str(which_sector)
                my2dprofile = ROOT.TProfile2D(name, name,  30, -130, 130, 30, -110, 110, -3, 3)
                my2dprofile.GetXaxis().SetTitle("local y")
                my2dprofile.GetYaxis().SetTitle("local x")
                my2dprofile.GetZaxis().SetTitle(which_residual)
                my2dprofile.Sumw2()
                #cut = which_residual + " < 100 && prop_location[1] == " + str(which_station) + " && prop_location[0] == " + str(which_wheel) + " && prop_location[2] == " + str(which_sector)
                # print("Cut is ", cut)
                myTree.Project(name, which_residual + ":prop_local_x:prop_local_y", which_residual + " < 100 && location[1] == " + str(which_station) + " && location[0] == " + str(which_wheel) + " && location[2] == " + str(which_sector) )
                #my2dprofile.SetStats(0)
                my2dprofile.Draw("colz")
                canvas.SaveAs("DTPlots/sectors/collision_MC_052522/ChamberLevel/" + str(which_residual) + "_W" + str(which_wheel) + "St" + str(which_station) + "Sec" + str(which_sector) + ".png")

# myHist.SetDirectory(0)
myFile.Close()
