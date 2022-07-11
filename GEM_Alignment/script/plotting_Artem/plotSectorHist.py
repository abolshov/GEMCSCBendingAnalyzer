import ROOT, array, math

myFile = ROOT.TFile.Open("DTPlots/DT_tbma_collisionMC_allRes.root", "READ")
#myFile = ROOT.TFile("CRUZET_march15.ROOT")
myTree = myFile.Get("DT_tbma/Inner_Prop_ChamberLevel")

canvas = ROOT.TCanvas("canvas","", 1200, 800)
canvas.SetWindowSize(500, 500)
canvas.SetRightMargin(0.15)
canvas.SetWindowSize(500, 500)

#for which_residual in ["res_dx", "res_dy"]:
which_residual = "res_dx"
# wheels = [i for i in range(-2, 3)]
# stations = [i for i in range(1, 5)]
# sectors = [i for i in range(1, 13)]
for which_wheel in [-2]:
    for which_station in [1]:
        for which_sector in [1]:
            name = which_residual + " " + str(which_wheel) + ", station" + str(which_station) + ", sector " + str(which_sector)
            title = which_residual + "_" + str(which_wheel) + "/" + str(which_station) + "/" + str(which_sector)
            myHist = ROOT.TH1D(name, title, 100, -2, 2)
            myHist.GetYaxis().SetTitle("events")
            myHist.GetXaxis().SetTitle(which_residual)
            myHist.Sumw2()
            cut = which_residual + " < 100 && location[1] == " + str(which_station) + " && location[0] == " + str(which_wheel) + " && location[2] == " + str(which_sector)
            print("Cut is: ", cut)
            # myTree.Project(name, which_residual + ":prop_local_x:prop_local_y", which_residual + " < 100 && location[1] == " + str(which_station) + " && location[0] == " + str(which_wheel) + " && location[2] == " + str(which_sector) )
            myTree.Project(name, which_residual, cut)
            myHist.Draw("h")
            #canvas.SaveAs("DT_cosmics/CRAFT_060822/sectors/ChamberLevel/Hsitograms/" + str(which_residual) + "_W" + str(which_wheel) + "St" + str(which_station) + "Sec" + str(which_sector) + ".png")
            canvas.SaveAs("DTPlots/sectors/collision_MC_052522/ChamberLevel/Hsitograms/" + title + ".png")

myHist.SetDirectory(0)
myFile.Close()
