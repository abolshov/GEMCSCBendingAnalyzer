import ROOT, sys, tdrstyle, os, array

f = ROOT.TFile("{}".format(sys.argv[1]))
event = f.Get("analyzer/ME11Seg_Prop")
RdPhi_Branch = "RdPhi_Corrected"
base_cut = "has_prop && abs(RdPhi) < 5 && n_ME11_segment == 1 && muon_pt > 3 && has_fidcut"
ROOT.gROOT.SetBatch(1)
tdrstyle.setTDRStyle()

if not os.path.exists("plots/"):
  os.makedirs("plots/")

do_1D = True
do_2D = False
do_Region = True
do_Chamber = True
do_Layer = True
do_Layer_Ratio = True

if do_1D:
  #1D Residual by Region
  if do_Region:
    plotdir = "plots/region_level/"
    if not os.path.exists(plotdir):
      os.makedirs(plotdir)
    if not os.path.exists(plotdir+"each/"):
      os.makedirs(plotdir+"each/")

    nX_plots = 2
    nY_plots = 1

    H_ref = 800*nY_plots
    W_ref = 800*nX_plots
    W = W_ref
    H = H_ref

    T = 0.12*H
    B = 0.16*H
    L = 0.16*W
    R = 0.08*W

    xbins = 100
    ybins = 100
    xlow  = -2
    xhigh = 2

    canvas = ROOT.TCanvas("c1", "c1", W, H)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetLeftMargin( L/W )
    canvas.SetRightMargin( R/W )
    canvas.SetTopMargin( T/H )
    canvas.SetBottomMargin( B/H )
    canvas.SetTickx(0)
    canvas.SetTicky(0)
    canvas.Divide(nX_plots, nY_plots)

    hlist = []
    for i, region in enumerate([-1, 1]):
      canvas.cd(i+1)
      ROOT.gPad.SetGridx()
      ROOT.gPad.SetGridy()
      cut = base_cut + " && prop_location[0] == {R} && prop_location[2] == ME11_location[3]".format(R = region)
      hlist.append(ROOT.TH1D("h{R}".format(R = region), "h{R}".format(R = region), xbins, xlow, xhigh))
      xAxis = hlist[i].GetXaxis()
      xAxis.SetTitleOffset(2)
      xAxis.SetTitleSize(0.04)
      xAxis.SetTitle("#DeltaR#phi [cm]")

      yAxis = hlist[i].GetYaxis()
      yAxis.SetTitleOffset(0)
      yAxis.SetTitleSize(0.04)
      yAxis.SetTitle("Entries")

      event.Project("h{R}".format(R = region), RdPhi_Branch, cut)

      hlist[i].SetLineWidth(3)
      hlist[i].Draw()

      latex = ROOT.TLatex()
      latex.SetNDC()
      latex.SetTextAngle(0)
      latex.SetTextColor(ROOT.kBlack)

      latex.SetTextFont(42)
      latex.SetTextSize(0.8*ROOT.gPad.GetTopMargin())

      latex.SetTextAlign(32)
      latex.DrawLatex(1-1.1*canvas.GetRightMargin(), 1 - 0.6*canvas.GetTopMargin(), "Entries: {entries}".format(entries = int(hlist[i].GetEntries())))
      latex.DrawLatex(1-1.1*canvas.GetRightMargin(), 1 - 1.0*canvas.GetTopMargin(), "Mean: {mean}".format(mean = round(hlist[i].GetMean(),3)))
      latex.DrawLatex(1-1.1*canvas.GetRightMargin(), 1 - 1.4*canvas.GetTopMargin(), "Std Dev: {stddev}".format(stddev = round(hlist[i].GetStdDev(),3)))
      latex.SetTextAlign(12)
      latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1 - 0.6*canvas.GetTopMargin(), "Region: {R}".format(R = region))

      latex.SetTextSize(0.8*ROOT.gPad.GetTopMargin())
      latex.SetTextFont(61)
      latex.DrawLatex(ROOT.gPad.GetLeftMargin(), 1 - 0.4*ROOT.gPad.GetTopMargin(), "CMS")
      latex.SetTextFont(52)
      latex.SetTextSize(0.6*ROOT.gPad.GetTopMargin())
      #latex.DrawLatex(1.9*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "Preliminary")
      latex.DrawLatex(1.9*ROOT.gPad.GetLeftMargin(), 1 - 0.45*ROOT.gPad.GetTopMargin(), "Work in Progress")

      latex.SetTextFont(42)
      latex.SetTextSize(0.5*ROOT.gPad.GetTopMargin())

      frame = canvas.GetFrame()
      frame.Draw()

      ROOT.gPad.SaveAs(plotdir+"each/res1D_R{R}.png".format(R = region))

    canvas.SaveAs(plotdir+"res1D.png")


  #1D Residual by Chamber
  if do_Chamber:
    plotdir = "plots/chamber_level/"
    if not os.path.exists(plotdir):
      os.makedirs(plotdir)
    if not os.path.exists(plotdir+"each/"):
      os.makedirs(plotdir+"each/")

    nX_plots = 9
    nY_plots = 4

    H_ref = 800
    W_ref = 800
    W = W_ref*nX_plots
    H = H_ref*nY_plots

    T = 0.12*H
    B = 0.16*H
    L = 0.16*W
    R = 0.08*W

    xbins = 100
    ybins = 100
    xlow  = -2
    xhigh = 2

    canvas = ROOT.TCanvas("c1", "c1", W, H)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetLeftMargin( L/W )
    canvas.SetRightMargin( R/W )
    canvas.SetTopMargin( T/H )
    canvas.SetBottomMargin( B/H )
    canvas.SetTickx(0)
    canvas.SetTicky(0)

    canvas.Divide(nX_plots, nY_plots)

    for region in [-1, 1]:
      hlist = []
      for i, chamber in enumerate(range(1, 37)):
        canvas.cd(i+1)
        ROOT.gPad.SetGridx()
        ROOT.gPad.SetGridy()
        cut = base_cut + " && prop_location[0] == {R} && prop_location[2] == {Ch} && prop_location[2] == ME11_location[3]".format(R = region, Ch = chamber)
        hlist.append(ROOT.TH1D("h{R}_{Ch}".format(R = region, Ch = chamber), "h{R}_{Ch}".format(R = region, Ch = chamber), xbins, xlow, xhigh))
        xAxis = hlist[i].GetXaxis()
        xAxis.SetTitleOffset(2)
        xAxis.SetTitleSize(0.04)
        xAxis.SetTitle("#DeltaR#phi [cm]")

        yAxis = hlist[i].GetYaxis()
        yAxis.SetTitleOffset(0)
        yAxis.SetTitleSize(0.04)
        yAxis.SetTitle("Entries")

        event.Project("h{R}_{Ch}".format(R = region, Ch = chamber), RdPhi_Branch, cut)

        hlist[i].SetLineWidth(3)
        hlist[i].Draw()

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextAngle(0)
        latex.SetTextColor(ROOT.kBlack)

        latex.SetTextFont(42)
        latex.SetTextSize(0.8*ROOT.gPad.GetTopMargin())

        latex.SetTextAlign(32)
        latex.DrawLatex(1-1.1*canvas.GetRightMargin(), 1 - 0.6*canvas.GetTopMargin(), "Entries: {entries}".format(entries = int(hlist[i].GetEntries())))
        latex.DrawLatex(1-1.1*canvas.GetRightMargin(), 1 - 1.0*canvas.GetTopMargin(), "Mean: {mean}".format(mean = round(hlist[i].GetMean(),3)))
        latex.DrawLatex(1-1.1*canvas.GetRightMargin(), 1 - 1.4*canvas.GetTopMargin(), "Std Dev: {stddev}".format(stddev = round(hlist[i].GetStdDev(),3)))
        latex.SetTextAlign(12)
        latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1 - 0.6*canvas.GetTopMargin(), "Region: {R}".format(R = region))
        latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1- 1.0*canvas.GetTopMargin(), "Chamber: {Ch}".format(Ch = chamber))

        latex.SetTextSize(0.8*ROOT.gPad.GetTopMargin())
        latex.SetTextFont(61)
        latex.DrawLatex(ROOT.gPad.GetLeftMargin(), 1 - 0.4*ROOT.gPad.GetTopMargin(), "CMS")
        latex.SetTextFont(52)
        latex.SetTextSize(0.6*ROOT.gPad.GetTopMargin())
        #latex.DrawLatex(1.9*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "Preliminary")
        latex.DrawLatex(1.9*ROOT.gPad.GetLeftMargin(), 1 - 0.45*ROOT.gPad.GetTopMargin(), "Work in Progress")

        latex.SetTextFont(42)
        latex.SetTextSize(0.5*ROOT.gPad.GetTopMargin())

        frame = canvas.GetFrame()
        frame.Draw()

        ROOT.gPad.SaveAs(plotdir+"each/res1D_R{R}_Ch{Ch}.png".format(R = region, Ch = chamber))

      canvas.SaveAs(plotdir+"res1D_Region{R}.png".format(R = region))


  #1D Residual by Chamber + Layer
  if do_Layer:
    plotdir = "plots/layer_level/"
    if not os.path.exists(plotdir):
      os.makedirs(plotdir)
    if not os.path.exists(plotdir+"each/"):
      os.makedirs(plotdir+"each/")

    nX_plots = 9
    nY_plots = 4

    H_ref = 800
    W_ref = 800
    W = W_ref*nX_plots
    H = H_ref*nY_plots

    T = 0.12*H
    B = 0.16*H
    L = 0.16*W
    R = 0.08*W

    xbins = 100
    ybins = 100
    xlow  = -2
    xhigh = 2

    canvas = ROOT.TCanvas("c1", "c1", W, H)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetLeftMargin( L/W )
    canvas.SetRightMargin( R/W )
    canvas.SetTopMargin( T/H )
    canvas.SetBottomMargin( B/H )
    canvas.SetTickx(0)
    canvas.SetTicky(0)

    canvas.Divide(nX_plots, nY_plots)

    for region in [-1, 1]:
      for layer in [1, 2]:
        hlist = []
        for i, chamber in enumerate(range(1, 37)):
          canvas.cd(i+1)
          ROOT.gPad.SetGridx()
          ROOT.gPad.SetGridy()
          cut = base_cut + " && prop_location[0] == {R} && prop_location[2] == {Ch} && prop_location[3] == {lay} && prop_location[2] == ME11_location[3]".format(R = region, Ch = chamber, lay = layer)
          hlist.append(ROOT.TH1D("h{R}_{Ch}_{lay}".format(R = region, Ch = chamber, lay = layer), "h{R}_{Ch}_{lay}".format(R = region, Ch = chamber, lay = layer), xbins, xlow, xhigh))
          xAxis = hlist[i].GetXaxis()
          xAxis.SetTitleOffset(2)
          xAxis.SetTitleSize(0.04)
          xAxis.SetTitle("#DeltaR#phi [cm]")

          yAxis = hlist[i].GetYaxis()
          yAxis.SetTitleOffset(0)
          yAxis.SetTitleSize(0.04)
          yAxis.SetTitle("Entries")

          event.Project("h{R}_{Ch}_{lay}".format(R = region, Ch = chamber, lay = layer), RdPhi_Branch, cut)

          hlist[i].SetLineWidth(3)
          hlist[i].Draw()

          latex = ROOT.TLatex()
          latex.SetNDC()
          latex.SetTextAngle(0)
          latex.SetTextColor(ROOT.kBlack)

          latex.SetTextFont(42)
          latex.SetTextSize(0.8*ROOT.gPad.GetTopMargin())

          latex.SetTextAlign(32)
          latex.DrawLatex(1-1.1*canvas.GetRightMargin(), 1 - 0.6*canvas.GetTopMargin(), "Entries: {entries}".format(entries = int(hlist[i].GetEntries())))
          latex.DrawLatex(1-1.1*canvas.GetRightMargin(), 1 - 1.0*canvas.GetTopMargin(), "Mean: {mean}".format(mean = round(hlist[i].GetMean(),3)))
          latex.DrawLatex(1-1.1*canvas.GetRightMargin(), 1 - 1.4*canvas.GetTopMargin(), "Std Dev: {stddev}".format(stddev = round(hlist[i].GetStdDev(),3)))
          latex.SetTextAlign(12)
          latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1 - 0.6*canvas.GetTopMargin(), "Region: {R}".format(R = region))
          latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1- 1.0*canvas.GetTopMargin(), "Chamber: {Ch}".format(Ch = chamber))
          latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1 - 1.4*canvas.GetTopMargin(), "Layer: {lay}".format(lay = layer))

          latex.SetTextSize(0.8*ROOT.gPad.GetTopMargin())
          latex.SetTextFont(61)
          latex.DrawLatex(ROOT.gPad.GetLeftMargin(), 1 - 0.4*ROOT.gPad.GetTopMargin(), "CMS")
          latex.SetTextFont(52)
          latex.SetTextSize(0.6*ROOT.gPad.GetTopMargin())
          #latex.DrawLatex(1.9*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "Preliminary")
          latex.DrawLatex(1.9*ROOT.gPad.GetLeftMargin(), 1 - 0.45*ROOT.gPad.GetTopMargin(), "Work in Progress")

          latex.SetTextFont(42)
          latex.SetTextSize(0.5*ROOT.gPad.GetTopMargin())

          frame = canvas.GetFrame()
          frame.Draw()

          canvas.SaveAs(plotdir+"each/res1D_Region{R}_Ch{Ch}_Layer{lay}.png".format(R = region, Ch = chamber, lay = layer))

        canvas.SaveAs(plotdir+"res1D_Region{R}_Layer{lay}.png".format(R = region, lay = layer))


  #1D Residual by Chamber + Layer RATIO
  if do_Layer_Ratio:
    plotdir = "plots/layer_ratio_level/"
    if not os.path.exists(plotdir):
      os.makedirs(plotdir)
    if not os.path.exists(plotdir+"each/"):
      os.makedirs(plotdir+"each/")

    nX_plots = 9
    nY_plots = 4

    H_ref = 800
    W_ref = 800
    W = W_ref*nX_plots
    H = H_ref*nY_plots

    T = 0.12*H
    B = 0.16*H
    L = 0.16*W
    R = 0.08*W

    xbins = 20
    ybins = 100
    xlow  = -0.5
    xhigh = 0.5

    canvas = ROOT.TCanvas("c1", "c1", W, H)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetLeftMargin( L/W )
    canvas.SetRightMargin( R/W )
    canvas.SetTopMargin( T/H )
    canvas.SetBottomMargin( B/H )
    canvas.SetTickx(0)
    canvas.SetTicky(0)

    canvas.Divide(nX_plots, nY_plots)




    for region in [-1, 1]:
      hlist_l1 = []
      hlist_l2 = []
      flist = []

      for i, chamber in enumerate(range(1, 37)):
        canvas.cd(i+1)
        ROOT.gPad.SetGridx()
        ROOT.gPad.SetGridy()
        cut = base_cut + " && prop_location[0] == {R} && prop_location[2] == {Ch} && prop_location[2] == ME11_location[3]".format(R = region, Ch = chamber)
        hlist_l1.append(ROOT.TH1D("h{R}_{Ch}_1_rat".format(R = region, Ch = chamber), "h{R}_{Ch}_1_rat".format(R = region, Ch = chamber), xbins, xlow, xhigh))
        hlist_l2.append(ROOT.TH1D("h{R}_{Ch}_2_rat".format(R = region, Ch = chamber), "h{R}_{Ch}_2_rat".format(R = region, Ch = chamber), xbins, xlow, xhigh))
        xAxis = hlist_l1[i].GetXaxis()
        xAxis.SetTitleOffset(2)
        xAxis.SetTitleSize(0.04)
        xAxis.SetTitle("#DeltaR#phi [cm]")

        yAxis = hlist_l1[i].GetYaxis()
        yAxis.SetTitleOffset(0)
        yAxis.SetTitleSize(0.04)
        yAxis.SetTitle("Entries")

        event.Project("h{R}_{Ch}_1_rat".format(R = region, Ch = chamber), RdPhi_Branch, cut+" && prop_location[3] == 1")
        event.Project("h{R}_{Ch}_2_rat".format(R = region, Ch = chamber), RdPhi_Branch, cut+" && prop_location[3] == 2")

        hlist_l1[i].SetLineWidth(3)
        hlist_l1[i].Divide(hlist_l2[i])
        hlist_l1[i].Draw()
        flist.append(ROOT.TF1("f{R}_{Ch}_rat".format(R = region, Ch = chamber), "[0]*x + [1]", -0.5, 0.5))
        flist[i].SetParameters(.1, 1.0)
        hlist_l1[i].Fit("f{R}_{Ch}_rat".format(R = region, Ch = chamber))
        flist[i].SetLineColor(ROOT.kRed)
        flist[i].Draw("same")


        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextAngle(0)
        #latex.SetTextColor(ROOT.kBlack)

        latex.SetTextFont(42)
        latex.SetTextSize(0.8*ROOT.gPad.GetTopMargin())

        latex.SetTextAlign(12)
        latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1 - 0.6*canvas.GetTopMargin(), "Region: {R}".format(R = region))
        latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1- 1.0*canvas.GetTopMargin(), "Chamber: {Ch}".format(Ch = chamber))

        latex.SetTextSize(0.8*ROOT.gPad.GetTopMargin())
        latex.SetTextFont(61)
        latex.DrawLatex(ROOT.gPad.GetLeftMargin(), 1 - 0.4*ROOT.gPad.GetTopMargin(), "CMS")
        latex.SetTextFont(52)
        latex.SetTextSize(0.6*ROOT.gPad.GetTopMargin())
        #latex.DrawLatex(1.9*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "Preliminary")
        latex.DrawLatex(1.9*ROOT.gPad.GetLeftMargin(), 1 - 0.45*ROOT.gPad.GetTopMargin(), "Work in Progress")

        latex.SetTextFont(42)
        latex.SetTextSize(0.8*ROOT.gPad.GetTopMargin())

        latex.DrawLatex(1.1*ROOT.gPad.GetLeftMargin(), 1 - 4*ROOT.gPad.GetTopMargin(), "{zer} x + {one}".format(zer = round(flist[i].GetParameter(0), 4), one = round(flist[i].GetParameter(1), 4)))


        #frame = canvas.GetFrame()
        #frame.Draw()

        #ROOT.gPad.SaveAs(plotdir+"each/ratio_Region{R}_Ch{Ch}.png".format(R = region, Ch = chamber))

      canvas.SaveAs(plotdir+"ratio_Region{R}.png".format(R = region))


if do_2D:
  #2D Residual by Region
  if do_Region:
    plotdir = "plots/region_level/"
    if not os.path.exists(plotdir):
      os.makedirs(plotdir)
    if not os.path.exists(plotdir+"each/"):
      os.makedirs(plotdir+"each/")

    nX_plots = 2
    nY_plots = 1

    H_ref = 800
    W_ref = 1000
    W = W_ref*nX_plots
    H = H_ref*nY_plots

    T = 0.08*H
    B = 0.16*H
    L = 0.16*W
    R = 0.20*W

    xbins = 100
    ybins = 100
    xlow  = -300
    xhigh = 300
    ylow  = -300
    yhigh = 300
    zlow  = -2
    zhigh = 2

    canvas = ROOT.TCanvas("c1", "c1", W, H)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetTickx(0)
    canvas.SetTicky(0)
    canvas.Divide(nX_plots, nY_plots)


    hlist = []
    for i, region in enumerate([-1, 1]):
      canvas.cd(i+1)
      ROOT.gPad.SetGridx()
      ROOT.gPad.SetGridy()
      ROOT.gPad.SetLeftMargin( L/W )
      ROOT.gPad.SetRightMargin( R/W )
      ROOT.gPad.SetTopMargin( T/H )
      ROOT.gPad.SetBottomMargin( B/H )
      cut = base_cut +" && prop_location[0] == {R}".format(R = region)
      hlist.append(ROOT.TProfile2D("h{R}".format(R = region), "h{R}".format(R = region), xbins, xlow, xhigh, ybins, ylow, yhigh))
      ablue = array.array("d", [1,1,0])
      ared = array.array("d", [0,1,1])
      agreen = array.array("d", [0,1,0])
      astop = array.array("d", [0,.5,1])
      myPalette = []
      fi = ROOT.TColor.CreateGradientColorTable(3, astop, ared, agreen, ablue, 100)
      for x in range(100):
        myPalette.append(fi+x)
      ROOT.gStyle.SetPalette(100, array.array("i", myPalette))

      xAxis = hlist[i].GetXaxis()
      xAxis.SetTitleOffset(2)
      xAxis.SetTitleSize(0.04)
      xAxis.SetTitle("Global x [cm]")

      yAxis = hlist[i].GetYaxis()
      yAxis.SetTitleOffset(0)
      yAxis.SetTitleSize(0.04)
      yAxis.SetTitle("Global y [cm]")

      zAxis = hlist[i].GetZaxis()
      zAxis.SetTitle("#DeltaR#phi [cm]")
      zAxis.SetTitleOffset(2)
      zAxis.SetTitleSize(0.04)
      zAxis.SetRangeUser(zlow, zhigh)

      event.Project("h{R}".format(R = region), RdPhi_Branch+":prop_GP[1]:prop_GP[0]", cut)

      hlist[i].Draw("colz")

      latex = ROOT.TLatex()
      latex.SetNDC()
      latex.SetTextAngle(0)
      latex.SetTextColor(ROOT.kBlack)

      latex.SetTextFont(42)
      latex.SetTextSize(0.6*ROOT.gPad.GetTopMargin())

      latex.SetTextAlign(32)
      latex.DrawLatex(1-1.1*ROOT.gPad.GetRightMargin(), 1 - 1.5*ROOT.gPad.GetTopMargin(), "Entries: {entries}".format(entries = int(hlist[i].GetEntries())))
      latex.SetTextAlign(12)
      latex.DrawLatex(0+1.1*ROOT.gPad.GetLeftMargin(), 1 - 1.5*ROOT.gPad.GetTopMargin(), "Region: {R}".format(R = region))

      latex.SetTextSize(0.8*ROOT.gPad.GetTopMargin())
      latex.SetTextFont(61)
      latex.DrawLatex(ROOT.gPad.GetLeftMargin(), 1 - 0.4*ROOT.gPad.GetTopMargin(), "CMS")
      latex.SetTextFont(52)
      latex.SetTextSize(0.6*ROOT.gPad.GetTopMargin())
      #latex.DrawLatex(1.9*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "Preliminary")
      latex.DrawLatex(1.9*ROOT.gPad.GetLeftMargin(), 1 - 0.45*ROOT.gPad.GetTopMargin(), "Work in Progress")

      latex.SetTextFont(42)
      latex.SetTextSize(0.5*ROOT.gPad.GetTopMargin())

      frame = canvas.GetFrame()
      frame.Draw()

      canvas.SaveAs(plotdir+"each/res2D_R{R}.png".format(R = region))


    canvas.SaveAs(plotdir+"res2D.png")


  #2D Residual by Chamber
  if do_Chamber:
    plotdir = "plots/chamber_level/"
    if not os.path.exists(plotdir):
      os.makedirs(plotdir)
    if not os.path.exists(plotdir+"each/"):
      os.makedirs(plotdir+"each/")

    nX_plots = 9
    nY_plots = 4

    H_ref = 800
    W_ref = 1000
    W = W_ref*nX_plots
    H = H_ref*nY_plots

    T = 0.08*H
    B = 0.16*H
    L = 0.20*W
    R = 0.20*W

    xbins = 100
    ybins = 100
    xlow  = -30
    xhigh = 30
    ylow  = -60
    yhigh = 60
    zlow  = -2
    zhigh = 2

    canvas = ROOT.TCanvas("c1", "c1", W, H)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetTickx(0)
    canvas.SetTicky(0)

    canvas.Divide(nX_plots, nY_plots)

    for region in [-1, 1]:
      hlist = []
      for i, chamber in enumerate(range(1, 37)):
        canvas.cd(i+1)
        ROOT.gPad.SetGridx()
        ROOT.gPad.SetGridy()
        ROOT.gPad.SetLeftMargin( L/W )
        ROOT.gPad.SetRightMargin( R/W )
        ROOT.gPad.SetTopMargin( T/H )
        ROOT.gPad.SetBottomMargin( B/H )
        cut = base_cut + " && prop_location[0] == {R} && prop_location[2] == {Ch}".format(R = region, Ch = chamber)
        hlist.append(ROOT.TProfile2D("h{R}_{Ch}".format(R = region, Ch = chamber), "h{R}_{Ch}".format(R = region, Ch = chamber), xbins, xlow, xhigh, ybins, ylow, yhigh))
        ablue = array.array("d", [1,1,0])
        ared = array.array("d", [0,1,1])
        agreen = array.array("d", [0,1,0])
        astop = array.array("d", [0,.5,1])
        myPalette = []
        fi = ROOT.TColor.CreateGradientColorTable(3, astop, ared, agreen, ablue, 100)
        for x in range(100):
          myPalette.append(fi+x)
        ROOT.gStyle.SetPalette(100, array.array("i", myPalette))

        xAxis = hlist[i].GetXaxis()
        xAxis.SetTitleOffset(2)
        xAxis.SetTitleSize(0.04)
        xAxis.SetTitle("#DeltaR#phi [cm]")

        yAxis = hlist[i].GetYaxis()
        yAxis.SetTitleOffset(0)
        yAxis.SetTitleSize(0.04)
        yAxis.SetTitle("Entries")

        zAxis = hlist[i].GetZaxis()
        zAxis.SetTitle("RdPhi [cm]")
        zAxis.SetTitleOffset(2)
        zAxis.SetTitleSize(0.04)
        zAxis.SetRangeUser(zlow, zhigh)

        event.Project("h{R}_{Ch}".format(R = region, Ch = chamber), RdPhi_Branch+":prop_LP[1]:prop_LP[0]", cut)

        hlist[i].Draw("colz")

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextAngle(0)
        latex.SetTextColor(ROOT.kBlack)

        latex.SetTextFont(42)
        latex.SetTextSize(0.6*ROOT.gPad.GetTopMargin())

        latex.SetTextAlign(32)
        latex.DrawLatex(1-1.1*ROOT.gPad.GetRightMargin(), 1 - 1.5*ROOT.gPad.GetTopMargin(), "Entries: {entries}".format(entries = int(hlist[i].GetEntries())))
        latex.SetTextAlign(12)
        latex.DrawLatex(0+1.1*ROOT.gPad.GetLeftMargin(), 1 - 1.5*ROOT.gPad.GetTopMargin(), "Region: {R}".format(R = region))
        latex.DrawLatex(0+1.1*ROOT.gPad.GetLeftMargin(), 1- 1.0*ROOT.gPad.GetTopMargin(), "Chamber: {Ch}".format(Ch = chamber))


        latex.SetTextSize(0.8*ROOT.gPad.GetTopMargin())
        latex.SetTextFont(61)
        latex.DrawLatex(ROOT.gPad.GetLeftMargin(), 1 - 0.4*ROOT.gPad.GetTopMargin(), "CMS")
        latex.SetTextFont(52)
        latex.SetTextSize(0.6*ROOT.gPad.GetTopMargin())
        #latex.DrawLatex(1.9*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "Preliminary")
        latex.DrawLatex(1.9*ROOT.gPad.GetLeftMargin(), 1 - 0.45*ROOT.gPad.GetTopMargin(), "Work in Progress")

        latex.SetTextFont(42)
        latex.SetTextSize(0.5*ROOT.gPad.GetTopMargin())

        frame = canvas.GetFrame()
        frame.Draw()

        canvas.SaveAs(plotdir+"each/res2D_Region{R}_Ch{Ch}.png".format(R = region, Ch = chamber))

      canvas.SaveAs(plotdir+"res2D_Region{R}.png".format(R = region))


  #2D Residual by Chamber + Layer
  if do_Layer:
    plotdir = "plots/layer_level/"
    if not os.path.exists(plotdir):
      os.makedirs(plotdir)
    if not os.path.exists(plotdir+"each/"):
      os.makedirs(plotdir+"each/")

    nX_plots = 9
    nY_plots = 4

    H_ref = 800
    W_ref = 1000
    W = W_ref*nX_plots
    H = H_ref*nY_plots

    T = 0.08*H
    B = 0.16*H
    L = 0.12*W
    R = 0.20*W

    xbins = 100
    ybins = 100
    xlow  = -30
    xhigh = 30
    ylow  = -70
    yhigh = 70
    zlow  = -2
    zhigh = 2

    canvas = ROOT.TCanvas("c1", "c1", W, H)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetTickx(0)
    canvas.SetTicky(0)

    canvas.Divide(nX_plots, nY_plots)

    for region in [-1, 1]:
      for layer in [1, 2]:
        hlist = []
        for i, chamber in enumerate(range(1, 37)):
          canvas.cd(i+1)
          ROOT.gPad.SetGridx()
          ROOT.gPad.SetGridy()
          ROOT.gPad.SetLeftMargin( L/W )
          ROOT.gPad.SetRightMargin( R/W )
          ROOT.gPad.SetTopMargin( T/H )
          ROOT.gPad.SetBottomMargin( B/H )
          cut = base_cut + " && prop_location[0] == {R} && prop_location[2] == {Ch} && prop_location[3] == {lay}".format(R = region, Ch = chamber, lay = layer)
          hlist.append(ROOT.TProfile2D("h{R}_{Ch}_{lay}".format(R = region, Ch = chamber, lay = layer), "h{R}_{Ch}_{lay}".format(R = region, Ch = chamber, lay = layer), xbins, xlow, xhigh, ybins, ylow, yhigh))
          ablue = array.array("d", [1,1,0])
          ared = array.array("d", [0,1,1])
          agreen = array.array("d", [0,1,0])
          astop = array.array("d", [0,.5,1])
          myPalette = []
          fi = ROOT.TColor.CreateGradientColorTable(3, astop, ared, agreen, ablue, 100)
          for x in range(100):
            myPalette.append(fi+x)
          ROOT.gStyle.SetPalette(100, array.array("i", myPalette))

          xAxis = hlist[i].GetXaxis()
          xAxis.SetTitleOffset(2)
          xAxis.SetTitleSize(0.04)
          xAxis.SetTitle("#DeltaR#phi [cm]")

          yAxis = hlist[i].GetYaxis()
          yAxis.SetTitleOffset(0)
          yAxis.SetTitleSize(0.04)
          yAxis.SetTitle("Entries")

          zAxis = hlist[i].GetZaxis()
          zAxis.SetTitle("RdPhi [cm]")
          zAxis.SetTitleOffset(2)
          zAxis.SetTitleSize(0.04)
          zAxis.SetRangeUser(zlow, zhigh)

          event.Project("h{R}_{Ch}_{lay}".format(R = region, Ch = chamber, lay = layer), RdPhi_Branch+":prop_LP[1]:prop_LP[0]", cut)

          hlist[i].Draw("colz")

          latex = ROOT.TLatex()
          latex.SetNDC()
          latex.SetTextAngle(0)
          latex.SetTextColor(ROOT.kBlack)

          latex.SetTextFont(42)
          latex.SetTextSize(0.6*ROOT.gPad.GetTopMargin())

          latex.SetTextAlign(32)
          latex.DrawLatex(1-1.1*ROOT.gPad.GetRightMargin(), 1 - 1.5*ROOT.gPad.GetTopMargin(), "Entries: {entries}".format(entries = int(hlist[i].GetEntries())))
          latex.SetTextAlign(12)
          latex.DrawLatex(0+1.1*ROOT.gPad.GetLeftMargin(), 1 - 1.4*ROOT.gPad.GetTopMargin(), "Region: {R}".format(R = region))
          latex.DrawLatex(0+1.1*ROOT.gPad.GetLeftMargin(), 1- 1.8*ROOT.gPad.GetTopMargin(), "Chamber: {Ch}".format(Ch = chamber))
          latex.DrawLatex(0+1.1*ROOT.gPad.GetLeftMargin(), 1 - 2.2*ROOT.gPad.GetTopMargin(), "Layer: {lay}".format(lay = layer))


          latex.SetTextSize(0.8*ROOT.gPad.GetTopMargin())
          latex.SetTextFont(61)
          latex.DrawLatex(ROOT.gPad.GetLeftMargin(), 1 - 0.4*ROOT.gPad.GetTopMargin(), "CMS")
          latex.SetTextFont(52)
          latex.SetTextSize(0.6*ROOT.gPad.GetTopMargin())
          #latex.DrawLatex(1.9*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "Preliminary")
          latex.DrawLatex(1.9*ROOT.gPad.GetLeftMargin(), 1 - 0.45*ROOT.gPad.GetTopMargin(), "Work in Progress")

          latex.SetTextFont(42)
          latex.SetTextSize(0.5*ROOT.gPad.GetTopMargin())

          frame = canvas.GetFrame()
          frame.Draw()

          canvas.SaveAs(plotdir+"each/res2D_Region{R}_Ch{Ch}_Layer{lay}.png".format(R = region, Ch = chamber, lay = layer))

        canvas.SaveAs(plotdir+"res2D_Region{R}_Layer{lay}.png".format(R = region, lay = layer))
