import ROOT, os
from make_plot import *

f = ROOT.TFile("out_file.root")

event = f.Get("analyzer/CSC_Prop")
event = f.Get("analyzer/Inner_Prop")
event = f.Get("analyzer/ME11Seg_Prop")

if os.path.exists("plots") == False:
  os.mkdir("plots/")


######################################__Branch Legend__#############################################
#Muons:                                                                                            #
#muon_charge	muon_pt		muon_eta 	muon_momentum                                      #
#evtNum		lumiBlock	muonIdx		runNum						   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  # # # #
#has_prop                                                                has successful propagation#
#prop_location[5]  prop_location[0] = region, [1] = station, [2] = chamber, [3] = layer, [4] = roll#
#prop_GP[3]                               Global propagation destination: [0] = x, [1] = y, [2] = z#
#prop_LP[3]                                Local propagation destination: [0] = x, [1] = y, [2] = z#
#prop_startingPoint_GP[3]                       Global starting position: [0] = x, [1] = y, [2] = z#
#prop_yroll		                                    fixed local y for full chamber plotting#
#prop_localphi_rad  	                                 Propagated position's local phi in radians#
#prop_localphi_deg	                                  Propagated position's localphi in degrees#
#has_fidcut	              current fidcut is |phi| < 4degrees, 5cm off top and bottom of chamber#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  # # # #
#Track:												   #
#track_chi2									 chi2 of track used#
#track_ndof									 ndof of track used#
#n_ME11_segment							   Number of ME11 segments in event#
#which_track							  Incoming (0) or Outgoing (1) muon#
#hasME11								      Includes ME11 segment#
#hasME11RecHit									 Hes rechit on ME11#
#hasME11A								     Includes ME11A segment#
#hasME11ARecHit								      Includes ME11A rechit#
#nCSCSeg								     Number of CSC segments#
#nDTSeg									      Number of DT segments#
#nME11RecHits								     Number of ME11 RecHits#
#ME11_BunchX									ME11 Bunch Crossing#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  # # # #
#Rechit:                                                                                           #
#has_rechit	                                                                                   #
#rechit_location[5]             [0] = region, [1] = station, [2] = chamber, [3] = layer, [4] = roll#
#rechit_GP[3]                                  Global position of rechit: [0] = x, [1] = y, [2] = z#
#rechit_LP[3]                                   Local position of rechit: [0] = x, [1] = y, [2] = z#
#rechit_yroll									   Local-y per roll#
#rechit_localphi_rad							       Local phi in radians#
#rechit_localphi_deg							       Local phi in degrees#
#rechit_first_strip                                                   		first cluster strip#
#rechit_CLS                                                                  	       cluster size#
#rechit_BunchX                                                             	     bunch crossing#
#RdPhi                                                          	      best matched residual#
#RdPhi_Corrected                      Residuals flipped for short/long (flips R+1 odd and R-1 even)#
#rechit_detId								Detector ID for matched hit#
#nRecHitsTot								 Number of rechits in event#
#nRecHits2								    Number with RdPhi < 2cm#
#nRecHits5								    Number with RdPhi < 5cm#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  # # # #
#Sim:												   #
#sim_GP[3]                                        Global simhit position: [0] = x, [1] = y, [2] = z#
#sim_LP[3]                                         local simhit position: [0] = x, [1] = y, [2] = z#
#simDy                                 prop_global_y - sim_global_y !!!should be local, will change#
#sim_yroll									   Local-y per roll#
#nSim                                                Counts simhits on same reg/sta/ch/lay/dRoll<=1#
####################################################################################################



#Example of 1D hist#
#plot1Dhist("branch", "Title of plot", event, [xbins, xlow, xhigh], "xaxis title", "yaxis title", "cuts1 && cuts2", "subdir", Save?, Logscale Y?)#
plot1Dhist("muon_pt", "Muon pt", event, [200, 0, 200], "muon pt [GeV]", "Entries", "", "pt", True, False)

#Example of 1D efficiency#
#plot1Defficiency("branch", "Title of plot", event, [xbins, xlow, xhigh], "xaxis title", "yaxis title", "total cut", "passed cut", "subdir", Save?)#
plot1Defficiency("prop_location[2]", "Chamber Efficiency", event, [37, 0, 37], "Chamber", "Eff", "has_prop == 1 && prop_location[4] < 10 && hasME11", "has_rechit == 1 && abs(RdPhi) < 5 && has_prop == 1 && prop_location[4] < 10 && hasME11", "makeEffi", True)

#Example of 2D hist / chamber&region#
#plot2Dhist("ybranch:xbranch", "Title of plot", event, [xbins, xlow, xhigh, ybins, ylow, yhigh], "x title", "y title", "cuts1 && cuts2", "subdir", Save?)
for region in [-1, 1]:
  for chamber in range(1, 37):
    plot2Dhist("prop_GP[1]:prop_GP[0]", "GEM Propagation Global region {i} chamber {j}".format(i = region, j = chamber), event, [100, -300, 300, 100, -300, 300], "prop x [cm]", "prop y [cm]", "prop_GP[1] < 1000 && prop_GP[0] < 1000 && prop_location[0] == {i} && prop_location[2] == {j}".format(i = region, j = chamber), "prop_2D", True)

#Example of 2D profile#
#plot2Dprofile("zbranch:ybranch:xbranch", "Title of plot", event, [xbins, xlow, xhigh, ybins, ylow, yhigh], "x title", "y title", "cuts1 && cuts2", "subdir", Save?)
plot2Dprofile("RdPhi:prop_GP[1]:prop_GP[0]", "Profile Test", event, [300, -300, 300, 300, -300, 300, -4, 4], "Global x", "Global y", "has_prop && abs(RdPhi) < 100", "example", True)
