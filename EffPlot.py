import ROOT
from ROOT import gROOT
from ROOT import gStyle
from ROOT import gPad 
from ROOT import TF1, TH1F
from ROOT import kRed, kBlue, kViolet, kBlack, kGreen, kYellow 
from ROOT import TPaveText 
from ROOT import TPaveStats
from ROOT import TLatex
from ROOT import TLegend 
from ROOT import TGraph 
from ROOT import TMultiGraph 
ROOT.gROOT.SetBatch(True)  # go into batch mode, stop trying to actively popup graphics on the screen
from subprocess import call
import string, re
from array import array
import sys

inputname = "EffStudy.root"
inputfile = ROOT.TFile.Open( inputname )

h_eff_phi = inputfile.Get("eff_phi")
h_eff_eta = inputfile.Get("eff_eta")

h_eff_phi_good = inputfile.Get("eff_phi_good")
h_eff_eta_good = inputfile.Get("eff_eta_good")


layers = []
mpv = []

for i in range(-16,17):
	for j in range(1,5):
		myhist = inputfile.Get("qpeak_sector_%d_layer_%d" % (i, j))
		if myhist == None: 
			layer= i+0.25*(j-1.0) 
			layers.append(float(layer))
			mpv.append(0.0) 
			continue
		gStyle.SetOptFit(11)
		gStyle.SetOptStat(0)

		layer= i+0.25*(j-1.0) 
		layers.append(float(layer))

		myc= ROOT.TCanvas("myc", "myc", 800, 600)
		myc.SetMargin(0.15,0.075,0.15,0.075)
		myhist.SetTitle("Sector_%d, Layer_%d" % (i,j))
		myhist.GetYaxis().SetTitle("Events")
		myhist.GetYaxis().SetTitleOffset(1.4)
		myhist.GetYaxis().SetTitleSize(0.04)
		myhist.GetXaxis().SetTitle("qpeak")
		myhist.Draw()

		fit_x = TF1("fit_x", "landau", 0., 4000000.)
#		fit_x = TF1("fit_x", "landau", 0., 6000000.)
		fit_x.SetLineColor(kRed)
		myhist.Fit("fit_x")

		mpv.append(fit_x.GetParameter(1)) 

		myc.Update()

		myc.SaveAs("effplots/qpeak/qpeak_sector_%d_layer_%d.eps" % (i,j))
		myc.SaveAs("effplots/qpeak/qpeak_sector_%d_layer_%d.png" % (i,j))

		myc.Close()

for k in range (0, len(layers)):
	print layers[k], mpv[k]

#print layers 
#print mpv
#print len(layers), len(mpv)




myc= ROOT.TCanvas("myc", "myc", 1200, 600)
myc.SetMargin(0.06,0.03,0.1,0.075)
h_eff_phi.SetStats(ROOT.kFALSE)
h_eff_phi.SetTitle("Before taking out stuck-bit groups")
h_eff_phi.SetMaximum(1.)
h_eff_phi.SetFillColor(ROOT.kRed-9)
h_eff_phi.SetLineColor(ROOT.kRed+2)
h_eff_phi.GetXaxis().SetTitle("Sector number")
h_eff_phi.GetXaxis().SetNdivisions(134)
h_eff_phi.GetYaxis().SetTitle("Phi Efficiency")
h_eff_phi.GetYaxis().SetTitleOffset(0.7)
h_eff_phi.GetYaxis().SetTitleSize(0.04)
h_eff_phi.Draw("HIST")
myc.Update()
error_eff_phi = ROOT.TGraphAsymmErrors(h_eff_phi)
error_eff_phi.SetFillStyle(3001)
error_eff_phi.SetFillColor(ROOT.kRed+2)
error_eff_phi.Draw("2 SAME")
myc.Update()
myc.SaveAs("effplots/h_eff_phi_before.eps")
myc.SaveAs("effplots/h_eff_phi_before.png")
myc.Close()



myc= ROOT.TCanvas("myc", "myc", 1200, 600)
myc.SetMargin(0.06,0.03,0.1,0.075)
h_eff_phi_good.SetStats(ROOT.kFALSE)
h_eff_phi_good.SetTitle("After taking out stuck-bit groups")
h_eff_phi_good.SetMaximum(1.)
h_eff_phi_good.SetFillColor(ROOT.kRed-9)
h_eff_phi_good.SetLineColor(ROOT.kRed+2)
h_eff_phi_good.GetXaxis().SetTitle("Sector number")
h_eff_phi_good.GetXaxis().SetNdivisions(134)
h_eff_phi_good.GetYaxis().SetTitle("Phi Efficiency")
h_eff_phi_good.GetYaxis().SetTitleOffset(0.7)
h_eff_phi_good.GetYaxis().SetTitleSize(0.04)
h_eff_phi_good.Draw("HIST")
myc.Update()
error_eff_phi_good = ROOT.TGraphAsymmErrors(h_eff_phi_good)
error_eff_phi_good.SetFillStyle(3001)
error_eff_phi_good.SetFillColor(ROOT.kRed+2)
error_eff_phi_good.Draw("2 SAME")
myc.Update()
myc.SaveAs("effplots/h_eff_phi_after.eps")
myc.SaveAs("effplots/h_eff_phi_after.png")
myc.Close()




myc= ROOT.TCanvas("myc", "myc", 1200, 600)
myc.SetMargin(0.06,0.03,0.1,0.075)
h_eff_eta.SetStats(ROOT.kFALSE)
h_eff_eta.SetTitle("Before taking out stuck-bit groups")
h_eff_eta.SetMaximum(1.)
h_eff_eta.SetFillColor(ROOT.kRed-9)
h_eff_eta.SetLineColor(ROOT.kRed+2)
h_eff_eta.GetXaxis().SetTitle("Sector number")
h_eff_eta.GetXaxis().SetNdivisions(134)
h_eff_eta.GetYaxis().SetTitle("Eta Efficiency")
h_eff_eta.GetYaxis().SetTitleOffset(0.7)
h_eff_eta.GetYaxis().SetTitleSize(0.04)
h_eff_eta.Draw("HIST")
myc.Update()
error_eff_eta = ROOT.TGraphAsymmErrors(h_eff_eta)
error_eff_eta.SetFillStyle(3001)
error_eff_eta.SetFillColor(ROOT.kRed+2)
error_eff_eta.Draw("2 SAME")
myc.Update()
myc.SaveAs("effplots/h_eff_eta_before.eps")
myc.SaveAs("effplots/h_eff_eta_before.png")
myc.Close()




myc= ROOT.TCanvas("myc", "myc", 1200, 600)
myc.SetMargin(0.06,0.03,0.1,0.075)
h_eff_eta_good.SetStats(ROOT.kFALSE)
h_eff_eta_good.SetTitle("After taking out stuck-bit groups")
h_eff_eta_good.SetMaximum(1.)
h_eff_eta_good.SetFillColor(ROOT.kRed-9)
h_eff_eta_good.SetLineColor(ROOT.kRed+2)
h_eff_eta_good.GetXaxis().SetTitle("Sector number")
h_eff_eta_good.GetXaxis().SetNdivisions(134)
h_eff_eta_good.GetYaxis().SetTitle("Eta Efficiency")
h_eff_eta_good.GetYaxis().SetTitleOffset(0.7)
h_eff_eta_good.GetYaxis().SetTitleSize(0.04)
h_eff_eta_good.Draw("HIST")
myc.Update()
error_eff_eta_good = ROOT.TGraphAsymmErrors(h_eff_eta_good)
error_eff_eta_good.SetFillStyle(3001)
error_eff_eta_good.SetFillColor(ROOT.kRed+2)
error_eff_eta_good.Draw("2 SAME")
myc.Update()
myc.SaveAs("effplots/h_eff_eta_after.eps")
myc.SaveAs("effplots/h_eff_eta_after.png")
myc.Close()



mpvNormal = array('d', [204142.70196, 232075.891083, 233767.716497, 226983.023572, 220583.10026, 234038.481075, 222878.26697, 211570.639714, 227815.302947, 226511.244185, 257635.727314, 238341.326044, 250239.099514, 213946.647816, 224093.744515, 219528.867857, 242460.582762, 231631.027117, 243019.406121, 258399.562603, 250839.166909, 214398.493148, 253661.366986, 239931.277569, 232918.205663, 246838.058269, 243905.963793, 265950.917208, 217303.461464, 240918.817347, 252295.92043, 218965.503996, 261904.433838, 225453.41466, 256461.15088, 205431.16761, 230005.441539, 232219.25105, 227586.685733, 260077.745487, 237404.387734, 276772.232085, 192661.332161, 271693.091782, 255217.426118, 293516.968766, 302337.539441, 304430.74213, 206718.75855, 222053.147352, 252629.238521, 213064.425141, 245501.834213, 231670.101926, 214212.312161, 276699.840263, 193287.382262, 229867.690543, 212896.881433, 235437.043852, 239200.098617, 249541.31511, 240636.290254, 240383.465401, 236481.25323, 238137.906988, 268381.288555, 216563.897021, 232886.277688, 266021.685389, 251256.813683, 235981.816833, 224276.890393, 267174.939711, 236422.015658, 218393.454567, 266642.196186, 232042.749888, 231501.697549, 259182.851929, 240153.605206, 212935.664383, 249393.564891, 230129.579333, 236962.146312, 241074.398461, 257696.742678, 230740.345222, 209741.092325, 245915.235351, 238243.67241, 251606.162434, 253285.347111, 238928.603339, 233894.749186, 280335.996639]) 



mpv_60= array('d', [137178.578153, 121772.951352, 115237.838921, 134976.140866, 142049.259244, 123700.204978, 129235.274052, 128886.174245, 111435.707786, 110347.155649, 117405.923667, 140354.440314, 141613.664392, 141479.659757, 129441.04817, 129883.313652, 126223.848157, 145503.447843, 122832.610102, 133399.073309, 178049.158722, 153139.255342, 133573.413511])



mpv_140 = array('d', [64033.1373138, 76756.4203126, 91256.6826855, 73329.6457697, 85088.4959123])

mpv_200= array('d', [54026.7802723])


phi_eff_normal = array('d', [0.955459, 0.955422, 0.952262, 0.949602, 0.949742, 0.946711, 0.95642, 0.953026, 0.950872, 0.945653, 0.966093, 0.942312, 0.963356, 0.957151, 0.948637, 0.950326, 0.948846, 0.953585, 0.955658, 0.952675, 0.947474, 0.946079, 0.953303, 0.948549, 0.950336, 0.952549, 0.950735, 0.943987, 0.947272, 0.950826, 0.950144, 0.948467, 0.955224, 0.950434, 0.946573, 0.92769, 0.952844, 0.943201, 0.949207, 0.952786, 0.950305, 0.943514, 0.958819, 0.961829, 0.949684, 0.952407, 0.950265, 0.944711, 0.950729, 0.952547, 0.952499, 0.932976, 0.954094, 0.946989, 0.938708, 0.936826, 0.94026, 0.938082, 0.953048, 0.953237, 0.948074, 0.943364, 0.946323, 0.952856, 0.914576, 0.949781, 0.945524, 0.942248, 0.950398, 0.950301, 0.93625, 0.959096, 0.948701, 0.955182, 0.947009, 0.952742, 0.9519, 0.956807, 0.953161, 0.951836, 0.925773, 0.948997, 0.952053, 0.902753, 0.957498, 0.954952, 0.953798, 0.943518, 0.945503, 0.947209, 0.965722, 0.934777, 0.954504, 0.954084, 0.955487, 0.951165])


phi_eff_lower60 = array('d', [0.936717, 0.917325, 0.91075, 0.672517, 0.950835, 0.936516, 0.923396, 0.917092, 0.908667, 0.866226, 0.975763, 0.943868, 0.934403, 0.931029, 0.907909, 0.931639, 0.574549, 0.950117, 0.925625, 0.929702, 0.93953, 0.888255, 0.918104])


phi_eff_lower140 = array('d', [0.718337, 0.829644, 0.921973, 0.86068, 0.930869])

phi_eff_lower200 = array('d', [0.883969])


myc= ROOT.TCanvas("myc", "myc", 1400, 700)
myc.SetMargin(0.06,0.05,0.1,0.05)
phi_mg = ROOT.TMultiGraph()

graphNormal = ROOT.TGraph(len(phi_eff_normal), mpvNormal, phi_eff_normal)
graphLower60 = ROOT.TGraph(len(phi_eff_lower60), mpv_60, phi_eff_lower60)
graphLower140 = ROOT.TGraph(len(phi_eff_lower140), mpv_140, phi_eff_lower140)
graphLower200 = ROOT.TGraph(len(phi_eff_lower200), mpv_200, phi_eff_lower200)


graphNormal.SetMarkerColor(kBlue)
graphLower60.SetMarkerColor(kRed)
graphLower140.SetMarkerColor(kBlack)
graphLower200.SetMarkerColor(kGreen)

graphNormal.SetMarkerStyle(20)
graphLower60.SetMarkerStyle(20) 
graphLower140.SetMarkerStyle(20)
graphLower200.SetMarkerStyle(20)

phi_mg.Add(graphNormal)
phi_mg.Add(graphLower60)
phi_mg.Add(graphLower140)
phi_mg.Add(graphLower200)

leg = ROOT.TLegend(0.78, 0.16, 0.94, 0.35)
leg.SetTextSize(0.04)
leg.AddEntry(graphNormal, "Normal layers", "p")
leg.AddEntry(graphLower60, "60V lower", "p")
leg.AddEntry(graphLower140, "140V lower", "p")
leg.AddEntry(graphLower200, "200V lower", "p")

phi_mg.Draw("AP") 
leg.Draw()

phi_mg.GetXaxis().SetTitle("Layer MPV")
phi_mg.GetXaxis().SetTitleSize(0.04)
#phi_mg.GetXaxis().SetLabelSize(0.04) 
#phi_mg.GetXaxis().SetNdivisions(66)
phi_mg.GetYaxis().SetTitle("Layer Phi Efficiency")
phi_mg.GetYaxis().SetTitleOffset(0.7)
phi_mg.GetYaxis().SetTitleSize(0.04)
#phi_mg.GetYaxis().SetRangeUser(0.295,1.15)
 

gPad.Modified()
myc.Update()
myc.SaveAs("effplots/PhiEff_vs_MPV.eps")
myc.SaveAs("effplots/PhiEff_vs_MPV.png")
myc.Close()





eta_eff_normal = array('d', [0.950853, 0.931425, 0.928411, 0.912871, 0.904825, 0.839082, 0.933462, 0.930141, 0.928498, 0.92244, 0.96966, 0.924042, 0.940922, 0.936956, 0.951617, 0.965532, 0.906228, 0.932996, 0.930195, 0.926594, 0.923357, 0.9139, 0.910195, 0.905394, 0.931471, 0.929911, 0.958221, 0.920842, 0.912269, 0.907867, 0.906219, 0.931982, 0.961356, 0.924742, 0.92322, 0.915787, 0.911915, 0.951156, 0.929931, 0.955112, 0.942612, 0.915413, 0.924038, 0.920117, 0.932585, 0.93077, 0.925688, 0.925034, 0.886812, 0.895882, 0.963704, 0.936612, 0.957621, 0.922513, 0.914672, 0.912268, 0.938693, 0.901227, 0.962448, 0.93221, 0.92799, 0.924147, 0.907771, 0.932595, 0.930309, 0.926525, 0.952705, 0.842632, 0.91417, 0.911654, 0.905415, 0.939066, 0.955549, 0.933656, 0.927823, 0.91482, 0.971368, 0.960057, 0.931804, 0.928637, 0.923777, 0.915959, 0.960666, 0.907046, 0.962657, 0.935108, 0.956719, 0.924566, 0.93248, 0.958909, 0.944821, 0.923457, 0.963448, 0.935503, 0.930294, 0.926822]) 

eta_eff_lower60 = array('d', [0.948467, 0.89454, 0.928609, 0.914254, 0.927189, 0.948872, 0.896766, 0.894675, 0.895273, 0.87356, 0.943305, 0.924731, 0.91872, 0.909179, 0.897387, 0.956109, 0.904342, 0.936113, 0.930517, 0.923921, 0.905354, 0.902059, 0.916959]) 

eta_eff_lower140 = array('d', [0.821693, 0.89313, 0.906436, 0.900503, 0.900479])

eta_eff_lower200 = array('d', [0.857279])




myc= ROOT.TCanvas("myc", "myc", 1400, 700)
myc.SetMargin(0.07,0.05,0.1,0.075)
eta_mg = ROOT.TMultiGraph()

graphNormalEta = ROOT.TGraph(len(eta_eff_normal), mpvNormal, eta_eff_normal)
graphLower60Eta = ROOT.TGraph(len(eta_eff_lower60), mpv_60, eta_eff_lower60)
graphLower140Eta = ROOT.TGraph(len(eta_eff_lower140), mpv_140, eta_eff_lower140)
graphLower200Eta = ROOT.TGraph(len(eta_eff_lower200), mpv_200, eta_eff_lower200)

graphNormalEta.SetMarkerColor(kBlue)
graphLower60Eta.SetMarkerColor(kRed)
graphLower140Eta.SetMarkerColor(kBlack)
graphLower200Eta.SetMarkerColor(kGreen)

graphNormalEta.SetMarkerStyle(20)
graphLower60Eta.SetMarkerStyle(20) 
graphLower140Eta.SetMarkerStyle(20)
graphLower200Eta.SetMarkerStyle(20)

eta_mg.Add(graphNormalEta)
eta_mg.Add(graphLower60Eta)
eta_mg.Add(graphLower140Eta)
eta_mg.Add(graphLower200Eta)

leg = ROOT.TLegend(0.78, 0.16, 0.94, 0.35)
leg.SetTextSize(0.04)
leg.AddEntry(graphNormalEta, "Normal layers", "p")
leg.AddEntry(graphLower60Eta, "60V lower", "p")
leg.AddEntry(graphLower140Eta, "140V lower", "p")
leg.AddEntry(graphLower200Eta, "200V lower", "p")

eta_mg.Draw("AP") 
leg.Draw()

eta_mg.SetTitle("All Layers")
eta_mg.GetXaxis().SetTitle("Layer MPV")
eta_mg.GetXaxis().SetTitleSize(0.04)
eta_mg.GetYaxis().SetTitle("Layer Eta Efficiency")
eta_mg.GetYaxis().SetTitleOffset(0.85)
eta_mg.GetYaxis().SetTitleSize(0.04)
eta_mg.GetYaxis().SetRangeUser(0.85,1.0)

gPad.Modified()
myc.Update()
myc.SaveAs("effplots/EtaEff_vs_MPV.eps")
myc.SaveAs("effplots/EtaEff_vs_MPV.png")
myc.Close()













