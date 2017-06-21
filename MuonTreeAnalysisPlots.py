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

inputname = "MuonTreeAnalysis.root"
inputfile = ROOT.TFile.Open( inputname) 

#h_hit_channel_phi = inputfile.Get("hit_channel_phi") 
h_hit_channel_phi_oneLayer = inputfile.Get("hit_channel_phi_oneLayer") 
#h2_hit_channel_phi_byLayers = inputfile.Get("hit_channel_phi_byLayers")
#h2_eff_channel = inputfile.Get("eff_channel") 

'''
myc= ROOT.TCanvas("myc", "myc", 800, 600)
myc.SetMargin(0.12,0.05,0.1,0.075)
h_hit_channel_phi.SetStats(ROOT.kFALSE)
h_hit_channel_phi.SetTitle("Muon Tree")
h_hit_channel_phi.GetXaxis().SetTitle("Phi Channel number")
h_hit_channel_phi.GetYaxis().SetTitle("Number of hits")
h_hit_channel_phi.GetYaxis().SetTitleOffset(1.4)
h_hit_channel_phi.GetYaxis().SetTitleSize(0.04)
h_hit_channel_phi.Draw()
myc.Update()
myc.SaveAs("muonPlots/h_hit_channel_phi.eps")
myc.SaveAs("muonPlots/h_hit_channel_phi.png")
myc.Close()
'''


myc= ROOT.TCanvas("myc", "myc", 800, 600)
myc.SetMargin(0.12,0.05,0.1,0.075)
h_hit_channel_phi_oneLayer.SetStats(ROOT.kFALSE)
h_hit_channel_phi_oneLayer.SetTitle("Muon Tree Sector 3, Layer 2")
h_hit_channel_phi_oneLayer.GetXaxis().SetTitle("Phi Channel number")
h_hit_channel_phi_oneLayer.GetYaxis().SetTitle("Number of hits")
h_hit_channel_phi_oneLayer.GetYaxis().SetTitleOffset(1.4)
h_hit_channel_phi_oneLayer.GetYaxis().SetTitleSize(0.04)
h_hit_channel_phi_oneLayer.Draw()
myc.Update()
myc.SaveAs("muonPlots/h_hit_channel_phi_oneLayer_sector3_layer2.eps")
myc.SaveAs("muonPlots/h_hit_channel_phi_oneLayer_sector3_layer2.png")
myc.Close()



'''
myc= ROOT.TCanvas("myc", "myc", 1000, 800)
myc.SetMargin(0.10,0.15,0.1,0.075)
h2_hit_channel_phi_byLayers.SetStats(ROOT.kFALSE)
h2_hit_channel_phi_byLayers.SetTitle("Muon Tree")
h2_hit_channel_phi_byLayers.GetXaxis().SetTitle("Phi Channel number")
h2_hit_channel_phi_byLayers.GetYaxis().SetTitle("Sector number")
h2_hit_channel_phi_byLayers.GetYaxis().SetTitleOffset(0.9)
h2_hit_channel_phi_byLayers.GetYaxis().SetTitleSize(0.04)
h2_hit_channel_phi_byLayers.GetZaxis().SetTitle("Hits")
h2_hit_channel_phi_byLayers.GetZaxis().SetTitleOffset(1.3)
h2_hit_channel_phi_byLayers.Draw("COLZ0")
myc.Update()
myc.SaveAs("muonPlots/h2_hit_channel_phi_byLayers.eps")
myc.SaveAs("muonPlots/h2_hit_channel_phi_byLayers.png")
myc.Close()
'''



'''
myc= ROOT.TCanvas("myc", "myc", 1000, 800)
myc.SetMargin(0.10,0.15,0.1,0.075)
h2_eff_channel.SetStats(ROOT.kFALSE)
h2_eff_channel.SetTitle("")
h2_eff_channel.GetXaxis().SetTitle("Phi Channel Number")
h2_eff_channel.GetYaxis().SetTitle("Sector number")
h2_eff_channel.GetYaxis().SetTitleOffset(0.9)
h2_eff_channel.GetYaxis().SetTitleSize(0.04)
h2_eff_channel.GetZaxis().SetTitle("Hits")
h2_eff_channel.GetZaxis().SetTitleOffset(1.3)
h2_eff_channel.Draw("COLZ0")
myc.Update()
myc.SaveAs("muonPlots/h2_eff_channel.eps")
myc.SaveAs("muonPlots/h2_eff_channel.png")
myc.Close()
'''







