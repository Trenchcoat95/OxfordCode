import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from array import array
import os
import sys
import ROOT
from ROOT import TVectorD, TMatrix, TMath, TVector3, TGraphErrors, TFile, TTree, gRandom, gPad, gROOT, gVirtualX, kTRUE, kRed, TProfile, gStyle,  TFile, gSystem
import sys 
import os
sys.path.append(os.path.abspath("/home/federico/Documents/Universita/Federico_2020-2021/Aliwork/fastMCKalman/fastMCKalman/MC/"))
#from fastSimulation import *


def Plot_prof_InRot(tree,ParamType1,ParamType2,Var,VarName,BinsRangeX,BinsRangeY,BinsRangeVar,legenda1,legenda2,RangeUserYsigma,RangeUserYmean,savesigma,savemean,extracond):
    tree.Draw("(fParamMC[0].GetP()-"+ParamType1+".GetP())/fParamMC[0].GetP():"+Var+">>hpkGAr("+BinsRangeX+","+BinsRangeY+")",ParamType1+".fP[4]!=0!=0&&fParamMC[0].fP[4]!=0"+extracond,"colz")
    hpkGAr = ROOT.gPad.GetPrimitive("hpkGAr")
    hpkGAr.FitSlicesY()
    hpkGAr_mean = ROOT.gDirectory.Get("hpkGAr_1")
    hpkGAr_sigma = ROOT.gDirectory.Get("hpkGAr_2")
    hpkGAr_mean.SetTitle(""+VarName+" #mu((p_{reco}-p_{MC})/p_{MC})")
    hpkGAr_sigma.SetTitle(""+VarName+" #sigma((p_{reco}-p_{MC})/p_{MC})")

    tree.Draw("(fParamMC[0].GetP()-"+ParamType2+".GetP())/fParamMC[0].GetP():"+Var+">>hpk("+BinsRangeX+","+BinsRangeY+")",ParamType2+".fP[4]!=0&&fParamMC[0].fP[4]!=0"+extracond,"colz")
    hpk = ROOT.gPad.GetPrimitive("hpk")
    hpk.FitSlicesY()
    hpk_mean = ROOT.gDirectory.Get("hpk_1")
    hpk_mean.SetTitle(""+VarName+" #mu((p_{reco}-p_{MC})/p_{MC}) (MeV/c)")
    hpk_sigma = ROOT.gDirectory.Get("hpk_2")
    hpk_sigma.SetTitle(""+VarName+" #sigma((p_{reco}-p_{MC})/p_{MC}) (MeV/c)")

    tree.Draw(Var+">>hlen("+BinsRangeVar+")",ParamType2+".fP[4]!=0&&fParamMC[0].fP[4]!=0"+extracond,"colz")
    hlen = ROOT.gPad.GetPrimitive("hlen")
    hlen.SetTitle(""+VarName+" n_{ev}")

    hpkGAr_sigma.SetMarkerStyle(20)
    hpk_sigma.SetMarkerStyle(21)
    hpk_sigma.SetMarkerColor(ROOT.kRed)
    hpkGAr_mean.SetMarkerStyle(20)
    hpk_mean.SetMarkerStyle(21)
    hpk_mean.SetMarkerColor(ROOT.kRed)

    ROOT.gStyle.SetOptStat(0)
    hqq = ROOT.TCanvas("hqq","hqq",800,800)
    hqq.Draw()
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.35, 1, 1.0)
    pad1.SetBottomMargin(0.015) # Upper and lower plot are joined
    pad1.SetRightMargin(0.05) # Upper and lower plot are joined
    pad1.SetLeftMargin(0.12) # Upper and lower plot are joined
    pad1.SetGridx()         # Vertical grid
    pad1.Draw()             # Draw the upper pad: pad1
    pad1.cd()               # pad1 becomes the current pad
    hpk_sigma.Draw()
    hpk_sigma.GetYaxis().SetLabelSize(0.03)
    hpk_sigma.GetXaxis().SetLabelSize(0.0)
    hpk_sigma.GetYaxis().SetRangeUser(RangeUserYsigma[0],RangeUserYsigma[1])
    axis = ROOT.TGaxis( -5, 20, -5, 220, 20,220,510,"")
    axis.SetLabelFont(43) # Absolute font size in pixel (precision 3)
    axis.SetLabelSize(15)
    axis.Draw()

    hpkGAr_sigma.Draw("same")
    legend = ROOT.TLegend(0.6,0.72,0.92,0.88)
    #legend.SetBorderSize(0)
    legend.AddEntry(hpkGAr_sigma,"ND-GAr KF Res","pl")
    legend.AddEntry(hpk_sigma,"fastMCKalman KF Res","pl")
    legend.Draw()
    hqq.Draw()

    hqq.cd()          # Go back to the main canvas before defining pad2
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0.05)
    pad2.SetBottomMargin(0.22)
    pad2.SetRightMargin(0.05) # Upper and lower plot are joined
    pad2.SetLeftMargin(0.12) # Upper and lower plot are joined
    pad2.SetGridx() # vertical grid
    pad2.Draw()
    pad2.cd()       # pad2 becomes the current pad

    hlen.GetYaxis().SetTitleSize(20)
    hlen.GetYaxis().SetTitleFont(43)
    hlen.GetYaxis().SetTitleOffset(1.55)
    hlen.GetYaxis().SetLabelFont(43) # Absolute font size in pixel (precision 3)
    hlen.GetYaxis().SetLabelSize(15)
    hlen.GetYaxis().SetNdivisions(505)

    hlen.GetXaxis().SetTitleSize(20)
    hlen.GetXaxis().SetTitleFont(43)
    hlen.GetXaxis().SetTitleOffset(1)
    hlen.GetXaxis().SetLabelFont(43) # Absolute font size in pixel (precision 3)
    hlen.GetXaxis().SetLabelSize(15)
    hlen.Draw()

    hqq.SaveAs(savesigma)

    ROOT.gStyle.SetOptStat(0)
    hqq2 = ROOT.TCanvas("hqq2","hqq2",800,800)
    pad1m = ROOT.TPad("pad1m", "pad1m", 0, 0.35, 1, 1.0)
    pad1m.SetBottomMargin(0.015) # Upper and lower plot are joined
    pad1m.SetRightMargin(0.05) # Upper and lower plot are joined
    pad1m.SetLeftMargin(0.12) # Upper and lower plot are joined
    pad1m.SetGridx()         # Vertical grid
    pad1m.Draw()             # Draw the upper pad: pad1
    pad1m.cd()               # pad1 becomes the current pad
    hpkGAr_mean.Draw()
    hpkGAr_mean.GetYaxis().SetRangeUser(RangeUserYmean[0],RangeUserYmean[1])
    hpkGAr_mean.GetYaxis().SetLabelSize(0.03)
    hpkGAr_mean.GetXaxis().SetLabelSize(0.0)
    axism = ROOT.TGaxis( -5, 20, -5, 220, 20,220,510,"")
    axism.SetLabelFont(43) # Absolute font size in pixel (precision 3)
    axism.SetLabelSize(15)
    axism.Draw()

    hpk_mean.Draw("same")
    legendm = ROOT.TLegend(0.6,0.72,0.92,0.88)
    #legend.SetBorderSize(0)
    legendm.AddEntry(hpkGAr_sigma,"ND-GAr KF Bias","pl")
    legendm.AddEntry(hpk_sigma,"fastMCKalman KF Bias","pl")
    legendm.Draw()
    hqq2.Draw()

    hqq2.cd()          # Go back to the main canvas before defining pad2
    pad2m = ROOT.TPad("pad2m", "pad2m", 0, 0.05, 1, 0.3)
    pad2m.SetTopMargin(0.05)
    pad2m.SetBottomMargin(0.22)
    pad2m.SetRightMargin(0.05) # Upper and lower plot are joined
    pad2m.SetLeftMargin(0.12) # Upper and lower plot are joined
    pad2m.SetGridx() # vertical grid
    pad2m.Draw()
    pad2m.cd()       # pad2 becomes the current pad

    hlen.Draw()

    hqq2.SaveAs(savemean)

def Plot_prof(tree1,tree2,ParamTypeMC,ParamType1,ParamType2,ParamFunc,ParamName,Var,VarName,BinsRangeX,BinsRangeY,BinsRangeVar,legenda1,legenda2,RangeUserYsigma,RangeUserYmean,savesigma,savemean,extracond,frac=True,extraform="",drawsecond = True):
    
    resInstr = "("+ParamTypeMC+ParamFunc+"-"+ParamType1+ParamFunc+")" 
    if(frac): resInstr +=  "/"+ParamTypeMC+ParamFunc
    tree1.Draw(resInstr+":"+Var+">>hpkGAr("+BinsRangeX+","+BinsRangeY+")",extracond,"colz")
    hpkGAr = ROOT.gPad.GetPrimitive("hpkGAr")
    hpkGAr.FitSlicesY()
    hpkGAr_mean = ROOT.gDirectory.Get("hpkGAr_1")
    hpkGAr_sigma = ROOT.gDirectory.Get("hpkGAr_2")
    namemean = ""+VarName+" #mu("+ParamName+"_{reco}-"+ParamName+"_{MC})"
    if(frac): namemean=""+VarName+" #mu(("+ParamName+"_{reco}-"+ParamName+"_{MC})/"+ParamName+"_{MC})"
    hpkGAr_mean.SetTitle(namemean)
    namesigma = ""+VarName+" #sigma("+ParamName+"_{reco}-"+ParamName+"_{MC})"
    if(frac): namesigma = ""+VarName+" #sigma(("+ParamName+"_{reco}-"+ParamName+"_{MC})/"+ParamName+"_{MC})"
    hpkGAr_sigma.SetTitle(namesigma)

    resInstr2 = "("+ParamTypeMC+ParamFunc+"-"+ParamType2+ParamFunc+")" 
    if(frac): resInstr2 +=  "/"+ParamTypeMC+ParamFunc+extraform
    tree2.Draw(resInstr2+":"+Var+">>hpk("+BinsRangeX+","+BinsRangeY+")",extracond,"colz")
    hpk = ROOT.gPad.GetPrimitive("hpk")
    hpk.FitSlicesY()
    hpk_mean = ROOT.gDirectory.Get("hpk_1")
    hpk_mean.SetTitle(namemean)
    hpk_sigma = ROOT.gDirectory.Get("hpk_2")
    hpk_sigma.SetTitle(namesigma)

    tree1.Draw(Var+">>hlen("+BinsRangeVar+")",extracond,"colz")
    hlen = ROOT.gPad.GetPrimitive("hlen")
    hlen.SetTitle(""+VarName+" n_{ev}")

    tree2.Draw(Var+">>hlen2("+BinsRangeVar+")",extracond,"colz")
    hlen2 = ROOT.gPad.GetPrimitive("hlen2")
    hlen2.SetTitle(""+VarName+" n_{ev}")

    hpkGAr_sigma.SetMarkerStyle(20)
    hpk_sigma.SetMarkerStyle(21)
    hpk_sigma.SetMarkerColor(ROOT.kRed)
    hpkGAr_mean.SetMarkerStyle(20)
    hpk_mean.SetMarkerStyle(21)
    hpk_mean.SetMarkerColor(ROOT.kRed)

    ROOT.gStyle.SetOptStat(0)
    hqq = ROOT.TCanvas("hqq","hqq",800,800)
    hqq.Draw()
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.35, 1, 1.0)
    pad1.SetBottomMargin(0.015) # Upper and lower plot are joined
    pad1.SetRightMargin(0.05) # Upper and lower plot are joined
    pad1.SetLeftMargin(0.12) # Upper and lower plot are joined
    pad1.SetGridx()         # Vertical grid
    pad1.Draw()             # Draw the upper pad: pad1
    pad1.cd()               # pad1 becomes the current pad
    hpk_sigma.Draw()
    hpk_sigma.GetYaxis().SetLabelSize(0.03)
    hpk_sigma.GetXaxis().SetLabelSize(0.0)
    hpk_sigma.GetYaxis().SetRangeUser(RangeUserYsigma[0],RangeUserYsigma[1])
    axis = ROOT.TGaxis( -5, 20, -5, 220, 20,220,510,"")
    axis.SetLabelFont(43) # Absolute font size in pixel (precision 3)
    axis.SetLabelSize(15)
    axis.Draw()

    if(drawsecond): hpkGAr_sigma.Draw("same")
    legend = ROOT.TLegend(0.6,0.72,0.92,0.88)
    #legend.SetBorderSize(0)
    if(drawsecond): legend.AddEntry(hpkGAr_sigma,legenda1+" Resolution","pl")
    legend.AddEntry(hpk_sigma,legenda2+" Resolution","pl")
    legend.Draw()
    hqq.Draw()

    hqq.cd()          # Go back to the main canvas before defining pad2
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0.05)
    pad2.SetBottomMargin(0.22)
    pad2.SetRightMargin(0.05) # Upper and lower plot are joined
    pad2.SetLeftMargin(0.12) # Upper and lower plot are joined
    pad2.SetGridx() # vertical grid
    pad2.Draw()
    pad2.cd()       # pad2 becomes the current pad

    hlen2.GetYaxis().SetTitleSize(20)
    hlen2.GetYaxis().SetTitleFont(43)
    hlen2.GetYaxis().SetTitleOffset(1.55)
    hlen2.GetYaxis().SetLabelFont(43) # Absolute font size in pixel (precision 3)
    hlen2.GetYaxis().SetLabelSize(15)
    hlen2.GetYaxis().SetNdivisions(505)

    hlen2.GetXaxis().SetTitleSize(20)
    hlen2.GetXaxis().SetTitleFont(43)
    hlen2.GetXaxis().SetTitleOffset(1)
    hlen2.GetXaxis().SetLabelFont(43) # Absolute font size in pixel (precision 3)
    hlen2.GetXaxis().SetLabelSize(15)
    hlen2.Draw()
    hlen2.SetLineColor(kRed)
    hlen.Draw("same")

    hqq.SaveAs(savesigma)

    ROOT.gStyle.SetOptStat(0)
    hqq2 = ROOT.TCanvas("hqq2","hqq2",800,800)
    pad1m = ROOT.TPad("pad1m", "pad1m", 0, 0.35, 1, 1.0)
    pad1m.SetBottomMargin(0.015) # Upper and lower plot are joined
    pad1m.SetRightMargin(0.05) # Upper and lower plot are joined
    pad1m.SetLeftMargin(0.12) # Upper and lower plot are joined
    pad1m.SetGridx()         # Vertical grid
    pad1m.Draw()             # Draw the upper pad: pad1
    pad1m.cd()               # pad1 becomes the current pad
    hpk_mean.Draw()
    hpk_mean.GetYaxis().SetRangeUser(RangeUserYmean[0],RangeUserYmean[1])
    hpk_mean.GetYaxis().SetLabelSize(0.03)
    hpk_mean.GetXaxis().SetLabelSize(0.0)
    axism = ROOT.TGaxis( -5, 20, -5, 220, 20,220,510,"")
    axism.SetLabelFont(43) # Absolute font size in pixel (precision 3)
    axism.SetLabelSize(15)
    axism.Draw()

    if(drawsecond): hpkGAr_mean.Draw("same")
    legendm = ROOT.TLegend(0.6,0.72,0.92,0.88)
    #legend.SetBorderSize(0)
    if(drawsecond): legendm.AddEntry(hpkGAr_sigma,legenda1+" Bias","pl")
    legendm.AddEntry(hpk_sigma,legenda2+" Bias","pl")
    legendm.Draw()
    #line = ROOT.TLine(0, 0, 700, 0)
    #line.Draw()
    hqq2.Draw()

    hqq2.cd()          # Go back to the main canvas before defining pad2
    pad2m = ROOT.TPad("pad2m", "pad2m", 0, 0.05, 1, 0.3)
    pad2m.SetTopMargin(0.05)
    pad2m.SetBottomMargin(0.22)
    pad2m.SetRightMargin(0.05) # Upper and lower plot are joined
    pad2m.SetLeftMargin(0.12) # Upper and lower plot are joined
    pad2m.SetGridx() # vertical grid
    pad2m.Draw()
    pad2m.cd()       # pad2 becomes the current pad

    hlen2.Draw()
    hlen2.SetLineColor(kRed)
    hlen.Draw("same")

    hqq2.SaveAs(savemean)


def Plot_prof3(tree1,tree2,tree3,ParamTypeMC,ParamType1,ParamType2,ParamType3,ParamFunc,ParamName,Var,VarName,BinsRangeX,BinsRangeY,BinsRangeVar,legenda1,legenda2,legenda3,RangeUserYsigma,RangeUserYmean,savesigma,savemean,extracond1,extracond2,extracond3,frac=True,extraform="",drawsecond = True):
    
    resInstr = "("+ParamTypeMC+ParamFunc+"-"+ParamType1+ParamFunc+")" 
    if(frac): resInstr +=  "/"+ParamTypeMC+ParamFunc
    tree1.Draw(resInstr+":"+Var+">>hpkGAr("+BinsRangeX+","+BinsRangeY+")",extracond1,"colz")
    hpkGAr = ROOT.gPad.GetPrimitive("hpkGAr")
    hpkGAr.FitSlicesY()
    hpkGAr_mean = ROOT.gDirectory.Get("hpkGAr_1")
    hpkGAr_sigma = ROOT.gDirectory.Get("hpkGAr_2")
    namemean = ""+VarName+" #mu("+ParamName+"_{reco}-"+ParamName+"_{MC})"
    if(frac): namemean=""+VarName+" #mu(("+ParamName+"_{reco}-"+ParamName+"_{MC})/"+ParamName+"_{MC})"
    hpkGAr_mean.SetTitle(namemean)
    namesigma = ""+VarName+" #sigma("+ParamName+"_{reco}-"+ParamName+"_{MC})"
    if(frac): namesigma = ""+VarName+" #sigma(("+ParamName+"_{reco}-"+ParamName+"_{MC})/"+ParamName+"_{MC})"
    hpkGAr_sigma.SetTitle(namesigma)

    resInstr2 = "("+ParamTypeMC+ParamFunc+"-"+ParamType2+ParamFunc+")" 
    if(frac): resInstr2 +=  "/"+ParamTypeMC+ParamFunc+extraform
    tree2.Draw(resInstr2+":"+Var+">>hpk("+BinsRangeX+","+BinsRangeY+")",extracond2,"colz")
    hpk = ROOT.gPad.GetPrimitive("hpk")
    hpk.FitSlicesY()
    hpk_mean = ROOT.gDirectory.Get("hpk_1")
    hpk_mean.SetTitle(namemean)
    hpk_sigma = ROOT.gDirectory.Get("hpk_2")
    hpk_sigma.SetTitle(namesigma)

    resInstr3 = "("+ParamTypeMC+ParamFunc+"-"+ParamType3+ParamFunc+")" 
    if(frac): resInstr3 +=  "/"+ParamTypeMC+ParamFunc+extraform
    tree3.Draw(resInstr3+":"+Var+">>hpk2("+BinsRangeX+","+BinsRangeY+")",extracond3,"colz")
    hpk2 = ROOT.gPad.GetPrimitive("hpk2")
    hpk2.FitSlicesY()
    hpk2_mean = ROOT.gDirectory.Get("hpk2_1")
    hpk2_mean.SetTitle(namemean)
    hpk2_sigma = ROOT.gDirectory.Get("hpk2_2")
    hpk2_sigma.SetTitle(namesigma)

    tree1.Draw(Var+">>hlen("+BinsRangeVar+")",extracond1,"colz")
    hlen = ROOT.gPad.GetPrimitive("hlen")
    hlen.SetTitle(""+VarName+" n_{ev}")

    tree2.Draw(Var+">>hlen2("+BinsRangeVar+")",extracond2,"colz")
    hlen2 = ROOT.gPad.GetPrimitive("hlen2")
    hlen2.SetTitle(""+VarName+" n_{ev}")

    tree3.Draw(Var+">>hlen3("+BinsRangeVar+")",extracond3,"colz")
    hlen3 = ROOT.gPad.GetPrimitive("hlen3")
    hlen3.SetTitle(""+VarName+" n_{ev}")

    

    hpkGAr_sigma.SetMarkerStyle(20)
    hpk_sigma.SetMarkerStyle(21)
    hpk_sigma.SetMarkerColor(ROOT.kRed)
    hpk2_sigma.SetMarkerStyle(22)
    hpk2_sigma.SetMarkerColor(ROOT.kBlue)

    hpkGAr_mean.SetMarkerStyle(20)
    hpk_mean.SetMarkerStyle(21)
    hpk_mean.SetMarkerColor(ROOT.kRed)
    hpk2_mean.SetMarkerStyle(22)
    hpk2_mean.SetMarkerColor(ROOT.kBlue)

    ROOT.gStyle.SetOptStat(0)
    hqq = ROOT.TCanvas("hqq","hqq",800,800)
    hqq.Draw()
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.35, 1, 1.0)
    pad1.SetBottomMargin(0.015) # Upper and lower plot are joined
    pad1.SetRightMargin(0.05) # Upper and lower plot are joined
    pad1.SetLeftMargin(0.12) # Upper and lower plot are joined
    pad1.SetGridx()         # Vertical grid
    pad1.Draw()             # Draw the upper pad: pad1
    pad1.cd()               # pad1 becomes the current pad
    hpk_sigma.Draw()
    hpk_sigma.GetYaxis().SetLabelSize(0.03)
    hpk_sigma.GetXaxis().SetLabelSize(0.0)
    hpk_sigma.GetYaxis().SetRangeUser(RangeUserYsigma[0],RangeUserYsigma[1])
    axis = ROOT.TGaxis( -5, 20, -5, 220, 20,220,510,"")
    axis.SetLabelFont(43) # Absolute font size in pixel (precision 3)
    axis.SetLabelSize(15)
    axis.Draw()

    if(drawsecond): hpkGAr_sigma.Draw("same")
    if(drawsecond): hpk2_sigma.Draw("same")
    legend = ROOT.TLegend(0.6,0.72,0.92,0.88)
    #legend.SetBorderSize(0)
    if(drawsecond): legend.AddEntry(hpkGAr_sigma,legenda1+" Resolution","pl")
    if(drawsecond): legend.AddEntry(hpk2_sigma,legenda3+" Resolution","pl")
    legend.AddEntry(hpk_sigma,legenda2+" Resolution","pl")
    legend.Draw()
    hqq.Draw()

    hqq.cd()          # Go back to the main canvas before defining pad2
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0.05)
    pad2.SetBottomMargin(0.22)
    pad2.SetRightMargin(0.05) # Upper and lower plot are joined
    pad2.SetLeftMargin(0.12) # Upper and lower plot are joined
    pad2.SetGridx() # vertical grid
    pad2.Draw()
    pad2.cd()       # pad2 becomes the current pad

    hlen2.GetYaxis().SetTitleSize(20)
    hlen2.GetYaxis().SetTitleFont(43)
    hlen2.GetYaxis().SetTitleOffset(1.55)
    hlen2.GetYaxis().SetLabelFont(43) # Absolute font size in pixel (precision 3)
    hlen2.GetYaxis().SetLabelSize(15)
    hlen2.GetYaxis().SetNdivisions(505)

    hlen2.GetXaxis().SetTitleSize(20)
    hlen2.GetXaxis().SetTitleFont(43)
    hlen2.GetXaxis().SetTitleOffset(1)
    hlen2.GetXaxis().SetLabelFont(43) # Absolute font size in pixel (precision 3)
    hlen2.GetXaxis().SetLabelSize(15)

    hlen2.SetLineColor(ROOT.kRed)
    hlen2.Draw()
    hlen3.SetLineColor(ROOT.kBlue)
    hlen3.Draw("same")
    hlen.SetLineColor(ROOT.kBlack)
    hlen.Draw("same")

    hqq.SaveAs(savesigma)

    ROOT.gStyle.SetOptStat(0)
    hqq2 = ROOT.TCanvas("hqq2","hqq2",800,800)
    pad1m = ROOT.TPad("pad1m", "pad1m", 0, 0.35, 1, 1.0)
    pad1m.SetBottomMargin(0.015) # Upper and lower plot are joined
    pad1m.SetRightMargin(0.05) # Upper and lower plot are joined
    pad1m.SetLeftMargin(0.12) # Upper and lower plot are joined
    pad1m.SetGridx()         # Vertical grid
    pad1m.Draw()             # Draw the upper pad: pad1
    pad1m.cd()               # pad1 becomes the current pad
    hpk_mean.Draw()
    hpk_mean.GetYaxis().SetRangeUser(RangeUserYmean[0],RangeUserYmean[1])
    hpk_mean.GetYaxis().SetLabelSize(0.03)
    hpk_mean.GetXaxis().SetLabelSize(0.0)
    axism = ROOT.TGaxis( -5, 20, -5, 220, 20,220,510,"")
    axism.SetLabelFont(43) # Absolute font size in pixel (precision 3)
    axism.SetLabelSize(15)
    axism.Draw()

    if(drawsecond): hpkGAr_mean.Draw("same")
    if(drawsecond): hpk2_mean.Draw("same")
    legendm = ROOT.TLegend(0.6,0.72,0.92,0.88)
    #legend.SetBorderSize(0)
    if(drawsecond): legendm.AddEntry(hpkGAr_sigma,legenda1+" Bias","pl")
    if(drawsecond): legendm.AddEntry(hpk2_sigma,legenda3+" Bias","pl")
    legendm.AddEntry(hpk_sigma,legenda2+" Bias","pl")
    legendm.Draw()
    #line = ROOT.TLine(0, 0, 700, 0)
    #line.Draw()
    hqq2.Draw()

    hqq2.cd()          # Go back to the main canvas before defining pad2
    pad2m = ROOT.TPad("pad2m", "pad2m", 0, 0.05, 1, 0.3)
    pad2m.SetTopMargin(0.05)
    pad2m.SetBottomMargin(0.22)
    pad2m.SetRightMargin(0.05) # Upper and lower plot are joined
    pad2m.SetLeftMargin(0.12) # Upper and lower plot are joined
    pad2m.SetGridx() # vertical grid
    pad2m.Draw()
    pad2m.cd()       # pad2 becomes the current pad

    
    hlen2.SetLineColor(ROOT.kRed)
    hlen2.Draw()
    hlen3.SetLineColor(ROOT.kBlue)
    hlen3.Draw("same")
    hlen.SetLineColor(ROOT.kBlack)
    hlen.Draw("same")


    hqq2.SaveAs(savemean)


def Plot_prof_from_histo(h,histoname,ParamName,VarName,legenda1,legenda2,RangeUserYsigma,RangeUserYmean,savesigma,savemean,dores=True,domean=False):
    
    #h.SetName("h")
    fran = ROOT.TF1("fran","([2]*[1]) / ( TMath::Pi() * ([1]*[1] + (x-[0])*(x-[0])) )",-1000,1000) 
    fran.SetParameters(0,5,100)
    #fran.SetParLimits(2,0,200)
    h.FitSlicesY(fran,0,-1,0,"QNRL")
    #ROOT.gDirectory.GetList().Print()
    h_mean = ROOT.gDirectory.Get(histoname+"_0")
    h_sigma = ROOT.gDirectory.Get(histoname+"_1")
    h_mean.SetTitle(""+VarName+" #mu("+ParamName+") (MeV/c)")
    h_sigma.SetTitle(""+VarName+" #sigma("+ParamName+") (MeV/c)")

    proj = h.ProjectionX()
    proj.SetTitle(""+VarName+" n_{ev}")

    h_sigma.SetMarkerStyle(20)
    h_mean.SetMarkerStyle(20)

    if(dores):
        ROOT.gStyle.SetOptStat(0)
        hqq = ROOT.TCanvas("hqq","hqq",800,800)
        hqq.Draw()
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.35, 1, 1.0)
        pad1.SetBottomMargin(0.015) # Upper and lower plot are joined
        pad1.SetRightMargin(0.05) # Upper and lower plot are joined
        pad1.SetLeftMargin(0.12) # Upper and lower plot are joined
        pad1.SetGridx()         # Vertical grid
        pad1.Draw()             # Draw the upper pad: pad1
        pad1.cd()               # pad1 becomes the current pad
        h_sigma.GetYaxis().SetLabelSize(0.03)
        h_sigma.GetXaxis().SetLabelSize(0.0)
        h_sigma.GetYaxis().SetRangeUser(RangeUserYsigma[0],RangeUserYsigma[1])
        h_sigma.Draw()

        legend = ROOT.TLegend(0.6,0.72,0.92,0.88)
        legend.AddEntry(h_sigma,legenda1,"pl")
        legend.Draw("same")

        hqq.cd()          # Go back to the main canvas before defining pad2
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
        pad2.SetTopMargin(0.05)
        pad2.SetBottomMargin(0.22)
        pad2.SetRightMargin(0.05) # Upper and lower plot are joined
        pad2.SetLeftMargin(0.12) # Upper and lower plot are joined
        pad2.SetGridx() # vertical grid
        pad2.Draw()
        pad2.cd()       # pad2 becomes the current pad

        proj.GetYaxis().SetTitleSize(20)
        proj.GetYaxis().SetTitleFont(43)
        proj.GetYaxis().SetTitleOffset(1.55)
        proj.GetYaxis().SetLabelFont(43) # Absolute font size in pixel (precision 3)
        proj.GetYaxis().SetLabelSize(15)
        proj.GetYaxis().SetNdivisions(505)

        proj.GetXaxis().SetTitleSize(20)
        proj.GetXaxis().SetTitleFont(43)
        proj.GetXaxis().SetTitleOffset(1)
        proj.GetXaxis().SetLabelFont(43) # Absolute font size in pixel (precision 3)
        proj.GetXaxis().SetLabelSize(15)
        proj.Draw()

        hqq.Draw()
        hqq.SaveAs(savesigma)

    if(domean):
        ROOT.gStyle.SetOptStat(0)
        hqq2 = ROOT.TCanvas("hqq2","hqq2",800,800)
        pad1m = ROOT.TPad("pad1m", "pad1m", 0, 0.35, 1, 1.0)
        pad1m.SetBottomMargin(0.015) # Upper and lower plot are joined
        pad1m.SetRightMargin(0.05) # Upper and lower plot are joined
        pad1m.SetLeftMargin(0.12) # Upper and lower plot are joined
        pad1m.SetGridx()         # Vertical grid
        pad1m.Draw()             # Draw the upper pad: pad1
        pad1m.cd()               # pad1 becomes the current pad
        h_mean.GetYaxis().SetRangeUser(RangeUserYmean[0],RangeUserYmean[1])
        h_mean.GetYaxis().SetLabelSize(0.03)
        h_mean.GetXaxis().SetLabelSize(0.0)
        h_mean.Draw()

        h_mean.Draw("same")
        legendm = ROOT.TLegend(0.6,0.72,0.92,0.88)
        #legend.SetBorderSize(0)
        legendm.AddEntry(h_sigma,legenda2,"pl")
        legendm.Draw()

        hqq2.cd()          # Go back to the main canvas before defining pad2
        pad2m = ROOT.TPad("pad2m", "pad2m", 0, 0.05, 1, 0.3)
        pad2m.SetTopMargin(0.05)
        pad2m.SetBottomMargin(0.22)
        pad2m.SetRightMargin(0.05) # Upper and lower plot are joined
        pad2m.SetLeftMargin(0.12) # Upper and lower plot are joined
        pad2m.SetGridx() # vertical grid
        pad2m.Draw()
        pad2m.cd()       # pad2 becomes the current pad

        proj.Draw()
        hqq2.Draw()

        hqq2.SaveAs(savemean)


        


def Plot_prof_Unit(tree,ParamType,ParamMCType,Var,VarName,BinsRangeX,BinsRangeY,BinsRangeVar,RangeUserYsigma,savesigma,extracond):
    tree.Draw("("+ParamType+".fP[4]-"+ParamMCType+".fP[4])/sqrt("+ParamType+".fC[14]):"+Var+">>huseed4la("+BinsRangeX+","+BinsRangeY+")",extracond,"colz")
    huseed4la = ROOT.gPad.GetPrimitive("huseed4la")
    huseed4la.FitSlicesY()
    huseed4la_sigma = ROOT.gDirectory.Get("huseed4la_2")
    huseed4la_mean = ROOT.gDirectory.Get("huseed4la_1")
    huseed4la_sigma.SetTitle(""+VarName+" #sigma(UnitSeed)")
    huseed4la_mean.SetTitle(""+VarName+" #mu(UnitSeed)")

    tree.Draw("("+ParamType+".fP[3]-"+ParamMCType+".fP[3])/sqrt("+ParamType+".fC[9]):"+Var+">>huseed3la("+BinsRangeX+","+BinsRangeY+")",extracond,"colz")
    huseed3la = ROOT.gPad.GetPrimitive("huseed3la")
    huseed3la.FitSlicesY()
    huseed3la_sigma = ROOT.gDirectory.Get("huseed3la_2")
    huseed3la_mean = ROOT.gDirectory.Get("huseed3la_1")
    huseed3la_sigma.SetTitle(""+VarName+" #sigma(UnitSeed)")
    huseed3la_mean.SetTitle(""+VarName+" #mu(UnitSeed)")

    tree.Draw("("+ParamType+".fP[2]-"+ParamMCType+".fP[2])/sqrt("+ParamType+".fC[5]):"+Var+">>huseed2la("+BinsRangeX+","+BinsRangeY+")",extracond,"colz")
    huseed2la = ROOT.gPad.GetPrimitive("huseed2la")
    huseed2la.FitSlicesY()
    huseed2la_sigma = ROOT.gDirectory.Get("huseed2la_2")
    huseed2la_mean = ROOT.gDirectory.Get("huseed2la_1")
    huseed2la_sigma.SetTitle(""+VarName+" #sigma(UnitSeed)")
    huseed2la_mean.SetTitle(""+VarName+" #mu(UnitSeed)")

    tree.Draw(Var+">>hula("+BinsRangeVar+")",extracond,"")
    hla = ROOT.gPad.GetPrimitive("hula")
    hla.SetTitle(""+VarName+" n_{ev}")


    huseed2la_sigma.SetMarkerStyle(20)

    huseed3la_sigma.SetMarkerStyle(21)
    huseed3la_sigma.SetMarkerColor(ROOT.kRed)

    huseed4la_sigma.SetMarkerStyle(22)
    huseed4la_sigma.SetMarkerColor(ROOT.kBlue)

    ROOT.gStyle.SetOptStat(0)
    hqq = ROOT.TCanvas("hqq","hqq",800,800)
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.35, 1, 1.0)
    pad1.SetBottomMargin(0.015) # Upper and lower plot are joined
    pad1.SetRightMargin(0.05) # Upper and lower plot are joined
    pad1.SetLeftMargin(0.12) # Upper and lower plot are joined
    pad1.SetGridx()         # Vertical grid
    pad1.Draw()             # Draw the upper pad: pad1
    pad1.cd()               # pad1 becomes the current pad
    huseed4la_sigma.Draw()
    huseed4la_sigma.GetYaxis().SetLabelSize(0.03)
    huseed4la_sigma.GetXaxis().SetLabelSize(0.0)
    huseed4la_sigma.GetYaxis().SetRangeUser(RangeUserYsigma[0],RangeUserYsigma[1])
    axis = ROOT.TGaxis( -5, 20, -5, 220, 20,220,510,"")
    axis.SetLabelFont(43) # Absolute font size in pixel (precision 3)
    axis.SetLabelSize(15)
    axis.Draw()

    huseed3la_sigma.Draw("same")
    huseed2la_sigma.Draw("same")
    legend = ROOT.TLegend(0.6,0.75,0.88,0.88)
    #legend.SetBorderSize(0)
    #legend.AddEntry(hplaGAr_sigma,"GArSoft Resolution","pl")
    legend.AddEntry(huseed4la_sigma,"Seed Pull Test p4","pl")
    legend.AddEntry(huseed3la_sigma,"Seed Pull Test p3","pl")
    legend.AddEntry(huseed2la_sigma,"Seed Pull Test p2","pl")
    legend.Draw()
    hqq.Draw()

    hqq.cd()          # Go bacla to the main canvas before defining pad2
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0.05)
    pad2.SetBottomMargin(0.22)
    pad2.SetRightMargin(0.05) # Upper and lower plot are joined
    pad2.SetLeftMargin(0.12) # Upper and lower plot are joined
    pad2.SetGridx() # vertical grid
    pad2.Draw()
    pad2.cd()       # pad2 becomes the current pad

    hla.GetYaxis().SetTitleSize(20)
    hla.GetYaxis().SetTitleFont(43)
    hla.GetYaxis().SetTitleOffset(1.55)
    hla.GetYaxis().SetLabelFont(43) # Absolute font size in pixel (precision 3)
    hla.GetYaxis().SetLabelSize(15)
    hla.GetYaxis().SetNdivisions(505)

    hla.GetXaxis().SetTitleSize(20)
    hla.GetXaxis().SetTitleFont(43)
    hla.GetXaxis().SetTitleOffset(1)
    hla.GetXaxis().SetLabelFont(43) # Absolute font size in pixel (precision 3)
    hla.GetXaxis().SetLabelSize(15)
    hla.Draw()
    hqq.SaveAs(savesigma)





def Plot_residuals(tree, ParamType,ParamMCType,save0,save1,save2,save3,save4,savep2gaus,savepgaus,ranges,extracon,plimits):
    hq0 = ROOT.TCanvas("hq0","hq0",800,600)
    tree.Draw("("+ParamMCType+".fP[0]-"+ParamType+".fP[0])/"+ParamMCType+".fP[0]>>htemp0(100,-0.1,0.1)",ParamType+".fP[4]!=0&&"+ParamMCType+".fP[4]!=0"+extracon)
    h0 = ROOT.gPad.GetPrimitive("htemp0")
    h0.GetYaxis().SetRangeUser(ranges[0],ranges[1])
    h0.SetTitle("param 0 residual(p0_{reco}-p0_{MC})/p0_{MC}n")
    hq0.cd()
    h0.Draw()
    hq0.Draw()
    hq0.SaveAs(save0)

    hq1 = ROOT.TCanvas("hq1","hq1",800,600)
    tree.Draw("("+ParamMCType+".fP[1]-"+ParamType+".fP[1])/"+ParamMCType+".fP[1]>>htemp1(100,-0.2,0.2)",ParamType+".fP[4]!=0&&"+ParamMCType+".fP[4]!=0"+extracon)
    h1 = ROOT.gPad.GetPrimitive("htemp1")
    h1.GetYaxis().SetRangeUser(ranges[2],ranges[3])
    h1.SetTitle("param 1 residual(p1_{reco}-p1_{MC})/p1_{MC}n")
    hq1.cd()
    hq1.Draw()
    hq1.SaveAs(save1)

    hq2 = ROOT.TCanvas("hq2","hq2",800,600)
    tree.Draw("("+ParamMCType+".fP[2]-"+ParamType+".fP[2])/"+ParamMCType+".fP[2]>>htemp2(100,-2,2)",ParamType+".fP[4]!=0&&"+ParamMCType+".fP[4]!=0"+extracon)
    h2 = ROOT.gPad.GetPrimitive("htemp2")
    h2.SetTitle("param 2 residual(p2_{reco}-p2_{MC})/p2_{MC}n")
    h2.GetYaxis().SetRangeUser(ranges[4],ranges[5])
    hq2.cd()
    hq2.Draw()
    hq2.SaveAs(save2)

    hq3 = ROOT.TCanvas("hq3","hq3",800,600)
    tree.Draw("("+ParamMCType+".fP[3]-"+ParamType+".fP[3])/"+ParamMCType+".fP[3]>>htemp3(100,-0.2,0.2)",ParamType+".fP[4]!=0&&"+ParamMCType+".fP[4]!=0"+extracon)
    h3 = ROOT.gPad.GetPrimitive("htemp3")
    h3.SetTitle("param 3 residual(p3_{reco}-p3_{MC})/p3_{MC}n")
    h3.GetYaxis().SetRangeUser(ranges[6],ranges[7])
    hq3.cd()
    hq3.Draw()
    hq3.SaveAs(save3)

    hq4 = ROOT.TCanvas("hq4","hq4",800,600)
    tree.Draw("("+ParamMCType+".fP[4]-"+ParamType+".fP[4])/"+ParamMCType+".fP[4]>>htemp4(100,-0.3,0.3)",ParamType+".fP[4]!=0&&"+ParamMCType+".fP[4]!=0"+extracon)
    h4 = ROOT.gPad.GetPrimitive("htemp4")
    h4.SetTitle("param 4 residual(p4_{reco}-p4_{MC})/p4_{MC}n")
    h4.GetYaxis().SetRangeUser(ranges[8],ranges[9])
    hq4.cd()
    hq4.Draw()
    hq4.SaveAs(save4)

    ROOT.gStyle.SetOptFit(kTRUE)
    ROOT.gStyle.SetOptStat(1010)
    hqp = ROOT.TCanvas("hqp","hqp",800,600)
    tree.Draw("("+ParamMCType+".GetP()-"+ParamType+".GetP())/"+ParamMCType+".GetP()>>htempp("+plimits+")",ParamType+".fP[4]!=0&&"+ParamMCType+".fP[4]!=0"+extracon,"E")
    hp = ROOT.gPad.GetPrimitive("htempp")
    hp.SetLineWidth(2)
    hp.SetMarkerSize(0.9)
    hp.SetMarkerStyle(8)
    hp.GetYaxis().SetRangeUser(ranges[10],ranges[11])
    Formula = "([0]/[2]*exp(-0.5*((x-[1])/[2])^2)+[3]/[5]*exp(-0.5*((x-[4])/[5])^2))"
    double_gauss = ROOT.TF1("double_gauss",Formula,-0.1,0.1)
    double_gauss.SetParNames("A_{core}","#mu_{core}","#sigma_{core}","A_{tail}","#mu_{tail}","#sigma_{tail}",)
    double_gauss.SetParameters(hp.GetEntries(),hp.GetMean(),hp.GetRMS(),0.5*hp.GetEntries(),hp.GetRMS(),hp.GetRMS())
    double_gauss.SetParLimits(0,0,1000)
    double_gauss.SetParLimits(1,-0.05,0.05)
    double_gauss.SetParLimits(4,-0.2,0.2)
    double_gauss.SetParLimits(3,0,1000)
    double_gauss.SetParLimits(2,0.01,0.5)
    double_gauss.SetParLimits(5,0.01,0.6)
    hp.Fit("double_gauss")
    hp.SetTitle("momentum p residual(p_{reco}-p_{MC})/p_{MC}n")
    hqp.cd()
    hqp.Draw()
    hqp.SaveAs(savep2gaus)

    ROOT.gStyle.SetOptFit(kTRUE)
    hqpg = ROOT.TCanvas("hqpg","hqpg",800,600)
    tree.Draw("("+ParamMCType+".GetP()-"+ParamType+".GetP())/"+ParamMCType+".GetP()>>htempp2(40,-0.1,0.1)",ParamType+".fP[4]!=0&&"+ParamMCType+".fP[4]!=0"+extracon,"E")
    hp2 = ROOT.gPad.GetPrimitive("htempp2")
    hp2.SetLineWidth(2)
    hp2.SetMarkerSize(0.9)
    hp2.SetMarkerStyle(8)
    hp2.GetYaxis().SetRangeUser(ranges[12],ranges[13]) #pgun
    hp2.Fit("gaus")
    hp2.SetTitle("momentum p residual(p_{reco}-p_{MC})/p_{MC}n")
    hqpg.cd()
    hqpg.Draw()
    hqpg.SaveAs(savepgaus)

def SetHisto(histo,title="",color=ROOT.kBlack,markstyle = 20,yMargin = [0,0.15],extra = True):
    fgkTextSize = 0.05
    fgkTitleSize = 0.05
    fgkMarkerSize = 1
    fgkLineWidth = 3
    fgkTextFont = 42
    fgkLabelOffset = 0.01
    fgkXTitleOffset = 1.25
    fgkYTitleOffset = 1.7
    fgkTickLength = 0.02

    if (extra): histo.SetMarkerSize(fgkMarkerSize)

    ax = histo.GetXaxis()
    ax.SetTickLength(fgkTickLength)
    ax.SetLabelFont(fgkTextFont)
    ax.SetLabelSize(fgkTextSize)
    ax.SetLabelOffset(fgkLabelOffset)
    ax.SetTitleFont(fgkTextFont)
    ax.SetTitleSize(fgkTitleSize)
    kcen = 1
    ax.CenterTitle(kcen)
    ax.SetNdivisions(505)

    ay = histo.GetYaxis()
    ay.SetTickLength(fgkTickLength)
    ay.SetLabelFont(fgkTextFont)
    ay.SetLabelSize(fgkTextSize)
    ay.SetLabelOffset(fgkLabelOffset)
    ay.SetTitleFont(fgkTextFont)
    ay.SetTitleSize(fgkTitleSize)
    kcen = 1
    ay.CenterTitle(kcen)
    ay.SetNdivisions(505)
    ay.SetRangeUser(yMargin[0],yMargin[1])

    histo.GetXaxis().SetTitleOffset(fgkXTitleOffset)
    histo.GetYaxis().SetTitleOffset(fgkYTitleOffset)

    if(extra):
        histo.SetMarkerColor(color)
        histo.SetMarkerStyle(markstyle)
        histo.SetMarkerSize(1)
        histo.SetLineColor(color)
        histo.SetLineWidth(fgkLineWidth)
    

    gStyle.SetOptTitle(1)
    gStyle.SetTitleX(0.6)
    gStyle.SetTitleW(0.8)
    gStyle.SetOptStat(0)
    histo.SetTitle(title)

def SetEff(histo,graph,title="",color=ROOT.kBlack,markstyle = 20,yMargin = [0,0.15],extra = True):
    fgkTextSize = 0.05
    fgkTitleSize = 0.05
    fgkMarkerSize = 1
    fgkLineWidth = 3
    fgkTextFont = 42
    fgkLabelOffset = 0.01
    fgkXTitleOffset = 1.25
    fgkYTitleOffset = 1.7
    fgkTickLength = 0.02

    if (extra): histo.SetMarkerSize(fgkMarkerSize)

    ax = graph.GetXaxis()
    ax.SetTickLength(fgkTickLength)
    ax.SetLabelFont(fgkTextFont)
    ax.SetLabelSize(fgkTextSize)
    ax.SetLabelOffset(fgkLabelOffset)
    ax.SetTitleFont(fgkTextFont)
    ax.SetTitleSize(fgkTitleSize)
    kcen = 1
    ax.CenterTitle(kcen)
    ax.SetNdivisions(505)

    ay = graph.GetYaxis()
    ay.SetTickLength(fgkTickLength)
    ay.SetLabelFont(fgkTextFont)
    ay.SetLabelSize(fgkTextSize)
    ay.SetLabelOffset(fgkLabelOffset)
    ay.SetTitleFont(fgkTextFont)
    ay.SetTitleSize(fgkTitleSize)
    kcen = 1
    ay.CenterTitle(kcen)
    ay.SetNdivisions(505)
    ay.SetRangeUser(yMargin[0],yMargin[1])

    graph.GetXaxis().SetTitleOffset(fgkXTitleOffset)
    graph.GetYaxis().SetTitleOffset(fgkYTitleOffset)

    if(extra):
        histo.SetMarkerColor(color)
        histo.SetMarkerStyle(markstyle)
        histo.SetMarkerSize(1)
        histo.SetLineColor(color)
        histo.SetLineWidth(fgkLineWidth)
    

    gStyle.SetOptTitle(1)
    gStyle.SetTitleX(0.6)
    gStyle.SetTitleW(0.8)
    gStyle.SetOptStat(0)
    histo.SetTitle(title)

def SetLegend(lg):
    lg.SetMargin(0.18)
    lg.SetFillStyle(-1)
    lg.SetBorderSize(-1)
    lg.SetTextFont(42)
    lg.SetTextSize(0.05*0.68)
    lg.SetTextAlign(12)

def SetCanvas(cc):
    cc.SetTicks(1,1)
    cc.SetFillColor(0)
    cc.SetLeftMargin(0.17)
    cc.SetRightMargin(0.035)
    cc.SetTopMargin(0.06)
    cc.SetBottomMargin(0.14)

def Draw3HistosRes(tree,Y,X,xrange,yrange,Xname,extracond,cc,lg,yrangeuser,whichfit=2):

    paramname = ""
    if(whichfit==2): paramname = "(a)"
    else: paramname = "(b)"

    paramsymbol = ""
    if(whichfit==2): paramsymbol= "#it{#sigma} ( #it{p}_{reco}/#it{p}_{true} -1 )"
    else: paramsymbol= "#it{#mu} ( #it{p}_{reco}/#it{p}_{true} -1 )"

    paramnumber = ""
    if(whichfit==2): paramnumber = "2"
    else: paramnumber = "1"


    tree.Draw(Y+":"+X+">>hpkGAr("+xrange+","+yrange+")",extracond[0],"colz")
    hpkGAr = ROOT.gPad.GetPrimitive("hpkGAr")
    hpkGAr.FitSlicesY()
    hpkGAr_sigma = ROOT.gDirectory.Get("hpkGAr_"+paramnumber)

    tree.Draw(Y+":"+X+">>hpkGAr13("+xrange+","+yrange+")",extracond[1],"colz")
    hpkGAr13 = ROOT.gPad.GetPrimitive("hpkGAr13")
    hpkGAr13.FitSlicesY()
    hpkGAr13_sigma = ROOT.gDirectory.Get("hpkGAr13_"+paramnumber)

    tree.Draw(Y+":"+X+">>hpkGAr211("+xrange+","+yrange+")",extracond[2],"colz")
    hpkGAr211 = ROOT.gPad.GetPrimitive("hpkGAr211")
    hpkGAr211.FitSlicesY()
    hpkGAr211_sigma = ROOT.gDirectory.Get("hpkGAr211_"+paramnumber)

    SetHisto(hpkGAr_sigma,";"+Xname+" ;"+paramsymbol,ROOT.kBlack,20,yrangeuser)
    SetHisto(hpkGAr13_sigma,";"+Xname+" "+paramsymbol,ROOT.kBlue,26,yrangeuser)
    SetHisto(hpkGAr211_sigma," ;" +Xname+""+paramsymbol,ROOT.kRed,32,yrangeuser)

    SetLegend(lg)
    lg.SetHeader(paramname)
    lg.AddEntry(hpkGAr_sigma, " protons", "pl")
    lg.AddEntry(hpkGAr13_sigma, " muons ", "pl")
    lg.AddEntry(hpkGAr211_sigma, " pions ", "pl")

    SetCanvas(cc)
    hpkGAr_sigma.Draw("E1")
    hpkGAr13_sigma.Draw("E1 same")
    hpkGAr211_sigma.Draw("E1 same")

def Draw2HistosRes_fromhisto(h1,h2,h1name,h2name,Xname,cc,lg,lg1_txt,lg2_txt,yrangeuser,whichfit=2):

    paramname = ""
    if(whichfit==2): paramname = "(a)"
    else: paramname = "(b)"

    paramsymbol = ""
    if(whichfit==2): paramsymbol= "#it{#sigma} ( #it{p}_{reco}/#it{p}_{true} -1 )"
    else: paramsymbol= "#it{#mu} ( #it{p}_{reco}/#it{p}_{true} -1 )"

    paramnumber = ""
    if(whichfit==2): paramnumber = "2"
    else: paramnumber = "1"


    h1.FitSlicesY()
    h1_sigma = ROOT.gDirectory.Get(h1name+"_"+paramnumber)

    h2.FitSlicesY()
    h2_sigma = ROOT.gDirectory.Get(h2name+"_"+paramnumber)


    SetHisto(h1_sigma,";"+Xname+" ;"+paramsymbol,ROOT.kBlack,20,yrangeuser)
    SetHisto(h2_sigma,";"+Xname+" "+paramsymbol,ROOT.kBlue,26,yrangeuser)

    SetLegend(lg)
    lg.SetHeader(paramname)
    lg.AddEntry(h1_sigma, lg1_txt, "pl")
    lg.AddEntry(h2_sigma, lg2_txt, "pl")

    SetCanvas(cc)
    h1_sigma.Draw("E1")
    h2_sigma.Draw("E1 same")

def Draw3HistosRes_fromhisto(h1,h2,h3,h1name,h2name,h3name,Xname,cc,lg,lg1_txt,lg2_txt,lg3_txt,yrangeuser,whichfit=2):

    paramname = ""
    if(whichfit==2): paramname = "(a)"
    else: paramname = "(b)"

    paramsymbol = ""
    if(whichfit==2): paramsymbol= "#it{#sigma} ( #it{p}_{reco}/#it{p}_{true} -1 )"
    else: paramsymbol= "#it{#mu} ( #it{p}_{reco}/#it{p}_{true} -1 )"

    paramnumber = ""
    if(whichfit==2): paramnumber = "2"
    else: paramnumber = "1"


    h1.FitSlicesY()
    h1_sigma = ROOT.gDirectory.Get(h1name+"_"+paramnumber)

    h2.FitSlicesY()
    h2_sigma = ROOT.gDirectory.Get(h2name+"_"+paramnumber)

    h3.FitSlicesY()
    h3_sigma = ROOT.gDirectory.Get(h3name+"_"+paramnumber)


    SetHisto(h1_sigma,";"+Xname+" ;"+paramsymbol,ROOT.kBlack,20,yrangeuser)
    SetHisto(h2_sigma,";"+Xname+" "+paramsymbol,ROOT.kBlue,26,yrangeuser)
    SetHisto(h3_sigma,";"+Xname+" "+paramsymbol,ROOT.kRed,32,yrangeuser)

    h1_sigma.GetYaxis().SetMaxDigits()
    h2_sigma.GetYaxis().SetLabelFont(42)
    h3_sigma.GetYaxis().SetLabelFont(42)

    SetLegend(lg)
    lg.SetHeader(paramname)
    lg.AddEntry(h1_sigma, lg1_txt, "pl")
    lg.AddEntry(h2_sigma, lg2_txt, "pl")
    lg.AddEntry(h3_sigma, lg3_txt, "pl")

    SetCanvas(cc)
    h1_sigma.Draw("E1")
    h2_sigma.Draw("E1 same")
    h3_sigma.Draw("E1 same")    






def SetGlobalStyle(lStat=0, kcolor=1) :
    fgkTextSize = 0.05
    fgkTitleSize = 0.05
    fgkMarkerSize = 1
    fgkTextFont = 42
    fgkLabelOffset = 0.01
    fgkXTitleOffset = 1.25
    fgkYTitleOffset = 1.6
    fgkTickLength = 0.02

    # Set gStyle
    # From plain

    #turn off figure info
    gErrorIgnoreLevel = ROOT.kWarning

    gStyle.SetFrameBorderMode(0)
    gStyle.SetFrameFillColor(0)
    gStyle.SetCanvasBorderMode(0)
    gStyle.SetPadBorderMode(0)
    gStyle.SetPadColor(10)
    gStyle.SetCanvasColor(10)
    gStyle.SetTitleFillColor(10)
    gStyle.SetTitleBorderSize(-1)
    gStyle.SetStatColor(10)
    gStyle.SetStatBorderSize(-1)
    gStyle.SetLegendBorderSize(-1)

    gStyle.SetDrawBorder(0)
    gStyle.SetTextFont(fgkTextFont)
    gStyle.SetStatFont(fgkTextFont)
    gStyle.SetStatFontSize(fgkTextSize)
    gStyle.SetStatX(0.97)
    gStyle.SetStatY(0.98)
    gStyle.SetStatH(0.03)
    gStyle.SetStatW(0.3)
    gStyle.SetTickLength(fgkTickLength,"xy")
    gStyle.SetEndErrorSize(3)
    gStyle.SetLabelSize(fgkTextSize,"xyz")
    gStyle.SetLabelFont(fgkTextFont,"xyz") 
    gStyle.SetLabelOffset(fgkLabelOffset,"xyz")
    gStyle.SetTitleFont(fgkTextFont,"xyz")  
    gStyle.SetTitleFont(fgkTextFont,"")  
    gStyle.SetTitleFontSize(fgkTitleSize)
    gStyle.SetTitleOffset(fgkXTitleOffset,"x")  
    gStyle.SetTitleOffset(fgkYTitleOffset,"y")  
    gStyle.SetTitleOffset(1.0,"z")  
    gStyle.SetTitleSize(fgkTitleSize,"xyz")  
    gStyle.SetTitleSize(fgkTitleSize,"")  
    gStyle.SetMarkerSize(fgkMarkerSize) 
    gStyle.SetPalette(1,0) 
    if (lStat):
        gStyle.SetOptTitle(1)
        gStyle.SetOptStat(1111)
        gStyle.SetOptFit(1111)

    else :
        gStyle.SetOptTitle(0)
        gStyle.SetOptStat(0)
        gStyle.SetOptFit(0)


    ROOT.TGaxis.SetMaxDigits(3)
    gStyle.SetTitleBorderSize(-1)

    if(kcolor):
        SetColor()


    gROOT.ForceStyle()


  
def SetColor():

  gStyle.SetHistFillColor(0)
  #gStyle.SetFillColor(0)//it conflicts with color palette
  gStyle.SetFrameFillColor(0)
  gStyle.SetPadColor(0)
  gStyle.SetCanvasColor(0)
  gStyle.SetStatColor(0)
  gStyle.SetTitleFillColor(0)


# kDeepSea=51,          kGreyScale=52,    kDarkBodyRadiator=53,
# kBlueYellow= 54,      kRainBow=55,      kInvertedDarkBodyRadiator=56,
# kBird=57,             kCubehelix=58,    kGreenRedViolet=59,
# kBlueRedYellow=60,    kOcean=61,        kColorPrintableOnGrey=62,
# kAlpine=63,           kAquamarine=64,   kArmy=65,
# kAtlantic=66,         kAurora=67,       kAvocado=68,
# kBeach=69,            kBlackBody=70,    kBlueGreenYellow=71,
# kBrownCyan=72,        kCMYK=73,         kCandy=74,
# kCherry=75,           kCoffee=76,       kDarkRainBow=77,
# kDarkTerrain=78,      kFall=79,         kFruitPunch=80,
# kFuchsia=81,          kGreyYellow=82,   kGreenBrownTerrain=83,
# kGreenPink=84,        kIsland=85,       kLake=86,
# kLightTemperature=87, kLightTerrain=88, kMint=89,
# kNeon=90,             kPastel=91,       kPearl=92,
# kPigeon=93,           kPlum=94,         kRedBlue=95,
# kRose=96,             kRust=97,         kSandyTerrain=98,
# kSienna=99,           kSolar=100,       kSouthWest=101,
# kStarryNight=102,     kSunset=103,      kTemperatureMap=104,
# kThermometer=105,     kValentine=106,   kVisibleSpectrum=107,
# kWaterMelon=108,      kCool=109,        kCopper=110,
# kGistEarth=111,       kViridis=112,     kCividis=113

  gStyle.SetPalette(56)#only 56 available

  return
  



  