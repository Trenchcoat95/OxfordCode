import numpy as np
from array import array
import os
import sys
import ROOT
from ROOT import TVectorD, TMatrix, TMath, TVector3, TGraphErrors, TFile, TTree, gRandom, gPad, gROOT, gVirtualX, kTRUE, kRed, TProfile, gStyle

edge = [0.0]
nue = []
numu = []
nuebar = []
numubar = []

nbins=0

Folder = "/home/federico/Documents/Universita/Federico_2020-2021/OxfordCode/PRISM-TKI/"

with open(Folder+"DUNE_ND_neutrinomode_Flux_2017.dat") as fp:
    for line in fp:
        wholeline = []
        wholeline.extend(float(item) for item in line.split())
        edge.append(wholeline[0])
        nue.append(wholeline[1])
        numu.append(wholeline[2])
        nuebar.append(wholeline[4])
        numubar.append(wholeline[5])
        nbins+=1


hnumu = ROOT.TH1F("histonumu","histonumu",nbins,array('d',edge))
hnue = ROOT.TH1F("histonue","histonue",nbins,array('d',edge))
hnumubar = ROOT.TH1F("histonumubar","histonumubar",nbins,array('d',edge))
hnuebar = ROOT.TH1F("histonuebar","histonuebar",nbins,array('d',edge))

for i in range(nbins-1):
    hnumu.SetBinContent(i+1,numu[i])
    hnue.SetBinContent(i+1,nue[i])
    hnumubar.SetBinContent(i+1,numubar[i])
    hnuebar.SetBinContent(i+1,nuebar[i])

#for i in range(hnumu.GetNbinsX()-1): hnumu.SetBinContent(i+1,hnumu.GetBinContent(i+1)/hnumu.GetBinWidth(i+1))  
#for i in range(hnue.GetNbinsX()-1): hnue.SetBinContent(i+1,hnue.GetBinContent(i+1)/hnue.GetBinWidth(i+1))  
#for i in range(hnuebar.GetNbinsX()-1): hnuebar.SetBinContent(i+1,hnuebar.GetBinContent(i+1)/hnuebar.GetBinWidth(i+1))  
#for i in range(hnumubar.GetNbinsX()-1): hnumubar.SetBinContent(i+1,hnumubar.GetBinContent(i+1)/hnumubar.GetBinWidth(i+1)) 

hnumu.GetXaxis().SetRangeUser(0,20)
hnue.GetXaxis().SetRangeUser(0,20)
hnumubar.GetXaxis().SetRangeUser(0,20)
hnuebar.GetXaxis().SetRangeUser(0,20)

cEnumu = ROOT.TCanvas("Enumu","Enumu",800,600)
hnumu.SetTitle("#nu_{#mu} flux;E(GeV);")
hnumu.Draw()
cEnumu.Draw()
save=Folder+"Enumu.png"
cEnumu.Print(save)

cEnue = ROOT.TCanvas("Enue","Enue",800,600)
hnue.SetTitle("#nu_{e} flux;E(GeV);")
hnue.Draw()
cEnue.Draw()
save=Folder+"Enue.png"
cEnue.Print(save)

cEnumub = ROOT.TCanvas("Enumub","Enumub",800,600)
hnumubar.SetTitle("#bar{#nu}_{#mu} flux;E(GeV);")
hnumubar.Draw()
cEnumub.Draw()
save=Folder+"Enumubar.png"
cEnumub.Print(save)

cEnueb = ROOT.TCanvas("Enueb","Enueb",800,600)
hnuebar.SetTitle("#bar{#nu}_{e} flux;E(GeV);")
hnuebar.Draw()
cEnueb.Draw()
save=Folder+"Enuebar.png"
cEnueb.Print(save)


filename="/home/federico/Documents/Universita/Federico_2020-2021/OxfordCode/PRISM-TKI/histos_g4lbne_v3r5p9_QGSP_BERT_OfficialEngDesignSept2021_neutrino_GAr_center.root"
file = ROOT.TFile(filename)


histonumu = file.Get("Unosc_numu_flux_DUNEPRISM_GAr_center") 
histonue = file.Get("Unosc_nue_flux_DUNEPRISM_GAr_center") 
histonumubar = file.Get("Unosc_numubar_flux_DUNEPRISM_GAr_center") 
histonuebar = file.Get("Unosc_nuebar_flux_DUNEPRISM_GAr_center") 



for i in range(6):

    offaxis= (i)*5 #distance offaxis in meters
    nbin = offaxis*20+81  #correspondent bin to off-axis position

    Folder = "Flux_NDGAr_PRISM_"+str(offaxis)+"m/"
    Filename = Folder+"Flux_NDGAr_PRISM_"+str(offaxis)+"m.dat"
    text_file = open(Filename, "w")

        
    histonumuoff = histonumu.ProjectionX("",nbin,nbin,"").Clone() 
    histonueoff = histonue.ProjectionX("",nbin,nbin,"").Clone()
    histonumubaroff = histonumubar.ProjectionX("",nbin,nbin,"").Clone()
    histonuebaroff = histonuebar.ProjectionX("",nbin,nbin,"").Clone()

    for iev in range(histonumuoff.GetNbinsX()):
        if iev==0 : continue
        width = histonumuoff.GetBinWidth(iev)   

        totnumu=histonumuoff.GetBinContent(iev)
        totnue=histonueoff.GetBinContent(iev)
        totnumubar=histonumubaroff.GetBinContent(iev)
        totnuebar=histonuebaroff.GetBinContent(iev)
        text_file.write(str(histonumuoff.GetBinLowEdge(iev+1))+" "+str(totnue)+" "+str(totnumu)+" "+str(0.0)+" "+str(totnuebar)+" "+str(totnumubar)+" "+str(0.0)+"\n")


    for i in range(histonumuoff.GetNbinsX()-1): histonumuoff.SetBinContent(i+1,histonumuoff.GetBinContent(i+1)/histonumuoff.GetBinWidth(i+1))  
    for i in range(histonueoff.GetNbinsX()-1): histonueoff.SetBinContent(i+1,histonueoff.GetBinContent(i+1)/histonueoff.GetBinWidth(i+1))  
    for i in range(histonumubaroff.GetNbinsX()-1): histonumubaroff.SetBinContent(i+1,histonumubaroff.GetBinContent(i+1)/histonumubaroff.GetBinWidth(i+1))  
    for i in range(histonuebaroff.GetNbinsX()-1): histonuebaroff.SetBinContent(i+1,histonuebaroff.GetBinContent(i+1)/histonuebaroff.GetBinWidth(i+1))  

    histonumuoff.GetXaxis().SetRangeUser(0,20)
    histonueoff.GetXaxis().SetRangeUser(0,20)
    histonumubaroff.GetXaxis().SetRangeUser(0,20)
    histonuebaroff.GetXaxis().SetRangeUser(0,20)

    #histonumuoff.Rebin(nbins,"histonumuoff",array('d',edge))
    #histonueoff.Rebin(nbins,"histonueoff",array('d',edge))
    #histonumubaroff.Rebin(nbins,"histonumubaroff",array('d',edge))
    #histonuebaroff.Rebin(nbins,"histonuebaroff",array('d',edge))

    Canvasnumu = "Enumu"+str(offaxis)
    cEnumu = ROOT.TCanvas(Canvasnumu,Canvasnumu,800,600)
    Title = "#nu_{#mu} flux off-axis "+str(offaxis)+"m"
    histonumuoff.SetTitle(Title)
    histonumuoff.Draw()
    cEnumu.Draw()
    save=Folder+"Enumu_"+str(offaxis)+"m.png"
    cEnumu.Print(save)

    Canvasnue = "Enue"+str(offaxis)
    cEnue = ROOT.TCanvas(Canvasnue,Canvasnue,800,600)
    Title = "#nu_{e} flux off-axis "+str(offaxis)+"m"
    histonueoff.SetTitle(Title)
    histonueoff.Draw()
    cEnue.Draw()
    save=Folder+"Enue_"+str(offaxis)+"m.png"
    cEnue.Print(save)

    Canvasnumubar = "Enumubar"+str(offaxis)
    cEnumubar = ROOT.TCanvas(Canvasnumubar,Canvasnumubar,800,600)
    Title = "#bar{#nu}_{#mu} flux off-axis "+str(offaxis)+"m"
    histonumubaroff.SetTitle(Title)
    histonumubaroff.Draw()
    cEnumubar.Draw()
    save=Folder+"Enumubar_"+str(offaxis)+"m.png"
    cEnumubar.Print(save)

    Canvasnuebar = "Enuebar"+str(offaxis)
    cEnuebar = ROOT.TCanvas(Canvasnuebar,Canvasnuebar,800,600)
    Title = "#bar{#nu}_{e} flux off-axis "+str(offaxis)+"m"
    histonuebaroff.SetTitle(Title)
    histonuebaroff.Draw()
    cEnuebar.Draw()
    save=Folder+"Enuebar_"+str(offaxis)+"m.png"
    cEnuebar.Print(save)



    

