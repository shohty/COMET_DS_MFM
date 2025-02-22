#!/usr/bin/env python3

import ROOT
import numpy as np
import os

Xaxlist, p1list, p1errlist = [], [], []
def plot(OPERAfile,MEASfile):
    file0= ROOT.TFile.Open(OPERAfile, "READ")#OPERAシミュレーションファイルの読み込み
    file1 = ROOT.TFile.Open(MEASfile, "READ")#OPERAシミュレーションファイルの読み込み
    axname =  os.path.splitext(os.path.basename(OPERAfile))[0]
    if not file0 or file0.IsZombie():
        print("Error : could not open the ROOT File")
    if not file1 or file1.IsZombie():
        print("Error : could not open the ROOT File")
    #treeを取得
    opetree = file0.Get("tree")
    meastree = file1.Get("tree")
    #OPERAファイルより座標と磁場の情報をリストへ
    xlist, ylist, zlist=[], [], []
    Bxlist, Bylist, Bzlist, Blist = [], [], [], []
    for i in range(opetree.GetEntries()):
        opetree.GetEntry(i)
        xlist.append(getattr(opetree, "X"))
        ylist.append(getattr(opetree, "Y"))
        zlist.append(getattr(opetree, "Z"))
        Bxlist.append(getattr(opetree, "Bx"))
        Bylist.append(getattr(opetree, "By"))
        Bzlist.append(getattr(opetree, "Bz"))
        Blist.append(getattr(opetree, "B"))
    #OPERAファイルより座標と磁場の情報をリストへ
    #この値は3/5倍する必要あり
    Bxstdlist, Bystdlist, Bzstdlist = [], [], []
    Bxmeaslist, Bymeaslist, Bzmeaslist = [], [], []
    for i  in range(meastree.GetEntries()):
        meastree.GetEntry(i)
        Bxmeaslist.append(getattr(meastree, "Xmean") * (-3/5))#Bxは正方向が逆
        Bymeaslist.append(getattr(meastree, "Ymean") * (3/5))
        Bzmeaslist.append(getattr(meastree, "Zmean") * (3/5))
        Bxstdlist.append(getattr(meastree, "Xstdev") * (3/5))
        Bystdlist.append(getattr(meastree, "Ystdev") * (3/5))
        Bzstdlist.append(getattr(meastree, "Zstdev") * (3/5))
    
    #座標、磁場のarray：listより高速
    X = np.array(xlist, dtype=np.float64)
    Y = np.array(ylist, dtype=np.float64)
    Z = np.array(zlist, dtype=np.float64)
    Bx = np.array(Bxlist, dtype=np.float64)
    By = np.array(Bylist, dtype=np.float64)
    Bz = np.array(Bzlist, dtype=np.float64)
    B = np.array(Blist, dtype=np.float64)
    Bxmeas = np.array(Bxmeaslist, dtype=np.float64)
    Bymeas = np.array(Bymeaslist, dtype=np.float64)
    Bzmeas = np.array(Bzmeaslist, dtype=np.float64)
    Bmeas = np.sqrt(Bxmeas ** 2 + Bymeas ** 2 + Bzmeas ** 2)
    Bxstd = np.array(Bxstdlist, dtype=np.float64)
    Bystd = np.array(Bystdlist, dtype=np.float64)
    Bzstd = np.array(Bzstdlist, dtype=np.float64)
    Bstd = np.sqrt(((Bx/B) ** 2) * (Bxstd ** 2) + ((By/B) ** 2) * (Bystd ** 2) + ((Bz/B) ** 2) * (Bxstd ** 2))
    
    theta_x = np.deg2rad(0.4)#Bxセンサーのxz平面での角度
    phi_x = np.deg2rad(90-0.5)#Bxセンサーのxy平面での角度
    theta_y = np.deg2rad(-0.2)#Byセンサーのyz平面での角度
    phi_y = np.deg2rad(0.2)#Byセンサーのxy平面での角度
    theta_z = np.deg2rad(90+1.8)#Bzセンサーのxz平面での角度
    phi_z = np.deg2rad(90-0.4)#Bzセンサーのyz平面での角度
    
    #センサーの角度補正をかける。角度の定義や計算式などログノート要チェック
    Bxcorr = Bx * np.cos(theta_x) * np.sin(phi_x) - By * np.cos(phi_x) - Bz * np.sin(theta_x)
    Bycorr = By * np.cos(theta_y) * np.cos(phi_y) - Bx * np.sin(phi_y) - Bz * np.sin(theta_y)
    Bzcorr = Bz * np.sin(theta_z) * np.sin(phi_z) - By * np.cos(phi_z) - Bx * np.cos(theta_z)
    Bcorr =  np.sqrt(Bxcorr ** 2 + Bycorr ** 2 + Bzcorr ** 2)
    
    #キャンバス描画
    cfit = ROOT.TCanvas("cfit", f"{axname}", 800, 600)# Z vs B Fitあり angle corrected
    call = ROOT.TCanvas("call", f"{axname}", 1600, 1200)# Z vs B Fit なし　測定 & OPERA & OPERA angle corrected 各成分
    
    cfit.cd()
    g0 = ROOT.TGraph(len(Z), Z, Bcorr)
    g0.SetMarkerStyle(20)
    g0.SetMarkerColor(ROOT.kBlack)
    g0.SetTitle(f"{axname} : sensorangle-corrected")
    g0.GetXaxis().SetTitle("z (mm)")
    g0.GetYaxis().SetTitle("B (T)")
    g0.Draw("AP")
    Bmax = np.max(Bcorr)
    Bmax_index = np.argmax(Bcorr)
    BmaxZ = Z[Bmax_index]
    quadfit = ROOT.TF1("quadfit", "[0]*(x - [1])^2 + [2]", BmaxZ-150, BmaxZ+150)
    quadfit.SetParameters(-1e-6, BmaxZ, 0.85)  # 初期パラメータはp0は0にすべきではないなぜなら関数系が変わりlocal minimumに引っかかりがち
    g0.Fit(quadfit, "R")  # "R" = 指定範囲でフィット
    quadfit.SetLineColor(ROOT.kCyan)
    quadfit.Draw("same")
    cfit.SetGridx()
    cfit.SetGridy()
    cfit.Draw()
    ROOT.gStyle.SetOptFit(1)
    cfit.Update()
    cfit.Show()  # オプション: Canvas を表示
    cfit_dir = "/Users/shohtatakami/physics/COMETDS/DS189A_944turns_Integration/Anglecorr/corr_fit/"
    outputFile = ROOT.TFile(f"{cfit_dir}angelcorr{axname}_out.root", "RECREATE")
    cfit.Write()
    outputFile.Close()
    print(f"Successfuly Saved : {outputFile}")
    
    #peakposition分布のプロット作成に必要
    Xaxlist.append(np.mean(X))
    p1list.append(quadfit.GetParameter(1))
    p1errlist.append(quadfit.GetParError(1))
    
    call.Divide(2,2)
    gOPERAx = ROOT.TGraph(len(Z), Z, Bx)
    gOPERAx.SetMarkerStyle(8)
    gOPERAx.SetMarkerColor(ROOT.kRed+2)
    gOPERAy = ROOT.TGraph(len(Z), Z, By)
    gOPERAy.SetMarkerStyle(8)
    gOPERAy.SetMarkerColor(ROOT.kRed+2)
    gOPERAz = ROOT.TGraph(len(Z), Z, Bz)
    gOPERAz.SetMarkerStyle(8)
    gOPERAz.SetMarkerColor(ROOT.kRed+2)
    gOPERAB = ROOT.TGraph(len(Z), Z, B)
    gOPERAB.SetMarkerStyle(8)
    gOPERAB.SetMarkerColor(ROOT.kRed+2)
    gmeasx = ROOT.TGraphErrors(len(Z), Z, Bxmeas, np.zeros(len(Z)), Bxstd)
    gmeasx.SetMarkerStyle(7)
    gmeasx.SetMarkerColor(ROOT.kBlue)
    gmeasy = ROOT.TGraphErrors(len(Z), Z, Bymeas, np.zeros(len(Z)), Bystd)
    gmeasy.SetMarkerStyle(8)
    gmeasy.SetMarkerColor(ROOT.kBlue)
    gmeasz = ROOT.TGraphErrors(len(Z), Z, Bzmeas, np.zeros(len(Z)), Bzstd)
    gmeasz.SetMarkerStyle(8)
    gmeasz.SetMarkerColor(ROOT.kBlue)
    gmeasB = ROOT.TGraphErrors(len(Z), Z, Bmeas, np.zeros(len(Z)), Bstd)
    gmeasB.SetMarkerStyle(8)
    gmeasB.SetMarkerColor(ROOT.kBlue)
    gcorrx = ROOT.TGraph(len(Z), Z, Bxcorr)
    gcorrx.SetMarkerStyle(8)
    gcorrx.SetMarkerColor(ROOT.kSpring)
    gcorry = ROOT.TGraph(len(Z), Z, Bycorr)
    gcorry.SetMarkerStyle(8)
    gcorry.SetMarkerColor(ROOT.kSpring)
    gcorrz = ROOT.TGraph(len(Z), Z, Bzcorr)
    gcorrz.SetMarkerStyle(8)
    gcorrz.SetMarkerColor(ROOT.kSpring)
    gcorrB = ROOT.TGraph(len(Z), Z, Bcorr)
    gcorrB.SetMarkerStyle(8)
    gcorrB.SetMarkerColor(ROOT.kSpring)
    #それぞれを分割したキャンバスにMultiGraphにして描画
    call.cd(1)
    mgBx = ROOT.TMultiGraph()
    mgBx.Add(gOPERAx)
    mgBx.Add(gmeasx)
    mgBx.Add(gcorrx)
    mgBx.SetTitle(f"{axname} : Bx")
    mgBx.GetYaxis().SetTitle("Bx (T)")
    mgBx.GetXaxis().SetTitle("Z (mm)")
    mgBx.Draw("AP")
    ROOT.gPad.SetGrid(1,1)
    legendBx = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
    legendBx.AddEntry(gOPERAx, "OPERA", "p")
    legendBx.AddEntry(gmeasx, "Measured", "p")
    legendBx.AddEntry(gcorrx, "OPERA corr", "p")
    legendBx.SetFillStyle(0)
    legendBx.Draw()
    call.cd(2)
    mgBy = ROOT.TMultiGraph()
    mgBy.Add(gOPERAy)
    mgBy.Add(gmeasy)
    mgBy.Add(gcorry)
    mgBy.SetTitle(f"{axname} : By")
    mgBy.GetYaxis().SetTitle("By (T)")
    mgBy.GetXaxis().SetTitle("Z (mm)")
    mgBy.Draw("AP")
    ROOT.gPad.SetGrid(1,1)
    legendBy = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
    legendBy.AddEntry(gOPERAy, "OPERA", "p")
    legendBy.AddEntry(gmeasy, "Measured", "p")
    legendBy.AddEntry(gcorry, "OPERA corr", "p")
    legendBy.SetFillStyle(0)
    legendBy.Draw()
    call.cd(3)
    mgBz = ROOT.TMultiGraph()
    mgBz.Add(gOPERAz)
    mgBz.Add(gmeasz)
    mgBz.Add(gcorrz)
    mgBz.SetTitle(f"{axname} : Bz")
    mgBz.GetYaxis().SetTitle("Bz (T)")
    mgBz.GetXaxis().SetTitle("Z (mm)")
    mgBz.Draw("AP")
    ROOT.gPad.SetGrid(1,1)
    legendBz = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
    legendBz.AddEntry(gOPERAz, "OPERA", "p")
    legendBz.AddEntry(gmeasz, "Measured", "p")
    legendBz.AddEntry(gcorrz, "OPERA corr", "p")
    legendBz.SetFillStyle(0)
    legendBz.Draw()
    call.cd(4)
    mgB = ROOT.TMultiGraph()
    mgB.Add(gOPERAB)
    mgB.Add(gmeasB)
    mgB.Add(gcorrB)
    mgB.SetTitle(f"{axname} : B")
    mgB.GetYaxis().SetTitle("B (T)")
    mgB.GetXaxis().SetTitle("Z (mm)")
    mgB.Draw("AP")
    ROOT.gPad.SetGrid(1,1)
    legendB = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
    legendB.AddEntry(gOPERAB, "OPERA", "p")
    legendB.AddEntry(gmeasB, "Measured", "p")
    legendB.AddEntry(gcorrB, "OPERA corr", "p")
    legendBx.SetFillStyle(0)
    legendB.Draw()
    call.Update()
    #キャンバスをオブジェクトとして保存
    call_dir = "/Users/shohtatakami/physics/COMETDS/DS189A_944turns_Integration/Anglecorr/comp/"
    call_outFile = ROOT.TFile(f"{call_dir}{axname}_out.root", "RECREATE")
    call.Write()
    call_outFile.Close()
    print(f"Successfuly Saved : {call_outFile}")
    
if __name__=='__main__':
    operafile_directory =  "/Users/shohtatakami/physics/COMETDS/DS189A_944turns_Integration/DS189A_944turns_FieldIntegration_LaserTrackerPos/root/"
    measfile_directory = "/Users/shohtatakami/physics/COMETDS/Scan/"
    file_pairs= [
        ("X_-200Y_0_Map.root","test20240627-3.root"),
        ("X_-500Y_0_Map.root","test20240627-4.root"),
        ("X_-650Y_0_Map.root","test20240627-7.root"),
        ("X_-760Y_0_Map.root","test20240627-6.root"),
        #("X_0Y_-200_Map.root","test20240710-4.root"),
        #("X_0Y_-729_Map.root","test20240628-2.root"),
        ("X_0Y_0_Map.root","test20240626-6.root"),
        #("X_0Y_200_Map.root","test20240710-1.root"),
        #("X_0Y_780_Map.root","test20240628-0_rev.root"),
        ("X_200Y_0_Map.root","test20240627-0.root"),
        ("X_500Y_0_Map.root","test20240627-1.root"),
        ("X_650Y_0_Map.root","test20240627-8.root"),
        ("X_760Y_0_revMap.root","test20240627-2rev.root"),#rev : the version deleted the last several strange lines 
    ]
    
    for operafilename, measfilename in file_pairs:
        OPERAfile = os.path.join(operafile_directory, operafilename)
        MEASfile = os.path.join(measfile_directory, measfilename)
        plot(OPERAfile, MEASfile)
    
    Xaxis = np.array(Xaxlist, dtype = np.float64)
    p1 = np.array(p1list, dtype = np.float64)
    p1err= np.array(p1errlist, dtype = np.float64)
    cpp = ROOT.TCanvas("cpp", "Peak Position", 800, 600)
    gpp = ROOT.TGraphErrors(len(Xaxis), Xaxis, p1, np.zeros(len(Xaxis)), p1err)
    gpp.SetMarkerColor(ROOT.kRed+2)
    gpp.SetMarkerStyle(8)
    gpp.GetXaxis().SetTitle("X (mm)")
    gpp.GetYaxis().SetTitle("Peak Position (mm)")
    gpp.Draw("AP")
    cpp.SetGrid(1,1)
    cpp.Update()
    cpp_dir = "/Users/shohtatakami/physics/COMETDS/DS189A_944turns_Integration/Anglecorr/"
    ppoutFile = ROOT.TFile(f"{cpp_dir}PeakPosition.root", "RECREATE")
    cpp.Write()
    ppoutFile.Close()
    print(f"Successfuly Saved : {ppoutFile}")
    
    #ROOT.gApplication.Run()
    
    
    
    
   
    