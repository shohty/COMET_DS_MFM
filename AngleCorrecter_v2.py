#!/usr/bin/env python3
import ROOT
import numpy as np
import os
#version2ではセンサー補正に回転行列を用いる＆測定や補正前OPERAデータでのPeakPosition分布も出力
#三次元空間でのベクトルの表記方法は(z,x,y)の順
Xaxlist, p1corr_list, p1correrr_list = [], [], [] #PeakPosition 角度補正後
p1meas_list, p1measerr_list = [], []#PeakPosition 測定
p1opera_list, p1operaerr_list = [], []#PeakPosition 角度補正前
def rotM_xz(deg):
    #xz平面における回転行列
    theta = np.deg2rad(deg)
    return np.array([[np.cos(theta), 0, -np.sin(theta)],
                     [0, 1, 0],
                     [np.sin(theta), 0, np.cos(theta)]])
def rotM_xy(deg):
    #xy平面における回転行列
    theta = np.deg2rad(deg)
    return np.array([[np.cos(theta), -np.sin(theta), 0], 
                     [np.sin(theta), np.cos(theta), 0],
                     [0, 0, 1]])
def rotM_yz(deg):
    #yz平面における回転行列
    theta = np.deg2rad(deg)
    return np.array([[1, 0, 0],
                     [0, np.cos(theta), -np.sin(theta)],
                     [0, np.sin(theta), np.cos(theta)]])
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
    Bstd = np.sqrt(((Bxmeas/Bmeas) ** 2) * (Bxstd ** 2) + ((Bymeas/Bmeas) ** 2) * (Bystd ** 2) + ((Bzmeas/Bmeas) ** 2) * (Bxstd ** 2))
    
    #各センサーの法線ベクトルを出す
    theta_Bx = -0.4#Bxセンサーのxz平面での角度
    phi_Bx = -0.5#Bxセンサーのxy平面での角度
    theta_By = 0.2#Byセンサーのyz平面での角度
    phi_By = 0.2#Byセンサーのxy平面での角度
    theta_Bz = -1.8#Bzセンサーのxz平面での角度
    phi_Bz = 0.4#Bzセンサーのyz平面での角度
    #theta_Bx = np.deg2rad(0)#Bxセンサーのxz平面での角度
    #phi_Bx = np.deg2rad(0)#Bxセンサーのxy平面での角度
    #theta_By = np.deg2rad(0)#Byセンサーのyz平面での角度
    #phi_By = np.deg2rad(0)#Byセンサーのxy平面での角度
    #theta_Bz = np.deg2rad(0)#Bzセンサーのxz平面での角度
    #phi_Bz = np.deg2rad(0)#Bzセンサーのyz平面での角度
    
    nBx =   rotM_xy(phi_Bx) @ rotM_xz(theta_Bx) @ np.array([[1],
                                                            [0],
                                                            [0]])
    print("nBx : ", nBx)
    #print("rotM_xz(30) : ", rotM_xz(theta_Bx))
    #print("rotM_xy(60) : ", rotM_xy(phi_Bx))
    #print("cos : ", np.cos(theta_Bx))
    
    nBy = rotM_xy(phi_By) @ rotM_yz(theta_By) @ np.array([[0],
                                                          [1],
                                                          [0]])
    print("nBy : ", nBy)
    
    nBz = rotM_yz(phi_Bz) @ rotM_xz(theta_Bz) @ np.array([[0],
                                                          [0],
                                                          [1]])
    print("nBz : ", nBz)
    #OPERAデータに対して角度補正をかける。各点の磁場ベクトルと各センサーの法線ベクトルの内積
    Bxcorr = np.zeros(len(Bx))
    Bycorr = np.zeros(len(Bx))
    Bzcorr = np.zeros(len(Bx))
    Bcorr = np.zeros(len(Bx))
    for i in range(len(Bx)):
        B_vec = np.array([[Bx[i]],
                     [By[i]],
                     [Bz[i]]])#磁場の三次元ベクトル
        Bxcorr[i] = np.dot(B_vec.T, nBx).item()#itemをつけないとあくまで1×1の行列として処理される
        Bycorr[i] = np.dot(B_vec.T, nBy).item()
        Bzcorr[i] = np.dot(B_vec.T, nBz).item()
    Bcorr =  np.sqrt(Bxcorr ** 2 + Bycorr ** 2 + Bzcorr ** 2)
    
    #キャンバス描画
    cBzcomp = ROOT.TCanvas("cBzcomp", f"{axname}", 1000, 1200)# Bz成分の分布比較
    cBzcomp.Divide(1,2)
    call = ROOT.TCanvas("call", f"{axname}", 1600, 1200)# Z vs B Fit なし　測定 & OPERA & OPERA angle corrected 各成分
    call.Divide(2,2)
    
    #peakposition分布のプロット作成に必要
    Xaxlist.append(np.mean(X))
    
    gOPERAx = ROOT.TGraph(len(Z), Z, Bx)
    gOPERAx.SetMarkerStyle(8)
    gOPERAx.SetMarkerColor(ROOT.kRed+2)
    gOPERAy = ROOT.TGraph(len(Z), Z, By)
    gOPERAy.SetMarkerStyle(8)
    gOPERAy.SetMarkerColor(ROOT.kRed+2)
    gOPERAz = ROOT.TGraph(len(Z), Z, Bz)
    gOPERAz.SetMarkerStyle(8)
    gOPERAz.SetMarkerColor(ROOT.kRed+2)
    gOPERAB = ROOT.TGraphErrors(len(Z), Z, B)
    #gOPERAB = ROOT.TGraphErrors(len(Z), Z, B, np.zeros(len(Z)), Bstd)
    gOPERAB.SetMarkerStyle(8)
    gOPERAB.SetMarkerColor(ROOT.kRed+2)
    '''#Fitting
    Bmaxope = np.max(B)
    BmaxopeZ = B[np.argmax(B)]
    qfopera = ROOT.TF1("qfopera", "[0]*(x - [1])^2 + [2]", BmaxopeZ-150, BmaxopeZ+150)
    qfopera.SetParameters(-1e-7, BmaxopeZ, Bmaxope)  # 初期パラメータはp0は0にすべきではないなぜなら関数系が変わりlocal minimumに引っかかりがち
    gOPERAB.Fit(qfopera, "R")  # "R" = 指定範囲でフィット
    qfopera.SetLineColor(ROOT.kCyan)
    qfopera.Draw("same")
    p1opera_list.append(qfopera.GetParameter(1))
    p1operaerr_list.append(qfopera.GetParError(1))
    '''
    #測定値グラフ定義
    #gmeasx = ROOT.TGraphErrors(len(Z), Z, Bxmeas, np.zeros(len(Z)), Bxstd)
    #gmeasx.SetMarkerStyle(7)
    #gmeasx.SetMarkerColor(ROOT.kBlue)
    #gmeasy = ROOT.TGraphErrors(len(Z), Z, Bymeas, np.zeros(len(Z)), Bystd)
    #gmeasy.SetMarkerStyle(8)
    #gmeasy.SetMarkerColor(ROOT.kBlue)
    #gmeasz = ROOT.TGraphErrors(len(Z), Z, Bzmeas, np.zeros(len(Z)), Bzstd)
    #gmeasz.SetMarkerStyle(8)
    #gmeasz.SetMarkerColor(ROOT.kBlue)
    #gmeasB = ROOT.TGraph(len(Z), Z, Bmeas)
    #gmeasB = ROOT.TGraphErrors(len(Z), Z, Bmeas, np.zeros(len(Z)), Bstd)
    #gmeasB.SetMarkerStyle(8)
    #gmeasB.SetMarkerColor(ROOT.kBlue)
    '''#Fitting
    Bmeasmax = np.max(Bmeas)
    BmeasmaxZ = Bmeas[np.argmax(Bmeas)]
    qfmeas = ROOT.TF1("qfmeas", "[0]*(x - [1])^2 + [2]", BmeasmaxZ-150, BmeasmaxZ+150)
    qfmeas.SetParameters(-1e-7, BmeasmaxZ, Bmeasmax)  # 初期パラメータはp0は0にすべきではないなぜなら関数系が変わりlocal minimumに引っかかりがち
    gmeasB.Fit(qfmeas, "R")  # "R" = 指定範囲でフィット
    qfmeas.SetLineColor(ROOT.kCyan)
    qfmeas.Draw("same")
    p1meas_list.append(qfmeas.GetParameter(1))
    p1measerr_list.append(qfmeas.GetParError(1))
    '''
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
    #gcorrB = ROOT.TGraphErrors(len(Z), Z, Bcorr, np.zeros(len(Z)), Bstd)
    gcorrB.SetMarkerStyle(8)
    gcorrB.SetMarkerColor(ROOT.kSpring)
    #それぞれを分割したキャンバスにMultiGraphにして描画
    call.cd(1)
    mgBx = ROOT.TMultiGraph()
    mgBx.Add(gOPERAx)
    #mgBx.Add(gmeasx)
    mgBx.Add(gcorrx)
    mgBx.SetTitle(f"{axname} : Bx")
    mgBx.GetYaxis().SetTitle("Bx (T)")
    mgBx.GetXaxis().SetTitle("Z (mm)")
    mgBx.Draw("AP")
    ROOT.gPad.SetGrid(1,1)
    legendBx = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
    legendBx.AddEntry(gOPERAx, "OPERA", "p")
    #legendBx.AddEntry(gmeasx, "Measured", "p")
    legendBx.AddEntry(gcorrx, "OPERA corr", "p")
    legendBx.SetFillStyle(0)
    legendBx.Draw()
    call.cd(2)
    mgBy = ROOT.TMultiGraph()
    mgBy.Add(gOPERAy)
    #mgBy.Add(gmeasy)
    mgBy.Add(gcorry)
    mgBy.SetTitle(f"{axname} : By")
    mgBy.GetYaxis().SetTitle("By (T)")
    mgBy.GetXaxis().SetTitle("Z (mm)")
    mgBy.Draw("AP")
    ROOT.gPad.SetGrid(1,1)
    legendBy = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
    legendBy.AddEntry(gOPERAy, "OPERA", "p")
    #legendBy.AddEntry(gmeasy, "Measured", "p")
    legendBy.AddEntry(gcorry, "OPERA corr", "p")
    legendBy.SetFillStyle(0)
    legendBy.Draw()
    call.cd(3)
    mgBz = ROOT.TMultiGraph()
    mgBz.Add(gOPERAz)
    #mgBz.Add(gmeasz)
    mgBz.Add(gcorrz)
    mgBz.SetTitle(f"{axname} : Bz")
    mgBz.GetYaxis().SetTitle("Bz (T)")
    mgBz.GetXaxis().SetTitle("Z (mm)")
    mgBz.Draw("AP")
    ROOT.gPad.SetGrid(1,1)
    legendBz = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
    legendBz.AddEntry(gOPERAz, "OPERA", "p")
    #legendBz.AddEntry(gmeasz, "Measured", "p")
    legendBz.AddEntry(gcorrz, "OPERA corr", "p")
    legendBz.SetFillStyle(0)
    legendBz.Draw()
    call.cd(4)
    mgB = ROOT.TMultiGraph()
    mgB.Add(gOPERAB)
    #mgB.Add(gmeasB)
    mgB.Add(gcorrB)
    mgB.SetTitle(f"{axname} : B")
    mgB.GetYaxis().SetTitle("B (T)")
    mgB.GetXaxis().SetTitle("Z (mm)")
    mgB.Draw("AP")
    ROOT.gPad.SetGrid(1,1)
    legendB = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
    legendB.AddEntry(gOPERAB, "OPERA", "p")
    #legendB.AddEntry(gmeasB, "Measured", "p")
    legendB.AddEntry(gcorrB, "OPERA corr", "p")
    legendBx.SetFillStyle(0)
    legendB.Draw()
    call.Update()
    #キャンバスをオブジェクトとして保存
    call_dir = "/Users/shohtatakami/physics/COMETDS/DS189A_944turns_Integration/Anglecorr/comp_10mm/"
    call_outFile = ROOT.TFile(f"{call_dir}{axname}.root", "RECREATE")
    call.Write()
    call_outFile.Close()
    print(f"Successfuly Saved : {call_outFile}")
    #Bz分布の比較OPERA＆OPERA corr
    #センサーは傾いているのでそれらの感じているフラックスは独立ではないなのでBzの分布がシフトしているかどうかを確認
    gcorrz = ROOT.TGraph(len(Z), Z, Bzcorr)
    gcorrz.SetMarkerStyle(8)
    gcorrz.SetMarkerColor(ROOT.kSpring)
    gOPERAz = ROOT.TGraph(len(Z), Z, Bz)
    gOPERAz.SetMarkerStyle(8)
    gOPERAz.SetMarkerColor(ROOT.kRed+2)
    #BzへのBxの漏れ込みのラフな計算
    Bxproj = Bx*np.sin(np.deg2rad(-theta_Bz))#BxのBzへの漏れ込みの量
    gBxproj = ROOT.TGraph(len(Z), Z, Bxproj) 
    gBxproj.SetMarkerStyle(8)
    gBxproj.SetMarkerColor(ROOT.kAzure+1)
    cBzcomp.cd(1)
    mgBzcomp = ROOT.TMultiGraph()
    mgBzcomp.Add(gcorrz)
    mgBzcomp.Add(gOPERAz)
    mgBzcomp.SetTitle(f"{axname} : Bz comp")
    mgBzcomp.GetXaxis().SetTitle("z (mm)")
    mgBzcomp.GetYaxis().SetTitle("Bz (T)")
    mgBzcomp.Draw("AP")
    ROOT.gPad.SetGrid(1,1)
    legendBzcomp = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
    legendBzcomp.AddEntry(gOPERAz, "OPERA", "p")
    legendBzcomp.AddEntry(gcorrz, "OPERA corr", "p")
    legendBzcomp.SetFillStyle(0)
    legendBzcomp.Draw()
    cBzcomp.cd(2)
    gBxproj.SetTitle("Projection from Bx to Bz")
    gBxproj.GetXaxis().SetTitle("z (mm)")
    gBxproj.GetYaxis().SetTitle("Bx projection (T)")
    gBxproj.Draw("AP")
    ROOT.gPad.SetGrid(1,1)
    cBzcomp.Update()
    '''
    Bmax = np.max(Bcorr)
    Bmax_index = np.argmax(Bcorr)
    BmaxZ = Z[Bmax_index]
    quadfit = ROOT.TF1("quadfit", "[0]*(x - [1])^2 + [2]", BmaxZ-150, BmaxZ+150)
    quadfit.SetParameters(-1e-7, BmaxZ, Bmax)  # 初期パラメータはp0は0にすべきではないなぜなら関数系が変わりlocal minimumに引っかかりがち
    gcorrB.Fit(quadfit, "R")  # "R" = 指定範囲でフィット
    quadfit.SetLineColor(ROOT.kCyan)
    quadfit.Draw("same")
    p1corr_list.append(quadfit.GetParameter(1))
    p1correrr_list.append(quadfit.GetParError(1))
    '''
    cBzcomp_dir = "/Users/shohtatakami/physics/COMETDS/DS189A_944turns_Integration/Anglecorr/Bzcomp_10mm/"
    outputFile = ROOT.TFile(f"{cBzcomp_dir}Bzcomp{axname}.root", "RECREATE")
    cBzcomp.Write()
    outputFile.Close()
    print(f"Successfuly Saved : {outputFile}")
    
if __name__=='__main__':
    operafile_directory =  "/Users/shohtatakami/physics/COMETDS/DS189A_944turns_Integration/DS189A_944turns_FieldIntegration/root/"
    measfile_directory = "/Users/shohtatakami/physics/COMETDS/Scan/"
    file_pairs= [
        ("DS189A_944turns_FieldIntegration_X-200_Y0.root","test20240627-3.root"),
        ("DS189A_944turns_FieldIntegration_X-500_Y0.root","test20240627-4.root"),
        ("DS189A_944turns_FieldIntegration_X-650_Y0.root","test20240627-7.root"),
        ("DS189A_944turns_FieldIntegration_X-760_Y0.root","test20240627-6.root"),
        #("X_0Y_-200_Map.root","test20240710-4.root"),
        #("X_0Y_-729_Map.root","test20240628-2.root"),
        ("DS189A_944turns_FieldIntegration_X0_Y0.root","test20240626-6.root"),
        #("X_0Y_200_Map.root","test20240710-1.root"),
        #("X_0Y_780_Map.root","test20240628-0_rev.root"),
        ("DS189A_944turns_FieldIntegration_X200_Y0.root","test20240627-0.root"),
        ("DS189A_944turns_FieldIntegration_X500_Y0.root","test20240627-1.root"),
        ("DS189A_944turns_FieldIntegration_X650_Y0.root","test20240627-8.root"),
        ("DS189A_944turns_FieldIntegration_X760_Y0.root","test20240627-2rev.root"),#rev : the version deleted the last several strange lines 
    ]
    
    for operafilename, measfilename in file_pairs:
        OPERAfile = os.path.join(operafile_directory, operafilename)
        MEASfile = os.path.join(measfile_directory, measfilename)
        plot(OPERAfile, MEASfile)
    '''
    Xaxis = np.array(Xaxlist, dtype = np.float64)
    p1meas = np.array(p1meas_list, dtype = np.float64)
    p1measerr = np.array(p1measerr_list, dtype = np.float64)
    p1opera = np.array(p1opera_list, dtype = np.float64)
    p1operaerr = np.array(p1operaerr_list, dtype = np.float64)
    p1corr = np.array(p1corr_list, dtype = np.float64)
    p1correrr = np.array(p1correrr_list, dtype = np.float64)
    cpp = ROOT.TCanvas("cpp", "Peak Position", 800, 600)
    gppmeas = ROOT.TGraphErrors(len(Xaxis), Xaxis, p1meas, np.zeros(len(Xaxis)), p1measerr)
    gppmeas.SetMarkerColor(ROOT.kBlue)
    gppmeas.SetMarkerStyle(8)
    gppopera = ROOT.TGraphErrors(len(Xaxis), Xaxis, p1opera, np.zeros(len(Xaxis)), p1operaerr)
    gppopera.SetMarkerColor(ROOT.kRed+2)
    gppopera.SetMarkerStyle(8)
    gppcorr = ROOT.TGraphErrors(len(Xaxis), Xaxis, p1corr, np.zeros(len(Xaxis)), p1correrr)
    gppcorr.SetMarkerColor(ROOT.kGreen)
    gppcorr.SetMarkerStyle(8)
    mgpp = ROOT.TMultiGraph()
    mgpp.Add(gppmeas)
    mgpp.Add(gppopera)
    mgpp.Add(gppcorr)
    mgpp.GetXaxis().SetTitle("X (mm)")
    mgpp.GetYaxis().SetTitle("Peak Position (mm)")
    mgpp.Draw("AP")
    ROOT.gPad.SetGrid(1,1)
    legendpp = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
    legendpp.AddEntry(gppmeas, "Measured", "p")
    legendpp.AddEntry(gppopera, "OPERA", "p")
    legendpp.AddEntry(gppcorr, "Angle Corr", "p")
    legendpp.SetFillStyle(0)
    legendpp.Draw()
    cpp.Update()

    cpp_dir = "/Users/shohtatakami/physics/COMETDS/DS189A_944turns_Integration/Anglecorr/"
    ppoutFile = ROOT.TFile(f"{cpp_dir}PeakPosition.root", "RECREATE")
    cpp.Write()
    ppoutFile.Close()
    print(f"Successfuly Saved : {ppoutFile}")
    '''
    
    #ROOT.gApplication.Run()
    
    
    
    
   
    