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
def plot_OPERAcorrected(OPERAfile,MEASfile):
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
        #OPERAファイルより座標と磁場の情報をリストへ
        opetree.GetEntry(i)
        xlist.append(getattr(opetree, "X"))
        ylist.append(getattr(opetree, "Y"))
        zlist.append(getattr(opetree, "Z"))
        Bxlist.append(getattr(opetree, "Bx"))
        Bylist.append(getattr(opetree, "By"))
        Bzlist.append(getattr(opetree, "Bz"))
        Blist.append(getattr(opetree, "B"))
    
    #この値は3/5倍する必要あり
    Bxstdlist, Bystdlist, Bzstdlist = [], [], []
    Bxmeaslist, Bymeaslist, Bzmeaslist = [], [], []
    Xlist, Ylist,Zlist = [], [], []
    for i  in range(meastree.GetEntries()):
        meastree.GetEntry(i)
        Bxmeaslist.append(getattr(meastree, "Xmean") * (-3/5))#Bxは正方向が逆
        Bymeaslist.append(getattr(meastree, "Ymean") * (3/5))
        Bzmeaslist.append(getattr(meastree, "Zmean") * (3/5))
        Bxstdlist.append(getattr(meastree, "Xstdev") * (3/5))
        Bystdlist.append(getattr(meastree, "Ystdev") * (3/5))
        Bzstdlist.append(getattr(meastree, "Zstdev") * (3/5))
        Xlist.append(getattr(meastree, "X"))#測定点の座標
        Ylist.append(getattr(meastree, "Y"))#測定点の座標
        Zlist.append(getattr(meastree, "Z"))#測定点の座標
    
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
    Xmeas = np.array(Xlist, dtype=np.float64)#測定点の座標
    Ymeas = np.array(Ylist, dtype=np.float64)#測定点の座標
    Zmeas = np.array(Zlist, dtype=np.float64)#測定点の座標
     
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
    cBzcomp = ROOT.TCanvas("cBzcomp", f"{axname}", 1000, 600)# Bz成分の分布比較
    cBzcomp.Divide(1,2)
    
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
    #Fit角度補正なし
    BzMin = np.min(Bz)
    BzMinZ = Z[np.argmin(Bz)]
    print(f"補正前 左:{Z[np.argmin(Bz) - 1]}")
    print(f"補正前 :{BzMinZ}")
    print(f"補正前 右:{Z[np.argmin(Bz) + 1]}")
    qf = ROOT.TF1("qf", "[0]*(x - [1])^2 + [2]", BzMinZ - 50, BzMinZ + 50)#Fitting範囲は頂点＋ー150mm
    #qf = ROOT.TF1("qf", "[0]*(x - [1])^2 + [2]", Z[np.argmin(Bz) - 1], Z[np.argmin(Bz) + 1])#Fitting範囲は頂点と両サイド3点
    qf.SetParameters(1e-7, BzMinZ, BzMin)  # 初期パラメータはp0は0にすべきではないなぜなら関数系が変わりlocal minimumに引っかかりがち
    qf.SetLineColor(ROOT.kRed+2)
    gOPERAz.Fit(qf, "R")
    p1opera_list.append(qf.GetParameter(1))
    p1operaerr_list.append(qf.GetParError(1))
    qf.Draw("same")
    gOPERAB = ROOT.TGraphErrors(len(Z), Z, B)
    #gOPERAB = ROOT.TGraphErrors(len(Z), Z, B, np.zeros(len(Z)), Bstd)
    gOPERAB.SetMarkerStyle(8)
    gOPERAB.SetMarkerColor(ROOT.kRed+2)
    #測定値グラフ定義
    #gmeasx = ROOT.TGraphErrors(len(Z), Z, Bxmeas, np.zeros(len(Z)), Bxstd)
    #gmeasx.SetMarkerStyle(7)
    #gmeasx.SetMarkerColor(ROOT.kBlue)
    #gmeasy = ROOT.TGraphErrors(len(Z), Z, Bymeas, np.zeros(len(Z)), Bystd)
    #gmeasy.SetMarkerStyle(8)
    #gmeasy.SetMarkerColor(ROOT.kBlue)
    gmeasz = ROOT.TGraphErrors(len(Zmeas), Zmeas, Bzmeas, np.zeros(len(Z)), Bzstd)
    gmeasz.SetMarkerStyle(8)
    gmeasz.SetMarkerColor(ROOT.kBlue)
    BzmeasMin = np.min(Bzmeas)
    BzmeasMinZ = Zmeas[np.argmin(Bzmeas)]
    qfmeas = ROOT.TF1("qfmeas", "[0]*(x - [1])^2 + [2]", BzmeasMinZ - 50, BzmeasMinZ + 50)
    qfmeas.SetParameters(1e-7, BzmeasMinZ, BzmeasMin)
    qfmeas.SetLineColor(ROOT.kBlue)
    gmeasz.Fit(qfmeas, "R")
    p1meas_list.append(qfmeas.GetParameter(1))
    p1measerr_list.append(qfmeas.GetParError(1))
    qfmeas.Draw("same")
    #gmeasB = ROOT.TGraph(len(Z), Z, Bmeas)
    #gmeasB = ROOT.TGraphErrors(len(Z), Z, Bmeas, np.zeros(len(Z)), Bstd)
    #gmeasB.SetMarkerStyle(8)
    #gmeasB.SetMarkerColor(ROOT.kBlue)
    gcorrx = ROOT.TGraph(len(Z), Z, Bxcorr)
    gcorrx.SetMarkerStyle(8)
    gcorrx.SetMarkerColor(ROOT.kTeal+4)
    gcorry = ROOT.TGraph(len(Z), Z, Bycorr)
    gcorry.SetMarkerStyle(8)
    gcorry.SetMarkerColor(ROOT.kTeal+4)
    gcorrz = ROOT.TGraph(len(Z), Z, Bzcorr)
    gcorrz.SetMarkerStyle(8)
    gcorrz.SetMarkerColor(ROOT.kTeal+4)
    #Fit角度補正
    BzcorrMin = np.min(Bzcorr)#Bzなので分布の谷底、最小値になる
    BzcorrMinZ = Z[np.argmin(Bzcorr)]
    print(f"補正後 左:{Z[np.argmin(Bzcorr) - 1]}")
    print(f"補正後 :{BzcorrMinZ}")
    print(f"補正後 右:{Z[np.argmin(Bzcorr) + 1]}")
    #qfcorr = ROOT.TF1("qfcorr", "[0]*(x - [1])^2 + [2]", Z[np.argmin(Bzcorr) - 1], Z[np.argmin(Bzcorr) + 1])#頂点と両サイド3点
    qfcorr = ROOT.TF1("qfcorr", "[0]*(x - [1])^2 + [2]", BzcorrMinZ - 50, BzcorrMinZ + 50)#Fitting範囲は頂点と両サイド3点
    qfcorr.SetParameters(1e-7, BzcorrMinZ, BzcorrMin)#初期パラメータはp0は0にすべきではないなぜなら関数系が変わりlocal minimumに引っかかりがち(あとp0＋ー気をつける)
    qfcorr.SetLineColor(ROOT.kTeal+4)
    gcorrz.Fit(qfcorr, "R")
    p1corr_list.append(qfcorr.GetParameter(1))
    p1correrr_list.append(qfcorr.GetParError(1))
    qfcorr.Draw("same")
    gcorrB = ROOT.TGraph(len(Z), Z, Bcorr)
    #gcorrB = ROOT.TGraphErrors(len(Z), Z, Bcorr, np.zeros(len(Z)), Bstd)
    gcorrB.SetMarkerStyle(8)
    gcorrB.SetMarkerColor(ROOT.kTeal+4)
    #Bz分布OPERA＆OPERA corr (Fittingあり)
    cBzcomp.cd()
    mgBzcomp = ROOT.TMultiGraph()
    mgBzcomp.Add(gcorrz)
    mgBzcomp.Add(gOPERAz)
    mgBzcomp.Add(gmeasz)
    mgBzcomp.SetTitle(f"{axname} : Bz comp")
    mgBzcomp.GetXaxis().SetTitle("Z (mm)")
    mgBzcomp.GetYaxis().SetTitle("Bz (T)")
    mgBzcomp.Draw("AP")
    #qf.Draw("same")
    ROOT.gPad.SetGrid(1,1)
    legendBzcomp = ROOT.TLegend(0.7, 0.1, 0.9, 0.25)
    legendBzcomp.AddEntry(gOPERAz, "OPERA", "p")
    legendBzcomp.AddEntry(gcorrz, "OPERA corr", "p")
    legendBzcomp.AddEntry(gmeasz, "Measured", "p")
    legendBzcomp.SetFillStyle(0)
    legendBzcomp.Draw()
    cBzcomp.Update()
    cBzcomp_dir = "/Users/shohtatakami/physics/COMETDS/DS189A_944turns_Integration/Anglecorr/Bzcomp_10mmMap/"
    outputFile = ROOT.TFile(f"{cBzcomp_dir}{axname}Fitpm50.root", "RECREATE")
    cBzcomp.Write()
    outputFile.Close()
    print(f"Successfuly Saved : {outputFile}")
def ployt_MEAScorrected(OPERAfile,MEASfile):
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
        #OPERAファイルより座標と磁場の情報をリストへ
        opetree.GetEntry(i)
        xlist.append(getattr(opetree, "X"))
        ylist.append(getattr(opetree, "Y"))
        zlist.append(getattr(opetree, "Z"))
        Bxlist.append(getattr(opetree, "Bx"))
        Bylist.append(getattr(opetree, "By"))
        Bzlist.append(getattr(opetree, "Bz"))
        Blist.append(getattr(opetree, "B"))
    
    #この値は3/5倍する必要あり
    Bxstdlist, Bystdlist, Bzstdlist = [], [], []
    Bxmeaslist, Bymeaslist, Bzmeaslist = [], [], []
    Xlist, Ylist,Zlist = [], [], []
    for i  in range(meastree.GetEntries()):
        meastree.GetEntry(i)
        Bxmeaslist.append(getattr(meastree, "Xmean") * (-3/5))#Bxは正方向が逆
        Bymeaslist.append(getattr(meastree, "Ymean") * (3/5))
        Bzmeaslist.append(getattr(meastree, "Zmean") * (3/5))
        Bxstdlist.append(getattr(meastree, "Xstdev") * (3/5))
        Bystdlist.append(getattr(meastree, "Ystdev") * (3/5))
        Bzstdlist.append(getattr(meastree, "Zstdev") * (3/5))
        Xlist.append(getattr(meastree, "X"))#測定点の座標
        Ylist.append(getattr(meastree, "Y"))#測定点の座標
        Zlist.append(getattr(meastree, "Z"))#測定点の座標
    
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
    Xmeas = np.array(Xlist, dtype=np.float64)#測定点の座標
    Ymeas = np.array(Ylist, dtype=np.float64)#測定点の座標
    Zmeas = np.array(Zlist, dtype=np.float64)#測定点の座標
     
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
    cBzcomp = ROOT.TCanvas("cBzcomp", f"{axname}", 1000, 600)# Bz成分の分布比較
    cBzcomp.Divide(1,2)
    
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
    #Fit角度補正なし
    BzMin = np.min(Bz)
    BzMinZ = Z[np.argmin(Bz)]
    print(f"補正前 左:{Z[np.argmin(Bz) - 1]}")
    print(f"補正前 :{BzMinZ}")
    print(f"補正前 右:{Z[np.argmin(Bz) + 1]}")
    qf = ROOT.TF1("qf", "[0]*(x - [1])^2 + [2]", BzMinZ - 50, BzMinZ + 50)#Fitting範囲は頂点＋ー150mm
    #qf = ROOT.TF1("qf", "[0]*(x - [1])^2 + [2]", Z[np.argmin(Bz) - 1], Z[np.argmin(Bz) + 1])#Fitting範囲は頂点と両サイド3点
    qf.SetParameters(1e-7, BzMinZ, BzMin)  # 初期パラメータはp0は0にすべきではないなぜなら関数系が変わりlocal minimumに引っかかりがち
    qf.SetLineColor(ROOT.kRed+2)
    gOPERAz.Fit(qf, "R")
    p1opera_list.append(qf.GetParameter(1))
    p1operaerr_list.append(qf.GetParError(1))
    qf.Draw("same")
    gOPERAB = ROOT.TGraphErrors(len(Z), Z, B)
    #gOPERAB = ROOT.TGraphErrors(len(Z), Z, B, np.zeros(len(Z)), Bstd)
    gOPERAB.SetMarkerStyle(8)
    gOPERAB.SetMarkerColor(ROOT.kRed+2)
    #測定値グラフ定義
    #gmeasx = ROOT.TGraphErrors(len(Z), Z, Bxmeas, np.zeros(len(Z)), Bxstd)
    #gmeasx.SetMarkerStyle(7)
    #gmeasx.SetMarkerColor(ROOT.kBlue)
    #gmeasy = ROOT.TGraphErrors(len(Z), Z, Bymeas, np.zeros(len(Z)), Bystd)
    #gmeasy.SetMarkerStyle(8)
    #gmeasy.SetMarkerColor(ROOT.kBlue)
    gmeasz = ROOT.TGraphErrors(len(Zmeas), Zmeas, Bzmeas, np.zeros(len(Z)), Bzstd)
    gmeasz.SetMarkerStyle(8)
    gmeasz.SetMarkerColor(ROOT.kBlue)
    BzmeasMin = np.min(Bzmeas)
    BzmeasMinZ = Zmeas[np.argmin(Bzmeas)]
    qfmeas = ROOT.TF1("qfmeas", "[0]*(x - [1])^2 + [2]", BzmeasMinZ - 50, BzmeasMinZ + 50)
    qfmeas.SetParameters(1e-7, BzmeasMinZ, BzmeasMin)
    qfmeas.SetLineColor(ROOT.kBlue)
    gmeasz.Fit(qfmeas, "R")
    p1meas_list.append(qfmeas.GetParameter(1))
    p1measerr_list.append(qfmeas.GetParError(1))
    qfmeas.Draw("same")
    #gmeasB = ROOT.TGraph(len(Z), Z, Bmeas)
    #gmeasB = ROOT.TGraphErrors(len(Z), Z, Bmeas, np.zeros(len(Z)), Bstd)
    #gmeasB.SetMarkerStyle(8)
    #gmeasB.SetMarkerColor(ROOT.kBlue)
    gcorrx = ROOT.TGraph(len(Z), Z, Bxcorr)
    gcorrx.SetMarkerStyle(8)
    gcorrx.SetMarkerColor(ROOT.kTeal+4)
    gcorry = ROOT.TGraph(len(Z), Z, Bycorr)
    gcorry.SetMarkerStyle(8)
    gcorry.SetMarkerColor(ROOT.kTeal+4)
    gcorrz = ROOT.TGraph(len(Z), Z, Bzcorr)
    gcorrz.SetMarkerStyle(8)
    gcorrz.SetMarkerColor(ROOT.kTeal+4)
    #Fit角度補正
    BzcorrMin = np.min(Bzcorr)#Bzなので分布の谷底、最小値になる
    BzcorrMinZ = Z[np.argmin(Bzcorr)]
    print(f"補正後 左:{Z[np.argmin(Bzcorr) - 1]}")
    print(f"補正後 :{BzcorrMinZ}")
    print(f"補正後 右:{Z[np.argmin(Bzcorr) + 1]}")
    #qfcorr = ROOT.TF1("qfcorr", "[0]*(x - [1])^2 + [2]", Z[np.argmin(Bzcorr) - 1], Z[np.argmin(Bzcorr) + 1])#頂点と両サイド3点
    qfcorr = ROOT.TF1("qfcorr", "[0]*(x - [1])^2 + [2]", BzcorrMinZ - 50, BzcorrMinZ + 50)#Fitting範囲は頂点と両サイド3点
    qfcorr.SetParameters(1e-7, BzcorrMinZ, BzcorrMin)#初期パラメータはp0は0にすべきではないなぜなら関数系が変わりlocal minimumに引っかかりがち(あとp0＋ー気をつける)
    qfcorr.SetLineColor(ROOT.kTeal+4)
    gcorrz.Fit(qfcorr, "R")
    p1corr_list.append(qfcorr.GetParameter(1))
    p1correrr_list.append(qfcorr.GetParError(1))
    qfcorr.Draw("same")
    gcorrB = ROOT.TGraph(len(Z), Z, Bcorr)
    #gcorrB = ROOT.TGraphErrors(len(Z), Z, Bcorr, np.zeros(len(Z)), Bstd)
    gcorrB.SetMarkerStyle(8)
    gcorrB.SetMarkerColor(ROOT.kTeal+4)
    #Bz分布OPERA＆OPERA corr (Fittingあり)
    cBzcomp.cd()
    mgBzcomp = ROOT.TMultiGraph()
    mgBzcomp.Add(gcorrz)
    mgBzcomp.Add(gOPERAz)
    mgBzcomp.Add(gmeasz)
    mgBzcomp.SetTitle(f"{axname} : Bz comp")
    mgBzcomp.GetXaxis().SetTitle("Z (mm)")
    mgBzcomp.GetYaxis().SetTitle("Bz (T)")
    mgBzcomp.Draw("AP")
    #qf.Draw("same")
    ROOT.gPad.SetGrid(1,1)
    legendBzcomp = ROOT.TLegend(0.7, 0.1, 0.9, 0.25)
    legendBzcomp.AddEntry(gOPERAz, "OPERA", "p")
    legendBzcomp.AddEntry(gcorrz, "OPERA corr", "p")
    legendBzcomp.AddEntry(gmeasz, "Measured", "p")
    legendBzcomp.SetFillStyle(0)
    legendBzcomp.Draw()
    cBzcomp.Update()
    cBzcomp_dir = "/Users/shohtatakami/physics/COMETDS/DS189A_944turns_Integration/Anglecorr/Bzcomp_10mmMap/"
    outputFile = ROOT.TFile(f"{cBzcomp_dir}{axname}Fitpm50.root", "RECREATE")
    cBzcomp.Write()
    outputFile.Close()
    print(f"Successfuly Saved : {outputFile}")
if __name__=='__main__':
    operafile_directory =  "/Users/shohtatakami/physics/COMETDS/DS189A_944turns_Integration/DS189A_944turns_FieldIntegration/root/"
    measfile_directory = "/Users/shohtatakami/physics/COMETDS/newScan/"#座標追加
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
        plot_OPERAcorrected(OPERAfile, MEASfile)
    
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
    '''
    linearmeas = ROOT.TF1("linearmeas", "[0] * x^3 + [1]", np.min(Xaxis), np.max(Xaxis))
    linearmeas.SetParameters(1e-08, 0)#初期パラメータはp0は0にすべきではないなぜなら関数系が変わりlocal minimumに引っかかりがち(あとp0＋ー気をつける)
    linearmeas.SetLineColor(ROOT.kBlue)
    gppmeas.Fit(linearmeas, "R")
    linearmeas.Draw("same")
    '''
    gppopera = ROOT.TGraphErrors(len(Xaxis), Xaxis, p1opera, np.zeros(len(Xaxis)), p1operaerr)
    gppopera.SetMarkerColor(ROOT.kRed+2)
    gppopera.SetMarkerStyle(8)
    gppcorr = ROOT.TGraphErrors(len(Xaxis), Xaxis, p1corr, np.zeros(len(Xaxis)), p1correrr)
    gppcorr.SetMarkerColor(ROOT.kTeal+4)
    gppcorr.SetMarkerStyle(8)
    '''
    linearcorr = ROOT.TF1("linearcorr", "[0] * x^3 + [1]", np.min(Xaxis), np.max(Xaxis))
    linearcorr.SetParameters(1e-08, 0)#初期パラメータはp0は0にすべきではないなぜなら関数系が変わりlocal minimumに引っかかりがち(あとp0＋ー気をつける)
    linearcorr.SetLineColor(ROOT.kSpring)
    gppcorr.Fit(linearcorr, "R")
    linearcorr.Draw("same")
    '''
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
    ppoutFile = ROOT.TFile(f"{cpp_dir}PeakPositionpm50.root", "RECREATE")
    cpp.Write()
    ppoutFile.Close()
    print(f"Successfuly Saved : {ppoutFile}")
    
    
    #ROOT.gApplication.Run()
    
    
    
    
   
    