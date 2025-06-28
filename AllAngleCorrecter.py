#!/usr/bin/env python3
import ROOT
import numpy as np
import os
import math
#三次元空間でのベクトルの表記方法は(x,y,z)の順
Axislist = []
p1corr_list, p1correrr_list = [], [] #PeakPosition 角度補正後(センサーのみ)
p1corrall_list, p1corrallerr_list = [], [] #PeakPosition 角度補正後(レールも)
p1meas_list, p1measerr_list = [], []#PeakPosition 測定
p1opera_list, p1operaerr_list = [], []#PeakPosition 角度補正前
#回転を考える時は各軸の右ねじの回転を考える。従って、xy平面,yz平面,zx平面
def rotM_zx(deg):
    #zx平面における回転行列
    theta = np.deg2rad(deg)
    return np.array([[np.cos(theta), 0, np.sin(theta)],
                     [0, 1, 0],
                     [-np.sin(theta), 0, np.cos(theta)]])
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
        
#レールの角度も
def plot_allanglecorrected(OPERAfile,MEASfile):
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
        Bzstdlist.append(getattr(meastree, "Zstdev") * (3/5) )
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
    
    #測定データに対して角度補正をかける。各点の測定磁場ベクトルと各センサーの法線ベクトルの内積
    Bxcorr = np.zeros(len(Bxmeas))#"測定"を補正したデータを格納する(センサーのみ)
    Bycorr = np.zeros(len(Bymeas))
    Bzcorr = np.zeros(len(Bzmeas))
    for i in range(len(Bxmeas)):
        B_vec = np.array([[Bxmeas[i]],
                     [Bymeas[i]],
                     [Bzmeas[i]]])#測定磁場の三次元ベクトル
        #レールの角度も補正
        theta_Bx = 0.4#Bxセンサーのzx平面での角度(レール込み)
        phi_Bx = -0.5#Bxセンサーのxy平面での角度
        theta_By = 0.2#Byセンサーのyz平面での角度
        phi_By = 0.2#Byセンサーのxy平面での角度
        theta_Bz = 1.8#Bzセンサーのzx平面での角度
        phi_Bz = 0.4#Bzセンサーのyz平面での角度
        nBx = rotM_xy(phi_Bx) @ rotM_zx(theta_Bx) @ np.array([[1],[0],[0]])
        nBy = rotM_xy(phi_By) @ rotM_yz(theta_By) @ np.array([[0],[1],[0]])
        nBz = rotM_yz(phi_Bz) @ rotM_zx(theta_Bz) @ np.array([[0],[0],[1]])
        rotM = np.hstack((nBx, nBy, nBz)).T #１行目nBxの横ベクトル、２行目nByの横ベクトル、3行目nBzの横ベクトル
        #逆行列を縦ベクトル(Bxmeas,Bymeas,Bzmeas)にかけて測定を補正
        invM = np.linalg.inv(rotM)
        Bcorr = invM @ B_vec
        Bxcorr[i] = Bcorr[0,0]
        Bycorr[i] = Bcorr[1,0]
        Bzcorr[i] = Bcorr[2,0]
    #測定データに対してレールの角度込みで角度補正をかける。各点の測定磁場ベクトルと各センサーの法線ベクトルの内積
    Bxcorr_all = np.zeros(len(Bxmeas))#"測定"を補正したデータを格納する(レールの角度も)
    Bycorr_all = np.zeros(len(Bymeas))
    Bzcorr_all = np.zeros(len(Bzmeas))
    Bangleerr = np.zeros(len(Bxmeas))
    
    #レールの角度
    #angle_xy = np.zeros(len(Bxmeas))#XY平面
    angle_zx = np.zeros(len(Bxmeas))#XZ平面
    angle_yz = np.zeros(len(Bxmeas))#YZ平面
    deltaX = np.zeros(len(Xmeas))
    deltaY = np.zeros(len(Ymeas))
    deltaZ = np.zeros(len(Zmeas))
    #レールの角度を計算
    for i in range(len(Xmeas)):
        if i==0:
            angle_zx[i] = 0
            angle_yz[i] = 0
            deltaX[i] = 0
            deltaY[i] = 0
            deltaZ[i] = 35.0
        if 0 < i < 212 :
            #レールの角度も補正(XY平面でのレールによる角度はトラック不可)
            deltaX[i] = Xmeas[i]-Xmeas[i-1]
            deltaY[i] = Ymeas[i]-Ymeas[i-1]
            deltaZ[i] = Zmeas[i]-Zmeas[i-1]
            
            #angle_zx[i] = 90+np.degrees(np.arctan2(deltaZ[i], deltaX[i]))
            angle_zx[i] = np.degrees(np.arctan(deltaX[i]/deltaZ[i]))
            angle_yz[i] = np.degrees(np.arctan(deltaY[i]/(-deltaZ[i])))
            #angle_yz[i] = 90+np.degrees(np.arctan2(deltaZ[i], deltaY[i]))
        if i ==212: #test2024062606のid=212と211の座標点が同じ->おそらく同じとこで2回座標測定した
            deltaX[i] = Xmeas[i]-Xmeas[i-1]
            deltaY[i] = Ymeas[i]-Ymeas[i-1]
            deltaZ[i] = 35.0
            angle_zx[i] = 0
            angle_yz[i] = 0
    #測定を補正
    for i in range(len(Bxmeas)):
        #print(f"{i} : {deltaZ[i]}")
        B_vec = np.array([[Bxmeas[i]],
                     [Bymeas[i]],
                     [Bzmeas[i]]])#測定磁場の三次元ベクトル
      #トータルでのセンサーの角度(センサー＋レール)
        theta_Bx = 0.4#Bxセンサーのzx平面での角度
        phi_Bx = -0.5 #Bxセンサーのxy平面での角度
        theta_By = 0.2#Byセンサーのyz平面での角度
        phi_By = 0.2 #Byセンサーのxy平面での角度
        theta_Bz = 1.8 #Bzセンサーのzx平面での角度
        phi_Bz = 0.4 #Bzセンサーのyz平面での角度
        nBx = rotM_zx(angle_zx[i]) @ rotM_xy(phi_Bx) @ rotM_zx(theta_Bx) @ np.array([[1],[0],[0]])
        nBy = rotM_yz(angle_yz[i]) @ rotM_xy(phi_By) @ rotM_yz(theta_By) @ np.array([[0],[1],[0]])
        nBz = rotM_yz(angle_yz[i]) @ rotM_zx(angle_zx[i]) @ rotM_yz(phi_Bz) @ rotM_zx(theta_Bz) @ np.array([[0],[0],[1]])
        rotM = np.hstack((nBx, nBy, nBz)).T #１行目nBxの横ベクトル、２行目nByの横ベクトル、3行目nBzの横ベクトル
        #逆行列を縦ベクトル(Bxmeas,Bymeas,Bzmeas)にかけて測定を補正
        invM = np.linalg.inv(rotM)
        identity_check = invM @ rotM
        if not np.allclose(identity_check, np.identity(3), atol=1e-10):
            raise ValueError(f"The inversed matrix might be WRONG !! i = {i}\n{identity_check}")
        Bcorr = invM @ B_vec
        Bxcorr_all[i] = Bcorr[0,0]
        Bycorr_all[i] = Bcorr[1,0]
        Bzcorr_all[i] = Bcorr[2,0]
        
    #print(f"{i}:{angle_zx[i]}")
    #print(f"{angle_yz[111]}")
    #print(f"{rotM_zx(angle_zx[111]) @ rotM_yz(angle_yz[111]) @ rotM_yz(phi_Bz) @ rotM_zx(theta_Bz) @ np.array([[0],[0],[1]])}")
    #print(f"{rotM_yz(angle_yz[111]) @ rotM_xy(phi_By) @ rotM_yz(theta_By) @ np.array([[0],[1],[0]])}")
    #print(f"{np.argmax(angle_zx)} : {np.max(angle_zx)}")
    #print(f"deltaX : {X[111]} - {X[110]} = {deltaX[111]}")
    #print(f"deltaY : {Y[111]} - {Y[110]} = {deltaY[111]}")
    #print(f"deltaZ : {Z[111]} - {Z[110]} = {deltaZ[111]}")     
    
    #peakposition分布のプロット作成に必要
    Axislist.append(np.mean(Y))
    
    #Draw Graphs from here
        #Graph configuration
    ROOT.gStyle.SetLabelSize(0.04, "X") 
    ROOT.gStyle.SetLabelSize(0.04, "Y")
    ROOT.gStyle.SetTitleSize(0.05, "X")
    ROOT.gStyle.SetTitleSize(0.05, "Y")
    ROOT.gStyle.SetLegendTextSize(0.03)
    #OPERA
    gOPERAx = ROOT.TGraph(len(Z), Z, Bx)
    gOPERAx.SetMarkerStyle(8)
    gOPERAx.SetMarkerColor(ROOT.kRed+2)
    gOPERAy = ROOT.TGraph(len(Z), Z, By)
    gOPERAy.SetMarkerStyle(8)
    gOPERAy.SetMarkerColor(ROOT.kRed+2)
    gOPERAz = ROOT.TGraphErrors(len(Zmeas), Zmeas, Bz, np.zeros(len(Z)), Bzstd)
    gOPERAz.SetMarkerStyle(8)
    gOPERAz.SetMarkerColor(ROOT.kRed+2)
    gOPERAz.SetLineWidth(2)
    gOPERAz.SetLineColor(ROOT.kRed+2)
        #Fitting
    BzMin = np.min(Bz)
    BzMinZ = Z[np.argmin(Bz)]
    qf = ROOT.TF1("qf", "[0]*(x - [1])^2 + [2]", BzMinZ - 50, BzMinZ + 50)#Fitting範囲は頂点＋ー150mm
    qf.SetParameters(1e-7, BzMinZ, BzMin)  # 初期パラメータはp0は0にすべきではないなぜなら関数系が変わりlocal minimumに引っかかりがち
    qf.SetLineColor(ROOT.kRed+2)
    gOPERAz.Fit(qf, "R")
    p1opera_list.append(qf.GetParameter(1))
    p1operaerr_list.append(qf.GetParError(1))
    qf.Draw("same")
    #測定値補正なし
    gmeasx = ROOT.TGraphErrors(len(Z), Z, Bxmeas, np.zeros(len(Z)), Bxstd)
    gmeasx.SetMarkerStyle(7)
    gmeasx.SetMarkerColor(ROOT.kBlue)
    gmeasy = ROOT.TGraphErrors(len(Z), Z, Bymeas, np.zeros(len(Z)), Bystd)
    gmeasy.SetMarkerStyle(8)
    gmeasy.SetMarkerColor(ROOT.kBlue)
    gmeasz = ROOT.TGraphErrors(len(Zmeas), Zmeas, Bzmeas, np.zeros(len(Z)), Bzstd)
    gmeasz.SetMarkerStyle(8)
    gmeasz.SetMarkerColor(ROOT.kBlue)
    gmeasz.SetLineWidth(2)
    gmeasz.SetLineColor(ROOT.kBlue)
        #Fitting
    BzmeasMin = np.min(Bzmeas)
    BzmeasMinZ = Zmeas[np.argmin(Bzmeas)]
    qfmeas = ROOT.TF1("qfmeas", "[0]*(x - [1])^2 + [2]", BzmeasMinZ - 50, BzmeasMinZ + 50)
    qfmeas.SetParameters(1e-7, BzmeasMinZ, BzmeasMin)
    qfmeas.SetLineColor(ROOT.kBlue)
    gmeasz.Fit(qfmeas, "R")
    p1meas_list.append(qfmeas.GetParameter(1))
    p1measerr_list.append(qfmeas.GetParError(1))
    qfmeas.Draw("same")
    #測定値補正あり(センサーの角度のみ)
    gcorrx = ROOT.TGraph(len(Zmeas), Zmeas, Bxcorr)
    gcorrx.SetMarkerStyle(8)
    gcorrx.SetMarkerColor(ROOT.kAzure+7)
    gcorry = ROOT.TGraph(len(Zmeas), Zmeas, Bycorr)
    gcorry.SetMarkerStyle(8)
    gcorry.SetMarkerColor(ROOT.kAzure+7)
    gcorrz = ROOT.TGraphErrors(len(Zmeas), Zmeas, Bzcorr, np.zeros(len(Zmeas)), Bzstd)#測定補正につけるエラーは測定のものを内挿
    gcorrz.SetMarkerStyle(8)
    gcorrz.SetMarkerColor(ROOT.kAzure+7)
    gcorrz.SetLineWidth(2)
    gcorrz.SetLineColor(ROOT.kAzure+7)
        #Fitting
    BzcorrMin = np.min(Bzcorr)#Bzなので分布の谷底、最小値になる
    BzcorrMinZ = Zmeas[np.argmin(Bzcorr)]
    #qfcorr = ROOT.TF1("qfcorr", "[0]*(x - [1])^2 + [2]", Z[np.argmin(Bzcorr) - 1], Z[np.argmin(Bzcorr) + 1])#頂点と両サイド3点
    qfcorr = ROOT.TF1("qfcorr", "[0]*(x - [1])^2 + [2]", BzcorrMinZ - 50, BzcorrMinZ + 50)#Fitting範囲は頂点と両サイド3点
    qfcorr.SetParameters(1e-7, BzcorrMinZ, BzcorrMin)#初期パラメータはp0は0にすべきではないなぜなら関数系が変わりlocal minimumに引っかかりがち(あとp0＋ー気をつける)
    qfcorr.SetLineColor(ROOT.kAzure+7)
    gcorrz.Fit(qfcorr, "R")
    p1corr_list.append(qfcorr.GetParameter(1))
    p1correrr_list.append(qfcorr.GetParError(1))
    qfcorr.Draw("same")
    #測定値補正あり(レールアングルも)
    gcorrallx = ROOT.TGraph(len(Zmeas), Zmeas, Bxcorr_all)
    gcorrallx.SetMarkerStyle(8)
    gcorrallx.SetMarkerColor(ROOT.kOrange+10)
    gcorrally = ROOT.TGraph(len(Zmeas), Zmeas, Bycorr_all)
    gcorrally.SetMarkerStyle(8)
    gcorrally.SetMarkerColor(ROOT.kOrange+10)
    gcorrallz = ROOT.TGraphErrors(len(Zmeas), Zmeas, Bzcorr_all, np.zeros(len(Zmeas)), Bzstd)#測定補正につけるエラーは測定のものを内挿
    gcorrallz.SetMarkerStyle(8)
    gcorrallz.SetMarkerColor(ROOT.kOrange+10)
    gcorrallz.SetLineWidth(2)
    gcorrallz.SetLineColor(ROOT.kOrange+10)
        #Fitting
    BzcorrallMin = np.min(Bzcorr_all)#Bzなので分布の谷底、最小値になる
    BzcorrallMinZ = Zmeas[np.argmin(Bzcorr)]
    qfcorrall = ROOT.TF1("qfcorr", "[0]*(x - [1])^2 + [2]", BzcorrallMinZ - 50, BzcorrallMinZ + 50)#Fitting範囲は頂点と両サイド3点
    qfcorrall.SetParameters(1e-7, BzcorrallMinZ, BzcorrallMin)#初期パラメータはp0は0にすべきではないなぜなら関数系が変わりlocal minimumに引っかかりがち(あとp0＋ー気をつける)
    qfcorrall.SetLineColor(ROOT.kOrange+10)
    gcorrallz.Fit(qfcorrall, "R")
    p1corrall_list.append(qfcorr.GetParameter(1))
    p1corrallerr_list.append(qfcorr.GetParError(1))
    qfcorrall.Draw("same")
    
    #キャンバス描画
    cBzcomp = ROOT.TCanvas("cBzcomp", f"{axname}_Bzcomp", 1200, 800)# Bz成分の分布比較
    
    #Bz分布OPERA＆OPERA corr (Fittingあり)
    cBzcomp.cd()
    mgBzcomp = ROOT.TMultiGraph()
    mgBzcomp.Add(gcorrz)
    mgBzcomp.Add(gcorrallz)
    mgBzcomp.Add(gOPERAz)
    mgBzcomp.Add(gmeasz)
    mgBzcomp.SetTitle(f"{axname} : Bz comp")
    mgBzcomp.GetXaxis().SetTitle("Z (mm)")
    mgBzcomp.GetYaxis().SetTitle("Bz (T)")
    mgBzcomp.Draw("AP")
    ROOT.gPad.SetGrid(1,1)
    legendBzcomp = ROOT.TLegend(0.7, 0.1, 0.9, 0.25)
    legendBzcomp.AddEntry(gOPERAz, "OPERA", "p")
    legendBzcomp.AddEntry(gcorrz, "MEAS corr", "p")
    legendBzcomp.AddEntry(gmeasz, "Measured", "p")
    legendBzcomp.AddEntry(gcorrallz, "Rail corr", "p")
    legendBzcomp.SetFillStyle(0)
    legendBzcomp.Draw()
    cBzcomp.Update()
    #cBzcomp_dir = "/Users/shohtatakami/physics/COMETDS/DS189A_944turns_Integration/Allcorr/Bzcomp/"
    #outputFile = ROOT.TFile(f"{cBzcomp_dir}{axname}.root", "RECREATE")
    #cBzcomp.Write()
    #outputFile.Close()
    #print(f"Successfully Saved : {outputFile}")

    #測定点分布
    '''
    crailangle = ROOT.TCanvas("railangle", f"{axname}_railangle", 1200, 800)#レール角度の分布
    crailangle.Divide(2,2)
    crailangle.cd(1)
    gX = ROOT.TGraph(len(Zmeas), Zmeas, Xmeas)
    gX.SetTitle(f"{axname} : ZX coordinate")
    gX.SetMarkerStyle(7)
    gX.SetMarkerColor(ROOT.kBlack)
    gX.GetXaxis().SetTitle("Z (mm)")
    gX.GetYaxis().SetTitle("X (mm)")
    gX.Draw("AP")
    gXnegedge = ROOT.TLine(-1660, gX.GetYaxis().GetXmin(), -1660, gX.GetYaxis().GetXmax())
    gXnegedge.SetLineStyle(2)          # 点線
    gXnegedge.SetLineColor(ROOT.kGray+2)  
    gXnegedge.SetLineWidth(2)
    gXnegedge.Draw("same")
    gXposedge = ROOT.TLine(1660, gX.GetYaxis().GetXmin(), 1660, gX.GetYaxis().GetXmax())
    gXposedge.SetLineStyle(2)          # 点線
    gXposedge.SetLineColor(ROOT.kGray+2)  
    gXposedge.SetLineWidth(2)
    gXposedge.Draw("same")
    crailangle.cd(2)
    gY = ROOT.TGraph(len(Zmeas), Zmeas, Ymeas)
    gY.SetTitle(f"{axname} : YZ coordinate")
    gY.SetMarkerStyle(7)
    gY.SetMarkerColor(ROOT.kBlack)
    gY.GetXaxis().SetTitle("Z (mm)")
    gY.GetYaxis().SetTitle("Y (mm)")
    gY.Draw("AP")
    gYnegedge = ROOT.TLine(-1660, gY.GetYaxis().GetXmin(), -1660, gY.GetYaxis().GetXmax())
    gYnegedge.SetLineStyle(2)          # 点線
    gYnegedge.SetLineColor(ROOT.kGray+2)  
    gYnegedge.SetLineWidth(2)
    gYnegedge.Draw("same")
    gYposedge = ROOT.TLine(1660, gY.GetYaxis().GetXmin(), 1660, gY.GetYaxis().GetXmax())
    gYposedge.SetLineStyle(2)          # 点線
    gYposedge.SetLineColor(ROOT.kGray+2)  
    gYposedge.SetLineWidth(2)
    gYposedge.Draw("same")
    
    #レール角度分布
    crailangle.cd(3)
    gangleXZ = ROOT.TGraph(len(Zmeas), Zmeas, angle_zx)
    gangleXZ.SetTitle(f"{axname} : ZX Angle")
    gangleXZ.SetMarkerStyle(7)
    gangleXZ.SetMarkerColor(ROOT.kBlack)
    gangleXZ.GetXaxis().SetTitle("Z (mm)")
    gangleXZ.GetYaxis().SetTitle("angle_ZX (#circ)")
    gangleXZ.Draw("AP")
    gangleXZnegedge = ROOT.TLine(-1660, gangleXZ.GetYaxis().GetXmin(), -1660, gangleXZ.GetYaxis().GetXmax())
    gangleXZnegedge.SetLineStyle(2)          # 点線
    gangleXZnegedge.SetLineColor(ROOT.kGray+2)  
    gangleXZnegedge.SetLineWidth(2)
    gangleXZnegedge.Draw("same")
    gangleXZposedge = ROOT.TLine(1660, gangleXZ.GetYaxis().GetXmin(), 1660, gangleXZ.GetYaxis().GetXmax())
    gangleXZposedge.SetLineStyle(2)          # 点線
    gangleXZposedge.SetLineColor(ROOT.kGray+2)  
    gangleXZposedge.SetLineWidth(2)
    gangleXZposedge.Draw("same")
    crailangle.cd(4)
    gangleYZ = ROOT.TGraph(len(Zmeas), Zmeas, angle_yz)
    gangleYZ.SetTitle(f"{axname} : YZ Angle")
    gangleYZ.SetMarkerStyle(7)
    gangleYZ.SetMarkerColor(ROOT.kBlack)
    gangleYZ.GetXaxis().SetTitle("Z (mm)")
    gangleYZ.GetYaxis().SetTitle("angle_YZ (#circ)")
    gangleYZ.Draw("AP")
    gangleYZnegedge = ROOT.TLine(-1660, gangleYZ.GetYaxis().GetXmin(), -1660, gangleYZ.GetYaxis().GetXmax())
    gangleYZnegedge.SetLineStyle(2)          # 点線
    gangleYZnegedge.SetLineColor(ROOT.kGray+2)  
    gangleYZnegedge.SetLineWidth(2)
    gangleYZnegedge.Draw("same")
    gangleYZposedge = ROOT.TLine(1660, gangleYZ.GetYaxis().GetXmin(), 1660, gangleYZ.GetYaxis().GetXmax())
    gangleYZposedge.SetLineStyle(2)          # 点線
    gangleYZposedge.SetLineColor(ROOT.kGray+2)  
    gangleYZposedge.SetLineWidth(2)
    gangleYZposedge.Draw("same")
    crailangle.Update()
    #crailangle_dir = "/Users/shohtatakami/physics/COMETDS/DS189A_944turns_Integration/AllAnglecorr/railangle/"
    crailangle_dir = "/Users/shohtatakami/physics/COMETDS/DS189A_944turns_Integration/Allcorr/railangle/"
    railangleFile = ROOT.TFile(f"{crailangle_dir}{axname}_railangle.root", "RECREATE")
    crailangle.Write()
    railangleFile.Close()
    print(f"Successfully saved : {railangleFile}")
    '''
if __name__=='__main__':
    operafile_directory =  "/Users/shohtatakami/physics/COMETDS/DS189A_944turns_Integration/DS189A_944turns_FieldIntegration_LaserTrackerPos/root/"
    measfile_directory = "/Users/shohtatakami/physics/COMETDS/newScan/"#座標追加
    #PeakPositionX分布
    file_pairs_Y0Plane= [
        ("X_-200Y_0_Map.root","test20240627-3.root"),
        ("X_-500Y_0_Map.root","test20240627-4.root"),
        ("X_-650Y_0_Map.root","test20240627-7.root"),
        ("X_-760Y_0_Map.root","test20240627-6.root"),
        ("X_0Y_0_Map.root","test20240626-6.root"),
        ("X_200Y_0_Map.root","test20240627-0.root"),
        ("X_500Y_0_Map.root","test20240627-1.root"),
        ("X_650Y_0_Map.root","test20240627-8.root"),
        ("X_760Y_0_Map.root","test20240627-2rev.root")#rev : the version deleted the last several strange lines 
    ]
    file_pairs_X0Plane= [
        ("X_0Y_-200_Map.root", "test20240710-4.root"),
        ("X_0Y_-729_Map.root", "test20240628-2.root"),
        ("X_0Y_200_Map.root", "test20240710-1.root"),
        ("X_0Y_780_Map.root", "test20240628-0rev.root"),#rev : the version deleted the last several strange lines 
    ]
    
    for operafilename, measfilename in file_pairs_X0Plane:
        OPERAfile = os.path.join(operafile_directory, operafilename)
        MEASfile = os.path.join(measfile_directory, measfilename)
        plot_allanglecorrected(OPERAfile, MEASfile)
        #factorypillar(OPERAfile, MEASfile)
    #peak position 取得
    
    Axis = np.array(Axislist, dtype = np.float64)
    p1meas = np.array(p1meas_list, dtype = np.float64)
    p1measerr = np.array(p1measerr_list, dtype = np.float64)
    p1opera = np.array(p1opera_list, dtype = np.float64)
    p1operaerr = np.array(p1operaerr_list, dtype = np.float64)
    p1corr = np.array(p1corr_list, dtype = np.float64)
    p1correrr = np.array(p1correrr_list, dtype = np.float64)
    p1corrall = np.array(p1corrall_list, dtype = np.float64)
    p1corrallerr = np.array(p1corrallerr_list, dtype = np.float64)
    cpp = ROOT.TCanvas("cpp", "Peak Position", 800, 600)
    gppmeas = ROOT.TGraphErrors(len(Axis), Axis, p1meas, np.zeros(len(Axis)), p1measerr)
    gppmeas.SetMarkerColor(ROOT.kBlue)
    gppmeas.SetMarkerStyle(8)
    gppmeas.SetLineWidth(2)
    gppmeas.SetLineColor(ROOT.kBlue)
    
    #linearmeas = ROOT.TF1("linearmeas", "[0] * x^3 + [1]", np.min(Xaxis), np.max(Xaxis))
    #linearmeas.SetParameters(1e-08, 0)#初期パラメータはp0は0にすべきではないなぜなら関数系が変わりlocal minimumに引っかかりがち(あとp0+ー気をつける)
    #linearmeas.SetLineColor(ROOT.kBlue)
    #gppmeas.Fit(linearmeas, "R")
    #linearmeas.Draw("same")
    
    gppopera = ROOT.TGraphErrors(len(Axis), Axis, p1opera, np.zeros(len(Axis)), p1operaerr)
    #gppopera = ROOT.TGraphErrors(len(Xaxis), Xaxis, p1opera, np.zeros(len(Xaxis)), np.zeros(len(Xaxis)))
    gppopera.SetMarkerColor(ROOT.kRed+2)
    gppopera.SetMarkerStyle(8)
    gppopera.SetLineWidth(2)
    gppopera.SetLineColor(ROOT.kRed+2)
    gppcorr = ROOT.TGraphErrors(len(Axis), Axis, p1corr, np.zeros(len(Axis)), p1correrr)
    gppcorr.SetMarkerColor(ROOT.kAzure+7)
    gppcorr.SetMarkerStyle(8)
    gppcorr.SetLineWidth(2)
    gppcorr.SetLineColor(ROOT.kAzure+7)
    gppcorrall = ROOT.TGraphErrors(len(Axis), Axis, p1corrall, np.zeros(len(Axis)), p1corrallerr)
    gppcorrall.SetMarkerColor(ROOT.kOrange+10)
    gppcorrall.SetMarkerStyle(8)
    gppcorrall.SetLineWidth(2)
    gppcorrall.SetLineColor(ROOT.kOrange+10)
    
    #linearcorr = ROOT.TF1("linearcorr", "[0] * x^3 + [1]", np.min(Xaxis), np.max(Xaxis))
    #linearcorr.SetParameters(1e-08, 0)#初期パラメータはp0は0にすべきではないなぜなら関数系が変わりlocal minimumに引っかかりがち(あとp0+ー気をつける)
    #linearcorr.SetLineColor(ROOT.kSpring)
    #gppcorr.Fit(linearcorr, "R")
    #linearcorr.Draw("same")
    
    mgpp = ROOT.TMultiGraph()
    mgpp.Add(gppmeas)
    mgpp.Add(gppopera)
    mgpp.Add(gppcorr)
    mgpp.Add(gppcorrall)
    mgpp.GetXaxis().SetTitle("Y (mm)")#peakposition X分布
    #mgpp.GetXaxis().SetTitle("Y (mm)")#peakposition Y分布
    mgpp.GetYaxis().SetTitle("Peak Position (mm)")
    mgpp.Draw("AP")
    ROOT.gPad.SetGrid(1,1)
    legendpp = ROOT.TLegend(0.7, 0.75, 0.9, 0.9)
    legendpp.AddEntry(gppmeas, "Measured", "p")
    legendpp.AddEntry(gppopera, "OPERA", "p")
    legendpp.AddEntry(gppcorr, "MEAS Corr", "p")
    legendpp.AddEntry(gppcorrall, "Rail Corr", "p")
    legendpp.SetFillStyle(0)
    legendpp.Draw()
    cpp.Update()

    #cpp_dir = "/Users/shohtatakami/physics/COMETDS/DS189A_944turns_Integration/AllAnglecorr/"
    cpp_dir = "/Users/shohtatakami/Desktop/"
    ppoutFile = ROOT.TFile(f"{cpp_dir}PeakPositionY.root", "RECREATE")
    cpp.Write()
    ppoutFile.Close()
    print(f"Successfuly Saved : {ppoutFile}")
    #ROOT.gApplication.Run()
    
    
    
    
   
    