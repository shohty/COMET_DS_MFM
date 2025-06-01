#!/usr/bin/env python3
import ROOT
import numpy as np
import os
import math
#三次元空間でのベクトルの表記方法は(x,y,z)の順
Xaxlist = []
#Yaxlist = []
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
def DeltaofAngle(OPERAfile,MEASfile):
    file0= ROOT.TFile.Open(OPERAfile, "READ")#OPERAシミュレーションファイルの読み込み
    file1 = ROOT.TFile.Open(MEASfile, "READ")#OPERAシミュレーションファイルの読み込み
    axname =  os.path.splitext(os.path.basename(OPERAfile))[0]
    print(f"{axname}")
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
    
    #測定データに対してレールの角度込みで角度補正をかける。各点の測定磁場ベクトルと各センサーの法線ベクトルの内積
    Bxcorr= np.zeros(len(Bxmeas))#"測定"を補正したデータを格納する(レールの角度も)
    Bycorr= np.zeros(len(Bymeas))
    Bzcorr= np.zeros(len(Bzmeas))
    deltaBx_pl = np.zeros(len(Bxmeas))#誤差プラス側のBxとBxcorrの差
    deltaBx_mi = np.zeros(len(Bxmeas))
    deltaBy_pl = np.zeros(len(Bxmeas))#誤差プラス側のByとBycorrの差
    deltaBy_mi = np.zeros(len(Bxmeas))
    deltaBz_pl = np.zeros(len(Bxmeas))#誤差プラス側のBzとBzcorrの差
    deltaBz_mi = np.zeros(len(Bxmeas))
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
        Bxcorr[i] = Bcorr[0,0]
        Bycorr[i] = Bcorr[1,0]
        Bzcorr[i] = Bcorr[2,0]
        
        #角度の誤差の磁場への寄与を計算
        #ZX平面(Byは影響なし)
        nBx_pl = rotM_zx(angle_zx[i]) @ rotM_xy(phi_Bx) @ rotM_zx(theta_Bx + 0.1) @ np.array([[1],[0],[0]]) #誤差のプラス側(Bxセンサー法線ベクトル)
        nBx_mi = rotM_zx(angle_zx[i]) @ rotM_xy(phi_Bx) @ rotM_zx(theta_Bx - 0.1) @ np.array([[1],[0],[0]]) #誤差のマイナス側
        nBy_pl = rotM_yz(angle_yz[i]) @ rotM_xy(phi_By + 0.1) @ rotM_yz(theta_By) @ np.array([[0],[1],[0]])
        nBy_mi = rotM_yz(angle_yz[i]) @ rotM_xy(phi_By - 0.1) @ rotM_yz(theta_By) @ np.array([[0],[1],[0]])
        nBz_pl = rotM_yz(angle_yz[i]) @ rotM_zx(angle_zx[i]) @ rotM_yz(phi_Bz ) @ rotM_zx(theta_Bz + 0.1) @ np.array([[0],[0],[1]])#誤差のプラス側(Bzセンサー法線ベクトル)
        nBz_mi = rotM_yz(angle_yz[i]) @ rotM_zx(angle_zx[i]) @ rotM_yz(phi_Bz ) @ rotM_zx(theta_Bz - 0.1) @ np.array([[0],[0],[1]])#誤差のマイナス側
        deltaBz_pl[i] = (nBz_pl.T @ B_vec).item() - Bzcorr[i].item()#誤差のプラス側のBxの値_____________.item()をつけないと1×1のarrayとして扱われるためwarningが出る
        deltaBz_mi[i] = (nBz_mi.T @ B_vec).item() - Bzcorr[i].item()#誤差のマイナス側のBxの値
        deltaBx_pl[i] = (nBx_pl.T @ B_vec).item() - Bxcorr[i].item()#誤差のプラス側のBzの値
        deltaBx_mi[i] = (nBx_mi.T @ B_vec).item() - Bxcorr[i].item()#誤差のマイナス側のBzの値
    #print(type(deltaBx_pl[i]))
    #print(deltaBx_pl[i].shape)
    #Draw Graphs form here
            #Graph configuration
    ROOT.gStyle.SetLabelSize(0.06, "X") 
    ROOT.gStyle.SetLabelSize(0.06, "Y")
    ROOT.gStyle.SetTitleSize(0.05, "X")
    ROOT.gStyle.SetTitleSize(0.05, "Y")
    ROOT.gStyle.SetLegendTextSize(0.03)
    
    c_Bdelta = ROOT.TCanvas(f"delta{axname}", f"delta{axname}",1200,900)
    c_Bdelta.Divide(1,2)
        
    g_Bxpldelta = ROOT.TGraph(len(Z), Z, deltaBx_pl)
    g_Bxpldelta.SetMarkerStyle(8)
    g_Bxpldelta.SetMarkerColor(ROOT.kRed)
    g_Bxmidelta = ROOT.TGraph(len(Z), Z, deltaBx_mi)
    g_Bxmidelta.SetMarkerStyle(8)
    g_Bxmidelta.SetMarkerColor(ROOT.kBlue)
    
    g_Bzpldelta = ROOT.TGraph(len(Z), Z, deltaBz_pl)
    g_Bzpldelta.SetMarkerStyle(8)
    g_Bzpldelta.SetMarkerColor(ROOT.kRed)
    g_Bzmidelta = ROOT.TGraph(len(Z), Z, deltaBz_mi)
    g_Bzmidelta.SetMarkerStyle(8)
    g_Bzmidelta.SetMarkerColor(ROOT.kBlue)
    
    c_Bdelta.cd(1)
    mg_Bxdelta = ROOT.TMultiGraph()
    mg_Bxdelta.Add(g_Bxpldelta)
    mg_Bxdelta.Add(g_Bxmidelta)
    mg_Bxdelta.SetTitle(f"{axname} : Bx delta")#
    mg_Bxdelta.GetXaxis().SetTitle("Z (mm)")
    mg_Bxdelta.GetYaxis().SetTitle("deltaB (T)")
    mg_Bxdelta.Draw("AP")
    ROOT.gPad.SetGrid(1,1)
    Bxdelta_legend = ROOT.TLegend(0.7, 0.1, 0.9, 0.25)
    Bxdelta_legend.AddEntry(g_Bxpldelta, "#theta + #delta#theta", "p")
    Bxdelta_legend.AddEntry(g_Bxmidelta, "#theta - #delta#theta", "p")
    Bxdelta_legend.SetFillStyle(0)
    Bxdelta_legend.Draw()
    c_Bdelta.cd(2)
    mg_Bzdelta = ROOT.TMultiGraph()
    mg_Bzdelta.Add(g_Bzpldelta)
    mg_Bzdelta.Add(g_Bzmidelta)
    mg_Bzdelta.SetTitle(f"{axname} : Bz delta")
    mg_Bzdelta.GetXaxis().SetTitle("Z (mm)")
    mg_Bzdelta.GetYaxis().SetTitle("deltaB (T)")
    mg_Bzdelta.Draw("AP")
    ROOT.gPad.SetGrid(1,1)
    Bzdelta_legend = ROOT.TLegend(0.7, 0.1, 0.9, 0.25)
    Bzdelta_legend.AddEntry(g_Bzpldelta, "#theta + #delta#theta", "p")
    Bzdelta_legend.AddEntry(g_Bzmidelta, "#theta - #delta#theta", "p")
    Bzdelta_legend.SetFillStyle(0)
    Bzdelta_legend.Draw()

    c_Bdelta.Update()
    
    c_Bdelta_dir = "/Users/shohtatakami/physics/COMETDS/ErrorBudget/delta_pitch/"#保存ディレクトリは平面に応じて分ける
    delta_File = ROOT.TFile(f"{c_Bdelta_dir}{axname}_sensor.root", "RECREATE")
    c_Bdelta.Write()
    delta_File.Close()
    print(f"Successfully saved : {delta_File}")  
    
if __name__=='__main__':
    operafile_directory =  "/Users/shohtatakami/physics/COMETDS/DS189A_944turns_Integration/DS189A_944turns_FieldIntegration_LaserTrackerPos/root/"
    measfile_directory = "/Users/shohtatakami/physics/COMETDS/newScan/"#座標追加
    #PeakPositionX分布
    file_pairs= [
        ("X_-200Y_0_Map.root","test20240627-3.root"),
        ("X_-500Y_0_Map.root","test20240627-4.root"),
        ("X_-650Y_0_Map.root","test20240627-7.root"),
        ("X_-760Y_0_Map.root","test20240627-6.root"),
        ("X_0Y_0_Map.root","test20240626-6.root"),
        ("X_200Y_0_Map.root","test20240627-0.root"),
        ("X_500Y_0_Map.root","test20240627-1.root"),
        ("X_650Y_0_Map.root","test20240627-8.root"),
        ("X_760Y_0_Map.root","test20240627-2rev.root"),#rev : the version deleted the last several strange lines 
    ]
    '''
    #鉄柱ありOPERAファイル
    file_pairs= [
        ("DS189A_944turns_FactoryPillar_FieldIntegration_X-200_Y0.root","test20240627-3.root"),
        ("DS189A_944turns_FactoryPillar_FieldIntegration_X-500_Y0.root","test20240627-4.root"),
        ("DS189A_944turns_FactoryPillar_FieldIntegration_X-650_Y0.root","test20240627-7.root"),
        ("DS189A_944turns_FactoryPillar_FieldIntegration_X-760_Y0.root","test20240627-6.root"),
        #("X_0Y_-200_Map.root","test20240710-4.root"),
        #("X_0Y_-729_Map.root","test20240628-2.root"),
        ("DS189A_944turns_FactoryPillar_FieldIntegration_X0_Y0.root","test20240626-6.root"),
        #("X_0Y_200_Map.root","test20240710-1.root"),
        #("X_0Y_780_Map.root","test20240628-0_rev.root"),
        ("DS189A_944turns_FactoryPillar_FieldIntegration_X200_Y0.root","test20240627-0.root"),
        ("DS189A_944turns_FactoryPillar_FieldIntegration_X500_Y0.root","test20240627-1.root"),
        ("DS189A_944turns_FactoryPillar_FieldIntegration_X650_Y0.root","test20240627-8.root"),
        ("DS189A_944turns_FactoryPillar_FieldIntegration_X760_Y0.root","test20240627-2rev.root"),#rev : the version deleted the last several strange lines 
    ]
    '''
    #Peak Position　Y分布　Laser TrackerPos
    '''
    file_pairs= [
        ("X_0Y_-200_Map.root","test20240710-4.root"),
        ("X_0Y_-729_Map.root","test20240628-2.root"),
        ("X_0Y_0_Map.root","test20240626-6.root"),
        ("X_0Y_200_Map.root","test20240710-1.root"),
        ("X_0Y_780_Map.root","test20240628-0_rev.root"),
    ]
    '''
    
    for operafilename, measfilename in file_pairs:
        OPERAfile = os.path.join(operafile_directory, operafilename)
        MEASfile = os.path.join(measfile_directory, measfilename)
        DeltaofAngle(OPERAfile, MEASfile)
        #factorypillar(OPERAfile, MEASfile)
    #peak position 取得
    '''
    Xaxis = np.array(Xaxlist, dtype = np.float64)
    #Yaxis = np.array(Yaxlist, dtype = np.float64)
    p1meas = np.array(p1meas_list, dtype = np.float64)
    p1measerr = np.array(p1measerr_list, dtype = np.float64)
    p1opera = np.array(p1opera_list, dtype = np.float64)
    p1operaerr = np.array(p1operaerr_list, dtype = np.float64)
    p1corr = np.array(p1corr_list, dtype = np.float64)
    p1correrr = np.array(p1correrr_list, dtype = np.float64)
    p1corrall = np.array(p1corrall_list, dtype = np.float64)
    p1corrallerr = np.array(p1corrallerr_list, dtype = np.float64)
    cpp = ROOT.TCanvas("cpp", "Peak Position", 800, 600)
    gppmeas = ROOT.TGraphErrors(len(Xaxis), Xaxis, p1meas, np.zeros(len(Xaxis)), p1measerr)
    gppmeas.SetMarkerColor(ROOT.kBlue)
    gppmeas.SetMarkerStyle(8)
    gppmeas.SetLineWidth(2)
    gppmeas.SetLineColor(ROOT.kBlue)
    
    #linearmeas = ROOT.TF1("linearmeas", "[0] * x^3 + [1]", np.min(Xaxis), np.max(Xaxis))
    #linearmeas.SetParameters(1e-08, 0)#初期パラメータはp0は0にすべきではないなぜなら関数系が変わりlocal minimumに引っかかりがち(あとp0+ー気をつける)
    #linearmeas.SetLineColor(ROOT.kBlue)
    #gppmeas.Fit(linearmeas, "R")
    #linearmeas.Draw("same")
    
    gppopera = ROOT.TGraphErrors(len(Xaxis), Xaxis, p1opera, np.zeros(len(Xaxis)), p1operaerr)
    #gppopera = ROOT.TGraphErrors(len(Xaxis), Xaxis, p1opera, np.zeros(len(Xaxis)), np.zeros(len(Xaxis)))
    gppopera.SetMarkerColor(ROOT.kRed+2)
    gppopera.SetMarkerStyle(8)
    gppopera.SetLineWidth(2)
    gppopera.SetLineColor(ROOT.kRed+2)
    gppcorr = ROOT.TGraphErrors(len(Xaxis), Xaxis, p1corr, np.zeros(len(Xaxis)), p1correrr)
    gppcorr.SetMarkerColor(ROOT.kAzure+7)
    gppcorr.SetMarkerStyle(8)
    gppcorr.SetLineWidth(2)
    gppcorr.SetLineColor(ROOT.kAzure+7)
    gppcorrall = ROOT.TGraphErrors(len(Xaxis), Xaxis, p1corrall, np.zeros(len(Xaxis)), p1corrallerr)
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
    mgpp.GetXaxis().SetTitle("X (mm)")#peakposition X分布
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
    cpp_dir = "/Users/shohtatakami/physics/COMETDS/DS189A_944turns_Integration/Allcorr/"
    ppoutFile = ROOT.TFile(f"{cpp_dir}PeakPositionX.root", "RECREATE")
    cpp.Write()
    ppoutFile.Close()
    print(f"Successfuly Saved : {ppoutFile}")
    '''
    #ROOT.gApplication.Run()
    
    
    
    
   
    