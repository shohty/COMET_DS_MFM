#!/usr/bin/env python3
import ROOT
import numpy as np
import os
import math
# the configuration of Graph by ROOT
ROOT.gStyle.SetLabelSize(0.06, "X") 
ROOT.gStyle.SetLabelSize(0.06, "Y")
ROOT.gStyle.SetTitleSize(0.05, "X")
ROOT.gStyle.SetTitleSize(0.05, "Y")
ROOT.gStyle.SetLegendTextSize(0.03)
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
def Norm(x, y, z):
    return np.sqrt(x ** 2 +y ** 2 + z ** 2)

def DeltaofAngle(OPERAfile,MEASfile, Rotation, RailErrors, SensorErrors):
    #it calculates effects on each component
    #3rd argument means rotation plane, 4th argument means a value of error
    err_str = f"{Error:.3f}".replace('.','p')
    print(f"{err_str}")
             
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
    Bxcorr_pl = np.zeros(len(Bxmeas))
    Bxcorr_mi = np.zeros(len(Bxmeas))
    Bycorr_pl = np.zeros(len(Bxmeas))
    Bycorr_mi = np.zeros(len(Bxmeas))
    Bzcorr_pl = np.zeros(len(Bxmeas))
    Bzcorr_mi = np.zeros(len(Bxmeas))
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
        if Rotation == "Role": # XY    
            #平面(Byは影響なし)
            nBx_pl = rotM_zx(angle_zx[i] + RailErrors[0]) @ rotM_xy(phi_Bx) @ rotM_zx(theta_Bx) @ np.array([[1],[0],[0]]) #誤差のプラス側(Bxセンサー法線ベクトル)
            nBy_pl = rotM_yz(angle_yz[i] + RailErrors[0]) @ rotM_xy(phi_By) @ rotM_yz(theta_By) @ np.array([[0],[1],[0]])
            nBx_mi = rotM_zx(angle_zx[i] - RailErrors[0]) @ rotM_xy(phi_Bx) @ rotM_zx(theta_Bx) @ np.array([[1],[0],[0]]) #誤差のマイナス側
            nBy_mi = rotM_yz(angle_yz[i] - RailErrors[0]) @ rotM_xy(phi_By) @ rotM_yz(theta_By) @ np.array([[0],[1],[0]])
            rotM_pl = np.hstack((nBx_pl, nBy_pl, nBz)).T
            invM_pl = np.linalg.inv(rotM_pl)
            Bcorr_pl = invM_pl @ B_vec
            rotM_mi = np.hstack((nBx_mi, nBy_mi, nBz)).T
            invM_mi = np.linalg.inv(rotM_mi)
            Bcorr_mi = invM_mi @ B_vec
            Bxcorr_pl[i] = Bcorr_pl[0,0]
            Bycorr_pl[i] = Bcorr_pl[1,0]
            Bzcorr_pl[i] = Bcorr_pl[2,0]
            Bxcorr_mi[i] = Bcorr_mi[0,0]
            Bycorr_mi[i] = Bcorr_mi[1,0]
            Bzcorr_mi[i] = Bcorr_mi[2,0]
        elif Rotation == "Pitch": #YZ
            nBy_pl = rotM_yz(angle_yz[i] + RailErrors[1]) @ rotM_xy(phi_By) @ rotM_yz(theta_By) @ np.array([[0],[1],[0]])
            nBy_mi = rotM_yz(angle_yz[i] + RailErrors[1]) @ rotM_xy(phi_By) @ rotM_yz(theta_By) @ np.array([[0],[1],[0]])
            nBz_pl = rotM_yz(angle_yz[i] - RailErrors[1]) @ rotM_zx(angle_zx[i]) @ rotM_yz(phi_Bz) @ rotM_zx(theta_Bz) @ np.array([[0],[0],[1]])#誤差のプラス側(Bzセンサー法線ベクトル)
            nBz_mi = rotM_yz(angle_yz[i] - RailErrors[1]) @ rotM_zx(angle_zx[i]) @ rotM_yz(phi_Bz) @ rotM_zx(theta_Bz) @ np.array([[0],[0],[1]])#誤差のマイナス側
            rotM_pl = np.hstack((nBx, nBy_pl, nBz_pl)).T
            invM_pl = np.linalg.inv(rotM_pl)
            Bcorr_pl = invM_pl @ B_vec
            rotM_mi = np.hstack((nBx, nBy_mi, nBz_mi)).T
            invM_mi = np.linalg.inv(rotM_mi)
            Bcorr_mi = invM_mi @ B_vec
            Bxcorr_pl[i] = Bcorr_pl[0,0]
            Bycorr_pl[i] = Bcorr_pl[1,0]
            Bzcorr_pl[i] = Bcorr_pl[2,0]
            Bxcorr_mi[i] = Bcorr_mi[0,0]
            Bycorr_mi[i] = Bcorr_mi[1,0]
            Bzcorr_mi[i] = Bcorr_mi[2,0]
        elif Rotation == "Yaw": #ZX
            nBx_pl = rotM_zx(angle_zx[i] + RailErrors[2]) @ rotM_xy(phi_Bx) @ rotM_zx(theta_Bx) @ np.array([[1],[0],[0]]) #誤差のプラス側(Bxセンサー法線ベクトル)
            nBx_mi = rotM_zx(angle_zx[i] + RailErrors[2]) @ rotM_xy(phi_Bx) @ rotM_zx(theta_Bx) @ np.array([[1],[0],[0]]) #誤差のマイナス側
            nBz_pl = rotM_yz(angle_yz[i] - RailErrors[2]) @ rotM_zx(angle_zx[i]) @ rotM_yz(phi_Bz ) @ rotM_zx(theta_Bz) @ np.array([[0],[0],[1]])#誤差のプラス側(Bzセンサー法線ベクトル)
            nBz_mi = rotM_yz(angle_yz[i] - RailErrors[2]) @ rotM_zx(angle_zx[i]) @ rotM_yz(phi_Bz ) @ rotM_zx(theta_Bz) @ np.array([[0],[0],[1]])#誤差のマイナス側
            rotM_pl = np.hstack((nBx_pl, nBy, nBz_pl)).T
            invM_pl = np.linalg.inv(rotM_pl)
            Bcorr_pl = invM_pl @ B_vec
            rotM_mi = np.hstack((nBx_mi, nBy, nBz_mi)).T
            invM_mi = np.linalg.inv(rotM_mi)
            Bcorr_mi = invM_mi @ B_vec
            Bxcorr_pl[i] = Bcorr_pl[0,0]
            Bycorr_pl[i] = Bcorr_pl[1,0]
            Bzcorr_pl[i] = Bcorr_pl[2,0]
            Bxcorr_mi[i] = Bcorr_mi[0,0]
            Bycorr_mi[i] = Bcorr_mi[1,0]
            Bzcorr_mi[i] = Bcorr_mi[2,0]    
        else:
            "seems something is wrong !"
    def DrawGraphProtocol(Plane):
    #general rule : right handed system. First components corresponds to first letter of plane name. Second is same as First.
        if Rotation == "XY":
            return{
                "DIR": f"delta_yaw/",
                "ID": f"YAW_",
                "Title1": "Bx", "B1corr": Bxcorr, "B1corr_pl": Bxcorr_pl, "B1corr_mi": Bxcorr_mi,
                "Title2": "By", "B2corr": Bycorr, "B2corr_pl": Bycorr_pl, "B2corr_mi": Bycorr_mi
            }
        elif Rotation == "YZ":
            return{
                "DIR": f"delta_role/",
                "ID": f"ROLE_",
                "Title1": "By", "B1corr": Bycorr, "B1corr_pl": Bycorr_pl, "B1corr_mi": Bycorr_mi,
                "Title2": "Bz", "B2corr": Bzcorr, "B2corr_pl": Bzcorr_pl, "B2corr_mi": Bzcorr_mi
            }
        elif Rotation == "ZX":
            return{
                "DIR": f"delta_pitch/",
                "ID": f"PITCH_",
                "Title1": "Bz", "B1corr": Bzcorr, "B1corr_pl": Bzcorr_pl, "B1corr_mi": Bzcorr_mi,
                "Title2": "Bx", "B2corr": Bxcorr, "B2corr_pl": Bxcorr_pl, "B2corr_mi": Bxcorr_mi
            }
        else:
            raise print("Wrong Input !!")              

    canvas = ROOT.TCanvas(f"{Rotation}AngleError_{axname}", f"{Rotation}Angle_Error{axname}",1200,900)
    canvas.Divide(2,2)
    config = DrawGraphProtocol(Rotation)
        
    #1はこの場合x,2はy
    gB1_pl = ROOT.TGraph(len(Z), Z, config["B1corr_pl"])
    gB1_pl.SetMarkerStyle(8)
    gB1_pl.SetMarkerColor(ROOT.kRed)
    gB1 = ROOT.TGraph(len(Z),Z, config["B1corr"])
    gB1.SetMarkerStyle(8)
    gB1.SetMarkerColor(ROOT.kSpring)
    gB1_mi = ROOT.TGraph(len(Z), Z, config["B1corr_mi"])
    gB1_mi.SetMarkerStyle(8)
    gB1_mi.SetMarkerColor(ROOT.kBlue)
    
    gB2_pl = ROOT.TGraph(len(Z), Z, config["B2corr_pl"])
    gB2_pl.SetMarkerStyle(8)
    gB2_pl.SetMarkerColor(ROOT.kRed)
    gB2 = ROOT.TGraph(len(Z),Z, config["B2corr"])
    gB2.SetMarkerStyle(8)
    gB2.SetMarkerColor(ROOT.kSpring)
    gB2_mi = ROOT.TGraph(len(Z), Z, config["B2corr_mi"])
    gB2_mi.SetMarkerStyle(8)
    gB2_mi.SetMarkerColor(ROOT.kBlue)
    
    gB1delta_pl = ROOT.TGraph(len(Z), Z, config["B1corr_pl"] - config["B1corr"])
    gB1delta_pl.SetMarkerStyle(8)
    gB1delta_pl.SetMarkerColor(ROOT.kRed)
    gB1delta_mi = ROOT.TGraph(len(Z), Z, config["B1corr_mi"] - config["B1corr"])
    gB1delta_mi.SetMarkerStyle(8)
    gB1delta_mi.SetMarkerColor(ROOT.kBlue)
    
    gB2delta_pl = ROOT.TGraph(len(Z), Z, config["B2corr_pl"] - config["B2corr"])
    gB2delta_pl.SetMarkerStyle(8)
    gB2delta_pl.SetMarkerColor(ROOT.kRed)
    gB2delta_mi = ROOT.TGraph(len(Z), Z, config["B2corr_mi"] - config["B2corr"])
    gB2delta_mi.SetMarkerStyle(8)
    gB2delta_mi.SetMarkerColor(ROOT.kBlue)
    canvas.cd(1)
    #\pm delta \thetaと元の補正値
    mg_B1 = ROOT.TMultiGraph()
    mg_B1.Add(gB1_pl)
    mg_B1.Add(gB1_mi)
    mg_B1.Add(gB1)
    mg_B1.SetTitle(f"{axname} : {config["Title1"]} {config["ID"]}_{err_str}")#
    mg_B1.GetXaxis().SetTitle("Z (mm)")
    mg_B1.GetYaxis().SetTitle(f"{config["Title1"]} (T)")
    mg_B1.Draw("AP")
    B1_legend = ROOT.TLegend(0.7, 0.1, 0.9, 0.25)
    B1_legend.AddEntry(gB1_pl, "#theta + #delta#theta", "p")
    B1_legend.AddEntry(gB1, "#theta ", "p")
    B1_legend.AddEntry(gB1_mi, "#theta - #delta#theta", "p")
    B1_legend.SetFillStyle(0)
    B1_legend.Draw()
    ROOT.gPad.SetGrid(1,1)
    canvas.cd(2)
    mg_B2 = ROOT.TMultiGraph()
    mg_B2.Add(gB2_pl)
    mg_B2.Add(gB2_mi)
    mg_B2.Add(gB2)
    mg_B2.SetTitle(f"{axname} : {config["Title2"]} {config["ID"]}_{err_str}")#
    mg_B2.GetXaxis().SetTitle("Z (mm)")
    mg_B2.GetYaxis().SetTitle(f"{config["Title2"]}(T)")
    mg_B2.Draw("AP")
    B2_legend = ROOT.TLegend(0.7, 0.1, 0.9, 0.25)
    B2_legend.AddEntry(gB2_pl, "#theta + #delta#theta", "p")
    B2_legend.AddEntry(gB2, "#theta ", "p")
    B2_legend.AddEntry(gB2_mi, "#theta - #delta#theta", "p")
    B2_legend.SetFillStyle(0)
    B2_legend.Draw()
    ROOT.gPad.SetGrid(1,1)
    canvas.cd(3)
    mg_B1delta = ROOT.TMultiGraph()
    mg_B1delta.Add(gB1delta_pl)
    mg_B1delta.Add(gB1delta_mi)
    mg_B1delta.SetTitle(f"{axname} : #delta {config["Title1"]} {config["ID"]}_{err_str}")#
    mg_B1delta.GetXaxis().SetTitle("Z (mm)")
    mg_B1delta.GetYaxis().SetTitle(f"#delta {config["Title1"]}(T)")
    mg_B1delta.Draw("AP")
    B1delta_legend = ROOT.TLegend(0.7, 0.1, 0.9, 0.25)
    B1delta_legend.AddEntry(gB1delta_pl, "#theta + #delta#theta", "p")
    B1delta_legend.AddEntry(gB1delta_mi, "#theta - #delta#theta", "p")
    B1delta_legend.SetFillStyle(0)
    B1delta_legend.Draw()
    ROOT.gPad.SetGrid(1,1)
    canvas.cd(4)
    mg_B2delta = ROOT.TMultiGraph()
    mg_B2delta.Add(gB2delta_pl)
    mg_B2delta.Add(gB2delta_mi)
    mg_B2delta.SetTitle(f"{axname} : #delta {config["Title2"]} {config["ID"]}_{err_str}")#
    mg_B2delta.GetXaxis().SetTitle("Z (mm)")
    mg_B2delta.GetYaxis().SetTitle(f"#delta {config["Title2"]}(T)")
    mg_B2delta.Draw("AP")
    B2delta_legend = ROOT.TLegend(0.7, 0.1, 0.9, 0.25)
    B2delta_legend.AddEntry(gB2delta_pl, "#theta + #delta#theta", "p")
    B2delta_legend.AddEntry(gB2delta_mi, "#theta - #delta#theta", "p")
    B2delta_legend.SetFillStyle(0)
    B2delta_legend.Draw()
    ROOT.gPad.SetGrid(1,1)
    canvas.Update()
    canvas_dir = "/Users/shohtatakami/physics/COMETDS/ErrorBudget/" + config["DIR"] + f"{err_str}/"#保存ディレクトリは平面に応じて分ける
    os.makedirs(canvas_dir, exist_ok = True)#ディレクトリがなければ作成
    delta_File = ROOT.TFile(f"{canvas_dir}{axname}.root", "RECREATE")
    canvas.Write()
    delta_File.Close()
    print(f"Successfully saved : {delta_File}") 
    
def DeltaofAngleAbs(OPERAfile,MEASfile,Plane, Sensor_err, Rail_err):
    #it calculates effects on |B|
    #3rd argument means rotation plane, 4th argument means a value of error
    err_str = f"{Error:.3f}".replace('.','p')
    print(f"{err_str}")
         
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
    Babscorr= np.zeros(len(Bxmeas))
    Bcorr= np.zeros(len(Bzmeas))
    Bxcorr_pl = np.zeros(len(Bxmeas))
    Bxcorr_mi = np.zeros(len(Bxmeas))
    Bycorr_pl = np.zeros(len(Bxmeas))
    Bycorr_mi = np.zeros(len(Bxmeas))
    Bzcorr_pl = np.zeros(len(Bxmeas))
    Bzcorr_mi = np.zeros(len(Bxmeas))
    Babscorr_pl = np.zeros(len(Bxmeas))
    Babscorr_mi = np.zeros(len(Bxmeas))
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
        Babscorr[i] = np.linalg.norm(Bcorr)
        #角度の誤差の磁場への寄与を計算
        if Plane == "XY":    
            #各センサーの回転
            nBx_pl = rotM_zx(angle_zx[i]) @ rotM_xy(phi_Bx + Error) @ rotM_zx(theta_Bx) @ np.array([[1],[0],[0]]) #誤差のプラス側(Bxセンサー法線ベクトル)
            nBy_pl = rotM_yz(angle_yz[i]) @ rotM_xy(phi_By + Error) @ rotM_yz(theta_By) @ np.array([[0],[1],[0]])
            nBx_mi = rotM_zx(angle_zx[i]) @ rotM_xy(phi_Bx - Error) @ rotM_zx(theta_Bx) @ np.array([[1],[0],[0]]) #誤差のマイナス側
            nBy_mi = rotM_yz(angle_yz[i]) @ rotM_xy(phi_By - Error) @ rotM_yz(theta_By) @ np.array([[0],[1],[0]])
            rotM_pl = np.hstack((nBx_pl, nBy_pl, nBz)).T #1行目からnBx,nBy,nBzとなるようにしている
            invM_pl = np.linalg.inv(rotM_pl)#それの逆行列
            Bcorr_pl = invM_pl @ B_vec#逆行列を測定座標系(各種回転が混じった)のBベクトルから世紀直交座標に補正
            rotM_mi = np.hstack((nBx_mi, nBy_mi, nBz)).T
            invM_mi = np.linalg.inv(rotM_mi)
            Bcorr_mi = invM_mi @ B_vec
            Bxcorr_pl[i] = Bcorr_pl[0,0]
            Bycorr_pl[i] = Bcorr_pl[1,0]
            Bxcorr_mi[i] = Bcorr_mi[0,0]
            Bycorr_mi[i] = Bcorr_mi[1,0]
            Babscorr_pl[i] = np.linalg.norm(Bcorr_pl)
            Babscorr_mi[i] = np.linalg.norm(Bcorr_mi)
        elif Plane == "YZ":
            nBy_pl = rotM_yz(angle_yz[i]) @ rotM_xy(phi_By) @ rotM_yz(theta_By + Error) @ np.array([[0],[1],[0]])
            nBy_mi = rotM_yz(angle_yz[i]) @ rotM_xy(phi_By) @ rotM_yz(theta_By - Error) @ np.array([[0],[1],[0]])
            nBz_pl = rotM_yz(angle_yz[i]) @ rotM_zx(angle_zx[i]) @ rotM_yz(phi_Bz + Error) @ rotM_zx(theta_Bz) @ np.array([[0],[0],[1]])#誤差のプラス側(Bzセンサー法線ベクトル)
            nBz_mi = rotM_yz(angle_yz[i]) @ rotM_zx(angle_zx[i]) @ rotM_yz(phi_Bz - Error) @ rotM_zx(theta_Bz) @ np.array([[0],[0],[1]])#誤差のマイナス側
            rotM_pl = np.hstack((nBx, nBy_pl, nBz_pl)).T
            invM_pl = np.linalg.inv(rotM_pl)
            Bcorr_pl = invM_pl @ B_vec
            rotM_mi = np.hstack((nBx, nBy_mi, nBz_mi)).T
            invM_mi = np.linalg.inv(rotM_mi)
            Bcorr_mi = invM_mi @ B_vec
            Bycorr_pl[i] = Bcorr_pl[1,0]
            Bzcorr_pl[i] = Bcorr_pl[2,0]
            Bycorr_mi[i] = Bcorr_mi[1,0]
            Bzcorr_mi[i] = Bcorr_mi[2,0]
            Babscorr_pl[i] = np.linalg.norm(Bcorr_pl)
            Babscorr_mi[i] = np.linalg.norm(Bcorr_mi)
              
        elif Plane == "ZX":
            nBx_pl = rotM_zx(angle_zx[i]) @ rotM_xy(phi_Bx) @ rotM_zx(theta_Bx + Error) @ np.array([[1],[0],[0]]) #誤差のプラス側(Bxセンサー法線ベクトル)
            nBx_mi = rotM_zx(angle_zx[i]) @ rotM_xy(phi_Bx) @ rotM_zx(theta_Bx - Error) @ np.array([[1],[0],[0]]) #誤差のマイナス側
            nBz_pl = rotM_yz(angle_yz[i]) @ rotM_zx(angle_zx[i]) @ rotM_yz(phi_Bz ) @ rotM_zx(theta_Bz + Error) @ np.array([[0],[0],[1]])#誤差のプラス側(Bzセンサー法線ベクトル)
            nBz_mi = rotM_yz(angle_yz[i]) @ rotM_zx(angle_zx[i]) @ rotM_yz(phi_Bz ) @ rotM_zx(theta_Bz - Error) @ np.array([[0],[0],[1]])#誤差のマイナス側
            rotM_pl = np.hstack((nBx_pl, nBy, nBz_pl)).T
            invM_pl = np.linalg.inv(rotM_pl)
            Bcorr_pl = invM_pl @ B_vec
            rotM_mi = np.hstack((nBx_mi, nBy, nBz_mi)).T
            invM_mi = np.linalg.inv(rotM_mi)
            Bcorr_mi = invM_mi @ B_vec
            Bzcorr_pl[i] = Bcorr_pl[2,0]
            Bxcorr_pl[i] = Bcorr_pl[0,0]
            Bzcorr_mi[i] = Bcorr_mi[2,0]
            Bxcorr_mi[i] = Bcorr_mi[0,0]
            Babscorr_pl[i] = np.linalg.norm(Bcorr_pl)
            Babscorr_mi[i] = np.linalg.norm(Bcorr_mi)
        else:
            "seems something is wrong !"
    def DrawGraphProtocol(Plane):
    #general rule : right handed system. First components corresponds to first letter of plane name. Second is same as First.
        if Plane == "XY":
            return{
                "DIR": f"delta_yaw/",
                "ID": f"YAW_",
                "Title": "Babs", "Babscorr": Babscorr, "Babscorr_pl": Babscorr_pl, "Babscorr_mi": Babscorr_mi,
            }
        elif Plane == "YZ":
            return{
                "DIR": f"delta_role/",
                "ID": f"ROLE_",
                "Title": "Babs", "Babscorr": Babscorr, "Babscorr_pl": Babscorr_pl, "Babscorr_mi": Babscorr_mi,
            }
        elif Plane == "ZX":
            return{
                "DIR": f"delta_pitch/",
                "ID": f"PITCH_",
                "Title": "Babs", "Babscorr": Babscorr, "Babscorr_pl": Babscorr_pl, "Babscorr_mi": Babscorr_mi,
            }
        else:
            raise print("Wrong Input !!")              

    canvas = ROOT.TCanvas(f"{Plane}AngleError_{axname}", f"{Plane}Angle_Error{axname}",1200,900)
    canvas.Divide(1,2)
    config = DrawGraphProtocol(Plane)
        
    #1はこの場合x,2はy
    gBabs_pl = ROOT.TGraph(len(Z), Z, config["Babscorr_pl"])
    gBabs_pl.SetMarkerStyle(8)
    gBabs_pl.SetMarkerColor(ROOT.kRed)
    gBabs = ROOT.TGraph(len(Z),Z, config["Babscorr"])
    gBabs.SetMarkerStyle(8)
    gBabs.SetMarkerColor(ROOT.kSpring)
    gBabs_mi = ROOT.TGraph(len(Z), Z, config["Babscorr_mi"])
    gBabs_mi.SetMarkerStyle(8)
    gBabs_mi.SetMarkerColor(ROOT.kBlue)
    
    gBabsdelta_pl = ROOT.TGraph(len(Z), Z, config["Babscorr_pl"] - config["Babscorr"])
    gBabsdelta_pl.SetMarkerStyle(8)
    gBabsdelta_pl.SetMarkerColor(ROOT.kRed)
    gBabsdelta_mi = ROOT.TGraph(len(Z), Z, config["Babscorr_mi"] - config["Babscorr"])
    gBabsdelta_mi.SetMarkerStyle(8)
    gBabsdelta_mi.SetMarkerColor(ROOT.kBlue)
    canvas.cd(1)
    #\pm delta \thetaと元の補正値
    mg_Babs = ROOT.TMultiGraph()
    mg_Babs.Add(gBabs_pl)
    mg_Babs.Add(gBabs_mi)
    mg_Babs.Add(gBabs)
    mg_Babs.SetTitle(f"{axname} : {config["Title"]} {config["ID"]}_{err_str}")#
    mg_Babs.GetXaxis().SetTitle("Z (mm)")
    mg_Babs.GetYaxis().SetTitle(f"{config["Title"]} (T)")
    mg_Babs.Draw("AP")
    Babs_legend = ROOT.TLegend(0.7, 0.1, 0.9, 0.25)
    Babs_legend.AddEntry(gBabs_pl, "#theta + #delta#theta", "p")
    Babs_legend.AddEntry(gBabs, "#theta ", "p")
    Babs_legend.AddEntry(gBabs_mi, "#theta - #delta#theta", "p")
    Babs_legend.SetFillStyle(0)
    Babs_legend.Draw()
    ROOT.gPad.SetGrid(1,1)
    canvas.cd(2)
    mg_Babsdelta = ROOT.TMultiGraph()
    mg_Babsdelta.Add(gBabsdelta_pl)
    mg_Babsdelta.Add(gBabsdelta_mi)
    mg_Babsdelta.SetTitle(f"{axname} : #delta {config["Title"]} {config["ID"]}_{err_str}")#
    mg_Babsdelta.GetXaxis().SetTitle("Z (mm)")
    mg_Babsdelta.GetYaxis().SetTitle(f"#delta {config["Title"]}(T)")
    mg_Babsdelta.Draw("AP")
    Babsdelta_legend = ROOT.TLegend(0.7, 0.1, 0.9, 0.25)
    Babsdelta_legend.AddEntry(gBabsdelta_pl, "#theta + #delta#theta", "p")
    Babsdelta_legend.AddEntry(gBabsdelta_mi, "#theta - #delta#theta", "p")
    Babsdelta_legend.SetFillStyle(0)
    Babsdelta_legend.Draw()
    ROOT.gPad.SetGrid(1,1)
    canvas.Update()
    canvas_dir = "/Users/shohtatakami/physics/COMETDS/ErrorBudget/" + config["DIR"] + f"{err_str}/abs/"#保存ディレクトリは平面に応じて分ける
    os.makedirs(canvas_dir, exist_ok = True)#ディレクトリがなければ作成
    delta_File = ROOT.TFile(f"{canvas_dir}{axname}.root", "RECREATE")
    canvas.Write()
    delta_File.Close()
    print(f"Successfully saved : {delta_File}") 
          
if __name__=='__main__':
    SensorErrors = np.array([0.1], #Error of Sensors
                       [0]) #Decoy
    RailErrors_ = np.array([[0.15], #Error of Rail_Role
                       [0.02], #Error of Rail_Pitch
                       [0.006], #Error of Rail_Yaw
                       [0]] ) #Decoy
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
        ("X_760Y_0_revMap.root","test20240627-2rev.root"),#rev : the version deleted the last several strange lines 
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
        DeltaofAngleAbs(OPERAfile, MEASfile,"XY",0.14)
        DeltaofAngleAbs(OPERAfile, MEASfile,"YZ",0.14)
        DeltaofAngleAbs(OPERAfile, MEASfile,"ZX",0.17)
        DeltaofAngleAbs(OPERAfile, MEASfile,"XY",0.1)
        DeltaofAngleAbs(OPERAfile, MEASfile,"YZ",0.1)
        DeltaofAngleAbs(OPERAfile, MEASfile,"ZX",0.1)
        DeltaofAngleAbs(OPERAfile, MEASfile,"XY",0.15)
        DeltaofAngleAbs(OPERAfile, MEASfile,"YZ",0.02)
        DeltaofAngleAbs(OPERAfile, MEASfile,"ZX",0.006)
        
   