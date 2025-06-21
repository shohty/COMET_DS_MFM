#!/usr/bin/env python3
import ROOT
import numpy as np
import os
import math
from tqdm import tqdm
# the configuration of Graph by ROOT
ROOT.gStyle.SetLabelSize(0.05, "X") 
ROOT.gStyle.SetLabelSize(0.05, "Y")
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
Errmat = np.array([               #.      Rail, Bx, By, Bz
            [0.15, 0.1, 0.1, 0.0],#Role.   
            [0.10, 0.0, 0.1, 0.1],#Pitch
            [0.10, 0.1, 0.0, 0.1] #Yaw
        ])
def ErrorMatrix(Error_Origin):
    if Error_Origin == 'Rail':
        return np.array([
            [0.15, 0],
            [0.10, 0],
            [0.10, 0]
        ])
    elif Error_Origin == 'BxSensor':
        return np.array([
            [0, 0.1],
            [0, 0.0],
            [0, 0.0]
        ])
    elif Error_Origin == 'BySensor':
        return np.array([
            [0, 0.0],
            [0, 0.1],
            [0, 0.0]
        ])
    elif Error_Origin == 'BzSensor':
        return np.array([
            [0, 0.0],
            [0, 0.0],
            [0, 0.1]
        ])
    elif Error_Origin == 'Sensors':
        return np.array([
            [0, 0.1],
            [0, 0.1],
            [0, 0.1]
        ])
    else:
        print("Wrong Input !!")
def DeltaofAngle(OPERAfile,MEASfile, Error_Origin, Rotation):
    #it calculates effects on each component
    #3rd argument means rotation plane, 4th argument means a value of error
    
    AngleErrors = ErrorMatrix(Error_Origin)
    print(f"{AngleErrors}")
    #err_str = f"{Error:.3f}".replace('.','p')
    #print(f"{err_str}")
             
    #file0= ROOT.TFile.Open(OPERAfile, "READ")#OPERAシミュレーションファイルの読み込み
    file1 = ROOT.TFile.Open(MEASfile, "READ")#OPERAシミュレーションファイルの読み込み
    axname =  os.path.splitext(os.path.basename(OPERAfile))[0]
    print(f"{axname}")
    #if not file0 or file0.IsZombie():
    #    print("Error : could not open the ROOT File")
    if not file1 or file1.IsZombie():
        print("Error : could not open the ROOT File")
    #treeを取得
    #opetree = file0.Get("tree")
    meastree = file1.Get("tree")
    #OPERAファイルより座標と磁場の情報をリストへ
    '''
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
    '''
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
    #X = np.array(xlist, dtype=np.float64)
    #Y = np.array(ylist, dtype=np.float64)
    #Z = np.array(zlist, dtype=np.float64)
    #Bx = np.array(Bxlist, dtype=np.float64)
    #By = np.array(Bylist, dtype=np.float64)
    #Bz = np.array(Bzlist, dtype=np.float64)
    #B = np.array(Blist, dtype=np.float64)
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
    Babscorr= np.zeros(len(Bzmeas))
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
        nBx = rotM_xy(0) @ rotM_zx(angle_zx[i]) @ rotM_xy(phi_Bx) @ rotM_zx(theta_Bx) @ np.array([[1],[0],[0]])
        nBy = rotM_xy(0) @ rotM_yz(angle_yz[i]) @ rotM_xy(phi_By) @ rotM_yz(theta_By) @ np.array([[0],[1],[0]])
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
        Babscorr[i] = Norm(Bcorr[0,0],Bcorr[1,0],Bcorr[2,0])
        err = 0 #誤差の角度（文字に変えて表示）
        #角度の誤差の磁場への寄与を計算
        if Rotation == "Role": # XY  
            #平面(Byは影響なし)
            nBx_pl = rotM_xy(0 + AngleErrors[0,0]) @ rotM_zx(angle_zx[i]) @ rotM_xy(phi_Bx + AngleErrors[0,1]) @ rotM_zx(theta_Bx) @ np.array([[1],[0],[0]]) #誤差のプラス側(Bxセンサー法線ベクトル)
            nBx_mi = rotM_xy(0 - AngleErrors[0,0]) @ rotM_zx(angle_zx[i]) @ rotM_xy(phi_Bx - AngleErrors[0,1]) @ rotM_zx(theta_Bx) @ np.array([[1],[0],[0]]) #誤差のマイナス側
            nBy_pl = rotM_xy(0 + AngleErrors[0,0]) @ rotM_yz(angle_yz[i]) @ rotM_xy(phi_By + AngleErrors[1,1]) @ rotM_yz(theta_By) @ np.array([[0],[1],[0]])
            nBy_mi = rotM_xy(0 - AngleErrors[0,0]) @ rotM_yz(angle_yz[i]) @ rotM_xy(phi_By - AngleErrors[1,1]) @ rotM_yz(theta_By) @ np.array([[0],[1],[0]])
            rotM_pl = np.hstack((nBx_pl, nBy_pl, nBz)).T
            invM_pl = np.linalg.inv(rotM_pl)
            Bcorr_pl = invM_pl @ B_vec
            rotM_mi = np.hstack((nBx_mi, nBy_mi, nBz)).T
            invM_mi = np.linalg.inv(rotM_mi)
            Bcorr_mi = invM_mi @ B_vec
            Bxcorr_pl[i] = Bcorr_pl[0,0]
            Bycorr_pl[i] = Bcorr_pl[1,0]
            Bzcorr_pl[i] = Bcorr_pl[2,0]
            Babscorr_pl[i] = Norm(Bcorr_pl[0,0],Bcorr_pl[1,0],Bcorr_pl[2,0])
            Bxcorr_mi[i] = Bcorr_mi[0,0]
            Bycorr_mi[i] = Bcorr_mi[1,0]
            Bzcorr_mi[i] = Bcorr_mi[2,0]
            Babscorr_mi[i] = Norm(Bcorr_mi[0,0],Bcorr_mi[1,0],Bcorr_mi[2,0])
        elif Rotation == "Pitch": #YZ
            nBy_pl = rotM_xy(0) @ rotM_yz(angle_yz[i] + AngleErrors[1,0]) @ rotM_xy(phi_By) @ rotM_yz(theta_By + AngleErrors[1,1]) @ np.array([[0],[1],[0]])
            nBy_mi = rotM_xy(0) @ rotM_yz(angle_yz[i] - AngleErrors[1,0]) @ rotM_xy(phi_By) @ rotM_yz(theta_By - AngleErrors[1,1]) @ np.array([[0],[1],[0]])
            nBz_pl = rotM_yz(angle_yz[i] + AngleErrors[1,0]) @ rotM_zx(angle_zx[i]) @ rotM_yz(phi_Bz + AngleErrors[2,1]) @ rotM_zx(theta_Bz) @ np.array([[0],[0],[1]])#誤差のプラス側(Bzセンサー法線ベクトル)
            nBz_mi = rotM_yz(angle_yz[i] - AngleErrors[1,0]) @ rotM_zx(angle_zx[i]) @ rotM_yz(phi_Bz - AngleErrors[2,1]) @ rotM_zx(theta_Bz) @ np.array([[0],[0],[1]])#誤差のマイナス側
            rotM_pl = np.hstack((nBx, nBy_pl, nBz_pl)).T
            invM_pl = np.linalg.inv(rotM_pl)
            Bcorr_pl = invM_pl @ B_vec
            rotM_mi = np.hstack((nBx, nBy_mi, nBz_mi)).T
            invM_mi = np.linalg.inv(rotM_mi)
            Bcorr_mi = invM_mi @ B_vec
            Bxcorr_pl[i] = Bcorr_pl[0,0]
            Bycorr_pl[i] = Bcorr_pl[1,0]
            Bzcorr_pl[i] = Bcorr_pl[2,0]
            Babscorr_pl[i] = Norm(Bcorr_pl[0,0],Bcorr_pl[1,0],Bcorr_pl[2,0])
            Bxcorr_mi[i] = Bcorr_mi[0,0]
            Bycorr_mi[i] = Bcorr_mi[1,0]
            Bzcorr_mi[i] = Bcorr_mi[2,0]
            Babscorr_mi[i] = Norm(Bcorr_mi[0,0],Bcorr_mi[1,0],Bcorr_mi[2,0])
        elif Rotation == "Yaw": #ZX
            nBx_pl = rotM_xy(0) @ rotM_zx(angle_zx[i] + AngleErrors[2,0]) @ rotM_xy(phi_Bx) @ rotM_zx(theta_Bx + AngleErrors[0,1]) @ np.array([[1],[0],[0]]) #誤差のプラス側(Bxセンサー法線ベクトル)
            nBx_mi = rotM_xy(0) @ rotM_zx(angle_zx[i] - AngleErrors[2,0]) @ rotM_xy(phi_Bx) @ rotM_zx(theta_Bx - AngleErrors[0,1]) @ np.array([[1],[0],[0]]) #誤差のマイナス側
            nBz_pl = rotM_yz(angle_yz[i]) @ rotM_zx(angle_zx[i] + AngleErrors[2,0]) @ rotM_yz(phi_Bz ) @ rotM_zx(theta_Bz + AngleErrors[2,1]) @ np.array([[0],[0],[1]])#誤差のプラス側(Bzセンサー法線ベクトル)
            nBz_mi = rotM_yz(angle_yz[i]) @ rotM_zx(angle_zx[i] - AngleErrors[2,0]) @ rotM_yz(phi_Bz ) @ rotM_zx(theta_Bz - AngleErrors[2,1]) @ np.array([[0],[0],[1]])#誤差のマイナス側
            rotM_pl = np.hstack((nBx_pl, nBy, nBz_pl)).T
            invM_pl = np.linalg.inv(rotM_pl)
            Bcorr_pl = invM_pl @ B_vec
            rotM_mi = np.hstack((nBx_mi, nBy, nBz_mi)).T
            invM_mi = np.linalg.inv(rotM_mi)
            Bcorr_mi = invM_mi @ B_vec
            Bxcorr_pl[i] = Bcorr_pl[0,0]
            Bycorr_pl[i] = Bcorr_pl[1,0]
            Bzcorr_pl[i] = Bcorr_pl[2,0]
            Babscorr_pl[i] = Norm(Bcorr_pl[0,0],Bcorr_pl[1,0],Bcorr_pl[2,0])
            Bxcorr_mi[i] = Bcorr_mi[0,0]
            Bycorr_mi[i] = Bcorr_mi[1,0]
            Bzcorr_mi[i] = Bcorr_mi[2,0]
            Babscorr_mi[i] = Norm(Bcorr_mi[0,0],Bcorr_mi[1,0],Bcorr_mi[2,0])
        else:
            "seems something is wrong !"
    canvas = ROOT.TCanvas(f"{Rotation}AngleError_{axname}", f"{Rotation}AngleError_{axname}",1200,1200)
    cdelta = ROOT.TCanvas(f"{Rotation}AngleErrorDelta_{axname}", f"{Rotation}AngleErrorDelta{axname}",1200,1200)
    canvas.Divide(2,2)
    cdelta.Divide(2,2)
    canvas.SetMargin(0.1,0.1,0.1,0.1)
    components = {
        'x':{'Title':'Bx', 'base': Bxcorr, 'plus':Bxcorr_pl, 'minus':Bxcorr_mi},
        'y':{'Title':'By','base': Bycorr, 'plus':Bycorr_pl, 'minus':Bycorr_mi},
        'z':{'Title':'Bz','base': Bzcorr, 'plus':Bzcorr_pl, 'minus':Bzcorr_mi},
        'abs':{'Title':'B','base':Babscorr, 'plus':Babscorr_pl, 'minus':Babscorr_mi},
    }
    multigraphs = {}
    delta_multigraphs = {}
    legends = {}
    delta_legends = {}
    for comp, data in components.items():
        g_pl = ROOT.TGraph(len(Zmeas), Zmeas, data['plus'])
        g_pl.SetMarkerStyle(8)
        g_pl.SetMarkerColor(ROOT.kRed)
        g = ROOT.TGraph(len(Zmeas), Zmeas, data['base'])
        g.SetMarkerStyle(8)
        g.SetMarkerColor(ROOT.kSpring)
        g_mi = ROOT.TGraph(len(Zmeas), Zmeas, data['minus'])
        g_mi.SetMarkerStyle(8)
        g_mi.SetMarkerColor(ROOT.kBlue)
        g_pl_delta = ROOT.TGraph(len(Zmeas), Zmeas, data['plus'] - data['base'])
        g_pl_delta.SetMarkerStyle(8)
        g_pl_delta.SetMarkerColor(ROOT.kRed)
        g_mi_delta = ROOT.TGraph(len(Zmeas), Zmeas, data['minus'] - data['base'])
        g_mi_delta.SetMarkerStyle(8)
        g_mi_delta.SetMarkerColor(ROOT.kBlue)
        mg = ROOT.TMultiGraph()
        mg.Add(g_pl)
        mg.Add(g)
        mg.Add(g_mi)
        mg.SetTitle(f"{axname} : {data['Title']} {Rotation}")#
        mg.GetXaxis().SetTitle("Z (mm)")
        mg.GetYaxis().SetTitle(f"{data["Title"]} (T)")
        mglegend = ROOT.TLegend(0.7, 0.1, 0.9, 0.25)
        mglegend.AddEntry(g_pl, "#theta + #delta#theta", "p")
        mglegend.AddEntry(g, "#theta", "p")
        mglegend.AddEntry(g_mi, "#theta - #delta#theta", "p")
        legends[data['Title']] = mglegend
        mg_delta = ROOT.TMultiGraph()
        mg_delta.Add(g_pl_delta)
        mg_delta.Add(g_mi_delta)
        mg_delta.SetTitle(f"{axname} : #Delta{data['Title']} {Rotation}")#
        mg_delta.GetXaxis().SetTitle("Z (mm)")
        mg_delta.GetYaxis().SetTitle(f"#Delta{data["Title"]} (T)")
        mg_deltalegend = ROOT.TLegend(0.7, 0.1, 0.9, 0.25)
        mg_deltalegend.AddEntry(g_pl_delta, "#theta + #delta#theta", "p")
        mg_deltalegend.AddEntry(g_mi_delta, "#theta - #delta#theta", "p")
        delta_legends[data['Title']] = mg_deltalegend
        #TGraphをdictで保管
        multigraphs[f'{data['Title']}'] = mg
        delta_multigraphs[f'{data['Title']}'] = mg_delta
    
    for i, label in enumerate(['Bx', 'By', 'Bz', 'B'], start =1):
        canvas.cd(i)
        multigraphs[f'{label}'].Draw('AP')
        legends[f'{label}'].Draw()
        ROOT.gPad.SetGrid(1,1)
        canvas.Update()
        cdelta.cd(i)
        delta_multigraphs[f'{label}'].Draw('AP')
        delta_legends[f'{label}'].Draw()
        ROOT.gPad.SetGrid(1,1)
        cdelta.Update()
    canvas_dir = "/Users/shohtatakami/physics/COMETDS/ErrorBudget/" + f"{Error_Origin}/" + f"{Rotation}/"#保存ディレクトリは平面に応じて分ける
    os.makedirs(canvas_dir, exist_ok = True)#ディレクトリがなければ作成
    canvasFile = ROOT.TFile(f"{canvas_dir}{axname}.root", "RECREATE")
    canvas.Write()
    canvasFile.Close()
    print(f"Successfully saved : {canvasFile}") 
    cdelta_dir = "/Users/shohtatakami/physics/COMETDS/ErrorBudget/" + f"{Error_Origin}/" + f"{Rotation}_Delta/"#保存ディレクトリは平面に応じて分ける
    os.makedirs(cdelta_dir, exist_ok = True)#ディレクトリがなければ作成
    cdeltaFile = ROOT.TFile(f"{cdelta_dir}{axname}_delta.root", "RECREATE")
    cdelta.Write()
    cdeltaFile.Close()
    print(f"Successfully saved : {cdeltaFile}") 
    
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
        ("X_760Y_0_revMap.root","test20240627-2rev.root"),#rev : the version deleted the last several strange lines 
    ]
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
    
    for operafilename, measfilename in tqdm (file_pairs):
        OPERAfile = os.path.join(operafile_directory, operafilename)
        MEASfile = os.path.join(measfile_directory, measfilename)
        DeltaofAngle(OPERAfile, MEASfile, "Rail", "Role")
        DeltaofAngle(OPERAfile, MEASfile, "Rail", "Pitch")
        DeltaofAngle(OPERAfile, MEASfile, "Rail", "Yaw")
        DeltaofAngle(OPERAfile, MEASfile, "BxSensor", "Role")
        DeltaofAngle(OPERAfile, MEASfile, "BxSensor", "Yaw")
        DeltaofAngle(OPERAfile, MEASfile, "BySensor", "Role")
        DeltaofAngle(OPERAfile, MEASfile, "BySensor", "Pitch")
        DeltaofAngle(OPERAfile, MEASfile,"BzSensor", "Pitch")
        DeltaofAngle(OPERAfile, MEASfile,"BzSensor", "Yaw")
        DeltaofAngle(OPERAfile, MEASfile,"Sensors", "Role")
        DeltaofAngle(OPERAfile, MEASfile,"Sensors", "Pitch")
        DeltaofAngle(OPERAfile, MEASfile,"Sensors", "Yaw")
        #センサーの誤差はセンサー個別に効くはずー＞せいぶんごとにやらないといけない（自動化したいが今は手動で消す）

        
   