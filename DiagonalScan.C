#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <TCanvas.h>
#include <TF1.h>        
#include <TGraphErrors.h>
#include <TGraph2D.h>
#include <TMultiGraph.h>
#include <TPolyLine3D.h>
#include <TApplication.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TGaxis.h>
#include <TFile.h>
#include <TTree.h>

using namespace std;
///
double linear_Interpolation(double left, double left_B, double right, double right_B, double x){
    //線形補間の関数
    // leftはdown,bottomに対応
    // rightはup,topに対応
    double l = right -left;
    double B = ((right - x)/l) * left_B + ((x - left)/l) * right_B;//内分計算
    if(right -x > 10){
        cout<<"right : WRONG !!"<<endl;
    }
    if(x - left > 10){
        cout<<"left : WRONG !!"<<endl;
    }
    return B;
}

struct weighted_ave_Param{
    double p0_x,p0_y,p0_z,p1_x,p1_y,p1_z,p2_x,p2_y,p2_z,p3_x,p3_y,p3_z;
    double p4_x,p4_y,p4_z,p5_x,p5_y,p5_z,p6_x,p6_y,p6_z,p7_x,p7_y,p7_z;
    double p0_B,p1_B,p2_B,p3_B,p4_B,p5_B,p6_B,p7_B;
    double p0_Bx,p1_Bx,p2_Bx,p3_Bx,p4_Bx,p5_Bx,p6_Bx,p7_Bx;
    double p0_By,p1_By,p2_By,p3_By,p4_By,p5_By,p6_By,p7_By;
    double p0_Bz,p1_Bz,p2_Bz,p3_Bz,p4_Bz,p5_Bz,p6_Bz,p7_Bz;
    double x, y, z; //任意の点の座標
};

tuple<double, double, double, double> weighted_average(const weighted_ave_Param& params){
    /*_____________________________________________________________________________
    任意の点(x,y,z)が所属する三次元の格子を考える
    DSをDAQシステム側、z軸向かって左側から見る
    (z,x,y)の右手系
    y軸を垂直上下方向とすると上から見たときに底面左下p0,右下p1,左上p2,右上p3
     　　　　　　　　　　　　　　　　　　　　　　上面左下p4,右下p5,左上p6,右上p7とする
    z正方向 right : x正方向 up : y正方向 top
    z負方向 left : x負方向 down : y負方向 bottom
    
    box_numは各軸方向の格子の番号
    それぞれ原点から正方向に0,1,2,3,,,,
                 負方向に1,2,3,,,,
    _______________________________________________________________________________*/
    
    //底面の下辺上辺の磁場を計算し、底面のz,xにおける磁場を計算
    double B_bottomdown = linear_Interpolation(params.p0_z, params.p0_B, params.p1_z, params.p1_B, params.z);
    double B_bottomup = linear_Interpolation(params.p2_z, params.p2_B, params.p3_z, params.p3_B, params.z);
    double B_bottom = linear_Interpolation(params.p0_x, B_bottomdown, params.p2_x, B_bottomup, params.x);
    //上面の下辺上辺の磁場を計算し、上面のz,xにおける磁場を計算
    double B_topdown = linear_Interpolation(params.p4_z, params.p4_B, params.p5_z, params.p5_B, params.z);
    double B_topup = linear_Interpolation(params.p6_z, params.p6_B, params.p7_z, params.p7_B, params.z);
    double B_top = linear_Interpolation(params.p4_x, B_topdown, params.p6_x, B_topup, params.x);
    //y方向の補間
    double B = linear_Interpolation(params.p0_y, B_bottom, params.p4_y, B_top, params.y);
    //Bxも同様
    //底面の下辺上辺の磁場を計算し、底面のz,xにおける磁場を計算
    double Bx_bottomdown = linear_Interpolation(params.p0_z, params.p0_Bx, params.p1_z, params.p1_Bx, params.z);
    double Bx_bottomup = linear_Interpolation(params.p2_z, params.p2_Bx, params.p3_z, params.p3_Bx, params.z);
    double Bx_bottom = linear_Interpolation(params.p0_x, Bx_bottomdown, params.p2_x, Bx_bottomup, params.x);
    //上面の下辺上辺の磁場を計算し、上面のz,xにおける磁場を計算
    double Bx_topdown = linear_Interpolation(params.p4_z, params.p4_Bx, params.p5_z, params.p5_Bx, params.z);
    double Bx_topup = linear_Interpolation(params.p6_z, params.p6_Bx, params.p7_z, params.p7_Bx, params.z);
    double Bx_top = linear_Interpolation(params.p4_x, Bx_topdown, params.p6_x, Bx_topup, params.x);
    //y方向の補間
    double Bx = linear_Interpolation(params.p0_y, Bx_bottom, params.p4_y, Bx_top, params.y);
    //Byも同様
    //底面の下辺上辺の磁場を計算し、底面のz,xにおける磁場を計算
    double By_bottomdown = linear_Interpolation(params.p0_z, params.p0_By, params.p1_z, params.p1_By, params.z);
    double By_bottomup = linear_Interpolation(params.p2_z, params.p2_By, params.p3_z, params.p3_By, params.z);
    double By_bottom = linear_Interpolation(params.p0_x, By_bottomdown, params.p2_x, By_bottomup, params.x);
    //上面の下辺上辺の磁場を計算し、上面のz,xにおける磁場を計算
    double By_topdown = linear_Interpolation(params.p4_z, params.p4_By, params.p5_z, params.p5_By, params.z);
    double By_topup = linear_Interpolation(params.p6_z, params.p6_By, params.p7_z, params.p7_By, params.z);
    double By_top = linear_Interpolation(params.p4_x, By_topdown, params.p6_x, By_topup, params.x);
    //y方向の補間
    double By = linear_Interpolation(params.p0_y, By_bottom, params.p4_y, By_top, params.y);
    //Bzも同様
    //底面の下辺上辺の磁場を計算し、底面のz,xにおける磁場を計算
    double Bz_bottomdown = linear_Interpolation(params.p0_z, params.p0_Bz, params.p1_z, params.p1_Bz, params.z);
    double Bz_bottomup = linear_Interpolation(params.p2_z, params.p2_Bz, params.p3_z, params.p3_Bz, params.z);
    double Bz_bottom = linear_Interpolation(params.p0_x, Bz_bottomdown, params.p2_x, Bz_bottomup, params.x);
    //上面の下辺上辺の磁場を計算し、上面のz,xにおける磁場を計算
    double Bz_topdown = linear_Interpolation(params.p4_z, params.p4_Bz, params.p5_z, params.p5_Bz, params.z);
    double Bz_topup = linear_Interpolation(params.p6_z, params.p6_Bz, params.p7_z, params.p7_Bz, params.z);
    double Bz_top = linear_Interpolation(params.p4_x, Bz_topdown, params.p6_x, Bz_topup, params.x);
    //y方向の補間
    double Bz = linear_Interpolation(params.p0_y, Bz_bottom, params.p4_y, Bz_top, params.y);

    return make_tuple(Bx,By,Bz,B);
}

tuple <double,double,double,double> interpolated_Bfield(double x, double y, double z){
    /*_____________________________________________________________________________
    box_numは各軸方向の格子の番号
    それぞれ原点から正方向に0,1,2,3,,,,
                 負方向に1,2,3,,,,
    _______________________________________________________________________________*/
    //p0のインデックスからp1~p7のインデックスを計算
    int line_numX = floor(x/10) + 400;
    int line_numY = floor(y/10) + 100;
    int line_numZ = floor(z/10) + 100;
    int p0line_num =  (line_numX * 201 * 801) + (line_numY * 801) + (line_numZ);
    int p1line_num = p0line_num +1;
    int p2line_num = p0line_num + 201*801;
    int p3line_num = p2line_num + 1;
    int p4line_num = p0line_num + 801;
    int p5line_num = p4line_num +1;
    int p6line_num = p4line_num + 201*801;
    int p7line_num = p6line_num +1;

    TFile *file = TFile::Open("/Users/shohtatakami/physics/COMETDS/DS189A_944turns/DS189A_944turns_Map.root");
    if(!file || file->IsZombie()){
        cerr<<"Failed to open ROOT file!"<< endl;
        return make_tuple(0.0, 0.0, 0.0, 0.0);
    }
    TTree *tree = (TTree*)file->Get("tree");
    if(!tree){
        cerr<<"Failed to get TTree from file!" << endl;
        return make_tuple(0.0, 0.0, 0.0, 0.0);
    }
    //make variables correspoding to Branches of Tree
    double tx,ty,tz,tBx,tBy,tBz,tB;
    tree->SetBranchAddress("x", &tx);
    tree->SetBranchAddress("y", &ty);
    tree->SetBranchAddress("z", &tz);
    tree->SetBranchAddress("Bx", &tBx);
    tree->SetBranchAddress("By", &tBy);
    tree->SetBranchAddress("Bz", &tBz);
    tree->SetBranchAddress("B", &tB);

    weighted_ave_Param params;

    //read data of each line and set to struct
    auto set_values = [&](int line_num, double &px, double &py, double &pz, double &pBx, double &pBy, double &pBz, double &pB){
        tree->GetEntry(line_num);
        px =tx; py = ty; pz = tz;
        pBx = tBx; pBy = tBy; pBz = tBz; pB = tB;
    };
    
    //store data to struct
    set_values(p0line_num, params.p0_x, params.p0_y, params.p0_z, params.p0_Bx, params.p0_By, params.p0_Bz, params.p0_B);
    set_values(p1line_num, params.p1_x, params.p1_y, params.p1_z, params.p1_Bx, params.p1_By, params.p1_Bz, params.p1_B);
    set_values(p2line_num, params.p2_x, params.p2_y, params.p2_z, params.p2_Bx, params.p2_By, params.p2_Bz, params.p2_B);
    set_values(p3line_num, params.p3_x, params.p3_y, params.p3_z, params.p3_Bx, params.p3_By, params.p3_Bz, params.p3_B);
    set_values(p4line_num, params.p4_x, params.p4_y, params.p4_z, params.p4_Bx, params.p4_By, params.p4_Bz, params.p4_B);
    set_values(p5line_num, params.p5_x, params.p5_y, params.p5_z, params.p5_Bx, params.p5_By, params.p5_Bz, params.p5_B);
    set_values(p6line_num, params.p6_x, params.p6_y, params.p6_z, params.p6_Bx, params.p6_By, params.p6_Bz, params.p6_B);
    set_values(p7line_num, params.p7_x, params.p7_y, params.p6_z, params.p7_Bx, params.p7_By, params.p7_Bz, params.p7_B);
    //arbitrary point
    params.x = x;
    params.y = y;
    params.z = z;

    auto[Bx, By, Bz, B] = weighted_average(params);

    file->Close();
    return make_tuple(Bx,By,Bz,B);
}

int main(int argc, char** argv) {
    TApplication app("app", &argc, argv);
    vector <double> x,y,z;
    vector <double> Bx,By,Bz,B;
    vector <double> x_neg760,y_neg760,z_neg760;
    vector <double> Bx_neg760,By_neg760,Bz_neg760,B_neg760;
    //x1,y1,x2,y2の値を設定
    double x1 = 5.0; double x2 = 5.0; 
    double y1 = 0.0; double y2 =0.0;
    double L_dash = sqrt(pow(3220,2) + pow((y1+y2),2+pow((x1+x2),2)))/2;

    for(int i = 0; i < 3690; i += 10){
        double t = double(i)/(2*L_dash);
        double x_value = x2 - t*x1 - t*x2;
        double y_value = x2 - t*y1 - t*y2;
        double z_value = -L_dash + t*2*L_dash;;
        auto[Bx_value, By_value, Bz_value, B_value] = interpolated_Bfield(x_value,y_value,z_value);
        x.push_back(x_value); y.push_back(y_value); z.push_back(z_value);
        Bx.push_back(Bx_value); By.push_back(By_value); Bz.push_back(Bz_value); B.push_back(B_value);
    }

//X=-760 mmの傾いたScan軸上の座標
    for(int i = 0; i < x.size(); ++i){
        // x方向に-760 mm平行移動
        double x_value = x[i] - 760;
        double y_value = y[i]; double z_value = z[i];
        auto[Bx_value, By_value, Bz_value, B_value] = interpolated_Bfield(x_value,y_value,z_value);
        x_neg760.push_back(x_value); y_neg760.push_back(y_value); z_neg760.push_back(z_value);
        Bx_neg760.push_back(Bx_value); By_neg760.push_back(By_value); 
        By_neg760.push_back(By_value); B_neg760.push_back(B_value); 
    }

    TCanvas *c0 = new TCanvas( "c0" , "center" , 800 , 600 );
    TCanvas *cxneg760 = new TCanvas( "cx-760" , "X=-760" , 800 , 600 );
    TF1 *quadfit = new TF1("quadfit", "[0]*(x-[1])^2+[2]",-150,150);
    quadfit->SetParameters(0,0.0,0.85);
    quadfit->SetLineColor(kCyan);
    quadfit->SetLineWidth(2);

    c0->cd();
    TGraph *g0 = new TGraph(z.size(), &z[0], &B[0]);
    g0->SetMarkerStyle(6); g0->GetXaxis()->SetTitle("z' (mm)");
    g0->SetMarkerColor(kBlue); g0->GetYaxis()->SetTitle("B (T)");
    g0->Draw("AP");
    c0->Update();
    cxneg760->cd();
    TGraph *gxneg760 = new TGraph(z_neg760.size(), &z_neg760[0], &B[0]);
    gxneg760->SetMarkerStyle(6); gxneg760->GetXaxis()->SetTitle("z' (mm)");
    gxneg760->SetMarkerColor(kRed); gxneg760->GetYaxis()->SetTitle("B (T)");
    gxneg760->Draw("AP");
    cxneg760->Update();
    app.Run();
    return 0;
}