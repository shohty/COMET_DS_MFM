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
#include <chrono>
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

/*
格子のxy平面の側面にzは固定する
*/
double linear_Interpolation(double left, double left_B, double right, double right_B, double x){
    //線形補間の関数
    // leftはdown,bottomに対応
    // rightはup,topに対応
    double l = right -left;
    double B = ((right - x)/l) * left_B + ((x - left)/l) * right_B;//内分計算
    if(abs(right -x) > 10){
        cout<<"WRONG !!"<<"| right : "<<right<<"| x : "<< x <<endl;
    }
    if(abs(x - left) > 10){
        cout<<"WRONG !!"<<"| left : "<<left<< "| x :"<< x <<endl;
    }
    return B;
}

struct weighted_ave_Param{
    double p0_x,p0_y,p0_z,p1_x,p1_y,p1_z,p2_x,p2_y,p2_z,p3_x,p3_y,p3_z;
    double p0_B,p1_B,p2_B,p3_B;
    double p0_Bx,p1_Bx,p2_Bx,p3_Bx;
    double p0_By,p1_By,p2_By,p3_By;
    double p0_Bz,p1_Bz,p2_Bz,p3_Bz;
    double x, y, z; //任意の点の座標
};

tuple<double, double, double, double> weighted_average(const weighted_ave_Param& params){
    /*_____________________________________________________________________________
    斜めのスキャン軸上の任意の点Pのうちz= -4000,-3990,,,,,3990,4000の点について考える
    たとえばz=1000の場合z=1000nの平面上でPが所属する格子(二次元)内でバイリニア補間
    任意の点(x,y,z)が所属する三次元の格子を考える
    DSをDAQシステム側、z軸向かって左側から見る
    (z,x,y)の右手系
    y軸を垂直上方向、x軸を右方向となるように見たときに
    左下p0,右下p1,左上p2,右上p3とする
    x正方向 right : y正方向 up
    x負方向 left  : y負方向 down
    _______________________________________________________________________________*/
    //Bx
    double Bx_down = linear_Interpolation(params.p0_x, params.p0_Bx, params.p1_x, params.p1_Bx, params.x);
    double Bx_up = linear_Interpolation(params.p2_x, params.p2_Bx, params.p3_x, params.p3_Bx, params.x);
    double Bx = linear_Interpolation(params.p0_y, Bx_down, params.p2_y, Bx_up, params.y);
    //By
    double By_down = linear_Interpolation(params.p0_x, params.p0_By, params.p1_x, params.p1_By, params.x);
    double By_up = linear_Interpolation(params.p2_x, params.p2_By, params.p3_x, params.p3_By, params.x);
    double By = linear_Interpolation(params.p0_y, By_down, params.p2_y, By_up, params.y);
    //Bz
    double Bz_down = linear_Interpolation(params.p0_x, params.p0_Bz, params.p1_x, params.p1_Bz, params.x);
    double Bz_up = linear_Interpolation(params.p2_x, params.p2_Bz, params.p3_x, params.p3_Bz, params.x);
    double Bz = linear_Interpolation(params.p0_y, Bz_down, params.p2_y, Bz_up, params.y);
    //B
    double B_down = linear_Interpolation(params.p0_x, params.p0_B, params.p1_x, params.p1_B, params.x);
    double B_up = linear_Interpolation(params.p2_x, params.p2_B, params.p3_x, params.p3_B, params.x);
    double B = linear_Interpolation(params.p0_y, B_down, params.p2_y, B_up, params.y);

    return make_tuple(Bx,By,Bz,B);
}

tuple <double,double,double,double> interpolated_Bfield(double x, double y, double z){
    /*_____________________________________________________________________________
    box_numは各軸方向の格子の番号
    それぞれ原点から正方向に0,1,2,3,,,,
                 負方向に1,2,3,,,,
    _______________________________________________________________________________*/
    //p0のインデックスからp1~p7のインデックスを計算
    //p0のx,y,zの番
    //各軸最も小さい値から0,1,2,,,,,と数えている
    int line_numX = floor(x/10) + 100;
    int line_numY = floor(y/10) + 100;
    int line_numZ = floor(z/10) + 400;
    int p0line_num = (line_numX * 201 * 801) + (line_numY * 801) + line_numZ;
    int p1line_num = ((line_numX + 1) * 201 * 801) + (line_numY * 801) + line_numZ;
    int p2line_num = (line_numX * 201 * 801) + ((line_numY + 1) * 801) + line_numZ;
    int p3line_num = ((line_numX + 1) * 201 * 801) + ((line_numY + 1) * 801) + line_numZ;
    


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
    //arbitrary point
    params.x = x;
    params.y = y;
    params.z = z;
    auto[Bx, By, Bz, B] = weighted_average(params);
    /*cout<<"x : "<<params.x<<"y : "<<params.y <<"z : "<<params.z<<"\n"<<
    "p0_X : "<< params.p0_x<<" | p0_Y : "<<params.p0_y<<" | p0_Z : "<<params.p0_z<<"\n"<<
    "p1_X : "<< params.p1_x<<" | p1_Y : "<<params.p1_y<<" | p1_Z : "<<params.p1_z<<"\n"<<
    "p2_X : "<< params.p2_x<<" | p2_Y : "<<params.p2_y<<" | p2_Z : "<<params.p2_z<<"\n"<<
    "p3_X : "<< params.p3_x<<" | p3_Y : "<<params.p3_y<<" | p3_Z : "<<params.p3_z<<"\n"<<endl;
    cout<<"x : "<<params.x<<" | y : "<<params.y <<" | z : "<<params.z<<endl;*/
    file->Close();
    return make_tuple(Bx,By,Bz,B);
}

int main(int argc, char** argv) {
    auto start = chrono::high_resolution_clock::now();
    TApplication app("app", &argc, argv);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    //x1,y1,x2,y2の値を設定
    /*傾いた軸上の任意の点
    ->  ->    ->  ->
    p = b + t(a - b)    
    x   x2 + t(x1-x2)
    y = y2 + t(y1-y2)
    z   -L + t*(L-(-L))   L=1660 mm 
    */
    double x1 = -10.0; double x2 = 10.0; 
    double y1 = 0.0; double y2 = 0.0;
    //double L_dash = sqrt(pow(3220,2) + pow((y1+y2),2)+pow((x1+x2),2))/2;
    vector <double> x,y,z;//斜めの中心軸においてzを10mmピッチで（-4000,-3990,・・・3990,4000）動かした時の各座標
    vector <double> Bx,By,Bz,B;//上のx,y,z座標においてinterpolateして計算した磁場
    vector <double> x_neg760,y_neg760,z_neg760;//x=-760mm斜め軸
    vector <double> Bx_neg760,By_neg760,Bz_neg760,B_neg760;//x=-760mm斜め軸上のinterpolateして計算した磁場
    //中心軸における斜めScan　zは10mmピッチで-4000から4000
    for(int i = -4000; i <4000; i +=10){//x,yが＋ー1000を超えたのは省くように改良するべき
        double z_value = i;
        double t_value = ((z_value/1660) + 1)/2;//(z+1660/3220)
        double x_value = x2 + t_value*(x1-x2);
        double y_value = y2 + t_value*(y1-y2);
        auto[Bx_value, By_value, Bz_value, B_value] = interpolated_Bfield(x_value, y_value, z_value);
        x.push_back(x_value);
        y.push_back(y_value); 
        z.push_back(z_value);
        Bx.push_back(Bx_value); 
        By.push_back(By_value); 
        Bz.push_back(Bz_value); 
        B.push_back(B_value);
    }
    //x=-760mmにおける斜めScan(ずれ方は同じとする)
    for(int i = 0; i < z.size(); ++i){
        double x_value = x[i] - 760;
        double y_value = y[i];
        double z_value = z[i];
        auto[Bx_value, By_value, Bz_value, B_value] = interpolated_Bfield(x_value, y_value, z_value);
        x_neg760.push_back(x_value); 
        y_neg760.push_back(y_value);
        z_neg760.push_back(z_value);
        Bx_neg760.push_back(Bx_value); 
        By_neg760.push_back(By_value); 
        Bz_neg760.push_back(Bz_value); 
        B_neg760.push_back(B_value);
    }
    TCanvas *c0 = new TCanvas( "c0" , "center" , 800 , 600 );
    //TCanvas *c1 = new TCanvas("c1", "x=760", 800, 600);
    TCanvas *cxneg760 = new TCanvas( "cx-760" , "X=-760" , 800 , 600 );
    TCanvas *c3D_center = new TCanvas("c3D" , "scatter plot of points @ center" , 2400 , 1800);
    TCanvas *c3D_neg760 = new TCanvas("c3D_neg760" , " scatter plot of points @ x=-760" , 2400 , 1800);
    TF1 *quadfit = new TF1("quadfit", "[0]*(x-[1])^2+[2]",-150,150);
    quadfit->SetParameters(0,0.0,0.85);
    quadfit->SetLineColor(kCyan);
    quadfit->SetLineWidth(2);

    c0->cd();
    TGraph *g0 = new TGraph(z.size(), &z[0], &B[0]);
    g0->SetTitle("center");
    g0->SetMarkerStyle(6); 
    g0->GetXaxis()->SetTitle("z' (mm)");
    g0->SetMarkerColor(kBlue); 
    g0->GetYaxis()->SetTitle("B (T)");
    g0->Draw("AP");
    g0->Fit(quadfit,"R");
    c0->SetGridx();
    c0->SetGridy();
    c0->Update();

    cxneg760->cd();
    TGraph *gxneg760 = new TGraph(z_neg760.size(), &z_neg760[0], &B[0]);
    string titleneg760 = "x=-760 : x1 = " + to_string(x1) + "x2 = " + to_string(x2);
    gxneg760->SetTitle(titleneg760.c_str());
    gxneg760->SetMarkerStyle(6); 
    gxneg760->SetMarkerColor(kRed); 
    gxneg760->GetXaxis()->SetTitle("z' (mm)");
    gxneg760->GetYaxis()->SetTitle("B (T)");
    gxneg760->Draw("AP");
    gxneg760->Fit(quadfit,"R");
    cxneg760->Update();

    // TGraph2Dを作成
    /*c3D_center->cd();
    TGraph2D *g3D = new TGraph2D(x.size(), &z[0], &x[0], &y[0]);
    g3D->SetMarkerColor(kAzure);
    g3D->SetMarkerStyle(6);
    g3D->SetTitle("3D Scatter Plot");
    g3D->GetXaxis()->SetTitle("z (mm)");
    g3D->GetYaxis()->SetTitle("x (mm)");
    g3D->GetZaxis()->SetTitle("y (mm)");
    g3D->Draw("P");  // "P0" はポイントのみ表示
    c3D_center->Update();

    c3D_neg760->cd();
    TGraph2D *g3D_neg760 = new TGraph2D(x.size(), &z_neg760[0], &x_neg760[0], &y_neg760[0]);
    g3D_neg760->SetMarkerColor(kAzure);
    g3D_neg760->SetMarkerStyle(6);
    g3D_neg760->SetTitle("3D Scatter Plot");
    g3D_neg760->GetXaxis()->SetTitle("z (mm)");
    g3D_neg760->GetYaxis()->SetTitle("x (mm)");
    g3D_neg760->GetZaxis()->SetTitle("y (mm)");
    g3D_neg760->Draw("P");  // "P0" はポイントのみ表示
    c3D_neg760->Update();*/
    app.Run();
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    cout << "Execution time: " << elapsed.count() << " seconds" << endl;
    return 0;
}