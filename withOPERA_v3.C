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
#include <TFile.h>
#include <TTree.h>
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
using namespace std;
/*________________________________________________
version2で
通常の二次関数Fit、Bmax,OPERAに測定データのエラーバーを載せたものの三種類でpeakpositonを比較する機能を追加
plot関数に測定とOPERAの測定の差と比を出すコードを追加
______________________________________________________*/
void plot(const char* meas_file="test20240626-6.root", const char* opera_file="X_0Y_0.csv_Map.csv.root"){
    
    string FilenameStr = opera_file;
    FilenameStr.erase(0,FilenameStr.find_last_of('/') + 1);
    string directory = "/Users/shohtatakami/physics/COMETDS/DS189A_944turns/Bdiffratio/";
    string Filename_B = directory + FilenameStr + "_B.pdf";
    string Filename_Bdiffratio = directory + FilenameStr + "_Bdiffratio .pdf";
    // 各ROOTファイルを開く
    TFile *measRootFile = TFile::Open(meas_file);
    if (!measRootFile || measRootFile->IsZombie()){
        cerr << "Failed to open measured ROOT File!" << endl;
        return;
    }
    TFile *operaRootFile = TFile::Open(opera_file);  
    if (!operaRootFile || operaRootFile->IsZombie()) {
        cerr << "Failed to open opera ROOT file!" << endl;
        return;
    }

    // TTreeを取得
    TTree *meastree = (TTree*)measRootFile->Get("tree");  // "tree" を実際のツリーネームに変更
    if (!meastree) {
        cerr << "Failed to get TTree from measured file!" << endl;
        return;
    }
    TTree *operatree = (TTree*)operaRootFile->Get("tree");  // "tree" を実際のツリーネームに変更
    if (!operatree) {
        cerr << "Failed to get TTree from opera file!" << endl;
        return;
    }
    //TTreeのBranchをセット(measured)
    double xmean, xstdev, ymean, ystdev, zmean, zstdev, xrange, yrange, zrange;
    meastree->SetBranchAddress("Xmean", &xmean);
    meastree->SetBranchAddress("Xstdev", &xstdev);
    meastree->SetBranchAddress("Ymean", &ymean);
    meastree->SetBranchAddress("Ystdev", &ystdev);
    meastree->SetBranchAddress("Zmean", &zmean);
    meastree->SetBranchAddress("Zstdev", &zstdev);
    meastree->SetBranchAddress("Xrange", &xrange);
    meastree->SetBranchAddress("Yrange", &yrange);
    meastree->SetBranchAddress("Zrange", &zrange);
    // TTreeのBranchをセット(OPERA)
    double x_val,y_val,z_val,operaBx_val,operaBy_val,operaBz_val,operaB_val; //必要な変数を定義
    operatree->SetBranchAddress("X", &x_val);
    operatree->SetBranchAddress("Y", &y_val);
    operatree->SetBranchAddress("Z", &z_val);
    operatree->SetBranchAddress("Bx", &operaBx_val);
    operatree->SetBranchAddress("By", &operaBy_val);
    operatree->SetBranchAddress("Bz", &operaBz_val);
    operatree->SetBranchAddress("B", &operaB_val);

    //データを読み込む(measured)
    vector<double> measBx,measBy,measBz,measB,measBxerr,measByerr,measBzerr,measBerr;
    Long64_t meas_nentries = meastree->GetEntries();
    for(Long64_t i=0; i < meas_nentries; ++i){
        meastree->GetEntry(i);
        double Bx = -xmean * xrange/5; 
        double By = ymean * yrange/5; 
        double Bz = zmean * zrange/5;
        double B = sqrt(pow(Bx,2)+pow(By,2)+pow(Bz,2));
        double Bxerr = xstdev * xrange/5;
        double Byerr = ystdev * yrange/5; 
        double Bzerr = zstdev * zrange/5;
        double Berr = sqrt(pow(Bx/B,2)*pow(Bxerr,2)+pow(By/B,2)*pow(Byerr,2)+pow(Bz/B,2)*pow(Bzerr,2));
        measBx.push_back(Bx);
        measBy.push_back(By);
        measBz.push_back(Bz);
        measB.push_back(B);
        measBxerr.push_back(Bxerr);
        measByerr.push_back(Byerr);
        measBzerr.push_back(Bzerr);
        measBerr.push_back(Berr);
    }

    // データを読み込む(OPERA)
    vector<double> x,y,z;
    vector<double> operaBx,operaBy,operaBz,operaB;
    Long64_t opera_nentries = operatree->GetEntries();
    for (Long64_t i = 0; i < opera_nentries; ++i){
        operatree->GetEntry(i);
        x.push_back(x_val);
        y.push_back(y_val);
        z.push_back(z_val);
        operaBx.push_back(operaBx_val);
        operaBy.push_back(operaBy_val);
        operaBz.push_back(operaBz_val);
        operaB.push_back(operaB_val);
    }
    //測定とデータの差分と比
    vector<double> Bdiff,Bratio;
    for(int i = 0; i < z.size(); ++i){
        double Bdiff_tmp =  measB[i] - operaB[i] ;
        double Bratio_tmp = measB[i]/operaB[i];
        Bdiff.push_back(Bdiff_tmp);
        Bratio.push_back(Bratio_tmp);
    }

    //描画
    TCanvas *cB = new TCanvas(FilenameStr.c_str(), FilenameStr.c_str(), 800, 600);
    TCanvas *c_diffratio = new TCanvas("diffratio", "Bdiffratio", 1000, 3000);
    cB->cd();
    TGraphErrors *g_measB = new TGraphErrors(z.size(), &z[0], &measB[0],0,&measBerr[0]);
    TGraph *g_operaB = new TGraph(z.size(), &z[0], &operaB[0]);
    g_measB->SetMarkerColor(kBlue);
    g_measB->SetMarkerStyle(20);
    g_operaB->SetMarkerColor(kRed+2);
    g_operaB->SetMarkerStyle(20);
    //二次関数で範囲を絞ってフィット
    //measured
    TF1* fit_meas = new TF1("fitFunc", "[0]*(x-[1])^2+[2]", -150, 150);
    fit_meas->SetParameters(0,0.0,0.85);
    fit_meas->SetLineColor(kBlue);
    fit_meas->SetLineWidth(2);
    g_measB->Fit(fit_meas, "R");
    //フィットパラメータを取得
    double p0_meas = fit_meas->GetParameter(0);
    double p1_meas = fit_meas->GetParameter(1);
    double p1_meas_err = fit_meas->GetParError(1);
    double p2_meas = fit_meas->GetParameter(2);
    double p2_meas_err = fit_meas->GetParError(2);
    double chi2_meas = fit_meas->GetChisquare();
    //ピーク位置を計算 (y = ax^2 + bx + c の頂点の位置 x = -b/(2a) を使う)
    double peakPosition_meas = p1_meas;
    //OPERA
    TF1* fit_opera = new TF1("fit_simu", "[0]*(x-[1])^2+[2]", -150, 150);
    fit_opera->SetParameters(0,0.0,0.85);
    fit_opera->SetLineColor(kRed+2);
    fit_opera->SetLineWidth(2);
    g_operaB->Fit(fit_opera, "R");
    //フィットパラメータを取得
    double p0_opera = fit_opera->GetParameter(0);
    double p1_opera = fit_opera->GetParameter(1);
    double p1_opera_err = fit_opera->GetParError(1);
    double p2_opera = fit_opera->GetParameter(2);
    double p2_opera_err = fit_opera->GetParError(2);
    double chi2_opera = fit_opera->GetChisquare();
    //ピーク位置を計算 (y = ax^2 + bx + c の頂点の位置 x = -b/(2a) を使う)
    double peakPosition_simu = p1_opera;

    TMultiGraph *mgB=new TMultiGraph();
    mgB->Add(g_measB);
    mgB->Add(g_operaB);
    mgB->SetTitle(FilenameStr.c_str());
    mgB->GetXaxis()->SetTitle("z (mm)");
    mgB->GetXaxis()->SetLimits(-1200.0,1200.0);
    mgB->GetYaxis()->SetTitle("B (T)");
    mgB->GetYaxis()->SetRangeUser(0.7,0.92);
    mgB->Draw("AP");
    
    cB->SetGridx();
    cB->SetGridy();
    // カスタムテキストをカンバス上に表示する
    TLatex latex;
    latex.SetNDC();  // カンバス上のNDC座標系で位置を指定
    latex.SetTextSize(0.03);  // 文字サイズを設定 (0〜1の範囲)
    latex.SetTextColor(kRed+2);

    // 測定データのχ²とp1を表示
    latex.DrawLatex(0.7, 0.9, Form("#chi^{2} (measured) = %.10e", chi2_meas));
    latex.DrawLatex(0.7, 0.85, Form("p1 (measured) = %.5f #pm %.5f", p1_meas,p1_meas_err));
    latex.DrawLatex(0.7, 0.80, Form("p2 (measured) = %.5f #pm %.5f", p2_meas,p2_meas_err));
    
    // OPERAデータのχ²とp1を表示
    latex.DrawLatex(0.7, 0.75, Form("#chi^{2} (OPERA) = %.10e", chi2_opera));
    latex.DrawLatex(0.7, 0.70, Form("p1 (OPERA) = %.5f #pm %.5f", p1_opera,p1_opera_err));
    latex.DrawLatex(0.7, 0.65, Form("p2 (OPERA) = %.5f #pm %.5f", p2_opera,p2_opera_err));

    cB->Update();
    c_diffratio->Divide(1,2);
    c_diffratio->cd(1);
    TGraph *g_Bdiff = new TGraph(z.size(), &z[0], &Bdiff[0]);
    g_Bdiff->SetMarkerColor(kBlue);
    g_Bdiff->SetMarkerStyle(6);
    g_Bdiff->SetTitle("measured - opera");
    g_Bdiff->GetXaxis()->SetTitle("z (mm)");
    g_Bdiff->GetYaxis()->SetTitle("Bdiff (T)");
    g_Bdiff->Draw("AP");
    c_diffratio->cd(1)->SetGridx();
    c_diffratio->cd(1)->SetGridy();
    c_diffratio->cd(2);
    TGraph *g_Bratio = new TGraph(z.size(), &z[0], &Bratio[0]);
    g_Bratio->SetMarkerColor(kRed+2);
    g_Bratio->SetMarkerStyle(6);
    g_Bratio->SetTitle("measured / opera");
    g_Bratio->GetXaxis()->SetTitle("z (mm)");
    g_Bratio->GetYaxis()->SetTitle("Bratio");
    g_Bratio->Draw("AP");
    c_diffratio->cd(2)->SetGridx();
    c_diffratio->cd(2)->SetGridy();
    c_diffratio->Update();
    c_diffratio->SaveAs(Filename_Bdiffratio.c_str());
    //cB->SaveAs(Filename_B.c_str());
    string rootplot_directory = "physics/COMETDS/DS189A_944turns/LaserTrackerPosition/rootoutput/";
    auto outFile = new TFile(Form("%s%s_out.root", rootplot_directory.c_str(), FilenameStr.c_str()),"RECRREATE"); //オブジェクトをrootで保存
    cB->Write();
    outFile->Close();
}

tuple<double,double,double,double,double,double,double,double,double>peakposi(const char* meas_file="test20240626-6.root", const char* opera_file="X_0Y_0.csv_Map.csv.root"){
    string FilenameStr = opera_file;
    FilenameStr.erase(0,FilenameStr.find_last_of('/') + 1);
    //string directory = "/Users/shohtatakami/physics/COMETDS/DS189A_FactoryPillar_LaserTrackerPos/";
    //string Filename_B = directory + FilenameStr + "_B.pdf"; 
    // 各ROOTファイルを開く
    TFile *measRootFile = TFile::Open(meas_file);
    if (!measRootFile || measRootFile->IsZombie()){
        cerr << "Failed to open measured ROOT File!" << endl;
        return make_tuple(0, 0, 0, 0, 0, 0, 0, 0, 0); // 早期リターン
    }
    TFile *operaRootFile = TFile::Open(opera_file);  
    if (!operaRootFile || operaRootFile->IsZombie()) {
        cerr << "Failed to open opera ROOT file!" << endl;
        return make_tuple(0, 0, 0, 0, 0, 0, 0, 0, 0); // 早期リターン
    }

    // TTreeを取得
    TTree *meastree = (TTree*)measRootFile->Get("tree");  // "tree" を実際のツリーネームに変更
    if (!meastree) {
        cerr << "Failed to get TTree from measured file!" << endl;
        return make_tuple(0, 0, 0, 0, 0, 0, 0, 0, 0); // 早期リターン
    }
    TTree *operatree = (TTree*)operaRootFile->Get("tree");  // "tree" を実際のツリーネームに変更
    if (!operatree) {
        cerr << "Failed to get TTree from opera file!" << endl;
        return make_tuple(0, 0, 0, 0, 0, 0, 0, 0, 0); // 早期リターン
    }
    //TTreeのBranchをセット(measured)
    double xmean, xstdev, ymean, ystdev, zmean, zstdev, xrange, yrange, zrange;
    meastree->SetBranchAddress("Xmean", &xmean);
    meastree->SetBranchAddress("Xstdev", &xstdev);
    meastree->SetBranchAddress("Ymean", &ymean);
    meastree->SetBranchAddress("Ystdev", &ystdev);
    meastree->SetBranchAddress("Zmean", &zmean);
    meastree->SetBranchAddress("Zstdev", &zstdev);
    meastree->SetBranchAddress("Xrange", &xrange);
    meastree->SetBranchAddress("Yrange", &yrange);
    meastree->SetBranchAddress("Zrange", &zrange);
    // TTreeのBranchをセット(OPERA)
    double x_val,y_val,z_val,operaBx_val,operaBy_val,operaBz_val,operaB_val; //必要な変数を定義
    operatree->SetBranchAddress("X", &x_val);
    operatree->SetBranchAddress("Y", &y_val);
    operatree->SetBranchAddress("Z", &z_val);
    operatree->SetBranchAddress("Bx", &operaBx_val);
    operatree->SetBranchAddress("By", &operaBy_val);
    operatree->SetBranchAddress("Bz", &operaBz_val);
    operatree->SetBranchAddress("B", &operaB_val);

    //データを読み込む(measured)
    vector<double> measBx,measBy,measBz,measB,measBxerr,measByerr,measBzerr,measBerr;
    Long64_t meas_nentries = meastree->GetEntries();
    for(Long64_t i=0; i < meas_nentries; ++i){
        meastree->GetEntry(i);
        double Bx = -xmean * xrange/5; 
        double By = ymean * yrange/5; 
        double Bz = zmean * zrange/5;
        double B = sqrt(pow(Bx,2)+pow(By,2)+pow(Bz,2));
        double Bxerr = -xstdev * xrange/5;
        double Byerr = ystdev * yrange/5; 
        double Bzerr = zstdev * zrange/5;
        double Berr = sqrt(pow(Bx/B,2)*pow(Bxerr,2)+pow(By/B,2)*pow(Byerr,2)+pow(Bz/B,2)*pow(Bzerr,2));
        measBx.push_back(Bx);
        measBy.push_back(By);
        measBz.push_back(Bz);
        measB.push_back(B);
        measBxerr.push_back(Bxerr);
        measByerr.push_back(Byerr);
        measBzerr.push_back(Bzerr);
        measBerr.push_back(Berr);
    }

    // データを読み込む(OPERA)
    vector<double> x,y,z;
    vector<double> operaBx,operaBy,operaBz,operaB;
    Long64_t opera_nentries = operatree->GetEntries();
    for (Long64_t i = 0; i < opera_nentries; ++i){
        operatree->GetEntry(i);
        x.push_back(x_val);
        y.push_back(y_val);
        z.push_back(z_val);
        operaBx.push_back(operaBx_val);
        operaBy.push_back(operaBy_val);
        operaBz.push_back(operaBz_val);
        operaB.push_back(operaB_val);
    }
    vector<double>::iterator Bmax_itr= max_element(begin(operaB),end(operaB));
    double Bmax_tmp;
    double Bmaxposi_tmp;
    int Bmax_indent;
    if(Bmax_itr != end(operaB)){
        Bmax_tmp = *Bmax_itr;
        Bmax_indent = distance(begin(operaB),Bmax_itr); 
        Bmaxposi_tmp = z[Bmax_indent];
        if (Bmax_indent >= 0 && Bmax_indent < z.size()) {  // **範囲チェック**
        double Bmaxposi_tmp = z[Bmax_indent];// **ここで範囲外アクセスの可能性**
        } else {
            cerr << "Error: Bmax_indent out of range!" << endl;
            }
        cout<<"Bmax : " << Bmax_tmp<<"\n" <<"position of Bmax : "<< Bmaxposi_tmp <<endl;
    }

    //描画
    TGraphErrors *g_measB = new TGraphErrors(z.size(), &z[0], &measB[0], 0, &measBerr[0]);
    TGraph *g_operaB = new TGraph(z.size(), &z[0], &operaB[0]);
    TGraphErrors *g_operaBerr = new TGraphErrors(z.size(), &z[0], &operaB[0], 0, &measBerr[0]);
    //OPERAデータにエラーがついていないことでFitがうまくいっていない?->測定と同じエラーをOPERAにもつけよう
    g_measB->SetMarkerColor(kBlue);
    g_measB->SetMarkerStyle(6);
    g_operaB->SetMarkerColor(kRed+2);
    g_operaB->SetMarkerStyle(6);
    g_operaBerr->SetMarkerColor(kGreen);
    g_operaBerr->SetMarkerStyle(6);
    //二次関数で範囲を絞ってフィット
    //measured
    TF1* fit_meas = new TF1("fitFunc", "[0]*(x-[1])^2+[2]", -150, 150);
    fit_meas->SetParameters(0,0.0,0.85);
    fit_meas->SetLineColor(kBlue);
    fit_meas->SetLineWidth(2);
    g_measB->Fit(fit_meas, "R");
    //フィットパラメータを取得
    double p0_meas = fit_meas->GetParameter(0);
    double p1_meas = fit_meas->GetParameter(1);
    double p1_meas_err = fit_meas->GetParError(1);
    double p2_meas = fit_meas->GetParameter(2);
    double p2_meas_err = fit_meas->GetParError(2);
    double chi2_meas = fit_meas->GetChisquare();
    //OPERA
    TF1* fit_opera = new TF1("fit_simu", "[0]*(x-[1])^2+[2]", -150, 150);
    fit_opera->SetParameters(0,0,0.85);
    fit_opera->SetLineColor(kRed);
    fit_opera->SetLineWidth(2);
    g_operaB->Fit(fit_opera, "R");
    //フィットパラメータを取得
    double p0_opera = fit_opera->GetParameter(0);
    double p1_opera = fit_opera->GetParameter(1);
    double p1_opera_err = fit_opera->GetParError(1);
    double p2_opera = fit_opera->GetParameter(2);
    double p2_opera_err = fit_opera->GetParError(2);
    double chi2_opera = fit_opera->GetChisquare();
    //OPERA with error　_fiterrはフィッティングエラー
    TF1* fit_operaerr = new TF1("fit_simu", "[0]*(x-[1])^2+[2]", -150, 150);
    fit_operaerr->SetParameters(0,0,0.85);
    fit_operaerr->SetLineColor(kGreen);
    fit_operaerr->SetLineWidth(2);
    g_operaBerr->Fit(fit_operaerr, "R");
    //フィットパラメータを取得
    double p0_operaerr = fit_operaerr->GetParameter(0);
    double p1_operaerr = fit_operaerr->GetParameter(1);
    double p1_operaerr_fiterr = fit_operaerr->GetParError(1);
    double p2_operaerr = fit_operaerr->GetParameter(2);
    double p2_operaerr_fiterr = fit_operaerr->GetParError(2);
    double chi2_operaerr = fit_operaerr->GetChisquare();

    double sumx = accumulate(x.begin(),x.end(),0.0);
    double avex = sumx/x.size();
    cout<< "successful excution of peakposi()" <<endl;
    cout<< FilenameStr.c_str() <<endl;
    return make_tuple(avex,p1_meas,p1_meas_err,p1_opera,p1_opera_err, p1_operaerr, p1_operaerr_fiterr, Bmax_tmp, Bmaxposi_tmp);
}

tuple<double,double> err_evaluate(const char* meas_file="test20240626-6.root", const char* opera_file = "X_0Y_0.csv.root"){
    //Fit範囲-150 mm < z <150 mmでのBのerrを各軸で平均して比較
    string directory = "/Users/shohtatakami/physics/COMETDS/DS189A_FactoryPillar_LaserTrackerPos/MagneticFieldPlot/";
    // ROOTファイルを開く
    TFile *measRootFile = TFile::Open(meas_file);
    if (!measRootFile || measRootFile->IsZombie()){
        cerr << "Failed to open measured ROOT File!" << endl;
    }
    TFile *operaRootFile = TFile::Open(opera_file);  
    if (!operaRootFile || operaRootFile->IsZombie()) {
        cerr << "Failed to open opera ROOT file!" << endl;
    }

    // TTreeを取得
    TTree *meastree = (TTree*)measRootFile->Get("tree");  // "tree" を実際のツリーネームに変更
    if (!meastree) {
        cerr << "Failed to get TTree from measured file!" << endl;
    }
    TTree *operatree = (TTree*)operaRootFile->Get("tree");  // "tree" を実際のツリーネームに変更
    if (!operatree) {
        cerr << "Failed to get TTree from opera file!" << endl;
    }
    //TTreeのBranchをセット(measured)
    double xmean, xstdev, ymean, ystdev, zmean, zstdev, xrange, yrange, zrange;
    meastree->SetBranchAddress("Xmean", &xmean);
    meastree->SetBranchAddress("Xstdev", &xstdev);
    meastree->SetBranchAddress("Ymean", &ymean);
    meastree->SetBranchAddress("Ystdev", &ystdev);
    meastree->SetBranchAddress("Zmean", &zmean);
    meastree->SetBranchAddress("Zstdev", &zstdev);
    meastree->SetBranchAddress("Xrange", &xrange);
    meastree->SetBranchAddress("Yrange", &yrange);
    meastree->SetBranchAddress("Zrange", &zrange);
    // TTreeのBranchをセット(OPERA)
    double x_val,y_val,z_val,operaBx_val,operaBy_val,operaBz_val,operaB_val; //必要な変数を定義
    operatree->SetBranchAddress("X", &x_val);
    operatree->SetBranchAddress("Y", &y_val);
    operatree->SetBranchAddress("Z", &z_val);
    operatree->SetBranchAddress("Bx", &operaBx_val);
    operatree->SetBranchAddress("By", &operaBy_val);
    operatree->SetBranchAddress("Bz", &operaBz_val);
    operatree->SetBranchAddress("B", &operaB_val);


    //データをvectorへ(measured)
    vector<double> measBx,measBy,measBz,measB,measBxerr,measByerr,measBzerr,measBerr;
    Long64_t meas_nentries = meastree->GetEntries();
    for(Long64_t i=0; i < meas_nentries; ++i){
        meastree->GetEntry(i);
        double Bx = -xmean * xrange/5; 
        double By = ymean * yrange/5; 
        double Bz = zmean * zrange/5;
        double B = sqrt(pow(Bx,2)+pow(By,2)+pow(Bz,2));
        double Bxerr = xstdev * xrange/5;
        double Byerr = ystdev * yrange/5; 
        double Bzerr = zstdev * zrange/5;
        double Berr = sqrt(pow(Bx/B,2)*pow(Bxerr,2)+pow(By/B,2)*pow(Byerr,2)+pow(Bz/B,2)*pow(Bzerr,2));
        measBx.push_back(Bx);
        measBy.push_back(By);
        measBz.push_back(Bz);
        measB.push_back(B);
        measBxerr.push_back(Bxerr);
        measByerr.push_back(Byerr);
        measBzerr.push_back(Bzerr);
        measBerr.push_back(Berr);
    }
    // データをvectorへ(OPERA)
    vector<double> x,y,z;
    vector<double> operaBx,operaBy,operaBz,operaB;
    Long64_t opera_nentries = operatree->GetEntries();
    for (Long64_t i = 0; i < opera_nentries; ++i){
        operatree->GetEntry(i);
        x.push_back(x_val);
        y.push_back(y_val);
        z.push_back(z_val);
        operaBx.push_back(operaBx_val);
        operaBy.push_back(operaBy_val);
        operaBz.push_back(operaBz_val);
        operaB.push_back(operaB_val);
    }
    double sum = 0;
    double count = 0;
    for(int i=0; i<meas_nentries; ++i){
        if(z[i] >= -150 && z[i] <=150){ 
            count = count + 1;
            sum = sum + measBerr[i];
        }
    }
    double sumx = accumulate(x.begin(),x.end(),0.0);
    double avex = sumx/x.size();
    double err_ave = sum/count;
    return make_tuple(avex, err_ave); 
}

int main(int argc, char** argv){
    TApplication app("app", &argc, argv);
    string meas_filedirectory = "/Users/shohtatakami/physics/COMETDS/Scan/";
    string opera_filedirectory = "/Users/shohtatakami/physics/COMETDS/DS189A_944turns/LaserTrackerPosition/root/";
    vector<pair<string,string>> MFMfile_pairs = {
        {"test20240626-6.root","X_0Y_0.root"},
        {"test20240627-0.root","X_200Y_0.root"},
        {"test20240627-1.root","X_500Y_0.root"},
        {"test20240627-2rev.root","X_760Y_0rev.root"},
        {"test20240627-3.root","X_-200Y_0.root"},
        {"test20240627-4.root","X_-500Y_0.root"},
        {"test20240627-6.root","X_-760Y_0.root"},
        {"test20240627-7.root","X_-650Y_0.root"},
        {"test20240627-8.root","X_650Y_0.root"},
        //{"test20240628-0_rev.root","X_0Y_780.csv_Map.csv.root"},//(test20240628-0_revはB値×5/3倍したもの)
        //{"test20240628-2.root","X_0Y_-729.csv_Map.csv.root"},
        //{"test20240710-1.root","X_0Y_200.csv_Map.csv.root"},
        //{"test20240710-4.root","X_0Y_-200.csv_Map.csv.root"}
    };
    
    for(const auto& pair : MFMfile_pairs){
        string meas_file = meas_filedirectory + pair.first;
        string opera_file = opera_filedirectory + pair.second;
        plot(meas_file.c_str(), opera_file.c_str());
    }

    //peakpositionのプロット比較
    /*
    vector<double> x;
    vector<double> pp_meas,pperr_meas,pp_opera,pperr_opera, pp_operaerr, pperr_operaerr,Bmax, Bmaxposi;//_operaerrはOPERAに測定データのエラーをつけた物へのFittingの情報
    for(const auto& pair : MFMfile_pairs){
        string meas_filename = meas_filedirectory + pair.first;
        string opera_filename = opera_filedirectory + pair.second;

        auto[avex,p1_meas,p1_meas_err,p1_opera,p1_opera_err,p1_operaerr, p1_operaerr_fiterr,Bmax_tmp,Bmaxposi_tmp] = peakposi(meas_filename.c_str(),opera_filename.c_str());
        x.push_back(avex);
        pp_meas.push_back(p1_meas);
        pperr_meas.push_back(p1_meas_err);
        pp_opera.push_back(p1_opera);
        pperr_opera.push_back(p1_opera_err);
        pp_operaerr.push_back(p1_operaerr);
        pperr_operaerr.push_back(p1_operaerr_fiterr);
        Bmax.push_back(Bmax_tmp);
        Bmaxposi.push_back(Bmaxposi_tmp);
    }

    TCanvas *cpp = new TCanvas("cpp", "peakposition", 800, 600);
    TCanvas *cppcomp = new TCanvas("cppcomp", "peakposition_comp", 800, 600);
    TGraphErrors *g_measpp = new TGraphErrors(x.size(), &x[0], &pp_meas[0], 0, &pperr_meas[0]);
    TGraphErrors *g_operapp = new TGraphErrors(x.size(), &x[0], &pp_opera[0], 0, &pperr_opera[0]);
    TGraphErrors *g_operaerr = new TGraphErrors(x.size(), &x[0], &pp_operaerr[0], 0, &pperr_operaerr[0]);
    TGraph *g_Bmax = new TGraph(x.size(), &x[0], &Bmaxposi[0]);
    g_measpp->SetName("Measured Data");
    g_operapp->SetTitle("Opera ");
    g_operaerr->SetTitle("Opera with Measured Error");
    g_Bmax->SetTitle("Bmax Position");
    
    g_measpp->SetMarkerColor(kBlue);
    g_measpp->SetMarkerStyle(8);
    g_operapp->SetMarkerColor(kRed+2);
    g_operapp->SetMarkerStyle(8);
    g_operaerr->SetMarkerStyle(8);
    g_operaerr->SetMarkerColor(kGreen);
    g_Bmax->SetMarkerStyle(8);
    g_Bmax->SetMarkerColor(kYellow);
    cpp->cd();
    TMultiGraph *mgpp=new TMultiGraph();
    mgpp->Add(g_measpp);
    mgpp->Add(g_operapp);
    mgpp->SetTitle("Peak Position with FactoryPillar");
    mgpp->GetXaxis()->SetTitle("x (mm)");
    mgpp->GetYaxis()->SetTitle("peakposition (mm)");
    mgpp->Draw("AP");
    cpp->SetGridx();
    cpp->SetGridy();
    cpp->Update();
    //cpp->SaveAs("peakposition.pdf");
    cppcomp->cd();
    TMultiGraph *mg_ppcomp = new TMultiGraph();
    mg_ppcomp->Add(g_operapp);
    if(g_operaerr) mg_ppcomp->Add(g_operaerr);
    if(g_Bmax) mg_ppcomp->Add(g_Bmax);
    mg_ppcomp->SetTitle("Peak Position Comp");
    mg_ppcomp->GetXaxis()->SetTitle("x (mm)");
    mg_ppcomp->GetYaxis()->SetTitle("peak position (mm)");
    mg_ppcomp->Draw("AP");
    cppcomp->SetGridx();
    cppcomp->SetGridy();
    cppcomp->Update();
    //cpp->SaveAs("peakposition.pdf");
    
    


    //Fitting範囲のBのエラーの評価
    vector<double> x, Berr;
    for(const auto& pair : MFMfile_pairs){
        string meas_filename = meas_filedirectory + pair.first;
        string opera_filename = opera_filedirectory + pair.second;

        auto [avex,Berr_val] = err_evaluate(meas_filename.c_str(),opera_filename.c_str());
        x.push_back(avex);
        Berr.push_back(Berr_val);
    }
    TCanvas *cph = new TCanvas("canvas", "peakheight", 800, 600);
    TGraph *g_peakheight = new TGraph(x.size(), &x[0], &Berr[0]);
    g_peakheight->SetMarkerColor(kSpring);
    g_peakheight->SetMarkerStyle(8);

    g_peakheight->SetTitle("Berr btw -150~150");
    g_peakheight->GetXaxis()->SetTitle("x (mm)");
    g_peakheight->GetYaxis()->SetTitle("Berr (T)");
    g_peakheight->Draw("AP");
    
    cph->SetGridx();
    cph->SetGridy();
    cph->Update();
    cph->SaveAs("Berr.pdf");*/
    cerr << "Before app.Run()" << endl;
    app.Run();
    cerr << "After app.Run()" << endl;
    return 0;
}


