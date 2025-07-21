// AngleError.cpp - C++/ROOT version of AngleError.py
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TMath.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TBrowser.h>
#include <TApplication.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;
string Axname; string ErrorOrigin; //string Rotation; 
vector<double> Bxmeas, Bymeas, Bzmeas, Bmeas, Xmeas, Ymeas, Zmeas, Bxstd, Bystd, Bzstd, Bstd;
vector<double> Bxcorr, Bycorr, Bzcorr, Bcorr;//sensorangle corrected
vector<double> Bxplus, Byplus, Bzplus, Bplus;//センサー角度に加えて角度誤差プラスした角度補正
vector<double> Bxdiff, Bydiff, Bzdiff, Bdiff;//センサー角度に加えて角度誤差プラスした角度補正
//string AxIdentify(const string& MEASfile){
//    string name = MEASfile;
//    string Axname;
//    string meas_dir = "/Users/shohtatakami/physics/COMETDS/newScan/";
//    if (name == meas_dir + "test20240626-6.root"){Axname = "X_0Y_0";}
//    else if (name == meas_dir + "test20240627-3.root"){Axname = "X_-200Y_0";}
//    else if (name == meas_dir + "test20240627-4.root"){Axname = "X_-500Y_0";}
//    else if (name == meas_dir + "test20240627-6.root"){Axname = "X_-760Y_0";}
//    else if (name == meas_dir + "test20240627-7.root"){Axname = "X_-650Y_0";}
//    else if (name == meas_dir + "test20240627-0.root"){Axname = "X_200Y_0";}
//    else if (name == meas_dir + "test20240627-1.root"){Axname = "X_500Y_0";}
//    else if (name == meas_dir + "test20240627-8.root"){Axname = "X_650Y_0";}
//    else if (name == meas_dir + "test20240627-2_rev.root"){Axname = "X_760Y_0";}
//    return Axname;
//}
string AxIdentify(const string& MEASfile) {
    string meas_dir = "/Users/shohtatakami/physics/COMETDS/newScan/";
    if (MEASfile == meas_dir + "test20240626-6.root") return "X_0Y_0";
    else if (MEASfile == meas_dir + "test20240627-3.root") return "X_-200Y_0";
    else if (MEASfile == meas_dir + "test20240627-4.root") return "X_-500Y_0";
    else if (MEASfile == meas_dir + "test20240627-6.root") return "X_-760Y_0";
    else if (MEASfile == meas_dir + "test20240627-7.root") return "X_-650Y_0";
    else if (MEASfile == meas_dir + "test20240627-0.root") return "X_200Y_0";
    else if (MEASfile == meas_dir + "test20240627-1.root") return "X_500Y_0";
    else if (MEASfile == meas_dir + "test20240627-8.root") return "X_650Y_0";
    else if (MEASfile == meas_dir + "test20240627-2rev.root") return "X_760Y_0";
    else return "UNKNOWN";
}
// 回転行列
TMatrixD rotM_yaw(double deg) {
    double theta = deg * TMath::DegToRad();
    TMatrixD mat(3,3);
    mat(0,0)=cos(theta); mat(0,1)=0; mat(0,2)=sin(theta);
    mat(1,0)=0;          mat(1,1)=1; mat(1,2)=0;
    mat(2,0)=-sin(theta);mat(2,1)=0; mat(2,2)=cos(theta);
    return mat;
}
//xy平面の回転
TMatrixD rotM_roll(double deg) {
    double theta = deg * TMath::DegToRad();
    TMatrixD mat(3,3);
    mat(0,0)=cos(theta); mat(0,1)=-sin(theta); mat(0,2)=0;
    mat(1,0)=sin(theta); mat(1,1)=cos(theta);  mat(1,2)=0;
    mat(2,0)=0;          mat(2,1)=0;           mat(2,2)=1;
    return mat;
}
//yz平面の回転
TMatrixD rotM_pitch(double deg) {
    double theta = deg * TMath::DegToRad();
    TMatrixD mat(3,3);
    mat(0,0)=1; mat(0,1)=0;          mat(0,2)=0;
    mat(1,0)=0; mat(1,1)=cos(theta); mat(1,2)=-sin(theta);
    mat(2,0)=0; mat(2,1)=sin(theta); mat(2,2)=cos(theta);
    return mat;
}
//take absolute value
double Norm(double x, double y, double z) {
    return sqrt(x*x + y*y + z*z);
}
TMatrixD Matrix_Form(const TMatrixD& v1,const TMatrixD& v2, const TMatrixD& v3){
    for (const auto& v : {&v1, &v2, &v3}){
        if(v->GetNrows() !=3 ||v->GetNcols() !=1){
            throw invalid_argument("All input vectors must be 3×1 TMatrixD Object !!");
        }
    }
    TMatrixD M(3,3);
    for(int i = 0; i<3;++i){
        M(0,i) = v1(i,0);
        M(1,i) = v2(i,0);
        M(2,i)= v3(i,0);
    }
    return M;
}
void InverseCheck(const TMatrixD& M, const TMatrixD& M_inv){
    // 正則か判定
    double det = M.Determinant();
    if (abs(det) < 1e-10) {
        cerr << "Error: Matrix is not invertible (det ≈ 0)." << endl;
        return;
    }
    TMatrixD E = M * M_inv;
    // 単位行列を作成
    TMatrixD I(M.GetNrows(), M.GetNcols());
    I.UnitMatrix();
    // 誤差を評価
    TMatrixD diff = E - I;
    double error = diff.Norm1();//行ごとの絶対値の総和
    // 出力
    //std::cout << "Diff Norm1 = " << error << std::endl;
    if (error < 1e-10) {
        //cout << "✓ Succeeded to check inversed matrix(M × M⁻¹ ≈ I)" << endl;
    } else {
        throw invalid_argument("Inverting Fault");
        //cout << "✗ Failed to check inversed matrix" << endl;
    }
}
void checkMatrixForm(){
    TMatrixD v1(3,1), v2(3,1), v3(3,1);
    v1(0,0) = 2; v1(1,0) = 0; v1(2,0) = 1;
    v2(0,0) = 1; v2(1,0) = 2; v2(2,0) = 2;
    v3(0,0) = 4; v3(1,0) = 0; v3(2,0) = 1;
    TMatrixD test = Matrix_Form(v1,v2,v3);
    TMatrixD test_inv = test;
    test_inv.Invert();
    InverseCheck(test, test_inv);
    test.Print();
}
//sensor angles
//double Bx_yaw = 0.4; double Bx_roll = -0.5;
//double By_pitch = 0.2; double By_roll = 0.2;
//double Bz_yaw = 1.8; double Bz_pitch = 0.4;
double Bx_yaw = 0; double Bx_roll = 0;
double By_pitch = 0; double By_roll = 0;
double Bz_yaw = 0; double Bz_pitch = 0;
TMatrixD ErrorMatrix(const string& origin) {
    TMatrixD mat(3,5); // line : rotation, row : rotation source (cart or sensors)
    mat.Zero();
    if (origin == "Rail") {
        mat(0,0) = 0.15;
        mat(1,0) = 0.10;
        mat(2,0) = 0.10;
    } else if (origin == "BxSensor") {
        mat(0,1) = 0;
        mat(2,1) = 0.1;
    } else if (origin == "BySensor") {
        mat(0,2) = 0.1;
        mat(1,2) = 0.1;
    } else if (origin == "BzSensor") {
        mat(1,3) = 0.1;
        mat(2,3) = 0.1;
    } else if (origin == "Sensors") {
        mat(0,4) = 0.1;
        mat(1,4) = 0.1;
        mat(2,4) = 0.1;
    } else {
        cerr << "Unknown Error Origin: " << origin << endl;
    }
    mat.Print();
    return mat;
}
void DrawTGraphs(){
    //The Configuretion of TGraphs
    gStyle->SetTitleFontSize(0.06);
    gStyle->SetLabelSize(0.05, "X");
    gStyle->SetLabelSize(0.05, "Y");
    gStyle->SetTitleSize(0.05, "X");
    gStyle->SetTitleSize(0.05, "Y");
    gStyle->SetLegendTextSize(0.03); // これは ROOT v6.24以降
    auto c1 = new TCanvas("c1", "c1", 1200, 800);
    c1->Divide(2,2);
    auto c2 = new TCanvas("c2", "c2" ,1200,800);
    c2->Divide(2,2);
    vector<string> components = {"Bx", "By", "Bz", "B"};
    vector<vector<double>*> corrVecs = {&Bxcorr, &Bycorr, &Bzcorr, &Bcorr};
    vector<vector<double>*> plusVecs = {&Bxplus, &Byplus, &Bzplus, &Bplus};
    vector<vector<double>*> stdVecs = {&Bxstd, &Bystd, &Bzstd, &Bstd};
    vector<vector<double>*> diffVecs = {&Bxdiff, &Bydiff, &Bzdiff, &Bdiff};
    for (size_t i = 0; i < components.size(); ++i){
        c1->cd(i+1);
        TGraphErrors *gCorr = new TGraphErrors(Zmeas.size(), &Zmeas[0], &(*corrVecs[i])[0], 0, &(*stdVecs[i])[0]);
        TGraphErrors *gPlus = new TGraphErrors(Zmeas.size(), &Zmeas[0], &(*plusVecs[i])[0], 0, &(*stdVecs[i])[0]);
        
        gCorr->SetMarkerStyle(8);
        gCorr->SetMarkerColor(kAzure+7);
        gCorr->SetLineColor(kAzure+7);
        
        gPlus->SetMarkerStyle(8);
        gPlus->SetMarkerColor(kGreen+1);
        gPlus->SetLineColor(kGreen+1);
        
        TMultiGraph *mg = new TMultiGraph();
        mg->Add(gCorr);
        mg->Add(gPlus);
        //mg->GetYaxis()->SetMaxDigits(3);
        mg->Draw("AP");
        mg->SetTitle((components[i]).c_str());

        TLegend *legend = new TLegend(0.7, 0.1, 0.9, 0.25);
        legend->AddEntry(gCorr, (components[i] + " corr").c_str(), "P");
        legend->AddEntry(gPlus, (components[i] + " plus").c_str(), "P");
        legend->Draw();
        c1->Update();
    }
    c1->Update();
    string outputdir = "/Users/shohtatakami/physics/COMETDS/ErrorBudgetv2/" + ErrorOrigin + "/";
    auto AngleErrorFile = new TFile(Form("%s%s.root", outputdir.c_str(), Axname.c_str()),"RECREATE");
    c1->Write();
    AngleErrorFile->Close();
    std::cout << "Successfully save " << outputdir << " : "<< Axname << std::endl;

    for(size_t i =0; i < components.size(); ++i){
        c2->cd(i+1);
        TGraph *gDelta = new TGraph(Zmeas.size(), &Zmeas[0], &(*diffVecs[i])[0]);
        gDelta->SetMarkerStyle(8);
        gDelta->SetMarkerColor(kOrange+10);
        gDelta->SetLineColor(kOrange+10);
        gDelta->GetYaxis()->SetMaxDigits(3);
        gDelta->SetTitle((components[i] + "Diff").c_str());
        gDelta->Draw("AP");
        c2->Update();
    }
    c2->Update();
    auto AngleErrorDiffFile = new TFile(Form("%s%s_diff.root", outputdir.c_str(), Axname.c_str()), "RECREATE");
    c2->Write();
    AngleErrorDiffFile->Close();
    std::cout << "Successfully save " << outputdir << " : "<< Axname << "Diff" << std::endl;    
}

void Make_Vector(const string& MEASfile, const string& Error_Origin) {
    Bxmeas.clear(); Bymeas.clear(); Bzmeas.clear(); Bmeas.clear();
    Bxstd.clear(); Bystd.clear(); Bzstd.clear(); Bstd.clear();
    Xmeas.clear(); Ymeas.clear(); Zmeas.clear();
    Bxcorr.clear(); Bycorr.clear(); Bzcorr.clear(); Bcorr.clear();
    Bxplus.clear(); Byplus.clear(); Bzplus.clear(); Bplus.clear();
    Bxdiff.clear(); Bydiff.clear(); Bzdiff.clear(); Bdiff.clear();
    //cout << "Processing: " << MEASfile << " with ErrorOrigin = " << Error_Origin << ", Rotation = " << Rotation << endl;
    Axname = AxIdentify(MEASfile);
    ErrorOrigin = Error_Origin;
    //Rotation = rotation;

    cout << "AxisName : " << Axname << endl;

    TMatrixD AngleErrors = ErrorMatrix(Error_Origin);
    TFile* file1 = TFile::Open(MEASfile.c_str(), "READ");

    if (!file1 || file1->IsZombie()) {
        cerr << "Error opening file: " << MEASfile << endl;
        return;
    }

    TTree* meastree = (TTree*)file1->Get("tree");
    if (!meastree) {
        cerr << "Tree 'tree' not found in file: " << MEASfile << endl;
        return;
    }

    // Set up branches
    Double_t Xmean, Ymean, Zmean, Xstdev, Ystdev, Zstdev, X, Y, Z;
    meastree->SetBranchAddress("Xmean", &Xmean);
    meastree->SetBranchAddress("Ymean", &Ymean);
    meastree->SetBranchAddress("Zmean", &Zmean);
    meastree->SetBranchAddress("Xstdev", &Xstdev);
    meastree->SetBranchAddress("Ystdev", &Ystdev);
    meastree->SetBranchAddress("Zstdev", &Zstdev);
    meastree->SetBranchAddress("X", &X);
    meastree->SetBranchAddress("Y", &Y);
    meastree->SetBranchAddress("Z", &Z);

    Long64_t nEntries = meastree->GetEntries();
    for(Long64_t i = 0; i < nEntries; ++i){
        meastree->GetEntry(i);
        double Bx = Xmean * (-3.0/ 5.0);
        double By = Ymean * (3.0 / 5.0);
        double Bz = Zmean * (3.0 / 5.0);
        double B = sqrt(pow(Bx,2) + pow(By,2) + pow(Bz,2));
        double Bxdev = Xstdev * (3.0/5.0);
        double Bydev = Ystdev * (3.0/5.0);
        double Bzdev = Zstdev * (3.0/5.0);
        double Bdev = sqrt(pow(Bx/B,2)*pow(Bxdev ,2) + pow(By/B,2)*pow(Bzdev ,2) + pow(Bz/B,2)*pow(Bzdev ,2));
        Bxmeas.push_back(Bx); Bymeas.push_back(By); Bzmeas.push_back(Bz);
        Bxstd.push_back(Bxdev); Bystd.push_back(Bydev); Bzstd.push_back(Bzdev); Bstd.push_back(Bdev);
        Bmeas.push_back(Norm(Bx, By, Bz));
        Xmeas.push_back(X); Ymeas.push_back(Y); Zmeas.push_back(Z);
    }
    //calculation of correction(Only Sensor Angle)
    for(Long64_t i = 0; i < nEntries; ++i){
        TMatrixD B_vec(3,1);
        B_vec(0,0) = Bxmeas[i];
        B_vec(1,0) = Bymeas[i];
        B_vec(2,0) = Bzmeas[i];
        //correction matrix
        TMatrixD nx(3,1); nx(0,0)=1; nx(1,0)=0; nx(2,0)=0;
        TMatrixD ny(3,1); ny(0,0)=0; ny(1,0)=1; ny(2,0)=0;
        TMatrixD nz(3,1); nz(0,0)=0; nz(0,0)=0; nz(2,0)=1;
        TMatrixD nBx = rotM_yaw(Bx_yaw) * rotM_roll(Bx_roll) * nx;
        TMatrixD nBy = rotM_pitch(By_pitch) * rotM_roll(By_roll) * ny;
        TMatrixD nBz = rotM_pitch(Bz_pitch) * rotM_yaw(Bz_yaw) * nz;
        TMatrixD M = Matrix_Form(nBx,nBy,nBz);
        try{TMatrixD M = Matrix_Form(nBx, nBy, nBz);
        }catch(const invalid_argument& error){
            cerr<<"Error:"<<error.what()<<endl;//入力が3×1のベクトルか確認
        }
        TMatrixD M_inv = M;
        M_inv.Invert();
        try{InverseCheck(M, M_inv);
        }catch(const invalid_argument& error){
            cerr<<"Error:"<<error.what()<<endl;//逆行列を計算できているか確認
        }
        TMatrixD Bcorr_vec = M.Invert() * B_vec;
        Bxcorr.push_back(Bcorr_vec(0,0)); //Bxcorr[i]=でアクセスするときサイズが確保されてないと未定義動作Seg faultになる
        Bycorr.push_back(Bcorr_vec(1,0));
        Bzcorr.push_back(Bcorr_vec(2,0));
        Bcorr.push_back(Norm(Bcorr_vec(0,0), Bcorr_vec(1,0), Bcorr_vec(2,0)));
    }
    //calculation of correction(Only Sensor Angle + δθ)
    TMatrixD Error = ErrorMatrix(Error_Origin);
    for(Long64_t i = 0; i < nEntries; ++i){
        TMatrixD B_vec(3,1);
        B_vec(0,0) = Bxmeas[i];
        B_vec(1,0) = Bymeas[i];
        B_vec(2,0) = Bzmeas[i];
        //correction matrix
        TMatrixD nx(3,1); nx(0,0)=1; nx(1,0)=0; nx(2,0)=0;
        TMatrixD ny(3,1); ny(0,0)=0; ny(1,0)=1; ny(2,0)=0;
        TMatrixD nz(3,1); nz(0,0)=0; nz(0,0)=0; nz(2,0)=1;
        TMatrixD nBx = rotM_roll(Error(0,0)) * rotM_yaw(Error(2,0)) * // Rail Error
                        rotM_roll(Error(0,1)) * rotM_yaw(Error(2,1)) * // sensor Error
                        rotM_roll(Bx_roll) * rotM_yaw(Bx_yaw) * nx; //sensor tilt angle correction
        
        TMatrixD nBy =  rotM_roll(Error(0,0)) * rotM_pitch(Error(1,0)) * // Rail Error
                        rotM_roll(Error(0,2)) * rotM_pitch(Error(1,2)) * //sensor Error
                        rotM_roll(By_roll) * rotM_pitch(By_pitch) * ny;
        
        TMatrixD nBz = rotM_pitch(Error(1,0)) * rotM_yaw(Error(2,0)) * //Rail
                        rotM_pitch(Error(1,3)) * rotM_yaw(Error(2,3)) * //sensor Error
                        rotM_pitch(Bz_pitch) * rotM_yaw(Bz_yaw) * nz;

        TMatrixD M = Matrix_Form(nBx,nBy,nBz);
        try{TMatrixD M = Matrix_Form(nBx, nBy, nBz);
        }catch(const invalid_argument& error){
            cerr<<"Error:"<<error.what()<<endl;
        }
        TMatrixD M_inv = M;
        M_inv.Invert();
        try{InverseCheck(M, M_inv);
        }catch(const invalid_argument& error){
            cerr<<"Error:"<<error.what()<<endl;
        }
        TMatrixD Bcorr_vec = M.Invert() * B_vec;
        Bxplus.push_back(Bcorr_vec(0,0));
        Byplus.push_back(Bcorr_vec(1,0));
        Bzplus.push_back(Bcorr_vec(2,0));
        Bplus.push_back(Norm(Bcorr_vec(0,0), Bcorr_vec(1,0), Bcorr_vec(2,0)));
    }

    //make diff vector
    for(int i = 0; i < nEntries; ++i){
        Bxdiff.push_back(Bxplus[i] - Bxcorr[i]);
        Bydiff.push_back(Byplus[i] - Bycorr[i]);
        Bzdiff.push_back(Bzplus[i] - Bzcorr[i]);
        Bdiff.push_back(Bplus[i] - Bcorr[i]);
    }
    DrawTGraphs();
}

int main(int argc, char** argv){
    TApplication app("App", &argc, argv);  
    vector<string> measfiles = {
        "test20240627-6.root",
        "test20240626-6.root",
        "test20240627-2rev.root"
    };

    string opera_dir = "/Users/shohtatakami/physics/COMETDS/DS189A_944turns_Integration/DS189A_944turns_FieldIntegration_LaserTrackerPos/root/";
    string meas_dir = "/Users/shohtatakami/physics/COMETDS/newScan/";

    for (const auto& measfile : measfiles) {
        Make_Vector(meas_dir + measfile, "Rail");
        Make_Vector(meas_dir + measfile, "BxSensor");
        Make_Vector(meas_dir + measfile, "BySensor");
        Make_Vector(meas_dir + measfile, "BzSensor");
    }
    new TBrowser("BigBrowser");
    app.Run();
    return 0;
}
