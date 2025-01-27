#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <TCanvas.h>
#include <TF1.h>        
#include <TGraphErrors.h>
#include <TGraph2D.h>
#include <TMultiGraph.h>
#include <TPolyLine3D.h>
#include <TApplication.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TStyle.h>
using namespace std;
/*_________________________________________
レーザートラッカーによる三次元座標と視点からの移動距離を示すencoderの値から
測定点における三次元座標を線形補間にて算出するプログラム
*/
//zの正方向はレーザートラッカーに合わせる（ビーム下流->上流が正）


//レーザートラッカーの座標の値
vector <double> x;
vector <double> y;
vector <double> y_temp;
vector <double> z_temp;//z軸調整前
vector <double> z;//z軸調整後
vector <double> l_laser;//レーザートラッカー上での移動距離


//エンコーダーのカウント
vector <double> cnt_z; //encoder count
vector <double> encoder_z; //encoder (mm)
vector <double> l_encoder; //encoder (mm)
vector <double> delta;//encoder(mm) - z(mm)



void Intovectors(const char* filename_encoder = "test20240624-2"){
    //reading data of Laser Tracker
    //lasertracker_rawdataの説明
    //用いる軸に応じてstarlineとendlineの数字を変える
    /*
    X=0,Y=0 204行目から415行目
    X=+200,Y=0 452行目から461行目
    X=+500,Y=0 463行目から472行目
    X=+650,Y=0 473行目から482行目
    X=+760,Y=0 483行目から492行目 //491行目まで10点目おかしい
    X=-200,Y=0 493行目から502行目
    X=-500,Y=0 503行目から512行目
    X=-650,Y=0 513行目から522行目
    X=-760,Y=0 523行目から532行目
    X=0,Y=-200 533行目から542行目
    X=0,Y=-500 558行目から567行目
    X=0,Y=-730 576行目から585行目 
    X=0,Y=+200 590行目から599行目
    X=0,Y=+780 658行目から667行目*/
    string filename_laser = "lasertracker_rawdata.csv"; 
    ifstream file1(filename_laser);
    string line_laser;
    int startline_laser = 658;
    int endline_laser = 667;
    int currentline_laser = 1; //1行目は1

    if (!file1.is_open()){
        cerr << "failed to open the file" << filename_laser <<endl;   
    }

    while(getline(file1,line_laser)){
        if(currentline_laser > endline_laser ){
            break;
        }

        if(currentline_laser >= startline_laser && currentline_laser <= endline_laser){
            stringstream ss(line_laser);
            string item_laser;
            vector <string> items_laser;

            while(getline(ss, item_laser,',')){
                items_laser.push_back(item_laser);
            }
            x.push_back(stod(items_laser[4]));
            y_temp.push_back(stod(items_laser[5]));
            z_temp.push_back(stod(items_laser[6]));
        }
        
        //cout<< "currentline_laser: " <<currentline_laser << "xのサイズ"<< x.size() <<endl;
        currentline_laser = currentline_laser +1;

    }

    // for(int i =0; i<x.size(); i++){
    //     cout << i << "番目: " << x.at(i) << endl;
    // }
    
    file1.close();
    for(int i = 0; i < z_temp.size(); ++i){
        //double z_value = - z_temp[i] + 1660 + 22 - 7.3;
        double z_value = z_temp[i] - (22 - 7.3);
        z.push_back(z_value);
    }
    for(int i = 0; i < y_temp.size(); ++i){
        double y_value = y_temp[i]- 35;
        //cout << z_value << endl;
        y.push_back(y_value);
    }
    //x補正なし

    //reading encoder's position
    //string filename_encoder = "../MFMJune/Scan/test20240624-2";
    ifstream file2(filename_encoder);
    string line_encoder;
    int startline_encoder = 2;
    int endline_encoder = 214;
    int currentline_encoder = 1 ;//1行目は1

    if (!file2.is_open()){
        cerr << "failed to open the file" << filename_encoder <<endl;   
    }

    while(getline(file2,line_encoder)){

        if(currentline_encoder > endline_encoder ){
            break;
        }

        //201行目はレーザートラッカーが抜けているためスキップ
        if(currentline_encoder == 201){
            currentline_encoder = currentline_encoder + 1;
            continue;
        }

        if(currentline_encoder >= startline_encoder && currentline_encoder <= endline_encoder){
        stringstream ss(line_encoder);
        string item_encoder;
        vector <string> items_encoder;

            while(getline(ss, item_encoder,'\t')){
                items_encoder.push_back(item_encoder);
            }

            cnt_z.push_back(stod(items_encoder[2]));
        }
        
        //cout<< "currentline_laser: " <<currentline_laser << "xのサイズ"<< x.size() <<endl;
        currentline_encoder = currentline_encoder +1;

    }
    file2.close();
    for(int i = 0; i < cnt_z.size(); ++i){
        //double encoder_value = cnt_z[i]*(14.0/4000.0) - 7.3;
        double encoder_value = cnt_z[i]*(14.0/4000.0);
        encoder_z.push_back(encoder_value);
    }
   

}

void Interpolation(const char* filename = "test20240626-6"){
    //磁場測定したエンコーダーカウント
    vector <double> cnt_zact; //encoder count
    vector <double> encoder_zact; //enocodercount(mm)

    //線形補間パラメーター
    vector <double> k_x;
    vector <double> intercept_x;
    vector <double> x_rev;
    vector <double> k_y;
    vector <double> intercept_y;
    vector <double> y_rev;
    vector <double> k_z;
    vector <double> intercept_z;
    vector <double> z_rev;

    //reading encoder's position for magnetic field measurement
    //string filename_encoderzact = "../MFMJune/Scan/test20240626-6";
    ifstream file3(filename);
    string line_encoderzact;
    int startline_encoderzact = 2;
    int endline_encoderzact = 214;
    int currentline_encoderzact = 1 ;//1行目は1

    if (!file3.is_open()){
        cerr << "failed to open the file" << filename <<endl;   
    }

    while(getline(file3,line_encoderzact)){

        if(currentline_encoderzact > endline_encoderzact ){
            break;
        }

        if(currentline_encoderzact >= startline_encoderzact && currentline_encoderzact <= endline_encoderzact){
            stringstream ss(line_encoderzact);
            string item_encoderzact;
            vector <string> items_encoderzact;
                while(getline(ss, item_encoderzact, '\t')){
                    items_encoderzact.push_back(item_encoderzact);
                }

            cnt_zact.push_back(stod(items_encoderzact[2]));
        }
        
    //cout<< "currentline_laser: " <<currentline_laser << "xのサイズ"<< x.size() <<endl;
    currentline_encoderzact = currentline_encoderzact +1;
    }
    file3.close();
    for(int i = 0; i < cnt_zact.size(); ++i){
        //double encoder_value = cnt_zact[i]*(14.0/4000.0) - 7.3;
        double encoder_value = cnt_zact[i]*(14.0/4000.0);
        encoder_zact.push_back(encoder_value);
    }

    // linear interpolation of x
    for(int i = 0; i < encoder_z.size() - 1; ++i){
        double k_val = (x[i+1] - x[i])/(encoder_z[i+1] - encoder_z[i]);
        k_x.push_back(k_val);
        double intercept_val = x[i]-(((x[i+1] - x[i])/(encoder_z[i+1] - encoder_z[i])) * encoder_z[i]) ;
        intercept_x.push_back(intercept_val);
        }

    for(int i = 0; i < encoder_zact.size(); ++i){
        //cout<< "磁場測定エンコーダ" << encoder_zact[i] <<endl;
        for(int j = 0; j < encoder_z.size(); ++j){
            if(encoder_z[j] <= encoder_zact[i] && encoder_z[j+1] > encoder_zact[i]){
                double x_revtmp = k_x[j] * encoder_zact[i] + intercept_x[j];
                x_rev.push_back(x_revtmp);    
                //cout<< j << k[j] <<":"<< intercept[j] << endl;
            }
        }
    }
    // linear interpolation of y
    for(int i = 0; i < encoder_z.size() - 1; ++i){
        double k_val = (y[i+1] - y[i])/(encoder_z[i+1] - encoder_z[i]);
        k_y.push_back(k_val);
        double intercept_val = y[i]-(((y[i+1] - y[i])/(encoder_z[i+1] - encoder_z[i])) * encoder_z[i]) ;
        intercept_y.push_back(intercept_val);
        }

    for(int i = 0; i < encoder_zact.size(); ++i){
        //cout<< "磁場測定エンコーダ" << encoder_zact[i] <<endl;
        for(int j = 0; j < encoder_z.size(); ++j){
            if(encoder_z[j] <= encoder_zact[i] && encoder_z[j+1] > encoder_zact[i]){
                double y_revtmp = k_y[j] * encoder_zact[i] + intercept_y[j];
                y_rev.push_back(y_revtmp);    
                //cout<< j << k[j] <<":"<< intercept[j] << endl;
            }
        }
    }
    // linear interpolation of z
    for(int i = 0; i < encoder_z.size() - 1; ++i){
        double k_val = (z[i+1] - z[i])/(encoder_z[i+1] - encoder_z[i]);
        k_z.push_back(k_val);
        double intercept_val = z[i]-(((z[i+1] - z[i])/(encoder_z[i+1] - encoder_z[i])) * encoder_z[i]);
        intercept_z.push_back(intercept_val);
        }

    for(int i = 0; i < encoder_zact.size(); ++i){
        //cout<< "磁場測定エンコーダ" << encoder_zact[i] <<endl;
        for(int j = 0; j < encoder_z.size(); ++j){
            if(encoder_z[j] <= encoder_zact[i] && encoder_z[j+1] > encoder_zact[i]){
                double z_revtmp = k_z[j] * encoder_zact[i] + intercept_z[j];
                z_rev.push_back(z_revtmp);    
                //cout<< j << k[j] <<":"<< intercept[j] << endl;
            }
        }
    }
   /* for(int i = 0; i < encoder_zact.size(); ++i){
        cout << "Non-Interpolated (mm)" << encoder_zact[i]<<" : " << "Interpolated(mm)" << z_rev[i] << endl;
    }*/
    for(int i = 0; i < encoder_zact.size(); ++i){
        cout<< "enocoder"<<encoder_zact[i]<<" : xの位置"<< x_rev[i] << ": yの位置"<<y_rev[i]<<": zの位置"<<z_rev[i]<<endl;
    }
    cout<<"enocder "<<encoder_zact[212]<<" : x "<<k_x[210]*encoder_zact[212]+intercept_x[210]<<" : y "<<k_y[210]*encoder_zact[212]+intercept_y[210]<<" : z "<<k_z[210]*encoder_zact[212]+intercept_z[210]<<endl;
    //making csv file
    string outputFileName;
    if (strcmp(filename, "Scan/test20240626-6") == 0) {
        outputFileName = "X_0Y_0.csv";
    } else if (strcmp(filename, "Scan/test20240627-0") == 0) {
        outputFileName = "X_200Y_0.csv";
    } else if (strcmp(filename, "Scan/test20240627-1") == 0) {
        outputFileName = "X_500Y_0.csv";
    } else if (strcmp(filename, "Scan/test20240627-2rev") == 0) {
        outputFileName = "X_760Y_0.csv";
    } else if (strcmp(filename, "Scan/test20240627-3") == 0) {
        outputFileName = "X_-200Y_0.csv";
    } else if (strcmp(filename, "Scan/test20240627-4") == 0) {
        outputFileName = "X_-500Y_0.csv";
    } else if (strcmp(filename, "Scan/test20240627-6") == 0) {
        outputFileName = "X_-760Y_0.csv";
    } else if (strcmp(filename, "Scan/test20240627-7") == 0) {
        outputFileName = "X_-650Y_0.csv";
    } else if (strcmp(filename, "Scan/test20240627-8") == 0) {
        outputFileName = "X_650Y_0.csv";
    } else if (strcmp(filename, "Scan/test20240628-0_rev") == 0) {
        outputFileName = "X_0Y_780.csv";
    } else if (strcmp(filename, "Scan/test20240628-2") == 0) {
        outputFileName = "X_0Y_-729.csv";
    } else if (strcmp(filename, "Scan/test20240710-1") == 0) {
        outputFileName = "X_0Y_200.csv";
    } else if (strcmp(filename, "Scan/test20240710-4") == 0) {
        outputFileName = "X_0Y_-200.csv";
    } else {
        cerr << "Unsupported filename: " << filename << endl;
        return; // サポートされていないファイル名の場合、処理を終了
    }

    // 出力ファイルを開く
    ofstream outputFile(outputFileName);
    if (!outputFile.is_open()) {
        cerr << "Failed to open file: " << outputFileName << endl;
        return;
    }
    
    // header
    outputFile<<"Encoder,x,y,z"<<endl;
    for(int i = 0;i < encoder_zact.size(); ++i){
        outputFile <<encoder_zact[i]<<","<<x_rev[i]<<","<<y_rev[i]<<","<<z_rev[i]<<endl;
    }
    outputFile.close();

    //encoder vs lasertracker & encoder vs real coordinate by Interpolation
    TCanvas *c1 = new TCanvas( "c1" , "Graphs" , 1200 , 800 );
    TGraph *g_x = new TGraph(encoder_z.size(), &encoder_z[0], &x[0]);
    TGraph *g_y = new TGraph(encoder_z.size(), &encoder_z[0], &y[0]);
    TGraph *g_z = new TGraph(encoder_z.size(), &encoder_z[0], &z[0]);
    TGraph *g_xact = new TGraph(encoder_zact.size(), &encoder_zact[0], &x_rev[0]);
    TGraph *g_yact = new TGraph(encoder_zact.size(), &encoder_zact[0], &y_rev[0]);
    TGraph *g_zact = new TGraph(encoder_zact.size(), &encoder_zact[0], &z_rev[0]);
    c1->Divide(1,3);

    g_x->SetMarkerStyle(20);
    g_x->SetMarkerColor(kRed+2);
    g_xact->SetMarkerStyle(6);
    g_xact->SetMarkerColor(kBlue);
    g_y->SetMarkerStyle(20);
    g_y->SetMarkerColor(kRed+2);
    g_yact->SetMarkerStyle(6);
    g_yact->SetMarkerColor(kBlue);
    g_z->SetMarkerStyle(20);
    g_z->SetMarkerColor(kRed+2);
    g_zact->SetMarkerStyle(6);
    g_zact->SetMarkerColor(kBlue);

    TLegend *legendx = new TLegend(0.9, 0.1, 1.0, 0.2); 
    TLegend *legendy = new TLegend(0.9, 0.1, 1.0, 0.2); 
    TLegend *legendz = new TLegend(0.9, 0.1, 1.0, 0.2); 
    legendx->AddEntry(g_x, "Laser Tracker", "p");  
    legendx->AddEntry(g_xact, "Measured", "p");  
    legendy->AddEntry(g_y, "Laser Tracker", "p");  
    legendy->AddEntry(g_yact, "Measured", "p"); 
    legendz->AddEntry(g_z, "Laser Tracker", "p");  
    legendz->AddEntry(g_zact, "Measured", "p"); 
    TMultiGraph *mg_x = new TMultiGraph();
    TMultiGraph *mg_y = new TMultiGraph();
    TMultiGraph *mg_z = new TMultiGraph();
    mg_x->Add(g_x);
    mg_x->Add(g_xact);
    mg_y->Add(g_y);
    mg_y->Add(g_yact);
    mg_z->Add(g_z);
    mg_z->Add(g_zact);

    c1->cd(1);
    mg_x->Draw("AP");
    mg_x->GetXaxis()->SetTitle("encoder");
    mg_x->GetXaxis()->SetTitleSize(0.06);
    mg_x->GetXaxis()->SetLabelSize(0.06);
    mg_x->GetYaxis()->SetTitle("x (mm)"); 
    mg_x->GetYaxis()->SetTitleSize(0.06);
    mg_x->GetYaxis()->SetLabelSize(0.06);
    gPad->SetBottomMargin(0.15);
    legendx->Draw();  // 凡例をカンバスに描画
    c1->cd(2);
    mg_y->Draw("AP");
    mg_y->GetXaxis()->SetTitle("encoder");
    mg_y->GetXaxis()->SetTitleSize(0.06);
    mg_y->GetXaxis()->SetLabelSize(0.06);
    mg_y->GetYaxis()->SetTitle("y (mm)"); 
    mg_y->GetYaxis()->SetTitleSize(0.06);
    mg_y->GetYaxis()->SetLabelSize(0.06);
    gPad->SetBottomMargin(0.15);
    legendy->Draw();  // 凡例をカンバスに描画
    c1->cd(3);
    mg_z->Draw("AP");
    mg_z->GetXaxis()->SetTitle("encoder");
    mg_z->GetXaxis()->SetTitleSize(0.06);
    mg_z->GetXaxis()->SetLabelSize(0.06);
    mg_z->GetYaxis()->SetTitle("z (mm)"); 
    mg_z->GetYaxis()->SetTitleSize(0.06);
    mg_z->GetYaxis()->SetLabelSize(0.06);
    gPad->SetBottomMargin(0.15);
    legendz->Draw();  // 凡例をカンバスに描画

    
    string filenameStr = outputFileName;
    filenameStr.erase(0,filenameStr.find_last_of('/') + 1);
    string directory = "InterpolatedPlot_distiguishing/";
    string Filename = directory + filenameStr +"_xyz.pdf";
    c1->SaveAs(Filename.c_str());


    
}
//////////////////

void MFmap(const char* measure_filename = "test20240626-6",const char* OPERAfilename = "Interpolated_OPERA/X_0&Y_0rev.csv_Map.csv"){
    vector<double> cnt; //value of encoder count
    vector<double> xmean; //voltage value corresponding to magnetic field
    vector<double> B_x;//the value of B_x calibllated by the value of x range
    vector<double> ymean; 
    vector<double> B_y;
    vector<double> zmean;
    vector<double> B_z; 
    vector<double> Stdxmean;//vector of standard deviation of B_x
    vector<double> StdB_x;
    vector<double> Stdymean;
    vector<double> StdB_y;
    vector<double> Stdzmean;
    vector<double> StdB_z;
    vector<double> xrange;
    vector<double> yrange;
    vector<double> zrange;
    vector<double> B;
    vector<double> StdB;

    ifstream file3(measure_filename);
    string line_measure;
    int startline_measure = 2;
    int endline_measure = 214;
    int currentline_measure = 1 ;//1行目は1

    if (!file3.is_open()){
        cerr << "failed to open the file" << measure_filename <<endl;   
    }

    while(getline(file3,line_measure)){

        if(currentline_measure> endline_measure){
            break;
        }

        if(currentline_measure >= startline_measure && currentline_measure <= endline_measure){
            stringstream ss(line_measure);
            string item_measure;
            vector <string> items_measure;
            while(getline(ss, item_measure, '\t')){
                    items_measure.push_back(item_measure);
            }
            xmean.push_back(std::stod(items_measure[3]));
            Stdxmean.push_back(std::stod(items_measure[4]));
            ymean.push_back(std::stod(items_measure[5]));
            Stdymean.push_back(std::stod(items_measure[6]));
            zmean.push_back(std::stod(items_measure[7]));
            Stdzmean.push_back(std::stod(items_measure[8]));
            xrange.push_back(std::stod(items_measure[11]));
            yrange.push_back(std::stod(items_measure[12]));
            zrange.push_back(std::stod(items_measure[13]));
        }
    //cout<< "currentline_laser: " <<currentline_laser << "xのサイズ"<< x.size() <<endl;
    currentline_measure = currentline_measure +1;
    }
    file3.close();

    for(int i = 0; i<xmean.size(); ++i){
        double revised_x =xmean[i]*(xrange[i]/5);
        B_x.push_back(revised_x);
    }

    for(int i = 0; i<ymean.size(); ++i){
        double revised_y =ymean[i]*(yrange[i]/5);
        B_y.push_back(revised_y);
    }

    for(int i = 0; i<zmean.size(); ++i){
        double revised_z =zmean[i]*(zrange[i]/5);
        B_z.push_back(revised_z);
    }

    for(int i = 0; i<cnt.size(); ++i){
        double revised_Stdx =Stdxmean[i]*(xrange[i]/5);
        StdB_x.push_back(revised_Stdx);
    }

    for(int i = 0; i<cnt.size(); ++i){
        double revised_Stdy =Stdymean[i]*(yrange[i]/5);
        StdB_y.push_back(revised_Stdy);
    }   

    for(int i = 0; i<cnt.size(); ++i){
        double revised_Stdz =Stdzmean[i]*(zrange[i]/5);
        StdB_z.push_back(revised_Stdz);
    }
    
    //calculating the Amplitude of ,agnetic field
    for(int i=0; i<xmean.size(); ++i){
        double B_value = sqrt(pow(B_x[i],2)+pow(B_y[i],2)+pow(B_z[i],2));
        B.push_back(B_value);
        }

    for(int i =0; i<xmean.size(); ++i){
        double StdB_value = sqrt(pow(B_x[i]/B[i] ,2)*pow(StdB_x[i],2)+pow(B_y[i]/B[i] ,2)*pow(StdB_y[i],2)+pow(B_z[i]/B[i] ,2)*pow(StdB_z[i],2));
        StdB.push_back(StdB_value);
    }

    //reading simulation file by OPERA 3D    
    ifstream file4(OPERAfilename);
    string line_OPERA;
    if (!file4.is_open()){
        cerr << "failed to open the file" << measure_filename <<endl;   
    }
    while(getline(file4, line_OPERA )){
        stringstream ss(line_OPERA);// initialize line as stringstream
        string item_OPERA;
        vector<string> items_OPERA;

        //split by conmma
        while(getline(ss, item_OPERA ,',')){//read data from ss and store it in the item
            //cout<<item<<endl;
            items_OPERA.push_back(item_OPERA);
        }

        B_x.push_back(stod(items_OPERA[4]));
        B_y.push_back(stod(items_OPERA[5]));
        B_z.push_back(stod(items_OPERA[6]));
        B.push_back(stod(items_OPERA[7]));
    }
    file4.close();

    // TGraph2Dオブジェクトを作成 (磁場値は色で表現)
    TGraph2D *gmap = new TGraph2D(x.size(), &x[0], &y[0], &z[0]);
    TCanvas *c1 = new TCanvas("c1", "3D Magnetic Field Map ", 800, 600);
    gmap->SetTitle("3D Magnetic Field Map with Field Errors; X; Y; Z");
    gmap->SetMarkerStyle(20);
    gmap->SetMarkerSize(1);
    gmap->Draw("PCOLZ");

    // 磁場の誤差をエラーバーとして手動で描画
    /*for (int i = 0; i < x.size() ; i++) {
        // 磁場のエラーバーを z 軸方向に表示する (z ± B_err)
        TPolyLine3D *line = new TPolyLine3D(2);
        line->SetPoint(0, x[i], y[i], z[i] - B_err[i]);  // 下側のエラーバーの端点
        line->SetPoint(1, x[i], y[i], z[i] + B_err[i]);  // 上側のエラーバーの端点
        line->SetLineColor(kRed);  // エラーバーの色を赤に設定
        line->Draw();  // エラーバーを描画
    }*/

    // カンバスを更新して表示
    //c1->Update();
}




void PositionCorrecter(){
    //lasertracker上での移動距離計算
    for(int i = 0; i < z_temp.size() -1 ; ++i){
        double l_value = sqrt(pow(x[i+1] - x[i], 2) + pow(y[i+1] - y[i], 2) + pow(z[i+1] - z[i], 2));
        l_laser.push_back(l_value);
    }


    for(int i = 0; i < encoder_z.size(); ++i){
        double delta_value = encoder_z[i] - z[i];
        delta.push_back(delta_value);
    }
    
    for(int i = 0; i < z_temp.size() -1 ; ++i){
        double l_value = encoder_z[i+1] - encoder_z[i];
        l_encoder.push_back(l_value);
    }

    // エンコーダーのピッチとレーザートラッカー上での移動距離計算

    for (int i = 0; i < l_laser.size(); ++i){
        std::cout<< "encoder pitch :" << l_encoder[i] <<" lasertracker pitch :" << l_laser[i] <<endl;
    }

    double ave_laserpitch = reduce(l_laser.begin(),l_laser.end())/l_laser.size();
    double dev_laserpitch = 0.0;
    for(int i =0; i < l_laser.size(); ++i){
        dev_laserpitch = dev_laserpitch +pow(l_laser[i] - ave_laserpitch, 2)/l_laser.size();
    }
    
    
    double ave_encoderpitch = reduce(l_encoder.begin(),l_encoder.end())/l_encoder.size();
    double dev_encoderpitch = 0.0;
    for(int i =0; i < l_encoder.size(); ++i){
        dev_encoderpitch = dev_encoderpitch + pow(l_encoder[i] - ave_encoderpitch, 2)/l_encoder.size();
    }

    cout<< "エンコーダーのピッチの平均"<<ave_encoderpitch<<"エンコーダーのピッチのCoefficient of Variation"<<(sqrt(dev_encoderpitch)/ave_encoderpitch)*100<<"%"<<endl;
    cout<< "レーザートラッカーのピッチの平均"<<ave_laserpitch<<"レーザートラッカのCoefficient of Variation"<<(sqrt(dev_laserpitch)/ave_laserpitch)*100<<"%"<<endl;


    
    TCanvas *c1 = new TCanvas( "c1" , "Graphs" , 1200 , 800 );
    //TF1 *firstfit = new TF1("firstfit", "pol1", -2000, 6000);
    TGraph *g_z = new TGraph(z.size(), &encoder_z[0], &z[0]);
    //firstfit->SetLineColor(kGreen);
    //g_center->Fit(firstfit);
    TGraph *g_x = new TGraph(encoder_z.size(), &encoder_z[0], &x[0]);
    TGraph *g_y = new TGraph(encoder_z.size(), &encoder_z[0], &y[0]);
    c1->Divide(1,3);

    c1->cd(3);
    g_z->SetTitle(" encoder : lasertracker_z ");
    g_z->GetXaxis()->SetTitle("encoder (mm)");
    g_z->GetYaxis()->SetTitle("lasertracker_z (mm)");
    g_z->SetMarkerStyle(6);
    g_z->SetMarkerColor(kRed+2);
    g_z->Draw("AP");

    c1->cd(1);
    g_x->SetTitle(" encoder : lasertracker_x ");
    g_x->GetXaxis()->SetTitle("encoder (mm)");
    g_x->GetYaxis()->SetTitle("lasertracker_x (mm)");
    g_x->SetMarkerStyle(6);
    g_x->SetMarkerColor(kRed+2);
    g_x->Draw("AP");

    c1->cd(2);
    g_y->SetTitle(" encoder : lasertracker_y ");
    g_y->GetXaxis()->SetTitle("encoder (mm)");
    g_y->GetYaxis()->SetTitle("lasertracker_y (mm)");
    g_y->SetMarkerStyle(6);
    g_y->SetMarkerColor(kRed+2);
    g_y->Draw("AP");

    
    
    
    /*g_delta->SetTitle(" encoder_z - lasertracker_z ");
    g_delta->GetXaxis()->SetTitle("encoder (mm)");
    g_delta->GetYaxis()->SetTitle("delta (mm)");
    g_delta->SetMarkerStyle(6);
    g_delta->SetMarkerColor(kRed+2);
    g_delta->Draw("AP");

    TCanvas *c2 = new TCanvas( "c2" , "Graphs" , 1000 , 800 );
    TGraph *gx = new TGraph(z.size(), &encoder_z[0], &x[0]);
    TGraph *gy = new TGraph(z.size(), &encoder_z[0], &y[0]);
    c2->Divide(1,2);
    
    c2->cd(1);
    gx->SetTitle(" encoder : lasertracker_x ");
    gx->GetXaxis()->SetTitle("encoder (mm)");
    gx->GetYaxis()->SetTitle("lasertracker_x (mm)");
    gx->SetMarkerStyle(6);
    gx->SetMarkerColor(kRed+2);
    gx->Draw("AP");

    c2->cd(2);
    gy->SetTitle(" encoder : lasertracker_y ");
    gy->GetXaxis()->SetTitle("encoder (mm)");
    gy->GetYaxis()->SetTitle("lasertracker_y (mm)");
    gy->SetMarkerStyle(6);
    gy->SetMarkerColor(kRed+2);
    gy->Draw("AP");*/


}
int main(){
    
}



