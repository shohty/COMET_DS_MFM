#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
using namespace std;
/*
version2ではcsvファイルの末尾についたcsvなどの文字をとり、軸の位置＋withPillar.rootとなるように修正
*/

void FillToRoot(string inputFile, string outputFile){
    //open the csv file
    ifstream file(inputFile);
    if(!file.is_open()){
        cerr<<"failed to open csv file !!"<<endl;
        return ;
    }

    //create ROOT file and TTree
    TFile *rootFile =  TFile::Open(outputFile.c_str(), "RECREATE");
    TTree *tree = new TTree("tree","csv data");

    // names of rows of csv file
    double encoder, X, Y, Z, Bx, By, Bz, B;
    tree->Branch("encoder", &encoder, "encoder/D");
    tree->Branch("X", &X, "X/D");
    tree->Branch("Y", &Y, "Y/D");
    tree->Branch("Z", &Z, "Z/D");
    tree->Branch("Bx", &Bx, "Bx/D");
    tree->Branch("By", &By, "By/D");
    tree->Branch("Bz", &Bz, "Bz/D");
    tree->Branch("B", &B, "B/D");

    // Skip the first lines because it is a header line
    string line;
    /*for (int i = 0; i < 1; ++i) {
        if (!getline(file, line)) {
            cerr << "Error: There's no header line !!" << endl;
            return ;
        }
    }*/

    //read the csv file by a line
    while(getline(file, line)){
        istringstream ss(line);
        string value;
        if(getline(ss, value, ',')) encoder = stod(value);
        if(getline(ss, value, ',')) X = stod(value);
        if(getline(ss, value, ',')) Y = stod(value);
        if(getline(ss, value, ',')) Z= stod(value);
        if(getline(ss, value, ',')) Bx = stod(value);
        if(getline(ss, value, ',')) By = stod(value);
        if(getline(ss, value, ',')) Bz = stod(value);
        if(getline(ss, value, ',')) B = stod(value);
    //fill values to Tree
    tree->Fill();
    }
    //save root file and close
    tree->Write();
    rootFile->Close();
    cout<<"converting " <<inputFile << " to " <<outputFile<< " successfully done!"<<endl;
    return ;
}

int main(int argc, char** argv){
    string csvFile_directory = "/Users/shohtatakami/physics/COMETDS/DS189A_944turns/LaserTrackerPosition/csv/";
    string rootFile_directory = "/Users/shohtatakami/physics/COMETDS/DS189A_944turns/LaserTrackerPosition/root/";
    vector <string> scanfile = {
        "X_-200Y_0.csv_Map.csv",
        "X_-500Y_0.csv_Map.csv",
        "X_-650Y_0.csv_Map.csv",
        "X_-760Y_0.csv_Map.csv",
        "X_0Y_-200.csv_Map.csv",
        "X_0Y_-729.csv_Map.csv",
        "X_0Y_0.csv_Map.csv",
        "X_0Y_200.csv_Map.csv",
        "X_0Y_780.csv_Map.csv",
        "X_200Y_0.csv_Map.csv",
        "X_500Y_0.csv_Map.csv",
        "X_650Y_0.csv_Map.csv",
        "X_760Y_0.csv_Map.csv",
        "X_760Y_0rev.csv_Map.csv"
    };

    for(const auto& file : scanfile){
        string inputFile = csvFile_directory + file;
        //出力ファイル名の設定
        string originalName = file;
        size_t pos = originalName.find(".csv_Map.csv");
        if(pos != string::npos){
            originalName.erase(pos, 12);
        }
        string outputFile = rootFile_directory + originalName + ".root";
        FillToRoot(inputFile, outputFile);
    }
    return 0;
}


