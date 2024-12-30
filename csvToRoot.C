#include <TFile.h>
#include <TTree.h>
using namespace std;

int main(){
    string inputFile = "/Users/shohtatakami/physics/COMETDS/DS189A_944turns/DS189A_944turns_Map.csv";
    string outputFile = "/Users/shohtatakami/physics/COMETDS/DS189A_944turns/DS189A_944turns_Map.root";

    //open the csv file
    ifstream file(inputFile);
    if(!file.is_open()){
        cerr<<"failed to open csv file !!"<< endl;
        return 1;
    }

    //create ROOT file and TTree
    TFile *rootFile =  TFile::Open(outputFile.c_str(), "RECREATE");
    TTree *tree = new TTree("tree","csv data");

    // names of rows of csv file
    double x,y,z,Bx,By,Bz,B;
    tree->Branch("x", &x, "x/D");
    tree->Branch("y", &y, "y/D");
    tree->Branch("z", &z, "z/D");
    tree->Branch("Bx", &Bx, "Bx/D");
    tree->Branch("By", &By, "By/D");
    tree->Branch("Bz", &Bz, "Bz/D");
    tree->Branch("B", &B, "B/D");

    // Skip the first 9 lines
    string line;
    for (int i = 0; i < 9; ++i) {
        if (!getline(file, line)) {
            cerr << "Error: Not enough lines in the file to skip 9 lines!" << endl;
            return 1;
        }
    }

    //read the csv file by a line
    while(getline(file, line)){
        istringstream ss(line);
        string value;

        if(getline(ss, value, '\t')) x = stod(value);
        if(getline(ss, value, '\t')) y = stod(value);
        if(getline(ss, value, '\t')) z = stod(value);
        if(getline(ss, value, '\t')) Bx = stod(value);
        if(getline(ss, value, '\t')) By = stod(value);
        if(getline(ss, value, '\t')) Bz = stod(value);
        if(getline(ss, value, '\t')) B = stod(value);
    
    //fill values to Tree
    tree->Fill();
    }
    //save root file and close
    tree->Write();
    rootFile->Close();
    cout<<"converting csv to root successfully done "<<endl;
    return 0;
}
