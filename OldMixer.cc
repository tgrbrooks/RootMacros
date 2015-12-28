// C++ Includes
#include <iostream>
#include <string>
#include <vector>

// Root Includes
#include <TH2>
#include <TFile>

void CalcMix(vector <string> FileNames, string histoname1, string histoname2, string histoname3,int MixFlav){

  TFile *MixFile = new TFile("RawProb.root");
  TH2F *ProbHisto1 = (TH2F*)MixFile->Get(histoname1.c_str());  
  TH2F *ProbHisto2 = (TH2F*)MixFile->Get(histoname2.c_str());
  TH2F *ProbHisto3 = (TH2F*)MixFile->Get(histoname3.c_str());

  ifstream File1, File2, File3;
  ofstream Output1;
  float Energy, CosT, Flux1, Flux2, Flux3, Junk, NewFlux(0);
  for(unsigned int i(0); i < 1; i ++){
    File1.open(FileNames[i*3].c_str()); File2.open(FileNames[1+i*3].c_str()); File3.open(FileNames[2+i*3].c_str());
    Output1.open(("New_"+FileNames[i*3+MixFlav-1]).c_str()); 
    Output1 << "#A LINE OF JUNK\n";
    std::cout << "\nIT GES HERE!";
    while( File1 >> Energy >> CosT >> Flux1 >> Junk >> Junk ){

      File2 >> Junk >> Junk >> Flux2 >> Junk >> Junk;
      File3 >> Junk >> Junk >> Flux3 >> Junk >> Junk;

      float Prob1 = ProbHisto1->GetBinContent((int)(Energy*30),(int)((CosT+1)*100));
      float Prob2 = ProbHisto2->GetBinContent((int)(Energy*30),(int)((CosT+1)*100));
      float Prob3 = ProbHisto3->GetBinContent((int)(Energy*30),(int)((CosT+1)*100));

      if(MixFlav == 1) NewFlux = Flux1 - Flux1*Prob1 - Flux1*Prob2 + Flux3*Prob3;
      if(MixFlav == 2) NewFlux = Flux2 - Flux2*Prob1 - Flux2*Prob2 + Flux1*Prob3;
      if(MixFlav == 3) NewFlux = Flux1*Prob1 + Flux2*Prob2;
      
      Output1 << Energy << " " << CosT << " " << NewFlux << " " << Junk << " " << Junk << "\n;
    }
    File1.close(); File2.close(); File3.close(); Output1.close();
  }
}

void Mixer(){

  //SNO_MinENu.txt    SNO_MinEbNu.txt   SNO_MinMuNu.txt   SNO_MinMubNu.txt  SNO_MinTNu.txt    SNO_MinTbNu.txt   
  string FilePath = "FluxFiles/SNO/UnMixed/";
  vector <string> FileNames;
  FileNames.push_back( FilePath + "SNO_MinENu.txt" );
  FileNames.push_back( FilePath + "SNO_MinMuNu.txt" );
  FileNames.push_back( FilePath + "SNO_MinTNu.txt" );
  FileNames.push_back( FilePath + "SNO_MinEbNu.txt" );
  FileNames.push_back( FilePath + "SNO_MinMubNu.txt" );
  FileNames.push_back( FilePath + "SNO_MinTbNu.txt" );

  CalcMix(FileNames,"NuEToNuMu3f","NuEToNuT3f","NuMuToNuE3f",1);
  CalcMix(FileNames,"NuMuToNuE3f","NuMuToNuT3f","NuEToNuMu3f",2);
  CalcMix(FileNames,"NuEToNuT3f","NuMuToNuT3f","NuMuToNuE3f",3);

  std::cout << "\nHello World!";
}
