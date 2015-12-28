// C++ Includes
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
// Root Includes
#include <TH2>
#include <TFile>

void CalcMix(string FileName, string histoname){

  string FilePath = "UnmixedFluxFiles/";
  TFile *MixFile = new TFile("RawProb.root"); // Get probability mixing root file
  // Get histograms of relevant mixing from file
  TH2F *ProbHisto = (TH2F*)MixFile->Get(histoname.c_str());  

  // Get X and Y axis for each histogram 
  TAxis *xaxis = ProbHisto->GetXaxis();
  TAxis *yaxis = ProbHisto->GetYaxis();

  ifstream File;
  File.open((FilePath+FileName).c_str());
  string Throw; std::getline(File,Throw); 

  ofstream Output; 
  float Energy, CosT, Flux, Junk;

  Output.open(("New_"+FileName).c_str()); 
  Output << "#A LINE OF JUNK\n";
  while( File >> Energy >> CosT >> Flux >> Junk >> Junk ){
 
    Int_t binx = xaxis->FindBin(Energy);
    Int_t biny = yaxis->FindBin(CosT);
    float Prob(0);
    if(Energy < 200) Prob = ProbHisto->GetBinContent(binx,biny); // If statement for high energy fluxes
    
    Output << Energy << " " << CosT << " " << Prob*Flux << " " << Junk << " " << Junk << "\n";

  }
  File.close(); Output.close();
}

void CalcEMuMix(string FileName1, string FileName2, string histoname1, string NewName){

  string FilePath = "UnmixedFluxFiles/";  
  TFile *MixFile = new TFile("RawProb.root"); // Get probability mixing root file
  // Get histograms of relevant mixing from file
  TH2F *ProbHisto = (TH2F*)MixFile->Get(histoname1.c_str());  

  // Get X and Y axis for each histogram 
  TAxis *xaxis = ProbHisto->GetXaxis();
  TAxis *yaxis = ProbHisto->GetYaxis();

  ifstream File1, File2;
  ofstream Output; string Throw;

  float Energy, CosT, Flux1, Flux2, Junk, NewFlux(0);

  File1.open((FilePath+FileName1).c_str()); File2.open(FileName2.c_str());
  std::getline(File1,Throw); std::getline(File2,Throw);
  Output.open(NewName.c_str());
  Output << "#A LINE OF JUNK!\n";
  while( File1 >> Energy >> CosT >> Flux1 >> Junk >> Junk ){
    Int_t binx = xaxis->FindBin(Energy);
    Int_t biny = yaxis->FindBin(CosT);

    float Prob = ProbHisto->GetBinContent(binx,biny);//(int)((CosT+1)*100));
 
    File2 >> Junk >> Junk >> Flux2 >> Junk >> Junk;

    Output << Energy << " " << CosT << " " << (Flux2 + Flux1*Prob) << " " << Junk << " " << Junk << "\n";
  }
  File1.close(); File2.close(); Output.close();  
}

void CalcTauMix(string FileName1, string FileName2, string histoname1, string histoname2,string NewName){


  string FilePath = "UnmixedFluxFiles/";  
  TFile *MixFile = new TFile("RawProb.root"); // Get probability mixing root file
  // Get histograms of relevant mixing from file
  TH2F *ProbHisto1 = (TH2F*)MixFile->Get(histoname1.c_str());  
  TH2F *ProbHisto2 = (TH2F*)MixFile->Get(histoname2.c_str());

  // Get X and Y axis for each histogram 
  TAxis *xaxis1 = ProbHisto1->GetXaxis();
  TAxis *xaxis2 = ProbHisto2->GetXaxis();
  TAxis *yaxis1 = ProbHisto1->GetYaxis();
  TAxis *yaxis2 = ProbHisto2->GetYaxis();

  ifstream File1, File2;
  ofstream Output; string Throw;

  float Energy1, CosT1, Energy2, CosT2, Flux1, Flux2, Junk, NewFlux(0);

  File1.open((FilePath+FileName1).c_str()); File2.open((FilePath+FileName2).c_str());
  std::getline(File1,Throw); std::getline(File2,Throw);
  Output.open(NewName.c_str());
  Output << "#A LINE OF JUNK\n";
  while( File1 >> Energy1 >> CosT1 >> Flux1 >> Junk >> Junk ){

    Int_t binx1 = xaxis1->FindBin(Energy1);
    Int_t biny1 = yaxis1->FindBin(CosT1);

    File2 >> Energy2 >> CosT2 >> Flux2 >> Junk >> Junk;
    Int_t binx2 = xaxis2->FindBin(Energy2);
    Int_t biny2 = yaxis2->FindBin(CosT2);

    float Prob1 = ProbHisto1->GetBinContent(binx1,biny1);//(int)((CosT+1)*100));
    float Prob2 = ProbHisto2->GetBinContent(binx2,biny2);//(int)((CosT+1)*100));
   
    // }else{ NewFlux = Flux1; }
    Output << Energy1 << " " << CosT1 << " " << (Flux1*Prob1 + Flux2*Prob2) << " " << Junk << " " << Junk << "\n";
  }
  File1.close(); File2.close(); Output.close();  
}

void Mixer(){

  //SNO_MinENu.txt    SNO_MinEbNu.txt   SNO_MinMuNu.txt   SNO_MinMubNu.txt  SNO_MinTNu.txt    SNO_MinTbNu.txt   
  
  // FOR MATTER NEUTRINOS

  // Call function to open and write each file, providing the names of the mixing histograms and the flavour being mixed
  // Matter Mixing
  CalcMix("KAM_MinENuHE.txt","NuEToNuE3f");
  CalcEMuMix("KAM_MinMuNuHE.txt","New_KAM_MinENuHE.txt","NuMuToNuE3f","New_New_KAM_MinENuHE.txt");
  CalcMix("KAM_MinMuNuHE.txt","NuMuToNuMu3f");
  CalcEMuMix("KAM_MinENuHE.txt","New_KAM_MinMuNuHE.txt","NuEToNuMu3f","New_New_KAM_MinMuNuHE.txt");
  CalcTauMix("KAM_MinENuHE.txt","KAM_MinMuNuHE.txt","NuEToNuTau3f","NuMuToNuTau3f","New_KAM_MinTNuHE.txt");

  // Anti-Matter Mixing
  CalcMix("KAM_MinEbNuHE.txt","NuEToNuE3f");
  CalcEMuMix("KAM_MinMubNuHE.txt","New_KAM_MinEbNuHE.txt","NuMuToNuE3f","New_New_KAM_MinEbNuHE.txt");
  CalcMix("KAM_MinMubNuHE.txt","NuMuToNuMu3f");
  CalcEMuMix("KAM_MinEbNuHE.txt","New_KAM_MinMubNuHE.txt","NuEToNuMu3f","New_New_KAM_MinMubNuHE.txt");
  CalcTauMix("KAM_MinEbNuHE.txt","KAM_MinMubNuHE.txt","NuEToNuTau3f","NuMuToNuTau3f","New_KAM_MinTbNuHE.txt");


  std::cout << "\nHello World!\n";
}
