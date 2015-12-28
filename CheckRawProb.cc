#include <iostream>

void CheckRawProb(string histoname){

  TFile *MixFile = new TFile("RawProb.root"); // Get probability mixing root file 
  TH2F *ProbHisto = (TH2F*)MixFile->Get(histoname.c_str());  

  // Get X and Y axis for each histogram 
  TAxis *xaxis = ProbHisto->GetXaxis();
  TAxis *yaxis = ProbHisto->GetYaxis();
  // CosT Varies between -0.95 and 0.95. Energy between 0.1 and 10GeV.
  float CosT(-0.95), Energy(0.5);
  for(unsigned int i(0); i < 40; i ++){
    for(unsigned int j(0); j < 10; j++){
    Int_t binx = xaxis->FindBin(Energy);
    Int_t biny = yaxis->FindBin(CosT);

    float Prob = ProbHisto->GetBinContent(binx,biny);//(int)((CosT+1)*100));

    std::cout << "\n" << i << " : CosT = " << CosT << ", Energy = " << Energy << " :: Probability = " << Prob;
    Energy += 1;
    }
    CosT += 0.05; Energy = 0.5;
  }
  std::cout << "\nHello World!\n";
}
