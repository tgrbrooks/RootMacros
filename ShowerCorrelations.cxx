// C++ Includes
#include <iostream>
#include <string>
#include <vector>
#include <sstream>

// Root Includes
#include "THStack.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"

float Cor(vector<vector<float>> Events, int var1, int var2){
  // var = 0-19 hitno(t,u,v,y), tdcstd(t,u,v,y), tdciqr(t,u,v,y), adcamp(t,u,v,y), wfint(t,u,v,y)
  float mean1 = CalcMean(Events[var1]);
  float mean2 = CalcMean(Events[var2]);
  float sigma1 = CalcStandardDev(Events[var1]);
  float sigma2 = CalcStandardDev(Events[var2]);
  vector<float> var1var2;
  for(int i=0; i<Events[var1].size(); i++){
    var1var2.push_back(Events[var1][i]*Events[var2][i]);
  }
  float mean12 = CalcMean(var1var2);
  float correlation = (mean12 - mean1*mean2)/(sigma1*sigma2);
  return correlation;
}

float CalcMean(const vector <float> &Vec){
  float Mean = 0.0;
  for (unsigned int i =0; i < Vec.size(); i++){
    Mean += Vec[i];
  }
  return (Mean/Vec.size());
}

float CalcStandardDev(const vector <float> &Vec){
  float StandardDeviation = 0.0;
  float Mean = CalcMean(Vec);
  float sum = 0.0;
  float var = 0.0;
       
  for (unsigned int i =0; i < Vec.size(); i++){
    sum += pow((Vec[i] - Mean),2);
  }
  var = sqrt(sum/Vec.size());
  return var;
}

// Function to read in the ttree's from the input files and store the info in vector of vector of floats
vector<vector<float>> ReadTree(string RootFileName, bool isNeut){

  // Initilise float variables, vector of floats and vector of vector of floats
  float hitNo, hitNoU, hitNoV, hitNoY;
  float TDCstd, TDCstdU, TDCstdV, TDCstdY;
  float TDCiqr, TDCiqrU, TDCiqrV, TDCiqrY;
  float ADCamp, ADCampU, ADCampV, ADCampY;
  float Meanamp, MeanampU, MeanampV, MeanampY;
  float WFint, WFintU, WFintV, WFintY;
  float Meanint, MeanintU, MeanintV, MeanintY;
  float LowDen, LowDenU, LowDenV, LowDenY;
  float HiDen, HiDenU, HiDenV, HiDenY;
  float MeanRMS, MeanRMSU, MeanRMSV, MeanRMSY;
  float MeanMult, MeanMultU, MeanMultV, MeanMultY;
  float Wirestd, WirestdU, WirestdV, WirestdY;
  float Wireiqr, WireiqrU, WireiqrV, WireiqrY;
  vector<float> Evec;
  vector < vector < float > > vv_ReturnVect;

  TFile *TRootFile = new TFile(RootFileName.c_str());  
  TTree *FileTree = (TTree*)TRootFile->Get("ch_tree");
 
  // Set Branch Addresses to floats
  FileTree->SetBranchAddress("hitNo",&hitNo);
  FileTree->SetBranchAddress("hitNoU",&hitNoU);
  FileTree->SetBranchAddress("hitNoV",&hitNoV);
  FileTree->SetBranchAddress("hitNoY",&hitNoY);

  FileTree->SetBranchAddress("TDCstd", &TDCstd);
  FileTree->SetBranchAddress("TDCstdU", &TDCstdU);
  FileTree->SetBranchAddress("TDCstdV", &TDCstdV);
  FileTree->SetBranchAddress("TDCstdY", &TDCstdY);

  FileTree->SetBranchAddress("TDCiqr", &TDCiqr);
  FileTree->SetBranchAddress("TDCiqrU", &TDCiqrU);
  FileTree->SetBranchAddress("TDCiqrV", &TDCiqrV);
  FileTree->SetBranchAddress("TDCiqrY", &TDCiqrY);

  FileTree->SetBranchAddress("ADCamp", &ADCamp);
  FileTree->SetBranchAddress("ADCampU", &ADCampU);
  FileTree->SetBranchAddress("ADCampV", &ADCampV);
  FileTree->SetBranchAddress("ADCampY", &ADCampY);

  FileTree->SetBranchAddress("Meanamp", &Meanamp);
  FileTree->SetBranchAddress("MeanampU", &MeanampU);
  FileTree->SetBranchAddress("MeanampV", &MeanampV);
  FileTree->SetBranchAddress("MeanampY", &MeanampY);

  FileTree->SetBranchAddress("WFint", &WFint);
  FileTree->SetBranchAddress("WFintU", &WFintU);
  FileTree->SetBranchAddress("WFintV", &WFintV);
  FileTree->SetBranchAddress("WFintY", &WFintY);

  FileTree->SetBranchAddress("Meanint", &Meanint);
  FileTree->SetBranchAddress("MeanintU", &MeanintU);
  FileTree->SetBranchAddress("MeanintV", &MeanintV);
  FileTree->SetBranchAddress("MeanintY", &MeanintY);

  FileTree->SetBranchAddress("LowDen", &LowDen);
  FileTree->SetBranchAddress("LowDenU", &LowDenU);
  FileTree->SetBranchAddress("LowDenV", &LowDenV);
  FileTree->SetBranchAddress("LowDenY", &LowDenY);

  FileTree->SetBranchAddress("HiDen", &HiDen);
  FileTree->SetBranchAddress("HiDenU", &HiDenU);
  FileTree->SetBranchAddress("HiDenV", &HiDenV);
  FileTree->SetBranchAddress("HiDenY", &HiDenY);

  FileTree->SetBranchAddress("MeanRMS", &MeanRMS);
  FileTree->SetBranchAddress("MeanRMSU", &MeanRMSU);
  FileTree->SetBranchAddress("MeanRMSV", &MeanRMSV);
  FileTree->SetBranchAddress("MeanRMSY", &MeanRMSY);

  FileTree->SetBranchAddress("MeanMult", &MeanMult);
  FileTree->SetBranchAddress("MeanMultU", &MeanMultU);
  FileTree->SetBranchAddress("MeanMultV", &MeanMultV);
  FileTree->SetBranchAddress("MeanMultY", &MeanMultY);

  FileTree->SetBranchAddress("Wirestd", &Wirestd);
  FileTree->SetBranchAddress("WirestdU", &WirestdU);
  FileTree->SetBranchAddress("WirestdV", &WirestdV);
  FileTree->SetBranchAddress("WirestdY", &WirestdY);

  FileTree->SetBranchAddress("Wireiqr", &Wireiqr);
  FileTree->SetBranchAddress("WireiqrU", &WireiqrU);
  FileTree->SetBranchAddress("WireiqrV", &WireiqrV);
  FileTree->SetBranchAddress("WireiqrY", &WireiqrY);

  FileTree->SetBranchAddress("Evec", &Evec);

  vector<float> vhitNo, vhitNoU, vhitNoV, vhitNoY;
  vector<float> vTDCstd, vTDCstdU, vTDCstdV, vTDCstdY;
  vector<float> vTDCiqr, vTDCiqrU, vTDCiqrV, vTDCiqrY;
  vector<float> vADCamp, vADCampU, vADCampV, vADCampY;
  vector<float> vMeanamp, vMeanampU, vMeanampV, vMeanampY;
  vector<float> vWFint, vWFintU, vWFintV, vWFintY;
  vector<float> vMeanint, vMeanintU, vMeanintV, vMeanintY;
  vector<float> vLowDen, vLowDenU, vLowDenV, vLowDenY;
  vector<float> vHiDen, vHiDenU, vHiDenV, vHiDenY;
  vector<float> vMeanRMS, vMeanRMSU, vMeanRMSV, vMeanRMSY;
  vector<float> vMeanMult, vMeanMultU, vMeanMultV, vMeanMultY;
  vector<float> vWirestd, vWirestdU, vWirestdV, vWirestdY;
  vector<float> vWireiqr, vWireiqrU, vWireiqrV, vWireiqrY;
  vector<float> vEnergy;

  Long64_t EntryNumber = FileTree->GetEntries();
  for(unsigned int j=0; j!=EntryNumber; ++j){ // Loop over tree entries
    // Fill vectors of floats with each branch address data member.
    FileTree->GetEntry(j);

    vhitNo.push_back(hitNo);
    vhitNoU.push_back(hitNoU);
    vhitNoV.push_back(hitNoV);
    vhitNoY.push_back(hitNoY);

    vTDCstd.push_back(TDCstd);
    vTDCstdU.push_back(TDCstdU);
    vTDCstdV.push_back(TDCstdV);
    vTDCstdY.push_back(TDCstdY);

    vTDCiqr.push_back(TDCiqr);
    vTDCiqrU.push_back(TDCiqrU);
    vTDCiqrV.push_back(TDCiqrV);
    vTDCiqrY.push_back(TDCiqrY);

    vADCamp.push_back(ADCamp);
    vADCampU.push_back(ADCampU);
    vADCampV.push_back(ADCampV);
    vADCampY.push_back(ADCampY);

    vMeanamp.push_back(Meanamp);
    vMeanampU.push_back(MeanampU);
    vMeanampV.push_back(MeanampV);
    vMeanampY.push_back(MeanampY);

    vWFint.push_back(WFint);
    vWFintU.push_back(WFintU);
    vWFintV.push_back(WFintV);
    vWFintY.push_back(WFintY);

    vMeanint.push_back(Meanint);
    vMeanintU.push_back(MeanintU);
    vMeanintV.push_back(MeanintV);
    vMeanintY.push_back(MeanintY);

    vLowDen.push_back(LowDen);
    vLowDenU.push_back(LowDenU);
    vLowDenV.push_back(LowDenV);
    vLowDenY.push_back(LowDenY);

    vHiDen.push_back(HiDen);
    vHiDenU.push_back(HiDenU);
    vHiDenV.push_back(HiDenV);
    vHiDenY.push_back(HiDenY);

    vMeanRMS.push_back(MeanRMS);
    vMeanRMSU.push_back(MeanRMSU);
    vMeanRMSV.push_back(MeanRMSV);
    vMeanRMSY.push_back(MeanRMSY);

    vMeanMult.push_back(MeanMult);
    vMeanMultU.push_back(MeanMultU);
    vMeanMultV.push_back(MeanMultV);
    vMeanMultY.push_back(MeanMultY);

    vWirestd.push_back(Wirestd);
    vWirestdU.push_back(WirestdU);
    vWirestdV.push_back(WirestdV);
    vWirestdY.push_back(WirestdY);

    vWireiqr.push_back(Wireiqr);
    vWireiqrU.push_back(WireiqrU);
    vWireiqrV.push_back(WireiqrV);
    vWireiqrY.push_back(WireiqrY);

    vEnergy.push_back(Evec[0]);
  }
  // Push back data read in from tree onto vector of vectors. Columns of 2D vector correspond to events 1,2,3,...,N 
  // While rows correspond to the different cut variables.
  vv_ReturnVect.push_back(vhitNo);   vv_ReturnVect.push_back(vhitNoU);   vv_ReturnVect.push_back(vhitNoV);   vv_ReturnVect.push_back(vhitNoY);
  vv_ReturnVect.push_back(vTDCstd); vv_ReturnVect.push_back(vTDCstdU); vv_ReturnVect.push_back(vTDCstdV); vv_ReturnVect.push_back(vTDCstdY);
  vv_ReturnVect.push_back(vTDCiqr); vv_ReturnVect.push_back(vTDCiqrU); vv_ReturnVect.push_back(vTDCiqrV); vv_ReturnVect.push_back(vTDCiqrY);
  vv_ReturnVect.push_back(vADCamp); vv_ReturnVect.push_back(vADCampU); vv_ReturnVect.push_back(vADCampV); vv_ReturnVect.push_back(vADCampY);
  vv_ReturnVect.push_back(vWFint);  vv_ReturnVect.push_back(vWFintU);  vv_ReturnVect.push_back(vWFintV);  vv_ReturnVect.push_back(vWFintY);
  vv_ReturnVect.push_back(vMeanamp); vv_ReturnVect.push_back(vMeanampU); vv_ReturnVect.push_back(vMeanampV); vv_ReturnVect.push_back(vMeanampY);
  vv_ReturnVect.push_back(vMeanint); vv_ReturnVect.push_back(vMeanintU); vv_ReturnVect.push_back(vMeanintV); vv_ReturnVect.push_back(vMeanintY);
  vv_ReturnVect.push_back(vMeanRMS); vv_ReturnVect.push_back(vMeanRMSU); vv_ReturnVect.push_back(vMeanRMSV); vv_ReturnVect.push_back(vMeanRMSY);
  vv_ReturnVect.push_back(vMeanMult); vv_ReturnVect.push_back(vMeanMultU); vv_ReturnVect.push_back(vMeanMultV); vv_ReturnVect.push_back(vMeanMultY);
  vv_ReturnVect.push_back(vLowDen); vv_ReturnVect.push_back(vLowDenU); vv_ReturnVect.push_back(vLowDenV); vv_ReturnVect.push_back(vLowDenY);
  vv_ReturnVect.push_back(vHiDen); vv_ReturnVect.push_back(vHiDenU); vv_ReturnVect.push_back(vHiDenV); vv_ReturnVect.push_back(vHiDenY);
  vv_ReturnVect.push_back(vWirestd); vv_ReturnVect.push_back(vWirestdU); vv_ReturnVect.push_back(vWirestdV); vv_ReturnVect.push_back(vWirestdY);
  vv_ReturnVect.push_back(vWireiqr); vv_ReturnVect.push_back(vWireiqrU); vv_ReturnVect.push_back(vWireiqrV); vv_ReturnVect.push_back(vWireiqrY);
  vv_ReturnVect.push_back(vEnergy);

  return vv_ReturnVect;

} // End of ReadTree


// Main Script
void CorrelationFunctions(){

  string CutVariables[20] = {
    "hitNo","hitNoU","hitNoV","hitNoY","TDCstd","TDCstdU","TDCstdV","TDCstdY","TDCiqr","TDCiqrU","TDCiqrV","TDCiqrY",
    "ADCamp","ADCampU","ADCampV","ADCampY","WFint","WFintU","WFintV","WFintY"};

  // List nnbar and Neutrino files to be merged
  const int NumNbarFiles = 100;  string nnbarFileNames[NumNbarFiles];
  const int NumNeutFiles = 100; string NeutFileNames[NumNeutFiles];
  
  // Define String stream and temporary string to hold the file number
  stringstream ss; string FileNumAsString;
  // Loop over the number of files and allocate the file names to string vectors.
  // Use string stream to get the file number as a string
  for(unsigned int i(0); i < NumNbarFiles; i ++){
    ss << i; ss >> FileNumAsString; // Convert integer i of loop to string
    nnbarFileNames[i] = "outputs_run2/Nnbar_output_run2_" + FileNumAsString + ".root";
    NeutFileNames[i] = "outputs_run2/Atmo_output_run2_" + FileNumAsString + ".root";
    ss.str(""); ss.clear(); // Clear the string stream
  }

  gROOT->ProcessLine(".x lhcbStyle.C"); // Using lhcb Style file
 
  std::ofstream myfile; // Create output file stream with name AnalyserDataOutput.txt
  myfile.open("CorrelationFunctions.txt");
  
  // Create two TChains, one for nnbar trees and one for Neutrino data trees
  TChain nnbarTrees("ch_tree");
  TChain NeutTrees("ch_tree");
  
  // Message to output file
  myfile << "||CorrelationFunctions.cc|| Begin.\n||Analyser.cc|| The Following " << NumNbarFiles << " NNbar Files were merged:\n";  

  // Loop over nnbar files and add them to TChain
  for(unsigned int i(0); i < NumNbarFiles; i ++){
    nnbarTrees.Add(nnbarFileNames[i].c_str()); 
    myfile << nnbarFileNames[i] << "\n"; // Print File names to outputfile
  }
  nnbarTrees.Merge("allnnbar.root"); // Merge all Chained trees into a single root file
   
  // Call function GetChainSTD to get the Standard Deviations and Means from the merged allnnbar.root
  // containing all nnbar data. Feed in dummy empty vector "Means" which is not used as second argument is false
  vector<vector<float>> NnbarEvents = ReadTree("allnnbar.root",false);
  
  myfile << "\n" << "Number of nnbar events: "<<NnbarEvents[0].size();

  // Message to output file
  myfile << "\n\n||Analyser.cc|| The Folllowing " << NumNeutFiles << " Neutrino Files were analysed:\n";

  // Loop over Neutrino files and add them to TChain NeutTrees
  for(unsigned int i(0); i < NumNeutFiles; i ++){
    NeutTrees.Add(NeutFileNames[i].c_str()); 
    myfile << NeutFileNames[i] << "\n"; // Print File names to outputfile
  }  
  NeutTrees.Merge("allneut.root"); // MAY NOT NEED THIS LINE ANYMORE
  vector<vector<float>> NeutEvents = ReadTree("allneut.root",true); 

  myfile << "\nCorrelation Functions:\n\n";
  for(int k = 0; k < 20; k++){
    for(int j = 0; j < 20; j++){
      float nnbarcorr = Cor(NnbarEvents,k,j);
      myfile << "\nnnbar: Cor(" << CutVariables[k] << ", " << CutVariables[j] << ") = "<< nnbarcorr << "\n";
      float neutcorr = Cor(NeutEvents,k,j);
      myfile << "\nneutrino: Cor(" << CutVariables[k] << ", " << CutVariables[j] << ") = "<< neutcorr << "\n";
    }
  }

  myfile.close();  

}
