// C++ Includes
#include <iostream>
#include <string>
#include <vector>
#include <sstream>

// Root Includes
#include "THStack.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMVA/Types.h"

// Main Script
void AnalyseDecisionTree(){

  const int NumOfVars = 52;
  string CutVariables[NumOfVars] = {
    "hitNo","hitNoU","hitNoV","hitNoY","TDCstd","TDCstdU","TDCstdV","TDCstdY","TDCiqr","TDCiqrU","TDCiqrV","TDCiqrY",
    "ADCamp","ADCampU","ADCampV","ADCampY","WFint","WFintU","WFintV","WFintY","Meanamp","MeanampU","MeanampV","MeanampY",
    "Meanint","MeanintU","MeanintV","MeanintY","LowDen","LowDenU","LowDenV","LowDenY","HiDen","HiDenU","HiDenV","HiDenY",
    "MeanRMS","MeanRMSU","MeanRMSV","MeanRMSY","MeanMult","MeanMultU","MeanMultV","MeanMultY",
    "Wirestd","WirestdU","WirestdV","WirestdY","Wireiqr","WireiqrU","WireiqrV","WireiqrY" };

 //  BNBTrees.Merge("allBNB.root");
 // EminusTrees.Merge("allEminus.root"); // Merge all Chained trees into a single root file
 
  // INITIATE TMVA::Reader METHOD

  // Input data file
  TFile* data = new TFile("../RootOutputs/ShowerAnaOutputs/ShowerAna_Eminus_output_143.root"); 
  TTree* dataTree = (TTree*)(data->Get("ch_tree"));

  // Output file
  TFile *target = new TFile("real_data-mva_output.root","RECREATE" );
  TTree *tree = new TTree("tree","treelibrated tree");

  TMVA::Reader* reader = new TMVA::Reader( "V:Color:!Silent" ); //V:Color:!Silent

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

  reader->AddVariable("hitNo",  &hitNo);
  reader->AddVariable("hitNoU", &hitNoU);
  reader->AddVariable("hitNoV", &hitNoV);
  reader->AddVariable("hitNoY", &hitNoY);

  reader->AddVariable("TDCstd", &TDCstd);
  reader->AddVariable("TDCstdU", &TDCstdU);
  reader->AddVariable("TDCstdV", &TDCstdV);
  reader->AddVariable("TDCstdY", &TDCstdY);

  reader->AddVariable("TDCiqr", &TDCiqr);
  reader->AddVariable("TDCiqrU", &TDCiqrU);
  reader->AddVariable("TDCiqrV", &TDCiqrV);
  reader->AddVariable("TDCiqrY", &TDCiqrY);

  reader->AddVariable("ADCamp", &ADCamp);
  reader->AddVariable("ADCampU", &ADCampU);
  reader->AddVariable("ADCampV", &ADCampV);
  reader->AddVariable("ADCampY", &ADCampY);

  reader->AddVariable("WFint", &WFint);
  reader->AddVariable("WFintU", &WFintU);
  reader->AddVariable("WFintV", &WFintV);
  reader->AddVariable("WFintY", &WFintY);

  reader->AddVariable("Meanamp", &Meanamp);
  reader->AddVariable("MeanampU", &MeanampU);
  reader->AddVariable("MeanampV", &MeanampV);
  reader->AddVariable("MeanampY", &MeanampY);

  reader->AddVariable("Meanint", &Meanint);
  reader->AddVariable("MeanintU", &MeanintU);
  reader->AddVariable("MeanintV", &MeanintV);
  reader->AddVariable("MeanintY", &MeanintY);

  reader->AddVariable("LowDen", &LowDen);
  reader->AddVariable("LowDenU", &LowDenU);
  reader->AddVariable("LowDenV", &LowDenV);
  reader->AddVariable("LowDenY", &LowDenY);

  reader->AddVariable("HiDen", &HiDen);
  reader->AddVariable("HiDenU", &HiDenU);
  reader->AddVariable("HiDenV", &HiDenV);
  reader->AddVariable("HiDenY", &HiDenY);

  reader->AddVariable("MeanRMS", &MeanRMS);
  reader->AddVariable("MeanRMSU", &MeanRMSU);
  reader->AddVariable("MeanRMSV", &MeanRMSV);
  reader->AddVariable("MeanRMSY", &MeanRMSY);

  reader->AddVariable("MeanMult", &MeanMult);
  reader->AddVariable("MeanMultU", &MeanMultU);
  reader->AddVariable("MeanMultV", &MeanMultV);
  reader->AddVariable("MeanMultY", &MeanMultY);

  reader->AddVariable("Wirestd", &Wirestd);
  reader->AddVariable("WirestdU", &WirestdU);
  reader->AddVariable("WirestdV", &WirestdV);
  reader->AddVariable("WirestdY", &WirestdY);

  reader->AddVariable("Wireiqr", &Wireiqr);
  reader->AddVariable("WireiqrU", &WireiqrU);
  reader->AddVariable("WireiqrV", &WireiqrV);
  reader->AddVariable("WireiqrY", &WireiqrY);

  // Add branch to output tree of BDT response
  Double_t BDT_response;
  tree->Branch("BDT_response",&BDT_response);
 
  reader->BookMVA( "BDT", "weights/JobName_BDT.weights.xml" );

  UInt_t userVar[3];
  dataTree->SetBranchAddress("hitNo",&userVar[0]);

  for (Long64_t ievt=0; ievt < dataTree->GetEntries(); ievt++) {

    if (ievt%100000 == 0) std::cout << "\n\n--- ... Processing event: " << ievt <<std::endl;
    
    dataTree->GetEntry(ievt);
    
    
    BDT_response = reader->EvaluateMVA("BDT");  // THIS LINE CAUSES AN ERROR

    // Error on classifier response - must be called after "EvaluateMVA"
    // (not available for all methods, returns -1 in that case)
    Double_t mvaErr   = reader->GetMVAError();
    if(mvaErr != 0)  std::cout << "\nMVA error status is: " << mvaErr;

    // Regression response for one target
    Double_t regValue = (reader->EvaluateRegression( "BDT" ))[0];
    Bool_t passed = reader->EvaluateMVA( "BDT",  0.01 );
//    Double_t pSig = reader->GetProba( "BDT", 0.01 );
//    Double_t rarity = reader->GetRarity( "BDT" );


    if( passed == true ) std::cout << "\nThe event passed the cut!";
    else std::cout << "\nThe event failed the cut!";
    std::cout << " RegValue = " << regValue;// << ", pSig = " << pSig << ", rarity = " << rarity;
    tree->Fill();
  }
 
  tree->Write();

  //Double_t pSig = reader->GetProba( "<YourClassifierName>", sigFrac );


  // factory->AddRegressionTree( regTree, weight );

  // END TMVA::Reader METHOD

}
