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
void TrainandTestDecisionTree(){

  const int NumOfVars = 52;
  string CutVariables[NumOfVars] = {
    "hitNo","hitNoU","hitNoV","hitNoY","TDCstd","TDCstdU","TDCstdV","TDCstdY","TDCiqr","TDCiqrU","TDCiqrV","TDCiqrY",
    "ADCamp","ADCampU","ADCampV","ADCampY","WFint","WFintU","WFintV","WFintY","Meanamp","MeanampU","MeanampV","MeanampY",
    "Meanint","MeanintU","MeanintV","MeanintY","LowDen","LowDenU","LowDenV","LowDenY","HiDen","HiDenU","HiDenV","HiDenY",
    "MeanRMS","MeanRMSU","MeanRMSV","MeanRMSY","MeanMult","MeanMultU","MeanMultV","MeanMultY",
    "Wirestd","WirestdU","WirestdV","WirestdY","Wireiqr","WireiqrU","WireiqrV","WireiqrY" };

  // List of Eminus and BNB files to be merged
  const int NumEminusFiles = 195;  string EminusFileNames[NumEminusFiles];
  const int NumBNBFiles = 325; string BNBFileNames[NumBNBFiles];
  
  // Define String stream and temporary string to hold the file number
  stringstream ss; string FileNumAsString;
  // Loop over the number of files and allocate the file names to string vectors.
  // Use string stream to get the file number as a string
  for(unsigned int i(0); i < NumBNBFiles; i ++){
    ss << i; ss >> FileNumAsString; // Convert integer i of loop to string
    if( i < 195 ){
      EminusFileNames[i] = "../RootOutputs/ShowerAnaOutputs/ShowerAna_Eminus_output_" + FileNumAsString + ".root";
    }
    BNBFileNames[i] = "../RootOutputs/BNBAnaOutputs/ShowerAna_BNB_output_" + FileNumAsString + ".root"; 
    ss.str(""); ss.clear(); // Clear the string stream
  }

  gROOT->ProcessLine(".x lhcbStyle.C"); // Using lhcb Style file

  std::ofstream myfile; // Create output file stream with name ShowerPlotter.txt
  myfile.open("TrainandTestDecisionTree.txt");

  // Create two TChains, one for eminus trees and one for BNB data trees
  TChain EminusTrees("ch_tree");
  TChain BNBTrees("ch_tree");

  // Message to output file
  myfile << "\n\n||TrainandTestDecisionTree.cc|| The Following " << NumBNBFiles << " BNB Files will be merged.\n";

  // Loop over BNB files and add them to TChainBNBTrees
  for(unsigned int i(0); i < NumBNBFiles; i ++){
    BNBTrees.Add(BNBFileNames[i].c_str()); 
    myfile << BNBFileNames[i] << "\n"; // Print File names to outputfile
  }  
  BNBTrees.Merge("allBNB.root");

  
  // Message to output file
  myfile << "\n||TrainandTestDecisionTree.cc|| Success!\n||TrainandTestDecisionTree.cc|| The Folllowing " << NumEminusFiles << " Eminus Files will be merged.\n";  
   // Loop over nnbar files and add them to TChain
  for(unsigned int i(0); i < NumEminusFiles; i ++){
    EminusTrees.Add(EminusFileNames[i].c_str()); 
    myfile << EminusFileNames[i] << "\n"; // Print File names to outputfile
  }
  EminusTrees.Merge("allEminus.root"); // Merge all Chained trees into a single root file
  
  std::cout << "\nEminus and BNB Files successfully read and merged. Now attempting to gather data from chained file.\n";
 
  // INITIATE TMVA::Factory METHOD

  TFile *NodeFile = new TFile("NodeFile.root","RECREATE");
  TMVA::Factory* factory = new TMVA::Factory( "JobName", NodeFile, "" ); 


  TFile *SrcSigTreeTrain = new TFile( "allEminus.root " );
  TFile *SrcBkgTreeTrain = new TFile( "allBNB.root"     );
  TFile *SrcSigTreeTest  = new TFile( "../RootOutputs/ShowerOutputs/ShowerAna_output_1.root"  );
  TFile *SrcBkgTreeTest  = new TFile( "../RootOutputs/BNBAnaOutputs/ShowerAna_BNB_output_0.root"        );

  // Get the signal and background trees from TFile source(s);
  // multiple trees can be registered with the Factory
  // Get the signal and background training and test trees from TFile source(s)
  TTree* sigTreeTrain = (TTree*)SrcSigTreeTrain->Get( "ch_tree" );
  TTree* bkgTreeTrain = (TTree*)SrcBkgTreeTrain->Get( "ch_tree" );
  TTree* sigTreeTest  = (TTree*)SrcSigTreeTest ->Get( "ch_tree" );
  TTree* bkgTreeTest  = (TTree*)SrcBkgTreeTest ->Get( "ch_tree" );

  // Set the event weights per tree (these weights are applied in
  // addition to individual event weights that can be specified)
  Double_t sigWeight  = 1.0; Double_t bkgWeight = 1.0;

  factory->AddSignalTree    ( sigTreeTrain, sigWeight, TMVA::Types::kTraining );
  factory->AddBackgroundTree( bkgTreeTrain, bkgWeight, TMVA::Types::kTraining );
  factory->AddSignalTree    ( sigTreeTest,  sigWeight, TMVA::Types::kTesting  );
  factory->AddBackgroundTree( bkgTreeTest,  bkgWeight, TMVA::Types::kTesting  );
  
  for(unsigned int i(0); i < 4;         i ++)  factory->AddVariable( CutVariables[i].c_str(), 'I' ); 
  for(unsigned int i(4); i < NumOfVars; i ++)  factory->AddVariable( CutVariables[i].c_str(), 'F' );
  
  //  factory->AddTarget ( "Evec" );
/*
  TCut preCut1 = (CutVariables[12].c_str() + ">245319");
  TCut preCut2 = (CutVariables[13].c_str() + ">81922.9");
  TCut preCut3 = (CutVariables[14].c_str() + ">91397");
  TCut preCut4 = (CutVariables[15].c_str() + ">76464.4");
  TCut preCut5 = (CutVariables[47].c_str() + ">604.51");
*/
  //  TCut preselectionCut = (preCut1 + ":" + preCut2 + ":" + preCut3 + ":" + preCut4 + ":" + preCut5).c_str()
  TCut preselectionCut = ""; //ADCamp>245319&&ADCampU>81922.9&&ADCampV>91397&&ADCampY>76464.4&&WirestdV>604.51";
  // factory->PrepareTrainingAndTestTree( preselectionCut, "<options>" );
  // <options> = "nTrainSignal=5000:nTrainBackground=5000:nTestSignal=4000:nTestBackground=5000"
  factory->PrepareTrainingAndTestTree( preselectionCut, "" );

  factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=10" ); // Book decision Tree method

  // Train, Test and evaluate Boosted Decision Tree
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

  // END TMVA::Factory METHOD

  myfile << "\n" << "Finished writing to file. TrainandTestDecisionTree.cxx completed successfully!\n\n ";
  myfile.close();  
  NodeFile->Write();
  NodeFile->Close(); 
  // if (!gROOT->IsBatch()) TMVAGui( "NodeFile.root" );
}
