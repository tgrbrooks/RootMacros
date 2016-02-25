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


// Root Macro Includes
#include "ReadTree.cc"
#include "TTreeOneDHist.cxx"
#include "TTreeVsEnergyHist.cxx"
#include "TTreeEnergyAndOneD.cxx"

// Main Script
void ShowerPlotter1D(){

  const int NumOfVars = 52;
  string CutVariables[NumOfVars] = {
    "hitNo","hitNoU","hitNoV","hitNoY","TDCstd","TDCstdU","TDCstdV","TDCstdY","TDCiqr","TDCiqrU","TDCiqrV","TDCiqrY",
    "ADCamp","ADCampU","ADCampV","ADCampY","WFint","WFintU","WFintV","WFintY","Meanamp","MeanampU","MeanampV","MeanampY",
    "Meanint","MeanintU","MeanintV","MeanintY","LowDen","LowDenU","LowDenV","LowDenY","HiDen","HiDenU","HiDenV","HiDenY",
    "MeanRMS","MeanRMSU","MeanRMSV","MeanRMSY","MeanMult","MeanMultU","MeanMultV","MeanMultY",
    "Wirestd","WirestdU","WirestdV","WirestdY","Wireiqr","WireiqrU","WireiqrV","WireiqrY" };

  // Create list to contain shower files to be merged
  const int NumShowerFiles = 100;  string ShowerFileNames[NumShowerFiles];
 
  // Define String stream and temporary string to hold the file number
  stringstream ss; string FileNumAsString;
  // Loop over the number of files and allocate the file names to string vectors.
  // Use string stream to get the file number as a string
  for(unsigned int i(0); i < NumShowerFiles; i ++){
    ss << i; ss >> FileNumAsString; // Convert integer i of loop to string
    ShowerFileNames[i] = " ../RootOutputs/ShowerOutputs/ShowerAna_output_" + FileNumAsString + ".root";
    ss.str(""); ss.clear(); // Clear the string stream
  }

  gROOT->ProcessLine(".x lhcbStyle.C"); // Using lhcb Style file
 
  std::ofstream myfile; // Create output file stream with name ShowerPlotter.txt
  myfile.open("ShowerPlotter1D.txt");
  
  // Create TChain to group root files
  TChain ShowerTrees("ch_tree");
 
  // Message to output file
  myfile << "||ShowerPlotter1D.cc|| Begin.\nThe Following " << NumShowerFiles << " Shower Files were merged:\n";  

  // Loop over shower files and add them to TChain
  for(unsigned int i(0); i < NumShowerFiles; i ++){
    ShowerTrees.Add(ShowerFileNames[i].c_str()); 
    myfile << ShowerFileNames[i] << "\n"; // Print File names to outputfile
  }
  ShowerTrees.Merge("allshower.root"); // Merge all Chained trees into a single root file
   
  // Message to output file
  myfile << "\n" << NumShowerFiles << " Shower Files successfully merged:\n";

  // Call function ReadTree with merged root file to read in the tree
  vector<vector<float>> ShowerEvents = ReadTree("allshower.root");
  
  myfile << "\n" << "Number of shower events: "<<ShowerEvents[0].size();

  myfile << "\n\n" << "Now plotting 1D Histograms of all cut variables.";
  for(unsigned int i(0); i < 6; i ++){ // < NumOfVars
    TTreeEnergyAndOneD(ShowerTrees,CutVariables[i]);
  }
  myfile << "\n" << "Plotting 1D Histograms completed.\n\n ";
  myfile << "\n" << "Finished writing to file. ShowerPlotter1D completed successfully!\n\n ";
  myfile.close();  
}
