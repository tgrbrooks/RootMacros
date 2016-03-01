// Analyser/plotter with optimisation for single and double plane cuts

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

// Root Macro Includes
#include "ReadTree.cc"
#include "CalcMean.cc"
#include "CalcStandardDev.cc"

// Function to cut on all of the variables independently
// Different standard deviations for optimisations, rpm, suv, suy, svy for output run
vector<vector<float>> SingleCut(const vector<vector<float>> & vv_ReturnVect, bool isNeut, const vector<float> & MEANS, 
                                const vector<float> & STDEVP1P, const vector<float> & STDEVP1M, const vector<float> & STDEVP2P, 
                                const vector<float> & STDEVP2M, vector<vector<float>> Rpm, vector<vector<float>> suv, 
                                      vector<vector<float>> suy, vector<vector<float>> svy){

  // Total cut
  vector<vector<float>> RemovedType;
  vector<float> EmptyVect; 
  for(unsigned int i(0); i < 9; i ++) EmptyVect.push_back(0);

  for(unsigned int i(0); i < 20; i ++){
    RemovedType.push_back(EmptyVect); // Removed Type After this loop will be of dimensions RemovedType[20][8]

    for(unsigned int j(0); j < vv_ReturnVect[0].size(); j ++){ // LOOP NUMBER DICTATED BY v_Hits.size()
      // Multiply sigma by optimised +/-
      if( vv_ReturnVect[i][j] > (MEANS[i] + STDEVP1P[i]*Rpm[0][i]) || vv_ReturnVect[i][j] < (MEANS[i] - STDEVP1M[i]*Rpm[1][i]) ){ 
        RemovedType[i][0] ++;
        if(isNeut == true){
          for(unsigned int k=0; k < 8; k ++){
            if(vv_ReturnVect[21][j]==k) RemovedType[i][k+1] ++;
          }
        }
      }
    }
  }
  
  for(unsigned int j(0); j < 20; j++) RemovedType.push_back(EmptyVect); // Removed Type After this loop will be of dimensions RemovedType[40][8]
  for(unsigned int p(0); p < 5; p ++){
    // Index == (1+p*4) takes values 1,5,9,13,17. These correspond to vector elements v_UHits, v_TDCstdU, v_TDCiqrU, v_ADCampU and v_WFintU
    // Index == (2+p*4) takes values 2,6,10,14,18. These correspond to vector elements v_VHits, v_TDCstdV, v_TDCiqrV, v_ADCampV and v_WFintV
    // Index == (3+p*4) takes values 3,7,11,15,19. These correspond to vector elements v_YHits, v_TDCstdY, v_TDCiqrY, v_ADCampY and v_WFintY
    unsigned int UPI = (1+p*4); // UPlaneIndex
    unsigned int VPI = (2+p*4); // VPlaneIndex
    unsigned int YPI = (3+p*4); // YPlaneIndex
    for(unsigned int j(0); j < vv_ReturnVect[0].size(); j ++){ // LOOP NUMBER DICTATED BY v_Hits.size()

      // UV Plane Cut -- multiply sigma+ by +/- ratio and optimised no of sigma for uv plane
      if(
        (vv_ReturnVect[UPI][j]>(MEANS[UPI] + STDEVP1P[UPI]*(Rpm[0][UPI]/Rpm[1][UPI])*suv[0][p])||vv_ReturnVect[UPI][j]<(MEANS[UPI] - STDEVP1M[UPI]*suv[0][p])) ||
        (vv_ReturnVect[VPI][j]>(MEANS[VPI] + STDEVP2P[VPI]*(Rpm[0][VPI]/Rpm[1][VPI])*suv[1][p])||vv_ReturnVect[VPI][j]<(MEANS[VPI] - STDEVP2M[VPI]*suv[1][p]))){
        RemovedType[20+p][0] ++; // A Neutrino was removed from UV plane cut
	if(isNeut == true){
          for(unsigned int k=0; k<8;k++){
            if(vv_ReturnVect[21][j]==k){RemovedType[20+p][k+1] ++;} // Find flavour of neutrino
          }
        }
      }
      
      // UY Plane Cut
      if(
        (vv_ReturnVect[UPI][j]>(MEANS[UPI] + STDEVP1P[UPI]*(Rpm[0][UPI]/Rpm[1][UPI])*suy[0][p])||vv_ReturnVect[UPI][j]<(MEANS[UPI] - STDEVP1M[UPI]*suy[0][p])) ||
        (vv_ReturnVect[YPI][j]>(MEANS[YPI] + STDEVP2P[YPI]*(Rpm[0][YPI]/Rpm[1][YPI])*suy[1][p])||vv_ReturnVect[YPI][j]<(MEANS[YPI] - STDEVP2M[YPI]*suy[1][p]))){ 
        RemovedType[25+p][0] ++; // A Neutrino was removed from UY plane cut
        if(isNeut == true){
	  for(unsigned int k=0; k<8;k++){
	    if(vv_ReturnVect[21][j]==k){RemovedType[25+p][k+1] ++;} // Find flavour of neutrino
          }
	}
      }
    
      // VY Plane Cut
      if(
        (vv_ReturnVect[VPI][j]>(MEANS[VPI] + STDEVP1P[VPI]*(Rpm[0][VPI]/Rpm[1][VPI])*svy[0][p])||vv_ReturnVect[VPI][j]<(MEANS[VPI] - STDEVP1M[VPI]*svy[0][p])) ||
        (vv_ReturnVect[YPI][j]>(MEANS[YPI] + STDEVP2P[YPI]*(Rpm[0][YPI]/Rpm[1][YPI])*svy[1][p])||vv_ReturnVect[YPI][j]<(MEANS[YPI] - STDEVP2M[YPI]*svy[1][p]))){
        RemovedType[30+p][0] ++; // A Neutrino was removed from VY plane cut
	if(isNeut == true){
          for(unsigned int k=0; k<8;k++){
	    if(vv_ReturnVect[21][j]==k){RemovedType[30+p][k+1] ++;} // Find flavour of neutrino
          }
	}
      }
     
      // UVY Plane Cut
      if(
        (vv_ReturnVect[UPI][j]>(MEANS[UPI] + STDEVP1P[UPI])||vv_ReturnVect[UPI][j]<(MEANS[UPI] - STDEVP1M[UPI])) ||
        (vv_ReturnVect[VPI][j]>(MEANS[VPI] + STDEVP2P[VPI])||vv_ReturnVect[VPI][j]<(MEANS[VPI] - STDEVP2M[VPI])) ||
        (vv_ReturnVect[YPI][j]>(MEANS[YPI] + STDEVP2P[YPI])||vv_ReturnVect[YPI][j]<(MEANS[YPI] - STDEVP2M[YPI]))){
        RemovedType[35+p][0] ++; // A Neutrino was removed from UVY plane cut
	if(isNeut == true){
          for(unsigned int k=0; k<8;k++){
            if(vv_ReturnVect[21][j]==k){RemovedType[35+p][k+1] ++;} // Find flavour of neutrino
          }
	}
      }
    }
  }
  if(isNeut == true){
    for(unsigned int i=0; i < vv_ReturnVect[21].size(); i++){
      for(unsigned int j=0;j<8;j++){
        if(vv_ReturnVect[21][i] == j) EmptyVect[j] ++; // Loop over all neutrino types. See vector initialisation above.
      }
    }
    RemovedType.push_back(EmptyVect);
  }
  return RemovedType;

} // End of SingleCut

// Main Script
void Analyser(){
  // String list of cut variables, same as histogram names
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
 
  std::ofstream myfile; // Create output file stream with name AnalyserDataOutput.txt
  myfile.open("ShowerCutOptimiserOutput.txt");
  
  // Create two TChains, one for eminus trees and one for BNB data trees
  TChain EminusTrees("ch_tree");
  TChain BNBTrees("ch_tree");
  
  std::cout << "\nBeginning Shower Cut Optimiser main script. Now reading in files (Check Output file)\n";
  
  // Message to output file
  myfile << "||ShowerCutOptimiser.cc|| Begin.\n||ShowerCutOptimiser.cc|| The Folllowing " << NumEminusFiles << " Eminus Files were merged:\n";  

  // Loop over nnbar files and add them to TChain
  for(unsigned int i(0); i < NumEminusFiles; i ++){
    EminusTrees.Add(EminusFileNames[i].c_str()); 
    myfile << EminusFileNames[i] << "\n"; // Print File names to outputfile
  }
  std::cout << "\nNow merging together Eminusfiles into single file allEminus.root.\n";
  EminusTrees.Merge("allEminus.root"); // Merge all Chained trees into a single root file
  
  std::cout << "\nEminus Files successfully read and merged. Now attempting to gather data from chained file.\n";
   
  vector <float> StandardDevs, Means;

  // Use ReadTree to return vector of vectors that contains all the information for all the events for nnbar
  // Only grandmaster luke understands the indexing
  vector<vector<float>> EminusEvents = ReadTree("allEminus.root");
  
  std::cout << "\nData successfully read. Now calculating Standard Deviations and Means.\n";

  // Output message to file
  myfile << "\n\nMean Values and standard Deviations of Merged Eminusfile for different cut variables\n";

  // Loop over data returned in EminusEvents given by ReadTree() to get the Means and StandardDeviations
  float tempfloat(0);
  for(unsigned int i(0); i < EminusEvents.size(); i ++){ 
    tempfloat = CalcStandardDev(EminusEvents[i]);
    StandardDevs.push_back( tempfloat ); 
    Means.push_back( CalcMean(EminusEvents[i]) );
    myfile << "\n" << CutVariables[i] << " : " << Means[i] << "+/-" << StandardDevs[i]; // Output values to the file
  }

  // Message to output file
  myfile << "\n\n||ShowerCutOptimiser.cc|| The Following " << NumBNBFiles << " BNB Files were analysed:\n";
  
  std::cout << "\nNow adding BNB files to chain and merginging into allBNB.root file.\n";

  // Loop over BNB files and add them to TChainBNBTrees
  for(unsigned int i(0); i < NumBNBFiles; i ++){
    BNBTrees.Add(BNBFileNames[i].c_str()); 
    myfile << BNBFileNames[i] << "\n"; // Print File names to outputfile
  }  
  BNBTrees.Merge("allBNB.root");
  
  // Use ReadTree to return vector of vectors that contains all the information for all the events for neutrinos
  // Only grandmaster luke understands the indexing
  vector<vector<float>> BNBEvents = ReadTree("allBNB.root",true); 

  myfile.close();  

  std::cout << "\nHello World!\n";
  StandardDevs.clear(); Means.clear(); EminusEvents.clear(); BNBEvents.clear();
}
