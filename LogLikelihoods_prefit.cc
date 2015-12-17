// C++ Includes
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <math.h>

// Root Includes
#include "THStack.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"

// Function to approximate simulated distributions as discrete probability density functions
vector<vector<Double_t>> GetDiscretePDFs(TChain &nnbarTree, TChain &NeutTree, const string CutVariable){

  // Plot nnbar distribution to temporary histogram
  NeutTree.Draw(CutVariable.c_str(),"","E1HIST",10000000000,0);

  // Get the number of bins and limits of histograms
  Int_t bins = htemp->GetSize();
  Double_t xmax = htemp->GetXaxis()->GetXmax();
  Double_t xmin = htemp->GetXaxis()->GetXmin();
  std::cout<<bins<<"  "<<xmax<<std::endl;
  // Get the bin content for each bin
  vector<Double_t> bincontent;
  for(int i=0; i<bins; i++){
    bincontent.push_back(htemp->GetBinContent(i));
  }

  // Plot neutrino distribution to temp hist
  TCanvas *c2 = new TCanvas("c2","c2",900,600);
  TH1D *htemp2 = new TH1D("htemp2","htemp2",bins-2,xmin,xmax);
  std::ostringstream s;
  s << CutVariable << ">>" << "htemp2";
  nnbarTree.Draw((s.str()).c_str(),"","E1HIST",10000000000,0);

  int bins2 = htemp2->GetSize();
  float xmax2 = htemp2->GetXaxis()->GetXmax();
  std::cout<<bins2<<"  "<<xmax2<<std::endl;
  vector<Double_t> bincontent2, pnnbar, pneut;
  for(int i=0; i<bins2; i++){
    bincontent2.push_back(htemp2->GetBinContent(i));
    // Normalize bins so that total area adds up to one
    pnnbar.push_back(bincontent2[i]/(10000));
    pneut.push_back(bincontent[i]/(10000));
    //std::cout<<pnnbar[i]<<"  "<<pneut[i]<<std::endl;
  }
  
  vector<vector<Double_t>> probabilities;
  probabilities.push_back(pnnbar);
  probabilities.push_back(pneut);
  return probabilities;

}

// Function that uses RooFit to create continuous PDFs

// Function to read in the ttree's from the input files and store the info in vector of vector of floats
vector<vector<float>> ReadTree(string RootFileName, bool isNeut){

  // Initilise float variables, vector of floats and vector of vector of floats
  float TDCstd, TDCstdU,TDCstdV,TDCstdY,TDCiqr,TDCiqrU,TDCiqrV,TDCiqrY,ADCamp,ADCampU,ADCampV,ADCampY,WFint,WFintU,WFintV,WFintY,Energy;
  vector <float> v_TDCstd, v_TDCstdU,v_TDCstdV,v_TDCstdY,v_TDCiqr,v_TDCiqrU,v_TDCiqrV,v_TDCiqrY,v_ADCamp,v_ADCampU,v_ADCampV,v_ADCampY,v_WFint,v_WFintU,v_WFintV,v_WFintY,v_Energ,v_Hits,v_UHits,v_VHits,v_YHits,v_Type,v_Energy;
  vector < vector < float > > vv_ReturnVect;

  UInt_t Hits,UHits,VHits,YHits,Type;

  TFile *TRootFile = new TFile(RootFileName.c_str());  
  TTree *FileTree = (TTree*)TRootFile->Get("ch_tree");
 
  // Set Branch Addresses to floats
  FileTree->SetBranchAddress("hitNo",&Hits);
  FileTree->SetBranchAddress("hitNoU",&UHits);
  FileTree->SetBranchAddress("hitNoV",&VHits);
  FileTree->SetBranchAddress("hitNoY",&YHits);
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
  FileTree->SetBranchAddress("WFint", &WFint);
  FileTree->SetBranchAddress("WFintU", &WFintU);
  FileTree->SetBranchAddress("WFintV", &WFintV);
  FileTree->SetBranchAddress("WFintY", &WFintY);
  if( isNeut == true ){ // Additional data if it's Neutrino
    FileTree->SetBranchAddress("Energy", &Energy);
    FileTree->SetBranchAddress("Type",&Type);
  }
  Long64_t EntryNumber = FileTree->GetEntries();
  for(unsigned int j=0; j!=EntryNumber; ++j){ // Loop over tree entries
    // Fill vectors of floats with each branch address data member.
    FileTree->GetEntry(j);
    v_Hits.push_back(Hits);
    v_UHits.push_back(UHits);
    v_VHits.push_back(VHits);
    v_YHits.push_back(YHits);
    v_TDCstd.push_back(TDCstd);
    v_TDCstdU.push_back(TDCstdU);
    v_TDCstdV.push_back(TDCstdV);
    v_TDCstdY.push_back(TDCstdY);
    v_TDCiqr.push_back(TDCiqr);
    v_TDCiqrU.push_back(TDCiqrU);
    v_TDCiqrV.push_back(TDCiqrV);
    v_TDCiqrY.push_back(TDCiqrY);
    v_ADCamp.push_back(ADCamp);
    v_ADCampU.push_back(ADCampU);
    v_ADCampV.push_back(ADCampV);
    v_ADCampY.push_back(ADCampY);
    v_WFint.push_back(WFint);
    v_WFintU.push_back(WFintU);
    v_WFintV.push_back(WFintV);
    v_WFintY.push_back(WFintY);

    if( isNeut == true ){ // Additional data if it's Neutrino
      v_Energy.push_back(Energy);
      v_Type.push_back(Type);
    }
  }
  // Push back data read in from tree onto vector of vectors. Columns of 2D vector correspond to events 1,2,3,...,N 
  // While rows correspond to the different cut variables.
  vv_ReturnVect.push_back(v_Hits);   vv_ReturnVect.push_back(v_UHits);   vv_ReturnVect.push_back(v_VHits);   vv_ReturnVect.push_back(v_YHits);
  vv_ReturnVect.push_back(v_TDCstd); vv_ReturnVect.push_back(v_TDCstdU); vv_ReturnVect.push_back(v_TDCstdV); vv_ReturnVect.push_back(v_TDCstdY);
  vv_ReturnVect.push_back(v_TDCiqr); vv_ReturnVect.push_back(v_TDCiqrU); vv_ReturnVect.push_back(v_TDCiqrV); vv_ReturnVect.push_back(v_TDCiqrY);
  vv_ReturnVect.push_back(v_ADCamp); vv_ReturnVect.push_back(v_ADCampU); vv_ReturnVect.push_back(v_ADCampV); vv_ReturnVect.push_back(v_ADCampY);
  vv_ReturnVect.push_back(v_WFint);  vv_ReturnVect.push_back(v_WFintU);  vv_ReturnVect.push_back(v_WFintV);  vv_ReturnVect.push_back(v_WFintY);
  if ( isNeut == true ){
    vv_ReturnVect.push_back(v_Energy); vv_ReturnVect.push_back(v_Type);
  }

  return vv_ReturnVect;

} // End of ReadTree

// Main Script
void LogLikelihoods_prefit(){

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
  myfile.open("LogLikelihoods.txt");
  
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
   
  // Message to output file
  myfile << "\n\n||Analyser.cc|| The Folllowing " << NumNeutFiles << " Neutrino Files were analysed:\n";

  // Loop over Neutrino files and add them to TChain NeutTrees
  for(unsigned int i(0); i < NumNeutFiles; i ++){
    NeutTrees.Add(NeutFileNames[i].c_str()); 
    myfile << NeutFileNames[i] << "\n"; // Print File names to outputfile
  }  
  NeutTrees.Merge("allneut.root"); // MAY NOT NEED THIS LINE ANYMORE
  
  myfile.close();  

  // Get discrete or continuous pdfs
  vector<vector<Double_t>> probabilities = GetDiscretePDFs(nnbarTrees,NeutTrees,CutVariables[13]);
  vector<Double_t> pnnbar = probabilities[0];
  vector<Double_t> pneut = probabilities[1];
  vector<vector<Double_t>> probabilities1 = GetDiscretePDFs(nnbarTrees,NeutTrees,CutVariables[14]);
  vector<Double_t> pnnbar1 = probabilities1[0];
  vector<Double_t> pneut1 = probabilities1[1];
  vector<vector<Double_t>> probabilities2 = GetDiscretePDFs(nnbarTrees,NeutTrees,CutVariables[15]);
  vector<Double_t> pnnbar2 = probabilities2[0];
  vector<Double_t> pneut2 = probabilities2[1];

  // Get data from files
  vector<vector<float>> NnbarEvents = ReadTree("allnnbar.root",false);
  vector<vector<float>> NeutEvents = ReadTree("allneut.root",true);
  
  Double_t _ll_nnbar, _ll_neut;
/*
  // Create tree for plotting/saving etc
  TTree* _t_ch = new TTree("ch_tree","");
  TTree* _t_ch2 = new TTree("ch_tree2","");
  _t_ch->Branch("ll_nnbar",&_ll_nnbar,"ll_nnbar/D");
  _t_ch2->Branch("ll_neut",&_ll_neut,"ll_neut/D");
  _t_ch->SetDirectory(0);
  _t_ch2->SetDirectory(0);
*/
  TH1D *h_NNBAR = new TH1D("h_NNBAR","h_NNBAR",80,-7,10);
  TH1D *h_NEUT = new TH1D("h_NEUT","h_NEUT",80,-7,10);

  int bin, bin2, bin3, bin4, bin5, bin6;
  //vector<Double_t> _ll_neut, _ll_nnbar; 
  // Loop over all events
  for(int k=0; k<NnbarEvents[13].size(); k++){
    _ll_nnbar = 0;
    _ll_neut = 0;
    // Size of each bin - should be calculated automatically
    float binSize = 470/102;
    // Loop over all bins
    for(int j=0; j<102; j++){
      float min = j*binSize;
      float max = min + binSize;
      // If the value for the nnbar or neutrino events corresponds to a bin record that number
      if(NnbarEvents[13][k]>=min&&NnbarEvents[13][k]<max){bin = j;}
      if(NeutEvents[13][k]>=min&&NeutEvents[13][k]<max){bin2 = j;}
    }
    // NOT WORKING PROPERLY - Issue with discrete probabilities often being 0 giving nans
    // Calculate the log likelihood for nnbar to be nnbar rather than neutrino
    if(pnnbar[bin]!=0&&pneut[bin]!=0){_ll_nnbar += log(pnnbar[bin]/(pneut[bin]));}
    // Calculate log likelihood for neutrino to be neutrino rather than nnbar
    if(pnnbar[bin2]!=0&&pneut[bin2]!=0){_ll_neut += log(pneut[bin2]/(pnnbar[bin2]));}

    float binSize = 550/102;
    for(int j=0; j<102; j++){
      float min = j*binSize;
      float max = min + binSize;
      // If the value for the nnbar or neutrino events corresponds to a bin record that number
      if(NnbarEvents[14][k]>=min&&NnbarEvents[14][k]<max){bin3 = j;}
      if(NeutEvents[14][k]>=min&&NeutEvents[14][k]<max){bin4 = j;}
    }
    // Calculate the log likelihood for nnbar to be nnbar rather than neutrino
    if(pnnbar1[bin3]!=0&&pneut1[bin3]!=0){_ll_nnbar += log(pnnbar1[bin3]/(pneut1[bin3]));}
    // Calculate log likelihood for neutrino to be neutrino rather than nnbar
    if(pnnbar1[bin4]!=0&&pneut1[bin4]!=0){_ll_neut += log(pneut1[bin4]/(pnnbar1[bin4]));}
    
    float binSize = 610/102;
    for(int j=0; j<102; j++){
      float min = j*binSize;
      float max = min + binSize;
      // If the value for the nnbar or neutrino events corresponds to a bin record that number
      if(NnbarEvents[15][k]>=min&&NnbarEvents[15][k]<max){bin5 = j;}
      if(NeutEvents[15][k]>=min&&NeutEvents[15][k]<max){bin6 = j;}
    }
    // Calculate the log likelihood for nnbar to be nnbar rather than neutrino
    if(pnnbar2[bin5]!=0&&pneut2[bin5]!=0){_ll_nnbar += log(pnnbar2[bin5]/(pneut2[bin5]));}
    // Calculate log likelihood for neutrino to be neutrino rather than nnbar
    if(pnnbar2[bin6]!=0&&pneut2[bin6]!=0){_ll_neut += log(pneut2[bin6]/(pnnbar2[bin6]));}

    if(_ll_neut!=0){h_NEUT->Fill(_ll_neut);}
    if(_ll_nnbar!=0){h_NNBAR->Fill(_ll_nnbar);}

  }

/*  // Write to root file
  TFile *TRootFile = new TFile("out.root","RECREATE"); 
  _t_ch->Write();*/

  // Draw Histograms
  TCanvas *c3 = new TCanvas("c3","c3",900,600);
  h_NEUT->SetLineColor(kRed);
  h_NNBAR->SetLineColor(kBlue);
  h_NNBAR->GetXaxis()->SetTitle("Log Likelihood");
  h_NNBAR->GetYaxis()->SetTitle("Events");
  h_NNBAR->Draw();
  h_NEUT->Draw("same");
  leg = new TLegend(0.3,0.3,0.5,0.5);
  TLegendEntry *le1 = leg->AddEntry("nnbarTree","n#bar{n}","l"); le1->SetLineColor(kBlue); le1->SetLineWidth(2);
  TLegendEntry *le2 = leg->AddEntry("NeutTree","Atmospheric","l"); le2->SetLineColor(kRed); le2->SetLineWidth(2);
  leg->Draw();
  /*_t_ch2->SetLineColor(kRed);
  _t_ch->SetLineColor(kBlue);
  _t_ch2->Draw("ll_neut>>h","","E1HIST",10000000000,0);
  //TCanvas *c4 = new TCanvas("c4","c4",900,600);
  _t_ch->Draw("ll_nnbar>>h","","E1HISTsame",10000000000,0);
  h->Draw();*/

}
