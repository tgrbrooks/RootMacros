#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooGaussian.h"
#include "RooConstVar.h"

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
#include "TMath.h"

using namespace RooFit;

// ====================== MISCELLANEOUS FUNCTIONS ============================ //

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

// ============================= RooFit FITTING FUNCTIONS ============================= //

// Function to fit a standard Gaussian PDF
vector<Double_t> fitGauss(RooRealVar variable, TTree *tree, Double_t meanGuess, Double_t sigGuess){

  TCanvas *c4 = new TCanvas("c4","",600,400);
  RooDataSet data("data","data",RooArgSet(variable),Import(*tree));
  RooRealVar mean("mean","mean",meanGuess,0,1000000);
  RooRealVar sigma("sigma","sigma",sigGuess,0.1,5000000);
  RooGaussian gauss("gauss","gauss",variable,mean,sigma);
  gauss.fitTo(data);
  
  RooPlot* frame1 = variable.frame(Title("Imported TH1 with Poisson error bars"));
  data.plotOn(frame1);
  gauss.plotOn(frame1);
  Double_t chi2 = 0;
  chi2 = frame1->chiSquare();
  frame1->Draw();

  std::cout<<"Gaussian chi^2 = "<<chi2<<std::endl;

  vector<Double_t> returnVec;
  returnVec.push_back(mean.getVal());
  returnVec.push_back(sigma.getVal());
  return returnVec;

} // End of fitGauss


// Function to fit a Landau PDF
vector<Double_t> fitLand(RooRealVar variable, TTree *tree, Double_t meanGuess, Double_t sigGuess){

  TCanvas *c5 = new TCanvas("c5","",600,400);
  RooDataSet data("data","data",RooArgSet(variable),Import(*tree));
  RooRealVar mean("mean","mean",meanGuess,10,10000);
  RooRealVar sigma("sigma","sigma",sigGuess,0.1,1000);
  RooLandau land("land","land",variable,mean,sigma);
  land.fitTo(data);

  RooPlot* frame2 = variable.frame(Title("Imported TH1 with Poisson error bars"));
  data.plotOn(frame2);
  land.plotOn(frame2);
  Double_t neutchi2 = 0;
  neutchi2 = frame2->chiSquare();
  frame2->Draw();
  std::cout<<"Landau chi^2 = "<<neutchi2<<std::endl;

  vector<Double_t> returnVec;
  returnVec.push_back(mean.getVal());
  returnVec.push_back(sigma.getVal());
  return returnVec;
} // End of fitLand

// Function to fit an exponential PDF with a negative coefficient
Double_t fitExp(RooRealVar variable, TTree *tree, Double_t cGuess){

  TCanvas *c6 = new TCanvas("c6","",600,400);
  RooDataSet data("data","data",RooArgSet(variable),Import(*tree));
  RooRealVar c("c","c",cGuess,-1.,0);
  RooExponential exp("exp","exp",variable,c);
  //exp.fitTo(data);

  RooPlot* frame3 = variable.frame();
  data.plotOn(frame3);
  //exp.plotOn(frame3);
  //Double_t chi2 = 0;
  //chi2 = frame3->chiSquare();
  frame3->Draw();

  //std::cout<<"Exponential Chi^2 = "<<chi2<<std::endl;

  return c.getVal();
}// End of fitExp

// Function to fit an exponential PDF with a negative coefficient
vector<Double_t> fitExpandGauss(RooRealVar variable, TTree *tree, Double_t cGuess, Double_t meanGuess, Double_t sigGuess){

  TCanvas *c6 = new TCanvas("c6","",600,400);
  RooDataSet data("data","data",RooArgSet(variable),Import(*tree));
  RooRealVar c("c","c",cGuess,-1.,0);
  RooExponential exp("exp","exp",variable,c);
  RooRealVar mean("mean","mean",meanGuess,-1000,1000);
  RooRealVar sigma("sigma","sigma",sigGuess,0,100);
  RooGaussian gauss("gauss","gauss",variable,mean,sigma);
  RooRealVar nsig("nsig","#signal events",150,0.,100000);
  RooRealVar nbkg("nbkg","#background events",1000,0.,100000);
  RooAddPdf model("model","g+a",RooArgList(gauss,exp),RooArgList(nsig,nbkg));
  model.fitTo(data);

  RooPlot* frame3 = variable.frame();
  data.plotOn(frame3);
  model.plotOn(frame3);
  model.plotOn(frame3,Components(gauss),LineStyle(kDashed));
  model.plotOn(frame3,Components(exp),LineStyle(kDashed));
  Double_t chi2 = 0;
  chi2 = frame3->chiSquare();
  frame3->Draw();

  std::cout<<"Exponential Chi^2 = "<<chi2<<std::endl;

  vector<Double_t> returnvec;
  returnvec.push_back(c.getVal());
  returnvec.push_back(mean.getVal());
  returnvec.push_back(sigma.getVal());

  return returnvec;
}// End of fitExpandGauss

// Function to fit a bifurcated Gaussian PDF
vector<Double_t> fitBifurGauss(RooRealVar variable, TTree *tree, Double_t meanGuess, Double_t sigLGuess, Double_t sigRGuess){

  TCanvas *c7 = new TCanvas("c7","",600,400);
  RooDataSet data("data","data",RooArgSet(variable),Import(*tree));
  RooRealVar mean("mean","mean",meanGuess,0,1000000);
  RooRealVar sigmaL("sigmaL","sigmaL",sigLGuess,0.1,5000000);
  RooRealVar sigmaR("sigmaR","sigmaR",sigRGuess,0.1,5000000);
  RooBifurGauss bgauss("bgauss","bgauss",variable,mean,sigmaL,sigmaR);
  bgauss.fitTo(data);
  
  RooPlot* frame1 = variable.frame(Title("Imported TH1 with Poisson error bars"));
  data.plotOn(frame1);
  bgauss.plotOn(frame1);
  Double_t chi2 = 0;
  chi2 = frame1->chiSquare();
  frame1->Draw();

  std::cout<<"Bifurcated Gaussian chi^2 = "<<chi2<<std::endl;

  vector<Double_t> returnVec;
  returnVec.push_back(mean.getVal());
  returnVec.push_back(sigmaL.getVal());
  returnVec.push_back(sigmaR.getVal());
  return returnVec;
}// End of fitBifurGauss

// ====================== LOG LIKELIHOOD CALCULATION FUCTIONS ============================= //

// Function to return the value of a normalised bifurcated Gaussian PDF - OBSOLETE
Double_t EvaluateBifurGauss(Double_t x, Double_t mean, Double_t sigmaL, Double_t sigmaR){

  Double_t arg = x - mean;

  Double_t coef(0.0);
  Double_t norm(0.0);

  if (arg < 0.0){
    if (TMath::Abs(sigmaL) > 1e-30) {
      coef = -0.5/(sigmaL*sigmaL);
      norm = 1/(TMath::Sqrt(TMath::TwoPi())*(sigmaL));
    }
  } else {
    if (TMath::Abs(sigmaR) > 1e-30) {
      coef = -0.5/(sigmaR*sigmaR);
      norm = 1/(TMath::Sqrt(TMath::TwoPi())*(sigmaR));
    }
  }
    
  return norm*exp(coef*arg*arg);

}

// Returns single value from bifurcated Gaussian PDF - ONLY USE THIS METHOD IF CAN'T IMPLEMENT MANUALLY - VERY SLOW
Double_t BifurGauss(Double_t x, Double_t mean, Double_t sigmaL, Double_t sigmaR){

  RooBifurGauss bgauss("bgauss","bgauss",RooConst(x),RooConst(mean),RooConst(sigmaL),RooConst(sigmaR));
  return bgauss.getVal();

}

// ========================== MAIN SCRIPT ============================ //

void LogLikelihoods(){

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

  // Open the merged files and retrieve the TTrees
  TFile *TRootFileNnbar = new TFile("allnnbar.root");  
  TTree *NnbarTree = (TTree*)TRootFileNnbar->Get("ch_tree");
  TFile *TRootFileNeut = new TFile("allneut.root");  
  TTree *NeutTree = (TTree*)TRootFileNeut->Get("ch_tree");

  // Create RooFit variable corresponding to all of the cut variables with the range over which to fit over
  RooRealVar hitNo("hitNo","hitNo",0,3000); // nnbar gaussian, neut half gaus/landau + peak at 0
  RooRealVar hitNoU("hitNoU","hitNoU",5,1000); // nnbar gaus/landau, neut exp + peak at 0
  RooRealVar hitNoV("hitNoV","hitNoV",5,1000); // same as above
  RooRealVar hitNoY("hitNoY","hitNoY",5,1500); // nnbar gaussian, neut exp
  RooRealVar TDCstd("TDCstd","TDCstd",20,1000); // nnbar gaussian, neut polynomial + peak at 0
  RooRealVar TDCstdU("TDCstdU","TDCstdU",20,1000); // same as above
  RooRealVar TDCstdV("TDCstdV","TDCstdV",5,1000); // same as above
  RooRealVar TDCstdY("TDCstdY","TDCstdY",10,1000); // nnbar gaussian, neut half landau + peak at 0
  RooRealVar TDCiqr("TDCiqr","TDCiqr",15,1500); // nnbar short gaussian, neut exp + peak at 0
  RooRealVar TDCiqrU("TDCiqrU","TDCiqrU",15,1500); // nnbar short gaussian, neut exp
  RooRealVar TDCiqrV("TDCiqrV","TDCiqrV",10,1500); // same as above
  RooRealVar TDCiqrY("TDCiqrY","TDCiqrY",10,1500); // same as above
  RooRealVar ADCamp("ADCamp","ADCamp",5,400); // nnbar gaussian, neut landau 
  RooRealVar ADCampU("ADCampU","ADCampU",5,400); // same as above
  RooRealVar ADCampV("ADCampV","ADCampV",5,400); // same as above
  RooRealVar ADCampY("ADCampY","ADCampY",5,400); // same as above
  RooRealVar WFint("WFint","WFint",5,1000000); // nnbar gauss, neut exp
  RooRealVar WFintU("WFintU","WFintU",5,500000); // same as above
  RooRealVar WFintV("WFintV","WFintV",5,500000); // same as above
  RooRealVar WFintY("WFintY","WFintY",5,750000); //same as above

  // Fit all of the hit number distributions
  //vector<Double_t> hitNoNnbarGaussVals = fitGauss(hitNo,NnbarTree,1385,283);
  /*Double_t hitNoNeutExpVal = fitExp(hitNo,NeutTree,-0.002);//,1,10);
  vector<Double_t> hitNoUNnbarBifurGaussVals = fitBifurGauss(hitNoU,NnbarTree,413,74,158);
  Double_t hitNoUNeutExpVal = fitExp(hitNoU,NeutTree,-0.006);
  vector<Double_t> hitNoVNnbarBifurGaussVals = fitBifurGauss(hitNoV,NnbarTree,341,70,146);
  Double_t hitNoVNeutExpVal = fitExp(hitNoV,NeutTree,-0.007);
  vector<Double_t> hitNoYNnbarBifurGaussVals = fitBifurGauss(hitNoY,NnbarTree,427,74,173);
  Double_t hitNoYNeutExpVal = fitExp(hitNoY,NeutTree,-0.005);

  // Fit all of the TDC standard deviation distributions
  vector<Double_t> TDCstdNnbarGaussVals = fitGauss(TDCstd,NnbarTree,1000,500);
  Double_t TDCstdNeutExpVal = fitExp(TDCstd,NeutTree,-0.000006);
  vector<Double_t> TDCstdUNnbarGaussVals = fitGauss(TDCstdU,NnbarTree,1000,500);
  Double_t TDCstdUNeutExpVal = fitExp(TDCstdU,NeutTree,-0.000006);
  vector<Double_t> TDCstdVNnbarGaussVals = fitGauss(TDCstdV,NnbarTree,1000,500);
  Double_t TDCstdVNeutExpVal = fitExp(TDCstdV,NeutTree,-0.000006);
  vector<Double_t> TDCstdYNnbarGaussVals = fitGauss(TDCstdY,NnbarTree,1000,500);
  Double_t TDCstdYNeutExpVal = fitExp(TDCstdY,NeutTree,-0.000006);

  // Fit all of the TDC interquartile range distributions
  vector<Double_t> TDCiqrNnbarBifurGaussVals = fitBifurGauss(TDCiqr,NnbarTree,1000,250,250);
  Double_t TDCiqrNeutExpVal = fitExp(TDCiqr,NeutTree,-0.000006);
  vector<Double_t> TDCiqrUNnbarBifurGaussVals = fitBifurGauss(TDCiqrU,NnbarTree,1000,250,250);
  Double_t TDCiqrUNeutExpVal = fitExp(TDCiqrU,NeutTree,-0.000006);
  vector<Double_t> TDCiqrVNnbarBifurGaussVals = fitBifurGauss(TDCiqrV,NnbarTree,1000,250,250);
  Double_t TDCiqrVNeutExpVal = fitExp(TDCiqrV,NeutTree,-0.000006);
  vector<Double_t> TDCiqrYNnbarBifurGaussVals = fitBifurGauss(TDCiqrY,NnbarTree,1000,250,250);
  Double_t TDCiqrYNeutExpVal = fitExp(TDCiqrY,NeutTree,-0.000006);*/

  // Fit all of the mean ADC amplitude distributions
  vector<Double_t> ADCampGaussvals = fitGauss(ADCamp,NnbarTree,25,500);
  vector<Double_t> ADCampLandvals = fitLand(ADCamp,NeutTree,30,500);
  vector<Double_t> ADCampUGaussvals = fitGauss(ADCampU,NnbarTree,25,500);
  vector<Double_t> ADCampULandvals = fitLand(ADCampU,NeutTree,30,500);
  vector<Double_t> ADCampVGaussvals = fitGauss(ADCampV,NnbarTree,25,500);
  vector<Double_t> ADCampVLandvals = fitLand(ADCampV,NeutTree,30,500);
  vector<Double_t> ADCampYGaussvals = fitGauss(ADCampY,NnbarTree,25,500);
  vector<Double_t> ADCampYLandvals = fitLand(ADCampY,NeutTree,30,500);

  // Fit all of the integrated waveform distributions
  /*Double_t WFintNeutExpVal = fitExp(WFint,NeutTree,-0.000006);
  vector<Double_t> WFintNnbarGaussVals = fitGauss(WFint,NnbarTree,550000,100000);
  Double_t WFintUNeutExpVal = fitExp(WFintU,NeutTree,-0.000006);
  vector<Double_t> WFintUNnbarGaussVals = fitGauss(WFintU,NnbarTree,200000,100000);
  Double_t WFintVNeutExpVal = fitExp(WFintV,NeutTree,-0.000006);
  vector<Double_t> WFintVNnbarGaussVals = fitGauss(WFintV,NnbarTree,200000,100000);
  Double_t WFintYNeutExpVal = fitExp(WFintY,NeutTree,-0.000006);
  vector<Double_t> WFintYNnbarGaussVals = fitGauss(WFintY,NnbarTree,200000,100000);*/

  // Get data from files
  vector<vector<float>> NnbarEvents = ReadTree("allnnbar.root",false);
  vector<vector<float>> NeutEvents = ReadTree("allneut.root",true);
  
  // Initialise the log likelihoods
  Double_t _ll_nnbar, _ll_neut;

  // Create histograms to plot the log likelihoods on
  TH1D *h_NNBAR = new TH1D("h_NNBAR","h_NNBAR",100,-30,10);
  TH1D *h_NEUT = new TH1D("h_NEUT","h_NEUT",100,-30,10);
  
  // Initialise removed counters
  int removedneut = 0;
  int removednnbar = 0; 

  // Loop over all of the events (10000) - MAKE SURE FIRST LL IS = INSTEAD OF +=
  for(int k=0; k<NnbarEvents[13].size(); k++){
/*
    // U plane hit number log likelihoods
    Double_t pnnbarnnbarhnu = BifurGauss(NnbarEvents[1][k],hitNoUNnbarBifurGaussVals[0],hitNoUNnbarBifurGaussVals[1],hitNoUNnbarBifurGaussVals[2]);
    Double_t pnnbarneuthnu = TMath::Exp(NnbarEvents[1][k]*hitNoUNeutExpVal);
    Double_t pneutneuthnu = TMath::Exp(NeutEvents[1][k]*hitNoUNeutExpVal);
    Double_t pneutnnbarhnu = BifurGauss(NeutEvents[1][k],hitNoUNnbarBifurGaussVals[0],hitNoUNnbarBifurGaussVals[1],hitNoUNnbarBifurGaussVals[2]);
    _ll_nnbar = log(pnnbarnnbarhnu/pnnbarneuthnu);
    _ll_neut = log(pneutnnbarhnu/pneutneuthnu);

    // V plane hit number log likelihoods
    Double_t pnnbarnnbarhnv = BifurGauss(NnbarEvents[2][k],hitNoVNnbarBifurGaussVals[0],hitNoVNnbarBifurGaussVals[1],hitNoVNnbarBifurGaussVals[2]);
    Double_t pnnbarneuthnv = TMath::Exp(NnbarEvents[2][k]*hitNoUNeutExpVal);
    Double_t pneutneuthnv = TMath::Exp(NeutEvents[2][k]*hitNoUNeutExpVal);
    Double_t pneutnnbarhnv = BifurGauss(NeutEvents[2][k],hitNoVNnbarBifurGaussVals[0],hitNoVNnbarBifurGaussVals[1],hitNoVNnbarBifurGaussVals[2]);
    _ll_nnbar += log(pnnbarnnbarhnv/pnnbarneuthnv);
    _ll_neut += log(pneutnnbarhnv/pneutneuthnv);

    // V plane hit number log likelihoods
    Double_t pnnbarnnbarhny = BifurGauss(NnbarEvents[3][k],hitNoYNnbarBifurGaussVals[0],hitNoYNnbarBifurGaussVals[1],hitNoYNnbarBifurGaussVals[2]);
    Double_t pnnbarneuthny = TMath::Exp(NnbarEvents[3][k]*hitNoUNeutExpVal);
    Double_t pneutneuthny = TMath::Exp(NeutEvents[3][k]*hitNoUNeutExpVal);
    Double_t pneutnnbarhny = BifurGauss(NnbarEvents[3][k],hitNoYNnbarBifurGaussVals[0],hitNoYNnbarBifurGaussVals[1],hitNoYNnbarBifurGaussVals[2]);
    _ll_nnbar += log(pnnbarnnbarhny/pnnbarneuthny);
    _ll_neut += log(pneutnnbarhny/pneutneuthny);

  
    // U plane TDC standard deviation log likelihoods
    Double_t pnnbarnnbartsu = TMath::Gaus(NnbarEvents[5][k], TDCstdUNnbarGaussVals[0], TDCstdUNnbarGaussVals[1]);
    Double_t pnnbarneuttsu = TMath::Exp(NnbarEvents[5][k]*TDCstdUNeutExpVal);
    Double_t pneutneuttsu = TMath::Exp(NeutEvents[5][k]*TDCstdUNeutExpVal);
    Double_t pneutnnbartsu = TMath::Gaus(NeutEvents[5][k], TDCstdUNnbarGaussVals[0], TDCstdUNnbarGaussVals[1]);
    _ll_nnbar += log(pnnbarnnbartsu/pnnbarneuttsu);
    _ll_neut += log(pneutnnbartsu/pneutneuttsu);

    // V plane TDC standard deviation log likelihoods
    Double_t pnnbarnnbartsv = TMath::Gaus(NnbarEvents[6][k], TDCstdVNnbarGaussVals[0], TDCstdVNnbarGaussVals[1]);
    Double_t pnnbarneuttsv = TMath::Exp(NnbarEvents[6][k]*TDCstdVNeutExpVal);
    Double_t pneutneuttsv = TMath::Exp(NeutEvents[6][k]*TDCstdVNeutExpVal);
    Double_t pneutnnbartsv = TMath::Gaus(NeutEvents[6][k], TDCstdVNnbarGaussVals[0], TDCstdVNnbarGaussVals[1]);
    _ll_nnbar += log(pnnbarnnbartsv/pnnbarneuttsv);
    _ll_neut += log(pneutnnbartsv/pneutneuttsv);

    // Y plane TDC standard deviation log likelihoods
    Double_t pnnbarnnbartsy = TMath::Gaus(NnbarEvents[7][k], TDCstdYNnbarGaussVals[0], TDCstdYNnbarGaussVals[1]);
    Double_t pnnbarneuttsy = TMath::Exp(NnbarEvents[7][k]*TDCstdYNeutExpVal);
    Double_t pneutneuttsy = TMath::Exp(NeutEvents[7][k]*TDCstdYNeutExpVal);
    Double_t pneutnnbartsy = TMath::Gaus(NeutEvents[7][k], TDCstdYNnbarGaussVals[0], TDCstdYNnbarGaussVals[1]);
    _ll_nnbar += log(pnnbarnnbartsy/pnnbarneuttsy);
    _ll_neut += log(pneutnnbartsy/pneutneuttsy);


    // U plane TDC interquartile range log likelihoods
    Double_t pnnbarnnbartiu = TMath::Gaus(NnbarEvents[9][k], TDCiqrUNnbarBifurGaussVals[0], TDCiqrUNnbarBifurGaussVals[1]);
    Double_t pnnbarneuttiu = TMath::Exp(NnbarEvents[9][k]*TDCiqrUNeutExpVal);
    Double_t pneutneuttiu = TMath::Exp(NeutEvents[9][k]*TDCiqrUNeutExpVal);
    Double_t pneutnnbartiu = TMath::Gaus(NeutEvents[9][k], TDCiqrUNnbarBifurGaussVals[0], TDCiqrUNnbarBifurGaussVals[1]);
    _ll_nnbar += log(pnnbarnnbartiu/pnnbarneuttiu);
    _ll_neut += log(pneutnnbartiu/pneutneuttiu);

    // V plane TDC interquartile range log likelihoods
    Double_t pnnbarnnbartiv = TMath::Gaus(NnbarEvents[10][k], TDCiqrVNnbarBifurGaussVals[0], TDCiqrVNnbarBifurGaussVals[1]);
    Double_t pnnbarneuttiv = TMath::Exp(NnbarEvents[10][k]*TDCiqrVNeutExpVal);
    Double_t pneutneuttiv = TMath::Exp(NeutEvents[10][k]*TDCiqrVNeutExpVal);
    Double_t pneutnnbartiv = TMath::Gaus(NeutEvents[10][k], TDCiqrVNnbarBifurGaussVals[0], TDCiqrVNnbarBifurGaussVals[1]);
    _ll_nnbar += log(pnnbarnnbartiv/pnnbarneuttiv);
    _ll_neut += log(pneutnnbartiv/pneutneuttiv);

    // Y plane TDC interquartile range log likelihoods
    Double_t pnnbarnnbartiy = TMath::Gaus(NnbarEvents[11][k], TDCiqrYNnbarBifurGaussVals[0], TDCiqrYNnbarBifurGaussVals[1]);
    Double_t pnnbarneuttiy = TMath::Exp(NnbarEvents[11][k]*TDCiqrYNeutExpVal);
    Double_t pneutneuttiy = TMath::Exp(NeutEvents[11][k]*TDCiqrYNeutExpVal);
    Double_t pneutnnbartiy = TMath::Gaus(NeutEvents[11][k], TDCiqrYNnbarBifurGaussVals[0], TDCiqrYNnbarBifurGaussVals[1]);
    _ll_nnbar += log(pnnbarnnbartiy/pnnbarneuttiy);
    _ll_neut += log(pneutnnbartiy/pneutneuttiy);*/


    // U plane mean ADC amplitude log likelihoods
    Double_t pnnbarnnbaraau = TMath::Gaus(NnbarEvents[13][k], ADCampUGaussvals[0], ADCampUGaussvals[1]);
    Double_t pnnbarneutaau = TMath::Landau(NnbarEvents[13][k], ADCampULandvals[0], ADCampULandvals[1]);
    Double_t pneutneutaau = TMath::Landau(NeutEvents[13][k], ADCampULandvals[0], ADCampULandvals[1]);
    Double_t pneutnnbaraau = TMath::Gaus(NeutEvents[13][k], ADCampUGaussvals[0], ADCampUGaussvals[1]);
    _ll_nnbar = log(pnnbarnnbaraau/pnnbarneutaau);
    _ll_neut = log(pneutnnbaraau/pneutneutaau);

    // V plane mean ADC amplitude log likelihoods
    Double_t pnnbarnnbaraav = TMath::Gaus(NnbarEvents[14][k], ADCampVGaussvals[0], ADCampVGaussvals[1]);
    Double_t pnnbarneutaav = TMath::Landau(NnbarEvents[14][k], ADCampVLandvals[0], ADCampVLandvals[1]);
    Double_t pneutneutaav = TMath::Landau(NeutEvents[14][k], ADCampVLandvals[0], ADCampVLandvals[1]);
    Double_t pneutnnbaraav = TMath::Gaus(NeutEvents[14][k], ADCampVGaussvals[0], ADCampVGaussvals[1]);
    _ll_nnbar += log(pnnbarnnbaraav/pnnbarneutaav);
    _ll_neut += log(pneutnnbaraav/pneutneutaav);

    // Y plane mean ADC amplitude log likelihoods
    Double_t pnnbarnnbaraay = TMath::Gaus(NnbarEvents[15][k], ADCampYGaussvals[0], ADCampYGaussvals[1]);
    Double_t pnnbarneutaay = TMath::Landau(NnbarEvents[15][k], ADCampYLandvals[0], ADCampYLandvals[1]);
    Double_t pneutneutaay = TMath::Landau(NeutEvents[15][k], ADCampYLandvals[0], ADCampYLandvals[1]);
    Double_t pneutnnbaraay = TMath::Gaus(NeutEvents[15][k], ADCampYGaussvals[0], ADCampYGaussvals[1]);
    _ll_nnbar += log(pnnbarnnbaraay/pnnbarneutaay);
    _ll_neut += log(pneutnnbaraay/pneutneutaay);

/*
    // U plane integrated waveform log likelihoods
    Double_t pnnbarnnbarwfu = TMath::Gaus(NnbarEvents[17][k], WFintUNnbarGaussVals[0], WFintUNnbarGaussVals[1]);
    Double_t pnnbarneutwfu = TMath::Exp(NnbarEvents[17][k]*WFintUNeutExpVal);
    Double_t pneutneutwfu = TMath::Exp(NeutEvents[17][k]*WFintUNeutExpVal);
    Double_t pneutnnbarwfu = TMath::Gaus(NeutEvents[17][k], WFintUNnbarGaussVals[0], WFintUNnbarGaussVals[1]);
    _ll_nnbar += log(pnnbarnnbarwfu/pnnbarneutwfu);
    _ll_neut += log(pneutnnbarwfu/pneutneutwfu);

    // V plane integrated waveform log likelihoods
    Double_t pnnbarnnbarwfv = TMath::Gaus(NnbarEvents[18][k], WFintVNnbarGaussVals[0], WFintVNnbarGaussVals[1]);
    Double_t pnnbarneutwfv = TMath::Exp(NnbarEvents[18][k]*WFintVNeutExpVal);
    Double_t pneutneutwfv = TMath::Exp(NeutEvents[18][k]*WFintVNeutExpVal);
    Double_t pneutnnbarwfv = TMath::Gaus(NeutEvents[18][k], WFintVNnbarGaussVals[0], WFintVNnbarGaussVals[1]);
    _ll_nnbar += log(pnnbarnnbarwfv/pnnbarneutwfv);
    _ll_neut += log(pneutnnbarwfv/pneutneutwfv);

    // Y plane integrated waveform log likelihoods
    Double_t pnnbarnnbarwfy = TMath::Gaus(NnbarEvents[19][k], WFintUNnbarGaussVals[0], WFintUNnbarGaussVals[1]);
    Double_t pnnbarneutwfy = TMath::Exp(NnbarEvents[19][k]*WFintYNeutExpVal);
    Double_t pneutneutwfy = TMath::Exp(NeutEvents[19][k]*WFintUNeutExpVal);
    Double_t pneutnnbarwfy = TMath::Gaus(NeutEvents[19][k], WFintUNnbarGaussVals[0], WFintUNnbarGaussVals[1]);
    _ll_nnbar += log(pnnbarnnbarwfy/(pnnbarneutwfy));
    _ll_neut += log(pneutnnbarwfy/pneutneutwfy);*/

    // Apply cuts on the log likelihoods - TO BE REPLACED BY LUKES OPTIMISER
    if(_ll_nnbar<1||_ll_nnbar>6){removednnbar += 1;}
    if(_ll_neut<1||_ll_neut>6){removedneut += 1;}
    // Fill histograms
    h_NNBAR->Fill(_ll_nnbar);
    h_NEUT->Fill(_ll_neut);
  }

std::cout<<removednnbar<<"  "<<removedneut<<std::endl;

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

} // End of LogLikelihoods
