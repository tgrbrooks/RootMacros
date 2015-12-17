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

// Function to remove TTree entries with zero hit number
void RemoveZeros(string RootFileName){

  // Get all data from old root file
  float TDCstd, TDCstdU,TDCstdV,TDCstdY,TDCiqr,TDCiqrU,TDCiqrV,TDCiqrY,ADCamp,ADCampU,ADCampV,ADCampY,WFint,WFintU,WFintV,WFintY,Energy;
  UInt_t Hits,UHits,VHits,YHits,Type;
  TFile *TRootFile = new TFile(RootFileName.c_str());
  TTree *FileTree = (TTree*)TRootFile->Get("ch_tree");
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
  FileTree->SetBranchAddress("Energy", &Energy);
  FileTree->SetBranchAddress("Type", &Type);
  Long64_t EntryNumber = FileTree->GetEntries();

  // Create new TTree - no way to just remove entries from old
  TTree* _t_ch;
  UInt_t _hitNo, _hitNoU, _hitNoV, _hitNoY, _Type; 
  float _TDCstd, _TDCstdU, _TDCstdV, _TDCstdY;
  float _TDCiqr, _TDCiqrU, _TDCiqrV, _TDCiqrY;
  float _ADCamp, _ADCampU, _ADCampV, _ADCampY;
  float _WFint, _WFintU, _WFintV, _WFintY, _Energy;

  _t_ch = new TTree("ch_tree","");
  _t_ch->Branch("hitNo",&_hitNo,"hitNo/i");
  _t_ch->Branch("hitNoU",&_hitNoU,"hitNoU/i");
  _t_ch->Branch("hitNoV",&_hitNoV,"hitNoV/i");
  _t_ch->Branch("hitNoY",&_hitNoY,"hitNoY/i");
  _t_ch->Branch("TDCstd",&_TDCstd,"TDCstd/F");
  _t_ch->Branch("TDCstdU",&_TDCstdU,"TDCstdU/F");
  _t_ch->Branch("TDCstdV",&_TDCstdV,"TDCstdV/F");
  _t_ch->Branch("TDCstdY",&_TDCstdY,"TDCstdY/F");
  _t_ch->Branch("TDCiqr",&_TDCiqr,"TDCiqr/F");
  _t_ch->Branch("TDCiqrU",&_TDCiqrU,"TDCiqrU/F");
  _t_ch->Branch("TDCiqrV",&_TDCiqrV,"TDCiqrV/F");
  _t_ch->Branch("TDCiqrY",&_TDCiqrY,"TDCiqrY/F");
  _t_ch->Branch("ADCamp",&_ADCamp,"ADCamp/F");
  _t_ch->Branch("ADCampU",&_ADCampU,"ADCampU/F");
  _t_ch->Branch("ADCampV",&_ADCampV,"ADCampV/F");
  _t_ch->Branch("ADCampY",&_ADCampY,"ADCampY/F");
  _t_ch->Branch("WFint",&_WFint,"WFint/F");
  _t_ch->Branch("WFintU",&_WFintU,"WFintU/F");
  _t_ch->Branch("WFintV",&_WFintV,"WFintV/F");
  _t_ch->Branch("WFintY",&_WFintY,"WFintY/F");
  _t_ch->Branch("Energy",&_Energy,"Energy/F");
  _t_ch->Branch("Type",&_Type,"Type/i");
  _t_ch->SetDirectory(0);

  for(unsigned int j=0; j!=EntryNumber; ++j){ // Loop over tree entries
    // Fill vectors of floats with each branch address data member.
    FileTree->GetEntry(j);
    if(Hits==0){}//Remove entry
    else{
      _hitNo = Hits;
      _hitNoU = UHits;
      _hitNoV = VHits;
      _hitNoY = YHits;
      _TDCstd = TDCstd;
      _TDCstdU = TDCstdU;
      _TDCstdV = TDCstdV;
      _TDCstdY = TDCstdY;
      _TDCiqr = TDCiqr;
      _TDCiqrU = TDCiqrU;
      _TDCiqrV = TDCiqrV;
      _TDCiqrY = TDCiqrY;
      _ADCamp = ADCamp;
      _ADCampU = ADCampU;
      _ADCampV = ADCampV;
      _ADCampY = ADCampY;
      _WFint = WFint;
      _WFintU = WFintU;
      _WFintV = WFintV;
      _WFintY = WFintY;
      _Energy = Energy;
      _Type = Type;
      _t_ch->Fill();
    }
  }
  // overwrite root file
  TRootFile->Close();
  TFile *OutFile = new TFile(RootFileName.c_str(),"RECREATE");
  _t_ch->Write();
  OutFile->Close();

}// End of RemoveZeros


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

  //TCanvas *c4 = new TCanvas("c4","",600,400);
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
  //frame1->Draw();

  std::cout<<"Gaussian chi^2 = "<<chi2<<std::endl;

  vector<Double_t> returnVec;
  returnVec.push_back(mean.getVal());
  returnVec.push_back(sigma.getVal());
  return returnVec;

} // End of fitGauss


// Function to fit a Landau PDF
vector<Double_t> fitLand(RooRealVar variable, TTree *tree, Double_t meanGuess, Double_t sigGuess){

  //TCanvas *c5 = new TCanvas("c5","",600,400);
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
  //frame2->Draw();
  std::cout<<"Landau chi^2 = "<<neutchi2<<std::endl;

  vector<Double_t> returnVec;
  returnVec.push_back(mean.getVal());
  returnVec.push_back(sigma.getVal());
  return returnVec;
} // End of fitLand

// Function to fit a Landau and Gaussian PDF
vector<Double_t> fitLandandGauss(RooRealVar variable, TTree *tree, Double_t LmeanGuess, Double_t LsigGuess, Double_t GmeanGuess, Double_t GsigGuess){

  //TCanvas *c5 = new TCanvas("c5","",600,400);
  RooDataSet data("data","data",RooArgSet(variable),Import(*tree));
  RooRealVar Lmean("Lmean","Lmean",LmeanGuess,10,10000);
  RooRealVar Lsigma("Lsigma","Lsigma",LsigGuess,0.1,1000);
  RooLandau land("land","land",variable,Lmean,Lsigma);
  RooRealVar Gmean("Gmean","Gmean",GmeanGuess,10,10000);
  RooRealVar Gsigma("Gsigma","Gsigma",GsigGuess,0.1,1000);
  RooGaussian gauss("gauss","gauss",variable,Gmean,Gsigma);
  RooRealVar ngauss("ngauss","#gauss events",150,0.,100000);
  RooRealVar nland("nland","#land events",1000,0.,100000);
  RooAddPdf model("model","g+a",RooArgList(land,gauss),RooArgList(nland,ngauss));
  model.fitTo(data);

  RooPlot* frame2 = variable.frame(Title("Imported TH1 with Poisson error bars"));
  data.plotOn(frame2);
  model.plotOn(frame2);
  model.plotOn(frame2,Components(gauss),LineStyle(kDashed));
  model.plotOn(frame2,Components(land),LineStyle(kDashed));
  Double_t neutchi2 = 0;
  neutchi2 = frame2->chiSquare();
  //frame2->Draw();
  std::cout<<"Landau chi^2 = "<<neutchi2<<std::endl;

  vector<Double_t> returnVec;
  returnVec.push_back(Lmean.getVal());
  returnVec.push_back(Lsigma.getVal());
  returnVec.push_back(Gmean.getVal());
  returnVec.push_back(Gsigma.getVal());
  returnVec.push_back(nland.getVal());
  returnVec.push_back(ngauss.getVal());
  return returnVec;
} // End of fitLandandGauss

// Function to fit an exponential PDF with a negative coefficient
Double_t fitExp(RooRealVar variable, TTree *tree, Double_t cGuess){

  //TCanvas *c6 = new TCanvas("c6","",600,400);
  RooDataSet data("data","data",RooArgSet(variable),Import(*tree));
  RooRealVar c("c","c",cGuess,-1.,0);
  RooExponential exp("exp","exp",variable,c);
  exp.fitTo(data);

  RooPlot* frame3 = variable.frame();
  data.plotOn(frame3);
  exp.plotOn(frame3);
  Double_t chi2 = 0;
  chi2 = frame3->chiSquare();
  //frame3->Draw();

  std::cout<<"Exponential Chi^2 = "<<chi2<<std::endl;

  return c.getVal();
}// End of fitExp

// Function to fit an exponential PDF with a negative coefficient
vector<Double_t> fitExpandGauss(RooRealVar variable, TTree *tree, Double_t cGuess, Double_t meanGuess, Double_t sigGuess){

  //TCanvas *c6 = new TCanvas("c6","",600,400);
  RooDataSet data("data","data",RooArgSet(variable),Import(*tree));
  RooRealVar c("c","c",cGuess,-1.,0);
  RooExponential exp("exp","exp",variable,c);
  RooRealVar mean("mean","mean",meanGuess,0,1000);
  RooRealVar sigma("sigma","sigma",sigGuess,0,100);
  RooGaussian gauss("gauss","gauss",variable,mean,sigma);
  RooRealVar ngauss("ngauss","#gauss events",150,0.,100000);
  RooRealVar nexp("nexp","#exp events",1000,0.,100000);
  RooAddPdf model("model","g+a",RooArgList(gauss,exp),RooArgList(ngauss,nexp));
  model.fitTo(data);

  RooPlot* frame3 = variable.frame();
  data.plotOn(frame3);
  model.plotOn(frame3);
  model.plotOn(frame3,Components(gauss),LineStyle(kDashed));
  model.plotOn(frame3,Components(exp),LineStyle(kDashed));
  Double_t chi2 = 0;
  chi2 = frame3->chiSquare();
  //frame3->Draw();

  std::cout<<"Exponential Chi^2 = "<<chi2<<std::endl;

  vector<Double_t> returnvec;
  returnvec.push_back(c.getVal());
  returnvec.push_back(mean.getVal());
  returnvec.push_back(sigma.getVal());
  returnvec.push_back(ngauss.getVal());
  returnvec.push_back(nexp.getVal());

  return returnvec;
}// End of fitExpandGauss

// Function to fit an exponential PDF with a negative coefficient
vector<Double_t> fitExpandGaussandGauss(RooRealVar variable, TTree *tree, Double_t cGuess, Double_t meanGuess, Double_t sigGuess, Double_t mean2Guess, Double_t sig2Guess){

  //TCanvas *c6 = new TCanvas("c6","",600,400);
  RooDataSet data("data","data",RooArgSet(variable),Import(*tree));
  RooRealVar c("c","c",cGuess,-1.,0);
  RooExponential exp("exp","exp",variable,c);
  RooRealVar mean("mean","mean",meanGuess,0,1000);
  RooRealVar sigma("sigma","sigma",sigGuess,0,100);
  RooGaussian gauss("gauss","gauss",variable,mean,sigma);
  RooRealVar mean2("mean2","mean2",mean2Guess,0,1000);
  RooRealVar sigma2("sigma2","sigma2",sig2Guess,0,100);
  RooGaussian gauss2("gauss2","gauss2",variable,mean2,sigma2);
  RooRealVar ngauss("ngauss","#gauss events",150,0.,100000);
  RooRealVar ngauss2("ngauss2","#gauss2 events",150,0.,100000);
  RooRealVar nexp("nexp","#exp events",1000,0.,100000);
  RooAddPdf model("model","g+a",RooArgList(gauss,gauss2,exp),RooArgList(ngauss,ngauss2,nexp));
  model.fitTo(data);

  RooPlot* frame3 = variable.frame();
  data.plotOn(frame3);
  model.plotOn(frame3);
  model.plotOn(frame3,Components(gauss),LineStyle(kDashed));
  model.plotOn(frame3,Components(gauss2),LineStyle(kDashed));
  model.plotOn(frame3,Components(exp),LineStyle(kDashed));
  Double_t chi2 = 0;
  chi2 = frame3->chiSquare();
  //frame3->Draw();

  std::cout<<"Exponential Chi^2 = "<<chi2<<std::endl;

  vector<Double_t> returnvec;
  returnvec.push_back(c.getVal());
  returnvec.push_back(mean.getVal());
  returnvec.push_back(sigma.getVal());
  returnvec.push_back(ngauss.getVal());
  returnvec.push_back(mean2.getVal());
  returnvec.push_back(sigma2.getVal());
  returnvec.push_back(ngauss2.getVal());
  returnvec.push_back(nexp.getVal());

  return returnvec;
}// End of fitExpandGauss

// Function to fit a bifurcated Gaussian PDF
vector<Double_t> fitBifurGauss(RooRealVar variable, TTree *tree, Double_t meanGuess, Double_t sigLGuess, Double_t sigRGuess){

  //TCanvas *c7 = new TCanvas("c7","",600,400);
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
  //frame1->Draw();

  std::cout<<"Bifurcated Gaussian chi^2 = "<<chi2<<std::endl;

  vector<Double_t> returnVec;
  returnVec.push_back(mean.getVal());
  returnVec.push_back(sigmaL.getVal());
  returnVec.push_back(sigmaR.getVal());
  return returnVec;
}// End of fitBifurGauss

vector<Double_t> fitSixPoly(RooRealVar variable, TTree *tree, Double_t a1Guess, Double_t a2Guess, Double_t a3Guess, Double_t a4Guess, Double_t a5Guess, Double_t a6Guess){

  //TCanvas *c8 = new TCanvas("c8","",600,400);
  RooDataSet data("data","data",RooArgSet(variable),Import(*tree));
  RooRealVar a1("a1","a1",a1Guess,-10000,10000);
  RooRealVar a2("a2","a2",a2Guess,-10000,10000);
  RooRealVar a3("a3","a3",a3Guess,-10000,10000);
  RooRealVar a4("a4","a4",a4Guess,-10000,10000);
  RooRealVar a5("a5","a5",a5Guess,-10000,10000);
  RooRealVar a6("a6","a6",a6Guess,-10000,10000);
  RooPolynomial poly("poly","poly",variable,RooArgList(a1,a2,a3,a4));
  poly.fitTo(data);

  RooPlot* frame1 = variable.frame(Title("Imported TH1 with Poisson error bars"));
  data.plotOn(frame1);
  poly.plotOn(frame1);
  Double_t chi2 = 0;
  chi2 = frame1->chiSquare();
  //frame1->Draw();

  std::cout<<"6-Polynomial chi^2 = "<<chi2<<std::endl;

  vector<Double_t> returnVec;
  returnVec.push_back(a1.getVal());
  returnVec.push_back(a2.getVal());
  returnVec.push_back(a3.getVal());
  returnVec.push_back(a4.getVal());
  returnVec.push_back(a5.getVal());
  returnVec.push_back(a6.getVal());
  return returnVec;

}

// ====================== LOG LIKELIHOOD CALCULATION FUCTIONS ============================= //

// Returns single value from bifurcated Gaussian PDF - ONLY USE THIS METHOD IF CAN'T IMPLEMENT MANUALLY - VERY SLOW
Double_t BifurGauss(Double_t x, Double_t mean, Double_t sigmaL, Double_t sigmaR){

  RooBifurGauss bgauss("bgauss","bgauss",RooConst(x),RooConst(mean),RooConst(sigmaL),RooConst(sigmaR));
  return bgauss.getVal();

}

Double_t ExpandGauss(Double_t x, Double_t c, Double_t mean, Double_t sigma, Double_t ngauss, Double_t nexp){

  Double_t expRatio = nexp/(nexp+ngauss);
  Double_t gaussRatio = ngauss/(nexp+ngauss);
  Double_t Value = expRatio*TMath::Exp(x*c) + gaussRatio*TMath::Gaus(x,mean,sigma);
  return Value;

}

Double_t LandandGauss(Double_t x, Double_t Lmean, Double_t Lsigma, Double_t Gmean, Double_t Gsigma, Double_t nland, Double_t ngauss){

  Double_t landRatio = nland/(nland+ngauss);
  Double_t gaussRatio = ngauss/(nland+ngauss);
  Double_t Value = landRatio*TMath::Landau(x,Lmean,Lsigma) + gaussRatio*TMath::Gaus(x,Gmean,Gsigma);
  return Value;

}

// ========================== MAIN SCRIPT ============================ //

void LogLikelihoods_run5(){

  // List nnbar and Neutrino files to be merged
  const int NumNbarFiles = 100;  string nnbarFileNames[NumNbarFiles];
  const int NumNeutFiles = 100; string NeutFileNames[NumNeutFiles];
  
  // Define String stream and temporary string to hold the file number
  stringstream ss; string FileNumAsString;
  // Loop over the number of files and allocate the file names to string vectors.
  // Use string stream to get the file number as a string
  for(unsigned int i(0); i < NumNeutFiles; i ++){
    ss << i; ss >> FileNumAsString; // Convert integer i of loop to string
    nnbarFileNames[i] = "outputs5_2/Nnbar_output_run5_" + FileNumAsString + ".root";
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
  RemoveZeros("allnnbar.root");
   
  // Message to output file
  myfile << "\n\n||Analyser.cc|| The Folllowing " << NumNeutFiles << " Neutrino Files were analysed:\n";

  // Loop over Neutrino files and add them to TChain NeutTrees
  for(unsigned int i(0); i < NumNeutFiles; i ++){
    NeutTrees.Add(NeutFileNames[i].c_str()); 
    myfile << NeutFileNames[i] << "\n"; // Print File names to outputfile
  }  
  NeutTrees.Merge("allneut.root"); // MAY NOT NEED THIS LINE ANYMORE
  RemoveZeros("allneut.root");
  
  myfile.close();  

  // Open the merged files and retrieve the TTrees
  TFile *TRootFileNnbar = new TFile("allnnbar.root");  
  TTree *NnbarTree = (TTree*)TRootFileNnbar->Get("ch_tree");
  TFile *TRootFileNeut = new TFile("allneut.root");  
  TTree *NeutTree = (TTree*)TRootFileNeut->Get("ch_tree");

  // Create RooFit variable corresponding to all of the cut variables with the range over which to fit over
  RooRealVar hitNonnbar("hitNo","hitNo",13,2500); // nnbar gaussian, neut half gaus/landau + peak at 0
  RooRealVar hitNoneut("hitNo","hitNo",0,1500); // nnbar gaussian, neut half gaus/landau + peak at 0
  RooRealVar hitNoUnnbar("hitNoU","hitNoU",0,1000); // nnbar gaus/landau, neut exp + peak at 0
  RooRealVar hitNoUneut("hitNoU","hitNoU",0,500); // nnbar gaus/landau, neut exp + peak at 0
  RooRealVar hitNoVnnbar("hitNoV","hitNoV",0,600); // same as above
  RooRealVar hitNoVneut("hitNoV","hitNoV",0,500); // same as above
  RooRealVar hitNoYnnbar("hitNoY","hitNoY",0,1200); // nnbar gaussian, neut exp
  RooRealVar hitNoYneut("hitNoY","hitNoY",0,500); // nnbar gaussian, neut exp

  RooRealVar TDCstd("TDCstd","TDCstd",0,1200); // nnbar gaussian, neut polynomial + peak at 0
  RooRealVar TDCstdU("TDCstdU","TDCstdU",0,1200); // same as above
  RooRealVar TDCstdVnnbar("TDCstdV","TDCstdV",0,1500); // same as above
  RooRealVar TDCstdVneut("TDCstdV","TDCstdV",0,300); // same as above
  RooRealVar TDCstdYnnbar("TDCstdY","TDCstdY",0,1500); // nnbar gaussian, neut half landau + peak at 0
  RooRealVar TDCstdYneut("TDCstdY","TDCstdY",0,500); // same as above

  RooRealVar TDCiqrnnbar("TDCiqr","TDCiqr",0,2500); // nnbar short gaussian, neut exp + peak at 0
  RooRealVar TDCiqrneut("TDCiqr","TDCiqr",0,1000); // nnbar short gaussian, neut exp + peak at 0
  RooRealVar TDCiqrUnnbar("TDCiqrU","TDCiqrU",0,2500); // nnbar short gaussian, neut exp
  RooRealVar TDCiqrUneut("TDCiqrU","TDCiqrU",1,400); // nnbar short gaussian, neut exp
  RooRealVar TDCiqrVnnbar("TDCiqrV","TDCiqrV",0,2500); // same as above
  RooRealVar TDCiqrVneut("TDCiqrV","TDCiqrV",1,1000); // same as above
  RooRealVar TDCiqrYnnbar("TDCiqrY","TDCiqrY",0,2500); // same as above
  RooRealVar TDCiqrYneut("TDCiqrY","TDCiqrY",1,1000); // same as above

  RooRealVar ADCampnnbar("ADCamp","ADCamp",14,38); // nnbar gaussian, neut landau 
  RooRealVar ADCampneut("ADCamp","ADCamp",0,150); // nnbar gaussian, neut landau
  RooRealVar ADCampUnnbar("ADCampU","ADCampU",14,38); // same as above
  RooRealVar ADCampUneut("ADCampU","ADCampU",0,150); // same as above
  RooRealVar ADCampVnnbar("ADCampV","ADCampV",10,38); // same as above
  RooRealVar ADCampVneut("ADCampV","ADCampV",0,150); // same as above
  RooRealVar ADCampYnnbar("ADCampY","ADCampY",0,80); // same as above
  RooRealVar ADCampYneut("ADCampY","ADCampY",0,200); // same as above

  RooRealVar WFintnnbar("WFint","WFint",0,800000); // nnbar gauss, neut exp
  RooRealVar WFintneut("WFint","WFint",0,400000); // nnbar gauss, neut exp
  RooRealVar WFintU("WFintU","WFintU",5,300000); // same as above
  RooRealVar WFintV("WFintV","WFintV",5,150000); // same as above
  RooRealVar WFintY("WFintY","WFintY",5,400000); //same as above

  // Fit all of the hit number distributions
  vector<Double_t> hitNoNnbarGaussVals = fitGauss(hitNonnbar,NnbarTree,1385,283);
  vector<Double_t> hitNoNeutExpVal = fitExpandGauss(hitNoneut,NeutTree,-0.002,5,10);
  vector<Double_t> hitNoUNnbarBifurGaussVals = fitBifurGauss(hitNoUnnbar,NnbarTree,413,74,158);
  vector<Double_t> hitNoUNeutExpVal = fitExpandGauss(hitNoUneut,NeutTree,-0.006,2,10);
  vector<Double_t> hitNoVNnbarBifurGaussVals = fitBifurGauss(hitNoVnnbar,NnbarTree,341,70,146);
  vector<Double_t> hitNoVNeutExpVal = fitExpandGauss(hitNoVneut,NeutTree,-0.007,1,10);
  vector<Double_t> hitNoYNnbarBifurGaussVals = fitBifurGauss(hitNoYnnbar,NnbarTree,427,74,173);
  vector<Double_t> hitNoYNeutExpVal = fitExpandGauss(hitNoYneut,NeutTree,-0.005,5,10);

  // Fit all of the TDC standard deviation distributions
  vector<Double_t> TDCstdNnbarGaussVals = fitGauss(TDCstd,NnbarTree,1000,500);
  vector<Double_t> TDCstdNeutExpVal = fitExpandGauss(TDCstd,NeutTree,-0.005,2,10);
  vector<Double_t> TDCstdUNnbarGaussVals = fitGauss(TDCstdU,NnbarTree,1000,500);
  vector<Double_t> TDCstdUNeutExpVal = fitExpandGauss(TDCstdU,NeutTree,-0.005,2,10);
  vector<Double_t> TDCstdVNnbarGaussVals = fitGauss(TDCstdVnnbar,NnbarTree,1000,500);
  vector<Double_t> TDCstdVNeutExpVal = fitExpandGauss(TDCstdVneut,NeutTree,-0.005,0,10);
  vector<Double_t> TDCstdYNnbarGaussVals = fitGauss(TDCstdYnnbar,NnbarTree,1000,500);
  vector<Double_t> TDCstdYNeutExpVal = fitExpandGauss(TDCstdYneut,NeutTree,-0.006,0,10);

  // Fit all of the TDC interquartile range distributions
  vector<Double_t> TDCiqrNnbarBifurGaussVals = fitBifurGauss(TDCiqrnnbar,NnbarTree,1000,250,250);
  vector<Double_t> TDCiqrNeutExpVal = fitExpandGauss(TDCiqrneut,NeutTree,-0.006,0,10);
  vector<Double_t> TDCiqrUNnbarBifurGaussVals = fitBifurGauss(TDCiqrUnnbar,NnbarTree,1000,250,250);
  vector<Double_t> TDCiqrUNeutExpVal = fitExpandGauss(TDCiqrUneut,NeutTree,-0.006,20,15);
  vector<Double_t> TDCiqrVNnbarBifurGaussVals = fitBifurGauss(TDCiqrVnnbar,NnbarTree,1000,250,250);
  vector<Double_t> TDCiqrVNeutExpVal = fitExpandGauss(TDCiqrVneut,NeutTree,-0.006,0,10);
  vector<Double_t> TDCiqrYNnbarBifurGaussVals = fitBifurGauss(TDCiqrYnnbar,NnbarTree,1000,250,250);
  vector<Double_t> TDCiqrYNeutExpVal = fitExpandGauss(TDCiqrYneut,NeutTree,-0.006,0,10);

  // Fit all of the mean ADC amplitude distributions
  vector<Double_t> ADCampGaussvals = fitBifurGauss(ADCampnnbar,NnbarTree,25,250,250);
  vector<Double_t> ADCampLandvals = fitLand(ADCampneut,NeutTree,30,500);
  vector<Double_t> ADCampUGaussvals = fitBifurGauss(ADCampUnnbar,NnbarTree,25,250,250);
  vector<Double_t> ADCampULandvals = fitLand(ADCampUneut,NeutTree,30,500);
  vector<Double_t> ADCampVGaussvals = fitLandandGauss(ADCampVnnbar,NnbarTree,17,150,20,150);
  vector<Double_t> ADCampVLandvals = fitLandandGauss(ADCampVneut,NeutTree,20,100,20,100);
  vector<Double_t> ADCampYGaussvals = fitLandandGauss(ADCampYnnbar,NnbarTree,30,100,30,100);
  vector<Double_t> ADCampYLandvals = fitLand(ADCampYneut,NeutTree,30,500);

  // Fit all of the integrated waveform distributions
  vector<Double_t> WFintNeutExpVal = fitExpandGauss(WFintneut,NeutTree,-0.000006,2000,1000);
  vector<Double_t> WFintNnbarGaussVals = fitGauss(WFintnnbar,NnbarTree,550000,100000);
  vector<Double_t> WFintUNeutExpVal = fitExpandGauss(WFintU,NeutTree,-0.000006,1000,1000);
  vector<Double_t> WFintUNnbarGaussVals = fitGauss(WFintU,NnbarTree,200000,100000);
  vector<Double_t> WFintVNeutExpVal = fitExpandGauss(WFintV,NeutTree,-0.000006,500,500);
  vector<Double_t> WFintVNnbarGaussVals = fitGauss(WFintV,NnbarTree,200000,100000);
  vector<Double_t> WFintYNeutExpVal = fitExpandGauss(WFintY,NeutTree,-0.000006,500,500);
  vector<Double_t> WFintYNnbarGaussVals = fitGauss(WFintY,NnbarTree,200000,100000);

  // Get data from files
  vector<vector<float>> NnbarEvents = ReadTree("allnnbar.root",false);
  vector<vector<float>> NeutEvents = ReadTree("allneut.root",true);
  
  TFile *OutFile2 = new TFile("NnbarLogLikelihoods.root","RECREATE");
  TTree* _ll_nnbar_tree;
  Double_t _ll_nnbarhnu, _ll_nnbarhnv, _ll_nnbarhny, _ll_nnbartsu, _ll_nnbartsv, _ll_nnbartsy, _ll_nnbartiu, _ll_nnbartiv, _ll_nnbartiy, _ll_nnbaraau, _ll_nnbaraav, _ll_nnbaraay, _ll_nnbarwfu, _ll_nnbarwfv, _ll_nnbarwfy;

  _ll_nnbar_tree = new TTree("ch_tree","");
  _ll_nnbar_tree->Branch("_ll_nnbarhnu",&_ll_nnbarhnu,"_ll_nnbarhnu/D");
  _ll_nnbar_tree->Branch("_ll_nnbarhnv",&_ll_nnbarhnv,"_ll_nnbarhnv/D");
  _ll_nnbar_tree->Branch("_ll_nnbarhny",&_ll_nnbarhny,"_ll_nnbarhny/D");
  _ll_nnbar_tree->Branch("_ll_nnbartsu",&_ll_nnbartsu,"_ll_nnbartsu/D");
  _ll_nnbar_tree->Branch("_ll_nnbartsv",&_ll_nnbartsv,"_ll_nnbartsv/D");
  _ll_nnbar_tree->Branch("_ll_nnbartsy",&_ll_nnbartsy,"_ll_nnbartsy/D");
  _ll_nnbar_tree->Branch("_ll_nnbartiu",&_ll_nnbartiu,"_ll_nnbartiu/D");
  _ll_nnbar_tree->Branch("_ll_nnbartiv",&_ll_nnbartiv,"_ll_nnbartiv/D");
  _ll_nnbar_tree->Branch("_ll_nnbartiy",&_ll_nnbartiy,"_ll_nnbartiy/D");
  _ll_nnbar_tree->Branch("_ll_nnbaraau",&_ll_nnbaraau,"_ll_nnbaraau/D");
  _ll_nnbar_tree->Branch("_ll_nnbaraav",&_ll_nnbaraav,"_ll_nnbaraav/D");
  _ll_nnbar_tree->Branch("_ll_nnbaraay",&_ll_nnbaraay,"_ll_nnbaraay/D");
  _ll_nnbar_tree->Branch("_ll_nnbarwfu",&_ll_nnbarwfu,"_ll_nnbarwfu/D");
  _ll_nnbar_tree->Branch("_ll_nnbarwfv",&_ll_nnbarwfv,"_ll_nnbarwfv/D");
  _ll_nnbar_tree->Branch("_ll_nnbarwfy",&_ll_nnbarwfy,"_ll_nnbarwfy/D");

  // Create histograms to plot the log likelihoods on
  TH1D *h_NNBAR = new TH1D("h_NNBAR","h_NNBAR",100,-20,20);
  TH1D *h_NEUT = new TH1D("h_NEUT","h_NEUT",100,-20,20);
  
  // Initialise removed counters
  int removedneut = 0;
  int removednnbar = 0; 

  // Loop over all of the events (10000) - MAKE SURE FIRST LL IS = INSTEAD OF +=
  for(int k=0; k<NnbarEvents[0].size(); k++){

    // U plane hit number log likelihoods
    Double_t pnnbarnnbarhnu = BifurGauss(NnbarEvents[1][k],hitNoUNnbarBifurGaussVals[0],hitNoUNnbarBifurGaussVals[1],hitNoUNnbarBifurGaussVals[2]);
    Double_t pnnbarneuthnu = ExpandGauss(NnbarEvents[1][k],hitNoUNeutExpVal[0],hitNoUNeutExpVal[1],hitNoUNeutExpVal[2],hitNoUNeutExpVal[3],hitNoUNeutExpVal[4]);
    _ll_nnbarhnu = log(pnnbarnnbarhnu/pnnbarneuthnu);

    // V plane hit number log likelihoods
    Double_t pnnbarnnbarhnv = BifurGauss(NnbarEvents[2][k],hitNoVNnbarBifurGaussVals[0],hitNoVNnbarBifurGaussVals[1],hitNoVNnbarBifurGaussVals[2]);
    Double_t pnnbarneuthnv = ExpandGauss(NnbarEvents[2][k],hitNoVNeutExpVal[0],hitNoVNeutExpVal[1],hitNoVNeutExpVal[2],hitNoVNeutExpVal[3],hitNoVNeutExpVal[4]);
    _ll_nnbarhnv = log(pnnbarnnbarhnv/pnnbarneuthnv);

    // V plane hit number log likelihoods
    Double_t pnnbarnnbarhny = BifurGauss(NnbarEvents[3][k],hitNoYNnbarBifurGaussVals[0],hitNoYNnbarBifurGaussVals[1],hitNoYNnbarBifurGaussVals[2]);
    Double_t pnnbarneuthny = ExpandGauss(NnbarEvents[3][k],hitNoYNeutExpVal[0],hitNoYNeutExpVal[1],hitNoYNeutExpVal[2],hitNoYNeutExpVal[3],hitNoYNeutExpVal[4]);
    _ll_nnbarhny = log(pnnbarnnbarhny/pnnbarneuthny);

  
    // U plane TDC standard deviation log likelihoods
    Double_t pnnbarnnbartsu = TMath::Gaus(NnbarEvents[5][k], TDCstdUNnbarGaussVals[0], TDCstdUNnbarGaussVals[1]);
    Double_t pnnbarneuttsu = ExpandGauss(NnbarEvents[5][k],TDCstdUNeutExpVal[0],TDCstdUNeutExpVal[1],TDCstdUNeutExpVal[2],TDCstdUNeutExpVal[3],TDCstdUNeutExpVal[4]);
    _ll_nnbartsu = log(pnnbarnnbartsu/pnnbarneuttsu);

    // V plane TDC standard deviation log likelihoods
    Double_t pnnbarnnbartsv = TMath::Gaus(NnbarEvents[6][k], TDCstdVNnbarGaussVals[0], TDCstdVNnbarGaussVals[1]);
    Double_t pnnbarneuttsv = ExpandGauss(NnbarEvents[6][k],TDCstdVNeutExpVal[0],TDCstdVNeutExpVal[1],TDCstdVNeutExpVal[2],TDCstdVNeutExpVal[3],TDCstdVNeutExpVal[4]);
    _ll_nnbartsv = log(pnnbarnnbartsv/pnnbarneuttsv);

    // Y plane TDC standard deviation log likelihoods
    Double_t pnnbarnnbartsy = TMath::Gaus(NnbarEvents[7][k], TDCstdYNnbarGaussVals[0], TDCstdYNnbarGaussVals[1]);
    Double_t pnnbarneuttsy = ExpandGauss(NnbarEvents[7][k],TDCstdYNeutExpVal[0],TDCstdYNeutExpVal[1],TDCstdYNeutExpVal[2],TDCstdYNeutExpVal[3],TDCstdYNeutExpVal[4]);
    _ll_nnbartsy = log(pnnbarnnbartsy/pnnbarneuttsy);


    // U plane TDC interquartile range log likelihoods
    Double_t pnnbarnnbartiu = BifurGauss(NnbarEvents[9][k], TDCiqrUNnbarBifurGaussVals[0], TDCiqrUNnbarBifurGaussVals[1],TDCiqrUNnbarBifurGaussVals[2]);
    Double_t pnnbarneuttiu = ExpandGauss(NnbarEvents[9][k],TDCiqrUNeutExpVal[0],TDCiqrUNeutExpVal[1],TDCiqrUNeutExpVal[2],TDCiqrUNeutExpVal[3],TDCiqrUNeutExpVal[4]);
    _ll_nnbartiu = log(pnnbarnnbartiu/pnnbarneuttiu);

    // V plane TDC interquartile range log likelihoods
    Double_t pnnbarnnbartiv = BifurGauss(NnbarEvents[10][k], TDCiqrVNnbarBifurGaussVals[0], TDCiqrVNnbarBifurGaussVals[1],TDCiqrVNnbarBifurGaussVals[2]);
    Double_t pnnbarneuttiv = ExpandGauss(NnbarEvents[10][k],TDCiqrVNeutExpVal[0],TDCiqrVNeutExpVal[1],TDCiqrVNeutExpVal[2],TDCiqrVNeutExpVal[3],TDCiqrVNeutExpVal[4]);
    _ll_nnbartiv = log(pnnbarnnbartiv/pnnbarneuttiv);

    // Y plane TDC interquartile range log likelihoods
    Double_t pnnbarnnbartiy = BifurGauss(NnbarEvents[11][k], TDCiqrYNnbarBifurGaussVals[0], TDCiqrYNnbarBifurGaussVals[1],TDCiqrYNnbarBifurGaussVals[2]);
    Double_t pnnbarneuttiy = ExpandGauss(NnbarEvents[11][k],TDCiqrYNeutExpVal[0],TDCiqrYNeutExpVal[1],TDCiqrYNeutExpVal[2],TDCiqrYNeutExpVal[3],TDCiqrYNeutExpVal[4]);
    _ll_nnbartiy = log(pnnbarnnbartiy/pnnbarneuttiy);


    // U plane mean ADC amplitude log likelihoods
    Double_t pnnbarnnbaraau = BifurGauss(NnbarEvents[13][k], ADCampUGaussvals[0], ADCampUGaussvals[1], ADCampUGaussvals[2]);
    Double_t pnnbarneutaau = TMath::Landau(NnbarEvents[13][k], ADCampULandvals[0], ADCampULandvals[1]);
    _ll_nnbaraau = log(pnnbarnnbaraau/pnnbarneutaau);

    // V plane mean ADC amplitude log likelihoods
    Double_t pnnbarnnbaraav = LandandGauss(NnbarEvents[14][k], ADCampVGaussvals[0], ADCampVGaussvals[1], ADCampVGaussvals[2], ADCampVGaussvals[3], ADCampVGaussvals[4], ADCampVGaussvals[5]);
    Double_t pnnbarneutaav = LandandGauss(NnbarEvents[14][k], ADCampVLandvals[0], ADCampVLandvals[1], ADCampVLandvals[2], ADCampVLandvals[3], ADCampVLandvals[4], ADCampVLandvals[5]);
    _ll_nnbaraav = log(pnnbarnnbaraav/pnnbarneutaav);

    // Y plane mean ADC amplitude log likelihoods
    Double_t pnnbarnnbaraay = LandandGauss(NnbarEvents[15][k], ADCampYGaussvals[0], ADCampYGaussvals[1], ADCampYGaussvals[2], ADCampYGaussvals[3], ADCampYGaussvals[4], ADCampYGaussvals[5]);
    Double_t pnnbarneutaay = TMath::Landau(NnbarEvents[15][k], ADCampYLandvals[0], ADCampYLandvals[1]);
    _ll_nnbaraay = log(pnnbarnnbaraay/pnnbarneutaay);


    // U plane integrated waveform log likelihoods
    Double_t pnnbarnnbarwfu = TMath::Gaus(NnbarEvents[17][k], WFintUNnbarGaussVals[0], WFintUNnbarGaussVals[1]);
    Double_t pnnbarneutwfu = ExpandGauss(NnbarEvents[17][k],WFintUNeutExpVal[0],WFintUNeutExpVal[1],WFintUNeutExpVal[2],WFintUNeutExpVal[3],WFintUNeutExpVal[4]);
    _ll_nnbarwfu = log(pnnbarnnbarwfu/pnnbarneutwfu);

    // V plane integrated waveform log likelihoods
    Double_t pnnbarnnbarwfv = TMath::Gaus(NnbarEvents[18][k], WFintVNnbarGaussVals[0], WFintVNnbarGaussVals[1]);
    Double_t pnnbarneutwfv = ExpandGauss(NnbarEvents[18][k],WFintVNeutExpVal[0],WFintVNeutExpVal[1],WFintVNeutExpVal[2],WFintVNeutExpVal[3],WFintVNeutExpVal[4]);
    _ll_nnbarwfv = log(pnnbarnnbarwfv/pnnbarneutwfv);

    // Y plane integrated waveform log likelihoods
    Double_t pnnbarnnbarwfy = TMath::Gaus(NnbarEvents[19][k], WFintYNnbarGaussVals[0], WFintYNnbarGaussVals[1]);
    Double_t pnnbarneutwfy = ExpandGauss(NnbarEvents[19][k],WFintYNeutExpVal[0],WFintYNeutExpVal[1],WFintYNeutExpVal[2],WFintYNeutExpVal[3],WFintYNeutExpVal[4]);
    _ll_nnbarwfy = log(pnnbarnnbarwfy/(pnnbarneutwfy));

    // Apply cuts on the log likelihoods - TO BE REPLACED BY LUKES OPTIMISER
    //if(_ll_nnbar<-5||_ll_nnbar>15){removednnbar += 1;}
    // Fill histograms
    //h_NNBAR->Fill(_ll_nnbar);
    _ll_nnbar_tree->Fill();
  }

  _ll_nnbar_tree->Write();
  OutFile2->Close();
  
  TFile *OutFile3 = new TFile("NeutLogLikelihoods.root","RECREATE");
  TTree* _ll_neut_tree;
  Double_t _ll_neuthnu, _ll_neuthnv, _ll_neuthny, _ll_neuttsu, _ll_neuttsv, _ll_neuttsy, _ll_neuttiu, _ll_neuttiv, _ll_neuttiy, _ll_neutaau, _ll_neutaav, _ll_neutaay, _ll_neutwfu, _ll_neutwfv, _ll_neutwfy;

  _ll_neut_tree = new TTree("ch_tree","");
  _ll_neut_tree->Branch("_ll_neuthnu",&_ll_neuthnu,"_ll_neuthnu/D");
  _ll_neut_tree->Branch("_ll_neuthnv",&_ll_neuthnv,"_ll_neuthnv/D");
  _ll_neut_tree->Branch("_ll_neuthny",&_ll_neuthny,"_ll_neuthny/D");
  _ll_neut_tree->Branch("_ll_neuttsu",&_ll_neuttsu,"_ll_neuttsu/D");
  _ll_neut_tree->Branch("_ll_neuttsv",&_ll_neuttsv,"_ll_neuttsv/D");
  _ll_neut_tree->Branch("_ll_neuttsy",&_ll_neuttsy,"_ll_neuttsy/D");
  _ll_neut_tree->Branch("_ll_neuttiu",&_ll_neuttiu,"_ll_neuttiu/D");
  _ll_neut_tree->Branch("_ll_neuttiv",&_ll_neuttiv,"_ll_neuttiv/D");
  _ll_neut_tree->Branch("_ll_neuttiy",&_ll_neuttiy,"_ll_neuttiy/D");
  _ll_neut_tree->Branch("_ll_neutaau",&_ll_neutaau,"_ll_neutaau/D");
  _ll_neut_tree->Branch("_ll_neutaav",&_ll_neutaav,"_ll_neutaav/D");
  _ll_neut_tree->Branch("_ll_neutaay",&_ll_neutaay,"_ll_neutaay/D");
  _ll_neut_tree->Branch("_ll_neutwfu",&_ll_neutwfu,"_ll_neutwfu/D");
  _ll_neut_tree->Branch("_ll_neutwfv",&_ll_neutwfv,"_ll_neutwfv/D");
  _ll_neut_tree->Branch("_ll_neutwfy",&_ll_neutwfy,"_ll_neutwfy/D");

  // Loop over all of the events (10000) - MAKE SURE FIRST LL IS = INSTEAD OF +=
  for(int k=0; k<NeutEvents[0].size(); k++){

    // U plane hit number log likelihoods
    Double_t pneutneuthnu = ExpandGauss(NeutEvents[1][k],hitNoUNeutExpVal[0],hitNoUNeutExpVal[1],hitNoUNeutExpVal[2],hitNoUNeutExpVal[3],hitNoUNeutExpVal[4]);
    Double_t pneutnnbarhnu = BifurGauss(NeutEvents[1][k],hitNoUNnbarBifurGaussVals[0],hitNoUNnbarBifurGaussVals[1],hitNoUNnbarBifurGaussVals[2]);
    _ll_neuthnu = log(pneutnnbarhnu/pneutneuthnu);

    // V plane hit number log likelihoods
    Double_t pneutneuthnv = ExpandGauss(NeutEvents[2][k],hitNoVNeutExpVal[0],hitNoVNeutExpVal[1],hitNoVNeutExpVal[2],hitNoVNeutExpVal[3],hitNoVNeutExpVal[4]);
    Double_t pneutnnbarhnv = BifurGauss(NeutEvents[2][k],hitNoVNnbarBifurGaussVals[0],hitNoVNnbarBifurGaussVals[1],hitNoVNnbarBifurGaussVals[2]);
    _ll_neuthnv = log(pneutnnbarhnv/pneutneuthnv);

    // Y plane hit number log likelihoods
    Double_t pneutneuthny = ExpandGauss(NeutEvents[3][k],hitNoYNeutExpVal[0],hitNoYNeutExpVal[1],hitNoYNeutExpVal[2],hitNoYNeutExpVal[3],hitNoYNeutExpVal[4]);
    Double_t pneutnnbarhny = BifurGauss(NeutEvents[3][k],hitNoYNnbarBifurGaussVals[0],hitNoYNnbarBifurGaussVals[1],hitNoYNnbarBifurGaussVals[2]);
    _ll_neuthny = log(pneutnnbarhny/pneutneuthny);

  
    // U plane TDC standard deviation log likelihoods
    Double_t pneutneuttsu = ExpandGauss(NeutEvents[5][k],TDCstdUNeutExpVal[0],TDCstdUNeutExpVal[1],TDCstdUNeutExpVal[2],TDCstdUNeutExpVal[3],TDCstdUNeutExpVal[4]);
    Double_t pneutnnbartsu = TMath::Gaus(NeutEvents[5][k], TDCstdUNnbarGaussVals[0], TDCstdUNnbarGaussVals[1]);
    _ll_neuttsu = log(pneutnnbartsu/pneutneuttsu);

    // V plane TDC standard deviation log likelihoods
    Double_t pneutneuttsv = ExpandGauss(NeutEvents[6][k],TDCstdVNeutExpVal[0],TDCstdVNeutExpVal[1],TDCstdVNeutExpVal[2],TDCstdVNeutExpVal[3],TDCstdVNeutExpVal[4]);
    Double_t pneutnnbartsv = TMath::Gaus(NeutEvents[6][k], TDCstdVNnbarGaussVals[0], TDCstdVNnbarGaussVals[1]);
    _ll_neuttsv = log(pneutnnbartsv/pneutneuttsv);

    // Y plane TDC standard deviation log likelihoods
    Double_t pneutneuttsy = ExpandGauss(NeutEvents[7][k],TDCstdYNeutExpVal[0],TDCstdYNeutExpVal[1],TDCstdYNeutExpVal[2],TDCstdYNeutExpVal[3],TDCstdYNeutExpVal[4]);
    Double_t pneutnnbartsy = TMath::Gaus(NeutEvents[7][k], TDCstdYNnbarGaussVals[0], TDCstdYNnbarGaussVals[1]);
    _ll_neuttsy = log(pneutnnbartsy/pneutneuttsy);


    // U plane TDC interquartile range log likelihoods
    Double_t pneutneuttiu = ExpandGauss(NeutEvents[9][k],TDCiqrUNeutExpVal[0],TDCiqrUNeutExpVal[1],TDCiqrUNeutExpVal[2],TDCiqrUNeutExpVal[3],TDCiqrUNeutExpVal[4]);
    Double_t pneutnnbartiu = BifurGauss(NeutEvents[9][k], TDCiqrUNnbarBifurGaussVals[0], TDCiqrUNnbarBifurGaussVals[1],TDCiqrUNnbarBifurGaussVals[2]);
    _ll_neuttiu = log(pneutnnbartiu/pneutneuttiu);

    // V plane TDC interquartile range log likelihoods
    Double_t pneutneuttiv = ExpandGauss(NeutEvents[10][k],TDCiqrVNeutExpVal[0],TDCiqrVNeutExpVal[1],TDCiqrVNeutExpVal[2],TDCiqrVNeutExpVal[3],TDCiqrVNeutExpVal[4]);
    Double_t pneutnnbartiv = BifurGauss(NeutEvents[10][k], TDCiqrVNnbarBifurGaussVals[0], TDCiqrVNnbarBifurGaussVals[1],TDCiqrVNnbarBifurGaussVals[2]);
    _ll_neuttiv = log(pneutnnbartiv/pneutneuttiv);

    // Y plane TDC interquartile range log likelihoods
    Double_t pneutneuttiy = ExpandGauss(NeutEvents[11][k],TDCiqrYNeutExpVal[0],TDCiqrYNeutExpVal[1],TDCiqrYNeutExpVal[2],TDCiqrYNeutExpVal[3],TDCiqrYNeutExpVal[4]);
    Double_t pneutnnbartiy = BifurGauss(NeutEvents[11][k], TDCiqrYNnbarBifurGaussVals[0], TDCiqrYNnbarBifurGaussVals[1],TDCiqrYNnbarBifurGaussVals[2]);
    _ll_neuttiy = log(pneutnnbartiy/pneutneuttiy);


    // U plane mean ADC amplitude log likelihoods
    Double_t pneutnnbaraau = BifurGauss(NeutEvents[13][k], ADCampUGaussvals[0], ADCampUGaussvals[1], ADCampUGaussvals[2]);
    Double_t pneutneutaau = TMath::Landau(NeutEvents[13][k], ADCampULandvals[0], ADCampULandvals[1]);
    if(pneutnnbaraau==0){pneutnnbaraau = 4.06461e-305;}
    _ll_neutaau = log(pneutnnbaraau/pneutneuttiy);

    // V plane mean ADC amplitude log likelihoods
    Double_t pneutneutaav = LandandGauss(NeutEvents[14][k], ADCampVLandvals[0], ADCampVLandvals[1], ADCampVLandvals[2], ADCampVLandvals[3], ADCampVLandvals[4], ADCampVLandvals[5]);
    Double_t pneutnnbaraav = LandandGauss(NeutEvents[14][k], ADCampVGaussvals[0], ADCampVGaussvals[1], ADCampVGaussvals[2], ADCampVGaussvals[3], ADCampVGaussvals[4], ADCampVGaussvals[5]);
    _ll_neutaav = log(pneutnnbaraav/pneutneutaav);

    // Y plane mean ADC amplitude log likelihoods
    Double_t pneutneutaay = TMath::Landau(NeutEvents[15][k], ADCampYLandvals[0], ADCampYLandvals[1]);
    Double_t pneutnnbaraay = LandandGauss(NeutEvents[15][k], ADCampYGaussvals[0], ADCampYGaussvals[1], ADCampYGaussvals[2], ADCampYGaussvals[3], ADCampYGaussvals[4], ADCampYGaussvals[5]);
    _ll_neutaay = log(pneutnnbaraay/pneutneutaay);


    // U plane integrated waveform log likelihoods
    Double_t pneutneutwfu = ExpandGauss(NeutEvents[17][k],WFintUNeutExpVal[0],WFintUNeutExpVal[1],WFintUNeutExpVal[2],WFintUNeutExpVal[3],WFintUNeutExpVal[4]);
    Double_t pneutnnbarwfu = TMath::Gaus(NeutEvents[17][k], WFintUNnbarGaussVals[0], WFintUNnbarGaussVals[1]);
    _ll_neutwfu = log(pneutnnbarwfu/pneutneutwfu);

    // V plane integrated waveform log likelihoods
    Double_t pneutneutwfv = ExpandGauss(NeutEvents[18][k],WFintVNeutExpVal[0],WFintVNeutExpVal[1],WFintVNeutExpVal[2],WFintVNeutExpVal[3],WFintVNeutExpVal[4]);
    Double_t pneutnnbarwfv = TMath::Gaus(NeutEvents[18][k], WFintVNnbarGaussVals[0], WFintVNnbarGaussVals[1]);
    _ll_neutwfv = log(pneutnnbarwfv/pneutneutwfv);

    // Y plane integrated waveform log likelihoods
    Double_t pneutneutwfy = ExpandGauss(NeutEvents[19][k],WFintYNeutExpVal[0],WFintYNeutExpVal[1],WFintYNeutExpVal[2],WFintYNeutExpVal[3],WFintYNeutExpVal[4]);
    Double_t pneutnnbarwfy = TMath::Gaus(NeutEvents[19][k], WFintYNnbarGaussVals[0], WFintYNnbarGaussVals[1]);
    _ll_neutwfy = log(pneutnnbarwfy/pneutneutwfy);

    // Apply cuts on the log likelihoods - TO BE REPLACED BY LUKES OPTIMISER
    //if(_ll_neut<-5||_ll_neut>15){removedneut += 1;}
    // Fill histograms
    //h_NEUT->Fill(_ll_neut);
    _ll_neut_tree->Fill();
  }

  _ll_neut_tree->Write();
  OutFile3->Close();

/*
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
  leg->Draw();*/

} // End of LogLikelihoods
