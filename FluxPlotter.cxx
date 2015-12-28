//example of macro illustrating how to superimpose two histograms

// Root Includes
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TH2F.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCollection.h"
#include "THStack.h"
#include "TAxis.h"

// C++ Includes
#include <string>
#include <iostream>
#include <fstream>

void FluxPlotter(string FileName)
{
   //example of macro illustrating how to superimpose two histograms
   //with different scales in the "same" pad.
   // To see the output of this macro, click begin_html <a href="gif/twoscales.gif" >here</a> end_html
   //Author: Rene Br

  gROOT->ProcessLine(".x lhcbStyle.C");


  string line;
  ifstream myfile (FileName.c_str());
  if (myfile.is_open()){
//    TFile *f1 = new TFile("Outputfileimgoingtobeusing.root");
    double Energy(0), OldCosTheta(0), CosTheta(0), Flux(0), Throw1(0), Throw2(0), TotalFlux(0);
    int i(0);
    // TH2F* h_EnergyCos = new TH2F("","",20,-1,1,20,0,20);
    TH2F* h_FluxEnergy = new TH2F("lolz","lolz",50,0,200,50,0,30); 
    TH1F* h_FluxCos = new TH1F("wtf","butsrsly",20,-1,1);
    TH1F* h_FluxEnergy1D = new TH1F("","Aname",20,0,10);

    while ( myfile >> Energy >> CosTheta >> Flux >> Throw1 >> Throw2 ){
      //std::cout << "Reading line " << i << "\n";
      i ++;
      if(i == 1) OldCosTheta = CosTheta;
      // h_EnergyCos->Fill(CosTheta,Energy);
      h_FluxEnergy->Fill(Flux,Energy);
      for(unsigned int j(0); j < Flux/40; j ++ ) h_FluxEnergy1D->Fill(Energy);
      if(OldCosTheta != CosTheta){
        for(unsigned int j(0); j < TotalFlux; j++) h_FluxCos->Fill(OldCosTheta);
        TotalFlux = Flux; OldCosTheta = CosTheta;
      }else{ TotalFlux += Flux; }
    }
    for(unsigned int j(0); j < TotalFlux; j++) h_FluxCos->Fill(OldCosTheta);

    myfile.close();
  }  else std::cout << "Unable to open file"; 
  /*
  TCanvas *c1 = new TCanvas("c1","Energy vs Cos",900,600);
  h_EnergyCos->Draw("box");
  h_EnergyCos->GetXaxis()->SetTitle("Cos(Theta) [arb]");
  h_EnergyCos->GetYaxis()->SetTitle("Energy [GeV]");
  */
  TCanvas *c2 = new TCanvas("c2","Flux vs Energy",900,600);
  h_FluxEnergy->Draw("box");
  h_FluxEnergy->GetXaxis()->SetTitle("Flux [Arb]");
  h_FluxEnergy->GetYaxis()->SetTitle("Energy [GeV]");
  
  TCanvas *c3 = new TCanvas("c3","Flux vs Cos",900,600);
  h_FluxCos->Draw();
  h_FluxCos->GetXaxis()->SetTitle("Cos(Theta) [arb]");
  h_FluxCos->GetYaxis()->SetTitle("Flux [arb]");

  TCanvas *c4 = new TCanvas("c4","Energy vs Flux",900,600);
  h_FluxEnergy1D->Draw();
  h_FluxEnergy1D->GetXaxis()->SetTitle("Energy [GeV]");
  h_FluxEnergy1D->GetYaxis()->SetTitle("Flux [arb]");

}
