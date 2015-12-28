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
#include <vector>

void MultipleFluxPlotter(string FileName,int Colour,TLegend* leg1,TLegend* leg2, TLegend* leg3)
{
   //example of macro illustrating how to superimpose two histograms
   //with different scales in the "same" pad.
   // To see the output of this macro, click begin_html <a href="gif/twoscales.gif" >here</a> end_html
   //Author: Rene Br

  gROOT->ProcessLine(".x lhcbStyle.C");
  
  ifstream myfile (FileName.c_str());
  if (myfile.is_open()){
//    TFile *f1 = new TFile("Outputfileimgoingtobeusing.root");
    double Energy(0), OldCosTheta(0), CosTheta(0), Flux(0), Throw1(0), Throw2(0), TotalFlux(0);
    int i(0), p(0), FileLines(0);
    // TH2F* h_EnergyCos = new TH2F("","",20,-1,1,20,0,20);
    TH2F* h_FluxEnergy = new TH2F(("FluxEnergy"+FileName).c_str(),"FluxEnergy",50,0,400,50,0,30); 
    TH1F* h_FluxCos = new TH1F(("FluxCos"+FileName).c_str(),"FluxCos",20,-1,1);
    TH1F* h_FluxEnergy1D = new TH1F(("FluxEnergy1D"+FileName).c_str(),"FluxEnergy1D",20,0,10);
    string Throw; std::getline(myfile,Throw); // THIS THROWS AWAY THE TOP LINE. CHECK FLUX FILE FORMAT
    while ( myfile >> Energy >> CosTheta >> Flux >> Throw1 >> Throw2 ){
      //std::cout << "\nIt gets here!";
      FileLines ++;
    }
    std::cout << "\nFileLines = " << FileLines << "\n";
    myfile.clear(); myfile.seekg(0, ios::beg);
    string Throw; std::getline(myfile,Throw); // THIS THROWS AWAY THE TOP LINE. CHECK FLUX FILE FORMAT
    while ( myfile >> Energy >> CosTheta >> Flux >> Throw1 >> Throw2 ){
      // std::cout << "Reading line " << i << "\n";
      i ++;
      if(i == 1) OldCosTheta = CosTheta;
      // Flux /= Energy; MAY NEED THIS LINE TO COUNTER LOG(E)
      // h_EnergyCos->Fill(CosTheta,Energy);
      h_FluxEnergy->Fill(Flux,Energy);
      for(unsigned int j(0); j < Flux; j ++ ) h_FluxEnergy1D->Fill(Energy);
      if(OldCosTheta != CosTheta){
        //std::cout << "\nTotalFlux = " << TotalFlux;
	//std::cout << "\nOldCosTheta = " << OldCosTheta;
        for(unsigned int k(0); k < TotalFlux / p; k++) h_FluxCos->Fill(OldCosTheta);
        TotalFlux = Flux*10000; OldCosTheta = CosTheta; p = 0;
      }else{ TotalFlux += Flux*10000; p ++; }
    }
    for(unsigned int j(0); j < TotalFlux / p; j++) h_FluxCos->Fill(OldCosTheta);
    myfile.close();
  }  else{ std::cout << "\nUnable to open file: " << FileName << "\n"; exit(1); }
  
  
  /*
  TCanvas *c1 = new TCanvas("c1","Energy vs Cos",900,600);
  h_EnergyCos->Draw("box");
  h_EnergyCos->GetXaxis()->SetTitle("Cos(Theta) [arb]");
  h_EnergyCos->GetYaxis()->SetTitle("Energy [GeV]");
  */

  c2->cd();
  h_FluxEnergy->Draw("samebox");
  h_FluxEnergy->SetLineColor(Colour);
  h_FluxEnergy->GetXaxis()->SetTitle("Flux [Arb]");
  h_FluxEnergy->GetYaxis()->SetTitle("Energy [GeV]");
  leg1->AddEntry(h_FluxEnergy,FileName.c_str(),"l");
  leg1->Draw();  

  c3->cd();
  h_FluxCos->Draw("same");
  h_FluxCos->SetLineColor(Colour);
  h_FluxCos->GetXaxis()->SetTitle("Cos(Theta) [arb]");
  h_FluxCos->GetYaxis()->SetTitle("1000*Flux [arb]");
  leg2->AddEntry(h_FluxCos,FileName.c_str(),"l");
  leg2->Draw();
  

  c4->cd();
  h_FluxEnergy1D->Draw("same");
  h_FluxEnergy1D->SetLineColor(Colour);
  h_FluxEnergy1D->GetXaxis()->SetTitle("Energy [GeV]");
  h_FluxEnergy1D->GetYaxis()->SetTitle("Flux [arb]");
  leg3->AddEntry(h_FluxEnergy1D,FileName.c_str(),"l");
  leg3->Draw("same");
}
