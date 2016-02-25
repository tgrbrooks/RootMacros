//example of macro illustrating how to superimpose two histograms
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TFile.h"
#include "TLegend.h"
#include "TChain.h"

// Function to plot a stacked histogram of all neutrino energies split by their flavours
void TTreeEnergyAndOneD(TChain &RootTree,string Variable){

  // Create new canvas
  TCanvas *c1 = new TCanvas(Variable.c_str(),"TTreeOneDHist Canvas",900,600);
  c1->Divide(1,2);
  c1->cd(1);
  RootTree.Draw(Variable.c_str(),"Evec<500","",10000000000,0); // Looking for E < 500MeV showers  

  // Draw histograms and adjust axis
  htemp->GetXaxis()->SetTitle(Variable.c_str());  htemp->GetXaxis()->SetTitleOffset(0.9);
  htemp->GetYaxis()->SetTitle("Frequency");       htemp->GetYaxis()->SetTitleOffset(0.9);
  htemp->Draw();

//  spad1 = new TPad(("spad1 " + Variable).c_str(),"The first subpad",.1,.1,.5,.5);
//  spad1->cd();
  c1->cd(2);
  RootTree.Draw(("Evec:"+Variable).c_str(),"Evec<500","colz", 5000, 0);

  // Draw stacked histograms and adjust axis
//  htemp2->GetXaxis()->SetTitle(Variable.c_str());  htemp2->GetXaxis()->SetTitleOffset(0.9);
//  htemp1->GetYaxis()->SetTitle("Evec (MeV)");      htemp1->GetYaxis()->SetTitleOffset(0.9);
//  htemp1->Draw("colz");
//  spad1->Draw("colz");

  c1->Update(); c1->Modified();
}
