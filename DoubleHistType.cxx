//example of macro illustrating how to superimpose two histograms
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCollection.h"
#include <iostream>
#include "THStack.h"
#include "TAxis.h"

void DoubleHistType()
{
   //example of macro illustrating how to superimpose two histograms
   //with different scales in the "same" pad.
   // To see the output of this macro, click begin_html <a href="gif/twoscales.gif" >here</a> end_html
   //Author: Rene Br

gROOT->ProcessLine(".x lhcbStyle.C");

   TCanvas *c1 = new TCanvas("c1","Hit Number",900,600);
   TFile *f1 = new TFile("TDCiqr-numi10-old.root");
   THStack *hs = new THStack("hs","Stacked 1D histograms");

   TH1D *hCCQE = (TH1D*)f1->Get("h_NCUTQCCQE");
   hCCQE->SetFillColor(kRed-4);

   TH1D *hNCQE = (TH1D*)f1->Get("h_NCUTQNCQE");
   hNCQE->SetFillColor(kYellow-4);

   TH1D *hCCRE = (TH1D*)f1->Get("h_NCUTQCCRE");
   hCCRE->SetFillColor(kGreen-4);

   TH1D *hNCRE = (TH1D*)f1->Get("h_NCUTQNCRE");
   hNCRE->SetFillColor(kCyan-4);

   TH1D *hCCDIS = (TH1D*)f1->Get("h_NCUTQCCDIS");
   hCCDIS->SetFillColor(kBlue-4);

   TH1D *hNCDIS = (TH1D*)f1->Get("h_NCUTQNCDIS");
   hNCDIS->SetFillColor(kMagenta-4);

   TH1D *hCCCO = (TH1D*)f1->Get("h_NCUTQCCCO");
   hCCCO->SetFillColor(kMagenta+2);

   TH1D *hNCCO = (TH1D*)f1->Get("h_NCUTQNCCO");
   hNCCO->SetFillColor(kMagenta+4);

   hs->Add(hCCQE);
   hs->Add(hNCQE);
   hs->Add(hCCRE);
   hs->Add(hNCRE);
   hs->Add(hCCDIS);
   hs->Add(hNCDIS);
   hs->Add(hCCCO);
   hs->Add(hNCCO);
   hs->Draw();
   hs->GetXaxis()->SetTitle("E_nu (GeV)");
   hs->GetYaxis()->SetTitle("Events (/0.2 GeV)");
   hs->GetXaxis()->SetTitleOffset(0.9);
   hs->GetYaxis()->SetTitleOffset(0.7);
   hs->Draw();
   leg = new TLegend(0.2,0.8,0.4,0.4);
   leg->AddEntry(hCCQE,"CCQE","f");
   leg->AddEntry(hNCQE,"NCQE","f");
   leg->AddEntry(hCCRE,"CCRE","f");
   leg->AddEntry(hNCRE,"NCRE","f");
   leg->AddEntry(hCCDIS,"CCDIS","f");
   leg->AddEntry(hNCDIS,"NCDIS","f");
   leg->AddEntry(hCCCO,"CCCO","f");
   leg->AddEntry(hNCCO,"NCCO","f");
   leg->Draw();
 
}
