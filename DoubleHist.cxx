//example of macro illustrating how to superimpose two histograms
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TFile.h"
#include "TLegend.h"

void DoubleHist()
{
   //example of macro illustrating how to superimpose two histograms
   //with different scales in the "same" pad.
   // To see the output of this macro, click begin_html <a href="gif/twoscales.gif" >here</a> end_html
   //Author: Rene Brun

gROOT->ProcessLine(".x lhcbStyle.C");
   TFile *f1 = new TFile("TDCiqr-nnbar10.root");
   TH1I *h1 = (TH1I*)f1->Get("h_UHITS");
   TCanvas *c1 = new TCanvas("c1","Hit Wires",900,600);
   TFile *f2 = new TFile("TDCiqr-numi10.root");
   TH1I *h2 = (TH1I*)f2->Get("h_UHITS");
   h2->SetLineColor(kRed);
   h1->GetXaxis()->SetTitle("Hit Number");
   h1->GetYaxis()->SetTitle("Events /(50)");
   h1->GetXaxis()->SetTitleOffset(0.8);
   h1->GetYaxis()->SetTitleOffset(0.6);
   h1->GetXaxis()->SetRange(0,1000);
   Double_t norm = 0.2;
   h2->Scale(norm);
   h1->Draw("E1HIST");
   h2->Draw("E1HISTsame");
   leg = new TLegend(0.3,0.3,0.5,0.5);
   leg->AddEntry(h1,"nnbar","l");
   leg->AddEntry(h2,"numi","l");
   leg->Draw();
}
