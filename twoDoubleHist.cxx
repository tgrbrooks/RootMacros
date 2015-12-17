//example of macro illustrating how to superimpose two histograms
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TFile.h"
#include "TLegend.h"

void twoDoubleHist()
{
   //example of macro illustrating how to superimpose two histograms
   //with different scales in the "same" pad.
   // To see the output of this macro, click begin_html <a href="gif/twoscales.gif" >here</a> end_html
   //Author: Rene Brun
gROOT->ProcessLine(".x lhcbStyle.C");
   TFile *f1 = new TFile("WFint-nnbar10-old.root");
   TH2I *h1 = (TH2I*)f1->Get("h_VYHITS");
   TCanvas *c1 = new TCanvas("c1","Hit Wires",900,600);
   TFile *f2 = new TFile("WFint-numi10-old.root");
   TH2I *h2 = (TH2I*)f2->Get("h_VYHITS");
   h2->SetLineColor(kRed);
   h2->GetXaxis()->SetTitle("V Integrated WF");
   h2->GetYaxis()->SetTitle("Y Integrated WF");
   h2->GetXaxis()->SetTitleOffset(0.75);
   h2->GetYaxis()->SetTitleOffset(0.9);
   h2->Draw("box");
   h1->Draw("boxsame");
   leg = new TLegend(0.3,0.3,0.5,0.5);
   leg->AddEntry(h1,"nnbar","l");
   leg->AddEntry(h2,"numi","l");
   leg->Draw();
}
