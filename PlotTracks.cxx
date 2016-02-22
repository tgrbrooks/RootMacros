#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TGraph.h"
#include <utility>

// Plots hit amplitude against time with shower times marked
void PlotTracks() {

  gROOT->ProcessLine(".x lhcbStyle.C");

  // Create canvas for drawing
  TCanvas *c4 = new TCanvas("c4","Amplitude Vs Channel - Y Tracks",200,10,1050,650);

  // Open file and get relevant info from TTree
  TFile *f = TFile::Open("FindShower_output.root","READ");

  TTree *_t_ch; f->GetObject("ch_tree",_t_ch);

  std::vector<std::vector<std::pair<float,float>>> *Tracks = 0;
  std::vector<std::pair<float,float>> Track;
  std::vector<float> TDC, Chan;

  TBranch *btracks = 0;
  _t_ch->SetBranchAddress("Tracks",&Tracks,&btracks);

  Long64_t EntryNumber = _t_ch->GetEntries();

  TGraph *gr[100];

  // Loop over number of events
  for(unsigned int j=0; j!=EntryNumber; ++j){
    btracks->GetEntry(j);

    int size = Tracks->size();
    std::cout<<"Number of tracks: "<<size<<std::endl;

    TMultiGraph *mg = new TMultiGraph();

    c4->cd();
    for(int k = 0; k < size; ++k){
      Track = (*Tracks)[k];
      int length = Track.size();
      std::cout<<"Track "<<k<<" Length "<<length<<std::endl;
      for(int m = 0; m < length; ++m){
        TDC.push_back(Track[m].first);
        Chan.push_back(Track[m].second);
      }
      std::cout<<"TDC length "<<TDC.size()<<" Chan length "<<Chan.size()<<std::endl;
      // Plot hit time against hit amplitude
      std::cout<<"s";
      gr[k] = new TGraph(length,&(TDC[0]),&(Chan[0]));
      std::cout<<"e";
      gr[k]->SetMarkerSize(0.4);
      gr[k]->SetMarkerColor(k+1);
      std::cout<<"g";
      mg->Add(gr[k]);
      TDC.clear();
      Chan.clear();
    }

    //mg->GetXaxis()->SetTitle("TDC Time");
    //mg->GetYaxis()->SetTitle("Channel Number");
    mg->Draw("AP");
    c4->Modified();
    c4->Update();
   
    for (int n = 0; n < size; ++n){
      gr[n] = 0;
    }    

    std::cout<<"Graph "<<j<<std::endl;
    std::cin.get();
    Track.clear();
  }
}
