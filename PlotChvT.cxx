#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TGraph.h"
#include <utility>

// Plots hit amplitude against time with shower times marked
void PlotChvT() {

  gROOT->ProcessLine(".x lhcbStyle.C");

  // Create canvas for drawing
  TCanvas *c2 = new TCanvas("c2","Channel Vs Time - U",200,10,1050,650);
  TCanvas *c3 = new TCanvas("c3","Channel Vs Time - V",200,10,1050,650);
  TCanvas *c4 = new TCanvas("c4","Channel Vs Time - Y",200,10,1050,650);

  // Open file and get relevant info from TTree
  TFile *f = TFile::Open("bnb/ShowerAna_output_bnb_pandora.root","READ");

  TTree *_t_ch; f->GetObject("ch_tree",_t_ch);

  std::vector<float> *UChvec = 0;
  std::vector<float> *VChvec = 0;
  std::vector<float> *YChvec = 0;

  std::vector<float> *UTDCvec = 0;
  std::vector<float> *VTDCvec = 0;
  std::vector<float> *YTDCvec = 0;

  std::vector<float> *Evec = 0;

  UInt_t *NoHits = 0;

  std::vector<std::pair<float,float>> *ShowerStartEnd = 0;

  TBranch *buch = 0;
  _t_ch->SetBranchAddress("UChvec",&UChvec,&buch);
  TBranch *bvch = 0;
  _t_ch->SetBranchAddress("VChvec",&VChvec,&bvch);
  TBranch *bych = 0;
  _t_ch->SetBranchAddress("YChvec",&YChvec,&bych);
  TBranch *butdc = 0;
  _t_ch->SetBranchAddress("UTDCvec",&UTDCvec,&butdc);
  TBranch *bvtdc = 0;
  _t_ch->SetBranchAddress("VTDCvec",&VTDCvec,&bvtdc);
  TBranch *bytdc = 0;
  _t_ch->SetBranchAddress("YTDCvec",&YTDCvec,&bytdc);
  TBranch *bsse = 0;
  _t_ch->SetBranchAddress("ShowerStartEnd",&ShowerStartEnd,&bsse);
  TBranch *ben = 0;
  _t_ch->SetBranchAddress("Evec",&Evec,&ben);
  TBranch *bnohits = 0;
  _t_ch->SetBranchAddress("hitNo",&NoHits,&bnohits);

  Long64_t EntryNumber = _t_ch->GetEntries();

  // Loop over number of events
  for(unsigned int j=0; j!=EntryNumber; ++j){
    buch->GetEntry(j);
    bvch->GetEntry(j);
    bych->GetEntry(j);
    butdc->GetEntry(j);
    bvtdc->GetEntry(j);
    bytdc->GetEntry(j);
    bsse->GetEntry(j);
    ben->GetEntry(j);
    bnohits->GetEntry(j);

    int usize = UChvec->size();
    std::cout<<"Number of U hits: "<<usize<<std::endl;

    c2->cd();
    // Plot hit time against hit amplitude
    TGraph *g2 = new TGraph(usize,&((*UTDCvec)[0]),&((*UChvec)[0]));
    g2->SetMarkerSize(0.4);
    g2->SetMarkerColor(kBlue);
    g2->GetXaxis()->SetTitle("TDC Time");
    g2->GetYaxis()->SetTitle("Channel Number");
    g2->Draw("AP");
    c2->Modified();
    c2->Update();

    int vsize = VChvec->size();
    std::cout<<"Number of V hits: "<<vsize<<std::endl;

    c3->cd();
    // Plot hit time against hit amplitude
    TGraph *g3 = new TGraph(vsize,&((*VTDCvec)[0]),&((*VChvec)[0]));
    g3->SetMarkerSize(0.4);
    g3->SetMarkerColor(kBlue);
    g3->GetXaxis()->SetTitle("TDC Time");
    g3->GetYaxis()->SetTitle("Channel Number");
    g3->Draw("AP");
    c3->Modified();
    c3->Update();

    int ysize = YChvec->size();
    std::cout<<"Number of Y hits: "<<ysize<<std::endl;

    c4->cd();
    // Plot hit time against hit amplitude
    TGraph *g4 = new TGraph(ysize,&((*YTDCvec)[0]),&((*YChvec)[0]));
    g4->SetMarkerSize(0.4);
    g4->SetMarkerColor(kBlue);
    g4->GetXaxis()->SetTitle("TDC Time");
    g4->GetYaxis()->SetTitle("Channel Number");
    g4->Draw("AP");
    c4->Modified();
    c4->Update();

    // Loop over number of showers
    for(int i = 0; i < ShowerStartEnd->size(); i++){
      std::cout<<"Shower!  Energy= "<<Evec->at(i)<<" MeV"<<std::endl;
      std::cout<<"Start: "<<(ShowerStartEnd->at(i)).first<<", End: "<<(ShowerStartEnd->at(i)).second<<std::endl;
    }
   
    std::cout<<"Graph "<<j<<std::endl;
    std::cin.get();
  }
}
