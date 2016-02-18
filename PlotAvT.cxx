#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TGraph.h"
#include <utility>

// Plots hit amplitude against time with shower times marked
void PlotAvT() {

  gROOT->ProcessLine(".x lhcbStyle.C");

  // Open file and get relevant info from TTree
  TFile *f = TFile::Open("ShowerAna_output_pandora.root","READ");

  TTree *_t_ch; f->GetObject("ch_tree",_t_ch);

  std::vector<float> *ADCvec = 0;

  std::vector<float> *TDCvec = 0;

  std::vector<float> *Evec = 0;

  UInt_t *NoHits = 0;

  std::vector<std::pair<float,float>> *ShowerStartEnd = 0;

  std::vector<float> ShowerADC, ShowerTDC;

  // Create canvas for drawing
  TCanvas *c1 = new TCanvas("c1","Amplitude Vs Time",200,10,1050,650);

  TBranch *badc = 0;
  _t_ch->SetBranchAddress("ADCvec",&ADCvec,&badc);
  TBranch *btdc = 0;
  _t_ch->SetBranchAddress("TDCvec",&TDCvec,&btdc);
  TBranch *bsse = 0;
  _t_ch->SetBranchAddress("ShowerStartEnd",&ShowerStartEnd,&bsse);
  TBranch *ben = 0;
  _t_ch->SetBranchAddress("Evec",&Evec,&ben);
  TBranch *bnohits = 0;
  _t_ch->SetBranchAddress("hitNo",&NoHits,&bnohits);

  Long64_t EntryNumber = _t_ch->GetEntries();

  // Loop over number of events
  for(unsigned int j=0; j!=EntryNumber; ++j){
    badc->GetEntry(j);
    btdc->GetEntry(j);
    bsse->GetEntry(j);
    ben->GetEntry(j);
    bnohits->GetEntry(j);
    int size = ADCvec->size();
    int no = NoHits;
    
    std::cout<<"Number of hits: "<<size<<std::endl;
    std::cout<<"Number of TDC: " <<TDCvec->size()<<"\n";

    // Plot hit time against hit amplitude
    TGraph *g = new TGraph(size,&((*TDCvec)[0]),&((*ADCvec)[0]));
    g->SetMarkerSize(0.7);
    g->GetXaxis()->SetTitle("TDC Time");
    g->GetYaxis()->SetTitle("ADC Amplitude");
    g->Draw("AP");
    c1->Modified();
    c1->Update();

    bool showerFlag = false;

    // Loop over number of showers
    for(int i = 0; i < ShowerStartEnd->size(); i++){
      std::cout<<"Shower!  Energy= "<<Evec->at(i)<<" MeV"<<std::endl;
      std::cout<<"Start: "<<(ShowerStartEnd->at(i)).first<<", End: "<<(ShowerStartEnd->at(i)).second<<std::endl;
      // For each shower loop over the hits
      for(int k = 0; k < ADCvec->size(); k++){
        // Check if hit occured within shower time window
        if(TDCvec->at(k) > (ShowerStartEnd->at(i)).first && TDCvec->at(k) < (ShowerStartEnd->at(i)).second ){
          ShowerADC.push_back(ADCvec->at(k));
          ShowerTDC.push_back(TDCvec->at(k));
          showerFlag = true;
        }
      }
    }
    
    // Mark the hits that occured within the shower time window
    if(showerFlag){
      TGraph *gshow = new TGraph(ShowerADC.size(),&(ShowerTDC[0]),&(ShowerADC[0]));
      gshow->SetMarkerSize(1.7);
      gshow->SetMarkerColor(kBlue);
      gshow->Draw("P");
      c1->Modified();
      c1->Update();
    }
    
    std::cout<<"Graph "<<j<<std::endl;
    std::cin.get();
    ShowerADC.clear(); ShowerTDC.clear();
  }
}
