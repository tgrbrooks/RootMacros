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

  // Create canvas for drawing
  TCanvas *c1 = new TCanvas("c1","Amplitude Vs Time - Total",200,10,1050,650);
  TCanvas *c2 = new TCanvas("c2","Amplitude Vs Time - U",200,10,1050,650);
  TCanvas *c3 = new TCanvas("c3","Amplitude Vs Time - V",200,10,1050,650);
  TCanvas *c4 = new TCanvas("c4","Amplitude Vs Time - Y",200,10,1050,650);

  // Open file and get relevant info from TTree
  TFile *f = TFile::Open("bnb/ShowerAna_output_bnb_pandora.root","READ");

  TTree *_t_ch; f->GetObject("ch_tree",_t_ch);

  std::vector<float> *ADCvec = 0;
  std::vector<float> *UADCvec = 0;
  std::vector<float> *VADCvec = 0;
  std::vector<float> *YADCvec = 0;

  std::vector<float> *TDCvec = 0;
  std::vector<float> *UTDCvec = 0;
  std::vector<float> *VTDCvec = 0;
  std::vector<float> *YTDCvec = 0;

  std::vector<float> *Evec = 0;

  UInt_t *NoHits = 0;

  std::vector<std::pair<float,float>> *ShowerStartEnd = 0;

  std::vector<float> ShowerADC, ShowerTDC;

  TBranch *badc = 0;
  _t_ch->SetBranchAddress("ADCvec",&ADCvec,&badc);
  TBranch *buadc = 0;
  _t_ch->SetBranchAddress("UADCvec",&UADCvec,&buadc);
  TBranch *bvadc = 0;
  _t_ch->SetBranchAddress("VADCvec",&VADCvec,&bvadc);
  TBranch *byadc = 0;
  _t_ch->SetBranchAddress("YADCvec",&YADCvec,&byadc);
  TBranch *btdc = 0;
  _t_ch->SetBranchAddress("TDCvec",&TDCvec,&btdc);
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
    badc->GetEntry(j);
    buadc->GetEntry(j);
    bvadc->GetEntry(j);
    byadc->GetEntry(j);
    btdc->GetEntry(j);
    butdc->GetEntry(j);
    bvtdc->GetEntry(j);
    bytdc->GetEntry(j);
    bsse->GetEntry(j);
    ben->GetEntry(j);
    bnohits->GetEntry(j);
    int size = ADCvec->size();
    int no = NoHits;
    
    std::cout<<"Number of total hits: "<<size<<std::endl;

    c1->cd();
    // Plot hit time against hit amplitude
    TGraph *g = new TGraph(size,&((*TDCvec)[0]),&((*ADCvec)[0]));
    g->SetMarkerSize(0.7);
    g->SetMarkerColor(kBlue);
    g->GetXaxis()->SetTitle("TDC Time");
    g->GetYaxis()->SetTitle("ADC Amplitude");
    g->Draw("AP");
    g->GetXaxis()->SetLimits(0,10000);
    g->Draw("AP");
    c1->Modified();
    c1->Update();

    int usize = UADCvec->size();
    std::cout<<"Number of U hits: "<<usize<<std::endl;

    c2->cd();
    // Plot hit time against hit amplitude
    TGraph *g2 = new TGraph(usize,&((*UTDCvec)[0]),&((*UADCvec)[0]));
    g2->SetMarkerSize(0.7);
    g2->SetMarkerColor(kBlue);
    g2->GetXaxis()->SetTitle("TDC Time");
    g2->GetYaxis()->SetTitle("ADC Amplitude");
    g2->Draw("AP");
    g2->GetXaxis()->SetLimits(0,10000);
    g2->Draw("AP");
    c2->Modified();
    c2->Update();

    int vsize = VADCvec->size();
    std::cout<<"Number of V hits: "<<vsize<<std::endl;

    c3->cd();
    // Plot hit time against hit amplitude
    TGraph *g3 = new TGraph(vsize,&((*VTDCvec)[0]),&((*VADCvec)[0]));
    g3->SetMarkerSize(0.7);
    g3->SetMarkerColor(kBlue);
    g3->GetXaxis()->SetTitle("TDC Time");
    g3->GetYaxis()->SetTitle("ADC Amplitude");
    g3->Draw("AP");
    g3->GetXaxis()->SetLimits(0,10000);
    g3->Draw("AP");
    c3->Modified();
    c3->Update();

    int ysize = YTDCvec->size();
    std::cout<<"Number of Y hits: "<<ysize<<std::endl;

    c4->cd();
    // Plot hit time against hit amplitude
    TGraph *g4 = new TGraph(ysize,&((*YTDCvec)[0]),&((*YADCvec)[0]));
    g4->SetMarkerSize(0.7);
    g4->SetMarkerColor(kBlue);
    g4->GetXaxis()->SetTitle("TDC Time");
    g4->GetYaxis()->SetTitle("ADC Amplitude");
    g4->Draw("AP");
    g4->GetXaxis()->SetLimits(0,10000);
    g4->Draw("AP");
    c4->Modified();
    c4->Update();

    bool showerFlag = false;

    // Loop over number of showers
    for(int i = 0; i < ShowerStartEnd->size(); i++){
      std::cout<<"Shower!  Energy= "<<Evec->at(i)<<" MeV"<<std::endl;
      std::cout<<"Start: "<<(ShowerStartEnd->at(i)).first<<", End: "<<(ShowerStartEnd->at(i)).second<<std::endl;
      // For each shower loop over the hits
      for(int k = 0; k < ADCvec->size(); k++){
        // Check if hit occured within shower time window
        if((ShowerStartEnd->at(i)).first<(ShowerStartEnd->at(i)).second){
          if(TDCvec->at(k) > (ShowerStartEnd->at(i)).first-1 && TDCvec->at(k) < (ShowerStartEnd->at(i)).second+1 ){
            ShowerADC.push_back(ADCvec->at(k));
            ShowerTDC.push_back(TDCvec->at(k));
            showerFlag = true;
          }
        } else if((ShowerStartEnd->at(i)).first>(ShowerStartEnd->at(i)).second){
          if(TDCvec->at(k) < (ShowerStartEnd->at(i)).first+1 && TDCvec->at(k) > (ShowerStartEnd->at(i)).second-1 ){
            ShowerADC.push_back(ADCvec->at(k));
            ShowerTDC.push_back(TDCvec->at(k));
            showerFlag = true;
          }
        }
      }
    }
/*    
    // Mark the hits that occured within the shower time window
    if(showerFlag){
      TGraph *gshow = new TGraph(ShowerADC.size(),&(ShowerTDC[0]),&(ShowerADC[0]));
      gshow->SetMarkerSize(0.7);
      gshow->SetMarkerColor(kBlue);
      gshow->Draw("P");
      c4->Modified();
      c4->Update();
    }
*/    
    std::cout<<"Graph "<<j<<std::endl;
    std::cin.get();
    ShowerADC.clear(); ShowerTDC.clear();
  }
}
