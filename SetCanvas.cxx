// Root Includes                                                                                                                                  
#include "TCanvas.h"
#include "TStyle"
#include "TLegend.h"
void SetCanvas(){
  
  TCanvas *c2 = new TCanvas("c2","Flux vs Energy",900,600);
  TCanvas *c3 = new TCanvas("c3","Flux vs Cos",900,600);
  TCanvas *c4 = new TCanvas("c4","Energy vs Flux",900,600);
  
}
