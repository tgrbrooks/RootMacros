#include <string>
#include "SetCanvas.cxx"
#include "MultipleFluxPlotter.cxx"
/*
KamMaxENu.txt   KamMinaMuNu.txt  SudMaxaENu.txt   SudMaxMuNu.txt  SuMaxENu.txt   SuMinaMuNu.txt
KamMaxaENu.txt   KamMaxMuNu.txt  KamMinENu.txt    SudMaxaMuNu.txt  SuMaxaENu.txt   SuMaxMuNu.txt  SuMinENu.txt
KamMaxaMuNu.txt  KamMinaENu.txt  KamMinMuNu.txt   SudMaxENu.txt    SuMaxaMuNu.txt  SuMinaENu.txt  SuMinMuNu.txt
*/
PlotFluxes(){
  NumOfFiles = 4;
  // string FluxFileNames[] = {"New_SNO_MinENu.txt", "New_SNO_MinMuNu.txt","New_SNO_MinTNu.txt", "SNO_MinENu.txt","SNO_MinMuNu.txt"}; 
  string FluxFileNames[] = { "UnmixedFluxFiles/KAM_MinENuHE.txt", "UnmixedFluxFiles/KAM_MinMuNuHE.txt","MixedFluxFiles/New_KAM_MinENuHE.txt","MixedFluxFiles/New_KAM_MinMuNuHE.txt","MixedFluxFiles/New_KAM_MinEbNuHE.txt", "UnmixedFluxFiles/KAM_MinMubNuHE.txt", "UnmixedFluxFiles/KAM_MinEbNuHE.txt","MixedFluxFiles/New_KAM_MinMubNuHE.txt"};
  //string FluxFileNames[] = { "New_SNO_MinENu.txt", "New_SNO_MinMuNu.txt","New_SNO_MinTNu.txt"};
  //"New_SNO_MinENu.txt", "New_SNO_MinMuNu.txt","New_SNO_MinTNu.txt", "SNO_MinENu.txt","SNO_MinMuNu.txt"}; 

  //"New_SNO_MinMubNu.txt","New_SNO_MinEbNu.txt","New_SNO_MinTbNu.txt"
  SetCanvas();
  
  TLegend* leg1 = new TLegend(0.3,0.3,0.5,0.5);
  TLegend* leg2 = new TLegend(0.3,0.3,0.5,0.5);
  TLegend* leg3 = new TLegend(0.3,0.3,0.5,0.5);

  // TLegend leg[] = {new TLegend(0.3,0.3,0.5,0.5),new TLegend(0.3,0.3,0.5,0.5),new TLegend(0.3,0.3,0.5,0.5), new TLegend(0.3,0.3,0.5,0.5) } 

  for(unsigned int i(0); i < NumOfFiles; i ++){
    MultipleFluxPlotter(FluxFileNames[i],i+1,leg1,leg2,leg3); 
  }
}
