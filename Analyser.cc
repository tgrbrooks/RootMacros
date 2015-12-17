// Analyser/plotter with optimisation for single and double plane cuts

// C++ Includes
#include <iostream>
#include <string>
#include <vector>
#include <sstream>

// Root Includes
#include "THStack.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"

// Print text to file showing number of cut neutrinos of each flavour
void PrintText(std::ofstream & fileptr, std::vector<float> const & removedno, std::vector<float> const & nutypeno,float size){
    fileptr<<removedno[0]/size<<"% Total, "<<removedno[1]<<"/"<<nutypeno[0]<<" ("<<(removedno[1]/nutypeno[0])*size<<"%) CCQE, "
    <<removedno[2]<<"/"<<nutypeno[1]<<" ("<<(removedno[2]/nutypeno[1])*size<<"%) NCQE, "<<removedno[3]<<"/"<<nutypeno[2]<<" ("<<(removedno[3]/nutypeno[2])*size<<"%) CCRE, "
    <<removedno[4]<<"/"<<nutypeno[3]<<" ("<<(removedno[4]/nutypeno[3])*size<<"%) NCRE, "<<removedno[5]<<"/"<<nutypeno[4]<<" ("<<(removedno[5]/nutypeno[4])*size<<"%) CCDIS, "
    <<removedno[6]<<"/"<<nutypeno[5]<<" ("<<(removedno[6]/nutypeno[5])*size<<"%) NCDIS, "<<removedno[7]<<"/"<<nutypeno[6]<<" ("<<(removedno[7]/nutypeno[6])*size<<"%) CCCO, "
    <<removedno[8]<<"/"<<nutypeno[7]<<" ("<<(removedno[8]/nutypeno[7])*size<<"%) NCCO";
}

// Function to plot a stacked histogram of all neutrino energies split by their flavours
void DoubleHistType(TChain &NeutTree){
 
  // Create new canvas
  TCanvas *c2 = new TCanvas("SOMETHING","A NAME",900,600);
  c2->cd(); 
 
  // Get numbers corresponding to types as strings, Colours for each stacked plot  and list different 
  string Numbers[] = {"0","1","2","3","4","5","6","7"};
  string LabelNames[] = {"CCQE (40.0%)","NCQE (25.4%)","CCRE (13.8%)","NCRE (7.5%)","CCDIS (9.0%)","NCDIS (3.8%)","CCCO (0.4%)","NCCO (0.3%)"};
  string TypeNames[] = {"CCQE","NCQE","CCRE","NCRE","CCDIS","NCDIS","CCCO","NCCO"};
  short Colours[] = {kRed-4, kYellow-4, kGreen-4, kCyan-4, kBlue-4, kMagenta-4, kMagenta+2, kMagenta+4};  
  string CutStrings [8];

  // Initialise stack and legend objects
  THStack *hs = new THStack("hs","Stacked 1D histograms");
  TLegend *leg = new TLegend(0.4,0.4,0.7,0.7);

  // Loop over all neutrino flavours and define the "Cut String" which tells the TTree::Draw() function
  // which neutrino flavour to plot for
  for(unsigned int i(0); i < 8; i ++){
    CutStrings[i] = "Type == " + Numbers[i];
    NeutTree.SetLineColor(kBlack);  NeutTree.SetFillColor(Colours[i]);  
    NeutTree.Draw(("Energy>>"+TypeNames[i]+"(50,0,10)").c_str(),CutStrings[i].c_str(),"",10000000000,0);
  }
   
  // Access each histogram, add a legend entry for it and add it to the stack of histograms hs.
  leg->AddEntry(CCQE,LabelNames[0].c_str(),"f");  hs->Add(CCQE);
  leg->AddEntry(NCQE,LabelNames[1].c_str(),"f");  hs->Add(NCQE);
  leg->AddEntry(CCRE,LabelNames[2].c_str(),"f");  hs->Add(CCRE);
  leg->AddEntry(NCRE,LabelNames[3].c_str(),"f");  hs->Add(NCRE);
  leg->AddEntry(CCDIS,LabelNames[4].c_str(),"f"); hs->Add(CCDIS);
  leg->AddEntry(NCDIS,LabelNames[5].c_str(),"f"); hs->Add(NCDIS);
  leg->AddEntry(CCCO,LabelNames[6].c_str(),"f");  hs->Add(CCCO);
  leg->AddEntry(NCCO,LabelNames[7].c_str(),"f");  hs->Add(NCCO);
 
  // Draw stacked histograms and adjust axis
  hs->Draw();
  hs->SetMaximum(2100);// CURRENTLY THIS IS SET MANUALLY!
  hs->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  hs->GetYaxis()->SetTitle("Events (/0.2 GeV)");
  hs->GetXaxis()->SetTitleOffset(0.9);
  hs->GetYaxis()->SetTitleOffset(0.9);
  hs->Draw();

  // Draw legend and update canvas
  leg->Draw();  c2->Update();  c2->Modified();
}

// Function to plot hist
void DoubleHist(TChain &nnbarTree, TChain &NeutTree, const string CutVariable, const string AxisLabel ){

  NeutTree.SetLineColor(kRed);
  nnbarTree.SetLineColor(kBlue);
  
  if(CutVariable == "ADCampV" || CutVariable == "ADCampY" || CutVariable == "ADCamp" || CutVariable == "ADCampU" ){ 

    nnbarTree.Draw(CutVariable.c_str(),"","E1HIST",10000000000,0);
    NeutTree.Draw(CutVariable.c_str(),"","E1HISTsame",10000000000,0);

  } else {
 
    NeutTree.Draw(CutVariable.c_str(),"","E1HIST",10000000000,0);
    nnbarTree.Draw(CutVariable.c_str(),"","E1HISTsame",10000000000,0);
  }

  htemp->GetXaxis()->SetTitle(AxisLabel.c_str());
  htemp->GetXaxis()->SetTitleOffset(0.8);
  htemp->GetYaxis()->SetTitleOffset(0.9);
  //htemp->GetXaxis()->SetRange(0,1000);
  double per = (htemp->GetXaxis()->GetXmax())/(htemp->GetSize());
  std::ostringstream s;
  s << "Events /(" << per << ")";
  htemp->GetYaxis()->SetTitle((s.str()).c_str());

  leg = new TLegend(0.3,0.3,0.5,0.5);
  TLegendEntry *le1 = leg->AddEntry("nnbarTree","n#bar{n}","l"); le1->SetLineColor(kBlue); le1->SetLineWidth(2);
  TLegendEntry *le2 = leg->AddEntry("NeutTree","NuMI Beam","l"); le2->SetLineColor(kRed); le2->SetLineWidth(2);
  leg->Draw();
}

// Function to plot 2D hist for two plane cuts
void twoDoubleHist(TChain &nnbarTree, TChain &NeutTree,const string CutVar1, const string CutVar2, const string AxisLabel1, const string AxisLabel2){

  nnbarTree.SetLineColor(kBlue);
  nnbarTree.Draw((CutVar1+":"+CutVar2).c_str(),"","box",10000000000,0);

  htemp->GetXaxis()->SetTitle(AxisLabel1.c_str());
  htemp->GetYaxis()->SetTitle(AxisLabel2.c_str());
  htemp->GetXaxis()->SetTitleOffset(0.75);
  htemp->GetYaxis()->SetTitleOffset(0.9);

  NeutTree.SetLineColor(kRed);
  NeutTree.Draw((CutVar1+":"+CutVar2).c_str(),"","boxsame",10000000000,0); // 
  
  leg = new TLegend(0.3,0.3,0.5,0.5);
  TLegendEntry *le1 = leg->AddEntry("nnbarTree","n#bar{n}","l"); le1->SetLineColor(kBlue); le1->SetLineWidth(2);
  TLegendEntry *le2 = leg->AddEntry("NeutTree","NuMI Beam","l"); le2->SetLineColor(kRed); le2->SetLineWidth(2);
  leg->SetFillStyle(0);
  leg->Draw();
}

float CalcMean(const vector <float> &Vec){
  float Mean = 0.0;
  for (unsigned int i =0; i < Vec.size(); i++){
    Mean += Vec[i];
  }
  return (Mean/Vec.size());
}

float CalcStandardDev(const vector <float> &Vec){
  float StandardDeviation = 0.0;
  float Mean = CalcMean(Vec);
  float sum = 0.0;
  float var = 0.0;
       
  for (unsigned int i =0; i < Vec.size(); i++){
    sum += pow((Vec[i] - Mean),2);
  }
  var = sqrt(sum/Vec.size());
  return var;
}

// Function to read in the ttree's from the input files and store the info in vector of vector of floats
vector<vector<float>> ReadTree(string RootFileName, bool isNeut){

  // Initilise float variables, vector of floats and vector of vector of floats
  float TDCstd, TDCstdU,TDCstdV,TDCstdY,TDCiqr,TDCiqrU,TDCiqrV,TDCiqrY,ADCamp,ADCampU,ADCampV,ADCampY,WFint,WFintU,WFintV,WFintY,Energy;
  vector <float> v_TDCstd, v_TDCstdU,v_TDCstdV,v_TDCstdY,v_TDCiqr,v_TDCiqrU,v_TDCiqrV,v_TDCiqrY,v_ADCamp,v_ADCampU,v_ADCampV,v_ADCampY,v_WFint,v_WFintU,v_WFintV,v_WFintY,v_Energ,v_Hits,v_UHits,v_VHits,v_YHits,v_Type,v_Energy;
  vector < vector < float > > vv_ReturnVect;

  UInt_t Hits,UHits,VHits,YHits,Type;

  TFile *TRootFile = new TFile(RootFileName.c_str());  
  TTree *FileTree = (TTree*)TRootFile->Get("ch_tree");
 
  // Set Branch Addresses to floats
  FileTree->SetBranchAddress("hitNo",&Hits);
  FileTree->SetBranchAddress("hitNoU",&UHits);
  FileTree->SetBranchAddress("hitNoV",&VHits);
  FileTree->SetBranchAddress("hitNoY",&YHits);
  FileTree->SetBranchAddress("TDCstd", &TDCstd);
  FileTree->SetBranchAddress("TDCstdU", &TDCstdU);
  FileTree->SetBranchAddress("TDCstdV", &TDCstdV);
  FileTree->SetBranchAddress("TDCstdY", &TDCstdY);
  FileTree->SetBranchAddress("TDCiqr", &TDCiqr);
  FileTree->SetBranchAddress("TDCiqrU", &TDCiqrU);
  FileTree->SetBranchAddress("TDCiqrV", &TDCiqrV);
  FileTree->SetBranchAddress("TDCiqrY", &TDCiqrY);
  FileTree->SetBranchAddress("ADCamp", &ADCamp);
  FileTree->SetBranchAddress("ADCampU", &ADCampU);
  FileTree->SetBranchAddress("ADCampV", &ADCampV);
  FileTree->SetBranchAddress("ADCampY", &ADCampY);
  FileTree->SetBranchAddress("WFint", &WFint);
  FileTree->SetBranchAddress("WFintU", &WFintU);
  FileTree->SetBranchAddress("WFintV", &WFintV);
  FileTree->SetBranchAddress("WFintY", &WFintY);
  if( isNeut == true ){ // Additional data if it's Neutrino
    FileTree->SetBranchAddress("Energy", &Energy);
    FileTree->SetBranchAddress("Type",&Type);
  }
  Long64_t EntryNumber = FileTree->GetEntries();
  for(unsigned int j=0; j!=EntryNumber; ++j){ // Loop over tree entries
    // Fill vectors of floats with each branch address data member.
    FileTree->GetEntry(j);
    v_Hits.push_back(Hits);
    v_UHits.push_back(UHits);
    v_VHits.push_back(VHits);
    v_YHits.push_back(YHits);
    v_TDCstd.push_back(TDCstd);
    v_TDCstdU.push_back(TDCstdU);
    v_TDCstdV.push_back(TDCstdV);
    v_TDCstdY.push_back(TDCstdY);
    v_TDCiqr.push_back(TDCiqr);
    v_TDCiqrU.push_back(TDCiqrU);
    v_TDCiqrV.push_back(TDCiqrV);
    v_TDCiqrY.push_back(TDCiqrY);
    v_ADCamp.push_back(ADCamp);
    v_ADCampU.push_back(ADCampU);
    v_ADCampV.push_back(ADCampV);
    v_ADCampY.push_back(ADCampY);
    v_WFint.push_back(WFint);
    v_WFintU.push_back(WFintU);
    v_WFintV.push_back(WFintV);
    v_WFintY.push_back(WFintY);

    if( isNeut == true ){ // Additional data if it's Neutrino
      v_Energy.push_back(Energy);
      v_Type.push_back(Type);
    }
  }
  // Push back data read in from tree onto vector of vectors. Columns of 2D vector correspond to events 1,2,3,...,N 
  // While rows correspond to the different cut variables.
  vv_ReturnVect.push_back(v_Hits);   vv_ReturnVect.push_back(v_UHits);   vv_ReturnVect.push_back(v_VHits);   vv_ReturnVect.push_back(v_YHits);
  vv_ReturnVect.push_back(v_TDCstd); vv_ReturnVect.push_back(v_TDCstdU); vv_ReturnVect.push_back(v_TDCstdV); vv_ReturnVect.push_back(v_TDCstdY);
  vv_ReturnVect.push_back(v_TDCiqr); vv_ReturnVect.push_back(v_TDCiqrU); vv_ReturnVect.push_back(v_TDCiqrV); vv_ReturnVect.push_back(v_TDCiqrY);
  vv_ReturnVect.push_back(v_ADCamp); vv_ReturnVect.push_back(v_ADCampU); vv_ReturnVect.push_back(v_ADCampV); vv_ReturnVect.push_back(v_ADCampY);
  vv_ReturnVect.push_back(v_WFint);  vv_ReturnVect.push_back(v_WFintU);  vv_ReturnVect.push_back(v_WFintV);  vv_ReturnVect.push_back(v_WFintY);
  if ( isNeut == true ){
    vv_ReturnVect.push_back(v_Energy); vv_ReturnVect.push_back(v_Type);
  }

  return vv_ReturnVect;

} // End of ReadTree

// Function to cut on all of the variables independently
// Different standard deviations for optimisations, rpm, suv, suy, svy for output run
vector<vector<float>> SingleCut(const vector<vector<float>> & vv_ReturnVect, bool isNeut, const vector<float> & MEANS, 
                                const vector<float> & STDEVP1P, const vector<float> & STDEVP1M, const vector<float> & STDEVP2P, 
                                const vector<float> & STDEVP2M, vector<vector<float>> Rpm, vector<vector<float>> suv, vector<vector<float>> suy, vector<vector<float>> svy){

  // Total cut
  vector<vector<float>> RemovedType;
  vector<float> EmptyVect; 
  for(unsigned int i(0); i < 9; i ++) EmptyVect.push_back(0);

  for(unsigned int i(0); i < 20; i ++){
    RemovedType.push_back(EmptyVect); // Removed Type After this loop will be of dimensions RemovedType[20][8]

    for(unsigned int j(0); j < vv_ReturnVect[0].size(); j ++){ // LOOP NUMBER DICTATED BY v_Hits.size()
      // Multiply sigma by optimised +/-
      if( vv_ReturnVect[i][j] > (MEANS[i] + STDEVP1P[i]*Rpm[0][i]) || vv_ReturnVect[i][j] < (MEANS[i] - STDEVP1M[i]*Rpm[1][i]) ){ 
        RemovedType[i][0] ++;
        if(isNeut == true){
          for(unsigned int k=0; k < 8; k ++){
            if(vv_ReturnVect[21][j]==k) RemovedType[i][k+1] ++;
          }
        }
      }
    }
  }
  
  for(unsigned int j(0); j < 20; j++) RemovedType.push_back(EmptyVect); // Removed Type After this loop will be of dimensions RemovedType[40][8]
  for(unsigned int p(0); p < 5; p ++){
    // Index == (1+p*4) takes values 1,5,9,13,17. These correspond to vector elements v_UHits, v_TDCstdU, v_TDCiqrU, v_ADCampU and v_WFintU
    // Index == (2+p*4) takes values 2,6,10,14,18. These correspond to vector elements v_VHits, v_TDCstdV, v_TDCiqrV, v_ADCampV and v_WFintV
    // Index == (3+p*4) takes values 3,7,11,15,19. These correspond to vector elements v_YHits, v_TDCstdY, v_TDCiqrY, v_ADCampY and v_WFintY
    unsigned int UPI = (1+p*4); // UPlaneIndex
    unsigned int VPI = (2+p*4); // VPlaneIndex
    unsigned int YPI = (3+p*4); // YPlaneIndex
    for(unsigned int j(0); j < vv_ReturnVect[0].size(); j ++){ // LOOP NUMBER DICTATED BY v_Hits.size()

      // UV Plane Cut -- multiply sigma+ by +/- ratio and optimised no of sigma for uv plane
      if(
        (vv_ReturnVect[UPI][j]>(MEANS[UPI] + STDEVP1P[UPI]*(Rpm[0][UPI]/Rpm[1][UPI])*suv[0][p])||vv_ReturnVect[UPI][j]<(MEANS[UPI] - STDEVP1M[UPI]*suv[0][p])) ||
        (vv_ReturnVect[VPI][j]>(MEANS[VPI] + STDEVP2P[VPI]*(Rpm[0][VPI]/Rpm[1][VPI])*suv[1][p])||vv_ReturnVect[VPI][j]<(MEANS[VPI] - STDEVP2M[VPI]*suv[1][p]))){
        RemovedType[20+p][0] ++; // A Neutrino was removed from UV plane cut
	if(isNeut == true){
          for(unsigned int k=0; k<8;k++){
            if(vv_ReturnVect[21][j]==k){RemovedType[20+p][k+1] ++;} // Find flavour of neutrino
          }
        }
      }
      
      // UY Plane Cut
      if(
        (vv_ReturnVect[UPI][j]>(MEANS[UPI] + STDEVP1P[UPI]*(Rpm[0][UPI]/Rpm[1][UPI])*suy[0][p])||vv_ReturnVect[UPI][j]<(MEANS[UPI] - STDEVP1M[UPI]*suy[0][p])) ||
        (vv_ReturnVect[YPI][j]>(MEANS[YPI] + STDEVP2P[YPI]*(Rpm[0][YPI]/Rpm[1][YPI])*suy[1][p])||vv_ReturnVect[YPI][j]<(MEANS[YPI] - STDEVP2M[YPI]*suy[1][p]))){ 
        RemovedType[25+p][0] ++; // A Neutrino was removed from UY plane cut
        if(isNeut == true){
	  for(unsigned int k=0; k<8;k++){
	    if(vv_ReturnVect[21][j]==k){RemovedType[25+p][k+1] ++;} // Find flavour of neutrino
          }
	}
      }
    
      // VY Plane Cut
      if(
        (vv_ReturnVect[VPI][j]>(MEANS[VPI] + STDEVP1P[VPI]*(Rpm[0][VPI]/Rpm[1][VPI])*svy[0][p])||vv_ReturnVect[VPI][j]<(MEANS[VPI] - STDEVP1M[VPI]*svy[0][p])) ||
        (vv_ReturnVect[YPI][j]>(MEANS[YPI] + STDEVP2P[YPI]*(Rpm[0][YPI]/Rpm[1][YPI])*svy[1][p])||vv_ReturnVect[YPI][j]<(MEANS[YPI] - STDEVP2M[YPI]*svy[1][p]))){
        RemovedType[30+p][0] ++; // A Neutrino was removed from VY plane cut
	if(isNeut == true){
          for(unsigned int k=0; k<8;k++){
	    if(vv_ReturnVect[21][j]==k){RemovedType[30+p][k+1] ++;} // Find flavour of neutrino
          }
	}
      }
     
      // UVY Plane Cut
      if(
        (vv_ReturnVect[UPI][j]>(MEANS[UPI] + STDEVP1P[UPI])||vv_ReturnVect[UPI][j]<(MEANS[UPI] - STDEVP1M[UPI])) ||
        (vv_ReturnVect[VPI][j]>(MEANS[VPI] + STDEVP2P[VPI])||vv_ReturnVect[VPI][j]<(MEANS[VPI] - STDEVP2M[VPI])) ||
        (vv_ReturnVect[YPI][j]>(MEANS[YPI] + STDEVP2P[YPI])||vv_ReturnVect[YPI][j]<(MEANS[YPI] - STDEVP2M[YPI]))){
        RemovedType[35+p][0] ++; // A Neutrino was removed from UVY plane cut
	if(isNeut == true){
          for(unsigned int k=0; k<8;k++){
            if(vv_ReturnVect[21][j]==k){RemovedType[35+p][k+1] ++;} // Find flavour of neutrino
          }
	}
      }
    }
  }
  if(isNeut == true){
    for(unsigned int i=0; i < vv_ReturnVect[21].size(); i++){
      for(unsigned int j=0;j<8;j++){
        if(vv_ReturnVect[21][i] == j) EmptyVect[j] ++; // Loop over all neutrino types. See vector initialisation above.
      }
    }
    RemovedType.push_back(EmptyVect);
  }
  return RemovedType;

} // End of SingleCut

// Fnction to combine two cuts together
vector<vector<float>> CombineTwoCuts(const vector<vector<float>> & vv_ReturnVect, bool isNeut, const vector<float> & MEAN, 
                                     const vector<float> & STDEV1, const vector<float> & STDEV2, int cut1, int cut2){//cut=0,4,8,12,16 (hn,ts,ti,aa,wi)

  // Count number of removed events outside of some tolerance of hit wires
  // Total cut
  vector<vector<float>> RemovedType;
  vector<float> EmptyVect; 
  for(unsigned int i(0); i < 9; i ++) EmptyVect.push_back(0);

  for(unsigned int i(0); i < 4; i ++){
    RemovedType.push_back(EmptyVect);

    for(unsigned int j(0); j < vv_ReturnVect[0].size(); j ++){ // LOOP NUMBER DICTATED BY v_Hits.size()

      if( vv_ReturnVect[i+cut1][j] > (MEAN[i+cut1] + STDEV1[i+cut1]) || vv_ReturnVect[i+cut1][j] < (MEAN[i+cut1] - STDEV1[i+cut1]) || 
          vv_ReturnVect[i+cut2][j] > (MEAN[i+cut2] + STDEV2[i+cut2]) || vv_ReturnVect[i+cut2][j] < (MEAN[i+cut2] - STDEV2[i+cut2]) ){ 
        RemovedType[i][0] ++;
	if(isNeut == true){
          for(unsigned int k=0; k < 8; k ++){
            if(vv_ReturnVect[21][j]==k) RemovedType[i][k+1] ++;
          }
	}
      }       
    }
  }

  RemovedType.push_back(EmptyVect);

  for(unsigned int j(0); j < vv_ReturnVect[0].size(); j ++){ // LOOP NUMBER DICTATED BY v_Hits.size()

    // UV Plane Cut
    if(
      (vv_ReturnVect[1+cut1][j]>(MEAN[1+cut1] + STDEV1[1+cut1])||vv_ReturnVect[1+cut1][j]<(MEAN[1+cut1] - STDEV1[1+cut1])) ||
      (vv_ReturnVect[2+cut1][j]>(MEAN[2+cut1] + STDEV1[2+cut1])||vv_ReturnVect[2+cut1][j]<(MEAN[2+cut1] - STDEV1[2+cut1])) ||
      (vv_ReturnVect[1+cut2][j]>(MEAN[1+cut2] + STDEV2[1+cut2])||vv_ReturnVect[1+cut2][j]<(MEAN[1+cut2] - STDEV2[1+cut2])) ||
      (vv_ReturnVect[2+cut2][j]>(MEAN[2+cut2] + STDEV2[2+cut2])||vv_ReturnVect[2+cut2][j]<(MEAN[2+cut2] - STDEV2[2+cut2]))){
      RemovedType[4][0] ++; // A Neutrino was removed from UV plane cut
      if(isNeut == true){
        for(unsigned int k=0; k<8;k++){
	  if(vv_ReturnVect[21][j]==k){RemovedType[4][k+1] ++;} // Find flavour of neutrino
        }
      }
    }
  }

  RemovedType.push_back(EmptyVect);

  for(unsigned int j(0); j < vv_ReturnVect[0].size(); j ++){ // LOOP NUMBER DICTATED BY v_Hits.size()

    // UY Plane Cut
    if(
      (vv_ReturnVect[1+cut1][j]>(MEAN[1+cut1] + STDEV1[1+cut1])||vv_ReturnVect[1+cut1][j]<(MEAN[1+cut1] - STDEV1[1+cut1])) ||
      (vv_ReturnVect[3+cut1][j]>(MEAN[3+cut1] + STDEV1[3+cut1])||vv_ReturnVect[3+cut1][j]<(MEAN[3+cut1] - STDEV1[3+cut1])) ||
      (vv_ReturnVect[1+cut2][j]>(MEAN[1+cut2] + STDEV2[1+cut2])||vv_ReturnVect[1+cut2][j]<(MEAN[1+cut2] - STDEV2[1+cut2])) ||
      (vv_ReturnVect[3+cut2][j]>(MEAN[3+cut2] + STDEV2[3+cut2])||vv_ReturnVect[3+cut2][j]<(MEAN[3+cut2] - STDEV2[3+cut2]))){
      RemovedType[5][0] ++; // A Neutrino was removed from UV plane cut
      if(isNeut == true){
        for(unsigned int k=0; k<8;k++){
	  if(vv_ReturnVect[21][j]==k){RemovedType[5][k+1] ++;} // Find flavour of neutrino
        }
      }
    }
  }

  RemovedType.push_back(EmptyVect);

  for(unsigned int j(0); j < vv_ReturnVect[0].size(); j ++){ // LOOP NUMBER DICTATED BY v_Hits.size()

    // VY Plane Cut
    if(
      (vv_ReturnVect[2+cut1][j]>(MEAN[2+cut1] + STDEV1[2+cut1])||vv_ReturnVect[2+cut1][j]<(MEAN[2+cut1] - STDEV1[2+cut1])) ||
      (vv_ReturnVect[3+cut1][j]>(MEAN[3+cut1] + STDEV1[3+cut1])||vv_ReturnVect[3+cut1][j]<(MEAN[3+cut1] - STDEV1[3+cut1])) ||
      (vv_ReturnVect[2+cut2][j]>(MEAN[2+cut2] + STDEV2[2+cut2])||vv_ReturnVect[2+cut2][j]<(MEAN[2+cut2] - STDEV2[2+cut2])) ||
      (vv_ReturnVect[3+cut2][j]>(MEAN[3+cut2] + STDEV2[3+cut2])||vv_ReturnVect[3+cut2][j]<(MEAN[3+cut2] - STDEV2[3+cut2]))){
      RemovedType[6][0] ++; // A Neutrino was removed from UV plane cut
      if(isNeut == true){
        for(unsigned int k=0; k<8;k++){
	  if(vv_ReturnVect[21][j]==k){RemovedType[6][k+1] ++;} // Find flavour of neutrino
        }
      }
    }
  }

  RemovedType.push_back(EmptyVect);

  for(unsigned int j(0); j < vv_ReturnVect[0].size(); j ++){ // LOOP NUMBER DICTATED BY v_Hits.size()

    // UVY Plane Cut
    if(
      (vv_ReturnVect[1+cut1][j]>(MEAN[1+cut1] + STDEV1[1+cut1])||vv_ReturnVect[1+cut1][j]<(MEAN[1+cut1] - STDEV1[1+cut1])) ||
      (vv_ReturnVect[2+cut1][j]>(MEAN[2+cut1] + STDEV1[2+cut1])||vv_ReturnVect[2+cut1][j]<(MEAN[2+cut1] - STDEV1[2+cut1])) ||
      (vv_ReturnVect[3+cut1][j]>(MEAN[3+cut1] + STDEV1[3+cut1])||vv_ReturnVect[3+cut1][j]<(MEAN[3+cut1] - STDEV1[3+cut1])) ||
      (vv_ReturnVect[1+cut2][j]>(MEAN[1+cut2] + STDEV2[1+cut2])||vv_ReturnVect[1+cut2][j]<(MEAN[1+cut2] - STDEV2[1+cut2])) ||
      (vv_ReturnVect[2+cut2][j]>(MEAN[2+cut2] + STDEV2[2+cut2])||vv_ReturnVect[2+cut2][j]<(MEAN[2+cut2] - STDEV2[2+cut2])) ||
      (vv_ReturnVect[3+cut2][j]>(MEAN[3+cut2] + STDEV2[3+cut2])||vv_ReturnVect[3+cut2][j]<(MEAN[3+cut2] - STDEV2[3+cut2]))){
      RemovedType[7][0] ++; // A Neutrino was removed from UV plane cut
      if(isNeut == true){
        for(unsigned int k=0; k<8;k++){
	  if(vv_ReturnVect[21][j]==k){RemovedType[7][k+1] ++;} // Find flavour of neutrino
        }
      }
    }
  }
  return RemovedType;

} //End of CombineTwoCuts

// Function to optimise the plot limits for single planes
vector<vector<float>> OptimisePlusMinus(const vector<vector<float>> & NeutEvents, const vector<vector<float>> & NnbarEvents, 
                                const vector<float> & MEAN, const vector<float> & STDEV, float limit){
// plane = 1(u), 2(v), 3(y)
// cut = 0(hn), 4(ts), 8(ti), 12(aa), 16(wf)
  vector<vector<float>> plusminus;
  
  int lim = 15;
  vector<float> SigmaP, SigmaM;
  // Not needed for optimisation but required for output run
  vector<float> dummyratio;
  vector<vector<float>> dummyplanes;
  // sP = number of sigma above mean, sM = number of sigma below mean, maxrem = max number of neutrino events removed
  vector<float> sP, sM, maxrem;
  for(unsigned int k(0); k<5; k ++){
    for(unsigned int n(0); n < 4; n ++){
      sP.push_back(0); sM.push_back(0); maxrem.push_back(0); dummyratio.push_back(1);
    }
  }
  dummyplanes.push_back(dummyratio); dummyplanes.push_back(dummyratio);

  // Fill the sigmas with zeros
  for(unsigned int i(0); i < STDEV.size(); i ++){ SigmaP.push_back(0.0); SigmaM.push_back(0.0); }
  // loop over different combinations of sigma above and below mean
  for(unsigned int l(5); l < lim; l ++){
    for(unsigned int m(5); m < lim; m ++){
      for(unsigned int i(0); i < STDEV.size(); i ++){ SigmaP[i] = STDEV[i]*(Double_t)(l)/5; SigmaM[i] = STDEV[i]*(Double_t)(m)/5; }
      // Calculate the number of different events removed
      vector<vector<float>> NeutRemoved = SingleCut(NeutEvents,true,MEAN,SigmaP,SigmaM,SigmaP,SigmaM,dummyplanes,dummyplanes,dummyplanes,dummyplanes);
      vector<vector<float>> NnbarRemoved = SingleCut(NnbarEvents,false,MEAN,SigmaP,SigmaM,SigmaP,SigmaM,dummyplanes,dummyplanes,dummyplanes,dummyplanes);
      std::cout<<(l-5)*10+(m-5)<<std::endl;
      for(unsigned int k(0); k<5; k ++){
        int cut = k*4;
        for(unsigned int n(0); n < 4; n ++){
          int plane = n;
          // Find combination of sP and sM that maximises the number of neutrinos removed while 95% of nnbar remain
          if (NnbarRemoved[plane+cut][0]/100 < limit && NeutRemoved[plane+cut][0]/100 > maxrem[plane+cut]){
            maxrem[plane+cut] = NeutRemoved[plane+cut][0]/100; sP[plane+cut] = (Double_t)(l)/5; sM[plane+cut] = (Double_t)(m)/5;}
        }
      }
    }
  }

  plusminus.push_back(sP);
  plusminus.push_back(sM);
  return plusminus; // returns a 2x20 vector of vector
  
} // End of OptimisePlusMinus

// Function to optimise two plane cuts
vector<vector<float>> OptimiseTwoPlanes(const vector<vector<float>> & NeutEvents, const vector<vector<float>> & NnbarEvents, int p1, int p2,
                                        const vector<float> & MEAN, const vector<float> & STDEV, const vector<float> & PMRatio, float limit){
  //p1/2 = 1(u), 2(v), 3(y)
  vector<vector<float>> plane1plane2;

  // Not needed for optimisation
  vector<float> dummyratio;
  vector<vector<float>> dummyplanes;
  for(unsigned int k(0); k<5; k ++){
    for(unsigned int n(0); n < 4; n ++){
      dummyratio.push_back(1);
    }
  }
  dummyplanes.push_back(dummyratio); dummyplanes.push_back(dummyratio);
  
  int lim = 15;
  vector<float> SigmaP1P, SigmaP1M, SigmaP2P, SigmaP2M;
  // Record the sigmas for each plane (+/- can be obtained using the ratio)
  vector<float> sP1, sP2, maxrem;
  for(unsigned int n(0); n < 5; n ++){
    sP1.push_back(0); sP2.push_back(0); maxrem.push_back(0);
  }

  // Fill the sigmas with zeros
  for(unsigned int i(0); i < STDEV.size(); i ++){ SigmaP1P.push_back(0.0); SigmaP1M.push_back(0.0); SigmaP2P.push_back(0.0); SigmaP2M.push_back(0.0); }
  // Loop over different combinations of sigma for the two different planes
  for(unsigned int l(5); l < lim; l ++){
    for(unsigned int m(5); m < lim; m ++){
      std::cout<<(l-5)*10+(m-5)<<std::endl;
      // Calculate the +/- using the ratios calculated from OptimisePlusMinus
      for(unsigned int i(0); i < STDEV.size(); i ++){SigmaP1M[i] = STDEV[i]*(Double_t)(l)/5; SigmaP1P[i] = PMRatio[i]*SigmaP1M[i]; 
                                                     SigmaP2M[i] = STDEV[i]*(Double_t)(m)/5; SigmaP2P[i] = PMRatio[i]*SigmaP2M[i]; }
      vector<vector<float>> NeutRemoved = SingleCut(NeutEvents,true,MEAN,SigmaP1P,SigmaP1M,SigmaP2P,SigmaP2M,dummyplanes,dummyplanes,dummyplanes,dummyplanes);
      vector<vector<float>> NnbarRemoved = SingleCut(NnbarEvents,false,MEAN,SigmaP1P,SigmaP1M,SigmaP2P,SigmaP2M,dummyplanes,dummyplanes,dummyplanes,dummyplanes);
      for(unsigned int k(0); k<5; k ++){
        int cut = k;
        int plane;
        if ((p1==1&&p2==2)||(p2==1&&p1==2)){plane=20;}else if ((p1==1&&p2==3)||(p2==1&&p1==3)){plane=25;}
        else if ((p1==2&&p2==3)||(p2==2&&p1==3)){plane=30;}else{cerr << "\n\nError! You didn't enter a valid number!\n"; exit(1);}
        // Maximise number removed keeping 95% nnbar
        if (NnbarRemoved[plane+cut][0]/100 < limit && NeutRemoved[plane+cut][0]/100 > maxrem[k]){
          maxrem[k] = NeutRemoved[plane+cut][0]/100; sP1[k] = (Double_t)(l)/5; sP2[k] = (Double_t)(m)/5;}
      }
    }
  }

  plane1plane2.push_back(sP1);
  plane1plane2.push_back(sP2);
  return plane1plane2;//2x5 vector
  
}

// Main Script
void Analyser(){
  // String list of cut variables, same as histogram names
  string AxisLabels[20] = {
    "Hit Number","U Hit Number","V Hit Number","Y Hit Number","TDC Std Dev (0.5 #mus)","U TDC Std Dev (0.5 #mus)",
    "V TDC Std Dev (0.5 #mus)","Y TDC Std Dev (0.5 #mus)","TDC IQRange (0.5 #mus)","U TDC IQRange (0.5 #mus)","V TDC IQRange (0.5 #mus)",
    "Y TDC IQRange (0.5 #mus)","ADC Amplitude (arb)","U ADC Amplitude (arb)","V ADC Amplitude (arb)","Y ADC Amplitude (arb)",
    "Integrated Waveform (arb)","U Integrated Waveform (arb)","V Integrated Waveform (arb)","Y Integrated WAveform (arb)"};

  string CutVariables[20] = {
    "hitNo","hitNoU","hitNoV","hitNoY","TDCstd","TDCstdU","TDCstdV","TDCstdY","TDCiqr","TDCiqrU","TDCiqrV","TDCiqrY",
    "ADCamp","ADCampU","ADCampV","ADCampY","WFint","WFintU","WFintV","WFintY"};

  // List nnbar and Neutrino files to be merged
  const int NumNbarFiles = 100;  string nnbarFileNames[NumNbarFiles];
  const int NumNeutFiles = 100; string NeutFileNames[NumNeutFiles];
  
  // Define String stream and temporary string to hold the file number
  stringstream ss; string FileNumAsString;
  // Loop over the number of files and allocate the file names to string vectors.
  // Use string stream to get the file number as a string
  for(unsigned int i(0); i < NumNbarFiles; i ++){
    ss << i; ss >> FileNumAsString; // Convert integer i of loop to string
    nnbarFileNames[i] = "outputs_run2/Nnbar_output_run2_" + FileNumAsString + ".root";
    NeutFileNames[i] = "outputs_run2/Atmo_output_run2_" + FileNumAsString + ".root";
    ss.str(""); ss.clear(); // Clear the string stream
  }

  gROOT->ProcessLine(".x lhcbStyle.C"); // Using lhcb Style file
 
  std::ofstream myfile; // Create output file stream with name AnalyserDataOutput.txt
  myfile.open("DataOutputAnalyser.txt");
  
  // Create two TChains, one for nnbar trees and one for Neutrino data trees
  TChain nnbarTrees("ch_tree");
  TChain NeutTrees("ch_tree");
  
  std::cout << "\nBeginning Analyser main script. Now reading in files (Check Output file)\n";
  
  // Message to output file
  myfile << "||Analyser.cc|| Begin.\n||Analyser.cc|| The Folllowing " << NumNbarFiles << " NNbar Files were merged:\n";  

  // Loop over nnbar files and add them to TChain
  for(unsigned int i(0); i < NumNbarFiles; i ++){
    nnbarTrees.Add(nnbarFileNames[i].c_str()); 
    myfile << nnbarFileNames[i] << "\n"; // Print File names to outputfile
  }
  std::cout << "\nNow merging together nnbarfiles into single file allnnbar.root.\n";
  nnbarTrees.Merge("allnnbar.root"); // Merge all Chained trees into a single root file
  
  std::cout << "\nnnbar Files successfully read and merged. Now attempting to gather data from chained file.\n";
   
  vector <float> StandardDevs, Means, Sigma, StandardDevP, StandardDevM;

  // Use ReadTree to return vector of vectors that contains all the information for all the events for nnbar
  // Only grandmaster luke understands the indexing
  vector<vector<float>> NnbarEvents = ReadTree("allnnbar.root",false);
  
  // Output message to file
  myfile << "\n\nMean Values and standard Deviations of Merged nnbarfile for different cut variables\n";

  std::cout << "\nData successfully read. Now calculating Standard Deviations and Means.\n";

  // Loop over data returned in NnbarEvents given by ReadTree() to get the Means and StandardDeviations
  for(unsigned int i(0); i < NnbarEvents.size(); i ++){ 
    float tempfloat = CalcStandardDev(NnbarEvents[i]);
    StandardDevs.push_back( tempfloat ); 
    Sigma.push_back( tempfloat );
    StandardDevP.push_back( tempfloat );
    StandardDevM.push_back( tempfloat );
    Means.push_back( CalcMean(NnbarEvents[i]) );
    myfile << "\n" << CutVariables[i] << " : " << Means[i] << "+/-" << StandardDevs[i]; // Output values to the file
  }

  // Message to output file
  myfile << "\n\n||Analyser.cc|| The Following " << NumNeutFiles << " Neutrino Files were analysed:\n";
  
  std::cout << "\nNow adding Neutrino files to chain and merginging into allneut.root file.\n";

  // Loop over Neutrino files and add them to TChain NeutTrees
  for(unsigned int i(0); i < NumNeutFiles; i ++){
    NeutTrees.Add(NeutFileNames[i].c_str()); 
    myfile << NeutFileNames[i] << "\n"; // Print File names to outputfile
  }  
  NeutTrees.Merge("allneut.root"); // MAY NOT NEED THIS LINE ANYMORE
  
  // Use ReadTree to return vector of vectors that contains all the information for all the events for neutrinos
  // Only grandmaster luke understands the indexing
  vector<vector<float>> NeutEvents = ReadTree("allneut.root",true); 

  // Calculate optimal sigma+/-
//  vector<vector<float>> sigPlusMinus = OptimisePlusMinus(NeutEvents, NnbarEvents, Means, Sigma, 5);
  // Hardcoded for speed, uncomment above and remove these lines to do it properlt
  /*vector<float> sigPlus;
  sigPlus.push_back(2.2);sigPlus.push_back(2.4);sigPlus.push_back(2.4);sigPlus.push_back(2.4);sigPlus.push_back(2.2);sigPlus.push_back(2.2);sigPlus.push_back(2.2);sigPlus.push_back(2.2);sigPlus.push_back(2.2);sigPlus.push_back(2.2);sigPlus.push_back(2.2);sigPlus.push_back(2.4);sigPlus.push_back(2.4);sigPlus.push_back(2.2);sigPlus.push_back(2.4);sigPlus.push_back(2.4);sigPlus.push_back(2.2);sigPlus.push_back(2.4);sigPlus.push_back(2.2);sigPlus.push_back(2.2);
  vector<float> sigMinus;
  sigMinus.push_back(1.8);sigMinus.push_back(1.6);sigMinus.push_back(1.6);sigMinus.push_back(1.6);sigMinus.push_back(1.8);sigMinus.push_back(1.8);sigMinus.push_back(1.8);sigMinus.push_back(1.8);sigMinus.push_back(1.6);sigMinus.push_back(1.6);sigMinus.push_back(1.6);sigMinus.push_back(1.6);sigMinus.push_back(1.6);sigMinus.push_back(1.8);sigMinus.push_back(1.6);sigMinus.push_back(1.6);sigMinus.push_back(1.8);sigMinus.push_back(1.8);sigMinus.push_back(1.8);sigMinus.push_back(1.8);
  sigPlusMinus.push_back(sigPlus); sigPlusMinus.push_back(sigMinus);*/
  // Work out the ratios for the second stage of optimisation
/*  vector<float> PlusMinusRatio;// = {1.2222, 1.5, 1.5, 1.5, 1.2222, 1.2222, 1.2222, 1.2222, 1.375, 1.375, 1.375, 1.5, 1.5, 1.2222, 1.5, 1.5, 1.2222, 1.3333, 1.2222, 1.2222};
  for(unsigned int k(0); k<5; k ++){
    for(unsigned int n(0); n < 4; n ++){
      PlusMinusRatio.push_back(sigPlusMinus[0][k*4+n]/sigPlusMinus[1][k*4+n]);
    }
  }
  // Optimise for each plane combination
  vector<vector<float>> sigUV = OptimiseTwoPlanes(NeutEvents, NnbarEvents, 1, 2, Means, Sigma, PlusMinusRatio, 5);
  vector<vector<float>> sigUY = OptimiseTwoPlanes(NeutEvents, NnbarEvents, 1, 3, Means, Sigma, PlusMinusRatio, 5);
  vector<vector<float>> sigVY = OptimiseTwoPlanes(NeutEvents, NnbarEvents, 2, 3, Means, Sigma, PlusMinusRatio, 5);

  for(unsigned int k(0); k<5; k ++){
    std::cout<<sigUV[0][k]<<"  "<<sigUV[1][k]<<std::endl;
    std::cout<<sigUY[0][k]<<"  "<<sigUY[1][k]<<std::endl;
    std::cout<<sigVY[0][k]<<"  "<<sigVY[1][k]<<std::endl;
  }

  unsigned int NSTDs = 1;
  std::cout << "\nMerge Successful. Now calculating number of merged events.\nPlease enter the number of standard deviations you want to cut.\n";
  if( std::cin >> NSTDs ){ // Check input is valid
    //for(unsigned int i(0); i < StandardDevs.size(); i ++){ StandardDevP[i] *= sigPlusMinus[0][i]; StandardDevM[i] *= sigPlusMinus[1][i]; }
    // Calculate numbers of removed neutrino and nnbar events for the optimised standard deviations
    vector<vector<float>> RemovedEvents = SingleCut(NeutEvents,true,Means,StandardDevP,StandardDevM,StandardDevP,StandardDevM,sigPlusMinus,sigUV,sigUY,sigVY);
    vector<vector<float>> RemovedNnbar = SingleCut(NnbarEvents,false,Means,StandardDevP,StandardDevM,StandardDevP,StandardDevM,sigPlusMinus,sigUV,sigUY,sigVY);
  }else{ cerr << "\n\nError! You didn't enter a valid number!\n"; exit(1); }

//-------------------------------------------------------------------------------------------------------\\
  // OUTPUT SECTION
  std::cout << "\nNow performing plots!\n";
  float size = 100;
  myfile << "\n\nTotal number of events removed:";

  // Print messages to file about individual cut variables
  for(unsigned int i(0); i < 20; i ++){ 
    myfile << "\n\n" << CutVariables[i] << "\n"<<"Sigma + = "<<sigPlusMinus[0][i]<<", Sigma - = "<<sigPlusMinus[1][i]<<"\n";  
    PrintText(myfile,RemovedEvents[i],RemovedEvents[40],size);
    myfile << "\n" << "Nnbar Removed = " << RemovedNnbar[i][0]/size << "%";
  }

  myfile << "\n\nNumber of events removed for combined planes:";
  
  // Print messages to file about combined plane cut variables 
  for(unsigned int i(0); i < 5; i ++){ 
    myfile << "\n\n" << "Cut Variable is " << CutVariables[i*4];
    myfile << "\n" << "For UV Plane, \n"<<"Sigma U = "<<sigUV[0][i]<<", Sigma V = "<<sigUV[1][i]<<"\n";  
    PrintText(myfile,RemovedEvents[20+i],RemovedEvents[40],size);  
    myfile << "\n" << "Nnbar Removed = " << RemovedNnbar[20+i][0]/size << "%";
    myfile << "\n" << "For UY Plane, \n"<<"Sigma U = "<<sigUY[0][i]<<", Sigma Y = "<<sigUY[1][i]<<"\n";  
    PrintText(myfile,RemovedEvents[25+i],RemovedEvents[40],size);
    myfile << "\n" << "Nnbar Removed = " << RemovedNnbar[25+i][0]/size << "%";
    myfile << "\n" << "For VY Plane, \n"<<"Sigma V = "<<sigVY[0][i]<<", Sigma Y = "<<sigVY[1][i]<<"\n";  
    PrintText(myfile,RemovedEvents[30+i],RemovedEvents[40],size);
    myfile << "\n" << "Nnbar Removed = " << RemovedNnbar[30+i][0]/size << "%";
    myfile << "\n" << "For UVY Plane, \n";  
    PrintText(myfile,RemovedEvents[35+i],RemovedEvents[40],size);
    myfile << "\n" << "Nnbar Removed = " << RemovedNnbar[35+i][0]/size << "%";
  }*/
/*
  myfile << "\n\nCombined cuts:";

  for(unsigned int l(0); l < 5; l ++){
    for(unsigned int m(l+1); m < 5; m ++){
      vector<vector<float>> ComboRemoved = CombineTwoCuts(NeutEvents,true,Means,StandardDevs,StandardDevs,l*4,m*4);
      vector<vector<float>> ComboRemovedNnbar = CombineTwoCuts(NnbarEvents,false,Means,StandardDevs,StandardDevs,l*4,m*4);
      myfile << "\n\n" << CutVariables[l*4] << " + " << CutVariables[m*4];
      myfile << "\n" << "Total, \n";  PrintText(myfile,ComboRemoved[0],RemovedEvents[40],size);  
      myfile << "\n" << "Nnbar Removed = " << ComboRemovedNnbar[0][0]/size << "%";
      myfile << "\n" << "For U Plane, \n";  PrintText(myfile,ComboRemoved[1],RemovedEvents[40],size);  
      myfile << "\n" << "Nnbar Removed = " << ComboRemovedNnbar[1][0]/size << "%";
      myfile << "\n" << "For V Plane, \n";  PrintText(myfile,ComboRemoved[2],RemovedEvents[40],size);  
      myfile << "\n" << "Nnbar Removed = " << ComboRemovedNnbar[2][0]/size << "%";
      myfile << "\n" << "For Y Plane, \n";  PrintText(myfile,ComboRemoved[3],RemovedEvents[40],size);  
      myfile << "\n" << "Nnbar Removed = " << ComboRemovedNnbar[3][0]/size << "%";
      myfile << "\n" << "For UV Plane, \n";  PrintText(myfile,ComboRemoved[4],RemovedEvents[40],size);  
      myfile << "\n" << "Nnbar Removed = " << ComboRemovedNnbar[4][0]/size << "%";
      myfile << "\n" << "For UY Plane, \n";  PrintText(myfile,ComboRemoved[5],RemovedEvents[40],size);
      myfile << "\n" << "Nnbar Removed = " << ComboRemovedNnbar[5][0]/size << "%";
      myfile << "\n" << "For VY Plane, \n";  PrintText(myfile,ComboRemoved[6],RemovedEvents[40],size);
      myfile << "\n" << "Nnbar Removed = " << ComboRemovedNnbar[6][0]/size << "%";
      myfile << "\n" << "For UVY Plane, \n";  PrintText(myfile,ComboRemoved[7],RemovedEvents[40],size);
      myfile << "\n" << "Nnbar Removed = " << ComboRemovedNnbar[7][0]/size << "%";
      ComboRemoved.clear();
      ComboRemovedNnbar.clear();
    }
  }
*/
  myfile.close();  
/*
  // Surface plot of combined cuts
  vector <float> Sigma1, Sigma2;
  float limit = 5;
  //vector <float> s1, s2, numiRem;
  //vector <bool> pass;
  TCanvas *c11 = new TCanvas("c11","Graph2D example",0,0,600,400);
  TGraph2D *dt = new TGraph2D();
  TGraph2D *dtn = new TGraph2D();
  TGraph2D *dt1 = new TGraph2D();
  TGraph2D *dtn1 = new TGraph2D();
  TGraph2D *dt2 = new TGraph2D();
  TGraph2D *dtn2 = new TGraph2D();
  TGraph2D *dt3 = new TGraph2D();
  TGraph2D *dtn3 = new TGraph2D();
  TGraph2D *dt4 = new TGraph2D();
  TGraph2D *dtn4 = new TGraph2D();
  TGraph2D *dt5 = new TGraph2D();
  TGraph2D *dtn5 = new TGraph2D();
  TGraph2D *dt6 = new TGraph2D();
  TGraph2D *dtn6 = new TGraph2D();
  TGraph2D *dt7 = new TGraph2D();
  TGraph2D *dtn7 = new TGraph2D();
  float maxremt(0), maxremu(0), maxremv(0), maxremy(0), maxremuv(0), maxremuy(0), maxremvy(0), maxremuvy(0);
  int lim = 20;
  // Fill the sigmas with zeros
  for(unsigned int i(0); i < Sigma.size(); i ++){ Sigma1.push_back(0.0); Sigma2.push_back(0.0); }
  for(unsigned int l(0); l < lim; l ++){
    for(unsigned int m(0); m < lim; m ++){
      for(unsigned int i(0); i < Sigma.size(); i ++){ Sigma1[i] = Sigma[i]*(Double_t)(l)/5; Sigma2[i] = Sigma[i]*(Double_t)(m)/5; }
      //s1.push_back(l/2);
      //s2.push_back(m/2);
      vector<vector<float>> ComboRemoved = CombineTwoCuts(NeutEvents,true,Means,Sigma1,Sigma2,0,16);
      vector<vector<float>> ComboRemovedNnbar = CombineTwoCuts(NnbarEvents,false,Means,Sigma1,Sigma2,0,16);
      //numiRem.push_back(ComboRemoved[0][0]);
      //if(ComboRemovedNnbar[0][0]<limit){dt->SetPoint((l*8+m),(Double_t)(l)/2,(Double_t)(m)/2,(Double_t)ComboRemoved[0][0]);}
      //else{dt->SetPoint(((l*8+m)-36),(Double_t)(l)/2,(Double_t)(m)/2,0.0);}
      dt->SetPoint(((l-0)*lim+(m-0)),(Double_t)(l)/5,(Double_t)(m)/5,(Double_t)(ComboRemoved[0][0]/100));
      dtn->SetPoint(((l-0)*lim+(m-0)),(Double_t)(l)/5,(Double_t)(m)/5,(Double_t)(ComboRemovedNnbar[0][0]/100));
      dt1->SetPoint(((l-0)*lim+(m-0)),(Double_t)(l)/5,(Double_t)(m)/5,(Double_t)(ComboRemoved[1][0]/100));
      dtn1->SetPoint(((l-0)*lim+(m-0)),(Double_t)(l)/5,(Double_t)(m)/5,(Double_t)(ComboRemovedNnbar[1][0]/100));
      dt2->SetPoint(((l-0)*lim+(m-0)),(Double_t)(l)/5,(Double_t)(m)/5,(Double_t)(ComboRemoved[2][0]/100));
      dtn2->SetPoint(((l-0)*lim+(m-0)),(Double_t)(l)/5,(Double_t)(m)/5,(Double_t)(ComboRemovedNnbar[2][0]/100));
      dt3->SetPoint(((l-0)*lim+(m-0)),(Double_t)(l)/5,(Double_t)(m)/5,(Double_t)(ComboRemoved[3][0]/100));
      dtn3->SetPoint(((l-0)*lim+(m-0)),(Double_t)(l)/5,(Double_t)(m)/5,(Double_t)(ComboRemovedNnbar[3][0]/100));
      dt4->SetPoint(((l-0)*lim+(m-0)),(Double_t)(l)/5,(Double_t)(m)/5,(Double_t)(ComboRemoved[4][0]/100));
      dtn4->SetPoint(((l-0)*lim+(m-0)),(Double_t)(l)/5,(Double_t)(m)/5,(Double_t)(ComboRemovedNnbar[4][0]/100));
      dt5->SetPoint(((l-0)*lim+(m-0)),(Double_t)(l)/5,(Double_t)(m)/5,(Double_t)(ComboRemoved[5][0]/100));
      dtn5->SetPoint(((l-0)*lim+(m-0)),(Double_t)(l)/5,(Double_t)(m)/5,(Double_t)(ComboRemovedNnbar[5][0]/100));
      dt6->SetPoint(((l-0)*lim+(m-0)),(Double_t)(l)/5,(Double_t)(m)/5,(Double_t)(ComboRemoved[6][0]/100));
      dtn6->SetPoint(((l-0)*lim+(m-0)),(Double_t)(l)/5,(Double_t)(m)/5,(Double_t)(ComboRemovedNnbar[6][0]/100));
      dt7->SetPoint(((l-0)*lim+(m-0)),(Double_t)(l)/5,(Double_t)(m)/5,(Double_t)(ComboRemoved[7][0]/100));
      dtn7->SetPoint(((l-0)*lim+(m-0)),(Double_t)(l)/5,(Double_t)(m)/5,(Double_t)(ComboRemovedNnbar[7][0]/100));
      if (ComboRemovedNnbar[0][0]/100 < limit && ComboRemoved[0][0]/100 > maxremt){maxremt = ComboRemoved[0][0]/100;}
      if (ComboRemovedNnbar[1][0]/100 < limit && ComboRemoved[1][0]/100 > maxremu){maxremu = ComboRemoved[1][0]/100;}
      if (ComboRemovedNnbar[2][0]/100 < limit && ComboRemoved[2][0]/100 > maxremv){maxremv = ComboRemoved[2][0]/100;}
      if (ComboRemovedNnbar[3][0]/100 < limit && ComboRemoved[3][0]/100 > maxremy){maxremy = ComboRemoved[3][0]/100;}
      if (ComboRemovedNnbar[4][0]/100 < limit && ComboRemoved[4][0]/100 > maxremuv){maxremuv = ComboRemoved[4][0]/100;}
      if (ComboRemovedNnbar[5][0]/100 < limit && ComboRemoved[5][0]/100 > maxremuy){maxremuy = ComboRemoved[5][0]/100;}
      if (ComboRemovedNnbar[6][0]/100 < limit && ComboRemoved[6][0]/100 > maxremvy){maxremvy = ComboRemoved[6][0]/100;}
      if (ComboRemovedNnbar[7][0]/100 < limit && ComboRemoved[7][0]/100 > maxremuvy){maxremuvy = ComboRemoved[7][0]/100;}
      std::cout<<(l-0)*lim+(m-0)<<"  "<<(Double_t)(l)/5<<"  "<<(Double_t)(m)/5<<"  "<<(Double_t)(ComboRemoved[0][0]/100)<<"  "<<(ComboRemovedNnbar[0][0]/100)<<std::endl;
    }
  }
  std::cout<<"t: "<<maxremt<<" u: "<<maxremu<<" v: "<<maxremv<<" y: "<<maxremy<<" uv: "<<maxremuv<<" uy: "<<maxremuy<<" vy: "<<maxremvy<<" uvy: "<<maxremuvy;
   gStyle->SetPalette(1);
   dtn->Draw("surf");
   dt->Draw("surf1same");
   TCanvas *c12 = new TCanvas("c12","Graph2D example",0,0,600,400);
   gStyle->SetPalette(1);
   dtn1->Draw("surf");
   dt1->Draw("surf1same");
   TCanvas *c13 = new TCanvas("c13","Graph2D example",0,0,600,400);
   gStyle->SetPalette(1);
   dtn2->Draw("surf");
   dt2->Draw("surf1same");
   TCanvas *c14 = new TCanvas("c14","Graph2D example",0,0,600,400);
   gStyle->SetPalette(1);
   dtn3->Draw("surf");
   dt3->Draw("surf1same");
   TCanvas *c15 = new TCanvas("c15","Graph2D example",0,0,600,400);
   gStyle->SetPalette(1);
   dtn4->Draw("surf");
   dt4->Draw("surf1same");
   TCanvas *c16 = new TCanvas("c16","Graph2D example",0,0,600,400);
   gStyle->SetPalette(1);
   dtn5->Draw("surf");
   dt5->Draw("surf1same");
   TCanvas *c17 = new TCanvas("c17","Graph2D example",0,0,600,400);
   gStyle->SetPalette(1);
   dtn6->Draw("surf");
   dt6->Draw("surf1same");
   TCanvas *c18 = new TCanvas("c18","Graph2D example",0,0,600,400);
   gStyle->SetPalette(1);
   dtn7->Draw("surf");
   dt7->Draw("surf1same");
*/
/*
  // Plot each cut variable for nnbar and Neutrino on the same plot.  
  for(unsigned int i(0); i < 20; i ++ ){
    TCanvas *c1 = new TCanvas(("c" + CutVariables[i]).c_str(),CutVariables[i].c_str(),900,600);

    DoubleHist(nnbarTrees,NeutTrees,CutVariables[i],AxisLabels[i]);
    c1->Update();
    TLine *MinSig = new TLine(Means[i]-StandardDevM[i]*sigPlusMinus[1][i],0,Means[i]-StandardDevM[i]*sigPlusMinus[1][i],c1->GetUymax());
    TLine *MaxSig = new TLine(Means[i]+StandardDevP[i]*sigPlusMinus[0][i],0,Means[i]+StandardDevP[i]*sigPlusMinus[0][i],c1->GetUymax());
    MinSig->Draw(); MaxSig->Draw();

    c1->Update();
    c1->Modified(); 
  }
  
  // Plot 2D histograms for cuts on UV, VY and UY planes for each cut variable
  for(unsigned int i(1); i < 20; i = i + 4){

    float sigUP = StandardDevs[i]*(sigPlusMinus[0][i]/sigPlusMinus[1][i])*sigUV[0][(i-1)/4];
    float sigUM = StandardDevs[i]*sigUV[0][(i-1)/4];
    float sigVP = StandardDevs[i+1]*(sigPlusMinus[0][i+1]/sigPlusMinus[1][i+1])*sigUV[1][(i-1)/4];
    float sigVM = StandardDevs[i+1]*sigUV[1][(i-1)/4];
    TCanvas *c2 = new TCanvas(CutVariables[i].c_str(),"The Title Of This Canvas",900,600);
    twoDoubleHist(nnbarTrees,NeutTrees,CutVariables[i],CutVariables[i+1],AxisLabels[i],AxisLabels[i+1]);
    TLine *XMinSig1 = new TLine(Means[i+1]-sigVM,Means[i]-sigUM,Means[i+1]-sigVM,Means[i]+sigUP);
    TLine *XMaxSig1 = new TLine(Means[i+1]+sigVP,Means[i]-sigUM,Means[i+1]+sigVP,Means[i]+sigUP);
    TLine *YMinSig1 = new TLine(Means[i+1]-sigVM,Means[i]-sigUM,Means[i+1]+sigVP,Means[i]-sigUM);
    TLine *YMaxSig1 = new TLine(Means[i+1]-sigVM,Means[i]+sigUP,Means[i+1]+sigVP,Means[i]+sigUP);
    XMinSig1->Draw(); XMaxSig1->Draw(); 
    YMinSig1->Draw(); YMaxSig1->Draw();

    c2->Update(); c2->Modified();
    
    float sigYP = StandardDevs[i+2]*(sigPlusMinus[0][i+2]/sigPlusMinus[1][i+2])*sigVY[1][(i-1)/4];
    float sigYM = StandardDevs[i+2]*sigVY[1][(i-1)/4];
    float sigVP1 = StandardDevs[i+1]*(sigPlusMinus[0][i+1]/sigPlusMinus[1][i+1])*sigVY[0][(i-1)/4];
    float sigVM1 = StandardDevs[i+1]*sigVY[0][(i-1)/4];
    TCanvas *c3 = new TCanvas(CutVariables[(i+1)].c_str(),"The Title Of This Canvas",900,600);
    twoDoubleHist(nnbarTrees,NeutTrees,CutVariables[i+1],CutVariables[i+2],AxisLabels[i+1],AxisLabels[i+2]);
    TLine *XMinSig2 = new TLine(Means[i+2]-sigYM,Means[i+1]-sigVM1,Means[i+2]-sigYM,Means[i+1]+sigVP1);
    TLine *XMaxSig2 = new TLine(Means[i+2]+sigYP,Means[i+1]-sigVM1,Means[i+2]+sigYP,Means[i+1]+sigVP1);
    TLine *YMinSig2 = new TLine(Means[i+2]-sigYM,Means[i+1]-sigVM1,Means[i+2]+sigYP,Means[i+1]-sigVM1);
    TLine *YMaxSig2 = new TLine(Means[i+2]-sigYM,Means[i+1]+sigVP1,Means[i+2]+sigYP,Means[i+1]+sigVP1);
    XMinSig2->Draw(); XMaxSig2->Draw();
    YMinSig2->Draw(); YMaxSig2->Draw();
    c3->Update(); c3->Modified();
   
    float sigYP1 = StandardDevs[i+2]*(sigPlusMinus[0][i+2]/sigPlusMinus[1][i+2])*sigUY[1][(i-1)/4];
    float sigYM1 = StandardDevs[i+2]*sigUY[1][(i-1)/4];
    float sigUP1 = StandardDevs[i]*(sigPlusMinus[0][i]/sigPlusMinus[1][i])*sigUY[0][(i-1)/4];
    float sigUM1 = StandardDevs[i]*sigUY[0][(i-1)/4];
    TCanvas *c4 = new TCanvas(CutVariables[(i+2)].c_str(),"The Title Of This Canvas",900,600);
    twoDoubleHist(nnbarTrees,NeutTrees,CutVariables[i],CutVariables[i+2],AxisLabels[i],AxisLabels[i+2]);
    TLine *XMinSig3 = new TLine(Means[i+2]-sigYM1,Means[i]-sigUM1,Means[i+2]-sigYM1,Means[i]+sigUP1);
    TLine *XMaxSig3 = new TLine(Means[i+2]+sigYP1,Means[i]-sigUM1,Means[i+2]+sigYP1,Means[i]+sigUP1);
    TLine *YMinSig3 = new TLine(Means[i+2]-sigYM1,Means[i]-sigUM1,Means[i+2]+sigYP1,Means[i]-sigUM1);
    TLine *YMaxSig3 = new TLine(Means[i+2]-sigYM1,Means[i]+sigUP1,Means[i+2]+sigYP1,Means[i]+sigUP1);
    XMinSig3->Draw(); XMaxSig3->Draw();
    YMinSig3->Draw(); YMaxSig3->Draw();
    c4->Update(); c4->Modified();
  }

*/
  // Plot number of Each Type on a stacked plot.
  DoubleHistType(NeutTrees);


  std::cout << "\nHello World!\n";
  StandardDevs.clear(); Means.clear(); NnbarEvents.clear();RemovedEvents.clear(); NeutEvents.clear();
}
