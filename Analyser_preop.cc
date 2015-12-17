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
  string LabelNames[] = {"CCQE (38.3%)","NCQE (18.2%)","CCRE (20.7%)","NCRE (8.6%)","CCDIS (9.7%)","NCDIS (3.9%)","CCCO (0.4%)","NCCO (0.2%)"};
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
  hs->SetMaximum(2000);// CURRENTLY THIS IS SET MANUALLY!
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

// Function designed to give me a headache
vector < vector <float> > GetChainSTD(string RootFileName, bool isNeut, bool cut, vector <float> MEANS, vector <float> &STDEV){

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

  if( cut == true ){ // For Neutrino Events
    // Count number of removed events outside of some tolerance of hit wires
    // Total cut
    vector < vector < float > > RemovedType;
    vector < float > EmptyVect; 
    for(unsigned int i(0); i < 9; i ++) EmptyVect.push_back(0);

    for(unsigned int i(0); i < 20; i ++){
      RemovedType.push_back(EmptyVect); // Removed Type After this loop will be of dimensions RemovedType[20][8]
      //std::cout << "\n\nIn loop " << i << "\nMeans[" << i << "] is: " << MEANS[i] << "\nSTDEV[" << i << "] is: " << STDEV[i];

      for(unsigned int j(0); j < vv_ReturnVect[0].size(); j ++){ // LOOP NUMBER DICTATED BY v_Hits.size()

        //std::cout << "\n\nLooping over Neutrino Data vv_ReturnVect[" << j << "].\nvv_ReturnVect[" << i << "]["<<j<<"] is: " << vv_ReturnVect[i][j];
        if( vv_ReturnVect[i][j] > (MEANS[i] + STDEV[i]) || vv_ReturnVect[i][j] < (MEANS[i] - STDEV[i]) ){ 
          RemovedType[i][0] ++;
	  if(isNeut == true){
          for(unsigned int k=0; k < 8; k ++){
            if(v_Type[j]==k) RemovedType[i][k+1] ++;
          }
	  }
        }
        //for(unsigned int z(0); z < 9; z++){ std::cout << "\nRemovedType["<<i<<"]["<<z<<"] = " << RemovedType[i][z]; }
      }
    }
    //std::cout << "\n\nNow Performing Cuts for Multiple planes!\n";
    for(unsigned int j(0); j < 20; j++) RemovedType.push_back(EmptyVect); // Removed Type After this loop will be of dimensions RemovedType[40][8]
    for(unsigned int p(0); p < 5; p ++){
      // Index == (1+p*4) takes values 1,5,9,13,17. These correspond to vector elements v_UHits, v_TDCstdU, v_TDCiqrU, v_ADCampU and v_WFintU
      // Index == (2+p*4) takes values 2,6,10,14,18. These correspond to vector elements v_VHits, v_TDCstdV, v_TDCiqrV, v_ADCampV and v_WFintV
      // Index == (3+p*4) takes values 3,7,11,15,19. These correspond to vector elements v_YHits, v_TDCstdY, v_TDCiqrY, v_ADCampY and v_WFintY
      unsigned int UPI = (1+p*4); // UPlaneIndex
      unsigned int VPI = (2+p*4); // VPlaneIndex
      unsigned int YPI = (3+p*4); // YPlaneIndex
      /*std::cout << "\n\nIn loop " << p << ".";
      std::cout << "\nUMeans[" << UPI << "] is: " << MEANS[UPI] << "\nUSTDEV[" << UPI << "] is: " << STDEV[UPI];
      std::cout << "\nVMeans[" << VPI << "] is: " << MEANS[VPI] << "\nVSTDEV[" << VPI << "] is: " << STDEV[VPI];
      std::cout << "\nYMeans[" << YPI << "] is: " << MEANS[YPI] << "\nYSTDEV[" << YPI << "] is: " << STDEV[YPI];
      std::cout << "\n\nLooping over UPlane, VPLANE and YPLANE Neutrino Data Vectors [" << UPI << ", " << VPI << ", " << YPI << "].";*/
      for(unsigned int j(0); j < vv_ReturnVect[0].size(); j ++){ // LOOP NUMBER DICTATED BY v_Hits.size()

        /*std::cout << "\nvv_ReturnVect[" << UPI << "]["<<j<<"] is: " << vv_ReturnVect[UPI][j];
        std::cout << "\nvv_ReturnVect[" << VPI << "]["<<j<<"] is: " << vv_ReturnVect[VPI][j];
        std::cout << "\nvv_ReturnVect[" << YPI << "]["<<j<<"] is: " << vv_ReturnVect[YPI][j];*/

        // UV Plane Cut
        if(
          (vv_ReturnVect[UPI][j]>(MEANS[UPI] + STDEV[UPI])||vv_ReturnVect[UPI][j]<(MEANS[UPI] - STDEV[UPI])) ||
          (vv_ReturnVect[VPI][j]>(MEANS[VPI] + STDEV[VPI])||vv_ReturnVect[VPI][j]<(MEANS[VPI] - STDEV[VPI]))){
          RemovedType[20+p][0] ++; // A Neutrino was removed from UV plane cut
	  if(isNeut == true){
          for(unsigned int k=0; k<8;k++){
            if(v_Type[j]==k){RemovedType[20+p][k+1] ++;} // Find flavour of neutrino
          }
	  }
        }
      
        // UY Plane Cut
        if(
          (vv_ReturnVect[UPI][j]>(MEANS[UPI] + STDEV[UPI])||vv_ReturnVect[UPI][j]<(MEANS[UPI] - STDEV[UPI])) ||
          (vv_ReturnVect[YPI][j]>(MEANS[YPI] + STDEV[YPI])||vv_ReturnVect[YPI][j]<(MEANS[YPI] - STDEV[YPI]))){ 
          RemovedType[25+p][0] ++; // A Neutrino was removed from UY plane cut
	  if(isNeut == true){
	  for(unsigned int k=0; k<8;k++){
	    if(v_Type[j]==k){RemovedType[25+p][k+1] ++;} // Find flavour of neutrino
          }
	  }
        }
    
        // VY Plane Cut
        if(
          (vv_ReturnVect[VPI][j]>(MEANS[VPI] + STDEV[VPI])||vv_ReturnVect[VPI][j]<(MEANS[VPI] - STDEV[VPI])) ||
          (vv_ReturnVect[YPI][j]>(MEANS[YPI] + STDEV[YPI])||vv_ReturnVect[YPI][j]<(MEANS[YPI] - STDEV[YPI]))){
          RemovedType[30+p][0] ++; // A Neutrino was removed from VY plane cut
	  if(isNeut == true){
          for(unsigned int k=0; k<8;k++){
	    if(v_Type[j]==k){RemovedType[30+p][k+1] ++;} // Find flavour of neutrino
          }
	  }
        }
     
        // UVY Plane Cut
        if(
          (vv_ReturnVect[UPI][j]>(MEANS[UPI] + STDEV[UPI])||vv_ReturnVect[UPI][j]<(MEANS[UPI] - STDEV[UPI])) ||
          (vv_ReturnVect[VPI][j]>(MEANS[VPI] + STDEV[VPI])||vv_ReturnVect[VPI][j]<(MEANS[VPI] - STDEV[VPI])) ||
          (vv_ReturnVect[YPI][j]>(MEANS[YPI] + STDEV[YPI])||vv_ReturnVect[YPI][j]<(MEANS[YPI] - STDEV[YPI]))){
          RemovedType[35+p][0] ++; // A Neutrino was removed from UVY plane cut
	  if(isNeut == true){
          for(unsigned int k=0; k<8;k++){
            if(v_Type[j]==k){RemovedType[35+p][k+1] ++;} // Find flavour of neutrino
          }
	  }
        }
      }
           /*
      for(unsigned int z(0); z < 9; z++){ 
        std::cout << "\nUVPlane Cut RemovedType["<<(20+p)<<"]["<<z<<"] = " << RemovedType[20+p][z]; 
        std::cout << "\nUYPlane Cut RemovedType["<<(25+p)<<"]["<<z<<"] = " << RemovedType[25+p][z];
        std::cout << "\nVYPlane Cut RemovedType["<<(30+p)<<"]["<<z<<"] = " << RemovedType[30+p][z];
        std::cout << "\nUVYPlane Cut RemovedType["<<(35+p)<<"]["<<z<<"] = " << RemovedType[35+p][z] << "\n";
      }*/
    }
    if(isNeut == true){
    for(unsigned int i=0; i < v_Type.size(); i++){
      for(unsigned int j=0;j<8;j++){
        if(v_Type[i] == j) EmptyVect[j] ++; // Loop over all neutrino types. See vector initialisation above.
      }
    }
    RemovedType.push_back(EmptyVect);
    }
    return RemovedType;
  } // End of if( cut == true )

  return vv_ReturnVect;
}

// FUnction to combine two cuts together
vector < vector <float> > CombineTwoCuts(string RootFileName, bool isNeut, const vector <float> & MEAN, const vector <float> & STDEV1, 
                                         const vector <float> & STDEV2, int cut1, int cut2){//cut=0,4,8,12,16 (hn,ts,ti,aa,wi)

  // Initilise float variables, vector of floats and vector of vector of floats
  float TDCstd, TDCstdU,TDCstdV,TDCstdY,TDCiqr,TDCiqrU,TDCiqrV,TDCiqrY,ADCamp,ADCampU,ADCampV,ADCampY,WFint,WFintU,WFintV,WFintY,Energy;
  vector <float> v_TDCstd, v_TDCstdU,v_TDCstdV,v_TDCstdY,v_TDCiqr,v_TDCiqrU,v_TDCiqrV,v_TDCiqrY,v_ADCamp,v_ADCampU,v_ADCampV,v_ADCampY,v_WFint,v_WFintU,v_WFintV,v_WFintY,v_Energ,v_Hits,v_UHits,v_VHits,v_YHits,v_Type,v_Energy;
  vector < vector < float > > vv_ReturnVect;

  UInt_t Hits,UHits,VHits,YHits,Type;

  TFile *TRootFile = new TFile(RootFileName.c_str());  
  TTree *FileTree = (TTree*)TRootFile->Get("ch_tree");
//delete TRootFile; 
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

    // Count number of removed events outside of some tolerance of hit wires
    // Total cut
    vector < vector < float > > RemovedType;
    vector < float > EmptyVect; 
    for(unsigned int i(0); i < 9; i ++) EmptyVect.push_back(0);

    for(unsigned int i(0); i < 4; i ++){
      RemovedType.push_back(EmptyVect);

      for(unsigned int j(0); j < vv_ReturnVect[0].size(); j ++){ // LOOP NUMBER DICTATED BY v_Hits.size()

        if( vv_ReturnVect[i+cut1][j] > (MEAN[i+cut1] + STDEV1[i+cut1]) || vv_ReturnVect[i+cut1][j] < (MEAN[i+cut1] - STDEV1[i+cut1]) || 
            vv_ReturnVect[i+cut2][j] > (MEAN[i+cut2] + STDEV2[i+cut2]) || vv_ReturnVect[i+cut2][j] < (MEAN[i+cut2] - STDEV2[i+cut2]) ){ 
          RemovedType[i][0] ++;
	  if(isNeut == true){
          for(unsigned int k=0; k < 8; k ++){
            if(v_Type[j]==k) RemovedType[i][k+1] ++;
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
	if(v_Type[j]==k){RemovedType[4][k+1] ++;} // Find flavour of neutrino
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
	if(v_Type[j]==k){RemovedType[5][k+1] ++;} // Find flavour of neutrino
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
	if(v_Type[j]==k){RemovedType[6][k+1] ++;} // Find flavour of neutrino
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
	if(v_Type[j]==k){RemovedType[7][k+1] ++;} // Find flavour of neutrino
      }
      }
    }
    }
    return RemovedType;
}

// Main Script
void Analyser_preop(){
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
    nnbarFileNames[i] = "outputs5_2/Nnbar_output_run5_" + FileNumAsString + ".root";
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
   
  vector <float> StandardDevs, Means, Sigma;
  // Call function GetChainSTD to get the Standard Deviations and Means from the merged allnnbar.root
  // containing all nnbar data. Feed in dummy empty vector "Means" which is not used as second argument is false
  vector < vector <float> > ChainSTDs = GetChainSTD("allnnbar.root",false,false,Means,Means); 
  
  // Output message to file
  myfile << "\n\nMean Values and standard Deviations of Merged nnbarfile for different cut variables\n";
  myfile << "\n" << "Number of nnbar events: "<<ChainSTDs[0].size();

  std::cout << "\nData successfully read. Now calculating Standard Deviations and Means.\n";

  // Loop over data returned in ChainSTDs given by GetChainSTD() to get the Means and StandardDeviations
  for(unsigned int i(0); i < ChainSTDs.size(); i ++){ 
    StandardDevs.push_back( CalcStandardDev(ChainSTDs[i]) ); 
    Sigma.push_back( CalcStandardDev(ChainSTDs[i]) );
    Means.push_back( CalcMean(ChainSTDs[i]) );
    myfile << "\n" << CutVariables[i] << " : " << Means[i] << "+/-" << StandardDevs[i]; // Output values to the file
  }

  // Message to output file
  myfile << "\n\n||Analyser.cc|| The Folllowing " << NumNeutFiles << " Neutrino Files were analysed:\n";
  
  std::cout << "\nNow adding Neutrino files to chain and merginging into allneut.root file.\n";

  // Loop over Neutrino files and add them to TChain NeutTrees
  for(unsigned int i(0); i < NumNeutFiles; i ++){
    NeutTrees.Add(NeutFileNames[i].c_str()); 
    myfile << NeutFileNames[i] << "\n"; // Print File names to outputfile
  }  
  NeutTrees.Merge("allneut.root"); // MAY NOT NEED THIS LINE ANYMORE

  unsigned int NSTDs = 1;
  std::cout << "\nMerge Successful. Now calculating number of merged events.\nPlease enter the number of standard deviations you want to cut.\n";
  if( std::cin >> NSTDs ){ // Check input is valid
    for(unsigned int i(0); i < StandardDevs.size(); i ++){ StandardDevs[i] *= NSTDs; }// CHECK SAME OUTSIDE IF
    // Function GetChainSTD now returns the number of removed Neutrino Events using the Means and Standard Deviations given.
    vector < vector <float> > RemovedEvents = GetChainSTD("allneut.root",true,true,Means,StandardDevs);
    vector < vector <float> > RemovedNnbar = GetChainSTD("allnnbar.root",false,true,Means,StandardDevs);
  }else{ cerr << "\n\nError! You didn't enter a valid number!\n"; exit(1); }

//-------------------------------------------------------------------------------------------------------\\
  // OUTPUT SECTION
  std::cout << "\nNow performing plots!\n";
  float size = 100;
  myfile << "\n\nTotal number of events removed:";

  // Print messages to file about individual cut variables
  for(unsigned int i(0); i < 20; i ++){ 
    myfile << "\n\n" << CutVariables[i] << "\n";  PrintText(myfile,RemovedEvents[i],RemovedEvents[40],size);
    myfile << "\n" << "Nnbar Removed = " << RemovedNnbar[i][0]/size << "%";
  }

  myfile << "\n\nNumber of events removed for combined planes:";
  
  // Print messages to file about combined plane cut variables 
  for(unsigned int i(0); i < 5; i ++){ 
    myfile << "\n\n" << "Cut Variable is " << CutVariables[i*4];
    myfile << "\n" << "For UV Plane, \n";  PrintText(myfile,RemovedEvents[20+i],RemovedEvents[40],size);  
    myfile << "\n" << "Nnbar Removed = " << RemovedNnbar[20+i][0]/size << "%";
    myfile << "\n" << "For UY Plane, \n";  PrintText(myfile,RemovedEvents[25+i],RemovedEvents[40],size);
    myfile << "\n" << "Nnbar Removed = " << RemovedNnbar[25+i][0]/size << "%";
    myfile << "\n" << "For VY Plane, \n";  PrintText(myfile,RemovedEvents[30+i],RemovedEvents[40],size);
    myfile << "\n" << "Nnbar Removed = " << RemovedNnbar[30+i][0]/size << "%";
    myfile << "\n" << "For UVY Plane, \n";  PrintText(myfile,RemovedEvents[35+i],RemovedEvents[40],size);
    myfile << "\n" << "Nnbar Removed = " << RemovedNnbar[35+i][0]/size << "%";
  }

  myfile << "\n\nCombined cuts:";

  for(unsigned int l(0); l < 5; l ++){
    for(unsigned int m(l+1); m < 5; m ++){
      vector < vector <float> > ComboRemoved = CombineTwoCuts("allneut.root",true,Means,StandardDevs,StandardDevs,l*4,m*4);
      vector < vector <float> > ComboRemovedNnbar = CombineTwoCuts("allnnbar.root",false,Means,StandardDevs,StandardDevs,l*4,m*4);
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
      vector < vector <float> > ComboRemoved = CombineTwoCuts("allneut.root",true,Means,Sigma1,Sigma2,0,16);
      vector < vector <float> > ComboRemovedNnbar = CombineTwoCuts("allnnbar.root",false,Means,Sigma1,Sigma2,0,16);
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

  // Plot each cut variable for nnbar and Neutrino on the same plot.  
  for(unsigned int i(0); i < 20; i ++ ){
    TCanvas *c1 = new TCanvas(("c" + CutVariables[i]).c_str(),CutVariables[i].c_str(),900,600);

    DoubleHist(nnbarTrees,NeutTrees,CutVariables[i],AxisLabels[i]);
    c1->Update();
    TLine *MinSig = new TLine(Means[i]-StandardDevs[i],0,Means[i]-StandardDevs[i],c1->GetUymax());
    TLine *MaxSig = new TLine(Means[i]+StandardDevs[i],0,Means[i]+StandardDevs[i],c1->GetUymax());
    MinSig->Draw(); MaxSig->Draw();

    c1->Update();
    c1->Modified(); 
  }
  
  // Plot 2D histograms for cuts on UV, VY and UY planes for each cut variable
  for(unsigned int i(1); i < 20; i = i + 4){
    TCanvas *c2 = new TCanvas(CutVariables[i].c_str(),"The Title Of This Canvas",900,600);
    twoDoubleHist(nnbarTrees,NeutTrees,CutVariables[i],CutVariables[i+1],AxisLabels[i],AxisLabels[i+1]);
    TLine *XMinSig1 = new TLine(Means[i+1]-StandardDevs[i+1],Means[i]-StandardDevs[i],Means[i+1]-StandardDevs[i+1],Means[i]+StandardDevs[i]);
    TLine *XMaxSig1 = new TLine(Means[i+1]+StandardDevs[i+1],Means[i]-StandardDevs[i],Means[i+1]+StandardDevs[i+1],Means[i]+StandardDevs[i]);
    TLine *YMinSig1 = new TLine(Means[i+1]-StandardDevs[i+1],Means[i]-StandardDevs[i],Means[i+1]+StandardDevs[i+1],Means[i]-StandardDevs[i]);
    TLine *YMaxSig1 = new TLine(Means[i+1]-StandardDevs[i+1],Means[i]+StandardDevs[i],Means[i+1]+StandardDevs[i+1],Means[i]+StandardDevs[i]);
    XMinSig1->Draw(); XMaxSig1->Draw(); 
    YMinSig1->Draw(); YMaxSig1->Draw();

    c2->Update(); c2->Modified();
    
    TCanvas *c3 = new TCanvas(CutVariables[(i+1)].c_str(),"The Title Of This Canvas",900,600);
    twoDoubleHist(nnbarTrees,NeutTrees,CutVariables[i+1],CutVariables[i+2],AxisLabels[i+1],AxisLabels[i+2]);
    TLine *XMinSig2 = new TLine(Means[i+2]-StandardDevs[i+2],Means[i+1]-StandardDevs[i+1],Means[i+2]-StandardDevs[i+2],Means[i+1]+StandardDevs[i+1]);
    TLine *XMaxSig2 = new TLine(Means[i+2]+StandardDevs[i+2],Means[i+1]-StandardDevs[i+1],Means[i+2]+StandardDevs[i+2],Means[i+1]+StandardDevs[i+1]);
    TLine *YMinSig2 = new TLine(Means[i+2]-StandardDevs[i+2],Means[i+1]-StandardDevs[i+1],Means[i+2]+StandardDevs[i+2],Means[i+1]-StandardDevs[i+1]);
    TLine *YMaxSig2 = new TLine(Means[i+2]-StandardDevs[i+2],Means[i+1]+StandardDevs[i+1],Means[i+2]+StandardDevs[i+2],Means[i+1]+StandardDevs[i+1]);
    XMinSig2->Draw(); XMaxSig2->Draw();
    YMinSig2->Draw(); YMaxSig2->Draw();
    c3->Update(); c3->Modified();
   
    TCanvas *c4 = new TCanvas(CutVariables[(i+2)].c_str(),"The Title Of This Canvas",900,600);
    twoDoubleHist(nnbarTrees,NeutTrees,CutVariables[i],CutVariables[i+2],AxisLabels[i],AxisLabels[i+2]);
    TLine *XMinSig3 = new TLine(Means[i+2]-StandardDevs[i+2],Means[i]-StandardDevs[i],Means[i+2]-StandardDevs[i+2],Means[i]+StandardDevs[i]);
    TLine *XMaxSig3 = new TLine(Means[i+2]+StandardDevs[i+2],Means[i]-StandardDevs[i],Means[i+2]+StandardDevs[i+2],Means[i]+StandardDevs[i]);
    TLine *YMinSig3 = new TLine(Means[i+2]-StandardDevs[i+2],Means[i]-StandardDevs[i],Means[i+2]+StandardDevs[i+2],Means[i]-StandardDevs[i]);
    TLine *YMaxSig3 = new TLine(Means[i+2]-StandardDevs[i+2],Means[i]+StandardDevs[i],Means[i+2]+StandardDevs[i+2],Means[i]+StandardDevs[i]);
    XMinSig3->Draw(); XMaxSig3->Draw();
    YMinSig3->Draw(); YMaxSig3->Draw();
    c4->Update(); c4->Modified();
  }

  // Plot number of Each Type on a stacked plot.
  DoubleHistType(NeutTrees);


  std::cout << "\nHello World!\n";
  StandardDevs.clear(); Means.clear(); ChainSTDs.clear();RemovedEvents.clear();
}
