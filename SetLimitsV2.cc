#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <queue>


// Cut on vector of vector data vv_ReturnVect based on two Limits LIM1 and LIM2. Returns the number
// of cut neutrinos in the vector of Double_t's called RemovedType
vector <Double_t> DoTheCut( vector <vector <Double_t> > vv_ReturnVect, vector <Double_t> LIM1, vector <Double_t> LIM2, Size_t VectSize ){
  vector <Double_t> RemovedType;
  
  //std::cout << "\nvv_ReturnVect[0].size() is : " << vv_ReturnVect[0].size(); 
  for(unsigned int i(0); i < vv_ReturnVect.size(); i ++){ // Loop over number of vectors in vector of vectors
    RemovedType.push_back(0); // Removed Type After this loop will be of dimensions RemovedType[14] Typically
    //std::cout << "\n\nIn loop " << i << "\nMeans[" << i << "] is: " << MEANS[i];
    //std::cout << "\nSTDEVP[" << i << "] is: " << STDEVP[i] << "\nSTDEVM[" << i << "] is: " << STDEVM[i];
    for(unsigned int j(0); j < VectSize; j ++){ // LOOP NUMBER DICTATED BY v_Hits.size()
      //std::cout << "\n\nLooping over Neutrino Data vv_ReturnVect[" << j << "].\nvv_ReturnVect[" << i << "]["<<j<<"] is: " << vv_ReturnVect[i][j];
      if( LIM1[i] < vv_ReturnVect[i][j] && vv_ReturnVect[i][j] < LIM2[i] ){ // If data is within limits
        RemovedType[i] ++; // Add to RemovedType[i]
      }
    }
  }
  return RemovedType;
}

// Read normal tree data variables. Boolean determines whether the tree is neutrino type data containing additional variables "Type" and "Energy"
vector <vector<float>> ReadTree(string RootFileName, bool isNeut){

  // Initilise float variables, vector of floats and vector of vector of floats
  float TDCstd, TDCstdU,TDCstdV,TDCstdY,TDCiqr,TDCiqrU,TDCiqrV,TDCiqrY,ADCamp,ADCampU,ADCampV,ADCampY,WFint,WFintU,WFintV,WFintY,Energy;
  vector <float> v_TDCstd, v_TDCstdU,v_TDCstdV,v_TDCstdY,v_TDCiqr,v_TDCiqrU,v_TDCiqrV,v_TDCiqrY,v_ADCamp,v_ADCampU,v_ADCampV,v_ADCampY,v_WFint,v_WFintU,v_WFintV,v_WFintY,v_Energ,v_Hits,v_UHits,v_VHits,v_YHits,v_Type,v_Energy;
  vector < vector < float > > vv_ReturnVect;
  UInt_t Hits,UHits,VHits,YHits,Type;

  // Get tree from File
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
  // Get all data from tree
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
  return vv_ReturnVect;
}

// Read LogLilihood TTrees with branch names "_ll_"+Type+"<VARIABLE>" where Type is either "nnbar" or "neut"
vector <vector<Double_t>> ReadTreeLL(string RootFileName,string Type){

  // Initilise float variables, vector of floats and vector of vector of floats
  Double_t TDCstdU,TDCstdV,TDCstdY,TDCiqrU,TDCiqrV,TDCiqrY,ADCampU,ADCampV,ADCampY,WFintU,WFintV,WFintY,UHits,VHits,YHits;
  vector <Double_t> v_TDCstdU,v_TDCstdV,v_TDCstdY,v_TDCiqrU,v_TDCiqrV,v_TDCiqrY,v_ADCampU,v_ADCampV,v_ADCampY,v_WFintU,v_WFintV,v_WFintY,v_UHits,v_VHits,v_YHits;
  vector < vector < Double_t > > vv_ReturnVect;

  // Get the TTree from the file
  TFile *TRootFile = new TFile(RootFileName.c_str());  
  TTree *FileTree = (TTree*)TRootFile->Get("ch_tree");
 
  // Set Branch Addresses to floats
  FileTree->SetBranchAddress(("_ll_"+Type+"hnu").c_str(),&UHits);
  FileTree->SetBranchAddress(("_ll_"+Type+"hnv").c_str(),&VHits);
  FileTree->SetBranchAddress(("_ll_"+Type+"hny").c_str(),&YHits);
  FileTree->SetBranchAddress(("_ll_"+Type+"tsu").c_str(), &TDCstdU);
  FileTree->SetBranchAddress(("_ll_"+Type+"tsv").c_str(), &TDCstdV);
  FileTree->SetBranchAddress(("_ll_"+Type+"tsy").c_str(), &TDCstdY);
  FileTree->SetBranchAddress(("_ll_"+Type+"tiu").c_str(), &TDCiqrU);
  FileTree->SetBranchAddress(("_ll_"+Type+"tiv").c_str(), &TDCiqrV);
  FileTree->SetBranchAddress(("_ll_"+Type+"tiy").c_str(), &TDCiqrY);
  FileTree->SetBranchAddress(("_ll_"+Type+"aau").c_str(), &ADCampU);
  FileTree->SetBranchAddress(("_ll_"+Type+"aav").c_str(), &ADCampV);
  FileTree->SetBranchAddress(("_ll_"+Type+"aay").c_str(), &ADCampY);
  FileTree->SetBranchAddress(("_ll_"+Type+"wfu").c_str(), &WFintU);
  FileTree->SetBranchAddress(("_ll_"+Type+"wfv").c_str(), &WFintV);
  FileTree->SetBranchAddress(("_ll_"+Type+"wfy").c_str(), &WFintY);
  // Get the data from the tree
  Long64_t EntryNumber = FileTree->GetEntries();
  for(unsigned int j=0; j!=EntryNumber; ++j){ // Loop over tree entries
    // Fill vectors of floats with each branch address data member.
    FileTree->GetEntry(j);
    v_UHits.push_back(UHits);
    v_VHits.push_back(VHits);
    v_YHits.push_back(YHits);
    v_TDCstdU.push_back(TDCstdU);
    v_TDCstdV.push_back(TDCstdV);
    v_TDCstdY.push_back(TDCstdY);
    v_TDCiqrU.push_back(TDCiqrU);
    v_TDCiqrV.push_back(TDCiqrV);
    v_TDCiqrY.push_back(TDCiqrY);
    v_ADCampU.push_back(ADCampU);
    v_ADCampV.push_back(ADCampV);
    v_ADCampY.push_back(ADCampY);
    v_WFintU.push_back(WFintU);
    v_WFintV.push_back(WFintV);
    v_WFintY.push_back(WFintY);
  }
  // Push back data read in from tree onto vector of vectors. Columns of 2D vector correspond to events 1,2,3,...,N 
  // While rows correspond to the different cut variables.
  vv_ReturnVect.push_back(v_UHits);   vv_ReturnVect.push_back(v_VHits);   vv_ReturnVect.push_back(v_YHits);
  vv_ReturnVect.push_back(v_TDCstdU); vv_ReturnVect.push_back(v_TDCstdV); vv_ReturnVect.push_back(v_TDCstdY);
  vv_ReturnVect.push_back(v_TDCiqrU); vv_ReturnVect.push_back(v_TDCiqrV); vv_ReturnVect.push_back(v_TDCiqrY);
  vv_ReturnVect.push_back(v_ADCampU); vv_ReturnVect.push_back(v_ADCampV); vv_ReturnVect.push_back(v_ADCampY);
  vv_ReturnVect.push_back(v_WFintU);  vv_ReturnVect.push_back(v_WFintV);  vv_ReturnVect.push_back(v_WFintY); 
  return vv_ReturnVect;
}

// Function to calculate the mean
Double_t CalcMean(const vector <Double_t> &Vec){
  Double_t Mean = 0.0;
  for (unsigned int i =0; i < Vec.size(); i++){
    Mean += Vec[i];
  }
  return (Mean/Vec.size());
}

/* REDUNDANT METHOD
// Return the sum of the array Data after having performed a cut on it determining it is within the limits Lim1 < Data[i] < Lim2
Float_t TheCut(Float_t *Data, Size_t VectSize, Float_t Lim1, Float_t Lim2){
  //std::cout << "\nIn Optimise\n Lim1 is : "<< Lim1 << ", Lim2 is: " << Lim2;
  //if( Lim2 < (Lim1+1) ){ std::cout << "\n\nERROR!! Limits Equal or Limits larger than vector size!\n"; exit(1); }

  Float_t DataSum(0.0);
  for(unsigned int i(0); i < VectSize; i ++){ 
    if( Lim1 < Data[i] && Data[i] < Lim2 ){ DataSum += abs(Data[i]); }
  }
  
  return DataSum; // Makes this a probability
};
*/

// Return the sum of the array Data after having performed a cut on it determining it is within the limits Lim1 < Data[i] < Lim2
Float_t TheCutV2(Float_t *Data, Size_t VectSize, Float_t Lim1, Float_t Lim2){
  //std::cout << "\nIn Optimise\n Lim1 is : "<< Lim1 << ", Lim2 is: " << Lim2;
  //if( Lim2 < (Lim1+1) ){ std::cout << "\n\nERROR!! Limits Equal or Limits larger than vector size!\n"; exit(1); }

  Float_t DataSum(0.0);
  for(unsigned int i(0); i < VectSize; i ++){ 
    if( Lim1 < Data[i] && Data[i] < Lim2 ){ DataSum ++; }
  }
  
  return DataSum; // Makes this a probability
};



// Optimises a cut on two arrays to remove as much of neutbins as possible whilst maintaining over a certain
// Confidence percentage of nnbarbins events.
Float_t *OptimiseArrays(Float_t *nnbarbins, Float_t *neutbins, Size_t NnbarVectSize, Size_t NeutVectSize,Float_t Confidence, Float_t Accuracy,Float_t BotLim, Float_t UpLim){
  Float_t B(0.0), Beta(0.0), Alpha(1.0), MinAlpha(1000000000000000000.0), LowLim(0.0), HighLim(0.0), NormFact(0.0);

  for(unsigned int i(0); i < NnbarVectSize; i ++){ NormFact += abs(nnbarbins[i]); }
  std::cout << "\n\nConfidence = " << Confidence << ", Accuracy = " << Accuracy << ", BotLim = " << BotLim << ", UpLim = " << UpLim;
  for(Float_t k(BotLim); k < UpLim; k += Accuracy){ // Loop over lower limit
    for(Float_t j(BotLim); j < UpLim; j += Accuracy){ // Loop over upper limit
      if (k != j && k < j){ // If limits aren't equal and lower limit is smaller than upper limit

        //if( k > BotLim) B = TheCutV2(nnbarbins,NnbarVectSize,BotLim,k); // Lower part, IS NNBar Outside OF Critical Region
        //Beta = (B + TheCutV2(nnbarbins,NnbarVectSize,j,UpLim))/NnbarVectSize; // Probability IS NNBar Outside OF Critical Region
        //std::cout << "\n" << k << " to " << j << ". 1-Beta is : " << 1-Beta;// << ", B is : " << B;

        if( TheCutV2(nnbarbins,NnbarVectSize,k,j)/NnbarVectSize > Confidence ){ 
          Alpha = TheCutV2(neutbins,NeutVectSize,k,j); // Arbritary Probability IS NEUTRINO in Critical Region = ( Is not Nnbar in Critical region)
          //std::cout << "\nAlpha is : " << Alpha;
          if( Alpha < MinAlpha ){
            MinAlpha = Alpha; LowLim = k; HighLim = j; // Store indices and set new MinAlpha
            std::cout << "\nMinAlpha is: " << MinAlpha << ", LowLim = " << LowLim << ", HighLim = " << HighLim;
          }
        }
      }
    }
    //cin.get();
  }
  //cin.get();
  Float_t *ReturnFloatArr = new Size_t[2]; // Set the return float 
  ReturnFloatArr[0] = LowLim; ReturnFloatArr[1] = HighLim; // Return array contains best cuts
  return ReturnFloatArr;
}

// Calculate the Total Sum of a vector of Double_t
float TotalSum(const vector <Double_t> &Vect){
  Float_t ReturnSum(0.0);
  for(unsigned int i(0); i < Vect.size(); i ++) ReturnSum += Vect[i];
  return ReturnSum;
}

// Add together two vectors of Double_t's of the same length
Float_t *AddVects(vector <Double_t> &Vect1,const vector <Double_t> &Vect2){
  assert(Vect1.size() == Vect2.size()); // Check they are the same length
  Float_t* ReturnVect = new Float_t[Vect1.size()];
  for(unsigned int i(0); i < Vect1.size(); i ++) ReturnVect[i] = Vect1[i] + Vect2[i];
  return ReturnVect;
}

Float_t *AddVectToArray(vector <Double_t> &Vect,Float_t *Array){
  Float_t* ReturnArray = new Float_t[Vect.size()];
  for(unsigned int i(0); i < Vect.size(); i ++) ReturnArray[i] = Vect[i] + Array[i];
  return ReturnArray;
}


void SetLimitsV2(){
  // Define string arrays of the different variables
  string NeutCutVariables[14] = { "_ll_neuthnu","_ll_neuthnv","_ll_neuthny","_ll_neuttsu","_ll_neuttsv","_ll_neuttsy","_ll_neuttiu","_ll_neuttiv","_ll_neuttiy","_ll_neutaau","_ll_neutaav","_ll_neutaay","_ll_neutwfu","_ll_neutwfv","_ll_neutwfy" };
  string NnbarCutVariables[14] = { "_ll_nnbarhnu","_ll_nnbarhnv","_ll_nnbarhny","_ll_nnbartsu","_ll_nnbartsv","_ll_nnbartsy","_ll_nnbartiu","_ll_nnbartiv","_ll_nnbartiy","_ll_nnbaraau","_ll_nnbaraav","_ll_nnbaraay","_ll_nnbarwfu","_ll_nnbarwfv","_ll_nnbarwfy" };
  
  gROOT->ProcessLine(".x lhcbStyle.C"); // Using lhcb Style file
  
  // Get pointers to the two TTree objects
  TFile *TRootNNbarFile = new TFile("NnbarLogLikelihoods_kam.root");
  TTree *NNbarTree = (TTree*)TRootNNbarFile->Get("ch_tree");
  TFile *TRootNeutFile = new TFile("NeutLogLikelihoods_kam.root");
  TTree *NeutTree = (TTree*)TRootNeutFile->Get("ch_tree");
  
  NeutTree->SetLineColor(kRed);
  NNbarTree->SetLineColor(kBlue);
  
  // Read in the data from the Logliklihood root files.
  vector <vector<Double_t>> LLNeut = ReadTreeLL("NeutLogLikelihoods_kam.root","neut"); 
  vector <vector<Double_t>> LLNnbar = ReadTreeLL("NnbarLogLikelihoods_kam.root","nnbar");
  
  Size_t NeutVectSize = LLNeut[0].size();
  Size_t NnbarVectSize = LLNnbar[0].size();
  Size_t *Limits = new Float_t[2];
  
  vector <Double_t> LowerLimits, UpperLimits;
  Float_t *NeutArray; Float_t *NnbarArray;
  NeutArray = new Float_t[NeutVectSize];
  NnbarArray = new Float_t[NnbarVectSize];
  
  for(unsigned int j(0); j < LLNeut.size(); j ++){
    // Set values of arrays to values in Logliklihood root files.
    for(unsigned int i(0); i < NeutVectSize; i ++){ NeutArray[i] = LLNeut[j][i]; } 
    for(unsigned int i(0); i < NnbarVectSize; i++){ NnbarArray[i] = LLNnbar[j][i]; }
    
    // Get the limits through optimisation function
    Limits = OptimiseArrays(NnbarArray,NeutArray,NnbarVectSize,NeutVectSize,0.95,0.2,-5.0,5.0);
    std::cout << "\nLimits are: " << Limits[0] << ", " << Limits[1];
    
    std::cout << "\nNeutrinos Removed = " << 1-TheCutV2(NeutArray,NeutVectSize,Limits[0],Limits[1])/NeutVectSize << "%";
    std::cout << "\nNnbar Removed = " << 1-TheCutV2(NnbarArray,NnbarVectSize,Limits[0],Limits[1])/NnbarVectSize << "%\n";
    
    // Push back lower and upper limits onto vectors of lower and upper limits
    LowerLimits.push_back(Limits[0]);
    UpperLimits.push_back(Limits[1]);
  }

  std::ofstream MyFile; // Create output file stream with name AnalyserDataOutput.txt
  MyFile.open("OutputSetLimits.txt");

  vector <Double_t> Means, RemovedNeutrinos, RemovedNnbars;
  MyFile << "\nNOW OUTPUTTING RAW DATA\n";

  for(unsigned int p(0); p < LLNnbar.size(); p ++){
    Means.push_back( CalcMean(LLNnbar[p]) );
  } 
    
  RemovedNnbars = DoTheCut(LLNnbar,LowerLimits,UpperLimits,NnbarVectSize);
  RemovedNeutrinos = DoTheCut(LLNeut,LowerLimits,UpperLimits,NeutVectSize);
  
  for(unsigned int p(0); p < LLNnbar.size(); p ++){
    MyFile << "\n\n" << NeutCutVariables[p];
    MyFile << "\nMeans = " << Means[p] << ", L = " << LowerLimits[p] << ", U = " << UpperLimits[p];    
    MyFile << "\nRemovedNeutrinos: " << RemovedNeutrinos[p] << " : " << 100*(1-RemovedNeutrinos[p]/NeutVectSize) << "%";
    MyFile << "\nRemovedNnbar: " << RemovedNnbars[p] << " : " << 100*(1-RemovedNnbars[p]/NnbarVectSize) << "%";

    TCanvas *c2 = new TCanvas((NnbarCutVariables[p]).c_str(),"A NAME",900,600);
    NNbarTree->Draw(NnbarCutVariables[p].c_str(),"","E1HIST",10000000000,0);
    NeutTree->Draw(NeutCutVariables[p].c_str(),"","E1HISTSAME",10000000000,0);

    TLine *XMinSig = new TLine(LowerLimits[p],0,LowerLimits[p],500); XMinSig->Draw(); // c2->GetUymax();
    TLine *XMaxSig = new TLine(UpperLimits[p],0,UpperLimits[p],500); XMaxSig->Draw();
    c2->Update(); c2->Modified();
  }
  
  // Find Max sum difference
  float SumDiff(0.0);
  vector < float > DiffVec,OldVec; vector <int> Indexes;
  for(unsigned int p(0); p < LLNeut.size(); p ++){
    if( p != (LLNeut.size()-2)){
      //std::cout << p << " : NeutSum = " << TotalSum(LLNeut[p]) << ", NnbarSum = " << TotalSum(LLNnbar[p]) << "\n";
      SumDiff = TotalSum(LLNnbar[p]) - TotalSum(LLNeut[p]);
      DiffVec.push_back( SumDiff );
      //std::cout << "Difference = " << SumDiff << "\n";      
    }
  }
  OldVec = DiffVec;
  sort(DiffVec.begin(),DiffVec.end());
  // Find indexes of the sorted differences.
  for(unsigned int i(0); i < DiffVec.size(); i ++){
    for(unsigned int j(0); j < DiffVec.size(); j ++){
      if( OldVec[j] == DiffVec[i] ) Indexes.push_back(j);
    }
    // std::cout << "\n" << Indexes[i] << " : " << DiffVec[i];
  }
  
  // At this point, LLNeut[Indexes[13]] and LLNnbar[Indexes[13]] Gives the best discrimination, followed by Indexes[12], Indexes[11] etc.

  // TOM - If you want to test adding multiple LogLiklihoods see below:
  // To just test the two best ones, set "i < 13" to "i < 0" and "j < 13" to "j < 1"
  for(unsigned int i(0); i < 13; i ++ ){
    for( unsigned int j = i; j < 13; j ++){
      std::cout << "\nIn Loop i = " << i << ", j = " << j << "\n";
      // Add together different log liklihoods, Last two elements of vector are best
      Float_t *DoubleLLneut = AddVects(LLNeut[Indexes[13-i]],LLNeut[Indexes[12-j]]);
      Float_t *DoubleLLnnbar = AddVects(LLNnbar[Indexes[13-i]],LLNnbar[Indexes[12-j]]);
 
      MyFile << "\n\n" << NeutCutVariables[Indexes[13-i]] << " AND " << NeutCutVariables[Indexes[12-j]];
  
      Limits = OptimiseArrays(DoubleLLnnbar,DoubleLLneut,NnbarVectSize,NeutVectSize,0.95,0.5,-5.0,8.0);
  
      MyFile << "\nLimits are: " << Limits[0] << ", " << Limits[1];
   
      MyFile << "\nNeutrinos Removed = " << 100*(1-TheCutV2(DoubleLLneut,NeutVectSize,Limits[0],Limits[1])/NeutVectSize) << "%";
      MyFile << "\nNnbar Removed = " << 100*(1-TheCutV2(DoubleLLnnbar,NeutVectSize,Limits[0],Limits[1])/NeutVectSize) << "%\n";
    }
  }
  // TO TEST THE 3 (or more) BEST VARIABLES:
  std::cout << "\nNow Testing 3 variables\n";
  std::cout << "\nWARNING: Array length assumed to be the same as Vector size";

  // Add Together Arrays and Vectors
  Float_t *ManyLLneut = AddVects(LLNeut[Indexes[13]],LLNeut[Indexes[12]]);
  Float_t *ManyLLneut = AddVectToArray( LLNeut[Indexes[11]],ManyLLneut);
  Float_t *ManyLLnnbar = AddVects(LLNnbar[Indexes[13]],LLNnbar[Indexes[12]]);
  Float_t *ManyLLnnbar = AddVectToArray(LLNnbar[Indexes[11]],ManyLLnnbar);

  MyFile << "\n\n" << NeutCutVariables[Indexes[13]] << " AND " << NeutCutVariables[Indexes[12]] << " AND " << NeutCutVariables[Indexes[11]];
 
  Limits = OptimiseArrays(ManyLLnnbar,ManyLLneut,NnbarVectSize,NeutVectSize,0.95,0.5,-5.0,8.0);
  
  MyFile << "\nLimits are: " << Limits[0] << ", " << Limits[1];
   
  MyFile << "\nNeutrinos Removed = " << 100*(1-TheCutV2(ManyLLneut,NeutVectSize,Limits[0],Limits[1])/NeutVectSize) << "%";
  MyFile << "\nNnbar Removed = " << 100*(1-TheCutV2(ManyLLnnbar,NnbarVectSize,Limits[0],Limits[1])/NnbarVectSize) << "%\n";
  

  
  // Male a plot of the most recently calculated combination of logliklihoods 
  // (Move this into double for loop and rename histograms to plot a histogram for each different variable.
  TCanvas *c2 = new TCanvas("Combined","A NAME",900,600);
  TH1F* NnbarHisto = new TH1F("neuthist","neuthist",50,-10,10);
  TH1F* NeutHisto = new TH1F("nnbarhist","nnbarhist",50,-10,10);
  for(unsigned int i(0); i < NnbarVectSize; i ++){ NnbarHisto->Fill(ManyLLnnbar[i]);}
  for(unsigned int i(0); i < NeutVectSize; i ++){ NeutHisto->Fill(ManyLLneut[i]); }
  NnbarHisto->Draw();
  NeutHisto->Draw("same");
  TLine *XMinSig = new TLine(Limits[0],0,Limits[0],2000); XMinSig->Draw(); // c2->GetUymax();
  TLine *XMaxSig = new TLine(Limits[1],0,Limits[1],2000); XMaxSig->Draw();
  
  std::cout << "\nLimits are: LowLim =" << Limits[0] << ", HighLim = " << Limits[1] <<"\n";


  c2->Update(); c2->Modified();
  MyFile.close(); LLNeut.clear(); LLNnbar.clear();
}
