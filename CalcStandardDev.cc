#include <vector>
#include "TMath.h"

float CalcStandardDev(const vector <float> &Vec){
  float Mean = CalcMean(Vec);
  float sum(0.0), var(0.0);
  int rem(0);
       
  for (unsigned int i =0; i < Vec.size(); i++){
    if(TMath::IsNaN(Vec[i])){ rem += 1;
    }else sum += pow((Vec[i] - Mean),2); 
  }
  var = sqrt(sum/(Vec.size()-rem));
  return var;
}

