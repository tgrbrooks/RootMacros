#include <vector>
#include "TMath.h"

float CalcMean(const vector <float> &Vec){
  float Mean = 0.0; int rem(0);

  for (unsigned int i =0; i < Vec.size(); i++){
    if( TMath::IsNaN(Vec[i]) ){ rem += 1;
    }else Mean += Vec[i];
  }
  return ( Mean/(Vec.size() - rem));
}

