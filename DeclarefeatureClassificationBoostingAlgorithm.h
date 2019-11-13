//
//  DeclarefeatureClassificationBoostingAlgorithm.h
//  Created by MANGESH Sonawane
//  Copyright © 2019 Signum Audio. All rights reserved.
//

#include <stdio.h>
#include <iostream>
using namespace std;
class DeclarefeatureClassificationBoostingAlgorithm
{
public:
   double Threshold[20]={5.69921,0.845003,5.79869,0.220335,1.456317, 5.871387, 0.221293,-10.8799,5.797844, 1.408993, 5.803905, 0.206998, 5.748194, 0.201943,0.862561,5.796515, 1.451334, 0.201946, 5.735963, 5.795558};
    double weightingCoefficient[20]={0.921095,-0.76701,0.545198,-0.54092,0.427891,0.395285, -0.36794,0.279396,0.289932,0.160819,0.121731, -0.13833,0.137475, -0.11858,-0.07838,0.06853,0.071846,-0.05518,0.063226, 0.048795};
    int featureVectorIndexi_j[20]={1,5,1,6,4,3,7,2,1,4,1,7,3,6,5,1,4,7,3,1};
     // the weighted “weak” classifier output for the jth interim “weak” classifier stage.
    double*weightedWeakClassifier;
//Constructor
    DeclarefeatureClassificationBoostingAlgorithm();
    //Destructor
    ~DeclarefeatureClassificationBoostingAlgorithm();
};
