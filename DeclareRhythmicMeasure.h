//
//  DeclareRhythmicMeasure.h
//  Created by MANGESH Sonawane
//  Copyright © 2019 Signum Audio. All rights reserved.
//
#include<iostream>
#include <stdio.h>
using namespace std;

class DeclareRhythmicMeasure
{
public:
       //variance of sample x in block m
    double *variance;
      //Zero mean Sequence
    double *zeroMeanSequence;
        //AutoCorrelation of ZeroMean Sequence
    double *autoCorrelationOfZeroMeanSequence;

    //Constructor
    DeclareRhythmicMeasure(int blockNumbers);
    
    //Destructor
    ~DeclareRhythmicMeasure();
    
};
