//
//  DeclareLongRhythmicMeasure.h
//  Created by MANGESH Sonawane
//  Copyright © 2019 Signum Audio. All rights reserved.
//
#include<iostream>
#include <stdio.h>
using namespace std;
class DeclareLongRhythmicMeasure
{
public:
    //Total spectral weights includes values of spectral weights from previous frame and current frame
    double *totalSpectralWeights;
    //AutoCorrelation ScoreValue
    double *autoCorrelationScoreValue;
//Constructor
    DeclareLongRhythmicMeasure(int blockNumbers);
    //Destructor
    ~DeclareLongRhythmicMeasure();
};
