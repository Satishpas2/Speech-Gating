//
//  DeclareLongRhythmicMeasure.cpp
//  Created by MANGESH Sonawane
//  Copyright © 2019 Signum Audio. All rights reserved.
//
#include "DeclareLongRhythmicMeasure.h"
#include<iostream>
#include <stdio.h>
using namespace std;

//Constructor
    DeclareLongRhythmicMeasure::DeclareLongRhythmicMeasure(int blockNumbers)
    {
        int doubleBlockNumber= 2*blockNumbers;
        totalSpectralWeights=new double[doubleBlockNumber];
        autoCorrelationScoreValue=new double[doubleBlockNumber];
    }
    //Destructor
   DeclareLongRhythmicMeasure::~DeclareLongRhythmicMeasure()
    {
        delete []totalSpectralWeights;
        delete []autoCorrelationScoreValue;
    }

