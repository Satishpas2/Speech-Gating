//
//  DeclareRhythmicMeasure.cpp
//  Created by MANGESH Sonawane
//  Copyright © 2019 Signum Audio. All rights reserved.
//
#include "DeclareRhythmicMeasure.h"
#include<iostream>
#include <stdio.h>
using namespace std;

    //Constructor
    DeclareRhythmicMeasure:: DeclareRhythmicMeasure(int blockNumbers)
    {
    variance= new double [blockNumbers];
    zeroMeanSequence=new double[blockNumbers];
    autoCorrelationOfZeroMeanSequence=new double[blockNumbers];
    }
    
    //Destructor
     DeclareRhythmicMeasure::~DeclareRhythmicMeasure()
    {
        delete []variance;
        delete []zeroMeanSequence;
        delete []autoCorrelationOfZeroMeanSequence;
    }
    
