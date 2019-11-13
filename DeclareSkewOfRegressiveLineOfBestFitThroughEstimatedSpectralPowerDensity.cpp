//
//  DeclareSkewOfRegressiveLineOfBestFitThroughEstimatedSpectralPowerDensity.cpp
//  Created by MANGESH Sonawane
//  Copyright © 2019 Signum Audio. All rights reserved.
//
#include "DeclareSkewOfRegressiveLineOfBestFitThroughEstimatedSpectralPowerDensity.h"
#include <iostream>
#include <stdio.h>
using namespace std;


//Constructor
    DeclareSkewOfRegressiveLineOfBestFitThroughEstimatedSpectralPowerDensity::DeclareSkewOfRegressiveLineOfBestFitThroughEstimatedSpectralPowerDensity(int blockNumbers)
    {
     Xdb=new double *[blockNumbers];
    regressiveCoefficient=new double[blockNumbers];
    }
    
    //Destructor
    DeclareSkewOfRegressiveLineOfBestFitThroughEstimatedSpectralPowerDensity::~DeclareSkewOfRegressiveLineOfBestFitThroughEstimatedSpectralPowerDensity()
    {
        delete []Xdb;
    }
    
