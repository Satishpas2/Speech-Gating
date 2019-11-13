//
//  DeclareSkewOfRegressiveLineOfBestFitThroughEstimatedSpectralPowerDensity.h
//  Created by MANGESH Sonawane
//  Copyright © 2019 Signum Audio. All rights reserved.
//
#include <iostream>
#include <stdio.h>
using namespace std;

class DeclareSkewOfRegressiveLineOfBestFitThroughEstimatedSpectralPowerDensity
{
public:
    double **Xdb; //in dB log-domain
    double *regressiveCoefficient;// it is the 'G', regressive coefficient
//Constructor
    DeclareSkewOfRegressiveLineOfBestFitThroughEstimatedSpectralPowerDensity(int blockNumbers);
    
    //Destructor
    ~DeclareSkewOfRegressiveLineOfBestFitThroughEstimatedSpectralPowerDensity();
};
