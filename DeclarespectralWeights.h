//
//  DeclarespectralWeights.h
//  Created by MANGESH Sonawane
//  Copyright © 2019 Signum Audio. All rights reserved.
//
#include <iostream>
#include <stdio.h>
using namespace std;

class DeclarespectralWeights
{
public:
    //spectral Power Density in dB log-domain
    double **Xdb;
     //Maximum log domain power spectrum value equation 16
    double *maximumLogDomainPowerSpectrum;
     //Spectral Weight
    double *spectralWeight;
       
    //Constructor
    DeclarespectralWeights(int blockNumbers);
    
    //Destructor
    ~DeclarespectralWeights();
};
