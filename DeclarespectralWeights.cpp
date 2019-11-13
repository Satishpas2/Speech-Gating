//
//  DeclarespectralWeights.cpp
//  Created by MANGESH Sonawane
//  Copyright © 2019 Signum Audio. All rights reserved.
//
#include "DeclarespectralWeights.h"
#include <iostream>
#include <stdio.h>
using namespace std;
//Constructor
    DeclarespectralWeights::DeclarespectralWeights(int blockNumbers)
    {
        Xdb=new double *[blockNumbers];
        maximumLogDomainPowerSpectrum=new double [blockNumbers];
        spectralWeight=new double[blockNumbers];
    }
    
    //Destructor
    DeclarespectralWeights::~DeclarespectralWeights()
    {
        delete Xdb;
        delete maximumLogDomainPowerSpectrum;
        delete spectralWeight;
    }
