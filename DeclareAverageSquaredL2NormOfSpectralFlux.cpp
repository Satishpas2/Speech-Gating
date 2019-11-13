//
//  DeclareAverageSquaredL2NormOfSpectralFlux.cpp
//  Created by MANGESH Sonawane
//  Copyright © 2019 Signum Audio. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include "DeclareAverageSquaredL2NormOfSpectralFlux.h"
using namespace std;
//Feature 1

//Constructor
 
  DeclareAverageSquaredL2NormOfSpectralFlux::DeclareAverageSquaredL2NormOfSpectralFlux(int blockNumbers)
    {
        Weight=new double[blockNumbers];
        l2Norm=new double[blockNumbers];
    }
    
    //Destructor
    DeclareAverageSquaredL2NormOfSpectralFlux::~DeclareAverageSquaredL2NormOfSpectralFlux()
    {
        delete []Weight;
        delete []l2Norm;
    }

