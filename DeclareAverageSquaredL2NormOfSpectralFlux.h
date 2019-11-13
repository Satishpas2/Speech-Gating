//
//  DeclareAverageSquaredL2NormOfSpectralFlux.h
//  Created by MANGESH Sonawane
//  Copyright © 2019 Signum Audio. All rights reserved.
//

#include <stdio.h>
#include <iostream>
using namespace std;
//Feature 1
class DeclareAverageSquaredL2NormOfSpectralFlux
{
public:
     // 'W' is  the weight for the current block m. (EQUATION 2)
    double * Weight;
    //'l2Norm'is l2-norm of the weighted spectral flux for block (EQUATION 3)
    double *l2Norm;
   
    //Constructor
 
    DeclareAverageSquaredL2NormOfSpectralFlux(int blockNumbers);
    
    //Destructor
    ~DeclareAverageSquaredL2NormOfSpectralFlux();
};
