//
//  DeclarefeatureClassificationBoostingAlgorithm.cpp
//  Created by MANGESH Sonawane
//  Copyright © 2019 Signum Audio. All rights reserved.
//
#include"DeclarefeatureClassificationBoostingAlgorithm.h"
#include <stdio.h>
#include <iostream>
using namespace std;
//Constructor
    DeclarefeatureClassificationBoostingAlgorithm::DeclarefeatureClassificationBoostingAlgorithm()
    {
    weightedWeakClassifier=new double [20];
    }
    //Destructor
    DeclarefeatureClassificationBoostingAlgorithm::~DeclarefeatureClassificationBoostingAlgorithm()
    {
        delete []weightedWeakClassifier;
    }
