//
//  DeclarePauseCount.cpp
//  Created by MANGESH Sonawane
//  Copyright © 2019 Signum Audio. All rights reserved.
//
#include "DeclarePauseCount.h"
#include <iostream>
#include <stdio.h>
using namespace std;

    //Constructor
   DeclarePauseCount:: DeclarePauseCount(int blockNumber)
    {
        powerOfEachBlock= new double[blockNumber];
    }
    //Destructor
    DeclarePauseCount::~DeclarePauseCount()
   {
        delete []powerOfEachBlock;
    }

